#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(__x86_64__) || defined(_M_X64)
#include <cpuid.h>
#include <immintrin.h>
#endif

// Next-nearest-neighbour exclusion engine.
//
// The hop rate across a bond is allowed to depend on the two sites flanking it
// as well as on the bond itself.  Writing a bond as the pair of sites (b, b+1),
// the rate of a hop across it is a function of the flanking pair
//
//     (a, d) = (n_{b-1}, n_{b+2}),
//
// which takes four values.  Every allowed bond therefore belongs to one of four
// rate classes per direction, and within a class all bonds share one rate --
// the property the whole selection scheme relies on.  The Katz-Lebowitz-Spohn
// model is the instance provided here, but any model whose rates depend on the
// next-nearest occupations is expressed by changing the eight class weights.
namespace monitored_exclusion {

// Selecting the n-th allowed bond of a mask word is the hottest single step of
// an event: the portable form walks past n set bits one at a time, a chain of
// dependent operations averaging about a quarter of the word's width at half
// filling.  BMI2's PDEP does it in one instruction, by depositing a single bit
// into the n-th set position of the mask.
//
// It cannot simply be assumed.  PDEP is fast on Intel from Haswell and on AMD
// from Zen 3, but on Zen 1 and Zen 2 it is microcoded and takes hundreds of
// cycles -- far worse than the loop it would replace.  So the instruction is
// compiled in behind a target attribute, and chosen once at start-up from what
// the processor reports.
namespace {

bool detect_usable_pdep() {
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
    if (!__builtin_cpu_supports("bmi2")) return false;
    unsigned eax = 0, ebx = 0, ecx = 0, edx = 0;
    if (!__get_cpuid(0, &eax, &ebx, &ecx, &edx)) return false;
    if (ebx != 0x68747541u) return true;          // "Auth", i.e. AuthenticAMD
    if (!__get_cpuid(1, &eax, &ebx, &ecx, &edx)) return false;
    const unsigned base_family = (eax >> 8) & 0xFu;
    const unsigned family = base_family == 0xFu ? base_family + ((eax >> 20) & 0xFFu) : base_family;
    return family >= 0x19u;                       // Zen 3 and later
#else
    return false;
#endif
}
// In the previous function, eax packs several fields into one word, the family among them:
    //
    //     bits 27:20  extended family      bits 11:8  family
    //     bits 19:16  extended model       bits  7:4  model
    //                                      bits  3:0  stepping
    //
// The family field is four bits wide, so it saturates at 15.  0xF is therefore an escape value meaning "the real family is larger, add the extended field to it".  Zen 3 is 0xF + 0xA = 0x19, Zen 2 is 0xF + 0x8 = 0x17, so the comparison below keeps Zen 3 and later and drops the microcoded Zen 1 and Zen 2.

const bool usable_pdep = detect_usable_pdep();

#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("bmi2"))) //This attribute is passed only to the declaration right below it
inline std::size_t select_in_word_pdep(std::uint64_t mask, std::size_t ordinal) {
    return static_cast<std::size_t>(std::countr_zero(_pdep_u64(std::uint64_t{1} << ordinal, mask)));
}
#endif

// Position of the (ordinal+1)-th set bit, counting from the least significant.
inline std::size_t select_in_word(std::uint64_t mask, std::size_t ordinal) {
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
    if (usable_pdep) return select_in_word_pdep(mask, ordinal);
#endif
    while (ordinal-- > 0) mask &= mask - 1;       // clear the lowest set bit
    return static_cast<std::size_t>(std::countr_zero(mask));
}

}  // namespace

// Flanking patterns, indexed as 2*a + d for (a, d) = (n_{b-1}, n_{b+2}).
inline constexpr int classes = 4;

enum class DynamicsKind {Heap, Uniformized, Direct};

// Every field below is overwritten by parse_parameters(), which reads all of
// them from the dictionary the Python wrapper builds and raises KeyError on a
// missing one, so none of these initializers is ever observed through the
// normal entry point.  They are kept so that a field added to the struct and
// forgotten in parse_parameters() degrades to a defined value instead of an
// indeterminate one, and they mirror the resolved defaults of the NNN
// constructor in Python_Cpp_Interface/NNN_class.py, which is the single source
// of truth for what an unspecified parameter means.
struct Parameters {
    std::size_t L = 20;
    std::size_t M = 100;
    double tmax = 20.0;
    double E = 0.0;
    double s = 0.0;
    double k = 0.0;

    std::string model = "KLS";
 
    // The KLS interaction, matching the constructor's KLS defaults -- which
    // apply only to a KLS run, the model this struct also defaults to.  A
    // WASEP or NNN run zeroes both before they reach here.
    double epsilon = 0.6;
    double delta_kls = 0.0;
    // Explicit class weights, four per direction indexed as 2*a + d.  Read
    // only when model = "NNN", which requires both; the other two families
    // ignore whatever is here.
    std::vector<double> right_weights;
    std::vector<double> left_weights;
    double filling = 0.5;
    double alpha = 1.0;
    double gamma = 1.0;
    double delta = 1.0;
    double beta = 1.0;
    // Supplying only target_max implies target_min as its two-thirds power,
    // so the pair below is 10 and 10^(2/3).
    double target_min = 4.641588833612778;
    double target_max = 10.0;
    // Infinity means "no clock-based check"; a finite value adds one on top
    // of whatever the mechanism schedules for itself.
    double cloning_interval = std::numeric_limits<double>::infinity();
    double record_interval = -1.0;
    std::uint64_t seed = 12345;
    std::string initial = "alternating";
    std::string boundary = "open";
    std::size_t progress = 0;
    DynamicsKind dynamics = DynamicsKind::Direct;

    bool rejection_free() const { return dynamics != DynamicsKind::Uniformized; }
};

class Model {
public:
    explicit Model(const Parameters& p): L(p.L), periodic(p.boundary == "periodic"), links(periodic ? p.L : p.L - 1), E(p.E), s(p.s), k(p.k), k_scaled(p.k / (static_cast<double>(p.L) * static_cast<double>(p.L))), alpha(p.alpha), gamma(p.gamma), delta(p.delta), beta(p.beta),dq(1.0 / static_cast<double>(links)),p0(std::exp(E / static_cast<double>(links))),q0(std::exp(-E / static_cast<double>(links))),ps(std::exp((E + s) / static_cast<double>(links))),qs(std::exp(-(E + s) / static_cast<double>(links))) {
        // The model alone decides which rate description is read: NNN takes
        // the weight vectors, KLS builds its own from epsilon and delta_kls,
        // WASEP uses neither.  Whatever a family does not read is ignored.
        const bool use_weights = (p.model == "NNN");
        uniform = (p.model == "WASEP");
        if (uniform) {
            for (int c = 0; c < classes; ++c)
                right_weight[c] = left_weight[c] = 1.0;
            rate_model = "wasep";
            heaviest_weight = 1.0;
            return;
        }
        if (!use_weights && p.model != "KLS")
            throw std::invalid_argument("model must be WASEP, KLS or NNN");
        if (use_weights && (p.right_weights.size() != static_cast<std::size_t>(classes) || p.left_weights.size() != static_cast<std::size_t>(classes)))
            throw std::invalid_argument( "model NNN needs four right_weights and four left_weights");
        for (int a = 0; a < 2; ++a) {
            for (int d = 0; d < 2; ++d) {
                const int cls = 2 * a + d;
                if (use_weights) {
                    right_weight[cls] = p.right_weights[cls];
                    left_weight[cls] = p.left_weights[cls];
                    continue;
                }
                const double sum = static_cast<double>(a + d);
                const double difference = static_cast<double>(a - d);
                right_weight[cls] = 1.0 + p.delta_kls * (1.0 - sum) + p.epsilon * difference;
                left_weight[cls] = 1.0 + p.delta_kls * (1.0 - sum) - p.epsilon * difference;
            }
        }
        custom_weights = use_weights;
        rate_model = use_weights ? "custom" : "kls";
        heaviest_weight = 0.0;
        for (int c = 0; c < classes; ++c)
            heaviest_weight = std::max({heaviest_weight, right_weight[c],left_weight[c]});
    }

    // ``activity_right`` and ``activity_left`` are the class counts already
    // weighted by their rates, so the counting-field potential keeps exactly
    // the WASEP form.  For a bond of class c the untilted rate is c*p0 and the
    // tilted one c*ps, and the weighted counts absorb the factor c.
    double counting_potential(double activity_right, double activity_left) const {
        return activity_right * (ps - p0) + activity_left * (qs - q0);
    }

    double widest_possible_potential_spread() const {
        return static_cast<double>(links) * heaviest_weight* (std::abs(ps - p0) + std::abs(qs - q0))+ 2.0 * k_scaled * static_cast<double>(L);
    }

    double broad_jump_rate_estimate() const {
        return static_cast<double>(links) * heaviest_weight * (ps + qs) / 4.0 + (periodic ? 0.0 : (alpha + gamma + delta + beta) / 2.0);
    }

    std::size_t L;
    bool periodic;
    std::size_t links;
    double E, s, k;
    // The measurement coupling enters every weight as k/L^2; see the
    // monitoring section of docs/theory_algorithm_implementation.tex.  The raw k
    // is kept for the metadata and the k == 0 tests.
    double k_scaled;
    double alpha, gamma, delta, beta;
    double dq;
    double p0, q0, ps, qs;
    double right_weight[classes], left_weight[classes];
    double heaviest_weight;
    bool custom_weights = false;
    // True when all eight class rates are equal by construction, which lets the
    // engine keep one pair of bond masks instead of four.
    bool uniform = false;
    std::string rate_model = "kls";
};

struct Walker {
    std::size_t nR[classes] = {0, 0, 0, 0};
    std::size_t nL[classes] = {0, 0, 0, 0};
    double activity_right = 0.0;
    double activity_left = 0.0;
    int overlap = 0;
    double potential = 0.0;
    double log_weight = 0.0;
    double last_weight_time = 0.0;
    double current = 0.0;

    void update_weight(double time) {
        log_weight += potential * (time - last_weight_time);
        last_weight_time = time;
    }
};

struct Record {
    double time = 0.0;
    double cgf_cloning = 0.0;
    double mean_potential = 0.0;
    double accepted_current_rate = 0.0;
    double reference_current_rate = 0.0;
    double potential_spread = 0.0;
    double effective_sample_size = 0.0;
    std::uint64_t resampling_steps = 0;
    std::uint64_t resampling_checks = 0;
    std::uint64_t accepted_walker_events = 0;
    std::uint64_t rejected_walker_events = 0;
    std::vector<double> density_profile;
};

struct Event {
    double time;
    std::size_t process;
};

// Orders events so that std::make_heap/push_heap/pop_heap keep the *earliest*
// event at the front: the "greatest" element under this comparator is the one
// that happens first.
struct LaterEvent {
    bool operator()(const Event& a, const Event& b) const {
        if (a.time != b.time) return a.time > b.time;
        return a.process > b.process;
    }
};

// One reference-process event, stored so that a walker evolved on its own can
// replay it later.  The occupations are the reference's *before* the flip,
// which is all a walker needs to update its overlap without the live reference
// configuration, and which the reference itself no longer holds once it has
// run ahead.
struct ReferenceEventRecord {
    double time;
    std::size_t changed_a;
    std::size_t changed_b;
    bool before_a;
    bool before_b;
};

class Simulation;

// One event mechanism. The simulation owns exactly one of these and asks it to
// carry the population from the present time to the next scheduled instant;
// everything else -- the configurations, the weights, the cloning, the records
// -- belongs to the Simulation and is reached through the reference each
// implementation holds.
//
// Every implementation supplies two loops. ``advance_weighted`` is the general
// one. ``advance_unbiased`` is the s = k = 0 case, where all weights stay
// exactly one: no weight is integrated, no certificate is accumulated and no
// cloning can be triggered, so the loop reduces to plain trajectory sampling
// and none of that machinery is touched.
class Dynamics {
public:
    virtual ~Dynamics() = default;

    // Returns false when the weight certificate tripped before ``end_time``,
    // leaving the clock at the event that tripped it. Implementations without
    // a running certificate always return true.
    virtual bool advance_weighted(double end_time) = 0;
    virtual void advance_unbiased(double end_time) = 0;

    // Called after the population has been resampled, for whatever per-process
    // scheduling state the implementation keeps.
    virtual void after_resampling() {}

    // Called after a checkpoint restore, when the per-process event times are
    // already back in place and only derived structures need rebuilding.
    virtual void restore_schedule() {}

    // The interval until the next cloning examination.  A non-positive value
    // leaves the schedule to the running certificate.
    // Infinity means the mechanism schedules nothing itself and leaves the
    // examination to the running weight certificate.
    virtual double step_size() const { return std::numeric_limits<double>::infinity(); }

    // Called once the simulation's derived constants exist, before any block.
    virtual void prepare() {}

    // Scheduling state a mechanism carries between blocks and must therefore
    // survive a checkpoint.  Only the adaptive step of the direct mechanism
    // needs it; the others schedule from the certificate alone.
    virtual double schedule_state() const { return 0.0; }
    virtual void set_schedule_state(double /*value*/) {}

    // Whether this mechanism judges a finished block before it counts.
    virtual bool reviews_blocks() const { return false; }

    // Called with the exact log-weight range the block produced and how long
    // it ran.  Returns true when the block must be discarded and redone with a
    // shorter step.
    virtual bool review_block(double /*log_weight_range*/,double /*duration*/) { return false; }
};

std::unique_ptr<Dynamics> make_dynamics(Simulation& simulation, DynamicsKind kind);

class Simulation {
public:
    friend class HeapDynamics;
    friend class UniformizedDynamics;
    friend class DirectDynamics;

    explicit Simulation(Parameters parameters) : parameters_(std::move(parameters)), model_(parameters_), unbiased_(parameters_.s == 0.0 && parameters_.k == 0.0), words_((parameters_.L + 63) / 64), final_word_mask_(final_mask(parameters_.L)), final_bond_mask_(final_mask( (parameters_.L - 1) - ((parameters_.L - 1) / 64) * 64)), walkers_(parameters_.M), configs_(parameters_.M * words_, 0), masks_(parameters_.M * mask_vectors * words_, 0), reference_(words_, 0), reference_masks_(mask_vectors * words_, 0), rng_(parameters_.seed), event_times_(parameters_.M + 1, std::numeric_limits<double>::infinity()), clone_log_weights_(unbiased_ ? 0 : parameters_.M), clone_weights_(unbiased_ ? 0 : parameters_.M), spare_walkers_(unbiased_ ? 0 : parameters_.M), spare_configs_(unbiased_ ? 0 : parameters_.M * words_, 0), spare_masks_(unbiased_ ? 0 : parameters_.M * mask_vectors * words_, 0) {
        validate_parameters();
        dynamics_ = make_dynamics(*this, parameters_.dynamics);
        // Only the single-clock mechanisms can maintain a running bound on the
        // weight spread; direct dynamics move each walker on its own clock and
        // schedule cloning in advance instead, so they must not pay for one.
        track_certificate_ = parameters_.dynamics != DynamicsKind::Direct;
        choose_intervals();
        dynamics_->prepare();
        initialize_configuration();
        initialize_state();
        dynamics_->after_resampling();
        // Direct dynamics set their interval from the potentials, which exist
        // only once the walkers do; heap and uniformized are unaffected because
        // the clock has not moved.
        schedule_next_cloning_check();
    }

    std::vector<Record> run() {
        std::vector<Record> rows;
        if (!initial_row_emitted_) {
            rows.push_back(output_record());
            reset_recording_window();
            initial_row_emitted_ = true;
        }

        constexpr double tolerance = 1e-12;
        if (dynamics_->reviews_blocks()) {
            save_block_state();
            block_start_time_ = time_;
        }
        while (time_ < parameters_.tmax - tolerance) {
            const double next_time = std::min( {next_clone_time_, next_output_time_, parameters_.tmax});
            if (!advance_dynamics_until(next_time)) {
                // The weights spread past the tolerated ratio before the next
                // scheduled instant.  Examine them now and resume; no row is
                // due.  The exact range may show that resampling is not yet
                // necessary because the certificate is an upper bound.
                clone_population(false);
                schedule_next_cloning_check();
                continue;
            }

            // A mechanism that predicts its own interval only finds out
            // afterwards whether the prediction held.  Judge the block here,
            // before anything is recorded or resampled: a block that let the
            // weights spread past target_max is discarded exactly and replayed
            // with a shorter step.
            if (dynamics_->reviews_blocks()) {
                const double range = population_log_weight_range();
                if (dynamics_->review_block(range, time_ - block_start_time_)) {
                    restore_block_state();
                    ++rejected_blocks_;
                    schedule_next_cloning_check();
                    save_block_state();
                    block_start_time_ = time_;
                    continue;
                }
            }

            const bool clone_due = std::abs(time_ - next_clone_time_) <= tolerance;
            const bool scheduled_output_due = std::abs(time_ - next_output_time_) <= tolerance;
            const bool final_time = std::abs(time_ - parameters_.tmax) <= tolerance;
            const bool output_due = scheduled_output_due || final_time;

            if (output_due) {
                // With a tilt, absorb the weights before taking unweighted
                // population averages. At s=k=0 all weights stay exactly one,
                // so resampling would only copy the population and redraw its
                // clocks without changing the ensemble.
                if (!unbiased_) {
                    clone_population(true);
                    schedule_next_cloning_check();
                }
                rows.push_back(output_record());
                ++recordings_completed_;
                reset_recording_window();
                if (scheduled_output_due) {
                    ++output_number_;
                    next_output_time_ = static_cast<double>(output_number_) * record_interval_;
                }
                if (parameters_.progress > 0 && recordings_completed_ % parameters_.progress == 0) {
                    std::cerr << "Progress: " << std::fixed << std::setprecision(1) << 100.0 * time_ / parameters_.tmax << "%\n";
                }
            } else if (clone_due) {
                // An explicit cloning interval is a manual request to inspect
                // the weights.  It is not enabled automatically.
                clone_population(false);
                schedule_next_cloning_check();
            }
            if (dynamics_->reviews_blocks()) {
                save_block_state();
                block_start_time_ = time_;
            }
        }
        return rows;
    }

    std::string checkpoint() const {
        std::string bytes;
        append(bytes, std::uint64_t{0x4d45524543485031ULL}); // MERECHP1
        append(bytes, std::uint32_t{9});
        append(bytes, static_cast<std::uint64_t>(parameters_.L));
        append(bytes, static_cast<std::uint64_t>(parameters_.M));
        append(bytes, static_cast<std::uint32_t>(parameters_.dynamics));
        append(bytes, model_.periodic);
        append(bytes, time_);
        append(bytes, next_clone_time_);
        append(bytes, next_output_time_);
        append(bytes, static_cast<std::uint64_t>(output_number_));
        append(bytes, static_cast<std::uint64_t>(recordings_completed_));
        append(bytes, initial_row_emitted_);
        append(bytes, reference_current_);
        append(bytes, log_normalization_);
        append(bytes, effective_sample_size_sum_);
        append(bytes, effective_sample_size_samples_);
        append_vector(bytes, reference_);
        append_vector(bytes, configs_);
        append_vector(bytes, event_times_);
        append(bytes, accepted_walker_events_);
        append(bytes, rejected_walker_events_);
        append(bytes, accepted_reference_events_);
        append(bytes, rejected_reference_events_);
        append(bytes, resampling_steps_);
        append(bytes, resampling_checks_);
        append(bytes, walkers_killed_);
        append(bytes, extra_clones_);
        for (const Walker& walker : walkers_) {
            append(bytes, walker.log_weight);
            append(bytes, walker.last_weight_time);
            append(bytes, walker.current);
        }
        append(bytes, dynamics_->schedule_state());
        std::ostringstream rng_stream;
        rng_stream << rng_;
        append_string(bytes, rng_stream.str());
        return bytes;
    }

    void restore(const std::string& bytes) {
        std::size_t offset = 0;
        if (read<std::uint64_t>(bytes, offset) != 0x4d45524543485031ULL || read<std::uint32_t>(bytes, offset) != 9) {
            throw std::invalid_argument("Unsupported C++ checkpoint format");
        }
        const auto saved_L = read<std::uint64_t>(bytes, offset);
        const auto saved_M = read<std::uint64_t>(bytes, offset);
        const auto saved_dynamics = static_cast<DynamicsKind>( read<std::uint32_t>(bytes, offset));
        const bool saved_periodic = read<bool>(bytes, offset);
        if (saved_L != parameters_.L || saved_M != parameters_.M) {
            throw std::invalid_argument("Checkpoint L or M does not match");
        }
        if (saved_dynamics != parameters_.dynamics)
            throw std::invalid_argument("Checkpoint dynamics mode does not match");
        if (saved_periodic != model_.periodic)
            throw std::invalid_argument("Checkpoint boundary condition does not match");
        time_ = read<double>(bytes, offset);
        next_clone_time_ = read<double>(bytes, offset);
        next_output_time_ = read<double>(bytes, offset);
        output_number_ = read<std::uint64_t>(bytes, offset);
        recordings_completed_ = read<std::uint64_t>(bytes, offset);
        initial_row_emitted_ = read<bool>(bytes, offset);
        reference_current_ = read<double>(bytes, offset);
        log_normalization_ = read<double>(bytes, offset);
        effective_sample_size_sum_ = read<double>(bytes, offset);
        effective_sample_size_samples_ = read<std::uint64_t>(bytes, offset);
        reference_ = read_vector<std::uint64_t>(bytes, offset);
        configs_ = read_vector<std::uint64_t>(bytes, offset);
        event_times_ = read_vector<double>(bytes, offset);
        accepted_walker_events_ = read<std::uint64_t>(bytes, offset);
        rejected_walker_events_ = read<std::uint64_t>(bytes, offset);
        accepted_reference_events_ = read<std::uint64_t>(bytes, offset);
        rejected_reference_events_ = read<std::uint64_t>(bytes, offset);
        resampling_steps_ = read<std::uint64_t>(bytes, offset);
        resampling_checks_ = read<std::uint64_t>(bytes, offset);
        walkers_killed_ = read<std::uint64_t>(bytes, offset);
        extra_clones_ = read<std::uint64_t>(bytes, offset);
        for (Walker& walker : walkers_) {
            walker.log_weight = read<double>(bytes, offset);
            walker.last_weight_time = read<double>(bytes, offset);
            walker.current = read<double>(bytes, offset);
        }
        dynamics_->set_schedule_state(read<double>(bytes, offset));
        const std::string rng_text = read_string(bytes, offset);
        std::istringstream rng_stream(rng_text);
        rng_stream >> rng_;
        if (!rng_stream || offset != bytes.size()) {
            throw std::invalid_argument("Corrupt C++ checkpoint");
        }
        refresh_reference();
        for (std::size_t i = 0; i < parameters_.M; ++i) refresh_walker(i, false);
        if (!unbiased_) {
            // Restored weights need not be uniform, so restart the certificate
            // from the range the checkpoint actually carries.
            double lowest = walkers_[0].log_weight;
            double highest = walkers_[0].log_weight;
            for (std::size_t i = 1; i < parameters_.M; ++i) {
                lowest = std::min(lowest, walkers_[i].log_weight);
                highest = std::max(highest, walkers_[i].log_weight);
            }
            reanchor_weight_bound(highest - lowest);
        }
        // The default recording grid is a fraction of tmax, and resuming a
        // checkpoint is exactly the case where tmax grows.  Re-derive the next
        // recording time from the restored clock so the grid stays consistent
        // with the interval this instance will actually use.
        if (parameters_.record_interval <= 0) {
            output_number_ = static_cast<std::uint64_t>( std::floor(time_ / record_interval_)) + 1;
            next_output_time_ = static_cast<double>(output_number_) * record_interval_;
        }
        dynamics_->restore_schedule();
    }

    double time() const { return time_; }
    // Derived on demand rather than stored: the value is whatever the
    // mechanism and any explicit interval currently imply, and zero is the
    // reported spelling of "no interval of its own".
    // The reference trajectory exists only to be measured against: at k = 0
    // nothing reads its configuration, so it is not evolved at all and its
    // reported current stays at zero.
    bool monitored() const { return model_.k != 0.0; }

    double record_interval() const { return record_interval_; }
    double estimated_jump_rate() const { return estimated_jump_rate_; }
    bool unbiased() const { return unbiased_; }
    std::uint64_t accepted_walker_events() const {
        return accepted_walker_events_;
    }
    std::uint64_t accepted_reference_events() const {
        return accepted_reference_events_;
    }
    std::uint64_t rejected_walker_events() const {
        return rejected_walker_events_;
    }
    std::uint64_t rejected_reference_events() const {
        return rejected_reference_events_;
    }
    DynamicsKind dynamics_kind() const { return parameters_.dynamics; }
    const char* dynamics_name() const {
        switch (parameters_.dynamics) {
            case DynamicsKind::Uniformized: return "uniformized";
            case DynamicsKind::Direct: return "direct";
            default: return "heap";
        }
    }
    bool periodic() const { return model_.periodic; }
    std::uint64_t resampling_steps() const { return resampling_steps_; }
    std::uint64_t resampling_checks() const { return resampling_checks_; }
    std::uint64_t rejected_blocks() const { return rejected_blocks_; }
    double target_min() const { return parameters_.target_min; }
    double target_max() const { return parameters_.target_max; }
    double log_target_min() const { return log_resample_range_; }
    double log_target_max() const { return log_target_max_; }
    bool custom_weights() const { return model_.custom_weights; }
    const std::string& rate_model() const { return model_.rate_model; }
    std::vector<double> right_weights() const {
        return {model_.right_weight, model_.right_weight + classes};
    }
    std::vector<double> left_weights() const {
        return {model_.left_weight, model_.left_weight + classes};
    }
    std::uint64_t walkers_killed() const { return walkers_killed_; }

private:
    static std::uint64_t final_mask(std::size_t bits) {
        if (bits == 0) return 0;
        const std::size_t remainder = bits % 64;
        return remainder == 0 ? ~std::uint64_t{0} : (std::uint64_t{1} << remainder) - 1;
    }

    void validate_parameters() const {
        if (parameters_.L < 2) throw std::invalid_argument("L must be at least 2");
        if (parameters_.M < 2) throw std::invalid_argument("M must be at least 2");
        if (parameters_.tmax <= 0) throw std::invalid_argument("tmax must be positive");
        if (std::min({model_.k, model_.alpha, model_.gamma, model_.delta, model_.beta}) < 0) {
            throw std::invalid_argument("k and boundary rates must be non-negative");
        }
        const std::size_t expected = static_cast<std::size_t>(classes);
        const std::size_t given_right = parameters_.right_weights.size();
        const std::size_t given_left = parameters_.left_weights.size();
        if ((given_right != 0 && given_right != expected) || (given_left != 0 && given_left != expected)) {
            throw std::invalid_argument( "right_weights and left_weights must each hold four rates");
        }
        if ((given_right == 0) != (given_left == 0)) {
            throw std::invalid_argument( "right_weights and left_weights must be given together");
        }
        for (int cls = 0; cls < classes; ++cls) {
            if (!(model_.right_weight[cls] > 0.0) || !(model_.left_weight[cls] > 0.0)) {
                throw std::invalid_argument("every hop rate must be positive");
            }
        }
        if (!(parameters_.target_min > 1))
            throw std::invalid_argument("target_min must exceed 1");
        if (!(parameters_.target_max > parameters_.target_min))
            throw std::invalid_argument("target_max must exceed target_min");
        if (!std::isfinite(parameters_.filling) || parameters_.filling < 0.0 || parameters_.filling > 1.0) {
            throw std::invalid_argument("filling must lie between 0 and 1");
        }
        if (parameters_.boundary != "open" && parameters_.boundary != "periodic") {
            throw std::invalid_argument("boundary must be open or periodic");
        }
        if (parameters_.cloning_interval <= 0 || parameters_.record_interval == 0) {
            throw std::invalid_argument("Intervals must be positive when supplied");
        }
    }

    void choose_intervals() {
        // The gap between the two targets is what keeps examinations from
        // crowding together: one that declines to resample re-anchors at the
        // exact range, which sits below target_min and therefore strictly
        // below the trigger at target_max.  Their historical relation was
        // target_min = target_max^(2/3), which remains the default pair.
        log_target_max_ = std::log(parameters_.target_max);
        log_resample_range_ = std::log(parameters_.target_min);
        estimated_jump_rate_ = model_.broad_jump_rate_estimate();
        // Recording is O(M*L) per row and every row is kept in memory until the
        // run returns, so the default is a fixed hundred rows over the whole
        // evolution rather than one row per cloning time.
        record_interval_ = parameters_.record_interval > 0 ? parameters_.record_interval : 0.01 * parameters_.tmax;
        schedule_next_cloning_check();
        next_output_time_ = record_interval_;
    }

    // Heap and uniformized dynamics leave this to the running certificate and
    // only honour an explicitly requested clock.  Direct dynamics have no
    // certificate, so they fix the next interval here from the spread of the
    // potentials: the log-weight ratio between two walkers grows as the
    // integral of their potential difference, so over an interval tau it
    // cannot exceed tau*(V_max - V_min), and log(target_max) divided by
    // that spread is the time at which the tolerated ratio is reached.
    void schedule_next_cloning_check() {
        const double interval = unbiased_ ? std::numeric_limits<double>::infinity() : std::min(dynamics_->step_size(), parameters_.cloning_interval);
        next_clone_time_ = time_ + interval;
    }

    std::uint64_t& config_word(std::size_t walker, std::size_t word) {
        return configs_[walker * words_ + word];
    }
    std::uint64_t config_word(std::size_t walker, std::size_t word) const {
        return configs_[walker * words_ + word];
    }
    bool occupied(const std::vector<std::uint64_t>& config, std::size_t site) const {
        return (config[site / 64] >> (site % 64)) & 1U;
    }
    bool occupied(const std::uint64_t* config, std::size_t site) const {
        return (config[site / 64] >> (site % 64)) & 1U;
    }
    bool occupied(std::size_t walker, std::size_t site) const {
        return (config_word(walker, site / 64) >> (site % 64)) & 1U;
    }
    void flip(std::vector<std::uint64_t>& config, std::size_t site) {
        config[site / 64] ^= std::uint64_t{1} << (site % 64);
    }
    void flip(std::size_t walker, std::size_t site) {
        config_word(walker, site / 64) ^= std::uint64_t{1} << (site % 64);
    }

    void initialize_configuration() {
        std::vector<std::uint64_t> initial(words_, 0);
        if (model_.periodic) {
            // A closed ring conserves particle number. Sample uniformly from
            // configurations with the closest integer to filling * L.
            const std::size_t particles = static_cast<std::size_t>(std::llround( parameters_.filling * static_cast<double>(parameters_.L)));
            std::vector<std::size_t> sites(parameters_.L);
            for (std::size_t site = 0; site < parameters_.L; ++site)
                sites[site] = site;
            std::shuffle(sites.begin(), sites.end(), rng_);
            for (std::size_t index = 0; index < particles; ++index) {
                const std::size_t site = sites[index];
                initial[site / 64] |= std::uint64_t{1} << (site % 64);
            }
        } else if (parameters_.initial == "alternating") {
            for (std::size_t site = 0; site < parameters_.L; site += 2)
                initial[site / 64] |= std::uint64_t{1} << (site % 64);
        } else if (parameters_.initial == "random") {
            for (std::uint64_t& word : initial) word = rng_();
            initial.back() &= final_word_mask_;
        } else if (parameters_.initial == "full") {
            for (std::uint64_t& word : initial) word = ~std::uint64_t{0};
            initial.back() &= final_word_mask_;
        } else if (parameters_.initial != "empty") {
            throw std::invalid_argument("Unknown initial condition");
        }
        reference_ = initial;
        for (std::size_t i = 0; i < parameters_.M; ++i)
            std::copy(initial.begin(), initial.end(), configs_.begin() + static_cast<std::ptrdiff_t>(i * words_));
    }

    void initialize_state() {
        refresh_reference();
        for (std::size_t i = 0; i < parameters_.M; ++i) {
            walkers_[i].last_weight_time = 0.0;
            refresh_walker(i, false);
        }
        if (!unbiased_) reanchor_weight_bound(0.0);
    }

    // Masks are stored per walker as 2*classes contiguous bit vectors of
    // words_ words: the right-allowed bonds of each flanking class first, then
    // the left-allowed ones.
    static constexpr std::size_t mask_vectors = 2 * classes;
    std::size_t mask_offset(std::size_t walker, int direction, int cls) const {
        return ((walker * mask_vectors) + static_cast<std::size_t>(direction * classes + cls)) * words_;
    }

    // Flanking class of a bond, 2*n_{b-1} + n_{b+2}.  On a ring both flanks
    // always exist.  With open boundaries the outermost bonds have a flank off
    // the lattice; that missing occupation is mirrored from the one that is
    // present, so the pair is equal and the rate reduces to unity whenever
    // delta vanishes.  This keeps the two ends particle-hole symmetric, which
    // matters for the half-filled KLS analysis, and is in any case an
    // order-1/L surface convention.
    int flanking_class(const std::uint64_t* config, std::size_t bond) const {
        const std::size_t L = parameters_.L;
        int behind, ahead;
        const bool has_behind = model_.periodic || bond >= 1;
        const bool has_ahead = model_.periodic || bond + 2 < L;
        if (has_behind)
            behind = occupied(config, (bond + L - 1) % L) ? 1 : 0;
        if (has_ahead)
            ahead = occupied(config, (bond + 2) % L) ? 1 : 0;
        if (!has_behind && !has_ahead) behind = ahead = 0;
        else if (!has_behind) behind = ahead;
        else if (!has_ahead) ahead = behind;
        return 2 * behind + ahead;
    }

    // Remove a bond from whichever class vector currently holds it.
    void clear_bond(std::uint64_t* masks, std::size_t* nR, std::size_t* nL, std::size_t bond) const {
        const std::size_t word = bond / 64;
        const std::uint64_t bit = std::uint64_t{1} << (bond % 64);
        for (int cls = 0; cls < classes; ++cls) {
            std::uint64_t& r = masks[static_cast<std::size_t>(cls) * words_ + word];
            if (r & bit) { r &= ~bit; --nR[cls]; }
            std::uint64_t& l = masks[static_cast<std::size_t>(classes + cls) * words_ + word];
            if (l & bit) { l &= ~bit; --nL[cls]; }
        }
    }

    // Classify a bond and file it under the right vector.
    void set_bond(const std::uint64_t* config, std::uint64_t* masks, std::size_t* nR, std::size_t* nL, std::size_t bond) const {
        const std::size_t L = parameters_.L;
        const bool occupied_left = occupied(config, bond);
        const bool occupied_right = occupied(config, (bond + 1) % L);
        if (occupied_left == occupied_right) return;
        const int cls = flanking_class(config, bond);
        const std::size_t word = bond / 64;
        const std::uint64_t bit = std::uint64_t{1} << (bond % 64);
        if (occupied_left) {
            masks[static_cast<std::size_t>(cls) * words_ + word] |= bit;
            ++nR[cls];
        } else {
            masks[static_cast<std::size_t>(classes + cls) * words_ + word] |= bit;
            ++nL[cls];
        }
    }

    // WASEP bookkeeping.  When every class carries the same rate the flanking
    // occupations no longer matter, so a bond is allowed purely by its own two
    // sites and all allowed bonds can be found a word at a time: shifting the
    // configuration down by one site puts n_{b+1} beside n_b, and the two
    // masks are then ``n_b & ~n_{b+1}`` and ``~n_b & n_{b+1}``.  Everything is
    // filed under class 0, which the selection code already handles because the
    // remaining class vectors stay empty and contribute nothing.
    void rebuild_masks_uniform(const std::uint64_t* config, std::uint64_t* masks, std::size_t* nR, std::size_t* nL) const {
        std::fill_n(masks, mask_vectors * words_, std::uint64_t{0});
        std::fill_n(nR, classes, std::size_t{0});
        std::fill_n(nL, classes, std::size_t{0});
        std::uint64_t* right = masks;
        std::uint64_t* left = masks + static_cast<std::size_t>(classes) * words_;
        std::size_t right_count = 0;
        std::size_t left_count = 0;
        for (std::size_t word = 0; word < words_; ++word) {
            const std::uint64_t current = config[word];
            const std::uint64_t next_low = word + 1 < words_ ? (config[word + 1] & 1U) : 0;
            const std::uint64_t shifted = (current >> 1) | (next_low << 63);
            std::uint64_t r = current & ~shifted;
            std::uint64_t l = ~current & shifted;
            if (word + 1 == words_) {
                r &= final_bond_mask_;
                l &= final_bond_mask_;
            }
            right[word] = r;
            left[word] = l;
            right_count += static_cast<std::size_t>(std::popcount(r));
            left_count += static_cast<std::size_t>(std::popcount(l));
        }
        if (model_.periodic) {
            // The shift above covers the bonds x -> x+1 with x < L-1; the
            // closing bond L-1 -> 0 wraps and has to be filed by hand.
            const std::size_t bond = parameters_.L - 1;
            const std::size_t word = bond / 64;
            const std::uint64_t bit = std::uint64_t{1} << (bond % 64);
            const bool last = (config[word] & bit) != 0;
            const bool first = (config[0] & 1U) != 0;
            right[word] &= ~bit;
            left[word] &= ~bit;
            if (last && !first) {
                right[word] |= bit;
                ++right_count;
            } else if (!last && first) {
                left[word] |= bit;
                ++left_count;
            }
        }
        nR[0] = right_count;
        nL[0] = left_count;
    }

    void build_masks(const std::uint64_t* config, std::uint64_t* masks, std::size_t* nR, std::size_t* nL) const {
        if (model_.uniform) {
            rebuild_masks_uniform(config, masks, nR, nL);
            return;
        }
        std::fill_n(masks, mask_vectors * words_, std::uint64_t{0});
        std::fill_n(nR, classes, std::size_t{0});
        std::fill_n(nL, classes, std::size_t{0});
        for (std::size_t bond = 0; bond < model_.links; ++bond)
            set_bond(config, masks, nR, nL, bond);
    }

    // A hop across bond b flips sites b and b+1.  Whether a bond permits a hop
    // depends on its own two sites, so only bonds b-1, b, b+1 can change in
    // that respect.  Its *rate class*, however, depends on the four-site window
    // (b'-1, b', b'+1, b'+2), and that window still contains site b or b+1 for
    // b' = b-2 and b' = b+2.  Those two bonds keep their occupation pattern but
    // migrate between class vectors, so the update has to span b-2 .. b+2.
    void update_masks_near(const std::uint64_t* config, std::uint64_t* masks, std::size_t* nR, std::size_t* nL, std::size_t bond) const {
        if (model_.uniform) {
            // Rebuilding costs one pass over ``words_`` words, which for any
            // lattice that fits in a handful of words is cheaper than
            // reclassifying five bonds across eight class vectors.
            rebuild_masks_uniform(config, masks, nR, nL);
            return;
        }
        const std::size_t links = model_.links;
        for (int offset = -2; offset <= 2; ++offset) {
            const std::ptrdiff_t raw = static_cast<std::ptrdiff_t>(bond) + offset;
            std::size_t target;
            if (model_.periodic) {
                target = static_cast<std::size_t>( (raw % static_cast<std::ptrdiff_t>(links) + static_cast<std::ptrdiff_t>(links)) % static_cast<std::ptrdiff_t>(links));
            } else {
                if (raw < 0 || raw >= static_cast<std::ptrdiff_t>(links))
                    continue;
                target = static_cast<std::size_t>(raw);
            }
            clear_bond(masks, nR, nL, target);
            set_bond(config, masks, nR, nL, target);
        }
    }

    double weighted_activity(const std::size_t* counts, const double* weights) const {
        double total = 0.0;
        for (int cls = 0; cls < classes; ++cls)
            total += static_cast<double>(counts[cls]) * weights[cls];
        return total;
    }

    // The reference configuration a walker's overlap is measured against.
    // Direct dynamics evolve the reference ahead of the walkers and point this
    // at a shadow copy that is replayed forward walker by walker, so that a
    // walker never compares itself with the reference's future.
    const std::uint64_t* overlap_reference() const {
        return overlap_reference_ != nullptr ? overlap_reference_ : reference_.data();
    }

    void set_overlap_reference(const std::uint64_t* configuration) {
        overlap_reference_ = configuration;
    }

    int overlap(std::size_t walker) const {
        const std::uint64_t* against = overlap_reference();
        std::size_t mismatches = 0;
        for (std::size_t word = 0; word < words_; ++word)
            mismatches += std::popcount(config_word(walker, word) ^ against[word]);
        return static_cast<int>(parameters_.L - 2 * mismatches);
    }

    void refresh_walker(std::size_t i, bool keep_overlap) {
        Walker& walker = walkers_[i];
        build_masks(configs_.data() + i * words_, masks_.data() + mask_offset(i, 0, 0), walker.nR, walker.nL);
        walker.activity_right = weighted_activity(walker.nR, model_.right_weight);
        walker.activity_left = weighted_activity(walker.nL, model_.left_weight);
        if (!keep_overlap)
            walker.overlap = model_.k == 0.0 ? 0 : overlap(i);
        if (!unbiased_)
            walker.potential = model_.counting_potential( walker.activity_right, walker.activity_left) + model_.k_scaled * static_cast<double>(walker.overlap);
    }

    // Same as refresh_walker, but only the five bonds whose rate class the hop
    // across ``bond`` can have disturbed are reclassified.
    void refresh_walker_near(std::size_t i, std::size_t bond) {
        Walker& walker = walkers_[i];
        update_masks_near(configs_.data() + i * words_, masks_.data() + mask_offset(i, 0, 0), walker.nR, walker.nL, bond);
        walker.activity_right = weighted_activity(walker.nR, model_.right_weight);
        walker.activity_left = weighted_activity(walker.nL, model_.left_weight);
        if (!unbiased_)
            walker.potential = model_.counting_potential( walker.activity_right, walker.activity_left) + model_.k_scaled * static_cast<double>(walker.overlap);
    }

    void refresh_reference_near(std::size_t bond) {
        update_masks_near(reference_.data(), reference_masks_.data(), reference_nR_, reference_nL_, bond);
        reference_activity_right_ = weighted_activity(reference_nR_, model_.right_weight);
        reference_activity_left_ = weighted_activity(reference_nL_, model_.left_weight);
    }

    void refresh_reference() {
        build_masks(reference_.data(), reference_masks_.data(), reference_nR_, reference_nL_);
        reference_activity_right_ = weighted_activity(reference_nR_, model_.right_weight);
        reference_activity_left_ = weighted_activity(reference_nL_, model_.left_weight);
    }

    // Advance the certificate to the current time using the potential spread
    // that was in force over the elapsed stretch, then widen that spread with a
    // freshly changed potential.  Both steps are O(1) and must happen on every
    // event: a potential that grows inside an interval raises the rate at which
    // the weights can separate, and a bound fixed in advance would miss it.
    void accumulate_weight_bound() {
        weight_range_bound_ += (potential_maximum_ - potential_minimum_)
                             * (time_ - weight_bound_time_);
        weight_bound_time_ = time_;
    }

    void widen_potential_extremes(double potential) {
        if (potential > potential_maximum_) potential_maximum_ = potential;
        if (potential < potential_minimum_) potential_minimum_ = potential;
    }

    bool weight_bound_reached() const {
        return weight_range_bound_ >= log_target_max_;
    }

    // Called only where the whole population is already being traversed.
    void reanchor_weight_bound(double exact_range) {
        weight_range_bound_ = exact_range;
        weight_bound_time_ = time_;
        potential_maximum_ = walkers_[0].potential;
        potential_minimum_ = walkers_[0].potential;
        for (std::size_t i = 1; i < parameters_.M; ++i)
            widen_potential_extremes(walkers_[i].potential);
    }

    double walker_rate(std::size_t i) const {
        const Walker& w = walkers_[i];
        double rate = w.activity_right * model_.ps + w.activity_left * model_.qs;
        if (!model_.periodic) {
            rate += occupied(i, 0) ? model_.gamma : model_.alpha;
            rate += occupied(i, parameters_.L - 1) ? model_.beta : model_.delta;
        }
        return rate;
    }

    double reference_rate() const {
        double rate = reference_activity_right_ * model_.p0 + reference_activity_left_ * model_.q0;
        if (!model_.periodic) {
            rate += occupied(reference_, 0) ? model_.gamma : model_.alpha;
            rate += occupied(reference_, parameters_.L - 1) ? model_.beta : model_.delta;
        }
        return rate;
    }

    double uniform_open() {
        // Bit-for-bit what std::generate_canonical<double,53> returns for a
        // 64-bit engine, without its per-call template machinery: one draw
        // already carries more bits than a double can hold, so the whole word
        // is rounded to nearest and scaled by an exact power of two.  The
        // clamp reproduces the guard against a draw that rounds up to 1.0.
        const double value = static_cast<double>(rng_()) * 0x1p-64;
        return value < 1.0 ? value : std::nextafter(1.0, 0.0);
    }

    double exponential_wait(double rate) {
        if (rate <= 0) return std::numeric_limits<double>::infinity();
        double u = uniform_open();
        while (u <= 0.0) u = uniform_open();
        return -std::log(u) / rate;
    }

    std::size_t select_set_bit(const std::uint64_t* masks, std::size_t ordinal) const {
        for (std::size_t word = 0; word < words_; ++word) {
            std::uint64_t mask = masks[word];
            const std::size_t count = std::popcount(mask);
            if (ordinal >= count) {
                ordinal -= count;
                continue;
            }
            return word * 64 + select_in_word(mask, ordinal);
        }
        throw std::logic_error("Allowed-transition mask selection failed");
    }

    void apply_walker_event(std::size_t i, std::size_t changed_a, std::size_t changed_b, double current_change) {
        Walker& walker = walkers_[i];
        if (!unbiased_) {
            walker.update_weight(time_);
            // Direct dynamics move each walker on its own clock, so a running
            // bound over the population has no meaning there and is not kept.
            if (track_certificate_) accumulate_weight_bound();
        }
        // For short configurations a handful of hardware POPCNT operations
        // beats the indexed tests below. Beyond eight words, update overlap
        // from only the flipped sites so the monitored path remains O(1).
        constexpr std::size_t overlap_scan_words = 8;
        int overlap_change = 0;
        if (model_.k != 0.0 && words_ > overlap_scan_words) {
            const std::uint64_t* against = overlap_reference();
            int old_mismatches = occupied(i, changed_a) != occupied(against, changed_a);
            if (changed_b < parameters_.L)
                old_mismatches += occupied(i, changed_b) != occupied(against, changed_b);
            const int changed_count = changed_b < parameters_.L ? 2 : 1;
            overlap_change = 4 * old_mismatches - 2 * changed_count;
        }
        flip(i, changed_a);
        if (changed_b < parameters_.L) flip(i, changed_b);
        walker.current += current_change;
        if (model_.k != 0.0) {
            walker.overlap = words_ <= overlap_scan_words ? overlap(i) : walker.overlap + overlap_change;
        }
        if (changed_b < parameters_.L) {
            refresh_walker_near(i, changed_a);
        } else {
            // A reservoir event flips a single boundary site, which disturbs the
            // bonds within two of it in the same way a hop does.
            refresh_walker_near(i, changed_a == 0 ? 0 : model_.links - 1);
        }
        if (!unbiased_ && track_certificate_)
            widen_potential_extremes(walker.potential);
        ++accepted_walker_events_;
    }

    // Walk the class vectors of one direction, consuming the uniform draw.
    // Every bond inside a class shares one rate, so once the draw lands in a
    // class the remainder divided by that rate is a uniform index into it --
    // the same trick the two-rate WASEP uses, applied to eight categories.
    bool select_within_direction(const std::uint64_t* masks, const std::size_t* counts, const double* weights, double base_rate, double& choice, std::size_t& bond) const {
        for (int cls = 0; cls < classes; ++cls) {
            if (counts[cls] == 0) continue;
            const double rate = weights[cls] * base_rate;
            const double block = static_cast<double>(counts[cls]) * rate;
            if (choice < block) {
                const std::size_t ordinal = std::min( static_cast<std::size_t>(choice / rate), counts[cls] - 1);
                bond = select_set_bit( masks + static_cast<std::size_t>(cls) * words_, ordinal);
                return true;
            }
            choice -= block;
        }
        return false;
    }

    void execute_walker_event(std::size_t i) {
        const Walker& walker = walkers_[i];
        const double right_total = walker.activity_right * model_.ps;
        const double left_total = walker.activity_left * model_.qs;
        const double left_boundary = model_.periodic ? 0.0 : (occupied(i, 0) ? model_.gamma : model_.alpha);
        const double right_boundary = model_.periodic ? 0.0 : (occupied(i, parameters_.L - 1) ? model_.beta : model_.delta);
        const double total = right_total + left_total + left_boundary + right_boundary;
        double choice = uniform_open() * total;
        std::size_t changed_a = 0;
        std::size_t changed_b = parameters_.L;
        double current_change = 0.0;
        std::uint64_t* walker_masks = masks_.data() + mask_offset(i, 0, 0);
        std::size_t bond = 0;
        if (choice < right_total && select_within_direction(walker_masks, walker.nR, model_.right_weight, model_.ps, choice, bond)) {
            changed_a = bond;
            changed_b = (bond + 1) % parameters_.L;
            current_change = model_.dq;
        } else if ((choice = std::max(choice - right_total, 0.0)) < left_total && select_within_direction( walker_masks + static_cast<std::size_t>(classes) * words_, walker.nL, model_.left_weight, model_.qs, choice, bond)) { 
            changed_a = bond; changed_b = (bond + 1) % parameters_.L; current_change = -model_.dq; 
        } else if ((choice = std::max(choice - left_total, 0.0)) < left_boundary) { 
            changed_a = 0; 
        } else { 
            changed_a = parameters_.L - 1; 
        } 
        apply_walker_event(i, changed_a, changed_b, current_change); 
    }

    void apply_reference_event(std::size_t changed_a, std::size_t changed_b, double current_change, bool update_walkers = true) {
        if (model_.k != 0.0 && update_walkers) {
            // Every walker's overlap with the reference shifts at once, so the
            // certificate has to close off the preceding stretch before any
            // potential moves.
            accumulate_weight_bound();
            for (std::size_t i = 0; i < parameters_.M; ++i) {
                Walker& walker = walkers_[i];
                walker.update_weight(time_);
                int old_mismatches = occupied(i, changed_a) != occupied(reference_, changed_a);
                if (changed_b < parameters_.L)
                    old_mismatches += occupied(i, changed_b) != occupied(reference_, changed_b);
                const int changed_count = changed_b < parameters_.L ? 2 : 1;
                const int overlap_change = 4 * old_mismatches - 2 * changed_count;
                walker.overlap += overlap_change;
                walker.potential += model_.k_scaled * overlap_change;
                widen_potential_extremes(walker.potential);
            }
        }
        flip(reference_, changed_a);
        if (changed_b < parameters_.L) flip(reference_, changed_b);
        reference_current_ += current_change;
        refresh_reference_near(changed_b < parameters_.L ? changed_a: (changed_a == 0 ? 0 : model_.links - 1));
        ++accepted_reference_events_;
    }

    void execute_reference_event( std::vector<ReferenceEventRecord>* log = nullptr) {
        const double right_total = reference_activity_right_ * model_.p0;
        const double left_total = reference_activity_left_ * model_.q0;
        const double left_boundary = model_.periodic ? 0.0 : (occupied(reference_, 0) ? model_.gamma : model_.alpha);
        const double right_boundary = model_.periodic ? 0.0 : (occupied(reference_, parameters_.L - 1) ? model_.beta : model_.delta);
        const double total = right_total + left_total + left_boundary + right_boundary;
        double choice = uniform_open() * total;
        std::size_t changed_a = 0;
        std::size_t changed_b = parameters_.L;
        double current_change = 0.0;
        if (choice < right_total) {
            std::size_t bond = 0;
            select_within_direction(reference_masks_.data(), reference_nR_, model_.right_weight, model_.p0, choice, bond);
            changed_a = bond;
            changed_b = (changed_a + 1) % parameters_.L;
            current_change = model_.dq;
        } else if ((choice = std::max(choice - right_total, 0.0)) < left_total) {
            std::size_t bond = 0;
            select_within_direction( reference_masks_.data() + static_cast<std::size_t>(classes) * words_, reference_nL_, model_.left_weight, model_.q0, choice, bond);
            changed_a = bond;
            changed_b = (changed_a + 1) % parameters_.L;
            current_change = -model_.dq;
        } else if ((choice = std::max(choice - left_total, 0.0)) < left_boundary) {
            changed_a = 0;
        } else {
            changed_a = parameters_.L - 1;
        }

        if (log != nullptr) {
            // Store the reference's occupations before the flip; that is all a
            // walker needs to update its overlap once the reference has run
            // ahead of it.  The walkers are left untouched here and replay this
            // record themselves.
            log->push_back({time_, changed_a, changed_b, occupied(reference_, changed_a), changed_b < parameters_.L ? occupied(reference_, changed_b) : false});
            apply_reference_event(changed_a, changed_b, current_change, false);
            return;
        }
        apply_reference_event(changed_a, changed_b, current_change);
    }

    // then rejects down to the true one, so a proposal is accepted with
    // probability weight(class)/heaviest_weight on top of the exclusion test.
    bool accept_class_rate(const std::uint64_t* config, std::size_t bond, const double* weights) {
        const double weight = weights[flanking_class(config, bond)];
        if (weight >= model_.heaviest_weight) return true;
        return uniform_open() * model_.heaviest_weight < weight;
    }

    void attempt_uniformized_walker(std::size_t i, double choice) {
        const double scale = model_.heaviest_weight;
        const double right_end = static_cast<double>(model_.links) * model_.ps * scale;
        const double left_end = right_end + static_cast<double>(model_.links) * model_.qs * scale;
        bool accepted = false;
        std::size_t changed_a = 0;
        std::size_t changed_b = parameters_.L;
        double current_change = 0.0;
        if (choice < right_end) {
            const std::size_t bond = std::min( static_cast<std::size_t>(choice / (model_.ps * scale)), model_.links - 1);
            const std::size_t neighbor = (bond + 1) % parameters_.L;
            accepted = occupied(i, bond) && !occupied(i, neighbor) && accept_class_rate(configs_.data() + i * words_, bond, model_.right_weight);
            changed_a = bond;
            changed_b = neighbor;
            current_change = model_.dq;
        } else if (choice < left_end) {
            const std::size_t bond = std::min( static_cast<std::size_t>((choice - right_end) / (model_.qs * scale)), model_.links - 1);
            const std::size_t neighbor = (bond + 1) % parameters_.L;
            accepted = !occupied(i, bond) && occupied(i, neighbor) && accept_class_rate(configs_.data() + i * words_, bond, model_.left_weight);
            changed_a = bond;
            changed_b = neighbor;
            current_change = -model_.dq;
        } else if (!model_.periodic) {
            // A ring has no reservoirs and its proposal rate excludes them, so
            // only a floating-point tie can reach this far.  Without the guard
            // such a tie would inject or remove a particle and break the
            // conservation law of the closed chain.
            choice -= left_end;
            if (choice < model_.alpha) {
                accepted = !occupied(i, 0);
                changed_a = 0;
            } else if ((choice -= model_.alpha) < model_.gamma) {
                accepted = occupied(i, 0);
                changed_a = 0;
            } else if ((choice -= model_.gamma) < model_.delta) {
                accepted = !occupied(i, parameters_.L - 1);
                changed_a = parameters_.L - 1;
            } else {
                accepted = occupied(i, parameters_.L - 1);
                changed_a = parameters_.L - 1;
            }
        }
        if (accepted)
            apply_walker_event(i, changed_a, changed_b, current_change);
        else ++rejected_walker_events_;
    }

    void attempt_uniformized_reference(double choice) {
        const double scale = model_.heaviest_weight;
        const double right_end = static_cast<double>(model_.links) * model_.p0 * scale;
        const double left_end = right_end + static_cast<double>(model_.links) * model_.q0 * scale;
        bool accepted = false;
        std::size_t changed_a = 0;
        std::size_t changed_b = parameters_.L;
        double current_change = 0.0;
        if (choice < right_end) {
            const std::size_t bond = std::min( static_cast<std::size_t>(choice / (model_.p0 * scale)), model_.links - 1);
            const std::size_t neighbor = (bond + 1) % parameters_.L;
            accepted = occupied(reference_, bond) && !occupied(reference_, neighbor) && accept_class_rate(reference_.data(), bond, model_.right_weight);
            changed_a = bond;
            changed_b = neighbor;
            current_change = model_.dq;
        } else if (choice < left_end) {
            const std::size_t bond = std::min( static_cast<std::size_t>((choice - right_end) / (model_.q0 * scale)), model_.links - 1);
            const std::size_t neighbor = (bond + 1) % parameters_.L;
            accepted = !occupied(reference_, bond) && occupied(reference_, neighbor) && accept_class_rate(reference_.data(), bond, model_.left_weight);
            changed_a = bond;
            changed_b = neighbor;
            current_change = -model_.dq;
        } else if (!model_.periodic) {
            // Same closed-chain guard as the walker proposals above.
            choice -= left_end;
            if (choice < model_.alpha) {
                accepted = !occupied(reference_, 0);
                changed_a = 0;
            } else if ((choice -= model_.alpha) < model_.gamma) {
                accepted = occupied(reference_, 0);
                changed_a = 0;
            } else if ((choice -= model_.gamma) < model_.delta) {
                accepted = !occupied(reference_, parameters_.L - 1);
                changed_a = parameters_.L - 1;
            } else {
                accepted = occupied(reference_, parameters_.L - 1);
                changed_a = parameters_.L - 1;
            }
        }
        if (accepted)
            apply_reference_event(changed_a, changed_b, current_change);
        else ++rejected_reference_events_;
    }

    bool advance_dynamics_until(double end_time) {
        if (unbiased_) {
            dynamics_->advance_unbiased(end_time);
            return true;
        }
        return dynamics_->advance_weighted(end_time);
    }

    // Draw a fresh pending event time for every process from the current
    // clock.  Both rejection-free mechanisms keep these; uniformization
    // proposes globally and ignores them.
    void draw_event_times() {
        for (std::size_t i = 0; i < parameters_.M; ++i)
            event_times_[i] = time_ + exponential_wait(walker_rate(i));
        event_times_[parameters_.M] = monitored() ? time_ + exponential_wait(reference_rate()) : std::numeric_limits<double>::infinity();
    }

    // Apply one stored reference event to a single walker.  This is the
    // per-walker half of apply_reference_event, used when the walkers are
    // evolved separately and the reference has already moved on.
    void replay_reference_event(std::size_t i, const ReferenceEventRecord& event) {
        Walker& walker = walkers_[i];
        walker.update_weight(time_);
        int old_mismatches = occupied(i, event.changed_a) != event.before_a;
        if (event.changed_b < parameters_.L)
            old_mismatches += occupied(i, event.changed_b) != event.before_b;
        const int changed_count = event.changed_b < parameters_.L ? 2 : 1;
        const int overlap_change = 4 * old_mismatches - 2 * changed_count;
        walker.overlap += overlap_change;
        walker.potential += model_.k_scaled * overlap_change;
    }

    // Everything a block can change.  Saved before a reviewed block so that a
    // block which overshot its budget can be discarded exactly, down to the
    // generator, and replayed with a shorter step.
    void save_block_state() {
        saved_configs_ = configs_;
        saved_masks_ = masks_;
        saved_walkers_ = walkers_;
        saved_event_times_ = event_times_;
        saved_reference_ = reference_;
        saved_reference_masks_ = reference_masks_;
        std::copy(std::begin(reference_nR_), std::end(reference_nR_), std::begin(saved_reference_nR_));
        std::copy(std::begin(reference_nL_), std::end(reference_nL_), std::begin(saved_reference_nL_));
        saved_reference_activity_right_ = reference_activity_right_;
        saved_reference_activity_left_ = reference_activity_left_;
        saved_reference_current_ = reference_current_;
        saved_time_ = time_;
        saved_accepted_walker_events_ = accepted_walker_events_;
        saved_rejected_walker_events_ = rejected_walker_events_;
        saved_accepted_reference_events_ = accepted_reference_events_;
        saved_rejected_reference_events_ = rejected_reference_events_;
        saved_rng_ = rng_;
    }

    void restore_block_state() {
        configs_ = saved_configs_;
        masks_ = saved_masks_;
        walkers_ = saved_walkers_;
        event_times_ = saved_event_times_;
        reference_ = saved_reference_;
        reference_masks_ = saved_reference_masks_;
        std::copy(std::begin(saved_reference_nR_), std::end(saved_reference_nR_), std::begin(reference_nR_));
        std::copy(std::begin(saved_reference_nL_), std::end(saved_reference_nL_), std::begin(reference_nL_));
        reference_activity_right_ = saved_reference_activity_right_;
        reference_activity_left_ = saved_reference_activity_left_;
        reference_current_ = saved_reference_current_;
        time_ = saved_time_;
        accepted_walker_events_ = saved_accepted_walker_events_;
        rejected_walker_events_ = saved_rejected_walker_events_;
        accepted_reference_events_ = saved_accepted_reference_events_;
        rejected_reference_events_ = saved_rejected_reference_events_;
        rng_ = saved_rng_;
    }

    // The exact spread of the population's log-weights at the present time.
    double population_log_weight_range() {
        double lowest = std::numeric_limits<double>::infinity();
        double highest = -std::numeric_limits<double>::infinity();
        for (Walker& walker : walkers_) {
            walker.update_weight(time_);
            lowest = std::min(lowest, walker.log_weight);
            highest = std::max(highest, walker.log_weight);
        }
        return highest - lowest;
    }

    void clone_population(bool force_resampling) {
        std::vector<double>& log_weights = clone_log_weights_;
        for (std::size_t i = 0; i < parameters_.M; ++i) {
            walkers_[i].update_weight(time_);
            log_weights[i] = walkers_[i].log_weight;
        }
        const auto extremes = std::minmax_element( log_weights.begin(), log_weights.end());
        const double shift = *extremes.second;
        const double exact_range = *extremes.second - *extremes.first;
        ++resampling_checks_;

        // A certificate-triggered or explicitly scheduled examination may
        // find that the exact range is still small.  Re-anchor the upper bound
        // and return before any exponentials or ESS work.  Recording-time
        // calls bypass this test because their population averages require an
        // actual resampling.
        reanchor_weight_bound(exact_range);
        if (!force_resampling && exact_range < log_resample_range_) return;

        // The effective sample size, the cloning normalization, and the
        // resampling weights are all built from the same shifted exponentials.
        // They are needed only when resampling will actually happen.
        std::vector<double>& weights = clone_weights_;
        double total_weight = 0.0;
        double squared_weight = 0.0;
        for (std::size_t i = 0; i < parameters_.M; ++i) {
            const double weight = std::exp(log_weights[i] - shift);
            weights[i] = weight;
            total_weight += weight;
            squared_weight += weight * weight;
        }
        // Accumulate the degeneracy of actual resampling populations only.
        effective_sample_size_sum_ += total_weight * total_weight / squared_weight;
        ++effective_sample_size_samples_;

        log_normalization_ += shift + std::log( total_weight / static_cast<double>(parameters_.M));

        const double spacing = total_weight / static_cast<double>(parameters_.M);
        const double start = uniform_open() * spacing;
        std::vector<Walker>& new_walkers = spare_walkers_;
        std::vector<std::uint64_t>& new_configs = spare_configs_;
        std::vector<std::uint64_t>& new_masks = spare_masks_;
        std::size_t parent = 0, previous_parent = parameters_.M;
        std::size_t unique_parents = 0;
        double cumulative = weights[0];
        for (std::size_t slot = 0; slot < parameters_.M; ++slot) {
            const double point = start + static_cast<double>(slot) * spacing;
            while (point > cumulative && parent + 1 < parameters_.M)
                cumulative += weights[++parent];
            if (parent != previous_parent) {
                ++unique_parents;
                previous_parent = parent;
            }
            // A child is a bit-identical copy of its parent, so the parent's
            // allowed-transition masks and its cached nR, nL, overlap and
            // potential are already correct for the child.  Copying the masks
            // is cheaper than rebuilding them from the configuration, and it
            // removes a second full pass over the population.
            new_walkers[slot] = walkers_[parent];
            new_walkers[slot].log_weight = 0.0;
            new_walkers[slot].last_weight_time = time_;
            // With monitoring the potential is advanced incrementally on every
            // reference event, so rebuild it from the cached counts here.  That
            // keeps the cloning time a synchronisation point at which rounding
            // drift accumulated since the previous one is discarded.
            new_walkers[slot].potential = model_.counting_potential( new_walkers[slot].activity_right, new_walkers[slot].activity_left) + model_.k_scaled * static_cast<double>(new_walkers[slot].overlap);
            const auto from = static_cast<std::ptrdiff_t>(parent * words_);
            const auto to = static_cast<std::ptrdiff_t>(slot * words_);
            std::copy_n(configs_.begin() + from, words_, new_configs.begin() + to);
            std::copy_n(masks_.begin() + static_cast<std::ptrdiff_t>(mask_offset(parent, 0, 0)), mask_vectors * words_, new_masks.begin() + static_cast<std::ptrdiff_t>( mask_offset(slot, 0, 0)));
        }
        const auto replaced = static_cast<std::uint64_t>( parameters_.M - unique_parents);
        walkers_killed_ += replaced;
        extra_clones_ += replaced;
        // Swapping leaves the previous generation's storage in the scratch
        // buffers, so the next cloning reuses it and never allocates.
        walkers_.swap(new_walkers);
        configs_.swap(new_configs);
        masks_.swap(new_masks);
        // Every weight was just reset to one, so the range restarts at zero.
        reanchor_weight_bound(0.0);
        ++resampling_steps_;

        dynamics_->after_resampling();
    }

    Record output_record() const {
        Record row;
        row.time = time_;
        double mean_current = 0.0;
        double mean_potential = 0.0;
        for (std::size_t i = 0; i < parameters_.M; ++i) {
            mean_current += walkers_[i].current;
            mean_potential += walkers_[i].potential;
        }
        mean_current /= static_cast<double>(parameters_.M);
        mean_potential /= static_cast<double>(parameters_.M);
        row.mean_potential = mean_potential;
        for (const Walker& walker : walkers_)
            row.potential_spread = std::max( row.potential_spread, std::abs(walker.potential - mean_potential));
        if (time_ > 0.0) {
            // Every output follows an unconditional resampling, so all pending
            // log-weights are exactly zero. clone_population() has already
            // folded their previous log mean into this normalization.
            row.cgf_cloning = log_normalization_ / time_;
            row.accepted_current_rate = mean_current / time_;
            row.reference_current_rate = reference_current_ / time_;
        }
        row.effective_sample_size = effective_sample_size_samples_ > 0 ? effective_sample_size_sum_ / effective_sample_size_samples_ : static_cast<double>(parameters_.M);
        row.resampling_steps = resampling_steps_;
        row.resampling_checks = resampling_checks_;
        row.accepted_walker_events = accepted_walker_events_;
        row.rejected_walker_events = rejected_walker_events_;
        // Population profile at this exact, post-cloning recording time.
        // The strided read below looks cache-hostile, but it is branch-free and
        // vectorizes; an occupied-site-only sweep over each walker's words was
        // measured to be slower because of its data-dependent inner loop.
        row.density_profile.assign(parameters_.L, 0.0);
        for (std::size_t site = 0; site < parameters_.L; ++site) {
            std::size_t count = 0;
            for (std::size_t i = 0; i < parameters_.M; ++i)
                count += occupied(i, site);
            row.density_profile[site] = static_cast<double>(count) / parameters_.M;
        }
        return row;
    }

    void reset_recording_window() {
        effective_sample_size_sum_ = 0.0;
        effective_sample_size_samples_ = 0;
    }

    template <typename T> static void append(std::string& bytes, const T& value) {
        const char* raw = reinterpret_cast<const char*>(&value);
        bytes.append(raw, sizeof(T));
    }
    template <typename T> static void append_vector(std::string& bytes, const std::vector<T>& values) {
        append(bytes, static_cast<std::uint64_t>(values.size()));
        if (!values.empty()) bytes.append( reinterpret_cast<const char*>(values.data()), values.size() * sizeof(T));
    }
    static void append_string(std::string& bytes, const std::string& value) {
        append(bytes, static_cast<std::uint64_t>(value.size()));
        bytes.append(value);
    }
    template <typename T> static T read(const std::string& bytes, std::size_t& offset) {
        if (offset + sizeof(T) > bytes.size())
            throw std::invalid_argument("Truncated C++ checkpoint");
        T value;
        std::memcpy(&value, bytes.data() + offset, sizeof(T));
        offset += sizeof(T);
        return value;
    }
    template <typename T> static std::vector<T> read_vector(const std::string& bytes, std::size_t& offset) {
        const auto size = read<std::uint64_t>(bytes, offset);
        if (size > (bytes.size() - offset) / sizeof(T))
            throw std::invalid_argument("Invalid C++ checkpoint vector");
        std::vector<T> values(static_cast<std::size_t>(size));
        if (size) std::memcpy(values.data(), bytes.data() + offset, static_cast<std::size_t>(size) * sizeof(T));
        offset += static_cast<std::size_t>(size) * sizeof(T);
        return values;
    }
    static std::string read_string(const std::string& bytes, std::size_t& offset) {
        const auto size = read<std::uint64_t>(bytes, offset);
        if (size > bytes.size() - offset)
            throw std::invalid_argument("Invalid C++ checkpoint string");
        std::string value(bytes.data() + offset, static_cast<std::size_t>(size));
        offset += static_cast<std::size_t>(size);
        return value;
    }

    Parameters parameters_;
    Model model_;
    const bool unbiased_;
    std::size_t words_;
    std::uint64_t final_word_mask_;
    std::uint64_t final_bond_mask_;
    std::vector<Walker> walkers_;
    std::vector<std::uint64_t> configs_, masks_;
    std::vector<std::uint64_t> reference_, reference_masks_;
    std::size_t reference_nR_[classes] = {0, 0, 0, 0};
    std::size_t reference_nL_[classes] = {0, 0, 0, 0};
    double reference_activity_right_ = 0.0, reference_activity_left_ = 0.0;
    std::mt19937_64 rng_;
    // Binary heap held directly, rather than behind std::priority_queue, so
    // that a rebuild can reuse this storage instead of allocating a fresh
    // container at every cloning time.
    const std::uint64_t* overlap_reference_ = nullptr;
    // Block snapshot, allocated only when the chosen mechanism reviews blocks.
    std::vector<std::uint64_t> saved_configs_, saved_masks_;
    std::vector<Walker> saved_walkers_;
    std::vector<double> saved_event_times_;
    std::vector<std::uint64_t> saved_reference_, saved_reference_masks_;
    std::size_t saved_reference_nR_[classes] = {0, 0, 0, 0};
    std::size_t saved_reference_nL_[classes] = {0, 0, 0, 0};
    double saved_reference_activity_right_ = 0.0;
    double saved_reference_activity_left_ = 0.0;
    double saved_reference_current_ = 0.0;
    double saved_time_ = 0.0;
    std::uint64_t saved_accepted_walker_events_ = 0;
    std::uint64_t saved_rejected_walker_events_ = 0;
    std::uint64_t saved_accepted_reference_events_ = 0;
    std::uint64_t saved_rejected_reference_events_ = 0;
    std::mt19937_64 saved_rng_;
    std::uint64_t rejected_blocks_ = 0;
    double block_start_time_ = 0.0;
    std::unique_ptr<Dynamics> dynamics_;
    bool track_certificate_ = true;
    std::vector<double> event_times_;
    std::vector<double> clone_log_weights_, clone_weights_;
    std::vector<Walker> spare_walkers_;
    std::vector<std::uint64_t> spare_configs_, spare_masks_;
    double time_ = 0.0;
    double record_interval_ = 0.0;
    double estimated_jump_rate_ = 0.0;
    double next_clone_time_ = 0.0, next_output_time_ = 0.0;
    std::uint64_t output_number_ = 1, recordings_completed_ = 0;
    bool initial_row_emitted_ = false;
    double reference_current_ = 0.0;
    double log_normalization_ = 0.0;
    // Running certificate that the log-weight range cannot yet have reached
    // log(target_max).  ``weight_range_bound_`` integrates the widest
    // potential difference present in the population, so it is always an upper
    // bound on the true range; it is re-anchored to the exact range whenever
    // the population is actually examined.
    double log_target_max_ = 0.0;
    double log_resample_range_ = 0.0;
    double weight_range_bound_ = 0.0;
    double weight_bound_time_ = 0.0;
    double potential_maximum_ = 0.0, potential_minimum_ = 0.0;
    double effective_sample_size_sum_ = 0.0;
    std::uint64_t effective_sample_size_samples_ = 0;
    std::uint64_t resampling_checks_ = 0;
    std::uint64_t accepted_walker_events_ = 0;
    std::uint64_t rejected_walker_events_ = 0;
    std::uint64_t accepted_reference_events_ = 0;
    std::uint64_t rejected_reference_events_ = 0;
    std::uint64_t resampling_steps_ = 0;
    std::uint64_t walkers_killed_ = 0, extra_clones_ = 0;
};

// ---------------------------------------------------------------------------
// Heap dynamics: rejection-free events on one clock
//
// Every process -- each walker and the reference -- carries one pending event
// time, and a binary heap keeps the earliest of them at the front.  The
// population therefore advances in true chronological order, which is what
// lets the weight certificate be evaluated continuously and cloning be
// triggered at the exact instant the bound is reached.
// ---------------------------------------------------------------------------
class HeapDynamics final : public Dynamics {
public:
    explicit HeapDynamics(Simulation& simulation) : sim_(simulation) {
        heap_.reserve(sim_.parameters_.M + 1);
    }

    bool advance_weighted(double end_time) override {
        while (!heap_.empty() && heap_.front().time <= end_time) {
            step(heap_.front());
            if (sim_.weight_bound_reached()) return false;
        }
        sim_.time_ = end_time;
        return true;
    }

    // Same ordering, with none of the weight machinery: at s = k = 0 every
    // weight stays one, so there is no certificate to reach and no reason to
    // test for it.
    void advance_unbiased(double end_time) override {
        while (!heap_.empty() && heap_.front().time <= end_time)
            step(heap_.front());
        sim_.time_ = end_time;
    }

    void after_resampling() override {
        sim_.draw_event_times();
        rebuild_from_event_times();
    }

    void restore_schedule() override {
        if (sim_.event_times_.size() != sim_.parameters_.M + 1)
            throw std::invalid_argument("Checkpoint event-time count mismatch");
        rebuild_from_event_times();
    }

private:
    void step(const Event event) {
        sim_.time_ = event.time;
        if (event.process == sim_.parameters_.M)
            sim_.execute_reference_event();
        else sim_.execute_walker_event(event.process);
        const double rate = event.process == sim_.parameters_.M ? sim_.reference_rate() : sim_.walker_rate(event.process);
        sim_.event_times_[event.process] = sim_.time_ + sim_.exponential_wait(rate);
        replace_root({sim_.event_times_[event.process], event.process});
    }

    void rebuild_from_event_times() {
        // One O(M) heapify, instead of M separate O(log M) insertions.  The
        // storage reserved in the constructor is reused, so the hundreds of
        // thousands of rebuilds in a run perform no allocation at all.
        heap_.clear();
        for (std::size_t process = 0; process < sim_.event_times_.size(); ++process)
            heap_.push_back({sim_.event_times_[process], process});
        std::make_heap(heap_.begin(), heap_.end(), LaterEvent{});
    }

    void replace_root(const Event replacement) {
        // The event just executed was the root and its replacement cannot be
        // earlier than the current time. Restore the fixed-size min-heap with
        // one downward sift instead of two logarithmic heap traversals.
        std::size_t hole = 0;
        while (true) {
            const std::size_t left = 2 * hole + 1;
            if (left >= heap_.size()) break;
            const std::size_t right = left + 1;
            std::size_t earlier_child = left;
            if (right < heap_.size() && LaterEvent{}(heap_[left], heap_[right])) {
                earlier_child = right;
            }
            if (!LaterEvent{}(replacement, heap_[earlier_child])) break;
            heap_[hole] = heap_[earlier_child];
            hole = earlier_child;
        }
        heap_[hole] = replacement;
    }

    Simulation& sim_;
    std::vector<Event> heap_;
};

// ---------------------------------------------------------------------------
// Uniformized dynamics: blind proposals, rejected down to the true rates
//
// One global clock ticks at the heaviest rate the population could possibly
// have; each tick picks a process and a transition uniformly and then rejects
// unless the exclusion test and the class rate allow it.  Slower than the heap
// by construction, and kept because it shares nothing with it but the model:
// agreement between the two is evidence about the model, not about one loop.
// ---------------------------------------------------------------------------
class UniformizedDynamics final : public Dynamics {
public:
    explicit UniformizedDynamics(Simulation& simulation) : sim_(simulation) {}

    bool advance_weighted(double end_time) override {
        const Rates rates = proposal_rates();
        while (true) {
            const double proposal_time = sim_.time_ + sim_.exponential_wait(rates.global);
            if (proposal_time > end_time) break;
            sim_.time_ = proposal_time;
            propose(rates);
            if (sim_.weight_bound_reached()) return false;
        }
        sim_.time_ = end_time;
        return true;
    }

    void advance_unbiased(double end_time) override {
        const Rates rates = proposal_rates();
        while (true) {
            const double proposal_time = sim_.time_ + sim_.exponential_wait(rates.global);
            if (proposal_time > end_time) break;
            sim_.time_ = proposal_time;
            propose(rates);
        }
        sim_.time_ = end_time;
    }

private:
    struct Rates {
        double walker;
        double walker_population;
        double global;
    };

    Rates proposal_rates() const {
        const Model& model = sim_.model_;
        const double boundary_rate = model.periodic ? 0.0 : model.alpha + model.gamma + model.delta + model.beta;
        const double walker = static_cast<double>(model.links)
                                  * model.heaviest_weight
                                  * (model.ps + model.qs)
                            + boundary_rate;
        const double reference = sim_.monitored()
                               ? static_cast<double>(model.links)
                                     * model.heaviest_weight
                                     * (model.p0 + model.q0)
                                 + boundary_rate
                               : 0.0;
        const double population = static_cast<double>(sim_.parameters_.M) * walker;
        return {walker, population, population + reference};
    }

    void propose(const Rates& rates) {
        double position = sim_.uniform_open() * rates.global;
        if (position < rates.walker_population) {
            const std::size_t walker = std::min( static_cast<std::size_t>(position / rates.walker), sim_.parameters_.M - 1);
            position -= static_cast<double>(walker) * rates.walker;
            sim_.attempt_uniformized_walker(walker, position);
        } else {
            sim_.attempt_uniformized_reference(position - rates.walker_population);
        }
    }

    Simulation& sim_;
};

// ---------------------------------------------------------------------------
// Direct dynamics: one walker at a time, run alone to the next stop
//
// Between two cloning or recording instants the walkers do not interact, so
// there is nothing to gain from interleaving them in global time order and a
// great deal to lose: the heap costs a sift per event and, worse, brings a
// different walker's configuration and masks into cache on every step.  Here
// each walker is carried alone from its own pending event time to the end of
// the block, so its few hundred bytes of state stay in L1 for thousands of
// consecutive events.
//
// The price is that no running certificate over the population exists -- while
// walker 3 is at t = 7, walker 800 is still at t = 2 -- so the cloning instants
// cannot be discovered as the run proceeds and are fixed in advance instead;
// see Simulation::schedule_next_cloning_check.
//
// With monitoring the walkers are coupled to the reference, which is why it is
// evolved first and its events are stored: each walker then replays that one
// realisation as it goes, merging it with its own events in time order.
// ---------------------------------------------------------------------------
class DirectDynamics final : public Dynamics {
public:
    explicit DirectDynamics(Simulation& simulation) : sim_(simulation) {}

    // The first block has no measurement behind it, so it uses the rigorous
    // bound: over any interval the log-weight ratio cannot open faster than
    // the widest potential spread the model permits, so log(target_max)
    // divided by that spread is safe whatever the configurations are.
    void prepare() override {
        const double widest = sim_.model_.widest_possible_potential_spread();
        step_ = widest > 0.0 ? sim_.log_target_max() / widest : std::numeric_limits<double>::infinity();
    }

    bool advance_weighted(double end_time) override {
        const bool monitored = sim_.monitored();
        reference_events_.clear();
        // The reference must be captured before it moves: each walker replays
        // the block against this snapshot, so that its overlap is always
        // measured against the reference at the walker's own time rather than
        // at the block's end.
        if (monitored) reference_at_block_start_ = sim_.reference_;
        if (monitored) advance_reference(end_time, monitored);
        for (std::size_t i = 0; i < sim_.parameters_.M; ++i) {
            if (monitored) {
                shadow_reference_ = reference_at_block_start_;
                sim_.set_overlap_reference(shadow_reference_.data());
                advance_walker_with_reference(i, end_time);
            } else {
                advance_walker(i, end_time);
            }
            // Close the walker's weight integral at the block boundary, so the
            // population is comparable the moment the block ends.
            sim_.walkers_[i].update_weight(end_time);
        }
        if (monitored) sim_.set_overlap_reference(nullptr);
        sim_.time_ = end_time;
        return true;
    }

    // At s = k = 0 there is no weight to integrate and the reference no longer
    // touches the walkers, so this is plain independent trajectory sampling.
    void advance_unbiased(double end_time) override {
        // Unbiased means s = k = 0, so the reference has no reader here.
        for (std::size_t i = 0; i < sim_.parameters_.M; ++i)
            advance_walker(i, end_time);
        sim_.time_ = end_time;
    }

    void after_resampling() override { sim_.draw_event_times(); }

    double step_size() const override { return step_; }
    bool reviews_blocks() const override { return true; }
    double schedule_state() const override { return step_; }
    void set_schedule_state(double value) override {
        if (value > 0.0) step_ = value;
    }

    // Adaptive step control.  A block is judged by the exact log-weight range
    // it produced.  Above target_max the prediction was wrong in the dangerous
    // direction and the block is discarded and replayed at half the step.
    // Otherwise the block stands and the step is re-estimated from the rate at
    // which the range actually opened, which -- unlike the range itself -- is
    // unaffected by a block cut short by a recording.
    //
    // Growth is capped at twice the accepted step so that one lucky block
    // cannot stretch the next one out of range: shrink fast, grow slowly.  A
    // floor stops the halving from chasing a pathological single event for
    // ever; a block that overshoots there is accepted and counted instead.
    bool review_block(double log_weight_range, double duration) override {
        if (log_weight_range > sim_.log_target_max() && step_ > minimum_step()) {
            step_ = std::max(0.5 * step_, minimum_step());
            return true;
        }
        const double rate = duration > 0.0 ? log_weight_range / duration : 0.0;
        const double predicted = rate > 0.0 ? sim_.log_target_min() / rate : std::numeric_limits<double>::infinity();
        step_ = std::min(predicted, 2.0 * step_);
        return false;
    }

private:
    // Each walker keeps its own pending event time across block boundaries, so
    // stopping at a recording instant costs nothing: the walker is left with
    // its next event beyond the block, exactly where the heap would leave it,
    // and the recorded snapshot is its true state at ``end_time``.
    void advance_walker(std::size_t i, double end_time) {
        double next = sim_.event_times_[i];
        while (next <= end_time) {
            sim_.time_ = next;
            sim_.execute_walker_event(i);
            next = sim_.time_ + sim_.exponential_wait(sim_.walker_rate(i));
        }
        sim_.event_times_[i] = next;
    }

    // The merge of two ordered streams: this walker's own events and the
    // reference events recorded for this block.
    void advance_walker_with_reference(std::size_t i, double end_time) {
        double next = sim_.event_times_[i];
        std::size_t replayed = 0;
        while (true) {
            const double reference_time = replayed < reference_events_.size() ? reference_events_[replayed].time : std::numeric_limits<double>::infinity();
            if (next <= end_time && next <= reference_time) {
                sim_.time_ = next;
                sim_.execute_walker_event(i);
                next = sim_.time_ + sim_.exponential_wait(sim_.walker_rate(i));
            } else if (reference_time <= end_time) {
                const ReferenceEventRecord& event = reference_events_[replayed];
                sim_.time_ = reference_time;
                sim_.replay_reference_event(i, event);
                // Advance the shadow with it, so the walker's next event sees
                // the reference as it stands at that moment.
                flip_shadow(event.changed_a);
                if (event.changed_b < sim_.parameters_.L)
                    flip_shadow(event.changed_b);
                ++replayed;
            } else {
                break;
            }
        }
        sim_.event_times_[i] = next;
    }

    void advance_reference(double end_time, bool record) {
        const std::size_t process = sim_.parameters_.M;
        double next = sim_.event_times_[process];
        while (next <= end_time) {
            sim_.time_ = next;
            if (record) sim_.execute_reference_event(&reference_events_);
            else sim_.execute_reference_event();
            next = sim_.time_ + sim_.exponential_wait(sim_.reference_rate());
        }
        sim_.event_times_[process] = next;
    }

    void flip_shadow(std::size_t site) {
        shadow_reference_[site / 64] ^= std::uint64_t{1} << (site % 64);
    }

    // Halving stops here: a few mean waiting times of a single walker, below
    // which a shorter block cannot separate the weights any less.
    double minimum_step() const {
        const double rate = sim_.model_.broad_jump_rate_estimate();
        return rate > 0.0 ? 1.0 / rate : 0.0;
    }

    Simulation& sim_;
    std::vector<ReferenceEventRecord> reference_events_;
    std::vector<std::uint64_t> reference_at_block_start_;
    std::vector<std::uint64_t> shadow_reference_;
    double step_ = 0.0;
};

std::unique_ptr<Dynamics> make_dynamics(Simulation& simulation, DynamicsKind kind) {
    switch (kind) {
        case DynamicsKind::Uniformized:
            return std::make_unique<UniformizedDynamics>(simulation);
        case DynamicsKind::Direct:
            return std::make_unique<DirectDynamics>(simulation);
        case DynamicsKind::Heap:
        default:
            return std::make_unique<HeapDynamics>(simulation);
    }
}

} // namespace monitored_exclusion

using monitored_exclusion::DynamicsKind;
using monitored_exclusion::Parameters;
using monitored_exclusion::Record;
using monitored_exclusion::Simulation;

static const char* capsule_name = "Python_Cpp_Interface.NNN";

static PyObject* dict_item(PyObject* dictionary, const char* name) {
    PyObject* value = PyDict_GetItemString(dictionary, name);
    if (!value) {
        PyErr_Format(PyExc_KeyError, "Missing parameter: %s", name);
        throw std::invalid_argument("missing parameter");
    }
    return value;
}

static double as_double(PyObject* dictionary, const char* name) {
    const double value = PyFloat_AsDouble(dict_item(dictionary, name));
    if (PyErr_Occurred()) throw std::invalid_argument("invalid float parameter");
    return value;
}

static std::uint64_t as_uint64(PyObject* dictionary, const char* name) {
    const auto value = PyLong_AsUnsignedLongLong(dict_item(dictionary, name));
    if (PyErr_Occurred()) throw std::invalid_argument("invalid integer parameter");
    return value;
}

static std::string as_string(PyObject* dictionary, const char* name) {
    const char* value = PyUnicode_AsUTF8(dict_item(dictionary, name));
    if (!value) throw std::invalid_argument("invalid string parameter");
    return value;
}

static std::vector<double> as_double_list(PyObject* dictionary, const char* name) {
    PyObject* value = dict_item(dictionary, name);
    if (value == Py_None) return {};
    PyObject* sequence = PySequence_Fast(value, "expected a sequence of rates");
    if (!sequence) throw std::invalid_argument("invalid rate list");
    std::vector<double> values;
    const Py_ssize_t size = PySequence_Fast_GET_SIZE(sequence);
    for (Py_ssize_t i = 0; i < size; ++i)
        values.push_back( PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence, i)));
    Py_DECREF(sequence);
    if (PyErr_Occurred()) throw std::invalid_argument("invalid rate value");
    return values;
}

static Parameters parse_parameters(PyObject* dictionary) {
    if (!PyDict_Check(dictionary))
        throw std::invalid_argument("parameters must be a dictionary");
    Parameters p;
    p.L = as_uint64(dictionary, "L");
    p.M = as_uint64(dictionary, "M");
    p.tmax = as_double(dictionary, "tmax");
    p.E = as_double(dictionary, "E");
    p.s = as_double(dictionary, "s");
    p.k = as_double(dictionary, "k");
    p.model = as_string(dictionary, "model");
    p.epsilon = as_double(dictionary, "epsilon");
    p.delta_kls = as_double(dictionary, "delta_kls");
    p.right_weights = as_double_list(dictionary, "right_weights");
    p.left_weights = as_double_list(dictionary, "left_weights");
    p.filling = as_double(dictionary, "filling");
    p.alpha = as_double(dictionary, "alpha");
    p.gamma = as_double(dictionary, "gamma");
    p.delta = as_double(dictionary, "delta");
    p.beta = as_double(dictionary, "beta");
    p.target_min = as_double(dictionary, "target_min");
    p.target_max = as_double(dictionary, "target_max");
    PyObject* clone = dict_item(dictionary, "cloning_interval");
    p.cloning_interval = clone == Py_None ? std::numeric_limits<double>::infinity() : PyFloat_AsDouble(clone);
    PyObject* record = dict_item(dictionary, "record_interval");
    p.record_interval = record == Py_None ? -1.0 : PyFloat_AsDouble(record);
    p.seed = as_uint64(dictionary, "seed");
    p.initial = as_string(dictionary, "initial");
    p.boundary = as_string(dictionary, "boundary");
    const std::string dynamics = as_string(dictionary, "dynamics");
    if (dynamics == "heap")
        p.dynamics = DynamicsKind::Heap;
    else if (dynamics == "uniformized")
        p.dynamics = DynamicsKind::Uniformized;
    else if (dynamics == "direct")
        p.dynamics = DynamicsKind::Direct;
    else throw std::invalid_argument( "dynamics must be heap, uniformized or direct");
    PyObject* progress = dict_item(dictionary, "progress");
    p.progress = progress == Py_None ? 0 : PyLong_AsUnsignedLongLong(progress);
    if (PyErr_Occurred()) throw std::invalid_argument("invalid optional parameter");
    return p;
}

static Simulation* get_simulation(PyObject* capsule) {
    auto* simulation = static_cast<Simulation*>( PyCapsule_GetPointer(capsule, capsule_name));
    if (!simulation) throw std::invalid_argument("invalid simulation capsule");
    return simulation;
}

static void destroy_capsule(PyObject* capsule) {
    void* pointer = PyCapsule_GetPointer(capsule, capsule_name);
    if (pointer) delete static_cast<Simulation*>(pointer);
    else PyErr_Clear();
}

static PyObject* py_create(PyObject*, PyObject* dictionary) {
    try {
        auto* simulation = new Simulation(parse_parameters(dictionary));
        return PyCapsule_New(simulation, capsule_name, destroy_capsule);
    } catch (const std::exception& error) {
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, error.what());
        return nullptr;
    }
}

static bool dict_set(PyObject* dict, const char* name, PyObject* value) {
    if (!value) return false;
    const int result = PyDict_SetItemString(dict, name, value);
    Py_DECREF(value);
    return result == 0;
}

static bool dict_set_vector(PyObject* dict, const char* name, const std::vector<double>& values) {
    PyObject* list = PyList_New(static_cast<Py_ssize_t>(values.size()));
    if (!list) return false;
    for (std::size_t i = 0; i < values.size(); ++i) {
        PyObject* value = PyFloat_FromDouble(values[i]);
        if (!value) {
            Py_DECREF(list);
            return false;
        }
        PyList_SET_ITEM(list, static_cast<Py_ssize_t>(i), value);
    }
    const int result = PyDict_SetItemString(dict, name, list);
    Py_DECREF(list);
    return result == 0;
}

static PyObject* record_to_python(const Record& row) {
    PyObject* result = PyDict_New();
    if (!result) return nullptr;
    const bool ok = dict_set(result, "Time", PyFloat_FromDouble(row.time)) && dict_set(result, "CGFCloningNormalization", PyFloat_FromDouble(row.cgf_cloning)) && dict_set(result, "MeanPopulationPotential", PyFloat_FromDouble(row.mean_potential)) && dict_set(result, "AcceptedCurrentRate", PyFloat_FromDouble(row.accepted_current_rate)) && dict_set(result, "ReferenceCurrentRate", PyFloat_FromDouble(row.reference_current_rate)) && dict_set(result, "PotentialSpread", PyFloat_FromDouble(row.potential_spread)) && dict_set(result, "EffectiveSampleSize", PyFloat_FromDouble(row.effective_sample_size)) && dict_set(result, "ResamplingSteps", PyLong_FromUnsignedLongLong(row.resampling_steps)) && dict_set(result, "ResamplingChecks", PyLong_FromUnsignedLongLong(row.resampling_checks)) && dict_set(result, "AcceptedWalkerEvents", PyLong_FromUnsignedLongLong(row.accepted_walker_events)) && dict_set(result, "RejectedWalkerEvents", PyLong_FromUnsignedLongLong(row.rejected_walker_events));
    if (!ok) { Py_DECREF(result); return nullptr; }
    if (!dict_set_vector(result, "DensityProfile", row.density_profile)) {
        Py_DECREF(result);
        return nullptr;
    }
    return result;
}

static PyObject* py_run(PyObject*, PyObject* capsule) {
    try {
        Simulation* simulation = get_simulation(capsule);
        std::vector<Record> rows;
        PyThreadState* thread_state = PyEval_SaveThread();
        try {
            rows = simulation->run();
        } catch (...) {
            PyEval_RestoreThread(thread_state);
            throw;
        }
        PyEval_RestoreThread(thread_state);
        PyObject* result = PyList_New(static_cast<Py_ssize_t>(rows.size()));
        if (!result) return nullptr;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            PyObject* item = record_to_python(rows[i]);
            if (!item) { Py_DECREF(result); return nullptr; }
            PyList_SET_ITEM(result, static_cast<Py_ssize_t>(i), item);
        }
        return result;
    } catch (const std::exception& error) {
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

static PyObject* py_checkpoint(PyObject*, PyObject* capsule) {
    try {
        const std::string bytes = get_simulation(capsule)->checkpoint();
        return PyBytes_FromStringAndSize(bytes.data(), static_cast<Py_ssize_t>(bytes.size()));
    } catch (const std::exception& error) {
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

static PyObject* py_restore(PyObject*, PyObject* args) {
    PyObject* capsule;
    const char* bytes;
    Py_ssize_t size;
    if (!PyArg_ParseTuple(args, "Oy#", &capsule, &bytes, &size)) return nullptr;
    try {
        get_simulation(capsule)->restore(std::string(bytes, size));
        Py_RETURN_NONE;
    } catch (const std::exception& error) {
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, error.what());
        return nullptr;
    }
}

static PyObject* py_info(PyObject*, PyObject* capsule) {
    try {
        Simulation* simulation = get_simulation(capsule);
        PyObject* result = PyDict_New();
        if (!result) return nullptr;
        const bool ok = dict_set(result, "time", PyFloat_FromDouble(simulation->time())) && dict_set(result, "record_interval", PyFloat_FromDouble(simulation->record_interval())) && dict_set(result, "estimated_jump_rate", PyFloat_FromDouble(simulation->estimated_jump_rate())) && dict_set(result, "unbiased", PyBool_FromLong(simulation->unbiased())) && dict_set(result, "dynamics", PyUnicode_FromString( simulation->dynamics_name())) && dict_set(result, "boundary", PyUnicode_FromString( simulation->periodic() ? "periodic" : "open")) && dict_set(result, "accepted_walker_events", PyLong_FromUnsignedLongLong(simulation->accepted_walker_events())) && dict_set(result, "rejected_walker_events", PyLong_FromUnsignedLongLong(simulation->rejected_walker_events())) && dict_set(result, "accepted_reference_events", PyLong_FromUnsignedLongLong(simulation->accepted_reference_events())) && dict_set(result, "rejected_reference_events", PyLong_FromUnsignedLongLong(simulation->rejected_reference_events())) && dict_set(result, "resampling_steps", PyLong_FromUnsignedLongLong(simulation->resampling_steps())) && dict_set(result, "resampling_checks", PyLong_FromUnsignedLongLong(simulation->resampling_checks())) && dict_set(result, "rejected_blocks", PyLong_FromUnsignedLongLong(simulation->rejected_blocks())) && dict_set(result, "target_min", PyFloat_FromDouble(simulation->target_min())) && dict_set(result, "target_max", PyFloat_FromDouble(simulation->target_max())) && dict_set(result, "rate_model", PyUnicode_FromString( simulation->rate_model().c_str())) && dict_set(result, "walkers_killed", PyLong_FromUnsignedLongLong(simulation->walkers_killed()));
        if (!ok) { Py_DECREF(result); return nullptr; }
        if (!dict_set_vector(result, "right_weights", simulation->right_weights()) || !dict_set_vector(result, "left_weights", simulation->left_weights())) {
            Py_DECREF(result);
            return nullptr;
        }
        return result;
    } catch (const std::exception& error) {
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

static PyMethodDef methods[] = {
    {"create", reinterpret_cast<PyCFunction>(py_create), METH_O, "Create a C++ rejection-free simulation."}, {"run", reinterpret_cast<PyCFunction>(py_run), METH_O, "Run the simulation and return new records."}, {"checkpoint", reinterpret_cast<PyCFunction>(py_checkpoint), METH_O, "Return the opaque C++ checkpoint bytes."}, {"restore", reinterpret_cast<PyCFunction>(py_restore), METH_VARARGS, "Restore opaque C++ checkpoint bytes."}, {"info", reinterpret_cast<PyCFunction>(py_info), METH_O, "Return current engine diagnostics."}, {nullptr, nullptr, 0, nullptr}
};

static PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_engine", "C++ rejection-free monitored exclusion simulation engine.", -1, methods };

PyMODINIT_FUNC PyInit__engine() { return PyModule_Create(&module); }
