# NNN — next-nearest-neighbour exclusion engine

**A Monte Carlo simulator for rare events in stochastic lattice models** — a C++ engine behind a Python interface, with a bit-parallel state representation, three interchangeable event-generation implementations, and checkpointed runs that resume bit-for-bit.

**The physics it is built for**: a classical lattice gas is watched continuously by an observer who, from the information collected, tries to describe the statistics of the particle current. The whole distribution is sought, down to exponentially rare values, which are encoded in the cumulant generating function of the current. The models considered are exclusion processes whose hop rates depend on the sites adjacent to a bond; the weakly asymmetric exclusion process (WASEP) and the Katz–Lebowitz–Spohn model (KLS, [Katz, Lebowitz & Spohn 1984](https://doi.org/10.1007/BF01018556)) are two particular instances.
The observer never sees the configuration itself, only the measurement record, and so maintains a probability distribution over configurations that represents their knowledge of the system, updating it by Bayesian conditioning as the record arrives. This is exactly the classical problem of estimating the hidden state of a Markov chain from noisy observations — the *filtering* problem ([Bain & Crisan 2009](https://doi.org/10.1007/978-0-387-76896-0); [Doucet & Johansen 2011](https://www.stats.ox.ac.uk/~doucet/doucet_johansen_tutorialPF2011.pdf)). The measurement procedure is standard, inheriting its structure from the quantum measurement formalism described in [Jacobs & Steck 2006](https://doi.org/10.1080/00107510601101934); the single-particle counterpart of what is simulated here is the monitored random walker of [Jin & Martin 2022](https://doi.org/10.1103/PhysRevLett.129.260603). See the notes [`theory_algorithm_implementation.tex`](theory_algorithm_implementation.tex) for further details.

The theory, the algorithm and the implementation are documented in [`theory_algorithm_implementation.tex`](theory_algorithm_implementation.tex), with its references in [`algorithm_refs.bib`](algorithm_refs.bib). This README covers simply how to build the code, run it and extract the data, discussing a few generalities about the implementation.

## What's in it

- A **C++20 simulation engine** (~2500 lines) compiled as a CPython extension module and driven entirely from Python.
- A **bit-parallel state representation**: configurations are packed into 64-bit words and the allowed transitions are kept as bitmasks, so counting them is a `popcount` per word. Choosing which transition occurs then amounts to locating the `n`-th set bit (the `n`-th bit equal to `1`) of a bitmask, with `n` drawn uniformly among the allowed transitions. That is precisely what a single `PDEP` instruction does.
- **Runtime CPU dispatch.** `PDEP` is used only on processors where it is actually fast — it is microcoded, and hundreds of cycles slower on AMD Zen 1 and Zen 2, for instance. Elsewhere the fallback is a portable loop that clears (sets back to `0`) the lowest `n` set bits and then counts the trailing zeros. The use of `PDEP` is decided once the CPython extension module is loaded. Both paths return the same bond for the same draw, so results are unaffected.
- **Compile-time CPU baseline.** Counting the allowed transitions rests on `std::popcount`. Compiled with the `-march=x86-64-v2` flag, it becomes a single hardware `POPCNT` instruction; without it, the compiler emits a software routine that counts the bits one by one, giving the same answer far more slowly. Unlike the choice of `PDEP` above, this one is fixed in the binary at compile time rather than made when it loads.
- **Three interchangeable event-generation algorithms** that implement the same process and are cross-validated against each other.
- **Resumable runs.** A checkpoint stores the random generator's state along with the population, so a finished run can be extended. The extended run is bitwise identical to one carried out from the start without any interruption.
- **Parameter sweeps from a TOML plan.** [`simulations.toml`](simulations.toml) allows several simulations to be launched at once. It declares the parameters shared by every run, the ones to sweep over, and any correlated sets of values; the launcher collects all runs to be submitted, checks that no two jobs would write the same file, and runs them through one single-threaded process per available core.

## Build and run

Build the simulator:
```bash
conda env create -f environment.yml
conda activate monitored-exclusion
cd Simulator && python setup_cpp.py build_ext --inplace && cd ..
```

[`Simulator/setup_cpp.py`](Simulator/setup_cpp.py) must be run from `Simulator/` because setuptools resolves the package it copies the C++ extension module into against the working directory. Everything else is independent from where it is started.

This build passes `-march=x86-64-v2`, so it assumes any non-prehistoric x86-64 machine (roughly, Intel from Nehalem in 2008 and AMD from Bulldozer in 2011 onwards). On an older one the resulting binary aborts with an illegal instruction, and on a non-x86 machine, such as an ARM laptop, the compiler will not accept the flag at all. Both cases are fixed by removing that single flag and rebuilding, as detailed in [`Simulator/docs/CPU_REQUIREMENTS.md`](Simulator/docs/CPU_REQUIREMENTS.md).

To launch the code and begin the simulation, edit [`simulations.toml`](simulations.toml) (check out the documentation in [`LAUNCHER.md`](LAUNCHER.md)) and run `bash run_NNN.sh` with the environment activated. 

A single run, without the launcher:

```bash
python Simulator/run_NNN_cpp.py -L 100 -M 2000 -t 200 -E 7 -s -3.35 --epsilon 0.6 --boundary open --alpha 1 --gamma 1 --delta 1 --beta 1
```

## Algorithm, in one paragraph

The code simulates a Markov process in which each trajectory — which we call a *walker* — carries a statistical weight, set either by the measurement outcomes or by the bias imposed to compute the current cumulant generating function (see [`theory_algorithm_implementation.tex`](theory_algorithm_implementation.tex), Section *Large deviations*, for a comprehensive introduction to the theory). The Markov dynamics are generated by the Gillespie algorithm, and the weights are handled by the cloning method of [Giardinà, Kurchan & Peliti 2006](https://doi.org/10.1103/PhysRevLett.96.120603) in the continuous-time formulation of [Lecomte & Tailleur 2007](https://doi.org/10.1088/1742-5468/2007/03/P03004): pure Markov evolution is interleaved with resampling instants at which each walker is cloned or killed according to its weight. When measurements are considered there is in addition a reference trajectory, representing the true evolution of the system, upon which the measurements are performed; it couples to the walkers and shapes the observer's statistical description, which is encoded in the distribution of the walkers.

In what follows, an *event* is a configuration jump induced by the Markov evolution. Three interchangeable implementations generate those events, selected by `dynamics`:

| `dynamics`    | how events are generated                                                           |  speed (system size `L = 100`)           |
| ------------- | -----------------------------------------------------------------------------------| -----------------------------------------|
| `direct`      | rejection-free; each walker runs independently between resampling stops            | fastest; the default choice              |
| `heap`        | rejection-free; all processes on one clock, ordered by a binary heap               | 1.2–2.8× slower than `direct`            |
| `uniformized` | blind proposals of any jump between adjacent sites, rejected down to the true rate | simple but 2.7–5.1× slower than `direct` |

All three produce the same physics. `uniformized` is kept because its implementation is simple enough to be obviously correct, which makes it a good benchmark for the other two; `heap` because it is a cute idea, and one that may be useful in other settings; and `direct` because it is the most efficient. The full account — the algorithm as a sequence of steps, the resampling schedule, the implementation, the measured cost of each mechanism and the comparison between them — is in [`theory_algorithm_implementation.tex`](theory_algorithm_implementation.tex), Sections *The algorithm* and *Implementation*.

## Validation

Two independent kinds of check, in increasing order of what they can catch.

**Internal consistency.** The three mechanisms implement the same process, so they must agree, and they are compared rather than assumed to. `heap` and `uniformized` reproduce the earlier single-model engine's output bitwise. `direct` consumes the random stream in its own order, so it agrees in distribution rather than per seed — within one standard error over eight seeds in the cases checked, monitored ones included — and reproduces its own runs exactly across a checkpoint. `KLS` and `NNN` agree event for event at equal rates and equal seed, and `model = "WASEP"` reproduces the retired WASEP engine bitwise.

**Exact results.** Two references exist for the cumulant generating function itself, that we use to test the engine against:

- *Derrida's open SSEP.* The exact analytical result to the open SSEP's CGF in the limit of infinite `L`. [`validation/make_convergence_figure.py`](validation/make_convergence_figure.py) plots the engine against it over a range of sizes: the difference between the two falls as `1/L`, from 11 per cent at `L = 8` to 0.6 per cent at `L = 128`, confirming convergence of the engine's results to the correct asymptotic ones.
- *Transfer matrix.* [`validation/exact/SSEP_TM.py`](validation/exact/SSEP_TM.py) diagonalises the tilted generator of a small open SSEP, giving the CGF with no statistical error at all for `L` up to about 14. All three mechanisms match it to better than `1e-4` at `L = 5, 6, 8`. See [`validation/README.md`](validation/README.md).


## Model selection in [`simulations.toml`](simulations.toml)

The `model` value chooses the exclusion process to be simulated:

| `model`   | model                                | bookkeeping                    | rates from                      | files      |
| --------- | ------------------------------------ | ------------------------------ | ------------------------------- | ---------- |
| `WASEP`   | weakly asymmetric exclusion process  | one pair of bond masks         | none; one rate per direction    | `WASEP_*`  |
| `KLS`     | Katz–Lebowitz–Spohn                  | four class masks per direction | `epsilon`, `delta_kls`          | `KLS_*`    |
| `NNN`     | generic next-nearest-neighbour model | four class masks per direction | `right_weights`, `left_weights` | `NNN_*`    |

One engine serves all three. `NNN` is the general case: hops occur only between adjacent sites, subject to exclusion, but their rate may depend on the two sites flanking the pair. Those two sites are either empty or occupied, giving `2 × 2 = 4` possibilities, and a hop may go left or right, so there are `4 × 2 = 8` rates in total. All eight must be inserted explicitly; [`LAUNCHER.md`](LAUNCHER.md) gives their ordering and how to write them into [`simulations.toml`](simulations.toml). 

WASEP and KLS are two instances of this family, and both have been important for my own work, so each gets its own flag: choosing one of them means supplying only that model's parameters rather than the eight rates. `KLS` and `NNN` share the same implementation, and at equal rates and equal seed they agree event for event. `WASEP`, on the other hand, is a genuinely different code path rather than the special case of the former which makes it about twice as fast. Why, and at what cost, is discussed in the notes.

## Boundaries

`boundary` selects between a closed ring and a chain coupled to two reservoirs:

| `boundary`   | what it means                                                        | relevant parameters              |
| ------------ | -------------------------------------------------------------------- | -------------------------------- |
| `"periodic"` | a ring of `L` bonds; the particle number is conserved                | `filling`                        |
| `"open"`     | a chain of `L-1` bonds with a reservoir at each end                  | `alpha`, `gamma`, `delta`, `beta`|

On a ring, `filling` fixes the initial density and is rounded to the nearest integer particle number; the reservoir rates are ignored.

On an open chain, four extra transitions act on the ends: `alpha` injects a particle into an empty leftmost site and `gamma` removes one from an occupied leftmost site, while `delta` and `beta`, respectively, do the same at the rightmost site. The reservoir densities they impose are `alpha/(alpha+gamma)` on the left and `delta/(delta+beta)` on the right, so a density gradient — and with it a steady current — is set by making the two unequal.

It is important to stress our convention in the open case. A bond's rate depends on the two sites flanking it, but the outermost bonds of an open chain have a flank that falls off the lattice. The missing occupation is mirrored from the one that is present, which makes the flanking pair equal. This keeps both ends particle–hole symmetric (symmetry of the model under the bit flip transformation `1->0` and `0->1`) and is in any case a surface effect of order `1/L`.

## Output

Each run writes two files into `data/`, sharing a name built from the model (`NNN`, `KLS` or `WASEP`), the boundary (`O` for an open chain, `C` for a closed ring), the mechanism (`direct`, `heap` and `uniformized`), the rates, the sizes, the initial condition and the seed:

- a **CSV** with one row per recording time — the CGF estimate, the mean potential, the currents, the density profile, the effective sample size and the resampling counters;
- a **JSON** with the metadata: every parameter the run was launched with, the diagnostics at the end, and a *checkpoint*.

The column-by-column meaning of the CSV and the exact filename grammar are in [`Simulator/docs/OUTPUT_CONVENTIONS.md`](Simulator/docs/OUTPUT_CONVENTIONS.md).

### The checkpoint

The checkpoint is a snapshot of everything the engine needs to carry on from where it stopped, rather than of the results it has produced so far. Concretely, the C++ side serialises into a byte string: the current time and the next scheduled cloning and recording instants, the reference trajectory, all `M` walker configurations, each walker's log-weight and integrated current, the accumulated `log Z`, every counter, and the internal state of the random generator. Those bytes are binary, and JSON can only hold text, so they are [Base64](https://en.wikipedia.org/wiki/Base64)-encoded — a standard way of writing arbitrary bytes as ASCII characters, at a cost of about a third more space — and stored as one long string field. Restoring is the exact inverse: decode, hand the bytes back to the engine, and it is in precisely the state it was in.

Because the generator state is included, a restored run draws the *same* random numbers it would have drawn had it never stopped. That is what makes an extended run bitwise identical to the same run done in one go.

### The caching rules

Before starting, a run looks in `data/` for a stored run it can continue instead. A stored run qualifies only if it matches on **all** of:

- every parameter that affects the model — sizes, rates, fields, boundary, seed, initial condition, and even the event generator mechanism (since the three mechanisms have incompatible restart states);
- the engine version and the checkpoint format version, so that results from an older, differently-behaved engine are never silently mixed in;
- the presence of both the CSV and the JSON.

Only `tmax` is allowed to differ, and the stored run is then reused:

- if its `tmax` already reaches what you asked for, nothing is simulated and the stored rows are returned;
- otherwise its checkpoint is restored and only the missing time is simulated. 

## Analysis notebook

[`notebooks/plot_cgf_and_density.ipynb`](notebooks/plot_cgf_and_density.ipynb) loads every `WASEP_*`, `KLS_*` and `NNN_*` run in `data/` and plots the CGF estimates, the recorded density profiles, a colour plot of the density over site and time, and other observables, for one closed and one open selection. Numerical results are appropriately compared against existing theoretical predictions.



## Layout

```text
NNN/
  environment.yml      the conda environment everything below is run in
  simulations.toml     the plan: what to simulate
  run_NNN.sh           runs that plan
  data/                results, written here by every run
  notebooks/           analysis
  validation/          exact small-system references
  theory_algorithm_implementation.tex   the notes: physics, algorithm, code
  Simulator/           the simulator itself
    cpp/engine.cpp       the simulation
    Python_Cpp_Interface/  the Python package wrapping it
    setup_cpp.py         builds the extension
    run_NNN_cpp.py       one simulation, from the command line
    launch_simulations.py  expands simulations.toml and runs it in parallel
    build/               compiler output
    docs/, tests/
```

## Documentation

- [`theory_algorithm_implementation.tex`](theory_algorithm_implementation.tex) — **the notes**: the monitored chain, the large-deviation formalism, the algorithm step by step, and how it is implemented.
- [`LAUNCHER.md`](LAUNCHER.md) — the TOML plan format and parallel execution.
- [`Simulator/docs/SIMULATION_FLOW.md`](Simulator/docs/SIMULATION_FLOW.md) —
  flow charts of the run loop, of a single event, and of what the three mechanisms do differently.
- [`Simulator/docs/OUTPUT_CONVENTIONS.md`](Simulator/docs/OUTPUT_CONVENTIONS.md)
  — how to read what the engine writes: the output columns, and the naming,
  caching and restart rules for the files.
- [`Simulator/docs/CPU_REQUIREMENTS.md`](Simulator/docs/CPU_REQUIREMENTS.md) —
  build flags, and the run-time choice between the `PDEP` and portable bond-selection paths.

