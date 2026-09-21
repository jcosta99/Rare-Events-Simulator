# How a run flows

Four views of the same computation, from the outside in. The diagrams render
directly on GitHub, so the labels are kept short and anything that needs a
sentence is written as one underneath.

Two companion documents pick up where this one stops:
[`OUTPUT_CONVENTIONS.md`](OUTPUT_CONVENTIONS.md) for what a run writes, how the
files are named, and how caching and restarts work;
[`CPU_REQUIREMENTS.md`](CPU_REQUIREMENTS.md) for the bit-selection step that
dominates the cost of an event.

## Vocabulary

- **Walker** — one trajectory in the population, biased towards rare currents.
- **Reference** — a single untilted trajectory, used for the monitoring overlap.
- **Potential** `V_i` — the instantaneous rate at which walker `i`'s weight
  grows or shrinks.
- **Stop** — an instant at which the whole population is brought to the same
  time: the next cloning check, the next recording, or `tmax`.
- **Block** — everything that happens between two consecutive stops.
- **Resampling** — replacing the weighted population by `M` equally weighted
  descendants, preserving the influence of the old weights in the normalization.
- **Certificate** — a running upper bound on how far the log-weights can have
  spread. Only a mechanism that keeps one global clock can watch it.

## 1. From a plan to a figure

```mermaid
flowchart LR
    toml["simulations.toml<br/>the plan"] --> launcher["launch_simulations.py<br/>expands it"]
    launcher -->|"one process per CPU"| cli["run_NNN_cpp.py<br/>one run"]
    cli --> wrapper["NNN_class.py<br/>validate, name, cache"]
    wrapper --> engine["engine.cpp<br/>the simulation"]
    engine --> store[("data/<br/>CSV rows<br/>JSON metadata")]
    store --> nb["notebooks/<br/>CGF, profiles, ESS"]
```

The wrapper is
[`Simulator/Python_Cpp_Interface/NNN_class.py`](../Python_Cpp_Interface/NNN_class.py)
and the engine [`Simulator/cpp/engine.cpp`](../cpp/engine.cpp). Before starting
the engine, the wrapper looks in `data/` for a result whose parameters match: an
equal one is returned as it stands, and a shorter one is resumed from its
checkpoint so that only the missing time is simulated. What counts as a match is
set out in [`OUTPUT_CONVENTIONS.md`](OUTPUT_CONVENTIONS.md).

## 2. The run loop, in all three mechanisms at once

One iteration per block. All three mechanisms plan the same stops and, once the
block is over, record and resample identically. They differ only in how the
population is carried to the stop, and in how the block is then judged.

```mermaid
flowchart TD
    start([run]) --> init["write the first row"]
    init --> plan["plan the stop"]

    subgraph advance["advance to that stop"]
        heap["heap<br/>pop the earliest event,<br/>execute, redraw, sift"]
        uni["uniformized<br/>propose at the top rate,<br/>accept or reject"]
        dir["direct<br/>the reference, then each<br/>walker replaying it"]
    end

    plan --> heap
    plan --> uni
    plan --> dir

    heap --> cert
    uni --> cert
    cert{"certificate<br/>tripped?"}
    cert -->|yes| early["examine the weights"]
    early --> plan
    cert -->|no| kind

    dir --> rev{"overshot<br/>target_max?"}
    rev -->|yes| back["restore, halve tau,<br/>replay"]
    back --> plan
    rev -->|no| grow["keep it, predict<br/>the next tau"]
    grow --> kind

    kind{"which stop?"}
    kind -->|recording| force["resample"]
    force --> row["write a row"]
    kind -->|cloning check| maybe["examine the weights"]

    row --> fin{"tmax reached?"}
    maybe --> fin
    fin -->|no| plan
    fin -->|yes| out([return the rows])
```

The stop is the earliest of the next cloning check, the next recording, and
`tmax`. `heap` and `uniformized` keep one global clock, so the certificate can
be watched continuously and the block is cut the instant the bound is reached;
`direct` advances each walker separately, so it has no running certificate and
predicts `tau` instead, judging the block only once it is over.

Examining the weights is not the same as resampling. The certificate is an
upper bound, so the exact range may still be within `target_min`, in which case
nothing is resampled. A recording is the one stop that resamples on its own
account — unless the run is untilted at `s = k = 0`, where every weight stays
exactly one and resampling would only copy the population.

`tau` shrinks fast and grows slowly: a block that overshoots `target_max` is
discarded and replayed at half the step, while a block that holds may at most
double it, so one lucky block cannot stretch the next one out of range.

## 3. One accepted event

The path every walker event takes, in the rejection-free mechanisms.

```mermaid
flowchart TD
    draw["draw a uniform over<br/>the walker's total rate"] --> dir{"which part<br/>of the rate?"}
    dir -->|"right hops"| cls["find the class"]
    dir -->|"left hops"| cls
    dir -->|"reservoir"| bnd["inject or remove<br/>at an end site"]
    cls --> ordinal["the remainder indexes<br/>that class's allowed bonds"]
    ordinal --> pick["select_set_bit:<br/>pick that bond"]
    pick --> apply
    bnd --> apply
    apply["apply_walker_event:<br/>close the weight integral,<br/>flip the sites, add the current"]
    apply --> masks["update the masks"]
    masks --> rate["recompute the activity<br/>and the potential"]
    rate --> wait["draw the next<br/>waiting time"]
```

A class's allowed bonds are held as set bits in a mask word, so choosing the
`n`-th of them is what selects the bond. `select_set_bit` does that either with
a single `PDEP` instruction or with a portable walk over the bits; which route
is taken is decided once, when the extension module loads, and why is the
subject of [`CPU_REQUIREMENTS.md`](CPU_REQUIREMENTS.md).

Updating the masks is local for `KLS` and `NNN` — a hop can only change the five
bonds around it — whereas `WASEP` rebuilds them for the whole lattice a word at
a time.

## 4. Why three mechanisms

All three implement the same process and are checked against each other, so the
choice is one of cost rather than of physics.

| mechanism | how the next event is found | rejections | cost | why it is kept |
| --- | --- | --- | --- | --- |
| `direct` | each walker runs alone between stops, replaying the reference's stored events | none | fastest; the default | it is the most efficient |
| `heap` | every process holds one pending time; a binary heap keeps the earliest at the front | none | 1.2–2.8× `direct` | true chronological order, and an idea that may be useful elsewhere |
| `uniformized` | one clock at the heaviest possible rate proposes any jump between adjacent sites | rejected down to the true rate | 2.7–5.1× `direct` | simple enough to be obviously correct, so it benchmarks the other two |

Only `uniformized` reports a non-zero `RejectedWalkerEvents`; the other two
never propose an event they do not take. `heap` and `uniformized` reproduce the
earlier single-model engine bitwise, while `direct` consumes the random stream
in its own order and so agrees in distribution rather than per seed. The full
account, with the measurements behind the cost column, is in
[`theory_algorithm_implementation.tex`](../../theory_algorithm_implementation.tex),
sections *The algorithm* and *Implementation*.

## Where to look in the code

| stage | function |
| --- | --- |
| the loop above | `Simulation::run` |
| the stop schedule | `Simulation::schedule_next_cloning_check` |
| advancing a block | `HeapDynamics::advance_weighted`, `UniformizedDynamics::advance_weighted`, `DirectDynamics::advance_weighted` |
| the block review | `DirectDynamics::review_block` |
| event selection | `Simulation::execute_walker_event`, `select_within_direction`, `select_set_bit` |
| applying an event | `Simulation::apply_walker_event` |
| mask maintenance | `update_masks_near`, `rebuild_masks_uniform` |
| the certificate | `accumulate_weight_bound`, `weight_bound_reached` |
| resampling and the estimator | `Simulation::clone_population` |
| the recorded row | `Simulation::output_record` |
