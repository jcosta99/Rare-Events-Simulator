# NNN — next-nearest-neighbour exclusion engine

A Monte Carlo simulator for rare-current statistics in monitored stochastic lattice gases. A C++ engine evolves populations of weighted trajectories; Python provides the command-line interface, parameter sweeps, caching, and analysis workflow.

The engine supports the weakly asymmetric exclusion process (WASEP), the Katz–Lebowitz–Spohn model (KLS), and a general next-nearest-neighbour (NNN) model whose hopping rates depend on the occupations flanking each bond. Runs may use an open chain coupled to reservoirs or a periodic ring.

For the physical model, monitoring protocol, large-deviation formalism, and derivation of the estimator, see [`theory_algorithm_implementation.pdf`](Simulator/docs/theory_algorithm_implementation.pdf).

## Main features

- Bit-parallel configurations and transition masks.
- Population dynamics based on continuous-time cloning.
- Three interchangeable event generation mechanisms: `direct`, `heap`, and
  `uniformized`.
- TOML plans for parameter sweeps and correlated parameter sets.
- Checkpointed results that can be extended without restarting a run.
- Exact and asymptotic reference calculations for validation.

## Quick start

Create the environment and build the extension:

```bash
conda env create -f environment.yml
conda activate monitored-exclusion
cd Simulator
python setup_cpp.py build_ext --inplace
cd ..
```

The extension must be built from `Simulator/`; the remaining commands may be
run from the repository root.

Edit [`simulations.toml`](simulations.toml), inspect the expanded plan, and run
it:

```bash
bash run_NNN.sh --dry-run
bash run_NNN.sh
```

The launcher starts at most one single-threaded simulation per available CPU. Use `bash run_NNN.sh --workers N` to cap the number of simulations running in parallel, when memory rather than cores is the limit. The plan format and all launcher options are described in [`LAUNCHER.md`](Simulator/docs/LAUNCHER.md).

A single simulation can also be launched directly:

```bash
python Simulator/run_NNN_cpp.py \
  -L 100 -M 2000 -t 200 -E 7 -s -3.35 \
  --epsilon 0.6 --boundary open \
  --alpha 1 --gamma 1 --delta 1 --beta 1
```

The standard build targets a CPU supporting `x86-64-v2`. For older x86 processors, non-x86 machines, or details of the optional BMI2 optimization, see [`Simulator/docs/CPU_REQUIREMENTS.md`](Simulator/docs/CPU_REQUIREMENTS.md).

## Choosing a model

Set `model` in [`simulations.toml`](simulations.toml):

| value | process | additional rate parameters |
| --- | --- | --- |
| `WASEP` | weakly asymmetric exclusion process | none |
| `KLS` | Katz–Lebowitz–Spohn model | `epsilon`, `delta_kls` |
| `NNN` | general flanking-dependent model | `right_weights`, `left_weights` |

Set `boundary = "periodic"` for a ring, where `filling` fixes the conserved
particle number. Set `boundary = "open"` for a chain coupled to reservoirs,
with rates `alpha`, `gamma`, `delta`, and `beta`. Parameter meanings, defaults,
and the ordering of the NNN weights are kept in
[`LAUNCHER.md`](Simulator/docs/LAUNCHER.md).

## Choosing an event mechanism

All three values of `dynamics` simulate the same process:

- `direct` evolves walkers independently between population stops. It is the
  fastest mechanism and the default.
- `heap` maintains all pending events on a shared clock using a binary heap.
- `uniformized` generates simple maximal-rate proposals and rejects them down
  to the true rate. Its independent construction makes it useful for
  validation.

The mechanisms differ in event scheduling and in how they control the interval
between population examinations. Their flow chart is illustrated in
[`Simulator/docs/SIMULATION_FLOW.md`](Simulator/docs/SIMULATION_FLOW.md); the
full algorithm and performance comparison are in the theory notes [`theory_algorithm_implementation.pdf`](Simulator/docs/theory_algorithm_implementation.pdf).

## Results and restarts

Runs write a CSV time series and a JSON metadata file under `data/`. The JSON also contains the checkpoint required to resume the simulation. If a compatible run already exists, the wrapper restores it and simulates only the missing time (if needed).

The variables and respective definitions of the CSV columns, filenames, cache matching, and restart behavior are in [`Simulator/docs/OUTPUT_CONVENTIONS.md`](Simulator/docs/OUTPUT_CONVENTIONS.md).

[`notebooks/plot_cgf_and_density.ipynb`](notebooks/plot_cgf_and_density.ipynb) loads the stored runs and plots the cumulant generating function estimates, density profiles, space-time density, currents and population diagnostics.

## Validation

The three event mechanisms are cross-checked against one another. The cumulant generating function of small open systems is also compared with that extracted from exact diagonalization of the tilted generator, while for larger SSEP systems it is compared with the known asymptotic result. The reference data and regeneration instructions are in [`validation/README.md`](validation/README.md); the numerical comparisons are discussed in the theory notes [`theory_algorithm_implementation.pdf`](Simulator/docs/theory_algorithm_implementation.pdf).

## Repository map

| Path | What it is |
| --- | --- |
| [`simulations.toml`](simulations.toml) | The simulation plan |
| [`run_NNN.sh`](run_NNN.sh) | Runs that plan |
| [`data/`](data) | Generated results |
| [`notebooks/`](notebooks) | Analysis notebooks |
| [`validation/`](validation) | Exact and asymptotic checks |
| [`theory_algorithm_implementation.pdf`](Simulator/docs/theory_algorithm_implementation.pdf) | Physics, algorithm and implementation |
| [`Simulator/cpp/engine.cpp`](Simulator/cpp/engine.cpp) | The C++ simulation engine |
| [`Simulator/Python_Cpp_Interface/`](Simulator/Python_Cpp_Interface) | Python wrapper package |
| [`Simulator/run_NNN_cpp.py`](Simulator/run_NNN_cpp.py) | A single run, from the command line |
| [`Simulator/launch_simulations.py`](Simulator/launch_simulations.py) | Expansion of the TOML plan and parallel launcher |
| [`Simulator/setup_cpp.py`](Simulator/setup_cpp.py) | Extension build script |
| [`Simulator/docs/`](Simulator/docs) | Focused implementation documentation |

## Documentation guide

| Question | Source of truth |
| --- | --- |
| What is being computed, and why does the estimators work? | [`theory_algorithm_implementation.pdf`](Simulator/docs/theory_algorithm_implementation.pdf) |
| How do I configure and launch a collection of runs? | [`LAUNCHER.md`](Simulator/docs/LAUNCHER.md) |
| How does a run move through the code? | [`Simulator/docs/SIMULATION_FLOW.md`](Simulator/docs/SIMULATION_FLOW.md) |
| What do the output files contain, and when can they be reused? | [`Simulator/docs/OUTPUT_CONVENTIONS.md`](Simulator/docs/OUTPUT_CONVENTIONS.md) |
| Which processors are supported? | [`Simulator/docs/CPU_REQUIREMENTS.md`](Simulator/docs/CPU_REQUIREMENTS.md) |
| How are the numerical checks generated? | [`validation/README.md`](validation/README.md) |

The bibliography for the theory notes is in [`algorithm_refs.bib`](Simulator/docs/algorithm_refs.bib).

## Citation

If this simulator supported published work, a citation is welcome — GitHub's Cite this repository button gives a ready BibTeX entry.
