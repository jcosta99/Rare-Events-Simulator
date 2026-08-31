# NNN/KLS simulation launcher

**Activate the environment first.** `run_NNN.sh` calls plain `python3` and does not resolve an environment itself, so whatever is active when you launch it is what runs. The environment supplies the libraries the scripts import, and the Python that the C++ engine is built against, so you should compile it and run it from the same environment.

Create the environment once with:
```bash
conda env create -f environment.yml
```

Then, to activate it,
```bash
conda activate monitored-exclusion     
```

Finally, edit `simulations.toml` and run:
```bash
bash run_NNN.sh
```

Note that if activating is inconvenient — inside a script, or from an editor — point `PYTHON_BIN` at the interpreter instead, which has the same effect without changing the active environment:
```bash
PYTHON_BIN=$(conda run -n monitored-exclusion which python) bash run_NNN.sh
```

The shell file `run_NNN.sh` is only a small entry point; `Simulator/launch_simulations.py` expands the plan (which can contain instructions for multiple runs) and starts up to one single-threaded simulation per available CPU. Use these controls when checking a plan or limiting memory:

```bash
bash run_NNN.sh --dry-run      # expand and validate the plan, print every command, run nothing
bash run_NNN.sh --workers 4    # run at most 4 simulations at once instead of one per CPU
```

`--dry-run` runs every check on the plan and prints the full command of each job without starting any of them, so it tells you the plan is sound before you commit the CPUs. It is also the way to see what you are actually submitting: how many simulations the plan expands to, and what every parameter of each one is. Parameters you forgot to set appear there too, filled in with their defaults, which makes them easy to spot.
`--workers` is about memory rather than speed — it sets the maximum amount of processes running at the same time. Each process holds its own population of walkers, so a plan with large `walkers` or `length` can exhaust RAM well before it runs out of cores.

Set `SIMULATION_CONFIG=/path/to/another.toml` to run a plan kept elsewhere.

## The three plan sections

A plan has three sections, and together they describe the list of jobs to run.

- `[fixed]` — the parameters that are the same for every run. Write one here once instead of repeating it in every job.

- `[sweep]` — the parameters you want to scan, each given as an array. This is the right section when the parameters vary independently and you 
do want all the pairings.

- `[[cases]]` — for when a sweep is the wrong shape, because two or more parameters have to change *together* and only some of their combinations are meaningful (crossing them in a sweep would produce the unwanted combinations too). Each `[[cases]]` table is one wanted combination, written out as a unit, and the plan runs one job per `[[cases]]` table.

Every case is combined with the full sweep, so the number of simulations is the number of cases times the product of the lengths of all the sweep arrays. 
A parameter belongs to exactly one section; naming it in two is an error, not an override.

The next example plan runs each bias and system size for both a KLS model and an explicitly weighted NNN model, giving eight simulations:

```toml
[fixed]
walkers = 2000
time = 20000.0
field = 5.0
measurement_strength = 0.0
boundary = "periodic"
filling = 0.3
record_interval = 10.0
dynamics = "direct"
seed = 12345

[sweep]
length = [100, 200]
bias = [-8.0, -5.0]

[[cases]]
model = "KLS"
epsilon = 0.6
delta_kls = 0.0


[[cases]]
model = "NNN"
right_weights = [1.0, 0.4, 1.6, 1.0]
left_weights = [1.0, 1.6, 0.4, 1.0]
```

This submits, in this order:

| job | `length` | `bias` | `model` | rate parameters                                                               |
| --- | -------- | ------ | ------- | ----------------------------------------------------------------------------- |
| 1   | `100`    | `-8.0` | `"KLS"` | `epsilon = 0.6`, `delta_kls = 0.0`                                            |
| 2   | `100`    | `-8.0` | `"NNN"` | `right_weights = [1.0, 0.4, 1.6, 1.0]`, `left_weights = [1.0, 1.6, 0.4, 1.0]` |
| 3   | `100`    | `-5.0` | `"KLS"` | `epsilon = 0.6`, `delta_kls = 0.0`                                            |
| 4   | `100`    | `-5.0` | `"NNN"` | `right_weights = [1.0, 0.4, 1.6, 1.0]`, `left_weights = [1.0, 1.6, 0.4, 1.0]` |
| 5   | `200`    | `-8.0` | `"KLS"` | `epsilon = 0.6`, `delta_kls = 0.0`                                            |
| 6   | `200`    | `-8.0` | `"NNN"` | `right_weights = [1.0, 0.4, 1.6, 1.0]`, `left_weights = [1.0, 1.6, 0.4, 1.0]` |
| 7   | `200`    | `-5.0` | `"KLS"` | `epsilon = 0.6`, `delta_kls = 0.0`                                            |                    
| 8   | `200`    | `-5.0` | `"NNN"` | `right_weights = [1.0, 0.4, 1.6, 1.0]`, `left_weights = [1.0, 1.6, 0.4, 1.0]` |

with every other parameter taken from `[fixed]`. The two models could not have been written as a sweep: `epsilon` belongs only to `KLS` and the weight vectors only to `NNN`, and pairing either with the wrong model is rejected.

Both `[sweep]` and `[[cases]]` are optional. If both are absent and a plan consists of `[fixed]` alone, only one simulation is launched.

Two jobs that would write the same file are rejected before any process starts, and the message names both job numbers. Only the parameters that appear in the output filename are compared, so two jobs differing solely in `cloning_interval`, `record_interval`, `target_min` or `target_max` count as the same file and are rejected; separate them with the seed, or with any other parameter that shows up in the name. 

Leaving a parameter out of the plan entirely is not an error: it takes the runner default, which `--dry-run` will show you.

## Common parameters

| Name                                    | Meaning                                                                                         |
| --------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `model`      *                          | Choose between the general option `NNN` and the specific models: `WASEP` and `KLS`              |
| `length`, `walkers`, `time`             | Lattice size `L`, population size `M`, and final time `t_{max}`.                                |
| `field`, `bias`, `measurement_strength` | External field `E`, current-counting field `s`, and monitoring coupling `k`.                    |
| `boundary`                              | Open chain: `"open"`; closed chain: `"periodic"`.                                               |
| `filling`                               | Periodic: average density (conserved quantity); ignored for an open chain.                      |
| `alpha`, `gamma`                        | Left injection and removal rates for an open chain, respectively.                               |
| `delta`, `beta`                         | Right injection and removal rates for an open chain, respectively.                              |
| `seed`                                  | Random seed; sweep it for independent repetitions.                                              |
| `record_interval`                       | Time between records (recording always forces resampling).                                      |
| `target_min`                            | Resample once the heaviest walker weight exceeds the lightest by this factor.                   |
| `target_max` **                         | The heaviest walker weight cannot exceed the lightest by `target_max`. Force resampling before. |
| `cloning_interval`                      | Force a maximum cloning interval for resampling to occur; normally omit it.                     |
| `initial`    ***                        | `"empty"` (no particles), `"alternating"`, `"random"`, or `"full"` (every site occupied).       |
| `progress`                              | Print progress every this many records.                                                         |
| `dynamics`   ****                       | `"direct"` (fastest), `"heap"` (nicest), or `"uniformized"` for validation.                     |
| `overwrite`                             | Default `false`. Replaces a file whose *name* collides but whose parameters differ; see Output. |
| `output`                                | Optional explicit CSV path; normally omit it.                                                   |

*   Three possible scenarios:
    -`"NNN"`: Receives as input two vectors with 4 entries: `right_weights` and `left_weights`. For a bond `(b,b+1)`, NNN rates are indexed by the flanking occupations `(n[b-1], n[b+2])` in the order `(0,0), (0,1), (1,0), (1,1)`. 
    -`"KLS"`: Receives two parameters as input, `epsilon` and `delta_KLS`, which completely specify the model's rates (this is done internally).
    Domain of the parameters: `|epsilon| < 1` and `|delta_kls| < 1`.
    -`WASEP`: No extra parameter is received as input. The model rates are uniform and equal to 1.0 .

    Providing the wrong parameters to the wrong model (or forgetting the necessary input) raises an error.
    Note that all these are bare interaction rates: they do not include the external or counting field in them. The engine then applies `exp((E+s)/L)` to right hops and `exp(-(E+s)/L)` to left hops.
    See README.md or theory_algorithm_implementation.tex for more details.

    **Default**: `KLS` with `epsilon = 0.6` and `delta_KLS=0`.

**  `Heap` and `uniformized` examine the population as soon as their running bound reaches `target_max`; `direct` replays any block that
    overshot `target_max` with a halved step.

*** The initial condition flag is only used for open boundary conditions. It is ignored on a ring, where `filling` fixes the particle number.

****See README.md or theory_algorithm_implementation.tex for more details about each implementation.

## Output

All simulations launched here write under the folder `data/` — sitting at the current folder, beside `simulations.toml` — so the NNN notebook can read them together. Large NNN populations carry eight transition mask families per walker; cap `--workers` if RAM, rather than CPU, is limiting.

Result names carry the parameters that identify a run: the model family and the geometry (`O` open, `C` closed), the event mechanism, the rate description, then `length`, `walkers`, `time`, `field`, `bias`, `measurement_strength`, the boundary rates or the filling, the initial condition, and the seed. Several parameters that genuinely affect the result are deliberately left out of the name — `cloning_interval`, `record_interval`, `target_min`, `target_max`, and the engine version — so a job that differs from an existing file only in those will try to store a file with the same name. That collision is the one situation `overwrite` exists for: if true, the stored result is replaced, and its metadata then carries the new signature; if false, an error is raised and no data is saved. The name is only checked when the finished run goes to write, so a job rejected this way has already spent its full runtime. The full grammar is in [`Simulator/docs/OUTPUT_CONVENTIONS.md`](Simulator/docs/OUTPUT_CONVENTIONS.md).

Note: the result names carry the event mechanism — `..._Heap_...`, `_Direct_`, `_Uniformized_` — because a stored run may only be extended by the mechanism that produced it: their restart states are not interchangeable. Two jobs that differ only in `dynamics` are therefore separate results, not a collision.