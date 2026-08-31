#!/usr/bin/env python3

"""Expand an NNN/KLS/WASEP TOML plan and run jobs in parallel batches."""

import argparse
import itertools
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tomllib


OPTION_FLAGS = {
    'length': '--length',
    'walkers': '--walkers',
    'time': '--time',
    'field': '--field',
    'bias': '--bias',
    'measurement_strength': '--measurement-strength',
    'epsilon': '--epsilon',
    'delta_kls': '--delta-kls',
    'right_weights': '--right-weights',
    'left_weights': '--left-weights',
    'filling': '--filling',
    'alpha': '--alpha',
    'gamma': '--gamma',
    'delta': '--delta',
    'beta': '--beta',
    'boundary': '--boundary',
    'target_min': '--target-min',
    'target_max': '--target-max',
    'record_interval': '--record-interval',
    'seed': '--seed',
    'initial': '--initial',
    'progress': '--progress',
    'dynamics': '--dynamics',
    'output': '--output',
    'overwrite': '--overwrite',
}
VECTOR_OPTIONS = {'right_weights', 'left_weights'}
BOOLEAN_OPTIONS = {'overwrite'}
RATE_OPTIONS = {'epsilon', 'delta_kls', 'right_weights', 'left_weights'}
MODELS = {'WASEP', 'KLS', 'NNN'}
# The launcher and the runner it starts live in Simulator/; the plan, the data
# folder and everything else a run is described by live one level up, in the
# project folder.
SIMULATOR_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SIMULATOR_DIR.parent
RUNNER_DEFAULTS = {
    'length': 20,
    'walkers': 100,
    'time': 20.0,
    'field': 0.0,
    'bias': 0.0,
    'measurement_strength': 0.0,
    'epsilon': 0.6,
    'delta_kls': 0.0,
    'right_weights': None,
    'left_weights': None,
    'filling': 0.5,
    'alpha': 1.0,
    'gamma': 1.0,
    'delta': 1.0,
    'beta': 1.0,
    'boundary': 'open',
    'initial': 'alternating',
    'seed': 12345,
    'dynamics': 'direct',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run a TOML simulation plan using the available CPU cores.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('config', nargs='?', type=Path,
                        default=PROJECT_DIR / 'simulations.toml')
    parser.add_argument(
        '--workers', type=int,
        help='maximum simultaneous simulations; defaults to available CPUs',
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='check the plan and print every expanded command without running it',
    )
    return parser.parse_args()


def load_plan(path):
    with path.open('rb') as stream:
        return tomllib.load(stream)


def expand_plan(plan):
    unknown_sections = set(plan) - {'fixed', 'sweep', 'cases'}
    if unknown_sections:
        names = ', '.join(sorted(unknown_sections))
        raise ValueError(f'unknown plan section(s): {names}')

    fixed = plan.get('fixed', {})
    sweep = plan.get('sweep', {})
    cases = plan.get('cases', [])
    if not isinstance(fixed, dict) or not isinstance(sweep, dict):
        raise ValueError('[fixed] and [sweep] must be TOML tables')
    if not isinstance(cases, list) or not all(isinstance(case, dict) for case in cases):
        raise ValueError('cases must be written as one or more [[cases]] tables')

    overlap = set(fixed) & set(sweep)
    if overlap:
        raise ValueError(_overlap_message(overlap, '[fixed]', '[sweep]'))
    for number, case in enumerate(cases, start=1):
        overlap = set(case) & (set(fixed) | set(sweep))
        if overlap:
            raise ValueError(_overlap_message(
                overlap, f'[[cases]] entry {number}', '[fixed] or [sweep]'))

    sweep_names = list(sweep)
    sweep_values = []
    for name in sweep_names:
        values = sweep[name]
        if not isinstance(values, list) or not values:
            raise ValueError(f'[sweep].{name} must be a non-empty array')
        sweep_values.append(values)

    combinations = itertools.product(*sweep_values) if sweep_names else [()]
    case_values = cases or [{}]
    jobs = []
    for values in combinations:
        repeated = dict(zip(sweep_names, values))
        for case in case_values:
            jobs.append(fixed | repeated | case)
    return jobs


def _overlap_message(names, first, second):
    joined = ', '.join(sorted(names))
    return (f'each parameter belongs to one section; {joined} occurs in '
            f'{first} and {second}')


def validate_job(job):
    unknown = set(job) - set(OPTION_FLAGS) - {'model'}
    if unknown:
        raise ValueError(f"unknown parameter(s): {', '.join(sorted(unknown))}")

    model = str(job.get('model', 'KLS')).upper()
    if model not in MODELS:
        allowed = ', '.join(sorted(MODELS))
        raise ValueError(f'NNN plans support model = {allowed}; got {model}')

    present_rates = RATE_OPTIONS & set(job)
    if model == 'WASEP' and present_rates:
        raise ValueError('WASEP jobs must not specify KLS/NNN rate parameters')
    if model == 'KLS' and present_rates & VECTOR_OPTIONS:
        raise ValueError('KLS jobs use epsilon/delta_kls, not explicit weights')
    if model == 'NNN':
        if present_rates & {'epsilon', 'delta_kls'}:
            raise ValueError('NNN jobs use explicit weights, not epsilon/delta_kls')
        if not VECTOR_OPTIONS <= set(job):
            raise ValueError('NNN jobs require right_weights and left_weights')
        for name in VECTOR_OPTIONS:
            values = job[name]
            if not isinstance(values, list) or len(values) != 4:
                raise ValueError(f'{name} must contain four rates')
    for name in BOOLEAN_OPTIONS:
        if name in job and not isinstance(job[name], bool):
            raise ValueError(f'{name} must be true or false')
    return model


def command_for_job(job, simulator_dir=SIMULATOR_DIR,
                    python_executable=sys.executable):
    model = validate_job(job)
    # One engine serves all three models: WASEP selects its single-mask path
    # rather than a separate program.
    runner = simulator_dir / 'run_NNN_cpp.py'
    command = [str(python_executable), str(runner), '--model', model]
    for name, flag in OPTION_FLAGS.items():
        if name not in job:
            continue
        value = job[name]
        if name in BOOLEAN_OPTIONS:
            if value:
                command.append(flag)
        elif name in VECTOR_OPTIONS:
            command.append(flag)
            command.extend(str(item) for item in value)
        else:
            command.extend((flag, str(value)))
    return command


def output_identity(job, project_dir=PROJECT_DIR):
    if 'output' in job:
        output = Path(job['output'])
        return ('output', str(output if output.is_absolute()
                              else project_dir / output))
    values = RUNNER_DEFAULTS | job
    model = str(values.get('model', 'KLS')).upper()
    boundary = values['boundary']
    identity = [
        model, values['dynamics'], values['length'], values['walkers'],
        values['time'], values['field'], values['bias'],
        values['measurement_strength'], boundary, values['seed'],
    ]
    if boundary == 'periodic':
        identity.append(values['filling'])
    else:
        # ``initial`` only reaches the filename on an open chain, because a
        # ring ignores it and fixes its particle number from the filling.
        identity.extend(values[name] for name in ('alpha', 'gamma', 'delta', 'beta'))
        identity.append(values['initial'])
    if model == 'KLS':
        identity.extend((values['epsilon'], values['delta_kls']))
    elif model == 'NNN':
        identity.extend((tuple(values['right_weights']), tuple(values['left_weights'])))
    return tuple(identity)


def reject_output_collisions(jobs, project_dir=PROJECT_DIR):
    seen = {}
    for index, job in enumerate(jobs, start=1):
        validate_job(job)
        identity = output_identity(job, project_dir)
        if identity in seen:
            raise ValueError(
                f'jobs {seen[identity]} and {index} would write the same output; '
                'vary a filename-defining parameter or seed'
            )
        seen[identity] = index


def available_cpu_count():
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def run_batches(commands, project_dir, workers):
    failures = []
    total = len(commands)
    for batch_start in range(0, total, workers):
        batch = commands[batch_start:batch_start + workers]
        processes = []
        try:
            for offset, command in enumerate(batch):
                number = batch_start + offset + 1
                print(f'[{number}/{total}] starting: {shlex.join(command)}', flush=True)
                processes.append((number, subprocess.Popen(command, cwd=project_dir)))
            for number, process in processes:
                returncode = process.wait()
                if returncode == 0:
                    print(f'[{number}/{total}] finished', flush=True)
                else:
                    failures.append((number, returncode))
                    print(f'[{number}/{total}] failed with exit code {returncode}',
                          file=sys.stderr, flush=True)
        except KeyboardInterrupt:
            print('\nStopping active simulations...', file=sys.stderr, flush=True)
            for _, process in processes:
                if process.poll() is None:
                    process.terminate()
            for _, process in processes:
                process.wait()
            raise
    return failures


def main():
    args = parse_args()
    config = args.config.resolve()
    project_dir = PROJECT_DIR
    jobs = expand_plan(load_plan(config))
    if not jobs:
        raise ValueError('the plan expands to no simulations')
    reject_output_collisions(jobs, project_dir)
    commands = [command_for_job(job) for job in jobs]
    missing_runners = sorted({command[1] for command in commands if not Path(command[1]).is_file()})
    if missing_runners:
        raise ValueError(f"runner not found: {', '.join(missing_runners)}")

    available = available_cpu_count()
    workers = args.workers if args.workers is not None else available
    if workers < 1:
        raise ValueError('--workers must be positive')
    workers = min(workers, available, len(commands))
    print(f'Expanded {len(commands)} simulation(s); running up to {workers} '
          f'at once on {available} available CPU(s).')
    if args.dry_run:
        for number, command in enumerate(commands, start=1):
            print(f'[{number}/{len(commands)}] {shlex.join(command)}')
        return

    failures = run_batches(commands, project_dir, workers)
    if failures:
        summary = ', '.join(f'{number} (exit {code})'
                            for number, code in failures)
        raise SystemExit(f'{len(failures)} simulation(s) failed: {summary}')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f'Launcher error: {error}') from None
