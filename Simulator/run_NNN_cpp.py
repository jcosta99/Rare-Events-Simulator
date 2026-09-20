#!/usr/bin/env python3

"""Command-line interface for the C++ rejection-free simulation engine."""

import argparse
from pathlib import Path

from Python_Cpp_Interface import NNN


# Results belong to the project, not to the simulator that produced them, so
# they are written next to simulations.toml and the notebooks rather than
# inside Simulator/ -- and to the same place whatever the working directory.
DATA_DIR = Path(__file__).resolve().parents[1] / 'data'


def parse_args():
    parser = argparse.ArgumentParser(
        description='C++ cloning simulation of monitored open or periodic NNN',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-L', '--length', type=int, default=20)
    parser.add_argument('-M', '--walkers', type=int, default=100)
    parser.add_argument('-t', '--time', type=float, default=20.0)
    parser.add_argument('-E', '--field', type=float, default=0.0)
    parser.add_argument('-s', '--bias', type=float, default=0.0)
    parser.add_argument('-k', '--measurement-strength', type=float, default=0.0)
    parser.add_argument(
        '--model', choices=['WASEP', 'KLS', 'NNN'],
        help=('rate family: WASEP keeps one rate per direction and uses the '
              'cheaper single-mask bookkeeping, KLS builds the eight class '
              'rates from --epsilon and --delta-kls, NNN takes them as given '
              'weights; omitted, it is read off the rates'),
    )
    # No default here: argparse would apply one whether or not the flag was
    # given, so an unset interaction would be indistinguishable from 0.6.  The
    # constructor resolves None per model.
    parser.add_argument(
        '--epsilon', type=float,
        help=('KLS next-nearest-neighbour coupling, used by --model KLS and '
              'ignored by the others; unset it is 0.6'),
    )
    parser.add_argument(
        '--delta-kls', type=float,
        help=('KLS asymmetry parameter, used by --model KLS and ignored by '
              'the others; unset it is 0, which keeps particle-hole symmetry'),
    )
    parser.add_argument(
        '--right-weights', type=float, nargs=4, metavar=('W00', 'W01', 'W10', 'W11'),
        help=('explicit right-hop rates for flanking occupations '
              '(n_b-1, n_b+2) = (0,0) (0,1) (1,0) (1,1); used by --model NNN, '
              'which requires it together with --left-weights'),
    )
    parser.add_argument(
        '--left-weights', type=float, nargs=4, metavar=('W00', 'W01', 'W10', 'W11'),
        help='explicit left-hop rates, same ordering',
    )
    parser.add_argument(
        '--filling', type=float, default=0.5,
        help=(
            'target closed-ring density; ignored for open boundaries and '
            'rounded to the nearest integer particle number'
        ),
    )
    parser.add_argument(
        '--alpha', type=float, default=1.0,
        help='open case: inject a particle into an empty leftmost site',
    )
    parser.add_argument(
        '--gamma', type=float, default=1.0,
        help='open case: remove a particle from an occupied leftmost site',
    )
    parser.add_argument(
        '--delta', type=float, default=1.0,
        help='open case: inject a particle into an empty rightmost site',
    )
    parser.add_argument(
        '--beta', type=float, default=1.0,
        help='open case: remove a particle from an occupied rightmost site',
    )
    parser.add_argument(
        '--boundary', choices=['open', 'periodic'], default='open',
        help='open reservoirs or a closed periodic ring',
    )
    parser.add_argument(
        '--target-min', type=float,
        help=('resample once the heaviest walker weight exceeds the lightest '
              'one by this factor'),
    )
    parser.add_argument(
        '--target-max', type=float,
        help=('never let that ratio pass this factor: the single-clock '
              'mechanisms examine the population as soon as their bound '
              'reaches it, and the direct one replays any block that '
              'overshot it with a shorter step'),
    )
    parser.add_argument('--cloning-interval', type=float)
    parser.add_argument('--record-interval', type=float)
    parser.add_argument('--seed', type=int, default=12345)
    parser.add_argument(
        '--initial', choices=['empty', 'alternating', 'random', 'full'],
        default='alternating',
        help=('open-chain initial configuration; ignored on a ring, where the '
              'particle number is fixed by --filling'),
    )
    parser.add_argument('--progress', type=int, metavar='N')
    parser.add_argument(
        '--dynamics', choices=['direct', 'heap', 'uniformized'],
        default='direct',
        help=('per-walker evolution between fixed cloning instants, '
              'accepted-event heap on one clock, or blind '
              'rejected-proposal baseline'),
    )
    parser.add_argument('-o', '--output', type=Path)
    parser.add_argument(
        '--overwrite', action='store_true',
        help=('replace a stored result of the same name that was produced with '
              'different parameters; running parameters that match a stored '
              'result returns it unchanged and writes nothing either way'),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    parameters = dict(
        L=args.length, M=args.walkers, tmax=args.time,
        E=args.field, s=args.bias, k=args.measurement_strength,
        model=args.model, epsilon=args.epsilon, delta_kls=args.delta_kls,
        right_weights=args.right_weights, left_weights=args.left_weights,
        filling=args.filling,
        alpha=args.alpha, gamma=args.gamma,
        delta=args.delta, beta=args.beta,
        target_min=args.target_min, target_max=args.target_max,
        cloning_interval=args.cloning_interval,
        record_interval=args.record_interval,
        seed=args.seed, initial=args.initial, progress=args.progress,
        dynamics=args.dynamics,
        boundary=args.boundary,
    )
    requested = NNN(**parameters)
    preferred = args.output
    if preferred is None:
        preferred, _ = requested.result_paths(DATA_DIR)
    simulation, rows, csv_path, json_path = NNN.run_cached(
        data_dir=preferred.parent,
        preferred_output=preferred,
        overwrite=args.overwrite,
        **parameters,
    )
    final = rows[-1]
    diagnostics = simulation.diagnostics()
    print(f'Finished C++ {args.dynamics} simulation')
    print(f'Output: {csv_path}')
    print(f'Metadata: {json_path}')
    print(f"Cache action: {simulation.cache_action}")
    print(f"Rate model: {diagnostics['rate_model']}  "
          f"right {diagnostics['right_weights']}  "
          f"left {diagnostics['left_weights']}")
    print(f"Resampling: {diagnostics['resampling_steps']} of "
          f"{diagnostics['resampling_checks']} examinations "
          f"(targets {diagnostics['target_min']:g} to "
          f"{diagnostics['target_max']:g})")
    print(f"Accepted walker events: {diagnostics['accepted_walker_events']}")
    print(f"Rejected walker events: {diagnostics['rejected_walker_events']}")
    print(f"CGF cloning-normalization estimate: "
          f"{final['CGFCloningNormalization']}")


if __name__ == '__main__':
    main()
