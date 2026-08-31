#!/usr/bin/env python3

"""Exact cumulant generating function of the open SSEP, by diagonalisation.

The tilted generator of a small chain is written out in full and its dominant
eigenvalue taken.  That eigenvalue is the scaled cumulant generating function
of the current, carrying no statistical error whatever, which is what makes it
the sharpest available check on the simulator.  The cost is exponential in the
lattice --- the matrix is ``2**L`` on a side --- so ``L`` beyond about fourteen
is out of reach.

The counting field is applied exactly as the engine applies it: spread
uniformly over the bonds, so that a right hop carries ``exp(+s/ell)`` and a
left hop ``exp(-s/ell)`` with ``ell`` the number of bonds, while the reservoir
moves carry no factor at all.  The conjugate observable is therefore the same
one the engine reports, the bond-averaged transported charge, and the two sets
of numbers can be placed side by side with no conversion between them.

Regenerate the reference data with:

    python validation/exact/SSEP_TM.py
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs


DATA_DIR = Path(__file__).resolve().parents[1] / 'data'


def tilted_generator(s, L, alpha=1.0, beta=1.0, gamma=0.0, delta=0.0,
                     boundary='open'):
    """The tilted Markov generator, as a sparse matrix acting on columns.

    Site ``j`` is bit ``L-1-j`` of the state index, so state ``0`` is the empty
    lattice and the leftmost site is the most significant bit.  Every allowed
    hop has rate one; the reservoirs inject and remove at the four boundary
    rates, ``alpha`` and ``gamma`` on the left and ``delta`` and ``beta`` on
    the right.  The diagonal carries the untilted escape rate, since the
    counting field reweights transitions and not the waiting times.
    """
    if boundary not in ('open', 'periodic'):
        raise ValueError("boundary must be 'open' or 'periodic'")

    bonds = L if boundary == 'periodic' else L - 1
    right_weight, left_weight = np.exp(s / bonds), np.exp(-s / bonds)

    states = 1 << L
    rows, cols, values = [], [], []
    escape = np.zeros(states)

    def occupied(state, site):
        return (state >> (L - 1 - site)) & 1

    for state in range(states):
        for site in range(L):
            if not occupied(state, site):
                continue
            for target, weight in ((site + 1, right_weight),
                                   (site - 1, left_weight)):
                if boundary == 'periodic':
                    target %= L
                elif not 0 <= target < L:
                    continue
                if occupied(state, target):
                    continue
                rows.append(state ^ (1 << (L - 1 - site))
                                  ^ (1 << (L - 1 - target)))
                cols.append(state)
                values.append(weight)
                escape[state] += 1.0

        if boundary != 'open':
            continue
        for site, inject, remove in ((0, alpha, gamma), (L - 1, delta, beta)):
            rate = remove if occupied(state, site) else inject
            if rate <= 0.0:
                continue
            rows.append(state ^ (1 << (L - 1 - site)))
            cols.append(state)
            values.append(rate)
            escape[state] += rate

    rows.extend(range(states))
    cols.extend(range(states))
    values.extend(-escape)
    return coo_matrix((values, (rows, cols)), shape=(states, states)).tocsc()


def cgf(s, L, **parameters):
    """The dominant eigenvalue of the tilted generator."""
    generator = tilted_generator(s, L, **parameters)
    if generator.shape[0] <= 256:                     # dense is faster here
        return float(np.max(np.linalg.eigvals(generator.toarray()).real))
    value = eigs(generator, k=1, which='LR', return_eigenvectors=False,
                 maxiter=100000, tol=1e-13)
    return float(value[0].real)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--lengths', type=int, nargs='+', default=(5, 6, 8, 10))
    parser.add_argument('--points', type=int, default=41)
    parser.add_argument('--range', type=float, default=4.0)
    parser.add_argument('--alpha', type=float, default=1.0)
    parser.add_argument('--beta', type=float, default=1.0)
    parser.add_argument('--gamma', type=float, default=0.0)
    parser.add_argument('--delta', type=float, default=0.0)
    args = parser.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    rates = dict(alpha=args.alpha, beta=args.beta,
                 gamma=args.gamma, delta=args.delta)
    fields = np.linspace(-args.range, args.range, args.points)

    for L in args.lengths:
        destination = DATA_DIR / f'SSEP_TM_L{L}.csv'
        with destination.open('w') as stream:
            stream.write('# open SSEP, exact diagonalisation of the tilted '
                         'generator\n')
            stream.write(f'# alpha={args.alpha:g} gamma={args.gamma:g} '
                         f'delta={args.delta:g} beta={args.beta:g}\n')
            stream.write('bias_s,cgf\n')
            for s in fields:
                stream.write(f'{s:.10g},{cgf(s, L, **rates):.12g}\n')
        print(f'L={L:>3}: wrote {destination.name} '
              f'({args.points} points, 2**{L} states)')


if __name__ == '__main__':
    main()
