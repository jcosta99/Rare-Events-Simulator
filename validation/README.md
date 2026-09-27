# Exact reference material

Small-system results the simulator can be checked against, and the figures drawn from them. 

```text
validation/
  exact/SSEP_TM.py              exact diagonalisation of the tilted generator
  data/SSEP_TM_L<L>.csv         its output, one row per counting field
  figures/                      the large-L benchmark
```

## Exact diagonalisation

[`exact/SSEP_TM.py`](exact/SSEP_TM.py) builds the tilted generator of an open SSEP on a small lattice and takes its dominant eigenvalue, which is the scaled cumulant generating function of the current carrying no statistical error of any kind. That makes it the sharpest check available. The matrix is `2**L` by `2**L`, so the method runs out somewhere around `L = 14`.

The counting field is applied the way the engine applies it: spread uniformly over the bonds, `exp(+s/ell)` on a right hop and `exp(-s/ell)` on a left one, with the reservoir moves carrying no factor. The conjugate observable is then the same bond-averaged transported charge the engine reports.

Regenerate the stored tables with

```bash
python validation/exact/SSEP_TM.py
```

which writes `data/SSEP_TM_L<L>.csv` for the reservoir rates named in each file's header. The defaults are `alpha = beta = 1` with `gamma = delta = 0`, holding the reservoirs at `rho_L = 1` and `rho_R = 0`.

## Large-L benchmark

[`figures/`](figures) contains a benchmark of the engine against the large-L closed form of the open SSEP cumulant generating function, which holds with no external field and equal reservoir densities, to leading order in `L` (see the theory notes [`theory_algorithm_implementation.pdf`](../Simulator/docs/theory_algorithm_implementation.pdf)). The runs behind it are not included in the repository.
