# Testing

## Configurations

The included `build_gfortran.sh` runs two supported configurations:

### Debug

```text
-O0 -g -std=f2018 -Wall -Wextra -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

### Release

```text
-O3 -std=f2018 -Wall -Wextra -Werror
```

Both link with LAPACK and BLAS.

## Test programs

### `test_linear_quadratic`

- equality-constrained quadratic-program reference `x=(0.25,0.75)`;
- objective reference `1.875`;
- inequality-constrained LP reference objective `-9`;
- `dqp`/`cps` compatibility workflow; and
- solution extractor behavior.

### `test_cones`

- the original two-cone SOCP unit-test problem;
- the original two-block SDP unit-test problem;
- Lorentz feasibility; and
- minimum PSD eigenvalue feasibility.

### `test_mixed_cones`

- the original LP example containing NNOC, SOCC, and PSDC blocks
  simultaneously.

### `test_nonlinear`

- analytical centering with a linear equality;
- callback objective gradients and Hessians;
- linear objective with reciprocal nonlinear constraint;
- nonlinear phase-I recovery from an infeasible but in-domain initial point; and
- compatibility names `dcp` and `dnl`.

### `test_special`

- L1 regression;
- equal-risk-contribution portfolio identities;
- geometric-programming status and positivity; and
- an independently solved SciPy/SLSQP reference
  `(1.2866774443, 0.5304618318)`.

### `test_cone_algebra`

- Lorentz and PSD identities;
- Jordan products and inverses;
- cone inner products; and
- maximum-step/interiority calculations.

## Original package coverage

The translated tests cover the numerical categories in the original RUnit suite:

- LP with inequalities;
- QP with equalities and inequalities;
- SOCP;
- SDP;
- mixed cone programs;
- nonlinear convex programs;
- geometric programming; and
- specialized L1 and risk-parity programs.

The original test and demo files remain under `original/` for provenance.
