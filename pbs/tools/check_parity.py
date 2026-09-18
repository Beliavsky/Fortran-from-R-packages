#!/usr/bin/env python3
"""Independent SciPy parity checks for the translated pbs basis."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline


def basis_matrix(knots: np.ndarray, degree: int, x: np.ndarray) -> np.ndarray:
    ncoef = knots.size - degree - 1
    eye = np.eye(ncoef)
    return np.column_stack([BSpline(knots, eye[j], degree, extrapolate=False)(x) for j in range(ncoef)])


def periodic_reference(x: np.ndarray) -> np.ndarray:
    degree = 3
    base = np.array([0.0, 1.0, 2.0, 3.0, 6.0])
    close = base.copy()
    for i in range(1, degree + 1):
        close = np.r_[close, base[-1] + base[i] - base[0]]
    for i in range(1, degree + 1):
        close = np.r_[base[0] - base[-1] + base[-1 - i], close]
    raw = basis_matrix(close, degree, x)
    ncoef = raw.shape[1]
    middle = raw[:, degree : ncoef - degree]
    folded = raw[:, :degree] + raw[:, ncoef - degree :]
    return np.column_stack([middle, folded])


def ordinary_reference(x: np.ndarray) -> np.ndarray:
    degree = 2
    boundary = np.array([0.0, 2.5])
    knots = np.r_[np.repeat(boundary[0], degree + 1), [1.0, 2.0], np.repeat(boundary[1], degree + 1)]
    ncoef = knots.size - degree - 1
    eye = np.eye(ncoef)
    ans = np.zeros((x.size, ncoef))
    for i, xx in enumerate(x):
        if boundary[0] <= xx <= boundary[1]:
            ans[i] = [BSpline(knots, eye[j], degree, extrapolate=False)(xx) for j in range(ncoef)]
        else:
            pivot = boundary[0] if xx < boundary[0] else boundary[1]
            delta = xx - pivot
            for r in range(degree + 1):
                ans[i] += delta**r / math.factorial(r) * np.array(
                    [BSpline(knots, eye[j], degree, extrapolate=False).derivative(r)(pivot) for j in range(ncoef)]
                )
    return ans


def main() -> int:
    path = Path(sys.argv[1] if len(sys.argv) > 1 else 'fortran_parity.csv')
    rows = list(csv.DictReader(path.open(newline='')))
    diffs = {}
    for case in ('periodic', 'ordinary'):
        sub = [r for r in rows if r['case'] == case]
        nr = max(int(r['row']) for r in sub)
        nc = max(int(r['column']) for r in sub)
        x = np.zeros(nr)
        got = np.zeros((nr, nc))
        for r in sub:
            i = int(r['row']) - 1
            j = int(r['column']) - 1
            x[i] = float(r['x'])
            got[i, j] = float(r['value'])
        ref = periodic_reference(x) if case == 'periodic' else ordinary_reference(x)
        diffs[case] = float(np.max(np.abs(got - ref)))
    print(f"periodic max abs diff: {diffs['periodic']:.3e}")
    print(f"ordinary max abs diff: {diffs['ordinary']:.3e}")
    if diffs['periodic'] > 5e-13 or diffs['ordinary'] > 5e-13:
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
