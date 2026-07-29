#!/usr/bin/env python3
"""Generate the fixed two-dimensional BEKK references used by test_core.f90."""

import numpy as np

DATA = np.array([
    [0.12, -0.08],
    [-0.03, 0.11],
    [0.07, -0.02],
    [-0.15, -0.05],
    [0.04, 0.09],
    [-0.06, -0.12],
    [0.10, 0.03],
    [-0.02, 0.07],
], dtype=float)

C = np.array([[0.20, 0.00], [0.05, 0.15]])
A = np.array([[0.25, -0.01], [0.02, 0.20]])
B = np.array([[0.10, 0.00], [0.00, 0.08]])
G = np.array([[0.80, 0.00], [0.01, 0.75]])
SIGNS = np.array([-1.0, -1.0])


def indicator(x: np.ndarray) -> float:
    return float(np.all(SIGNS * x >= 0.0))


def loglike(model: str, asymmetric: bool = False) -> float:
    h = DATA.T @ DATA / DATA.shape[0]
    total = 0.0
    for row in DATA:
        sign, logdet = np.linalg.slogdet(h)
        assert sign > 0.0
        total += -0.5 * (DATA.shape[1] * np.log(2.0 * np.pi) + logdet + row @ np.linalg.solve(h, row))
        if model == "full":
            h = C @ C.T + A.T @ np.outer(row, row) @ A + G.T @ h @ G
            if asymmetric:
                h += indicator(row) * B.T @ np.outer(row, row) @ B
        elif model == "diagonal":
            ad = np.diag(np.diag(A))
            gd = np.diag(np.diag(G))
            h = C @ C.T + ad.T @ np.outer(row, row) @ ad + gd.T @ h @ gd
        elif model == "scalar":
            a = 0.10
            g = 0.75 if asymmetric else 0.80
            b = 0.05 if asymmetric else 0.0
            h = C @ C.T + (a + indicator(row) * b) * np.outer(row, row) + g * h
        h = 0.5 * (h + h.T)
    return total


print("full", loglike("full"))
print("asymmetric full", loglike("full", True))
print("diagonal", loglike("diagonal"))
print("scalar", loglike("scalar"))
print("asymmetric scalar", loglike("scalar", True))
