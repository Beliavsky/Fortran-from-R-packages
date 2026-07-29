#!/usr/bin/env python3
"""Generate independent NumPy references used by test_moments.f90."""
import numpy as np

r = np.array([
    [0.012, -0.004, 0.009],
    [-0.006, 0.011, 0.003],
    [0.018, 0.002, -0.005],
    [-0.009, -0.007, 0.014],
    [0.004, 0.015, -0.002],
    [0.021, -0.003, 0.006],
], dtype=float)
w = np.array([0.40, 0.35, 0.25], dtype=float)

z = r - r.mean(axis=0)
t, n = z.shape
m2 = z.T @ z / (t - 1)
m3 = np.empty((n, n**2))
m4 = np.empty((n, n**3))
for a in range(n):
    for b in range(n):
        for c in range(n):
            m3[a, b*n + c] = np.mean(z[:, a] * z[:, b] * z[:, c])
            for d in range(n):
                m4[a, (b*n + c)*n + d] = np.mean(
                    z[:, a] * z[:, b] * z[:, c] * z[:, d]
                )

ww = np.kron(w, w)
www = np.kron(ww, w)
pm2 = w @ m2 @ w
pm3 = w @ m3 @ ww
pm4 = w @ m4 @ www
dm2 = 2 * m2 @ w
dm3 = 3 * m3 @ ww
dm4 = 4 * m4 @ www
port_skew = pm3 / pm2**1.5
port_kurt = pm4 / pm2**2
port_skew_deriv = (
    pm2**1.5 * dm3 - pm3 * np.sqrt(pm2) * dm2
) / pm2**3
port_kurt_deriv = (pm2 * dm4 - pm4 * dm2) / (2 * pm2**3)

for name, value in {
    "M2": m2,
    "M3": m3,
    "M4": m4,
    "pm2": pm2,
    "pm3": pm3,
    "pm4": pm4,
    "dm2": dm2,
    "dm3": dm3,
    "dm4": dm4,
    "PortSkew": port_skew,
    "PortKurt": port_kurt,
    "PortSkewDeriv": port_skew_deriv,
    "PortKurtDeriv": port_kurt_deriv,
}.items():
    print(f"{name}:\n{value!r}\n")
