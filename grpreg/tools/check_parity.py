#!/usr/bin/env python3
"""Independent NumPy/SciPy checks for the translated grpreg numerical core."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.special import gammaln

path = Path(sys.argv[1] if len(sys.argv) > 1 else "fortran_parity.csv")
out = pd.read_csv(path).set_index("case")
n = 40
i = np.arange(1, n + 1)
x = np.column_stack(((i - 20.5) / 12.0, np.sin(0.37 * i), np.cos(0.23 * i), np.sin(0.11 * i * i)))
yg = 1.25 + 1.7*x[:,0] - 0.9*x[:,1] + 0.6*x[:,2] - 0.3*x[:,3]
yb = ((i % 5 == 0) | (i % 7 == 0) | (x[:,1] > 0.65)).astype(float)
yp = (i % 4).astype(float)
time = ((17*i) % 41 + 1).astype(float) + 0.01*i
event = (((5*i + 2) % 7) < 4).astype(float)

# Independent block-coordinate solution of the group-lasso objective.
center = x.mean(axis=0)
scale = np.sqrt(((x-center)**2).mean(axis=0))
xstd = (x-center)/scale
basis = np.zeros((4,4))
col = 0
for inds in ([0,1], [2,3]):
    gram = xstd[:,inds].T @ xstd[:,inds]
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    for value, vector in zip(values[order], vectors[:,order].T):
        if value > 1e-20:
            basis[np.asarray(inds),col] = vector*np.sqrt(n/value)
            col += 1
basis = basis[:,:col]
xw = xstd @ basis
b = np.zeros(col)
r = yg-yg.mean()
for _ in range(100000):
    change = 0.0
    for g, inds in enumerate(([0,1], [2,3])):
        inds = np.asarray(inds)
        old = b[inds].copy()
        z = xw[:,inds].T @ r/n + old
        zn = np.linalg.norm(z)
        threshold = 0.15*np.sqrt(2.0)
        new = np.maximum(zn-threshold,0.0)/zn*z if zn > 0 else np.zeros_like(z)
        b[inds] = new
        r -= xw[:,inds] @ (new-old)
        change = max(change, np.max(np.abs(new-old)))
    if change < 1e-13:
        break
raw_beta = basis @ b / scale
raw_intercept = yg.mean() - center @ raw_beta
rss = np.sum((yg-(raw_intercept+x@raw_beta))**2)
row = out.loc["gaussian_group_lasso"]
fref = np.r_[row.intercept,row.beta1,row.beta2,row.beta3,row.beta4,row.deviance]
pref = np.r_[raw_intercept,raw_beta,rss]
gaussian_diff = np.max(np.abs(fref-pref))

# Near-zero lambda should agree with independently optimized GLM likelihoods.
def bin_nll(theta):
    eta = theta[0] + x @ theta[1:]
    return np.sum(np.logaddexp(0.0,eta)-yb*eta)
bin_fit = minimize(bin_nll,np.zeros(5),method="BFGS",options={"gtol":1e-10,"maxiter":10000})
row = out.loc["binomial_near_mle"]
bin_fortran = np.r_[row.intercept,row.beta1,row.beta2,row.beta3,row.beta4]
binomial_diff = np.max(np.abs(bin_fit.x-bin_fortran))

def pois_nll(theta):
    eta = theta[0] + x @ theta[1:]
    return np.sum(np.exp(np.clip(eta,-50,50))-yp*eta+gammaln(yp+1.0))
pois_fit = minimize(pois_nll,np.zeros(5),method="BFGS",options={"gtol":1e-10,"maxiter":10000})
row = out.loc["poisson_near_mle"]
pois_fortran = np.r_[row.intercept,row.beta1,row.beta2,row.beta3,row.beta4]
poisson_diff = np.max(np.abs(pois_fit.x-pois_fortran))

order = np.argsort(time,kind="stable")
xs = x[order]
ds = event[order]
def cox_nll(beta):
    eta = xs @ beta
    risk = np.cumsum(np.exp(eta[::-1]))[::-1]
    return -np.sum(ds*(eta-np.log(risk)))
cox_fit = minimize(cox_nll,np.zeros(4),method="BFGS",options={"gtol":1e-10,"maxiter":10000})
row = out.loc["cox_near_mle"]
cox_fortran = np.r_[row.beta1,row.beta2,row.beta3,row.beta4]
cox_diff = np.max(np.abs(cox_fit.x-cox_fortran))
cox_dev_diff = abs(2.0*cox_fit.fun-row.deviance)

# B-spline expansion against SciPy's independent design-matrix implementation.
x1 = x[:,0]
boundary = (x1.min(),x1.max())
interior = np.quantile(x1,[0.5],method="linear")
knots = np.r_[np.repeat(boundary[0],4),interior,np.repeat(boundary[1],4)]
bs = BSpline.design_matrix(x1,knots,k=3,extrapolate=False).toarray()[:,1:]
row = out.loc["spline_row"]
spline_fortran = np.r_[row.beta1,row.beta2,row.beta3,row.beta4]
spline_diff = np.max(np.abs(bs[12]-spline_fortran))

limits = {
    "gaussian_group_lasso": (gaussian_diff, 5e-9),
    "binomial_near_mle": (binomial_diff, 5e-6),
    "poisson_near_mle": (poisson_diff, 5e-6),
    "cox_near_mle": (cox_diff, 5e-6),
    "cox_deviance": (cox_dev_diff, 5e-6),
    "bspline": (spline_diff, 5e-12),
}
for name, (diff, limit) in limits.items():
    print(f"{name}: max_abs_diff={diff:.3e} limit={limit:.1e}")
    if not np.isfinite(diff) or diff > limit:
        raise SystemExit(f"Parity check failed for {name}")
print("Independent grpreg NumPy/SciPy parity checks passed.")
