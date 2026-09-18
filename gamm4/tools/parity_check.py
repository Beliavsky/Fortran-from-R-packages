from pathlib import Path
import csv
import numpy as np
from scipy.optimize import minimize

csv_path = Path(__file__).resolve().parents[1] / "fortran_parity.csv"
if not csv_path.exists():
    csv_path = Path("fortran_parity.csv")
values = {row["name"]: float(row["value"]) for row in csv.DictReader(csv_path.open())}

n = 18
x = np.linspace(-1.0, 1.0, n)
idx = np.arange(n)
ii = np.arange(1, n + 1)
group_effect = np.where(idx < n // 2, -0.28, 0.28)
y = 0.9 + 1.1*x + 0.75*x*x - 0.55*x**3 + group_effect + 0.025*np.sin(1.31*ii)
x_fixed = np.column_stack((np.ones(n), x))
z_smooth = np.column_stack((x*x, x**3))
z_group = np.column_stack((idx < n//2, idx >= n//2)).astype(float)
z = np.column_stack((z_smooth, z_group))

def evaluate(eta, details=False):
    g = np.diag([np.exp(2*eta[0])]*2 + [np.exp(2*eta[1])]*2)
    v = np.eye(n) + z @ g @ z.T
    sign, logdet = np.linalg.slogdet(v)
    if sign <= 0:
        return 1e100
    vinv = np.linalg.inv(v)
    beta = np.linalg.solve(x_fixed.T @ vinv @ x_fixed, x_fixed.T @ vinv @ y)
    residual = y - x_fixed @ beta
    scale = (residual @ vinv @ residual) / n
    objective = n*(np.log(2*np.pi) + 1.0 + np.log(scale)) + logdet
    u = g @ z.T @ vinv @ residual
    if details:
        return objective, scale, beta, u
    return objective

opt = minimize(evaluate, np.log([0.5, 0.5]), method="L-BFGS-B", bounds=[(-8, 4), (-8, 4)],
               options={"ftol": 1e-14, "gtol": 1e-11, "maxiter": 1000})
objective, scale, beta, u = evaluate(opt.x, True)
coef_ref = np.r_[beta, u[:2]]
coef_ftn = np.array([values[f"coef{k}"] for k in range(1, 5)])
eta_ftn = np.array([values["eta1"], values["eta2"]])

# Reconstruct upstream getVb at the *Fortran* optimum, so this check isolates
# covariance assembly from the small optimizer-tolerance difference.
sp = np.exp(-2*eta_ftn[0])
scale_ftn = values["scale"]
x_orig = np.column_stack((np.ones(n), x, x*x, x**3))
phi = np.eye(2) * np.exp(2*eta_ftn[1])
v = np.eye(n)*scale_ftn + scale_ftn*z_group @ phi @ z_group.T
vinv = np.linalg.inv(v)
xvx = x_orig.T @ vinv @ x_orig
penalty = np.diag([0.0, 0.0, sp/scale_ftn, sp/scale_ftn])
vcov = np.linalg.inv(xvx + penalty)
vcov_diag_ftn = np.array([values[f"vcov_diag{k}"] for k in range(1, 5)])

checks = {
    "objective": abs(values["objective"] - objective),
    "scale": abs(values["scale"] - scale),
    "coefficients": np.max(np.abs(coef_ftn - coef_ref)),
    "eta": np.max(np.abs(eta_ftn - opt.x)),
    "vcov_diag": np.max(np.abs(vcov_diag_ftn - np.diag(vcov))),
}
for name, difference in checks.items():
    print(f"{name} max/abs difference: {difference:.3e}")
assert checks["objective"] < 2e-5
assert checks["scale"] < 2e-6
assert checks["coefficients"] < 2e-4
assert checks["eta"] < 2e-3
assert checks["vcov_diag"] < 2e-10
print("Independent NumPy/SciPy parity checks passed.")
