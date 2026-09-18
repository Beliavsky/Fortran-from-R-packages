#!/usr/bin/env python3
"""Independent NumPy/statsmodels checks for selected translated FlexMix kernels."""
from pathlib import Path
import csv
import math
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad
from scipy.stats import norm
from sklearn.linear_model import Lasso

here = Path(__file__).resolve().parent
path = Path.cwd() / "fortran_parity.csv"
if not path.exists():
    path = here / "fortran_parity.csv"
with path.open(newline="") as f:
    got = {row["metric"]: float(row["value"]) for row in csv.DictReader(f)}

x = np.column_stack([np.ones(5), np.arange(-2.0, 3.0)])
yg = np.array([1.2, 1.8, 3.1, 4.2, 5.1])
yp = np.array([1.0, 1.0, 2.0, 4.0, 7.0])
success = np.array([1.0, 2.0, 4.0, 7.0, 9.0])
trials = np.full(5, 10.0)
ym = np.column_stack([
    np.array([-1.0, 0.0, 1.0, 2.0, 3.0]),
    np.array([2.0, 1.0, 4.0, 3.0, 5.0]),
])

beta_g = np.linalg.lstsq(x, yg, rcond=None)[0]
resid = yg - x @ beta_g
sigma_g = math.sqrt(np.dot(resid, resid) / (len(yg) - x.shape[1]))
poisson = sm.GLM(yp, x, family=sm.families.Poisson()).fit()
binomial = sm.GLM(success / trials, x, family=sm.families.Binomial(), freq_weights=trials).fit()
cov = np.cov(ym, rowvar=False, ddof=1)
center = ym.mean(axis=0)

y_gamma = np.array([1.0, 2.0, 4.0, 8.0])
x_gamma = np.ones((len(y_gamma), 1))
gamma_model = sm.GLM(
    y_gamma,
    x_gamma,
    family=sm.families.Gamma(link=sm.families.links.InversePower()),
).fit()
mu_gamma = gamma_model.fittedvalues
gamma_deviance = np.sum(2.0 * ((y_gamma - mu_gamma) / mu_gamma - np.log(y_gamma / mu_gamma)))
gamma_shape = len(y_gamma) / gamma_deviance

y_offset = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
exposure = np.array([1.0, 1.0, 2.0, 2.0, 4.0])
poisson_offset = sm.GLM(
    y_offset,
    np.ones((len(y_offset), 1)),
    family=sm.families.Poisson(),
    offset=np.log(exposure),
).fit()

class_count = np.array([2.0, 3.0, 5.0])
class_probability = class_count / class_count.sum()
multinomial_beta = np.log(class_probability[1:] / class_probability[0])

tfix = np.tile(np.arange(10.0) / 9.0, 2)
yfix = np.concatenate([
    1.0 + 2.0 * tfix[:10] + 0.02 * np.array([(-1.0) ** i for i in range(1, 11)]),
    8.0 + 2.0 * tfix[10:] + 0.02 * np.array([(-1.0) ** i for i in range(11, 21)]),
])
zfix = np.column_stack([
    np.r_[np.ones(10), np.zeros(10)],
    np.r_[np.zeros(10), np.ones(10)],
    tfix,
])
beta_fix = np.linalg.lstsq(zfix, yfix, rcond=None)[0]

idx = np.arange(1.0, 61.0)
z1 = np.sin(0.37 * idx) + 0.3 * np.cos(0.11 * idx)
z2 = np.cos(0.29 * idx)
y_factor = np.column_stack([
    2.0 + z1 + 0.18 * np.sin(1.7 * idx),
    -1.0 + 0.85 * z1 + 0.25 * np.cos(1.3 * idx),
    0.5 + 0.65 * z1 + 0.35 * np.sin(1.1 * idx),
    3.0 + 0.45 * z1 + 0.55 * z2,
])
factor_corr = np.corrcoef(y_factor, rowvar=False)
factor_inv = np.linalg.inv(factor_corr)
factor_start = np.clip((1.0 - 0.5 / 4.0) / np.diag(factor_inv), 0.005, 1.0)

def factor_objective_gradient(uniqueness):
    scale = 1.0 / np.sqrt(uniqueness)
    scaled = factor_corr * np.outer(scale, scale)
    eigenvalues, eigenvectors = np.linalg.eigh(scaled)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    objective = -(np.sum(np.log(eigenvalues[1:]) - eigenvalues[1:])) + (1.0 - 4.0)
    loadings = (
        eigenvectors[:, :1]
        * np.sqrt(np.maximum(eigenvalues[:1] - 1.0, 0.0))[None, :]
        * np.sqrt(uniqueness)[:, None]
    )
    gradient = (np.diag(loadings @ loadings.T) + uniqueness - 1.0) / uniqueness**2
    return objective, gradient


x_cond = np.tile(np.array([-1.0, 0.0, 1.0]), 12)[:, None]
y_cond = np.zeros(36)
selected_cond = [0] * 3 + [1] * 4 + [2] * 5
for group, choice in enumerate(selected_cond):
    y_cond[3 * group + choice] = 1.0

def conditional_logit_nll(beta):
    eta = x_cond[:, 0] * beta[0]
    value = 0.0
    for group in range(12):
        rows = slice(3 * group, 3 * group + 3)
        e = eta[rows]
        value -= e[np.argmax(y_cond[rows])] - np.log(np.exp(e - e.max()).sum()) - e.max()
    return value

cond_fit = minimize_scalar(
    lambda value: conditional_logit_nll(np.array([value])),
    bounds=(-10.0, 10.0),
    method="bounded",
    options={"xatol": 1.0e-13},
)
if not cond_fit.success:
    raise SystemExit(f"SciPy conditional-logit reference optimization failed: {cond_fit.message}")

factor_fit = minimize(
    lambda u: factor_objective_gradient(u)[0],
    factor_start,
    jac=lambda u: factor_objective_gradient(u)[1],
    method="L-BFGS-B",
    bounds=[(0.005, 1.0)] * 4,
    options={"ftol": 1.0e-13, "gtol": 1.0e-10, "maxiter": 1000},
)
if not factor_fit.success:
    raise SystemExit(f"SciPy factor-analysis reference optimization failed: {factor_fit.message}")


# Fixed-lambda Gaussian lasso/elastic-net reference for FLXMRglmnet.
idx_pen = np.arange(1.0, 21.0)
x_pen = np.column_stack([
    np.ones(20),
    (idx_pen - 10.5) / 5.0,
    np.sin(idx_pen),
    np.cos(0.3 * idx_pen),
])
y_pen = 2.0 + 3.0 * x_pen[:, 1] + 0.02 * np.sin(2.0 * idx_pen)
lasso = Lasso(alpha=0.05, fit_intercept=True, max_iter=100000, tol=1.0e-13, selection="cyclic")
lasso.fit(x_pen[:, 1:], y_pen)

# Fixed-smoothing-parameter references for the translated FLXMRmgcv numerical core.
smooth_i = np.arange(1.0, 11.0)
smooth_t = -1.0 + 2.0 * np.arange(10.0) / 9.0
x_smooth = np.column_stack([np.ones(10), smooth_t, smooth_t**2])
y_smooth = 1.0 + 2.0*smooth_t + 0.5*smooth_t**2 + 0.05*np.sin(smooth_i)
y_smooth_pois = np.mod(smooth_i, 4.0) + 1.0
smooth_success = np.mod(smooth_i, 3.0)
smooth_trials = np.full(10, 2.0)
smooth_penalty = np.diag([0.0, 0.0, 1.0])
smooth_lambda = 5.0
smooth_system = x_smooth.T @ x_smooth + smooth_lambda*smooth_penalty
smooth_gbeta = np.linalg.solve(smooth_system, x_smooth.T @ y_smooth)
smooth_hat_df = np.trace(np.linalg.solve(smooth_system, x_smooth.T @ x_smooth))
smooth_resid = y_smooth - x_smooth @ smooth_gbeta
smooth_sigma = np.sqrt(np.dot(smooth_resid, smooth_resid) / (10.0 - smooth_hat_df))
smooth_edf = smooth_hat_df + 1.0

def smooth_poisson_objective(beta):
    eta = np.clip(x_smooth @ beta, -30.0, 30.0)
    mu = np.exp(eta)
    return np.sum(mu - y_smooth_pois*eta) + 0.5*smooth_lambda*beta @ smooth_penalty @ beta

def smooth_poisson_gradient(beta):
    eta = np.clip(x_smooth @ beta, -30.0, 30.0)
    mu = np.exp(eta)
    return x_smooth.T @ (mu-y_smooth_pois) + smooth_lambda*smooth_penalty @ beta

smooth_pois_fit = minimize(
    smooth_poisson_objective,
    np.array([np.log(y_smooth_pois.mean()), 0.0, 0.0]),
    jac=smooth_poisson_gradient,
    method="BFGS",
    options={"gtol": 1.0e-12, "maxiter": 1000},
)
if not smooth_pois_fit.success and np.max(np.abs(smooth_pois_fit.jac)) > 1.0e-8:
    raise SystemExit(f"SciPy smooth-Poisson reference optimization failed: {smooth_pois_fit.message}")

def smooth_binomial_objective(beta):
    eta = x_smooth @ beta
    value = np.sum(smooth_trials*np.logaddexp(0.0, eta) - smooth_success*eta)
    return value + 0.5*smooth_lambda*beta @ smooth_penalty @ beta

def smooth_binomial_gradient(beta):
    eta = x_smooth @ beta
    prob = np.where(eta >= 0.0, 1.0/(1.0+np.exp(-eta)), np.exp(eta)/(1.0+np.exp(eta)))
    return x_smooth.T @ (smooth_trials*prob-smooth_success) + smooth_lambda*smooth_penalty @ beta

smooth_binom_fit = minimize(
    smooth_binomial_objective,
    np.zeros(3),
    jac=smooth_binomial_gradient,
    method="BFGS",
    options={"gtol": 1.0e-12, "maxiter": 1000},
)
if not smooth_binom_fit.success and np.max(np.abs(smooth_binom_fit.jac)) > 1.0e-8:
    raise SystemExit(f"SciPy smooth-binomial reference optimization failed: {smooth_binom_fit.message}")

# One-component random-intercept mixed model reference for FLXMRlmm/lmer.
x_lmm = []
y_lmm = []
group_lmm = []
for group in range(1, 13):
    random_effect = 0.8 * np.sin(0.7 * group)
    for rep in range(1, 5):
        xr = (rep - 1.0) / 3.0
        x_lmm.append([1.0, xr])
        group_lmm.append(group)
        y_lmm.append(2.0 + 1.5 * xr + random_effect + 0.15 * np.cos(3.0 * group + rep))
x_lmm = np.asarray(x_lmm)
y_lmm = np.asarray(y_lmm)
mixed = sm.MixedLM(y_lmm, x_lmm, groups=np.asarray(group_lmm)).fit(reml=False, method="lbfgs", disp=False)

# Left-censored random-intercept mixed model. The first row of each group is
# censored from above at a threshold 0.22 above the generated latent value.
y_lmmc = y_lmm.copy()
censored_lmm = np.zeros(len(y_lmmc), dtype=bool)
censored_lmm[::4] = True
y_lmmc[censored_lmm] += 0.22
group_lmm_array = np.asarray(group_lmm)

def lmmc_loglik(theta):
    beta = theta[:2]
    random_var = math.exp(theta[2])
    residual_var = math.exp(theta[3])
    value = 0.0
    covariance = random_var * np.ones((4, 4)) + residual_var * np.eye(4)
    for group_value in range(1, 13):
        idx = np.flatnonzero(group_lmm_array == group_value)
        mean = x_lmm[idx] @ beta
        observed = np.array([1, 2, 3])
        residual = y_lmmc[idx][observed] - mean[observed]
        cov_obs = covariance[np.ix_(observed, observed)]
        sign, logdet = np.linalg.slogdet(cov_obs)
        inv_obs = np.linalg.inv(cov_obs)
        value += -0.5 * (3.0 * np.log(2.0 * np.pi) + logdet + residual @ inv_obs @ residual)
        cross = covariance[0, observed]
        cond_mean = mean[0] + cross @ inv_obs @ residual
        cond_var = covariance[0, 0] - cross @ inv_obs @ cross
        value += norm.logcdf((y_lmmc[idx][0] - cond_mean) / math.sqrt(cond_var))
    return value

lmmc_start = np.array([
    got["lmmc_beta0"],
    got["lmmc_beta1"],
    math.log(got["lmmc_random_var"]),
    math.log(got["lmmc_residual_var"]),
])
lmmc_opt = minimize(
    lambda theta: -lmmc_loglik(theta),
    lmmc_start,
    method="Nelder-Mead",
    options={"maxiter": 4000, "xatol": 1.0e-11, "fatol": 1.0e-11},
)
if not lmmc_opt.success:
    raise SystemExit(f"SciPy censored mixed-model reference optimization failed: {lmmc_opt.message}")

def lmc_loglik(theta):
    beta = theta[:2]
    residual_var = math.exp(theta[2])
    sigma = math.sqrt(residual_var)
    mean = x_lmm @ beta
    value = np.sum(
        -0.5 * (np.log(2.0 * np.pi * residual_var) + ((y_lmmc[~censored_lmm] - mean[~censored_lmm]) ** 2) / residual_var)
    )
    value += np.sum(norm.logcdf((y_lmmc[censored_lmm] - mean[censored_lmm]) / sigma))
    return value

lmc_start = np.array([got["lmc_beta0"], got["lmc_beta1"], math.log(got["lmc_residual_var"])])
lmc_opt = minimize(
    lambda theta: -lmc_loglik(theta),
    lmc_start,
    method="Nelder-Mead",
    options={"maxiter": 3000, "xatol": 1.0e-12, "fatol": 1.0e-12},
)
if not lmc_opt.success:
    raise SystemExit(f"SciPy censored regression reference optimization failed: {lmc_opt.message}")

expected = {
    "gaussian_beta0": beta_g[0],
    "gaussian_beta1": beta_g[1],
    "gaussian_sigma": sigma_g,
    "poisson_beta0": poisson.params[0],
    "poisson_beta1": poisson.params[1],
    "binomial_beta0": binomial.params[0],
    "binomial_beta1": binomial.params[1],
    "mvnorm_center1": center[0],
    "mvnorm_center2": center[1],
    "mvnorm_cov11": cov[0, 0],
    "mvnorm_cov12": cov[0, 1],
    "mvnorm_cov22": cov[1, 1],
    "gamma_inverse_beta0": gamma_model.params[0],
    "gamma_shape": gamma_shape,
    "poisson_offset_beta0": poisson_offset.params[0],
    "multinomial_beta_class2": multinomial_beta[0],
    "multinomial_beta_class3": multinomial_beta[1],
    "glmfix_intercept1": beta_fix[0],
    "glmfix_intercept2": beta_fix[1],
    "glmfix_shared_slope": beta_fix[2],
    "conditional_logit_beta1": cond_fit.x,
    "factanal_uniqueness1": factor_fit.x[0],
    "factanal_uniqueness2": factor_fit.x[1],
    "factanal_uniqueness3": factor_fit.x[2],
    "factanal_uniqueness4": factor_fit.x[3],
    "glmnet_gaussian_intercept": lasso.intercept_,
    "glmnet_gaussian_slope1": lasso.coef_[0],
    "glmnet_gaussian_slope2": lasso.coef_[1],
    "glmnet_gaussian_slope3": lasso.coef_[2],
    "lmm_beta0": mixed.fe_params[0],
    "lmm_beta1": mixed.fe_params[1],
    "lmm_random_var": mixed.cov_re[0, 0],
    "lmm_residual_var": mixed.scale,
    "lmm_loglik": mixed.llf,
    "lmmc_beta0": lmmc_opt.x[0],
    "lmmc_beta1": lmmc_opt.x[1],
    "lmmc_random_var": math.exp(lmmc_opt.x[2]),
    "lmmc_residual_var": math.exp(lmmc_opt.x[3]),
    "lmmc_loglik": -lmmc_opt.fun,
    "lmc_beta0": lmc_opt.x[0],
    "lmc_beta1": lmc_opt.x[1],
    "lmc_residual_var": math.exp(lmc_opt.x[2]),
    "lmc_loglik": -lmc_opt.fun,
    "mgcv_gaussian_beta0": smooth_gbeta[0],
    "mgcv_gaussian_beta1": smooth_gbeta[1],
    "mgcv_gaussian_beta2": smooth_gbeta[2],
    "mgcv_gaussian_sigma": smooth_sigma,
    "mgcv_gaussian_edf": smooth_edf,
    "mgcv_poisson_beta0": smooth_pois_fit.x[0],
    "mgcv_poisson_beta1": smooth_pois_fit.x[1],
    "mgcv_poisson_beta2": smooth_pois_fit.x[2],
    "mgcv_binomial_beta0": smooth_binom_fit.x[0],
    "mgcv_binomial_beta1": smooth_binom_fit.x[1],
    "mgcv_binomial_beta2": smooth_binom_fit.x[2],
}

for key, target in expected.items():
    err = abs(got[key] - target)
    tol = 2e-8 * (1.0 + abs(target))
    if key == "lmm_random_var":
        tol = 7.0e-5
    elif key == "lmm_residual_var":
        tol = 1.0e-6
    elif key == "lmm_loglik":
        tol = 5.0e-6
    elif key == "lmmc_beta0":
        tol = 8.0e-4
    elif key == "lmmc_beta1":
        tol = 3.0e-5
    elif key == "lmmc_random_var":
        tol = 3.0e-5
    elif key == "lmmc_residual_var":
        tol = 3.0e-6
    elif key == "lmmc_loglik":
        tol = 3.0e-5
    elif key == "lmc_beta0" or key == "lmc_beta1":
        tol = 3.0e-6
    elif key == "lmc_residual_var":
        tol = 5.0e-7
    elif key == "lmc_loglik":
        tol = 2.0e-7
    elif key.startswith("glmnet_gaussian"):
        tol = 2.0e-6 * (1.0 + abs(target))
    elif key.startswith("mgcv_poisson") or key.startswith("mgcv_binomial"):
        tol = 2.0e-8 * (1.0 + abs(target))
    if err > tol:
        raise SystemExit(
            f"parity failure {key}: Fortran={got[key]:.16g} "
            f"expected={target:.16g} error={err:.3g}"
        )
# Validate the deterministic multivariate-truncation likelihood at the exact
# parameters returned by the short Fortran multi-censor run. This checks the
# Genz-style QMC integration without requiring the deliberately short run to
# reach its final EM optimum.
def bivariate_normal_upper_probability(mean, covariance, upper):
    sd1 = math.sqrt(covariance[0, 0])
    sd2 = math.sqrt(covariance[1, 1])
    rho = covariance[0, 1] / (sd1 * sd2)
    rho = float(np.clip(rho, -0.999999999999, 0.999999999999))
    denom = math.sqrt(1.0 - rho * rho)
    z_upper = (upper[0] - mean[0]) / sd1

    def integrand(z):
        cond = (upper[1] - mean[1] - rho * sd2 * z) / (sd2 * denom)
        return math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi) * norm.cdf(cond)

    probability, _ = quad(integrand, -np.inf, z_upper, epsabs=2.0e-12, epsrel=2.0e-12, limit=200)
    return probability

x_multi = []
y_multi = []
group_multi = []
censor_multi = []
for group_value in range(1, 5):
    random_effect = 0.6 * np.sin(0.5 * group_value)
    for rep in range(1, 4):
        xr = (rep - 1.0) / 2.0
        response = 1.3 + 1.1 * xr + random_effect + 0.08 * np.cos(2.0 * group_value + rep)
        is_censored = rep <= 2
        if is_censored:
            response += 0.15
        x_multi.append([1.0, xr])
        y_multi.append(response)
        group_multi.append(group_value)
        censor_multi.append(is_censored)
x_multi = np.asarray(x_multi)
y_multi = np.asarray(y_multi)
group_multi = np.asarray(group_multi)
censor_multi = np.asarray(censor_multi)
beta_multi = np.array([got["lmmc_multi_beta0"], got["lmmc_multi_beta1"]])
random_var_multi = got["lmmc_multi_random_var"]
residual_var_multi = got["lmmc_multi_residual_var"]
exact_multi_loglik = 0.0
for group_value in range(1, 5):
    idx = np.flatnonzero(group_multi == group_value)
    mean = x_multi[idx] @ beta_multi
    covariance = random_var_multi * np.ones((3, 3)) + residual_var_multi * np.eye(3)
    cens_idx = np.array([0, 1])
    obs_idx = np.array([2])
    residual = y_multi[idx][obs_idx] - mean[obs_idx]
    cov_obs = covariance[np.ix_(obs_idx, obs_idx)]
    inv_obs = np.linalg.inv(cov_obs)
    sign, logdet = np.linalg.slogdet(cov_obs)
    exact_multi_loglik += -0.5 * (np.log(2.0 * np.pi) + logdet + residual @ inv_obs @ residual)
    cross = covariance[np.ix_(cens_idx, obs_idx)]
    cond_mean = mean[cens_idx] + (cross @ inv_obs @ residual).ravel()
    cond_cov = covariance[np.ix_(cens_idx, cens_idx)] - cross @ inv_obs @ cross.T
    probability = bivariate_normal_upper_probability(cond_mean, cond_cov, y_multi[idx][cens_idx])
    exact_multi_loglik += math.log(probability)
if abs(got["lmmc_multi_loglik"] - exact_multi_loglik) > 8.0e-4:
    raise SystemExit(
        "parity failure lmmc_multi_loglik: "
        f"Fortran={got['lmmc_multi_loglik']:.16g} exact={exact_multi_loglik:.16g}"
    )

print("Independent NumPy/SciPy/statsmodels parity checks passed.")
