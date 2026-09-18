from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]


def logistic(z):
    z = np.asarray(z, float)
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


def logit(p):
    return np.log(p / (1.0 - p))


def parse_fortran(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        out[parts[0]] = np.array([float(x) for x in parts[1:]], dtype=float)
    return out


def mixture_cell(probs, share, cell):
    val = 0.0
    for r in range(2):
        pr = share[r]
        for j in range(3):
            pr *= probs[r, cell[j] - 1, j]
        val += pr
    return val


def nocov_data():
    n = 2400
    probs = np.empty((2, 2, 3))
    probs[0, :, 0] = [0.90, 0.10]
    probs[0, :, 1] = [0.80, 0.20]
    probs[0, :, 2] = [0.85, 0.15]
    probs[1, :, 0] = [0.15, 0.85]
    probs[1, :, 1] = [0.25, 0.75]
    probs[1, :, 2] = [0.20, 0.80]
    share = np.array([0.4, 0.6])
    rows = []
    for i in range(8):
        cell = np.array([1 + (i & 1), 1 + ((i >> 1) & 1), 1 + ((i >> 2) & 1)])
        count = int(np.rint(mixture_cell(probs, share, cell) * n))
        rows.extend([cell.copy()] * count)
    rows = rows[:n]
    while len(rows) < n:
        rows.append(np.array([2, 2, 2]))
    return np.array(rows, dtype=int), probs, share


def deterministic_uniform(n, offset):
    state = (123457 + 104729 * int(offset)) % 2147483646 + 1
    u = np.empty(n)
    for i in range(n):
        state = (48271 * state) % 2147483647
        u[i] = state / 2147483647.0
    return u


def rmulti(prob, u):
    cs = np.cumsum(prob, axis=1)
    return 1 + np.sum(u[:, None] > cs[:, :-1], axis=1)


def regression_data():
    n = 1200
    probs = np.empty((2, 2, 3))
    probs[0, :, 0] = [0.90, 0.10]
    probs[0, :, 1] = [0.85, 0.15]
    probs[0, :, 2] = [0.80, 0.20]
    probs[1, :, 0] = [0.10, 0.90]
    probs[1, :, 1] = [0.20, 0.80]
    probs[1, :, 2] = [0.15, 0.85]
    beta = np.array([-0.25, 1.1])
    x = np.column_stack([np.ones(n), np.linspace(-1.5, 1.5, n)])
    p2 = logistic(x @ beta)
    prior = np.column_stack([1.0 - p2, p2])
    cls = rmulti(prior, deterministic_uniform(n, 17))
    y = np.empty((n, 3), dtype=int)
    for j in range(3):
        pmat = probs[cls - 1, :, j]
        y[:, j] = rmulti(pmat, deterministic_uniform(n, 41 + j))
    return y, x, probs, beta


def nocov_objective(theta, y):
    p2 = logistic(theta[0])
    share = np.array([1.0 - p2, p2])
    q = logistic(theta[1:]).reshape(2, 3)
    ll = np.zeros(y.shape[0])
    for r in range(2):
        pr = np.full(y.shape[0], share[r])
        for j in range(3):
            pr *= np.where(y[:, j] == 2, q[r, j], 1.0 - q[r, j])
        ll += pr
    return -np.sum(np.log(ll))


def reg_objective(theta, y, x):
    beta = theta[:2]
    q = logistic(theta[2:]).reshape(2, 3)
    p2 = logistic(x @ beta)
    ll = np.zeros(y.shape[0])
    for r in range(2):
        pr = p2 if r == 1 else 1.0 - p2
        pr = pr.copy()
        for j in range(3):
            pr *= np.where(y[:, j] == 2, q[r, j], 1.0 - q[r, j])
        ll += pr
    return -np.sum(np.log(ll))



def posterior_from_params(y, prior, q):
    n = y.shape[0]
    post = np.empty((n, 2))
    for r in range(2):
        val = prior[:, r].copy()
        for j in range(3):
            val *= np.where(y[:, j] == 2, q[r, j], 1.0 - q[r, j])
        post[:, r] = val
    post /= post.sum(axis=1, keepdims=True)
    return post


def polca_se_reference(y, x, prior, q):
    post = posterior_from_params(y, prior, q)
    score_cols = []
    for r in range(2):
        for j in range(3):
            score_cols.append(post[:, r] * ((y[:, j] == 2).astype(float) - q[r, j]))
    for r in range(1, 2):
        for l in range(x.shape[1]):
            score_cols.append(x[:, l] * (post[:, r] - prior[:, r]))
    score = np.column_stack(score_cols)
    vce = np.linalg.pinv(score.T @ score, rcond=np.sqrt(np.finfo(float).eps))
    qprob = 6
    jac = np.zeros((12, 6))
    rpos = 0
    cpos = 0
    for r in range(2):
        for j in range(3):
            p1 = 1.0 - q[r, j]
            p2 = q[r, j]
            jac[rpos + 0, cpos] = -p1 * p2
            jac[rpos + 1, cpos] = p2 * (1.0 - p2)
            rpos += 2
            cpos += 1
    vprob = jac @ vce[:qprob, :qprob] @ jac.T
    seprob = np.sqrt(np.clip(np.diag(vprob), 0, None)).reshape(2, 3, 2).transpose(0, 2, 1)
    qbeta = x.shape[1]
    vbeta = vce[qprob:qprob+qbeta, qprob:qprob+qbeta]
    jmix = np.zeros((2, qbeta))
    for l in range(qbeta):
        pp = prior[:, 0] * prior[:, 1] * x[:, l]
        jmix[0, l] = -pp.mean()
        jmix[1, l] = pp.mean()
    vmix = jmix @ vbeta @ jmix.T
    return seprob, np.sqrt(np.clip(np.diag(vmix), 0, None)), np.sqrt(np.clip(np.diag(vbeta), 0, None))

def assert_close(name, actual, expected, atol):
    err = np.max(np.abs(np.asarray(actual) - np.asarray(expected)))
    if not np.isfinite(err) or err > atol:
        raise AssertionError(f"{name}: max abs error {err:.6g} exceeds {atol}; actual={actual}; expected={expected}")
    print(f"PASS {name}: max_abs_error={err:.3g}")


def main():
    parity_path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / '.parity' / 'fortran.txt'
    f = parse_fortran(parity_path)

    y, p0, s0 = nocov_data()
    theta0 = np.r_[logit(s0[1]), logit(p0[:, 1, :].ravel())]
    opt = minimize(nocov_objective, theta0, args=(y,), method='BFGS', options={'gtol': 1e-9, 'maxiter': 3000})
    if not opt.success and np.linalg.norm(opt.jac) > 2e-4:
        raise RuntimeError(f"no-cov optimizer failed: {opt.message}; |grad|={np.linalg.norm(opt.jac)}")
    q = logistic(opt.x[1:]).reshape(2, 3)
    share2 = logistic(opt.x[0])
    share = np.array([1.0 - share2, share2])
    assert_close('no-cov loglik', f['nocov_loglik'][0], -opt.fun, 2e-6)
    assert_close('no-cov shares', f['nocov_share'], share, 2e-5)
    assert_close('no-cov response probabilities class 1', f['nocov_p2_c1'], q[0], 3e-5)
    assert_close('no-cov response probabilities class 2', f['nocov_p2_c2'], q[1], 3e-5)
    prior = np.tile(share, (y.shape[0], 1))
    seprob, share_se, _ = polca_se_reference(y, np.ones((y.shape[0], 1)), prior, q)
    assert_close('no-cov class-share SE', f['nocov_share_se'], share_se, 3e-7)
    assert_close('no-cov response SE class 1', f['nocov_p2_se_c1'], seprob[0, 1, :], 3e-7)
    assert_close('no-cov response SE class 2', f['nocov_p2_se_c2'], seprob[1, 1, :], 3e-7)

    y, x, p0, beta0 = regression_data()
    theta0 = np.r_[beta0, logit(p0[:, 1, :].ravel())]
    opt = minimize(reg_objective, theta0, args=(y, x), method='BFGS', options={'gtol': 1e-8, 'maxiter': 5000})
    if not opt.success and np.linalg.norm(opt.jac) > 5e-4:
        raise RuntimeError(f"reg optimizer failed: {opt.message}; |grad|={np.linalg.norm(opt.jac)}")
    q = logistic(opt.x[2:]).reshape(2, 3)
    assert_close('regression loglik', f['reg_loglik'][0], -opt.fun, 2e-5)
    assert_close('regression beta', f['reg_beta'], opt.x[:2], 5e-4)
    assert_close('regression response probabilities class 1', f['reg_p2_c1'], q[0], 7e-4)
    assert_close('regression response probabilities class 2', f['reg_p2_c2'], q[1], 7e-4)
    p2 = logistic(x @ opt.x[:2])
    prior = np.column_stack([1.0 - p2, p2])
    _, share_se, beta_se = polca_se_reference(y, x, prior, q)
    assert_close('regression beta SE', f['reg_beta_se'], beta_se, 3e-6)
    assert_close('regression class-share SE', f['reg_share_se'], share_se, 3e-6)

    print('Independent SciPy parity checks passed.')


if __name__ == '__main__':
    main()
