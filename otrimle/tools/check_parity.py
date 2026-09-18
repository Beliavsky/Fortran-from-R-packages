#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2, multivariate_normal


def read_records(path):
    rows = []
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            rows.append((r['key'], int(r['i']), int(r['j']), float(r['value'])))
    return rows


def scalar(rows, key):
    return next(v for k, i, j, v in rows if k == key and i == 0 and j == 0)


def values(rows, key):
    return {(i, j): v for k, i, j, v in rows if k == key}


def bw_nrd0(x):
    x = np.asarray(x, float)
    n = len(x)
    sd = x.std(ddof=1)
    q25, q75 = np.quantile(x, [0.25, 0.75], method='linear')
    scale = min(sd, (q75 - q25) / 1.34)
    if scale <= 0:
        scale = sd
    if scale <= 0:
        scale = abs(x[0])
    if scale <= 0:
        scale = 1.0
    return 0.9 * scale * n ** (-0.2)


def r_density_gaussian(x, weights, lower, upper, n_user, bandwidth):
    n = max(n_user, 512)
    if n > 512:
        n = 2 ** math.ceil(math.log2(n))
    lo = lower - 4 * bandwidth
    up = upper + 4 * bandwidth
    y = np.zeros(2 * n)
    xdelta = (up - lo) / (n - 1)
    for xi, wi in zip(x, weights):
        xpos = (xi - lo) / xdelta
        ix = math.floor(xpos)
        fx = xpos - ix
        if 0 <= ix <= n - 2:
            y[ix] += (1 - fx) * wi
            y[ix + 1] += fx * wi
        elif ix == -1:
            y[0] += fx * wi
        elif ix == n - 1:
            y[ix] += (1 - fx) * wi
    kords = np.linspace(0, 2 * (up - lo), 2 * n)
    kords[n + 1:] = -kords[n - 1:0:-1]
    kernel = np.exp(-0.5 * (kords / bandwidth) ** 2) / (bandwidth * math.sqrt(2 * math.pi))
    conv = np.fft.ifft(np.fft.fft(y) * np.conj(np.fft.fft(kernel))).real[:n]
    conv = np.maximum(0, conv)
    return np.interp(np.linspace(lower, upper, n_user), np.linspace(lo, up, n), conv)


def kd_measure_1d(x, weights, qmax=3.2905267314919255, nk=20):
    w = np.maximum(np.asarray(weights, float), 0)
    w = w / w.sum()
    z = np.asarray(x, float) - np.sum(w * x)
    h = bw_nrd0(z)
    dens = r_density_gaussian(z, w, -qmax, qmax, nk, h)
    ys = np.sort(dens)[::-1]
    cp = 0.5 * (ys[0::2] + ys[1::2])
    m = nk // 2
    return math.sqrt((np.sum((dens[:m][::-1] - cp) ** 2) +
                      np.sum((dens[m:] - cp) ** 2)) / nk)


def scale_tau2_reference(x):
    x = np.asarray(x, float)
    mu0 = np.median(x)
    xabs = np.abs(x - mu0)
    s0 = np.median(xabs)
    if s0 <= 0:
        return mu0, 0.0
    c1 = 4.5
    c2 = 3.0
    q75 = 0.6744897501960817
    b = c2 * q75
    pnorm = 0.5 * (1 + math.erf(b / math.sqrt(2)))
    dnorm = math.exp(-0.5 * b * b) / math.sqrt(2 * math.pi)
    es2 = 2 * ((1 - b * b) * pnorm - b * dnorm + b * b) - 1
    weights = np.maximum(0, 1 - (xabs / (s0 * c1)) ** 2) ** 2
    mu = np.sum(x * weights) / np.sum(weights)
    rho = np.minimum(((x - mu) / s0) ** 2, c2 * c2)
    scale = s0 * math.sqrt(np.sum(rho) / (len(x) * es2))
    return mu, scale


def kmeanfun(n):
    return math.exp(-1.83084) * n ** (-0.55323) if n < 20 else math.exp(-2.4562342) * n ** (-0.3524332)


def ksdfun(n):
    return math.exp(-2.1482546) * n ** (-0.4639833) if n < 12.5 else math.exp(-2.6575388) * n ** (-0.4547216)


def cluster_density_measure(x, pi, means, covs, tau):
    ddpm = []
    for j in range(2):
        vals, vecs = np.linalg.eigh(covs[j])
        scores = (x - means[j]) @ vecs
        cm = []
        nw = tau[:, j + 1].sum()
        for q in range(2):
            z = scores[:, q] / math.sqrt(max(vals[q], np.finfo(float).tiny))
            cm.append(kd_measure_1d(z, tau[:, j + 1]))
        cm = np.asarray(cm)
        st = (cm - kmeanfun(nw)) / ksdfun(nw)
        ddpm.append(np.sum(np.maximum(st, 0) ** 2) / 2)
    ddpm = np.asarray(ddpm)
    return math.sqrt(np.sum(pi[1:] * ddpm ** 2))


def main(path):
    rows = read_records(path)
    code = scalar(rows, 'code')
    if int(round(code)) <= 0:
        raise AssertionError('Fortran fit failed')
    x = np.array([
        [-2.2, -2.0], [-1.8, -2.1], [-2.0, -1.7], [-2.3, -2.2], [-1.7, -1.9],
        [2.0, 2.1], [2.2, 1.8], [1.7, 2.0], [2.3, 2.2], [1.9, 1.7],
        [8.0, -7.0], [-8.0, 7.0],
    ])
    pi_map = values(rows, 'pi')
    pi = np.array([pi_map[(i, 0)] for i in range(1, 4)])
    mean_map = values(rows, 'mean')
    means = np.array([[mean_map[(i, j)] for i in range(1, 3)] for j in range(1, 3)])
    col1 = values(rows, 'cov11row')
    col2 = values(rows, 'cov12row')
    covs = []
    for j in range(1, 3):
        covs.append(np.array([[col1[(1, j)], col2[(1, j)]],
                              [col1[(2, j)], col2[(2, j)]]]))
    covs = np.asarray(covs)
    tau_map = values(rows, 'tau')
    tau = np.array([[tau_map[(i, j)] for j in range(1, 4)] for i in range(1, 13)])
    smd_map = values(rows, 'smd')
    smd = np.array([[smd_map[(i, j)] for j in range(1, 3)] for i in range(1, 13)])

    icd = math.exp(scalar(rows, 'logicd'))
    densities = np.full(len(x), pi[0] * icd)
    for j in range(2):
        densities += pi[j + 1] * multivariate_normal.pdf(x, mean=means[j], cov=covs[j])
    loglik = np.log(densities).sum()
    reported_ll = scalar(rows, 'iloglik')
    if abs(loglik - reported_ll) > 2e-4:
        raise AssertionError(f'log-likelihood mismatch: python={loglik} fortran={reported_ll}')

    post = np.empty_like(tau)
    post[:, 0] = pi[0] * icd
    for j in range(2):
        post[:, j + 1] = pi[j + 1] * multivariate_normal.pdf(x, mean=means[j], cov=covs[j])
    post /= post.sum(axis=1, keepdims=True)
    if np.max(np.abs(post - tau)) > 3e-4:
        raise AssertionError(f'posterior mismatch: max abs {np.max(np.abs(post-tau))}')

    kd = []
    for j in range(2):
        order = np.argsort(smd[:, j])
        c = np.cumsum(tau[order, j + 1]) / tau[:, j + 1].sum()
        kd.append(np.max(np.abs(c - chi2.cdf(smd[order, j], df=2))))
    criterion = np.dot(kd, pi[1:]) / pi[1:].sum()
    reported_criterion = scalar(rows, 'criterion')
    if abs(criterion - reported_criterion) > 2e-7:
        raise AssertionError(f'criterion mismatch: python={criterion} fortran={reported_criterion}')

    eig = np.concatenate([np.linalg.eigvalsh(c) for c in covs])
    if eig.max() / eig.min() > 20.000001:
        raise AssertionError('eigenratio constraint violated')
    if tau[:, 0].mean() > 0.300001:
        raise AssertionError('noise posterior constraint violated')

    kd_total = cluster_density_measure(x, pi, means, covs, tau)
    reported_kd = scalar(rows, 'kdmeasure')
    if abs(kd_total - reported_kd) > 5e-10:
        raise AssertionError(f'kernel-density diagnostic mismatch: python={kd_total} fortran={reported_kd}')

    tau_mu, tau_sd = scale_tau2_reference([-4.0, -1.5, -0.3, 0.0, 0.2, 0.7, 1.1, 2.0, 8.0])
    if abs(tau_mu - scalar(rows, 'tau_location')) > 3e-14:
        raise AssertionError('scaleTau2 robust-location mismatch')
    if abs(tau_sd - scalar(rows, 'tau_scale')) > 3e-14:
        raise AssertionError('scaleTau2 robust-scale mismatch')

    print('Independent NumPy/SciPy parity checks passed.')
    print(f'  log likelihood difference: {abs(loglik-reported_ll):.3e}')
    print(f'  posterior max difference:  {np.max(np.abs(post-tau)):.3e}')
    print(f'  criterion difference:      {abs(criterion-reported_criterion):.3e}')
    print(f'  density measure difference:{abs(kd_total-reported_kd):.3e}')
    print(f'  scaleTau2 location diff:   {abs(tau_mu-scalar(rows, "tau_location")):.3e}')
    print(f'  scaleTau2 scale diff:      {abs(tau_sd-scalar(rows, "tau_scale")):.3e}')


if __name__ == '__main__':
    main(Path(sys.argv[1] if len(sys.argv) > 1 else 'fortran_parity.csv'))
