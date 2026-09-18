#!/usr/bin/env python3
import csv
import math
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import multivariate_normal

path = sys.argv[1] if len(sys.argv) > 1 else 'fortran_parity.csv'
rows = defaultdict(list)
with open(path, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        rows[row['name']].append((int(row['i']), int(row['j']), float(row['value'])))

def scalar(name):
    return rows[name][0][2]

def close(actual, expected, tol, name):
    if not math.isfinite(actual) or abs(actual - expected) > tol:
        raise SystemExit(f'{name}: got {actual:.17g}, expected {expected:.17g}, tol {tol:g}')

# Independent agreement-index formulas.
T = np.array([[1, 1, 0], [1, 2, 1], [0, 0, 4]], dtype=float)
n = T.sum()
r = T.sum(axis=1)
c = T.sum(axis=0)
nis = np.sum(r*r)
njs = np.sum(c*c)
tot = n*(n-1)/2
t2 = np.sum(T*T)
t3 = (nis+njs)/2
nc = (n*(n*n+1) - (n+1)*nis - (n+1)*njs + 2*nis*njs/n)/(2*(n-1))
A = tot+t2-t3
D = -t2+t3
close(scalar('rand_adjusted'), (A-nc)/(tot-nc), 2e-13, 'adjusted Rand')
close(scalar('rand_raw'), A/tot, 2e-13, 'Rand')
Tk = np.sum(T*T)-n
Pk = np.sum(r*r)-n
Qk = np.sum(c*c)-n
Bk = Tk/np.sqrt(Pk*Qk)
EBk = np.sqrt(Pk*Qk)/(n*(n-1))
close(scalar('fm_raw'), Bk, 2e-13, 'Fowlkes-Mallows')
close(scalar('fm_adjusted'), (Bk-EBk)/(1-EBk), 2e-13, 'adjusted Fowlkes-Mallows')

# Recreate deterministic clustering data.
x = np.empty((30, 2), float)
for i in range(1, 13):
    x[i-1] = [-5 + .15*(i-6), -4 + .10*((i % 4)-2)]
for i in range(13, 25):
    x[i-1] = [5 + .12*(i-18), 4 + .08*((i % 5)-2)]
for i in range(25, 31):
    x[i-1] = [30+i, -30-2*i]

labels = np.array([round(v) for _, _, v in sorted(rows['tclust_label'])], dtype=int)
weights = {i: v for i, _, v in rows['tclust_weight']}
centers = np.zeros((2,2))
for i,j,v in rows['tclust_center']:
    centers[i-1,j-1] = v
covs = np.zeros((2,2,2))
for i,j,v in rows['tclust_cov']:
    covs[i-1,i-1,j-1] = v
for _,j,v in rows['tclust_cov12']:
    covs[0,1,j-1] = covs[1,0,j-1] = v
for j in (1,2):
    assigned = x[labels == j]
    close(centers[0,j-1], assigned[:,0].mean(), 2e-12, f'tclust center x {j}')
    close(centers[1,j-1], assigned[:,1].mean(), 2e-12, f'tclust center y {j}')
# Objective independently evaluated from returned parameters and classifications.
obj = 0.0
for i,lab in enumerate(labels):
    if lab > 0:
        obj += math.log(weights[lab]) + multivariate_normal.logpdf(x[i], mean=centers[:,lab-1], cov=covs[:,:,lab-1])
close(scalar('tclust_obj'), obj, 2e-9, 'tclust HARD objective')
eigs = np.concatenate([np.linalg.eigvalsh(covs[:,:,j]) for j in range(2)])
if eigs.max()/eigs.min() > 20*(1+2e-10):
    raise SystemExit('tclust covariance eigenvalue restriction violated')

# Trimmed k-means objective from the returned centers and labels.
kmlab = np.array([round(v) for _, _, v in sorted(rows['tkmeans_label'])], dtype=int)
kmcent = np.zeros((2,2))
for i,j,v in rows['tkmeans_center']:
    kmcent[i-1,j-1] = v
ss = 0.0
for i,lab in enumerate(kmlab):
    if lab > 0:
        ss += np.sum((x[i]-kmcent[:,lab-1])**2)
close(scalar('tkmeans_obj'), ss/np.sum(kmlab > 0), 2e-12, 'trimmed k-means objective')

# Robust linear grouping: independently fit rank-1 PCA lines to each returned group.
y = np.empty((36,2), float)
for i in range(1,16):
    y[i-1] = [-4+.4*i, -2+.8*i]
for i in range(16,31):
    y[i-1] = [-5+.45*(i-15), 8-.65*(i-15)]
for i in range(31,37):
    y[i-1] = [30+i, -20-i]
rlab = np.array([round(v) for _,_,v in sorted(rows['rlg_label'])], dtype=int)
rss = 0.0
for lab in (1,2):
    z = y[rlab == lab]
    mu = z.mean(axis=0)
    _, _, vh = np.linalg.svd(z-mu, full_matrices=False)
    u = vh[0]
    resid = z-mu - np.outer((z-mu)@u, u)
    rss += np.sum(resid*resid)
close(scalar('rlg_obj'), rss, 2e-9, 'rlg affine residual objective')

print('Independent NumPy/SciPy parity checks passed.')
