"""Compare the deterministic Fortran validation case with scikit-learn.

Labels are compared as partitions because cluster numbers are arbitrary.  The
HDBSCAN membership probability is checked against the R dbscan package's
core-distance formula, not sklearn's different probability definition.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from sklearn.cluster import DBSCAN, HDBSCAN, OPTICS
from sklearn.neighbors import LocalOutlierFactor, NearestNeighbors

X = np.array([
    [0.00, 0.00], [0.10, 0.00], [0.00, 0.12], [0.12, 0.10], [0.20, 0.05],
    [3.00, 3.00], [3.10, 3.00], [3.00, 3.12], [3.12, 3.10], [3.20, 3.05],
    [1.50, 1.50], [6.00, 0.00],
])


def same_partition(a: np.ndarray, b: np.ndarray) -> bool:
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False
    # Fortran uses 0 for noise; sklearn uses -1.
    anoise = a == 0
    bnoise = b == -1
    if not np.array_equal(anoise, bnoise):
        return False
    for i in range(len(a)):
        for j in range(len(a)):
            if (a[i] == a[j]) != (b[i] == b[j]):
                return False
    return True


def load(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    return {
        "dbscan": np.array([int(r["dbscan"]) for r in rows]),
        "lof": np.array([float(r["lof"]) for r in rows]),
        "hdbscan": np.array([int(r["hdbscan"]) for r in rows]),
        "hprob": np.array([float(r["hprob"]) for r in rows]),
        "reach": np.array([float(r["reach"]) for r in rows]),
        "core": np.array([float(r["core"]) for r in rows]),
        "predecessor": np.array([int(r["predecessor"]) for r in rows]),
    }


def main() -> None:
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "fortran_parity.csv")
    got = load(path)

    db = DBSCAN(eps=0.25, min_samples=3).fit_predict(X)
    assert same_partition(got["dbscan"], db)

    lof = -LocalOutlierFactor(n_neighbors=2).fit_predict(X)  # force fitting
    del lof
    lof_model = LocalOutlierFactor(n_neighbors=2).fit(X)
    expected_lof = -lof_model.negative_outlier_factor_
    np.testing.assert_allclose(got["lof"], expected_lof, rtol=2e-9, atol=2e-9)

    hd = HDBSCAN(min_cluster_size=3, min_samples=3, copy=True).fit(X)
    assert same_partition(got["hdbscan"], hd.labels_)

    core = NearestNeighbors(n_neighbors=3).fit(X).kneighbors(X)[0][:, 2]
    expected_prob = np.zeros(len(X))
    for label in np.unique(got["hdbscan"]):
        if label == 0:
            continue
        mask = got["hdbscan"] == label
        max_core = core[mask].max()
        expected_prob[mask] = (max_core - core[mask]) / max_core if max_core > 0 else 1.0
    np.testing.assert_allclose(got["hprob"], expected_prob, rtol=2e-10, atol=2e-10)

    op = OPTICS(min_samples=3, max_eps=0.5).fit(X)
    np.testing.assert_allclose(got["reach"], op.reachability_, rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(got["core"], op.core_distances_, rtol=2e-10, atol=2e-10)
    # sklearn and dbscan use different priority tie-breaking, so ordering is not asserted.

    print("dbscan/scikit-learn parity checks passed.")


if __name__ == "__main__":
    main()
