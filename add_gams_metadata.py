#!/usr/bin/env python3
"""Add conservative top-level NIST GAMS classifications to FPM manifests.

The classifications are inferred from the existing package name, description,
categories, and keywords.  Run with ``--check`` to verify that every manifest
has the expected metadata without changing files.
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parent
GAMS_ORDER = "ABCDEFGHIJKLMNOPQRSZ"

# These are intentionally top-level classifications.  The translated packages
# often expose several algorithms, so assigning a narrow leaf code solely from
# upstream package metadata would imply more precision than we have verified.
RULES: dict[str, tuple[str, ...]] = {
    "A": (
        "error analysis", "floating point", "arbitrary precision",
        "interval arithmetic",
    ),
    "B": ("number theory", "combinatorics", "combinatorial enumeration"),
    "C": (
        "special functions", "polynomial", "polynomials", "hypergeometric", "bessel",
        "gamma function", "beta function",
    ),
    "D": (
        "linear algebra", "matrix computations", "sparse matrix", "sparse array",
        "singular value", "svd", "eigenvalue", "eigendecomposition", "cholesky",
        "principal component", "matrix functions", "pca",
    ),
    "E": ("interpolation", "spline", "splines", "kriging"),
    "F": (
        "nonlinear equations", "root finding", "root solver", "zero finding",
        "steady state",
    ),
    "G": (
        "optimization", "optimisation", "linear programming",
        "quadratic programming", "nonlinear programming", "integer programming",
        "conic programming", "least squares", "genetic algorithm",
        "differential evolution", "particle swarm", "simulated annealing",
        "feature selection", "subset selection", "bfgs", "lbfgs",
    ),
    "H": (
        "numerical integration", "quadrature", "numerical differentiation",
        "automatic differentiation", "finite difference",
    ),
    "I": (
        "differential equations", "differential equation", "integral equations",
        "integral equation", "ordinary differential", "partial differential",
        "differential algebraic", "ode", "dae",
    ),
    "J": (
        "integral transforms", "fourier", "fft", "wavelet", "wavelets", "hilbert transform",
        "laplace transform",
    ),
    "K": ("approximation theory", "function approximation", "approximation"),
    "L": (
        "statistics", "statistical", "probability", "distribution", "distributions", "regression",
        "time series", "econometrics", "machine learning", "bayesian", "clustering",
        "classification", "density estimation", "survival analysis", "bootstrap",
        "finance", "financial", "portfolio", "risk", "actuarial", "insurance",
        "copula", "garch", "volatility", "forecasting", "maximum likelihood",
        "estimation", "covariance", "correlation", "moments", "quantile",
        "goodness of fit", "hypothesis test", "kalman", "state space",
        "random number", "rng", "mcmc",
    ),
    "M": (
        "simulation", "monte carlo", "stochastic modeling", "stochastic modelling",
        "random number", "rng", "mcmc", "markov chain", "bootstrap",
        "resampling",
    ),
    "N": ("data handling", "data manipulation", "data structure", "file format"),
    "O": ("symbolic computation", "computer algebra", "symbolic algebra"),
    "P": (
        "computational geometry", "nearest neighbor", "spatial geometry",
        "convex hull", "delaunay", "voronoi",
    ),
    "Q": ("graphics", "visualization", "visualisation", "plotting"),
    "R": ("service routines", "utilities", "utility routines"),
    "S": (
        "software development", "code generation", "compiler", "testing framework",
    ),
}

# Some older manifests intentionally have generic tags.  These small overrides
# record classifications that are clear from the translated public API but
# cannot be recovered from that metadata alone.
PROJECT_OVERRIDES: dict[str, tuple[str, ...]] = {
    "contfrac": ("K",),
    "expm": ("D",),
    "nilde": ("B", "G"),
    "orthopolynom": ("C", "K"),
    "PSDistr": ("L", "M"),
    "quadform": ("D",),
    "splines": ("E", "K", "L"),
}


def normalize(value: str) -> str:
    """Normalize prose and tags for phrase matching."""
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def phrase_present(text: str, phrase: str) -> bool:
    return re.search(rf"(?:^| )({re.escape(normalize(phrase))})(?: |$)", text) is not None


def load_manifest(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def searchable_text(data: dict[str, object]) -> str:
    values: list[str] = []
    for field in ("name", "description"):
        value = data.get(field)
        if isinstance(value, str):
            values.append(value)
    for field in ("categories", "keywords"):
        value = data.get(field)
        if isinstance(value, list):
            values.extend(item for item in value if isinstance(item, str))
    return normalize(" ".join(values))


def classify(data: dict[str, object], project: str) -> list[str]:
    text = searchable_text(data)
    codes = {
        code
        for code in GAMS_ORDER
        if code in RULES and any(phrase_present(text, phrase) for phrase in RULES[code])
    }
    codes.update(PROJECT_OVERRIDES.get(project, ()))
    return [code for code in GAMS_ORDER if code in codes] or ["Z"]


def format_gams(codes: list[str]) -> str:
    return "gams = [" + ", ".join(f'\"{code}\"' for code in codes) + "]"


def update_text(text: str, codes: list[str]) -> str:
    newline = "\r\n" if "\r\n" in text else "\n"
    replacement = format_gams(codes)
    if re.search(r"(?m)^gams\s*=.*$", text):
        return re.sub(r"(?m)^gams\s*=.*$", replacement, text, count=1)

    lines = text.splitlines(keepends=True)
    insert_after = -1
    for index, line in enumerate(lines):
        if re.match(r"^(categories|keywords)\s*=", line):
            insert_after = index
            if line.startswith("categories"):
                break
    if insert_after < 0:
        for index, line in enumerate(lines):
            if re.match(r"^description\s*=", line):
                insert_after = index
                break

    insertion = replacement + newline
    lines.insert(insert_after + 1, insertion)
    return "".join(lines)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="report missing or outdated GAMS metadata without modifying files",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    manifests = sorted(REPOSITORY_ROOT.glob("*/fpm.toml"))
    changed: list[tuple[str, list[str]]] = []
    counts = {code: 0 for code in GAMS_ORDER}

    for manifest in manifests:
        try:
            data = load_manifest(manifest)
        except (OSError, tomllib.TOMLDecodeError) as exc:
            print(f"ERROR: {manifest}: {exc}", file=sys.stderr)
            return 1
        codes = classify(data, manifest.parent.name)
        for code in codes:
            counts[code] += 1

        with manifest.open("r", encoding="utf-8", newline="") as stream:
            current = stream.read()
        expected = update_text(current, codes)
        if current != expected:
            changed.append((manifest.parent.name, codes))
            if not args.check:
                with manifest.open("w", encoding="utf-8", newline="") as stream:
                    stream.write(expected)

    action = "need updates" if args.check else "updated"
    print(f"{len(changed)} of {len(manifests)} manifests {action}.")
    used = ", ".join(f"{code}={count}" for code, count in counts.items() if count)
    print(f"Assignments: {used}")
    if changed and args.check:
        for name, codes in changed:
            print(f"  {name}: {', '.join(codes)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
