#!/usr/bin/env python3
from __future__ import annotations
import csv, math, shutil, statistics, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / "comparisons"
SUITES = (
    ("rugarch", "rugarch"),
    ("fracdiff", "fracdiff"),
    ("performanceanalytics", "PerformanceAnalytics"),
    ("FinTS", "FinTS"),
    ("corpcor", "corpcor"),
    ("cluster", "cluster"),
    ("pracma", "pracma"),
    ("boot", "boot"),
)

def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)

def read_rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return {row["case"]: row for row in csv.DictReader(f)}

def timing_stats(values):
    values = [x for x in values if x is not None and math.isfinite(x) and x > 0]
    if not values:
        return None
    return {
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "geomean": math.exp(statistics.fmean(math.log(x) for x in values)),
        "min": min(values),
        "max": max(values),
    }

def print_timing_summary(results):
    print("\nTiming summaries (positive finite measurements only)")
    scopes = [("overall", results)]
    scopes += [(name, [r for r in results if r["package"] == name])
               for name in sorted({r["package"] for r in results})]
    for name, rows in scopes:
        valid = [r for r in rows if r["speedup"] is not None and
                 math.isfinite(r["speedup"]) and r["speedup"] > 0]
        fortran_wins = sum(r["speedup"] > 1.05 for r in valid)
        r_wins = sum(r["speedup"] < 0.95 for r in valid)
        ties = len(valid) - fortran_wins - r_wins
        unresolved = len(rows) - len(valid)
        print(f"\n{name}: cases={len(rows)}, valid={len(valid)}, "
              f"Fortran faster={fortran_wins}, R faster={r_wins}, "
              f"ties={ties}, unresolved={unresolved}")
        print("  measure       mean       median     geomean    min        max")
        measures = (
            ("R seconds", [r["r_seconds"] for r in valid]),
            ("F seconds", [r["fortran_seconds"] for r in valid]),
            ("R/F", [r["speedup"] for r in valid]),
        )
        for label, values in measures:
            stats = timing_stats(values)
            if stats is None:
                print(f"  {label:<12} N/A")
            else:
                print(f"  {label:<12} {stats['mean']:<10.4g} {stats['median']:<10.4g} "
                      f"{stats['geomean']:<10.4g} {stats['min']:<10.4g} {stats['max']:<10.4g}")
        if valid:
            fastest_f = max(valid, key=lambda r: r["speedup"])
            fastest_r = min(valid, key=lambda r: r["speedup"])
            print(f"  largest R/F: {fastest_f['package']}/{fastest_f['case']} "
                  f"({fastest_f['speedup']:.2f})")
            print(f"  smallest R/F: {fastest_r['package']}/{fastest_r['case']} "
                  f"({fastest_r['speedup']:.2f})")

def main():
    if not shutil.which("Rscript") or not shutil.which("fpm"):
        print("ERROR: Rscript and fpm must be on PATH.", file=sys.stderr); return 2
    results, failures = [], 0
    for name, r_package in SUITES:
        suite = BASE / name
        p = run(["Rscript", "-e", f"quit(status=!requireNamespace('{r_package}',quietly=TRUE))"], ROOT)
        if p.returncode:
            print(f"SKIP {name}: R package is not installed"); continue
        r_csv, f_csv = suite/"r_results.csv", suite/"fortran_results.csv"
        rr = run(["Rscript", "reference.R", str(r_csv)], suite)
        if rr.returncode:
            print(f"FAIL {name} R:\n{rr.stdout}{rr.stderr}", file=sys.stderr); failures += 1; continue
        fr = run(["fpm", "run", "--profile", "release", "--", str(f_csv)], suite)
        if fr.returncode:
            print(f"FAIL {name} Fortran:\n{fr.stdout}{fr.stderr}", file=sys.stderr); failures += 1; continue
        rrows, frows = read_rows(r_csv), read_rows(f_csv)
        for case, rv in rrows.items():
            if case not in frows:
                print(f"FAIL {name}/{case}: missing Fortran result"); failures += 1; continue
            fv=frows[case]; a,b=float(rv["value"]),float(fv["value"])
            atol,rtol=float(rv["abs_tol"]),float(rv["rel_tol"]); err=abs(a-b)
            passed=math.isfinite(a) and math.isfinite(b) and err <= atol+rtol*abs(a)
            rt,ft=float(rv["seconds"]),float(fv["seconds"])
            speedup=rt/ft if rt>0 and ft>0 else None
            status="PASS" if passed else "FAIL"
            ratio=f"{speedup:.2f}" if speedup is not None else "N/A"
            print(f"{status:4} {name}/{case}: error={err:.3g}, R={rt:.4g}s, F={ft:.4g}s, R/F={ratio}")
            results.append(dict(package=name,case=case,status=status,r_value=a,fortran_value=b,
                abs_error=err,abs_tol=atol,rel_tol=rtol,r_seconds=rt,fortran_seconds=ft,speedup=speedup))
            failures += not passed
    out=BASE/"results.csv"
    if results:
        with out.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(results[0])); w.writeheader(); w.writerows(results)
    passed=sum(x["status"]=="PASS" for x in results)
    print(f"\n{passed}/{len(results)} comparisons passed. Results: {out}")
    if results:
        print_timing_summary(results)
    return 1 if failures else 0

if __name__ == "__main__": raise SystemExit(main())
