# Actual R reference outputs

These text files were produced by the official installed CRAN CreditRisk
0.1.7, not by evaluating the translation. `R-reference-session.txt` records
the R runtime/platform and package version. From the package root, with
CreditRisk 0.1.7 installed, regenerate using:

```sh
Rscript reference/generate_reference.R reference
python3 reference/embed_reference.py
fpm test
```

All four calibration fixtures use real `stats::optim` minimizations at
nontrivial synthetic market quotes. BlackCox records its unmodified public
return separately from the optimizer's actual final coordinates because the
upstream wrapper returns starting L/K. A temporary read-only R trace captures
that optimizer return; the source and objective are unchanged. Simulation
fixtures include actual R output and the corresponding standardized normal
innovations so Fortran path calculations can be tested independently of RNG
choice. No R or Python is needed for ordinary FPM tests.

The R generation script deliberately uses cases on which the upstream
optimizer completes. Some other finite quote sets drive R BlackCox to
invalid K/negative survival during unconstrained trial evaluations and abort
with a nonfinite-objective error. The port assigns those trials a finite
penalty; it does not claim identical optimizer behavior in that situation.
