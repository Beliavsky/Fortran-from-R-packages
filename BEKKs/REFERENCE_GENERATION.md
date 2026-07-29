# Reference generation

The fixed likelihood values in `test/test_core.f90` were calculated with an independent NumPy implementation in:

```text
scripts/generate_references.py
```

Run:

```bash
python scripts/generate_references.py
```

The script independently implements:

- Gaussian multivariate log likelihood;
- initial second moment `data' data / T`;
- full, diagonal, and scalar BEKK recursions;
- the asymmetric sign indicator.

It does not import or call the Fortran library.

The simulation checks use fixed innovation matrices so the covariance and observation recursions can be compared independently of random-number generators.
