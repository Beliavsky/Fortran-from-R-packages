# Upstream provenance

- Project: ARPACK-NG
- Repository: <https://github.com/opencollab/arpack-ng>
- Release: `3.9.1`
- Commit: `40329031ae8deb7c1e26baf8353fa384fc37c251`
- License: BSD-3-Clause; see `LICENSE`

The double-precision symmetric and nonsymmetric sources and their required
utilities were mechanically converted to free source form with
`mnormt/tools/fixed_to_free.py`. The numerical algorithms were not translated
or rewritten. `arpack_backend_compat.f90` and `rfortran_arpack.f90` are local
integration files.
