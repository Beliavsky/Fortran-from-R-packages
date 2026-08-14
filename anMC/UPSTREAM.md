# Upstream provenance

Translation source:

- R package: `anMC`
- version: `0.2.5`
- author/copyright holder listed by upstream: Dario Azzimonti
- upstream license declaration: `GPL-3`
- source archive SHA-256:
  `f6e3a2532f3c37400a104afe249e762b4e07c1b2f4a4c6988d8351624e9c8a1b`

The supplied upstream `DESCRIPTION`, `NAMESPACE`, and `NEWS.md` are retained in
`upstream/`.

A previously supplied `mvtnorm-fortran` archive was inspected during this port:

- archive SHA-256:
  `8f84333b086f2525f226c881ab2a2101546b38d19588930a8cc766a507c5f9f1`
- declared port license: `GPL-2.0-only`.

It was used only as an independent numerical reference during validation and is
not included or linked.  The Gaussian rectangle integration in this release is
a native implementation in the GPL-3 anMC source tree.

## References retained from upstream

D. Azzimonti and D. Ginsbourger (2018), "Estimating orthant probabilities of
high dimensional Gaussian vectors with an application to set estimation",
Journal of Computational and Graphical Statistics 27(2), 255-267.

The upstream source also cites work by Genz on multivariate-normal probability
computation, Dickmann and Schweizer on nested conditional Monte Carlo, and
Chevalier's thesis on Gaussian-process uncertainty reduction.
