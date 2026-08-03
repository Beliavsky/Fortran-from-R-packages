# Translation coverage

The 33 computational exports in the upstream namespace are represented:

- 9 kernel constructors: `rbfdot`, `laplacedot`, `besseldot`, `polydot`,
  `tanhdot`, `vanilladot`, `anovadot`, `splinedot`, `stringdot`.
- 4 kernel operations: `kernelMatrix`, `kernelMult`, `kernelPol`, `kernelFast`.
- 13 main methods: `kmmd`, `kpca`, `kcca`, `kha`, `specc`, `kkmeans`, `ksvm`,
  `rvm`, `gausspr`, `ranking`, `csi`, `lssvm`, `kqr`.
- 4 utilities: `ipop`, `inchol`, `couple`, `sigest`.
- 3 incremental/feature methods: `onlearn`, `inlearn`, `kfa`.

The non-exported Fourier numeric kernel is also included.

Not translated as compiled APIs: S4 accessors/classes, formula/model-frame
wrappers, plotting/show methods, R cross-validation orchestration, bundled
binary datasets, and package/vignette infrastructure.
