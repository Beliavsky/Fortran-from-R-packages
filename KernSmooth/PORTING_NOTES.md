# Porting notes

The R package combines R wrappers, FFT-based convolution, and legacy fixed-form Fortran kernels. This port exposes a single modern module, `kernsmooth_mod`.

Key choices:

1. Binning follows the same linear interpolation idea as upstream.
2. KDE convolution is direct rather than FFT-based, eliminating an external FFT dependency.
3. Local polynomial regression is solved directly by weighted normal equations with pivoted Gaussian elimination.
4. `dpih` and `dpik` retain the Wand-Jones plug-in formulas; the current implementation covers the common/default plug-in levels and level 0 exactly. Higher-level recursive selectors can be added later if bit-for-bit R parity is required.
5. `dpill` uses a modern direct local-polynomial pilot estimate of curvature and residual variance instead of the upstream blocked quartic/Mallows-Cp implementation.
6. Plotting and R object/interface behavior do not exist in KernSmooth's numerical core and are not ported.
