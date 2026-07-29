# Porting notes

## Translated surface

The original package exports four computational functions. Their mappings are:

| R | Fortran |
|---|---|
| `varcast()` | `varcast()` |
| `trafftest()` | `trafftest()` |
| `covtest()` | `covtest()` |
| `lossfunc()` | `lossfunc()` |

The hidden R helpers `arfilt`, `quant`, and `dens` are translated as explicit numerical procedures.

## Model integration

The supplied Fortran `rugarch` modules are used for sGARCH, eGARCH, APARCH, and FIGARCH. The fit is performed on the in-sample standardized series only. A fixed in-sample variance backcast is then used while filtering the combined series, allowing actual out-of-sample returns to update each subsequent one-step forecast without allowing the out-of-sample observations to affect parameter estimation or initialization.

Log-GARCH uses the translated `smoots` conditional-Gaussian ARMA estimator. FI-Log-GARCH uses the translated Haslett-Raftery `fracdiff` estimator.

## Semiparametric scale estimation

For sGARCH, eGARCH, APARCH, and Log-GARCH, the port calls the translated `smoots::msmooth` computational core with the ufRisk default polynomial order of 3 and default algorithm `A`.

For FIGARCH and FI-Log-GARCH, ufRisk calls `esemifar::tsmoothlm`. The required workflow is embedded here rather than introducing an external R dependency:

1. local-polynomial trend smoothing of log-squared centered returns;
2. ARFIMA candidate estimation with `fracdiff`;
3. BIC order selection;
4. zero-frequency spectrum calculation;
5. fractional-memory bandwidth inflation;
6. derivative-based plug-in bandwidth update;
7. final scale normalization.

The Fortran implementation uses the equivalent-kernel constants already validated in the translated `smoots` library. Consequently, long-memory bandwidths can differ slightly from those produced by private lookup constants in a particular `esemifar` release, although the model structure and iterative algorithm are retained.

## Intentional corrections

### Christoffersen transition likelihood

The R source uses `p01^n10` in one branch of the Markov likelihood. The transition count attached to `p01` is `n01`; the Fortran implementation uses the conventional and documented `p01^n01` expression.

### Student-t fitting for semiparametric Log-GARCH

The R Log-GARCH branch divides already scale-standardized returns by a total volatility that includes the scale a second time when estimating Student-t degrees of freedom. The Fortran implementation estimates degrees of freedom from centered returns divided by total volatility once.

### Student-t upper bound

The supplied rugarch parameter transformation permits extremely large degrees of freedom. ufRisk only needs the normal limit, so the retained forecast result caps an estimated Student-t degree of freedom at 100. This avoids cancellation in gamma-function differences while being numerically indistinguishable from a normal innovation for the risk levels used here.

### Zero returns

The R implementation evaluates `log(return^2)` directly. The Fortran port floors squared returns at `tiny(1.0_dp)` before taking logarithms, preventing negative infinity and preserving a finite optimization problem.

## Numerical differences

- The supplied rugarch translation uses finite-difference BFGS and may not reproduce R's optimizer path exactly.
- The translated smoots ARMA estimator uses conditional Gaussian least squares rather than R's Kalman-likelihood `arima` implementation.
- RNG streams are irrelevant to the principal ufRisk calculations and are not intended to reproduce R streams.
- Numerical tolerances and optimizer iteration limits are explicit options.

## Omitted R infrastructure

- S3 print and plot methods
- startup messages
- R lists, data frames, and `ts` attributes
- bundled WMT and EURO STOXX 50 datasets
- graphics

These omissions do not remove a numerical algorithm exported by ufRisk.
