# API map

The upstream `mcmc` 0.9-8 package exports seven unique computational APIs.

| R API | Fortran mapping |
|---|---|
| `metrop` | `metrop` with a typed log-density callback and `metrop_result` |
| `temper` | `temper_serial` and `temper_parallel` |
| `morph` | `morph_create` returning `morph_transform` |
| `morph.identity` | `morph_identity` |
| `morph.metrop` | `morph_metrop` |
| `initseq` | `initseq` returning `initseq_result` |
| `olbm` | `olbm` |

The R S3 continuation methods (`metrop.metropolis`, `temper.tempering`, and
`morph.metrop.morph.metropolis`) are represented by explicit result objects:
use the previous result's `final` state and options in a new call. This keeps
the numerical state visible and avoids emulating R class dispatch.

## Proposal scaling

R's three Metropolis proposal forms map to:

- scalar: `scale_constant(s)`;
- element-wise vector: `scale_diagonal(v)`;
- matrix: `scale_full(A)`.

Serial and parallel tempering accept an array of `mcmc_scale` objects, one per
component, or a length-one array reused for all components.

## Output callbacks

`metrop` and serial tempering accept the native Fortran `output_callback`
interface. Parallel tempering has `parallel_output_callback`. If omitted,
default batch output matches the continuous state components described by the
upstream package.
