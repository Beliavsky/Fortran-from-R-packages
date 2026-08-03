# Translation coverage

## Exported Fingraph routines

| R routine | Fortran routine | Status |
|---|---|---|
| `learn_connected_graph` | `learn_connected_graph` | Translated |
| `learn_regular_heavytail_graph` | `learn_regular_heavytail_graph` | Translated |
| `learn_kcomp_heavytail_graph` | `learn_kcomp_heavytail_graph` | Translated |

## Non-exported computational routines

| R routine | Fortran routine | Status |
|---|---|---|
| `compute_student_weights` | `compute_student_weight` | Translated |
| `compute_augmented_lagrangian_ht` | `compute_augmented_lagrangian_ht` | Translated |
| `compute_augmented_lagrangian_kcomp_ht` | `compute_augmented_lagrangian_kcomp_ht` | Translated |

The singular Fortran name `compute_student_weight` reflects that the helper
returns one scalar weight for one observation.

## Dependency functionality retained

The required spectralGraphTopology functionality was adapted into the project:

- graph operators and inverses
- graph metrics and block-diagonal construction
- symmetric pseudoinverse and eigendecomposition
- nonnegative quadratic initialization

Portable Gaussian and Student-t generators adapted from fitHeavyTail are
included for reproducible tests and examples.

## Omitted noncomputational material

- progress-bar display
- plotting/network visualization
- R S3/package infrastructure
- RDS data files and external data-download notebooks
- generated HTML, images, and presentation PDFs

Selected original metadata, documentation, and all three computational R files
are retained under `original/`.
