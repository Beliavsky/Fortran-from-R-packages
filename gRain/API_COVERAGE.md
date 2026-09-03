# API and computational coverage

This document maps major upstream gRain 1.4.6 computational areas to the native
Fortran API. It is intentionally explicit about what is and is not parity.

| Upstream area | Fortran API / status | Notes |
| --- | --- | --- |
| CPT creation and normalization | `make_cpt` | Child variable is first; normalized over first dimension. |
| Logical AND/OR CPT helpers | `logical_and_cpt`, `logical_or_cpt` | Binary deterministic CPTs. |
| Mendelian segregation helper | `mendel_cpt`, `mendel_probability` | Unordered diploid genotype states for an arbitrary allele count. |
| DAG construction/validation from CPTs | `compile_network` | Parent-to-child DAG inferred from CPT domains; integer nodes replace R names. |
| Moralization | `compile_network` via gRbase | Reuses sibling gRbase graph kernel. |
| Root-node completion | optional `root=` in `compile_network` | Completes the requested root set before triangulation. |
| Weighted triangulation | `compile_network` via gRbase | Uses node cardinalities as elimination weights. |
| RIP/junction-tree construction | `compile_network` via gRbase | Reuses gRbase running-intersection representation. |
| Clique-potential construction/CPT insertion | `compile_network`, `rebuild_potentials` | Unity clique tables followed by CPT multiplication into host cliques. |
| Initialize from clique potentials | `initialize_network_from_potentials` | For already constructed RIP/potential systems. |
| Lauritzen-Spiegelhalter propagation | `propagate_ls` | Collect/distribute translation of upstream C++ algorithm. |
| Probability of evidence | `probability_of_evidence` | Normalizing constant from propagated potentials. |
| Hard evidence | `set_hard_evidence` | Optional overwrite and propagation controls. |
| Soft/likelihood evidence | `set_soft_evidence` | Nonnegative state weights. |
| Evidence retraction | `clear_evidence` | One node or all nodes. |
| Evidence absorption | `absorb_evidence` | Makes current evidence-updated potentials the new baseline and clears evidence markers. |
| Evidence presence test | `has_evidence` | Node-specific or whole-network state. |
| Single-node marginal query | `query_marginal` | Calibrated posterior marginal. |
| All one-node marginals | `query_all_marginals` | Returns one table per node. |
| Joint query | `query_joint` | Direct clique marginal when possible; otherwise reconstructs the calibrated global joint. |
| Conditional query | `query_conditional` | First requested node is normalized conditionally over remaining requested nodes. |
| Full calibrated joint reconstruction | `reconstructed_joint` | Product of clique marginals divided by separator marginals. Intended for moderate state spaces. |
| Forward simulation with supplied randomness | `simulate_from_uniforms` | Uses calibrated clique tree and supplied `U(0,1)` values. |
| Forward simulation with portable RNG | `simulate_network` | Park-Miller stream; not R-RNG compatible. |
| Empirical categorical table | `empirical_table` | Integer data, optional smoothing/normalization. |
| CPT estimate from data | `estimate_cpt` | Frequency estimate with optional additive smoothing. |
| All DAG CPT estimates | `estimate_cpts_from_data` | Uses supplied DAG and cardinalities. |
| Clique marginals from data | `estimate_clique_marginals` | Frequency tables for RIP cliques. |
| Clique marginals -> potentials | `marginals_to_potentials` | Junction-tree factor conversion. |
| Potentials -> calibrated clique marginals | `potentials_to_marginals` | Propagates a temporary network then returns marginals. |
| Replace CPT values | `replace_network_cpt` | Does not retriangulate when the CPT domain/graph structure is unchanged. |
| R formula/name parsing | Not translated | R-specific interface layer. |
| S3/S4 methods, print/summary, broom/predict wrappers | Not translated | R-specific object/UI layer. |
| Factor labels/dimnames | Adapted | Integer node/state IDs replace R names; dimensions and column-major value ordering are retained. |
| Hugin `.net` parser/writer | Not translated | File-format/interface functionality rather than numerical kernel. |
| Plotting/graph display | Not translated | Explicitly excluded interactive/plotting functionality. |
| R RNG bit-for-bit behavior | Not translated | `simulate_from_uniforms` allows an external RNG; built-in seeded stream is portable but different. |

## Public types

### `grain_network_t`

Stores the DAG, triangulated undirected graph, node cardinalities, RIP/junction
tree, original/current/equilibrated clique potentials, CPTs, evidence state,
and the current evidence probability.

### `table_t` and `rip_order_t`

These gRbase types are re-exported from the public `grain` module so callers can
construct or inspect compatible probability tables and running-intersection
orders without duplicating type definitions.

## Numerical parity boundary

The translated package covers the main discrete Bayesian-network computation:
network compilation to a junction tree, exact probability propagation,
evidence, queries, simulation, and categorical maximum-likelihood/frequency
estimation. It does not reproduce R object semantics or convenience APIs.

The largest deliberate computational limitation is explicit full-joint
reconstruction for joint queries whose variables do not coexist in a single
clique. This is exact, but its memory use grows with the product of all node
cardinalities. Upstream gRain similarly relies on clique-tree structure to avoid
forming the full joint for ordinary local queries. Callers should prefer clique-
local or low-dimensional queries for large networks.
