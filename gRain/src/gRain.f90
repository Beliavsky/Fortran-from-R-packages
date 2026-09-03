module grain
  use r_kinds, only : dp
  use grbase_types, only : table_t, rip_order_t
  use grain_types, only : grain_network_t
  use grain_compile, only : grain_ok, grain_invalid_input, grain_inconsistent_cardinality
  use grain_compile, only : grain_not_dag, grain_impossible_evidence, grain_not_compiled
  use grain_compile, only : compile_network, initialize_network_from_potentials, rebuild_potentials
  use grain_cpt, only : make_cpt, logical_and_cpt, logical_or_cpt
  use grain_cpt, only : mendel_cpt, mendel_probability
  use grain_propagation, only : propagate_ls, probability_of_evidence
  use grain_evidence, only : set_hard_evidence, set_soft_evidence, clear_evidence
  use grain_evidence, only : absorb_evidence, has_evidence
  use grain_query, only : query_marginal, query_all_marginals, query_joint
  use grain_query, only : query_conditional, reconstructed_joint
  use grain_simulation, only : simulate_from_uniforms, simulate_network
  use grain_data, only : empirical_table, estimate_cpt, estimate_cpts_from_data
  use grain_data, only : estimate_clique_marginals, marginals_to_potentials
  use grain_data, only : potentials_to_marginals
  use grain_update, only : replace_network_cpt
  implicit none
  private

  public :: dp
  public :: table_t
  public :: rip_order_t
  public :: grain_network_t
  public :: grain_ok
  public :: grain_invalid_input
  public :: grain_inconsistent_cardinality
  public :: grain_not_dag
  public :: grain_impossible_evidence
  public :: grain_not_compiled
  public :: make_cpt
  public :: logical_and_cpt
  public :: logical_or_cpt
  public :: mendel_cpt
  public :: mendel_probability
  public :: compile_network
  public :: initialize_network_from_potentials
  public :: rebuild_potentials
  public :: propagate_ls
  public :: probability_of_evidence
  public :: set_hard_evidence
  public :: set_soft_evidence
  public :: clear_evidence
  public :: absorb_evidence
  public :: has_evidence
  public :: query_marginal
  public :: query_all_marginals
  public :: query_joint
  public :: query_conditional
  public :: reconstructed_joint
  public :: simulate_from_uniforms
  public :: simulate_network
  public :: empirical_table
  public :: estimate_cpt
  public :: estimate_cpts_from_data
  public :: estimate_clique_marginals
  public :: marginals_to_potentials
  public :: potentials_to_marginals
  public :: replace_network_cpt

end module grain
