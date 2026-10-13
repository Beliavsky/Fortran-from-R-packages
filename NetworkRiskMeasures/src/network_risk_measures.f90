! SPDX-License-Identifier: GPL-3.0-only
! Public facade for the NetworkRiskMeasures 0.1.7 numerical translation.
module network_risk_measures
   use network_risk_base, only: dp, nr_ok, nr_shape, nr_domain, nr_iterations, nr_infeasible, &
      nr_random_exhausted, nr_numeric, risk_matrix, impact_matrix, vulnerability_matrix, criticality
   use network_communicability, only: communicability_matrix, impact_susceptibility, impact_fluidity, impact_diffusion
   use network_contagion, only: contagion, contagion_result, summarize_contagion
   use network_reconstruction, only: max_ent, min_dens, matrix_estimation
   implicit none
   public
end module network_risk_measures
