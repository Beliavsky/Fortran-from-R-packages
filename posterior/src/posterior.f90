! SPDX-License-Identifier: BSD-3-Clause
module posterior
   !! Public aggregation module; arrays store iterations/draws on their first axis.
   use posterior_diagnostics
   use posterior_draws
   use posterior_pareto
   use posterior_predictive
   use posterior_transforms
   implicit none
   public
end module posterior
