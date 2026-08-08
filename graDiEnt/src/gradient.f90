! SPDX-License-Identifier: MIT
module gradient
  use gradient_kinds, only : dp
  use gradient_rng, only : seed_rng
  use gradient_types, only : sqgde_options, sqgde_result, objective_fn, get_algo_params, &
                             validate_options, SQGDE_RAND, SQGDE_CURRENT, SQGDE_BEST, &
                             CONVERGE_STDEV, CONVERGE_PERCENT
  use gradient_sqgde, only : optim_sqgde, adapt_sqgde_particle, purify_population
  implicit none
  public
end module gradient
