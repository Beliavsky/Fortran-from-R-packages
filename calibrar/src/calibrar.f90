! SPDX-License-Identifier: GPL-2.0-only
module calibrar
  use calibrar_kinds, only : dp, pi_dp
  use calibrar_interfaces, only : scalar_objective, gradient_callback, vector_objective
  use calibrar_gradient
  use calibrar_fitness
  use calibrar_random
  use calibrar_splines
  use calibrar_objective
  use calibrar_stopping
  use calibrar_optimization
  use calibrar_test_functions
  implicit none
  public
end module calibrar
