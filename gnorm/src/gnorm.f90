! SPDX-License-Identifier: GPL-2.0-or-later
module gnorm
   use gnorm_kinds, only : dp, i8
   use gnorm_status, only : gnorm_success, gnorm_invalid_argument, &
      gnorm_numerical_failure
   use gnorm_rng, only : gnorm_rng_state
   use gnorm_special, only : regularized_gamma_p, regularized_gamma_q, &
      inverse_regularized_gamma_p
   use gnorm_distribution, only : dgnorm, pgnorm, qgnorm, rgnorm, &
      rgnorm_fill, gnorm_mean, gnorm_variance, gnorm_parameters_valid
   implicit none
   public
end module gnorm
