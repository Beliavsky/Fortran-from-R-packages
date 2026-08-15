! Public umbrella module for distr-fortran.
! SPDX-License-Identifier: LGPL-3.0-only
module distr
   use distr_kinds, only : dp, pi
   use distr_special, only : digamma_value, inverse_digamma
   use distr_rng, only : seed_rng
   use distr_core
   use distr_matrix
   implicit none
   public
end module distr
