! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_status
   implicit none
   private
   integer, parameter, public :: ht_success = 0
   integer, parameter, public :: ht_invalid_argument = 1
   integer, parameter, public :: ht_too_few_observations = 2
   integer, parameter, public :: ht_singular_matrix = 3
   integer, parameter, public :: ht_no_convergence = 4
   integer, parameter, public :: ht_numerical_error = 5
end module fitheavytail_status
