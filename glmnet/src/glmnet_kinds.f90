! SPDX-License-Identifier: GPL-2.0-only
module glmnet_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: glmnet_eps = epsilon(1.0_dp)
   real(dp), parameter, public :: glmnet_tiny = tiny(1.0_dp)
   real(dp), parameter, public :: glmnet_huge = huge(1.0_dp)
end module glmnet_kinds
