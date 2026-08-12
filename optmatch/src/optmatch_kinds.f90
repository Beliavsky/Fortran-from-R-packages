! SPDX-License-Identifier: MIT
module optmatch_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: optmatch_inf = huge(1.0_dp) / 16.0_dp
end module optmatch_kinds
