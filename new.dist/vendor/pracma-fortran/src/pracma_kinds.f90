! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: i8 = selected_int_kind(18)
   real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
   real(dp), parameter, public :: two_pi_dp = 2.0_dp*pi_dp
   real(dp), parameter, public :: half_pi_dp = 0.5_dp*pi_dp
   real(dp), parameter, public :: sqrt_two_pi_dp = sqrt(2.0_dp*pi_dp)
   real(dp), parameter, public :: eps_dp = epsilon(1.0_dp)
   real(dp), parameter, public :: tiny_dp = tiny(1.0_dp)
   real(dp), parameter, public :: huge_dp = huge(1.0_dp)
end module pracma_kinds
