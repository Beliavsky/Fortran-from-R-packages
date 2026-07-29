! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Dmitriy Mayorov
module vasicekfit_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter, public :: log_two_pi = log(2.0_dp * pi)
end module vasicekfit_kinds
