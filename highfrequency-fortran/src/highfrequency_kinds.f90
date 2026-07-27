! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_kinds
  implicit none
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module highfrequency_kinds
