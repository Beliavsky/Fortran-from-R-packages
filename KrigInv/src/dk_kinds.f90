! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! Vendored in KrigInv-fortran under the GPL-3 option; see
! licenses/DiceKriging/LICENSE-GPL-3 and UPSTREAM.md.
module dk_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
end module dk_kinds
