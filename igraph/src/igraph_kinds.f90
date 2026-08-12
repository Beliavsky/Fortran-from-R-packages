module igraph_kinds
  use iso_fortran_env, only : real64, int64
  implicit none
  private
  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: igraph_inf = huge(1.0_dp) / 16.0_dp
end module igraph_kinds
