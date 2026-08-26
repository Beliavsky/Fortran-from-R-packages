! cmaes-fortran - GPL-2.0-only
module cmaes_linalg
  use cmaes_kinds, only : dp
  use r_linalg, only : symmetric_eigen
  implicit none
  private
  public :: symmetric_eigen_descending

contains
  subroutine symmetric_eigen_descending(a, values, vectors, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out) :: info
    call symmetric_eigen(a, values, vectors, info, descending=.true.)
  end subroutine symmetric_eigen_descending
end module cmaes_linalg
