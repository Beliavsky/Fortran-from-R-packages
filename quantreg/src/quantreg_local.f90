! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_local
  use quantreg_kinds, only : dp
  use quantreg_types, only : rq_result, lprq_result
  use quantreg_dense, only : rq_wfit_fnb
  implicit none
  private
  public :: lprq
contains

  subroutine lprq(x, y, bandwidth, tau, ngrid, result)
    real(dp), intent(in) :: x(:), y(:), bandwidth, tau
    integer, intent(in) :: ngrid
    type(lprq_result), intent(out) :: result
    real(dp), allocatable :: design(:,:), weights(:), z(:)
    type(rq_result) :: fit
    real(dp) :: xmin, xmax
    integer :: i, n

    n = size(x)
    if (size(y) /= n .or. bandwidth <= 0.0_dp .or. ngrid < 2) then
      result%info = -1
      return
    end if
    allocate(result%x(ngrid), result%fitted(ngrid), result%derivative(ngrid))
    allocate(design(n,2), weights(n), z(n))
    xmin = minval(x)
    xmax = maxval(x)
    do i = 1, ngrid
      result%x(i) = xmin + real(i-1,dp)*(xmax-xmin)/real(ngrid-1,dp)
      z = x - result%x(i)
      design(:,1) = 1.0_dp
      design(:,2) = z
      weights = exp(-0.5_dp*(z/bandwidth)**2)
      call rq_wfit_fnb(design, y, weights, tau, fit)
      if (fit%info /= 0) then
        result%info = fit%info
        return
      end if
      result%fitted(i) = fit%coefficients(1)
      result%derivative(i) = fit%coefficients(2)
    end do
    result%info = 0
  end subroutine lprq

end module quantreg_local
