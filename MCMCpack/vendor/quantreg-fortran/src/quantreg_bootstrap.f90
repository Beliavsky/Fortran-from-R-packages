! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_bootstrap
  use quantreg_kinds, only : dp
  use quantreg_types, only : rq_result
  use quantreg_dense, only : rq_fit_fnb
  implicit none
  private
  public :: rq_bootstrap_xy, seed_rng
contains

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i, huge(1)-1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  subroutine rq_bootstrap_xy(x, y, tau, nrep, coefficients, info, seed)
    real(dp), intent(in) :: x(:,:), y(:), tau
    integer, intent(in) :: nrep
    real(dp), intent(out) :: coefficients(:,:)
    integer, intent(out) :: info
    integer, intent(in), optional :: seed
    real(dp), allocatable :: xb(:,:), yb(:)
    type(rq_result) :: fit
    real(dp) :: u
    integer :: i, j, k, n, p

    n = size(x,1)
    p = size(x,2)
    if (size(coefficients,1) /= p .or. size(coefficients,2) /= nrep) then
      info = -1
      return
    end if
    if (present(seed)) call seed_rng(seed)
    allocate(xb(n,p), yb(n))
    do k = 1, nrep
      do i = 1, n
        call random_number(u)
        j = min(n, int(u*real(n,dp))+1)
        xb(i,:) = x(j,:)
        yb(i) = y(j)
      end do
      call rq_fit_fnb(xb,yb,tau,fit)
      if (fit%info /= 0) then
        info = fit%info
        return
      end if
      coefficients(:,k) = fit%coefficients
    end do
    info = 0
  end subroutine rq_bootstrap_xy

end module quantreg_bootstrap
