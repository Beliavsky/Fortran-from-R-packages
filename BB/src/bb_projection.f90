! SPDX-License-Identifier: GPL-3.0-only
module bb_projection
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bb_kinds, only: dp
  use quadprog, only: qp_result, solve_qp
  implicit none
  private

  public :: project_box, project_linear

contains

  subroutine project_box(x, lower, upper, projected, ok)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: lower(:)
    real(dp), intent(in) :: upper(:)
    real(dp), intent(out) :: projected(:)
    logical, intent(out) :: ok

    if (size(lower) /= size(x) .or. size(upper) /= size(x) .or. &
        size(projected) /= size(x)) then
      projected = 0.0_dp
      ok = .false.
      return
    end if
    if (any(lower > upper)) then
      projected = 0.0_dp
      ok = .false.
      return
    end if

    projected = min(max(x, lower), upper)
    ok = all(ieee_is_finite(projected))
  end subroutine project_box

  subroutine project_linear(x, a, b, meq, projected, ok)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:)
    integer, intent(in) :: meq
    real(dp), intent(out) :: projected(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: dmat(:, :), dvec(:), amat(:, :), rhs(:)
    type(qp_result) :: fit
    integer :: n, m, i
    logical :: need_projection

    n = size(x)
    m = size(a, 1)
    projected = x
    ok = .false.

    if (size(a, 2) /= n .or. size(b) /= m .or. size(projected) /= n) return
    if (meq < 0 .or. meq > m) return
    if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(a)) .or. &
        .not. all(ieee_is_finite(b))) return

    rhs = b - matmul(a, x)
    need_projection = meq > 0
    if (m > meq) need_projection = need_projection .or. any(rhs(meq + 1:m) > 0.0_dp)
    if (.not. need_projection) then
      ok = .true.
      return
    end if

    allocate(dmat(n, n), dvec(n), amat(n, m))
    dmat = 0.0_dp
    do i = 1, n
      dmat(i, i) = 1.0_dp
    end do
    dvec = 0.0_dp
    amat = transpose(a)

    ! Identity is its own inverse Cholesky factor, matching BB's
    ! solve.QP(..., Dmat=diag(1,n), factorized=TRUE) call.
    fit = solve_qp(dmat, dvec, amat, rhs, meq=meq, factorized=.true.)
    if (.not. fit%succeeded()) return

    projected = x + fit%solution
    ok = all(ieee_is_finite(projected))
  end subroutine project_linear

end module bb_projection
