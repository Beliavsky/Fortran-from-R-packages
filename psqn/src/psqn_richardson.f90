! SPDX-License-Identifier: Apache-2.0
module psqn_richardson
  use psqn_types, only : dp
  implicit none
  private
  public :: richardson_vector_derivative

  abstract interface
    subroutine vector_fun(x, out)
      import dp
      real(dp), intent(in) :: x
      real(dp), intent(out) :: out(:)
    end subroutine vector_fun
  end interface

contains

  subroutine richardson_vector_derivative(fun, x, out, eps, scale, tol, order)
    procedure(vector_fun) :: fun
    real(dp), intent(in) :: x
    real(dp), intent(out) :: out(:)
    real(dp), intent(in), optional :: eps, scale, tol
    integer, intent(in), optional :: order
    real(dp) :: eps_use, scale_use, tol_use, step, delta, scale_sq, mult, err_est
    integer :: order_use, n_vars, i, j, k
    real(dp), allocatable :: fplus(:), fminus(:), threshold(:), apprx(:,:)
    logical, allocatable :: converged(:)
    logical :: passed

    eps_use = 1.0e-3_dp
    if (present(eps)) eps_use = eps
    scale_use = 2.0_dp
    if (present(scale)) scale_use = scale
    tol_use = 1.0e-9_dp
    if (present(tol)) tol_use = tol
    order_use = 6
    if (present(order)) order_use = order

    if (eps_use <= 0.0_dp) error stop "richardson: eps <= 0"
    if (scale_use <= 1.0_dp) error stop "richardson: scale <= 1"
    if (tol_use <= 0.0_dp) error stop "richardson: tol <= 0"
    if (order_use < 0) error stop "richardson: order < 0"

    n_vars = size(out)
    allocate(fplus(n_vars), fminus(n_vars), threshold(n_vars))
    allocate(apprx(n_vars, 0:order_use))
    allocate(converged(n_vars))
    converged = .false.

    step = max(eps_use, abs(x) * eps_use)
    delta = step
    call central(delta, apprx(:,0))

    if (order_use > 0) then
      call fun(x, threshold)
      threshold = max(tol_use, abs(threshold) * tol_use)
    end if

    scale_sq = scale_use * scale_use
    do i = 0, order_use - 1
      delta = delta / scale_use
      call central(delta, apprx(:,i+1))

      mult = 1.0_dp
      do j = i, 0, -1
        mult = mult * scale_sq
        do k = 1, n_vars
          if (.not. converged(k)) then
            apprx(k,j) = apprx(k,j+1) + (apprx(k,j+1) - apprx(k,j)) / (mult - 1.0_dp)
          end if
        end do
      end do

      passed = i > 0
      do k = 1, n_vars
        if (.not. converged(k)) then
          err_est = (apprx(k,0) - apprx(k,1)) * mult / (mult - 1.0_dp)
          converged(k) = abs(err_est) < threshold(k)
          passed = passed .and. converged(k)
        end if
      end do
      if (passed) exit

      mult = mult * scale_sq
      do k = 1, n_vars
        if (.not. converged(k)) then
          apprx(k,0) = apprx(k,1) + (apprx(k,1) - apprx(k,0)) / (mult - 1.0_dp)
        end if
      end do
    end do

    out = apprx(:,0)

  contains
    subroutine central(h, deriv)
      real(dp), intent(in) :: h
      real(dp), intent(out) :: deriv(:)
      call fun(x + h, fplus)
      call fun(x - h, fminus)
      deriv = (fplus - fminus) / (2.0_dp * h)
    end subroutine central
  end subroutine richardson_vector_derivative

end module psqn_richardson
