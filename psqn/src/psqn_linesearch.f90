! SPDX-License-Identifier: Apache-2.0
module psqn_linesearch
  use psqn_types, only : dp, psqn_objective_eval
  use psqn_interpolation, only : intrapolate_type
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: wolfe_line_search

contains

  logical function wolfe_line_search(eval_obj, f0, x0, gr0, direction, fnew, &
                                      c1, c2, strong_wolfe, n_eval, n_grad, trace) result(success)
    procedure(psqn_objective_eval) :: eval_obj
    real(dp), intent(in) :: f0
    real(dp), intent(inout) :: x0(:), gr0(:)
    real(dp), intent(in) :: direction(:)
    real(dp), intent(inout) :: fnew
    real(dp), intent(in) :: c1, c2
    logical, intent(in) :: strong_wolfe
    integer, intent(inout) :: n_eval, n_grad
    integer, intent(in), optional :: trace

    integer, parameter :: max_ls = 12
    real(dp), allocatable :: x_mem(:), g_dummy(:)
    real(dp) :: forg, dpsi_zero, fold, a_prev, ai, fi, dpsi_i, test_val, mult
    integer :: i, trace_use
    logical :: found_ok_prev, failed_once
    type(intrapolate_type) :: inter

    allocate(x_mem(size(x0)), g_dummy(size(x0)))
    trace_use = 0
    if (present(trace)) trace_use = trace
    forg = fnew
    dpsi_zero = dot_product(gr0, direction)
    if (dpsi_zero > 0.0_dp) then
      success = .false.
      return
    end if
    ! With a numerically zero directional derivative, the current point is
    ! already stationary to working precision. Treat the zero step as a
    ! successful line search; this also avoids roundoff-only Armijo failures.
    if (abs(dpsi_zero) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(f0))) then
      success = .true.
      return
    end if

    fold = f0
    a_prev = 0.0_dp
    ai = 0.25_dp
    found_ok_prev = .false.
    failed_once = .false.
    mult = 4.0_dp

    do i = 1, max_ls
      ai = ai * mult
      fi = psi(ai)
      call report_inner(a_prev, ai, fi, .false.)

      if (.not. ieee_is_finite(fi)) then
        failed_once = .true.
        mult = 0.5_dp
        if (.not. found_ok_prev) cycle
        fi = fold
        ai = a_prev
      end if

      if (fi > f0 + c1 * ai * dpsi_zero .or. (found_ok_prev .and. fi > fold)) then
        call inter%init(f0, dpsi_zero, ai, fi)
        success = zoom(a_prev, ai, inter)
        if (success .or. (ieee_is_finite(fnew) .and. fnew < forg)) then
          x0 = x_mem
        else
          fnew = forg
        end if
        return
      end if

      dpsi_i = dpsi(ai)
      call report_inner_deriv(a_prev, ai, fi, dpsi_i, .false.)
      if (strong_wolfe) then
        test_val = abs(dpsi_i)
      else
        test_val = -dpsi_i
      end if
      if (test_val <= -c2 * dpsi_zero) then
        x0 = x_mem
        success = .true.
        return
      end if

      if (failed_once .and. fi < f0) then
        x0 = x_mem
        success = .false.
        return
      end if

      if (dpsi_i >= 0.0_dp) then
        if (found_ok_prev) then
          call inter%init(f0, dpsi_zero, a_prev, fold)
          call inter%update(ai, fi)
        else
          call inter%init(f0, dpsi_zero, ai, fi)
        end if
        success = zoom(ai, a_prev, inter)
        if (success .or. (ieee_is_finite(fnew) .and. fnew < forg)) then
          x0 = x_mem
        else
          fnew = forg
        end if
        return
      end if

      found_ok_prev = .true.
      a_prev = ai
      fold = fi
    end do

    success = .false.

  contains

    real(dp) function psi(alpha) result(f)
      real(dp), intent(in) :: alpha
      x_mem = x0 + alpha * direction
      call eval_obj(x_mem, f, g_dummy, .false.)
      n_eval = n_eval + 1
    end function psi

    real(dp) function dpsi(alpha) result(df)
      real(dp), intent(in) :: alpha
      x_mem = x0 + alpha * direction
      call eval_obj(x_mem, fnew, gr0, .true.)
      n_grad = n_grad + 1
      df = dot_product(gr0, direction)
    end function dpsi

    logical function zoom(a_low_in, a_high_in, interp) result(ok)
      real(dp), intent(in) :: a_low_in, a_high_in
      type(intrapolate_type), intent(inout) :: interp
      real(dp) :: a_low, a_high, f_low, a_try, f_try, d_try, tval
      integer :: iz

      a_low = a_low_in
      a_high = a_high_in
      f_low = psi(a_low)
      do iz = 1, max_ls
        a_try = interp%get_value(a_low, a_high)
        f_try = psi(a_try)
        call report_inner(a_low, a_try, f_try, .true.)
        if (.not. ieee_is_finite(f_try)) then
          if (a_low < a_high) then
            a_high = a_try
          else
            a_low = a_try
          end if
          cycle
        end if

        call interp%update(a_try, f_try)
        if (f_try > f0 + c1 * a_try * dpsi_zero .or. f_try >= f_low) then
          a_high = a_try
          cycle
        end if

        d_try = dpsi(a_try)
        call report_inner_deriv(a_low, a_try, f_try, d_try, .true.)
        if (strong_wolfe) then
          tval = abs(d_try)
        else
          tval = -d_try
        end if
        if (tval <= -c2 * dpsi_zero) then
          ok = .true.
          return
        end if

        if (d_try * (a_high - a_low) >= 0.0_dp) a_high = a_low
        a_low = a_try
        f_low = f_try
      end do
      ok = .false.
    end function zoom

    subroutine report_inner(a_old, a_new, f_new, is_zoom)
      real(dp), intent(in) :: a_old, a_new, f_new
      logical, intent(in) :: is_zoom
      if (trace_use >= 4) then
        write(*,'(a,l1,3(a,es14.6))') 'line-search zoom=', is_zoom, &
          ' a_old=', a_old, ' a_new=', a_new, ' f=', f_new
      end if
    end subroutine report_inner

    subroutine report_inner_deriv(a_old, a_new, f_new, deriv, is_zoom)
      real(dp), intent(in) :: a_old, a_new, f_new, deriv
      logical, intent(in) :: is_zoom
      if (trace_use >= 5) then
        write(*,'(a,l1,4(a,es14.6))') 'line-search zoom=', is_zoom, &
          ' a_old=', a_old, ' a_new=', a_new, ' f=', f_new, ' d=', deriv
      end if
    end subroutine report_inner_deriv

  end function wolfe_line_search

end module psqn_linesearch
