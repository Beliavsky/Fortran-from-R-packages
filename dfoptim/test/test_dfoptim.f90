! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

program test_dfoptim
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
   use dfoptim, only : dp, dfoptim_result_t, hj_control_t, nmk_control_t, mads_control_t, &
      hjk, hjkb, nmk, nmkb, mads, mads_poll_full, dfoptim_success, &
      dfoptim_max_evaluations, dfoptim_invalid_input
   implicit none

   type :: shift_data_t
      real(dp) :: target(2) = [2.0_dp, -3.0_dp]
   end type shift_data_t

   integer :: failures
   failures = 0

   call test_hjk(failures)
   call test_hjkb(failures)
   call test_nmk(failures)
   call test_nmkb(failures)
   call test_mads(failures)
   call test_maximize(failures)
   call test_statuses(failures)
   call test_one_sided_bounds(failures)
   call test_user_data_and_cancel(failures)

   if (failures > 0) then
      write(*, '(a,i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   write(*, '(a)') 'All dfoptim tests passed.'

contains

   subroutine test_hjk(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(hj_control_t) :: c
      real(dp) :: x0(2)
      x0 = [-1.2_dp, 1.0_dp]
      c%tol = 1.0e-7_dp
      c%maxfeval = 20000
      c%seed = 123
      r = hjk(x0, rosenbrock, c)
      call check(r%convergence == dfoptim_success, 'hjk convergence', failures)
      call check(r%value < 1.0e-8_dp, 'hjk Rosenbrock value', failures)
      call check(maxval(abs(r%x - 1.0_dp)) < 5.0e-4_dp, 'hjk Rosenbrock parameters', failures)
   end subroutine test_hjk

   subroutine test_hjkb(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(hj_control_t) :: c
      real(dp) :: x0(3), lo(3), up(3)
      x0 = 0.0_dp
      lo = -2.0_dp
      up = 0.5_dp
      c%tol = 1.0e-7_dp
      c%maxfeval = 20000
      c%seed = 321
      r = hjkb(x0, extended_rosenbrock, lo, up, c)
      call check(r%convergence == dfoptim_success, 'hjkb convergence', failures)
      call check(all(r%x >= lo .and. r%x <= up), 'hjkb bounds', failures)
      call check(abs(r%x(1) - 0.5_dp) < 1.0e-8_dp, 'hjkb active upper bound', failures)
      call check(r%value < 0.85_dp, 'hjkb bounded objective', failures)
   end subroutine test_hjkb

   subroutine test_nmk(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(nmk_control_t) :: c
      real(dp) :: x0(2)
      x0 = [-1.2_dp, 1.0_dp]
      c%tol = 1.0e-10_dp
      c%maxfeval = 10000
      c%max_restarts = 8
      r = nmk(x0, rosenbrock, c)
      call check(r%convergence == dfoptim_success, 'nmk convergence', failures)
      call check(r%value < 1.0e-10_dp, 'nmk Rosenbrock value', failures)
      call check(maxval(abs(r%x - 1.0_dp)) < 1.0e-4_dp, 'nmk Rosenbrock parameters', failures)
   end subroutine test_nmk

   subroutine test_nmkb(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(nmk_control_t) :: c
      real(dp) :: x0(2), lo(2), up(2)
      x0 = [0.0_dp, 0.0_dp]
      lo = [-0.5_dp, -2.0_dp]
      up = [0.5_dp, 2.0_dp]
      c%tol = 1.0e-10_dp
      c%maxfeval = 10000
      c%max_restarts = 8
      r = nmkb(x0, shifted_sphere, lo, up, c)
      call check(r%convergence == dfoptim_success, 'nmkb convergence', failures)
      call check(all(r%x >= lo .and. r%x <= up), 'nmkb bounds', failures)
      call check(abs(r%x(1) - 0.5_dp) < 1.0e-5_dp, 'nmkb active bound', failures)
      call check(abs(r%x(2) + 1.0_dp) < 1.0e-5_dp, 'nmkb free parameter', failures)
   end subroutine test_nmkb

   subroutine test_mads(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(mads_control_t) :: c
      real(dp) :: x0(2), lo(2), up(2), sc(2)
      x0 = [0.0_dp, 0.0_dp]
      lo = [-2.0_dp, -2.0_dp]
      up = [2.0_dp, 2.0_dp]
      sc = [0.5_dp, 0.5_dp]
      c%trace = .false.
      c%tol = 1.0e-5_dp
      c%maxfeval = 20000
      c%poll_style = mads_poll_full
      c%delta_init = 0.05_dp
      c%seed = 1138
      r = mads(x0, nonsmooth, lo, up, sc, c)
      call check(r%convergence == dfoptim_success, 'mads convergence', failures)
      call check(r%value < 1.0e-4_dp, 'mads nonsmooth value', failures)
      call check(maxval(abs(r%x - [0.75_dp, -0.25_dp])) < 2.0e-2_dp, 'mads parameters', failures)
      call check(size(r%log%value) == r%niter + 1, 'mads iteration log', failures)
   end subroutine test_mads

   subroutine test_maximize(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(hj_control_t) :: c
      real(dp) :: x0(2)
      x0 = [0.0_dp, 0.0_dp]
      c%maximize = .true.
      c%tol = 1.0e-7_dp
      c%maxfeval = 10000
      r = hjk(x0, concave, c)
      call check(r%value > 2.999999_dp, 'maximize value', failures)
      call check(maxval(abs(r%x - [1.0_dp, -2.0_dp])) < 1.0e-4_dp, 'maximize parameters', failures)
   end subroutine test_maximize

   subroutine test_statuses(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(hj_control_t) :: hc
      type(nmk_control_t) :: nc
      real(dp) :: x0(2), lo(2), up(2)
      x0 = [-1.2_dp, 1.0_dp]
      hc%maxfeval = 2
      r = hjk(x0, rosenbrock, hc)
      call check(r%convergence == dfoptim_max_evaluations, 'evaluation limit status', failures)

      lo = [0.0_dp, 0.0_dp]
      up = [1.0_dp, 1.0_dp]
      nc%maxfeval = 100
      r = nmkb([0.0_dp, 0.5_dp], sphere, lo, up, nc)
      call check(r%convergence == dfoptim_invalid_input, 'nmkb boundary validation', failures)
   end subroutine test_statuses


   subroutine test_one_sided_bounds(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(nmk_control_t) :: c
      real(dp) :: x0(2), lo(2), up(2)
      x0 = [0.5_dp, 0.0_dp]
      lo = [0.0_dp, ieee_value(1.0_dp, ieee_negative_inf)]
      up = [ieee_value(1.0_dp, ieee_positive_inf), 1.0_dp]
      c%tol = 1.0e-10_dp
      c%maxfeval = 10000
      c%max_restarts = 8
      r = nmkb(x0, one_sided_objective, lo, up, c)
      call check(r%convergence == dfoptim_success, 'nmkb one-sided convergence', failures)
      call check(maxval(abs(r%x - [1.5_dp, -0.5_dp])) < 1.0e-4_dp, &
         'nmkb one-sided parameters', failures)
   end subroutine test_one_sided_bounds

   subroutine test_user_data_and_cancel(failures)
      integer, intent(inout) :: failures
      type(dfoptim_result_t) :: r
      type(nmk_control_t) :: nc
      type(hj_control_t) :: hc
      type(shift_data_t) :: data
      real(dp) :: x0(2)
      x0 = 0.0_dp
      nc%tol = 1.0e-10_dp
      nc%maxfeval = 5000
      nc%max_restarts = 8
      r = nmk(x0, data_objective, nc, data)
      call check(maxval(abs(r%x - data%target)) < 1.0e-5_dp, 'user data callback', failures)

      hc%maxfeval = 10000
      r = hjk(x0, sphere, hc, monitor=cancel_monitor)
      call check(r%convergence /= dfoptim_success, 'monitor cancellation', failures)
      call check(r%niter == 1, 'monitor cancellation iteration', failures)
   end subroutine test_user_data_and_cancel

   subroutine check(condition, name, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // trim(name)
      end if
   end subroutine check

   function rosenbrock(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      f = 100.0_dp * (x(2) - x(1) ** 2) ** 2 + (1.0_dp - x(1)) ** 2
   end function rosenbrock

   function extended_rosenbrock(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      integer :: i
      f = 0.0_dp
      do i = 1, size(x) - 1
         f = f + 100.0_dp * (x(i) ** 2 - x(i + 1)) ** 2 + (x(i) - 1.0_dp) ** 2
      end do
   end function extended_rosenbrock

   function sphere(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      f = sum(x * x)
   end function sphere

   function shifted_sphere(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      f = (x(1) - 1.0_dp) ** 2 + (x(2) + 1.0_dp) ** 2
   end function shifted_sphere

   function nonsmooth(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      f = abs(x(1) - 0.75_dp) + 2.0_dp * abs(x(2) + 0.25_dp)
   end function nonsmooth

   function concave(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      f = 3.0_dp - (x(1) - 1.0_dp) ** 2 - (x(2) + 2.0_dp) ** 2
   end function concave

   function one_sided_objective(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      f = (x(1) - 1.5_dp) ** 2 + (x(2) + 0.5_dp) ** 2
   end function one_sided_objective

   function data_objective(x, user_data) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
      f = huge(1.0_dp)
      if (present(user_data)) then
         select type (user_data)
         type is (shift_data_t)
            f = sum((x - user_data%target) ** 2)
         end select
      end if
   end function data_objective

   subroutine cancel_monitor(x, value, iteration, evaluations, stop, user_data)
      real(dp), intent(in) :: x(:), value
      integer, intent(in) :: iteration, evaluations
      logical, intent(out) :: stop
      class(*), intent(inout), optional :: user_data
      stop = iteration >= 1
   end subroutine cancel_monitor

end program test_dfoptim
