program test_minqa
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use minqa_module, only : dp, minqa_control_t, minqa_result_t, bobyqa, newuoa, uobyqa
   implicit none

   integer :: failures

   failures = 0
   call test_quadratic()
   call test_bounded_rosenbrock()
   call test_maxfun()
   call test_invalid_bounds()
   call test_nonfinite_objective()
   call test_chebyquad()

   if (failures /= 0) then
      write(*, '("minqa tests failed: ",i0)') failures
      error stop 1
   end if
   write(*, '("All minqa tests passed.")')

contains

   subroutine test_quadratic()
      real(dp) :: x(4)
      type(minqa_control_t) :: control
      type(minqa_result_t) :: result

      control%rhobeg = 0.5_dp
      control%rhoend = 1.0e-8_dp
      control%maxfun = 10000

      x = acos(-1.0_dp)
      call newuoa(shifted_quadratic, x, result, control)
      call check(result%status == 0, 'NEWUOA normal status')
      call check(maxval(abs(x - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])) < 5.0e-6_dp, &
         'NEWUOA quadratic minimizer')
      call check(abs(result%fval + 10.0_dp) < 1.0e-10_dp, 'NEWUOA objective')

      x = acos(-1.0_dp)
      call uobyqa(shifted_quadratic, x, result, control)
      call check(result%status == 0, 'UOBYQA normal status')
      call check(maxval(abs(x - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])) < 5.0e-7_dp, &
         'UOBYQA quadratic minimizer')
      call check(abs(result%fval + 10.0_dp) < 1.0e-10_dp, 'UOBYQA objective')

      x = acos(-1.0_dp)
      call bobyqa(shifted_quadratic, x, result, -10.0_dp, 10.0_dp, control)
      call check(result%status == 0, 'BOBYQA normal status')
      call check(maxval(abs(x - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])) < 5.0e-6_dp, &
         'BOBYQA quadratic minimizer')
      call check(abs(result%fval + 10.0_dp) < 1.0e-10_dp, 'BOBYQA objective')
   end subroutine test_quadratic

   subroutine test_bounded_rosenbrock()
      real(dp) :: x(2), lower(2), upper(2)
      type(minqa_control_t) :: control
      type(minqa_result_t) :: result

      x = [-1.2_dp, 1.0_dp]
      lower = [-2.0_dp, -1.0_dp]
      upper = [ 2.0_dp,  3.0_dp]
      control%rhobeg = 0.25_dp
      control%rhoend = 1.0e-7_dp
      control%maxfun = 5000
      control%npt = 5

      call bobyqa(rosenbrock, x, result, lower, upper, control)
      call check(result%status == 0, 'bounded Rosenbrock status')
      call check(maxval(abs(x - 1.0_dp)) < 2.0e-4_dp, 'bounded Rosenbrock minimizer')
      call check(result%fval < 1.0e-10_dp, 'bounded Rosenbrock objective')
      call check(all(x >= lower .and. x <= upper), 'BOBYQA respects bounds')
   end subroutine test_bounded_rosenbrock

   subroutine test_maxfun()
      real(dp) :: x(4)
      type(minqa_control_t) :: control
      type(minqa_result_t) :: result

      x = acos(-1.0_dp)
      control%rhobeg = 0.5_dp
      control%rhoend = 1.0e-12_dp
      control%maxfun = 10
      call newuoa(shifted_quadratic, x, result, control)
      call check(result%status == 1, 'maximum-evaluation status mapping')
      call check(result%raw_status == 390, 'maximum-evaluation raw status')
      call check(result%evaluations <= control%maxfun, 'maximum-evaluation count')
   end subroutine test_maxfun

   subroutine test_invalid_bounds()
      real(dp) :: x(2), lower(2), upper(2)
      type(minqa_result_t) :: result

      x = [0.0_dp, 0.0_dp]
      lower = [1.0_dp, -1.0_dp]
      upper = [0.0_dp,  1.0_dp]
      call bobyqa(rosenbrock, x, result, lower, upper)
      call check(result%status == 6, 'invalid-bound status')
      call check(result%evaluations == 0, 'invalid bounds avoid evaluations')
   end subroutine test_invalid_bounds

   subroutine test_nonfinite_objective()
      real(dp) :: x(2)
      type(minqa_control_t) :: control
      type(minqa_result_t) :: result

      x = [0.5_dp, 0.5_dp]
      control%rhobeg = 0.1_dp
      control%rhoend = 1.0e-5_dp
      control%maxfun = 200
      call bobyqa(partly_nan, x, result, 0.0_dp, 1.0_dp, control)
      call check(result%evaluations > 0, 'nonfinite objective handled')
      call check(all(x >= 0.0_dp .and. x <= 1.0_dp), 'nonfinite objective preserves bounds')
   end subroutine test_nonfinite_objective


   subroutine test_chebyquad()
      real(dp) :: x(6)
      type(minqa_control_t) :: control
      type(minqa_result_t) :: result

      control%rhobeg = 0.2_dp
      control%rhoend = 1.0e-7_dp
      control%maxfun = 25000

      x = 1.0_dp
      call newuoa(chebyquad, x, result, control)
      call check(result%status == 0, 'NEWUOA Chebyquad status')
      call check(result%fval < 1.0e-10_dp, 'NEWUOA Chebyquad objective')

      x = 1.0_dp
      call uobyqa(chebyquad, x, result, control)
      call check(result%status == 0, 'UOBYQA Chebyquad status')
      call check(result%fval < 1.0e-12_dp, 'UOBYQA Chebyquad objective')

      x = 1.0_dp
      call bobyqa(chebyquad, x, result, -10.0_dp, 10.0_dp, control)
      call check(result%status == 0, 'BOBYQA Chebyquad status')
      call check(result%fval < 1.0e-8_dp, 'BOBYQA Chebyquad objective')
   end subroutine test_chebyquad

   function shifted_quadratic(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      integer :: i

      f = sum((x - [(real(i, dp), i = 1, size(x))])**2) - 10.0_dp
   end function shifted_quadratic

   function rosenbrock(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f

      f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
   end function rosenbrock


   function chebyquad(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      real(dp) :: residual(size(x))
      real(dp) :: rr, z2, z6, z7, z8
      integer :: i, j, k, n

      n = size(x)
      residual = 0.0_dp
      do i = 1, n
         rr = 0.0_dp
         do k = 1, n
            z7 = 1.0_dp
            z2 = 2.0_dp * x(k) - 1.0_dp
            z8 = z2
            j = 1
            do while (j < i)
               z6 = z7
               z7 = z8
               z8 = 2.0_dp * z2 * z7 - z6
               j = j + 1
            end do
            rr = rr + z8
         end do
         rr = rr / real(n, dp)
         if (mod(i, 2) == 0) rr = rr + 1.0_dp / real(i * i - 1, dp)
         residual(i) = rr
      end do
      f = sum(residual**2)
   end function chebyquad

   function partly_nan(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f

      if (x(1) > 0.9_dp) then
         f = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         f = sum((x - 0.25_dp)**2)
      end if
   end function partly_nan

   subroutine check(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name

      if (.not. condition) then
         failures = failures + 1
         write(*, '("FAIL: ",a)') name
      end if
   end subroutine check

end program test_minqa
