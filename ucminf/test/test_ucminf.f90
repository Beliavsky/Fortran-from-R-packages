module test_ucminf_context_support
   use ucminf, only : dp
   implicit none

   type :: quadratic_context
      real(dp) :: center(2) = [3.0_dp, -2.0_dp]
   end type quadratic_context

contains

   real(dp) function contextual_quadratic(x, raw_context) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout) :: raw_context

      select type (context => raw_context)
      type is (quadratic_context)
         f = (x(1) - context%center(1))**2 + 2.0_dp * (x(2) - context%center(2))**2
      class default
         f = huge(1.0_dp)
      end select
   end function contextual_quadratic

end module test_ucminf_context_support

program test_ucminf
   use ucminf, only : dp, UCMINF_GRAD_CENTRAL, ucminf_options, ucminf_result, &
      ucminf_gradient_check, ucminf_minimize, ucminf_minimize_context, ucminf_check_gradient
   use test_ucminf_context_support, only : quadratic_context, contextual_quadratic
   implicit none

   call test_rosenbrock_analytic()
   call test_rosenbrock_forward()
   call test_rosenbrock_central()
   call test_quadratic_initial_hessian()
   call test_context_callback()
   call test_gradient_checker()
   call test_invalid_options()
   write(*,'(a)') "All ucminf tests passed."

contains

   subroutine test_rosenbrock_analytic()
      type(ucminf_result) :: r
      real(dp) :: x0(2)
      x0 = [2.0_dp, 0.5_dp]
      call ucminf_minimize(x0, rosenbrock, r, rosenbrock_gradient)
      call assert_true(r%convergence == 1 .or. r%convergence == 2, "analytic Rosenbrock convergence")
      call assert_close(r%par(1), 1.0_dp, 5.0e-6_dp, "analytic Rosenbrock x1")
      call assert_close(r%par(2), 1.0_dp, 5.0e-6_dp, "analytic Rosenbrock x2")
      call assert_true(r%value < 1.0e-10_dp, "analytic Rosenbrock objective")
      call assert_true(r%neval == 21, "analytic Rosenbrock upstream evaluation count")
   end subroutine test_rosenbrock_analytic


   subroutine test_rosenbrock_forward()
      type(ucminf_result) :: r
      real(dp) :: x0(2)
      x0 = [2.0_dp, 0.5_dp]
      call ucminf_minimize(x0, rosenbrock, r)
      call assert_true(r%convergence == 4, "forward Rosenbrock upstream convergence code")
      call assert_close(r%par(1), 9.9969996210797296e-1_dp, 1.0e-9_dp, "forward Rosenbrock x1")
      call assert_close(r%par(2), 9.9939932383899688e-1_dp, 1.0e-9_dp, "forward Rosenbrock x2")
      call assert_true(r%neval == 19, "forward Rosenbrock upstream evaluation count")
   end subroutine test_rosenbrock_forward

   subroutine test_rosenbrock_central()
      type(ucminf_options) :: opt
      type(ucminf_result) :: r
      real(dp) :: x0(2)
      x0 = [2.0_dp, 0.5_dp]
      opt%grad_method = UCMINF_GRAD_CENTRAL
      call ucminf_minimize(x0, rosenbrock, r, options=opt)
      call assert_true(r%convergence == 1 .or. r%convergence == 2, "central Rosenbrock convergence")
      call assert_close(r%par(1), 1.0_dp, 2.0e-5_dp, "central Rosenbrock x1")
      call assert_close(r%par(2), 1.0_dp, 2.0e-5_dp, "central Rosenbrock x2")
      call assert_true(r%neval == 21, "central Rosenbrock upstream evaluation count")
   end subroutine test_rosenbrock_central

   subroutine test_quadratic_initial_hessian()
      type(ucminf_options) :: opt
      type(ucminf_result) :: r
      real(dp) :: x0(2)
      x0 = [8.0_dp, -4.0_dp]
      allocate(opt%invhessian_lt(3))
      opt%invhessian_lt = [0.5_dp, 0.0_dp, 0.25_dp]
      call ucminf_minimize(x0, quadratic, r, quadratic_gradient, opt)
      call assert_true(r%convergence == 1 .or. r%convergence == 2, "quadratic convergence")
      call assert_close(r%par(1), 3.0_dp, 1.0e-9_dp, "quadratic x1")
      call assert_close(r%par(2), -2.0_dp, 1.0e-9_dp, "quadratic x2")
      call assert_true(size(r%invhessian_lt) == 3, "inverse Hessian packed size")
   end subroutine test_quadratic_initial_hessian

   subroutine test_context_callback()
      type(quadratic_context) :: context
      type(ucminf_result) :: r
      real(dp) :: x0(2)

      x0 = [8.0_dp, -4.0_dp]
      call ucminf_minimize_context(x0, contextual_quadratic, context, r)
      call assert_true(r%convergence > 0, "context callback convergence")
      call assert_close(r%par(1), context%center(1), 1.0e-4_dp, "context callback x1")
      call assert_close(r%par(2), context%center(2), 1.0e-4_dp, "context callback x2")
   end subroutine test_context_callback

   subroutine test_gradient_checker()
      type(ucminf_gradient_check) :: report
      real(dp) :: x(2)
      x = [0.7_dp, 0.8_dp]
      call ucminf_check_gradient(x, rosenbrock, rosenbrock_gradient, 1.0e-5_dp, report)
      call assert_true(report%success, "gradient checker success")
      call assert_true(abs(report%max_extrapolated_error) < 1.0e-6_dp, "gradient checker extrapolated error")
   end subroutine test_gradient_checker

   subroutine test_invalid_options()
      type(ucminf_options) :: opt
      type(ucminf_result) :: r
      real(dp) :: x0(1)
      x0 = [0.0_dp]
      opt%stepmax = 0.0_dp
      call ucminf_minimize(x0, scalar_quadratic, r, options=opt)
      call assert_true(r%convergence == -4, "invalid stepmax code")
   end subroutine test_invalid_options

   function rosenbrock(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2
   end function rosenbrock

   subroutine rosenbrock_gradient(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - 2.0_dp * (1.0_dp - x(1))
      g(2) = 200.0_dp * (x(2) - x(1)**2)
   end subroutine rosenbrock_gradient

   function quadratic(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = (x(1) - 3.0_dp)**2 + 2.0_dp * (x(2) + 2.0_dp)**2
   end function quadratic

   subroutine quadratic_gradient(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g(1) = 2.0_dp * (x(1) - 3.0_dp)
      g(2) = 4.0_dp * (x(2) + 2.0_dp)
   end subroutine quadratic_gradient

   function scalar_quadratic(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = x(1)**2
   end function scalar_quadratic

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') "FAILED: " // trim(label)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual, expected, tol
      character(len=*), intent(in) :: label
      call assert_true(abs(actual - expected) <= tol, label)
   end subroutine assert_close

end program test_ucminf
