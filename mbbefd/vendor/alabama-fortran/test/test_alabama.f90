program test_alabama
   use alabama
   implicit none
   real(dp) :: x(2)
   type(alabama_result_t) :: res
   type(alabama_outer_control_t) :: oc
   type(alabama_inner_control_t) :: ic
   integer :: failures

   failures = 0
   ic%max_iterations = 300
   ic%reltol = 1.0e-10_dp

   ! New augmented Lagrangian: both constraint types, infeasible start.
   x = [-0.2_dp, 2.6_dp]
   oc = alabama_outer_control_t()
   oc%trace = .false.
   call auglag3(x, quad_obj, positive_xy, sum_two, res, quad_grad, &
                positive_xy_jac, sum_two_jac, oc, ic)
   call check_vec('auglag3 analytic', x, [0.5_dp, 1.5_dp], 2.0e-4_dp, failures)
   call check_true('auglag3 equality', abs(sum(x)-2.0_dp) < 2.0e-6_dp, failures)
   call check_true('auglag3 inequalities', minval(x) > -2.0e-6_dp, failures)

   ! Same problem with all derivatives from numDeriv.
   x = [-0.1_dp, 2.4_dp]
   call auglag3(x, quad_obj, positive_xy, sum_two, res, &
                control_outer=oc, control_inner=ic)
   call check_vec('auglag3 numerical derivatives', x, [0.5_dp, 1.5_dp], &
                  8.0e-4_dp, failures)

   ! Three-variable example from the R package documentation.
   block
      real(dp) :: y(3)
      y = [0.3_dp, 0.4_dp, 0.5_dp]
      call auglag3(y, documented_obj, documented_hin, documented_heq, res, &
                   documented_grad, documented_hin_jac, documented_heq_jac, oc, ic)
      call check_vec('documented auglag example', y, [0.0_dp, 0.0_dp, 1.0_dp], &
                     2.0e-3_dp, failures)
   end block

   ! Equality-only path.
   x = [0.0_dp, 0.0_dp]
   call auglag1(x, quad_obj, sum_two, res, quad_grad, sum_two_jac, oc, ic)
   call check_vec('auglag1', x, [0.5_dp, 1.5_dp], 2.0e-4_dp, failures)

   ! Inequality-only path with active boundary x(1) >= 0.
   x = [-2.0_dp, 3.0_dp]
   call auglag2(x, boundary_obj, x_nonnegative, res, boundary_grad, &
                x_nonnegative_jac, oc, ic)
   call check_vec('auglag2 boundary', x, [0.0_dp, 2.0_dp], 5.0e-4_dp, failures)

   ! Legacy constrOptim.nl path requires strict inequality feasibility.
   x = [0.8_dp, 1.2_dp]
   oc = alabama_outer_control_t(mu0=0.01_dp, sig0=10.0_dp, eps=1.0e-7_dp, &
                                itmax=50, trace=.false.)
   call constr_optim_nl(x, quad_obj, res, quad_grad, positive_xy, &
                        positive_xy_jac, sum_two, sum_two_jac, oc, ic)
   call check_vec('constr_optim_nl', x, [0.5_dp, 1.5_dp], 2.0e-3_dp, failures)

   ! Equality-only legacy augmented penalty.
   x = [0.2_dp, 0.2_dp]
   call augpen(x, quad_obj, sum_two, res, quad_grad, sum_two_jac, oc, ic)
   call check_vec('augpen', x, [0.5_dp, 1.5_dp], 2.0e-3_dp, failures)

   ! Inequality-only adaptive barrier with an interior unconstrained solution.
   x = [0.2_dp, 0.3_dp]
   call adpbar(x, quad_obj, positive_xy, res, quad_grad, positive_xy_jac, oc, ic)
   call check_vec('adpbar', x, [1.0_dp, 2.0_dp], 3.0e-3_dp, failures)

   ! Maximization with an equality constraint.
   x = [0.1_dp, 0.9_dp]
   oc = alabama_outer_control_t()
   oc%maximize = .true.
   oc%trace = .false.
   call auglag1(x, concave_obj, sum_one, res, concave_grad, sum_one_jac, oc, ic)
   call check_vec('auglag maximize', x, [0.5_dp, 0.5_dp], 5.0e-4_dp, failures)

   ! Strict-feasibility error for the legacy barrier.
   x = [-0.1_dp, 1.0_dp]
   call adpbar(x, quad_obj, positive_xy, res, control_outer=oc, control_inner=ic)
   call check_true('adpbar rejects infeasible start', res%convergence == al_invalid_input, failures)

   if (failures /= 0) then
      write(*,'(a,i0)') 'FAILURES: ', failures
      error stop 1
   end if
   print *, 'All alabama tests passed.'

contains

   function quad_obj(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = (x(1)-1.0_dp)**2 + (x(2)-2.0_dp)**2
   end function quad_obj

   subroutine quad_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g = [2.0_dp*(x(1)-1.0_dp), 2.0_dp*(x(2)-2.0_dp)]
   end subroutine quad_grad

   function boundary_obj(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = (x(1)+1.0_dp)**2 + (x(2)-2.0_dp)**2
   end function boundary_obj

   subroutine boundary_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g = [2.0_dp*(x(1)+1.0_dp), 2.0_dp*(x(2)-2.0_dp)]
   end subroutine boundary_grad

   function concave_obj(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = -(x(1)-1.0_dp)**2 - (x(2)-1.0_dp)**2
   end function concave_obj

   subroutine concave_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g = [-2.0_dp*(x(1)-1.0_dp), -2.0_dp*(x(2)-1.0_dp)]
   end subroutine concave_grad

   function documented_obj(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = (x(1)+3.0_dp*x(2)+x(3))**2 + 4.0_dp*(x(1)-x(2))**2
   end function documented_obj

   subroutine documented_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g(1) = 2.0_dp*(x(1)+3.0_dp*x(2)+x(3)) + 8.0_dp*(x(1)-x(2))
      g(2) = 6.0_dp*(x(1)+3.0_dp*x(2)+x(3)) - 8.0_dp*(x(1)-x(2))
      g(3) = 2.0_dp*(x(1)+3.0_dp*x(2)+x(3))
   end subroutine documented_grad

   function documented_heq(x) result(e)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: e(:)
      allocate(e(1)); e(1)=sum(x)-1.0_dp
   end function documented_heq

   function documented_heq_jac(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(1,size(x))); j=1.0_dp
   end function documented_heq_jac

   function documented_hin(x) result(h)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: h(:)
      allocate(h(4))
      h = [6.0_dp*x(2)+4.0_dp*x(3)-x(1)**3-3.0_dp, x(1), x(2), x(3)]
   end function documented_hin

   function documented_hin_jac(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(4,size(x))); j=0.0_dp
      j(1,:)=[-3.0_dp*x(1)**2,6.0_dp,4.0_dp]
      j(2,1)=1.0_dp; j(3,2)=1.0_dp; j(4,3)=1.0_dp
   end function documented_hin_jac

   function positive_xy(x) result(h)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: h(:)
      allocate(h(2)); h = x(1:2)
   end function positive_xy

   function positive_xy_jac(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(2,size(x))); j = 0.0_dp
      j(1,1)=1.0_dp; j(2,2)=1.0_dp
   end function positive_xy_jac

   function x_nonnegative(x) result(h)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: h(:)
      allocate(h(1)); h(1)=x(1)
   end function x_nonnegative

   function x_nonnegative_jac(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(1,size(x))); j=0.0_dp; j(1,1)=1.0_dp
   end function x_nonnegative_jac

   function sum_two(x) result(e)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: e(:)
      allocate(e(1)); e(1)=sum(x)-2.0_dp
   end function sum_two

   function sum_two_jac(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(1,size(x))); j=1.0_dp
   end function sum_two_jac

   function sum_one(x) result(e)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: e(:)
      allocate(e(1)); e(1)=sum(x)-1.0_dp
   end function sum_one

   function sum_one_jac(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(1,size(x))); j=1.0_dp
   end function sum_one_jac

   subroutine check_vec(name, got, expected, tol, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: got(:), expected(:), tol
      integer, intent(inout) :: failures
      if (maxval(abs(got-expected)) > tol) then
         write(*,'(a,2x,*(es15.7,1x))') 'FAIL '//trim(name)//':', got
         failures = failures + 1
      else
         write(*,'(a)') 'PASS '//trim(name)
      end if
   end subroutine check_vec

   subroutine check_true(name, ok, failures)
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      integer, intent(inout) :: failures
      if (.not. ok) then
         write(*,'(a)') 'FAIL '//trim(name)
         failures = failures + 1
      else
         write(*,'(a)') 'PASS '//trim(name)
      end if
   end subroutine check_true

end program test_alabama
