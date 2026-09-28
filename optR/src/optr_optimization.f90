module optr_optimization
   use optr_kinds, only : dp
   use optr_status, only : optr_success, optr_invalid_shape, optr_no_convergence, &
      optr_invalid_argument
   use optr_linear, only : gaussian_solve
   implicit none
   private

   public :: finite_difference_jacobian
   public :: newton_raphson
   public :: vector_function

   abstract interface
      pure subroutine vector_evaluate_iface(x, fx)
         import :: dp
         real(dp), intent(in) :: x(:) !! Point at which the vector-valued function is evaluated.
         real(dp), intent(out) :: fx(:) !! Function values; caller supplies the expected output shape.
      end subroutine vector_evaluate_iface
   end interface

   type :: vector_function
      procedure(vector_evaluate_iface), pointer, nopass :: evaluate => null()
   end type vector_function

contains

   pure subroutine finite_difference_jacobian(f, x, jacobian_matrix, f0, status, step)
      type(vector_function), intent(in) :: f !! Vector-valued function wrapper with an associated pure evaluator.
      real(dp), intent(in) :: x(:) !! Point at which the Jacobian is estimated; output dimension is assumed equal to `size(x)`.
      real(dp), allocatable, intent(out) :: jacobian_matrix(:, :) !! Forward-difference Jacobian with columns for perturbed inputs.
      real(dp), allocatable, intent(out) :: f0(:) !! Function value at the unperturbed point.
      integer, intent(out) :: status !! Zero on success; nonzero if the callback is absent or the input is empty.
      real(dp), optional, intent(in) :: step !! Positive forward-difference step; defaults to `1.0e-4_dp`.

      real(dp), allocatable :: f1(:)
      real(dp), allocatable :: x_work(:)
      real(dp) :: h
      integer :: i
      integer :: n

      if (.not. associated(f%evaluate)) then
         status = optr_invalid_argument
         allocate(jacobian_matrix(0, 0), f0(0))
         return
      end if

      n = size(x)
      if (n < 1) then
         status = optr_invalid_shape
         allocate(jacobian_matrix(0, 0), f0(0))
         return
      end if
      h = 1.0e-4_dp
      if (present(step)) h = step
      if (h <= 0.0_dp) then
         status = optr_invalid_argument
         allocate(jacobian_matrix(0, 0), f0(0))
         return
      end if

      allocate(jacobian_matrix(n, n), f0(n), f1(n), x_work(n))
      call f%evaluate(x, f0)
      x_work = x
      do i = 1, n
         x_work(i) = x(i) + h
         call f%evaluate(x_work, f1)
         x_work(i) = x(i)
         jacobian_matrix(:, i) = (f1 - f0) / h
      end do
      status = optr_success
   end subroutine finite_difference_jacobian

   pure subroutine newton_raphson(f, x0, x, status, iterations, max_iter, tol, step)
      type(vector_function), intent(in) :: f !! Vector-valued function wrapper defining the nonlinear system `f(x)=0`.
      real(dp), intent(in) :: x0(:) !! Initial root estimate.
      real(dp), allocatable, intent(out) :: x(:) !! Final root estimate.
      integer, intent(out) :: status !! Zero on convergence; nonzero for callback/linear-solve failure or iteration limit.
      integer, intent(out) :: iterations !! Number of Newton iterations performed.
      integer, optional, intent(in) :: max_iter !! Maximum Newton iterations; defaults to 30.
      real(dp), optional, intent(in) :: tol !! Residual and step tolerance; defaults to `1.0e-9_dp`.
      real(dp), optional, intent(in) :: step !! Forward-difference Jacobian step; defaults to `1.0e-4_dp`.

      real(dp), allocatable :: dx(:)
      real(dp), allocatable :: f0(:)
      real(dp), allocatable :: jac(:, :)
      real(dp) :: h
      real(dp) :: local_tol
      real(dp) :: residual_rms
      real(dp) :: scale
      integer :: i
      integer :: linear_status
      integer :: local_max_iter
      integer :: jac_status
      integer :: n

      iterations = 0
      n = size(x0)
      if (n < 1) then
         status = optr_invalid_shape
         allocate(x(0))
         return
      end if
      if (.not. associated(f%evaluate)) then
         status = optr_invalid_argument
         allocate(x(0))
         return
      end if

      local_max_iter = 30
      if (present(max_iter)) local_max_iter = max_iter
      if (local_max_iter < 1) then
         status = optr_invalid_argument
         allocate(x(0))
         return
      end if
      local_tol = 1.0e-9_dp
      if (present(tol)) local_tol = max(tol, 0.0_dp)
      h = 1.0e-4_dp
      if (present(step)) h = step
      if (h <= 0.0_dp) then
         status = optr_invalid_argument
         allocate(x(0))
         return
      end if

      allocate(x(n))
      x = x0
      do i = 1, local_max_iter
         call finite_difference_jacobian(f, x, jac, f0, jac_status, h)
         if (jac_status /= optr_success) then
            status = jac_status
            iterations = i
            return
         end if

         residual_rms = sqrt(dot_product(f0, f0) / real(n, dp))
         if (residual_rms < local_tol) then
            status = optr_success
            iterations = i - 1
            return
         end if

         call gaussian_solve(jac, -f0, dx, linear_status, 1.0e-12_dp)
         if (linear_status /= optr_success) then
            status = linear_status
            iterations = i
            return
         end if
         x = x + dx
         iterations = i
         scale = max(1.0_dp, maxval(abs(x)))
         if (sqrt(dot_product(dx, dx)) < local_tol * scale) then
            status = optr_success
            return
         end if
      end do

      status = optr_no_convergence
   end subroutine newton_raphson

end module optr_optimization
