module optr_iterative
   use optr_kinds, only : dp
   use optr_status, only : optr_success, optr_invalid_shape, optr_singular, &
      optr_no_convergence, optr_invalid_argument
   implicit none
   private

   public :: conjugate_gradient
   public :: gauss_seidel
   public :: non_diag_multiplication

contains

   pure function non_diag_multiplication(i, a, beta) result(asum)
      integer, intent(in) :: i !! One-based row index whose off-diagonal contribution is requested.
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: beta(:) !! Current iterate vector with one element per matrix row.
      real(dp) :: asum

      integer :: j

      asum = 0.0_dp
      do j = 1, size(beta)
         if (j /= i) asum = asum + a(i, j) * beta(j)
      end do
   end function non_diag_multiplication

   pure subroutine gauss_seidel(a, b, x, status, iterations, tol, max_iter, relaxation, relax_after, x0)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A x = b`.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), allocatable, intent(out) :: x(:) !! Final iterate, converged when `status` is zero.
      integer, intent(out) :: status !! Zero on convergence; nonzero for bad input, singular diagonal, or iteration limit.
      integer, intent(out) :: iterations !! Number of iterations actually performed.
      real(dp), optional, intent(in) :: tol !! Euclidean step tolerance; defaults to `1.0e-7_dp`.
      integer, optional, intent(in) :: max_iter !! Maximum iterations; defaults to 500.
      real(dp), optional, intent(in) :: relaxation !! Upstream multiplicative relaxation factor `w`; defaults to 1.
      integer, optional, intent(in) :: relax_after !! Iteration used to estimate a new relaxation factor; defaults to `max_iter`.
      real(dp), optional, intent(in) :: x0(:) !! Initial iterate; defaults to a zero vector.

      real(dp), allocatable :: x_old(:)
      real(dp) :: dx
      real(dp) :: dx1
      real(dp) :: dx2
      real(dp) :: factor
      real(dp) :: local_tol
      real(dp) :: ratio
      integer :: i
      integer :: k
      integer :: local_max_iter
      integer :: local_relax_after
      integer :: n

      status = optr_success
      iterations = 0
      if (size(a, 1) /= size(a, 2) .or. size(b) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0))
         return
      end if

      n = size(b)
      if (present(x0)) then
         if (size(x0) /= n) then
            status = optr_invalid_shape
            allocate(x(0))
            return
         end if
      end if

      local_tol = 1.0e-7_dp
      if (present(tol)) local_tol = max(tol, 0.0_dp)
      local_max_iter = 500
      if (present(max_iter)) local_max_iter = max_iter
      if (local_max_iter < 1) then
         status = optr_invalid_argument
         allocate(x(0))
         return
      end if

      factor = 1.0_dp
      if (present(relaxation)) factor = relaxation
      local_relax_after = local_max_iter
      if (present(relax_after)) local_relax_after = relax_after
      if (local_relax_after < 1) then
         status = optr_invalid_argument
         allocate(x(0))
         return
      end if

      allocate(x(n), x_old(n))
      x = 0.0_dp
      if (present(x0)) x = x0
      dx1 = 0.0_dp

      do i = 1, local_max_iter
         x_old = x
         do k = 1, n
            if (abs(a(k, k)) <= tiny(1.0_dp)) then
               status = optr_singular
               iterations = i
               return
            end if
            x(k) = (factor / a(k, k)) * (b(k) - non_diag_multiplication(k, a, x))
         end do
         dx = sqrt(sum((x - x_old) ** 2))
         iterations = i
         if (dx <= local_tol) then
            status = optr_success
            return
         end if

         if (i == local_relax_after) dx1 = dx
         if (i == local_relax_after + 1 .and. dx1 > 0.0_dp) then
            dx2 = dx
            ratio = dx2 / dx1
            if (ratio >= 0.0_dp .and. ratio <= 1.0_dp) then
               factor = 2.0_dp / (1.0_dp + sqrt(1.0_dp - ratio))
            end if
         end if
      end do

      status = optr_no_convergence
   end subroutine gauss_seidel

   pure subroutine conjugate_gradient(a, b, x, status, iterations, convergence, tol, max_iter, x0)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix, normally symmetric positive definite.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), allocatable, intent(out) :: x(:) !! Final iterate, converged when `status` is zero.
      integer, intent(out) :: status !! Zero on convergence; nonzero for bad input, breakdown, or iteration limit.
      integer, intent(out) :: iterations !! Number of iterations actually performed.
      real(dp), allocatable, intent(out) :: convergence(:) !! Euclidean step norms for the iterations that were completed.
      real(dp), optional, intent(in) :: tol !! Euclidean step tolerance; defaults to `1.0e-7_dp`.
      integer, optional, intent(in) :: max_iter !! Maximum iterations; defaults to 500.
      real(dp), optional, intent(in) :: x0(:) !! Initial iterate; defaults to a zero vector.

      real(dp), allocatable :: a_direction(:)
      real(dp), allocatable :: history(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: search_direction(:)
      real(dp), allocatable :: x_old(:)
      real(dp) :: alpha
      real(dp) :: alpha_parameter
      real(dp) :: beta
      real(dp) :: denominator
      real(dp) :: dx
      real(dp) :: local_tol
      integer :: i
      integer :: local_max_iter
      integer :: n

      status = optr_success
      iterations = 0
      if (size(a, 1) /= size(a, 2) .or. size(b) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0), convergence(0))
         return
      end if

      n = size(b)
      if (present(x0)) then
         if (size(x0) /= n) then
            status = optr_invalid_shape
            allocate(x(0), convergence(0))
            return
         end if
      end if

      local_tol = 1.0e-7_dp
      if (present(tol)) local_tol = max(tol, 0.0_dp)
      local_max_iter = 500
      if (present(max_iter)) local_max_iter = max_iter
      if (local_max_iter < 1) then
         status = optr_invalid_argument
         allocate(x(0), convergence(0))
         return
      end if

      allocate(x(n), x_old(n), residual(n), search_direction(n), a_direction(n), history(local_max_iter))
      x = 0.0_dp
      if (present(x0)) x = x0
      search_direction = 0.0_dp
      history = 0.0_dp
      beta = 0.0_dp

      do i = 1, local_max_iter
         x_old = x
         residual = b - matmul(a, x)
         if (i > 1) then
            a_direction = matmul(a, search_direction)
            denominator = dot_product(search_direction, a_direction)
            if (abs(denominator) <= tiny(1.0_dp)) then
               status = optr_singular
               iterations = i
               allocate(convergence(i - 1))
               if (i > 1) convergence = history(1:i - 1)
               return
            end if
            beta = -dot_product(residual, a_direction) / denominator
         end if

         search_direction = residual + beta * search_direction
         alpha_parameter = dot_product(search_direction, residual)
         if (abs(alpha_parameter) <= tiny(1.0_dp) .and. i > 1) then
            iterations = i
            allocate(convergence(i - 1))
            convergence = history(1:i - 1)
            status = optr_success
            return
         end if

         a_direction = matmul(a, search_direction)
         denominator = dot_product(search_direction, a_direction)
         if (abs(denominator) <= tiny(1.0_dp)) then
            status = optr_singular
            iterations = i
            allocate(convergence(i - 1))
            if (i > 1) convergence = history(1:i - 1)
            return
         end if
         alpha = alpha_parameter / denominator
         x = x + alpha * search_direction
         dx = sqrt(sum((x - x_old) ** 2))
         history(i) = dx
         iterations = i
         if (dx <= local_tol) then
            allocate(convergence(i))
            convergence = history(1:i)
            status = optr_success
            return
         end if
      end do

      allocate(convergence(local_max_iter))
      convergence = history
      status = optr_no_convergence
   end subroutine conjugate_gradient

end module optr_iterative
