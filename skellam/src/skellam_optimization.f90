! SPDX-License-Identifier: GPL-2.0-or-later
module skellam_optimization
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use skellam_kinds, only : dp
   implicit none
   private

   type, public :: optimizer_result
      real(dp), allocatable :: parameters(:)
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 0
      logical :: converged = .false.
   end type optimizer_result

   abstract interface
      function objective_function(parameters) result(value)
         import :: dp
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value
      end function objective_function
   end interface

   public :: minimize_bfgs, minimize_nelder_mead
   public :: numerical_gradient, numerical_hessian
   public :: invert_matrix, solve_linear_system

contains

   subroutine minimize_bfgs(objective, initial, result, tolerance, max_iterations)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: initial(:)
      type(optimizer_result), intent(out) :: result
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: x(:), x_new(:), gradient(:), gradient_new(:)
      real(dp), allocatable :: inverse_hessian(:,:), direction(:), step_vector(:), y(:)
      real(dp), allocatable :: identity(:,:), outer_sy(:,:), outer_ys(:,:)
      real(dp) :: f, f_new, step, armijo, ys, rho, tol, slope
      integer :: n, iteration, max_iter, i, line_iteration
      logical :: accepted

      n = size(initial)
      tol = 1.0e-7_dp
      max_iter = 500
      if (present(tolerance)) tol = tolerance
      if (present(max_iterations)) max_iter = max_iterations

      allocate(x(n), x_new(n), gradient(n), gradient_new(n), inverse_hessian(n,n))
      allocate(direction(n), step_vector(n), y(n), identity(n,n), outer_sy(n,n), outer_ys(n,n))
      x = initial
      identity = 0.0_dp
      do i = 1, n
         identity(i,i) = 1.0_dp
      end do
      inverse_hessian = identity
      f = objective(x)
      result%evaluations = 1
      call numerical_gradient(objective, x, gradient, result%evaluations)

      do iteration = 1, max_iter
         result%iterations = iteration
         if (maxval(abs(gradient)) <= tol*(1.0_dp + abs(f))) then
            result%converged = .true.
            result%status = 0
            exit
         end if

         direction = -matmul(inverse_hessian, gradient)
         slope = dot_product(gradient, direction)
         if (slope >= -sqrt(epsilon(1.0_dp))*max(1.0_dp, norm2(gradient)*norm2(direction))) then
            direction = -gradient
            inverse_hessian = identity
            slope = -dot_product(gradient, gradient)
         end if

         step = 1.0_dp
         armijo = 1.0e-4_dp
         accepted = .false.
         do line_iteration = 1, 40
            x_new = x + step*direction
            f_new = objective(x_new)
            result%evaluations = result%evaluations + 1
            if (ieee_is_finite(f_new)) then
               if (f_new <= f + armijo*step*slope) then
                  accepted = .true.
                  exit
               end if
            end if
            step = 0.5_dp*step
         end do
         if (.not. accepted) then
            result%status = 2
            exit
         end if

         call numerical_gradient(objective, x_new, gradient_new, result%evaluations)
         step_vector = x_new - x
         y = gradient_new - gradient
         ys = dot_product(y, step_vector)
         if (ys > sqrt(epsilon(1.0_dp))*norm2(y)*norm2(step_vector)) then
            rho = 1.0_dp/ys
            outer_sy = outer_product(step_vector, y)
            outer_ys = outer_product(y, step_vector)
            inverse_hessian = matmul(identity - rho*outer_sy, &
               matmul(inverse_hessian, identity - rho*outer_ys)) &
               + rho*outer_product(step_vector, step_vector)
            inverse_hessian = 0.5_dp*(inverse_hessian + transpose(inverse_hessian))
         else
            inverse_hessian = identity
         end if

         x = x_new
         f = f_new
         gradient = gradient_new
         if (norm2(step_vector) <= tol*(1.0_dp + norm2(x))) then
            result%converged = .true.
            result%status = 0
            exit
         end if
      end do

      if (.not. result%converged .and. result%status == 0) result%status = 1
      result%parameters = x
      result%objective = f
   end subroutine minimize_bfgs

   subroutine minimize_nelder_mead(objective, initial, result, tolerance, max_iterations)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: initial(:)
      type(optimizer_result), intent(out) :: result
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: simplex(:,:), values(:), centroid(:), reflected(:), expanded(:), contracted(:)
      real(dp) :: tol, scale, reflected_value, expanded_value, contracted_value
      integer :: n, i, iteration, max_iter

      n = size(initial)
      tol = 1.0e-8_dp
      max_iter = 2000
      if (present(tolerance)) tol = tolerance
      if (present(max_iterations)) max_iter = max_iterations
      allocate(simplex(n,n + 1), values(n + 1), centroid(n), reflected(n), expanded(n), contracted(n))

      simplex(:,1) = initial
      do i = 1, n
         simplex(:,i + 1) = initial
         scale = 0.05_dp*max(1.0_dp, abs(initial(i)))
         simplex(i,i + 1) = simplex(i,i + 1) + scale
      end do
      do i = 1, n + 1
         values(i) = objective(simplex(:,i))
      end do
      result%evaluations = n + 1

      do iteration = 1, max_iter
         result%iterations = iteration
         call sort_simplex(simplex, values)
         if (maxval(abs(values - values(1))) <= tol*(1.0_dp + abs(values(1))) .and. &
             maxval(abs(simplex - spread(simplex(:,1), 2, n + 1))) <= tol*(1.0_dp + maxval(abs(simplex(:,1))))) then
            result%converged = .true.
            result%status = 0
            exit
         end if

         centroid = sum(simplex(:,1:n), dim=2)/real(n, dp)
         reflected = centroid + (centroid - simplex(:,n + 1))
         reflected_value = objective(reflected)
         result%evaluations = result%evaluations + 1

         if (reflected_value < values(1)) then
            expanded = centroid + 2.0_dp*(reflected - centroid)
            expanded_value = objective(expanded)
            result%evaluations = result%evaluations + 1
            if (expanded_value < reflected_value) then
               simplex(:,n + 1) = expanded
               values(n + 1) = expanded_value
            else
               simplex(:,n + 1) = reflected
               values(n + 1) = reflected_value
            end if
         else if (reflected_value < values(n)) then
            simplex(:,n + 1) = reflected
            values(n + 1) = reflected_value
         else
            if (reflected_value < values(n + 1)) then
               contracted = centroid + 0.5_dp*(reflected - centroid)
            else
               contracted = centroid + 0.5_dp*(simplex(:,n + 1) - centroid)
            end if
            contracted_value = objective(contracted)
            result%evaluations = result%evaluations + 1
            if (contracted_value < min(reflected_value, values(n + 1))) then
               simplex(:,n + 1) = contracted
               values(n + 1) = contracted_value
            else
               do i = 2, n + 1
                  simplex(:,i) = simplex(:,1) + 0.5_dp*(simplex(:,i) - simplex(:,1))
                  values(i) = objective(simplex(:,i))
               end do
               result%evaluations = result%evaluations + n
            end if
         end if
      end do

      call sort_simplex(simplex, values)
      if (.not. result%converged) result%status = 1
      result%parameters = simplex(:,1)
      result%objective = values(1)
   end subroutine minimize_nelder_mead

   subroutine numerical_gradient(objective, parameters, gradient, evaluations)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: parameters(:)
      real(dp), intent(out) :: gradient(:)
      integer, intent(inout), optional :: evaluations
      real(dp), allocatable :: plus(:), minus(:)
      real(dp) :: step, fplus, fminus
      integer :: i

      allocate(plus(size(parameters)), minus(size(parameters)))
      do i = 1, size(parameters)
         step = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(parameters(i)))
         plus = parameters
         minus = parameters
         plus(i) = plus(i) + step
         minus(i) = minus(i) - step
         fplus = objective(plus)
         fminus = objective(minus)
         gradient(i) = (fplus - fminus)/(2.0_dp*step)
         if (present(evaluations)) evaluations = evaluations + 2
      end do
   end subroutine numerical_gradient

   subroutine numerical_hessian(objective, parameters, hessian, status)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: parameters(:)
      real(dp), intent(out) :: hessian(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: plus(:), minus(:), gradient_plus(:), gradient_minus(:)
      real(dp) :: step
      integer :: i

      if (present(status)) status = 0
      if (size(hessian,1) /= size(parameters) .or. size(hessian,2) /= size(parameters)) then
         if (present(status)) status = 1
         return
      end if
      allocate(plus(size(parameters)), minus(size(parameters)))
      allocate(gradient_plus(size(parameters)), gradient_minus(size(parameters)))
      do i = 1, size(parameters)
         step = epsilon(1.0_dp)**0.25_dp*max(1.0_dp, abs(parameters(i)))
         plus = parameters
         minus = parameters
         plus(i) = plus(i) + step
         minus(i) = minus(i) - step
         call numerical_gradient(objective, plus, gradient_plus)
         call numerical_gradient(objective, minus, gradient_minus)
         hessian(:,i) = (gradient_plus - gradient_minus)/(2.0_dp*step)
      end do
      hessian = 0.5_dp*(hessian + transpose(hessian))
   end subroutine numerical_hessian

   subroutine solve_linear_system(matrix, rhs, solution, status)
      real(dp), intent(in) :: matrix(:,:), rhs(:)
      real(dp), intent(out) :: solution(:)
      integer, intent(out) :: status
      real(dp), allocatable :: augmented(:,:)
      real(dp) :: pivot_value, factor
      integer :: n, i, j, pivot_row

      n = size(rhs)
      status = 0
      if (size(matrix,1) /= n .or. size(matrix,2) /= n .or. size(solution) /= n) then
         status = 1
         return
      end if
      allocate(augmented(n,n + 1))
      augmented(:,1:n) = matrix
      augmented(:,n + 1) = rhs

      do i = 1, n
         pivot_row = i - 1 + maxloc(abs(augmented(i:n,i)), dim=1)
         if (abs(augmented(pivot_row,i)) <= epsilon(1.0_dp)*max(1.0_dp, maxval(abs(matrix)))) then
            status = 2
            return
         end if
         if (pivot_row /= i) call swap_rows(augmented, i, pivot_row)
         pivot_value = augmented(i,i)
         augmented(i,i:n + 1) = augmented(i,i:n + 1)/pivot_value
         do j = 1, n
            if (j == i) cycle
            factor = augmented(j,i)
            augmented(j,i:n + 1) = augmented(j,i:n + 1) - factor*augmented(i,i:n + 1)
         end do
      end do
      solution = augmented(:,n + 1)
   end subroutine solve_linear_system

   subroutine invert_matrix(matrix, inverse, status)
      real(dp), intent(in) :: matrix(:,:)
      real(dp), intent(out) :: inverse(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: augmented(:,:)
      real(dp) :: pivot_value, factor
      integer :: n, i, j, pivot_row

      n = size(matrix,1)
      status = 0
      if (size(matrix,2) /= n .or. size(inverse,1) /= n .or. size(inverse,2) /= n) then
         status = 1
         return
      end if
      allocate(augmented(n,2*n))
      augmented = 0.0_dp
      augmented(:,1:n) = matrix
      do i = 1, n
         augmented(i,n + i) = 1.0_dp
      end do

      do i = 1, n
         pivot_row = i - 1 + maxloc(abs(augmented(i:n,i)), dim=1)
         if (abs(augmented(pivot_row,i)) <= epsilon(1.0_dp)*max(1.0_dp, maxval(abs(matrix)))) then
            status = 2
            inverse = 0.0_dp
            return
         end if
         if (pivot_row /= i) call swap_rows(augmented, i, pivot_row)
         pivot_value = augmented(i,i)
         augmented(i,:) = augmented(i,:)/pivot_value
         do j = 1, n
            if (j == i) cycle
            factor = augmented(j,i)
            augmented(j,:) = augmented(j,:) - factor*augmented(i,:)
         end do
      end do
      inverse = augmented(:,n + 1:2*n)
   end subroutine invert_matrix

   pure function outer_product(left, right) result(matrix)
      real(dp), intent(in) :: left(:), right(:)
      real(dp) :: matrix(size(left),size(right))
      matrix = spread(left, 2, size(right))*spread(right, 1, size(left))
   end function outer_product

   subroutine sort_simplex(simplex, values)
      real(dp), intent(inout) :: simplex(:,:), values(:)
      real(dp) :: value_temp
      real(dp), allocatable :: column_temp(:)
      integer :: i, best

      allocate(column_temp(size(simplex,1)))
      do i = 1, size(values) - 1
         best = i - 1 + minloc(values(i:), dim=1)
         if (best /= i) then
            value_temp = values(i)
            values(i) = values(best)
            values(best) = value_temp
            column_temp = simplex(:,i)
            simplex(:,i) = simplex(:,best)
            simplex(:,best) = column_temp
         end if
      end do
   end subroutine sort_simplex

   subroutine swap_rows(matrix, row1, row2)
      real(dp), intent(inout) :: matrix(:,:)
      integer, intent(in) :: row1, row2
      real(dp), allocatable :: temporary(:)

      allocate(temporary(size(matrix,2)))
      temporary = matrix(row1,:)
      matrix(row1,:) = matrix(row2,:)
      matrix(row2,:) = temporary
   end subroutine swap_rows

end module skellam_optimization
