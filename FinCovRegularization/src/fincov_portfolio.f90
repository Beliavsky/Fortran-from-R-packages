! SPDX-License-Identifier: GPL-2.0-only
module fincov_portfolio
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fincov_kinds, only : dp
   use fincov_status, only : fincov_ok, fincov_invalid_input, fincov_size_mismatch, &
      fincov_singular_matrix, fincov_no_convergence
   use fincov_linalg, only : solve_linear_system, matrix_is_symmetric
   use fincov_utils, only : vector_norm2
   implicit none
   private

   public :: gmvp, risk_parity, risk_parity_objective
contains
   function gmvp(covariance, allow_short, status, tolerance, max_iterations) result(weights)
      real(dp), intent(in) :: covariance(:,:)
      logical, intent(in), optional :: allow_short
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp) :: weights(size(covariance,1))
      logical :: shorts
      logical, allocatable :: active(:)
      integer, allocatable :: index(:)
      real(dp), allocatable :: rhs(:), solution(:), submatrix(:,:), active_weights(:), gradient(:)
      real(dp) :: denominator, tol, lambda, violation, most_negative
      integer :: n, m, i, j, iteration, max_iter, local_status, remove_position, add_index

      n = size(covariance,1)
      weights = ieee_value(0.0_dp, ieee_quiet_nan)
      shorts = .true.
      if (present(allow_short)) shorts = allow_short
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = max(tolerance, 100.0_dp*epsilon(1.0_dp))
      max_iter = max(20, 10*n + 20)
      if (present(max_iterations)) max_iter = max(1, max_iterations)

      if (n < 1 .or. size(covariance,2) /= n) then
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      if (.not. matrix_is_symmetric(covariance)) then
         if (present(status)) status = fincov_invalid_input
         return
      end if

      allocate(rhs(n))
      rhs = 1.0_dp
      if (shorts) then
         call solve_linear_system(covariance, rhs, solution, local_status)
         if (local_status /= fincov_ok) then
            if (present(status)) status = local_status
            return
         end if
         denominator = sum(solution)
         if (abs(denominator) <= 100.0_dp*epsilon(1.0_dp)) then
            if (present(status)) status = fincov_singular_matrix
            return
         end if
         weights = solution/denominator
         if (present(status)) status = fincov_ok
         return
      end if

      allocate(active(n), gradient(n))
      active = .true.
      weights = 1.0_dp/real(n,dp)
      do iteration = 1, max_iter
         m = count(active)
         if (m < 1) then
            weights = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = fincov_no_convergence
            return
         end if
         if (allocated(index)) deallocate(index)
         if (allocated(submatrix)) deallocate(submatrix)
         if (allocated(active_weights)) deallocate(active_weights)
         if (allocated(solution)) deallocate(solution)
         if (allocated(rhs)) deallocate(rhs)
         allocate(index(m), submatrix(m,m), active_weights(m), rhs(m))
         j = 0
         do i = 1, n
            if (active(i)) then
               j = j + 1
               index(j) = i
            end if
         end do
         do i = 1, m
            do j = 1, m
               submatrix(i,j) = covariance(index(i),index(j))
            end do
         end do
         rhs = 1.0_dp
         call solve_linear_system(submatrix, rhs, solution, local_status)
         if (local_status /= fincov_ok) then
            weights = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = local_status
            return
         end if
         denominator = sum(solution)
         if (abs(denominator) <= 100.0_dp*epsilon(1.0_dp)) then
            weights = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = fincov_singular_matrix
            return
         end if
         active_weights = solution/denominator

         most_negative = minval(active_weights)
         if (most_negative < -tol .and. m > 1) then
            remove_position = minloc(active_weights, dim=1)
            active(index(remove_position)) = .false.
            cycle
         end if

         weights = 0.0_dp
         do i = 1, m
            weights(index(i)) = max(active_weights(i), 0.0_dp)
         end do
         if (sum(weights) > 0.0_dp) weights = weights/sum(weights)
         gradient = matmul(covariance, weights)
         lambda = dot_product(weights, gradient)
         violation = 0.0_dp
         add_index = 0
         do i = 1, n
            if (.not. active(i)) then
               if (gradient(i) - lambda < violation) then
                  violation = gradient(i) - lambda
                  add_index = i
               end if
            end if
         end do
         if (add_index == 0 .or. violation >= -tol) then
            if (present(status)) status = fincov_ok
            return
         end if
         active(add_index) = .true.
      end do

      if (present(status)) status = fincov_no_convergence
   end function gmvp

   function risk_parity(covariance, status, tolerance, max_iterations, objective_value) result(weights)
      real(dp), intent(in) :: covariance(:,:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), intent(out), optional :: objective_value
      real(dp) :: weights(size(covariance,1))
      real(dp), allocatable :: simplex(:,:), values(:), centroid(:), reflected(:), expanded(:), contracted(:)
      real(dp), allocatable :: temp_point(:), x(:)
      real(dp) :: tol, alpha, gamma, rho, shrink, reflected_value, expanded_value, contracted_value
      real(dp) :: value_spread, point_spread, tmp_value
      integer :: n, m, i, j, best, iteration, max_iter

      n = size(covariance,1)
      weights = ieee_value(0.0_dp, ieee_quiet_nan)
      if (n < 1 .or. size(covariance,2) /= n) then
         if (present(status)) status = fincov_size_mismatch
         if (present(objective_value)) objective_value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (.not. matrix_is_symmetric(covariance)) then
         if (present(status)) status = fincov_invalid_input
         if (present(objective_value)) objective_value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (n == 1) then
         weights(1) = 1.0_dp
         if (present(status)) status = fincov_ok
         if (present(objective_value)) objective_value = 0.0_dp
         return
      end if

      m = n - 1
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = max(tolerance, 100.0_dp*epsilon(1.0_dp))
      max_iter = max(1000, 300*m)
      if (present(max_iterations)) max_iter = max(1, max_iterations)
      alpha = 1.0_dp
      gamma = 2.0_dp
      rho = 0.5_dp
      shrink = 0.5_dp

      allocate(simplex(m,m+1), values(m+1), centroid(m), reflected(m), expanded(m), contracted(m), &
         temp_point(m), x(m))
      simplex = 0.0_dp
      do i = 1, m
         simplex(i,i+1) = 0.05_dp
      end do
      do i = 1, m + 1
         values(i) = objective_from_reduced(simplex(:,i), covariance)
      end do

      do iteration = 1, max_iter
         do i = 1, m
            best = i
            do j = i + 1, m + 1
               if (values(j) < values(best)) best = j
            end do
            if (best /= i) then
               tmp_value = values(i)
               values(i) = values(best)
               values(best) = tmp_value
               temp_point = simplex(:,i)
               simplex(:,i) = simplex(:,best)
               simplex(:,best) = temp_point
            end if
         end do

         value_spread = maxval(abs(values - values(1)))
         point_spread = 0.0_dp
         do i = 2, m + 1
            point_spread = max(point_spread, vector_norm2(simplex(:,i) - simplex(:,1)))
         end do
         if (value_spread <= tol*max(1.0_dp,abs(values(1))) .and. point_spread <= max(sqrt(epsilon(1.0_dp)), 10.0_dp*tol)) exit

         centroid = sum(simplex(:,1:m), dim=2)/real(m,dp)
         reflected = centroid + alpha*(centroid - simplex(:,m+1))
         reflected_value = objective_from_reduced(reflected, covariance)

         if (reflected_value < values(1)) then
            expanded = centroid + gamma*(reflected - centroid)
            expanded_value = objective_from_reduced(expanded, covariance)
            if (expanded_value < reflected_value) then
               simplex(:,m+1) = expanded
               values(m+1) = expanded_value
            else
               simplex(:,m+1) = reflected
               values(m+1) = reflected_value
            end if
         else if (reflected_value < values(m)) then
            simplex(:,m+1) = reflected
            values(m+1) = reflected_value
         else
            if (reflected_value < values(m+1)) then
               contracted = centroid + rho*(reflected - centroid)
            else
               contracted = centroid + rho*(simplex(:,m+1) - centroid)
            end if
            contracted_value = objective_from_reduced(contracted, covariance)
            if (contracted_value < min(reflected_value, values(m+1))) then
               simplex(:,m+1) = contracted
               values(m+1) = contracted_value
            else
               do i = 2, m + 1
                  simplex(:,i) = simplex(:,1) + shrink*(simplex(:,i) - simplex(:,1))
                  values(i) = objective_from_reduced(simplex(:,i), covariance)
               end do
            end if
         end if
      end do

      best = minloc(values, dim=1)
      x = simplex(:,best)
      weights(1:m) = x
      weights(n) = 1.0_dp - sum(x)
      if (present(objective_value)) objective_value = risk_parity_objective(weights, covariance)
      if (present(status)) then
         if (iteration > max_iter) then
            status = fincov_no_convergence
         else
            status = fincov_ok
         end if
      end if
   end function risk_parity

   pure function risk_parity_objective(weights, covariance) result(value)
      real(dp), intent(in) :: weights(:), covariance(:,:)
      real(dp) :: value
      real(dp) :: contributions(size(weights))
      integer :: n

      n = size(weights)
      if (n < 1 .or. size(covariance,1) /= n .or. size(covariance,2) /= n) then
         value = huge(1.0_dp)
         return
      end if
      contributions = weights*matmul(covariance, weights)
      value = 2.0_dp*real(n,dp)*sum(contributions*contributions) - 2.0_dp*sum(contributions)**2
   end function risk_parity_objective

   pure function objective_from_reduced(x, covariance) result(value)
      real(dp), intent(in) :: x(:), covariance(:,:)
      real(dp) :: value
      real(dp) :: weights(size(x)+1)

      weights(1:size(x)) = x
      weights(size(weights)) = 1.0_dp - sum(x)
      value = risk_parity_objective(weights, covariance)
   end function objective_from_reduced
end module fincov_portfolio
