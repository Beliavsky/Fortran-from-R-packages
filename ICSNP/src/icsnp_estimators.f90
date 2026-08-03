! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp_estimators
   use icsnp_kinds, only : dp
   use icsnp_status, only : icsnp_ok, icsnp_invalid_input, icsnp_singular, &
      icsnp_iteration_limit, icsnp_numerical_error
   use icsnp_types, only : location_scatter_result, spatial_sign_result
   use icsnp_linalg, only : invert_matrix, determinant, matrix_sqrt, matrix_inv_sqrt, &
      covariance_matrix, mahalanobis_squared, frobenius_norm, identity_matrix, sample_mean
   use icsnp_special, only : median_value, rank_average, normal_quantile, chi_square_quantile, &
      chi_square_survival
   use icsnp_pairs, only : pair_diff, pair_sum, pair_prod
   implicit none
   private
   public :: spatial_median, spatial_sign, tyler_shape, duembgen_shape
   public :: duembgen_shape_wt, symm_huber, symm_huber_wt
   public :: HR_Mest, HP1_shape, hl_loc, vdw_loc

contains

   subroutine spatial_median(x, center, status, iterations, init, maxiter, tolerance)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: center(:)
      integer, intent(out) :: status
      integer, intent(out), optional :: iterations
      real(dp), intent(in), optional :: init(:), tolerance
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: current(:), updated(:), distances(:), direction_sum(:)
      real(dp) :: tol, difference, inverse_sum, residual_norm, eta_ratio
      integer :: n, p, limit, iter, i, nonzero, coincident
      logical, allocatable :: active(:)

      n = size(x, 1)
      p = size(x, 2)
      allocate(center(0))
      status = icsnp_invalid_input
      if (present(iterations)) iterations = 0
      if (n < 1 .or. p < 1) return
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      limit = 500
      if (present(maxiter)) limit = maxiter
      allocate(current(p), updated(p), distances(n), active(n), direction_sum(p))
      if (present(init)) then
         if (size(init) /= p) return
         current = init
      else
         do i = 1, p
            current(i) = median_value(x(:, i))
         end do
      end if

      status = icsnp_iteration_limit
      do iter = 1, limit
         do i = 1, n
            distances(i) = sqrt(dot_product(x(i, :) - current, x(i, :) - current))
         end do
         active = distances > sqrt(tiny(1.0_dp))
         nonzero = count(active)
         coincident = n - nonzero
         if (nonzero == 0) then
            status = icsnp_numerical_error
            exit
         end if
         inverse_sum = sum(1.0_dp / pack(distances, active))
         updated = 0.0_dp
         direction_sum = 0.0_dp
         do i = 1, n
            if (.not. active(i)) cycle
            updated = updated + x(i, :) / distances(i)
            direction_sum = direction_sum + (x(i, :) - current) / distances(i)
         end do
         updated = updated / inverse_sum
         if (coincident > 0) then
            residual_norm = sqrt(dot_product(direction_sum, direction_sum))
            if (residual_norm > sqrt(tiny(1.0_dp))) then
               eta_ratio = real(coincident, dp) / residual_norm
               updated = max(0.0_dp, 1.0_dp - eta_ratio) * updated + &
                  min(1.0_dp, eta_ratio) * current
            else
               updated = current
            end if
         end if
         difference = sqrt(dot_product(updated - current, updated - current))
         current = updated
         if (difference <= tol) then
            status = icsnp_ok
            exit
         end if
      end do
      deallocate(center)
      allocate(center(p))
      center = current
      if (present(iterations)) iterations = min(iter, limit)
   end subroutine spatial_median

   subroutine tyler_shape(x, shape, status, iterations, location, init, maxiter, &
      tolerance, steps)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: shape(:,:)
      integer, intent(out) :: status
      integer, intent(out), optional :: iterations
      real(dp), intent(in), optional :: location(:), init(:,:), tolerance
      integer, intent(in), optional :: maxiter, steps
      real(dp), allocatable :: centered_all(:,:), centered(:,:), current(:,:), updated(:,:)
      real(dp), allocatable :: inverse(:,:), covariance(:,:)
      real(dp) :: center(size(x, 2)), q, tol, difference
      integer :: n, p, i, used, limit, stop_steps, iter, local_status

      n = size(x, 1)
      p = size(x, 2)
      allocate(shape(0, 0))
      status = icsnp_invalid_input
      if (present(iterations)) iterations = 0
      if (n < p + 1 .or. p < 2) return
      if (present(location)) then
         if (size(location) /= p) return
         center = location
      else
         center = sample_mean(x)
      end if
      allocate(centered_all(n, p))
      centered_all = x - spread(center, 1, n)
      used = 0
      do i = 1, n
         if (dot_product(centered_all(i, :), centered_all(i, :)) > tiny(1.0_dp)) used = used + 1
      end do
      if (used <= p) then
         status = icsnp_singular
         return
      end if
      allocate(centered(used, p))
      used = 0
      do i = 1, n
         if (dot_product(centered_all(i, :), centered_all(i, :)) <= tiny(1.0_dp)) cycle
         used = used + 1
         centered(used, :) = centered_all(i, :)
      end do

      if (present(init)) then
         if (size(init, 1) /= p .or. size(init, 2) /= p) return
         allocate(current(p, p))
         current = 0.5_dp * (init + transpose(init))
      else
         call covariance_matrix(centered, covariance, local_status, center=spread(0.0_dp, 1, p))
         if (local_status /= icsnp_ok) then
            status = local_status
            return
         end if
         allocate(current(p, p))
         current = covariance
      end if
      call normalize_det(current, status)
      if (status /= icsnp_ok) return

      tol = 1.0e-6_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      limit = 100
      if (present(maxiter)) limit = maxiter
      stop_steps = huge(1)
      if (present(steps)) stop_steps = max(0, steps)
      allocate(updated(p, p))
      status = icsnp_iteration_limit
      do iter = 1, limit
         call invert_matrix(current, inverse, local_status)
         if (local_status /= icsnp_ok) then
            status = local_status
            exit
         end if
         updated = 0.0_dp
         do i = 1, used
            q = dot_product(centered(i, :), matmul(inverse, centered(i, :)))
            if (q <= tiny(1.0_dp)) cycle
            updated = updated + spread(centered(i, :), 2, p) * &
               spread(centered(i, :), 1, p) / q
         end do
         updated = real(p, dp) * updated / real(used, dp)
         updated = 0.5_dp * (updated + transpose(updated))
         call normalize_det(updated, local_status)
         if (local_status /= icsnp_ok) then
            status = local_status
            exit
         end if
         difference = frobenius_norm(updated - current)
         current = updated
         if (difference <= tol .or. iter >= stop_steps) then
            status = icsnp_ok
            exit
         end if
      end do
      if (present(iterations)) iterations = min(iter, limit)
      if (status == icsnp_ok) then
         deallocate(shape)
         allocate(shape(p, p))
         shape = current
      end if
   end subroutine tyler_shape

   subroutine normalize_det(a, status)
      real(dp), intent(inout) :: a(:,:)
      integer, intent(out) :: status
      real(dp) :: det_value
      integer :: p
      p = size(a, 1)
      call determinant(a, det_value, status)
      if (status /= icsnp_ok .or. det_value <= 0.0_dp) then
         status = icsnp_singular
         return
      end if
      a = a / det_value**(1.0_dp / real(p, dp))
      status = icsnp_ok
   end subroutine normalize_det

   subroutine duembgen_shape(x, shape, status, iterations, init, maxiter, tolerance, steps)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: shape(:,:)
      integer, intent(out) :: status
      integer, intent(out), optional :: iterations
      real(dp), intent(in), optional :: init(:,:), tolerance
      integer, intent(in), optional :: maxiter, steps
      real(dp), allocatable :: differences(:,:)
      real(dp) :: zero(size(x, 2))
      zero = 0.0_dp
      call pair_diff(x, differences, status)
      if (status /= icsnp_ok) then
         allocate(shape(0, 0))
         return
      end if
      if (present(init)) then
         if (present(tolerance) .and. present(maxiter) .and. present(steps)) then
            call tyler_shape(differences, shape, status, iterations, zero, init, maxiter, tolerance, steps)
         else if (present(tolerance) .and. present(maxiter)) then
            call tyler_shape(differences, shape, status, iterations, zero, init, maxiter, tolerance)
         else
            call tyler_shape(differences, shape, status, iterations, location=zero, init=init)
         end if
      else
         if (present(tolerance) .and. present(maxiter) .and. present(steps)) then
            call tyler_shape(differences, shape, status, iterations, location=zero, &
               maxiter=maxiter, tolerance=tolerance, steps=steps)
         else if (present(tolerance) .and. present(maxiter)) then
            call tyler_shape(differences, shape, status, iterations, location=zero, &
               maxiter=maxiter, tolerance=tolerance)
         else
            call tyler_shape(differences, shape, status, iterations, location=zero)
         end if
      end if
   end subroutine duembgen_shape

   subroutine duembgen_shape_wt(x, weights, shape, status, iterations, init, maxiter, tolerance)
      real(dp), intent(in) :: x(:,:), weights(:)
      real(dp), allocatable, intent(out) :: shape(:,:)
      integer, intent(out) :: status
      integer, intent(out), optional :: iterations
      real(dp), intent(in), optional :: init(:,:), tolerance
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: differences(:,:), pair_weights_matrix(:,:), pair_weights(:)
      real(dp), allocatable :: current(:,:), updated(:,:), inverse(:,:), covariance(:,:)
      real(dp) :: q, tol, difference, weight_sum
      integer :: p, m, i, used, limit, iter, local_status
      logical, allocatable :: active(:)

      p = size(x, 2)
      allocate(shape(0, 0))
      status = icsnp_invalid_input
      if (present(iterations)) iterations = 0
      if (size(x, 1) /= size(weights) .or. p < 2) return
      if (any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) return
      call pair_diff(x, differences, status)
      if (status /= icsnp_ok) return
      call pair_prod(reshape(weights, [size(weights), 1]), pair_weights_matrix, status)
      if (status /= icsnp_ok) return
      allocate(pair_weights(size(pair_weights_matrix, 1)))
      pair_weights = pair_weights_matrix(:, 1)
      m = size(differences, 1)
      allocate(active(m))
      do i = 1, m
         active(i) = pair_weights(i) > 0.0_dp .and. &
            dot_product(differences(i, :), differences(i, :)) > tiny(1.0_dp)
      end do
      used = count(active)
      if (used <= p) then
         status = icsnp_singular
         return
      end if
      differences = pack_rows(differences, active)
      pair_weights = pack(pair_weights, active)
      weight_sum = sum(pair_weights)

      if (present(init)) then
         if (size(init, 1) /= p .or. size(init, 2) /= p) return
         allocate(current(p, p))
         current = init
      else
         allocate(covariance(p, p))
         covariance = 0.0_dp
         do i = 1, used
            covariance = covariance + pair_weights(i) * spread(differences(i, :), 2, p) * &
               spread(differences(i, :), 1, p)
         end do
         current = covariance / weight_sum
      end if
      call normalize_det(current, status)
      if (status /= icsnp_ok) return
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      limit = 100
      if (present(maxiter)) limit = maxiter
      allocate(updated(p, p))
      status = icsnp_iteration_limit
      do iter = 1, limit
         call invert_matrix(current, inverse, local_status)
         if (local_status /= icsnp_ok) then
            status = local_status
            exit
         end if
         updated = 0.0_dp
         do i = 1, used
            q = dot_product(differences(i, :), matmul(inverse, differences(i, :)))
            if (q <= tiny(1.0_dp)) cycle
            updated = updated + pair_weights(i) * spread(differences(i, :), 2, p) * &
               spread(differences(i, :), 1, p) / q
         end do
         updated = real(p, dp) * updated / weight_sum
         call normalize_det(updated, local_status)
         if (local_status /= icsnp_ok) then
            status = local_status
            exit
         end if
         difference = frobenius_norm(updated - current)
         current = updated
         if (difference <= tol) then
            status = icsnp_ok
            exit
         end if
      end do
      if (present(iterations)) iterations = min(iter, limit)
      if (status == icsnp_ok) then
         deallocate(shape)
         allocate(shape(p, p))
         shape = current
      end if
   end subroutine duembgen_shape_wt

   function pack_rows(x, mask) result(y)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in) :: mask(:)
      real(dp), allocatable :: y(:,:)
      integer :: i, j
      allocate(y(count(mask), size(x, 2)))
      j = 0
      do i = 1, size(mask)
         if (.not. mask(i)) cycle
         j = j + 1
         y(j, :) = x(i, :)
      end do
   end function pack_rows

   subroutine symm_huber(x, scatter, status, iterations, qg, init, maxiter, tolerance)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: scatter(:,:)
      integer, intent(out) :: status
      integer, intent(out), optional :: iterations
      real(dp), intent(in), optional :: qg, init(:,:), tolerance
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: differences(:,:), current(:,:), updated(:,:), inverse(:,:)
      real(dp) :: probability, c_square, sigma_square, distance, weight, tol, difference
      integer :: p, m, i, limit, iter, local_status

      p = size(x, 2)
      allocate(scatter(0, 0))
      status = icsnp_invalid_input
      if (present(iterations)) iterations = 0
      if (size(x, 1) < 2 .or. p < 2) return
      probability = 0.9_dp
      if (present(qg)) probability = qg
      if (probability <= 0.0_dp .or. probability >= 1.0_dp) return
      c_square = 2.0_dp * chi_square_quantile(probability, real(p, dp))
      sigma_square = 2.0_dp * (1.0_dp - chi_square_survival(0.5_dp * c_square, &
         real(p + 2, dp))) + (c_square / real(p, dp)) * (1.0_dp - probability)
      call pair_diff(x, differences, status)
      if (status /= icsnp_ok) return
      m = size(differences, 1)
      if (present(init)) then
         if (size(init, 1) /= p .or. size(init, 2) /= p) return
         allocate(current(p, p))
         current = init
      else
         call covariance_matrix(x, current, status)
         if (status /= icsnp_ok) return
      end if
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      limit = 100
      if (present(maxiter)) limit = maxiter
      allocate(updated(p, p))
      status = icsnp_iteration_limit
      do iter = 1, limit
         call invert_matrix(current, inverse, local_status)
         if (local_status /= icsnp_ok) then
            status = local_status
            exit
         end if
         updated = 0.0_dp
         do i = 1, m
            distance = dot_product(differences(i, :), matmul(inverse, differences(i, :)))
            if (distance <= c_square) then
               weight = 1.0_dp / sigma_square
            else
               weight = c_square / (distance * sigma_square)
            end if
            updated = updated + weight * spread(differences(i, :), 2, p) * &
               spread(differences(i, :), 1, p)
         end do
         updated = updated / real(m, dp)
         updated = 0.5_dp * (updated + transpose(updated))
         difference = frobenius_norm(updated - current)
         current = updated
         if (difference <= tol) then
            status = icsnp_ok
            exit
         end if
      end do
      if (present(iterations)) iterations = min(iter, limit)
      if (status == icsnp_ok) then
         deallocate(scatter)
         allocate(scatter(p, p))
         scatter = current
      end if
   end subroutine symm_huber

   subroutine symm_huber_wt(x, weights, scatter, status, iterations, qg, init, maxiter, tolerance)
      real(dp), intent(in) :: x(:,:), weights(:)
      real(dp), allocatable, intent(out) :: scatter(:,:)
      integer, intent(out) :: status
      integer, intent(out), optional :: iterations
      real(dp), intent(in), optional :: qg, init(:,:), tolerance
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: differences(:,:), pair_weights_matrix(:,:), pair_weights(:)
      real(dp), allocatable :: current(:,:), updated(:,:), inverse(:,:)
      real(dp) :: probability, c_square, sigma_square, distance, robust_weight
      real(dp) :: tol, difference, weight_sum
      integer :: p, m, i, limit, iter, local_status

      p = size(x, 2)
      allocate(scatter(0, 0))
      status = icsnp_invalid_input
      if (present(iterations)) iterations = 0
      if (size(x, 1) /= size(weights) .or. p < 2) return
      if (any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) return
      probability = 0.9_dp
      if (present(qg)) probability = qg
      if (probability <= 0.0_dp .or. probability >= 1.0_dp) return
      c_square = 2.0_dp * chi_square_quantile(probability, real(p, dp))
      sigma_square = 2.0_dp * (1.0_dp - chi_square_survival(0.5_dp * c_square, &
         real(p + 2, dp))) + (c_square / real(p, dp)) * (1.0_dp - probability)
      call pair_diff(x, differences, status)
      if (status /= icsnp_ok) return
      call pair_prod(reshape(weights, [size(weights), 1]), pair_weights_matrix, status)
      if (status /= icsnp_ok) return
      allocate(pair_weights(size(pair_weights_matrix, 1)))
      pair_weights = pair_weights_matrix(:, 1)
      weight_sum = sum(pair_weights)
      if (weight_sum <= 0.0_dp) return
      m = size(differences, 1)
      if (present(init)) then
         if (size(init, 1) /= p .or. size(init, 2) /= p) return
         allocate(current(p, p))
         current = init
      else
         allocate(current(p, p))
         current = 0.0_dp
         do i = 1, m
            current = current + pair_weights(i) * spread(differences(i, :), 2, p) * &
               spread(differences(i, :), 1, p)
         end do
         current = current / weight_sum
      end if
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      limit = 100
      if (present(maxiter)) limit = maxiter
      allocate(updated(p, p))
      status = icsnp_iteration_limit
      do iter = 1, limit
         call invert_matrix(current, inverse, local_status)
         if (local_status /= icsnp_ok) then
            status = local_status
            exit
         end if
         updated = 0.0_dp
         do i = 1, m
            distance = dot_product(differences(i, :), matmul(inverse, differences(i, :)))
            if (distance <= c_square) then
               robust_weight = 1.0_dp / sigma_square
            else
               robust_weight = c_square / (distance * sigma_square)
            end if
            updated = updated + pair_weights(i) * robust_weight * &
               spread(differences(i, :), 2, p) * spread(differences(i, :), 1, p)
         end do
         updated = updated / weight_sum
         updated = 0.5_dp * (updated + transpose(updated))
         difference = frobenius_norm(updated - current)
         current = updated
         if (difference <= tol) then
            status = icsnp_ok
            exit
         end if
      end do
      if (present(iterations)) iterations = min(iter, limit)
      if (status == icsnp_ok) then
         deallocate(scatter)
         allocate(scatter(p, p))
         scatter = current
      end if
   end subroutine symm_huber_wt

   subroutine HR_Mest(x, result, maxiter, eps_scale, eps_center)
      real(dp), intent(in) :: x(:,:)
      type(location_scatter_result), intent(out) :: result
      integer, intent(in), optional :: maxiter
      real(dp), intent(in), optional :: eps_scale, eps_center
      real(dp), allocatable :: theta1(:), theta2(:), shape1(:,:), shape2(:,:)
      real(dp), allocatable :: inv_sqrt(:,:), transformed(:,:), transformed_center(:)
      real(dp) :: tol_scale, tol_center, difference
      integer :: p, limit, iter, status, center_iterations

      p = size(x, 2)
      result%status = icsnp_invalid_input
      result%iterations = 0
      allocate(result%center(0), result%scatter(0, 0))
      if (size(x, 1) <= p .or. p < 2) return
      limit = 100
      if (present(maxiter)) limit = maxiter
      tol_scale = 1.0e-6_dp
      tol_center = 1.0e-6_dp
      if (present(eps_scale)) tol_scale = eps_scale
      if (present(eps_center)) tol_center = eps_center
      allocate(theta2(p))
      do iter = 1, p
         theta2(iter) = median_value(x(:, iter))
      end do
      result%status = icsnp_iteration_limit
      do iter = 1, limit
         call tyler_shape(x, shape1, status, location=theta2, maxiter=limit, tolerance=tol_scale)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         call matrix_inv_sqrt(shape1, inv_sqrt, status)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         transformed = matmul(x, inv_sqrt)
         call spatial_median(transformed, transformed_center, status, center_iterations, &
            maxiter=limit, tolerance=tol_center)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         theta1 = matmul(transformed_center, matrix_sqrt_value(shape1, status))
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         call tyler_shape(x, shape2, status, location=theta1, maxiter=limit, tolerance=tol_scale)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         call matrix_inv_sqrt(shape2, inv_sqrt, status)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         transformed = matmul(x, inv_sqrt)
         call spatial_median(transformed, transformed_center, status, center_iterations, &
            maxiter=limit, tolerance=tol_center)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         theta2 = matmul(transformed_center, matrix_sqrt_value(shape2, status))
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         difference = sqrt(dot_product(theta1 - theta2, theta1 - theta2))
         if (difference <= tol_scale) then
            result%status = icsnp_ok
            exit
         end if
      end do
      result%iterations = min(iter, limit)
      if (result%status == icsnp_ok) then
         deallocate(result%center, result%scatter)
         allocate(result%center(p), result%scatter(p, p))
         result%center = theta2
         result%scatter = shape2
      end if
   end subroutine HR_Mest

   function matrix_sqrt_value(a, status) result(root)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: root(:,:)
      call matrix_sqrt(a, root, status)
   end function matrix_sqrt_value

   subroutine spatial_sign(x, result, center, shape, estimate_center, estimate_shape, &
      maxiter, tolerance)
      real(dp), intent(in) :: x(:,:)
      type(spatial_sign_result), intent(out) :: result
      real(dp), intent(in), optional :: center(:), shape(:,:)
      logical, intent(in), optional :: estimate_center, estimate_shape
      integer, intent(in), optional :: maxiter
      real(dp), intent(in), optional :: tolerance
      type(location_scatter_result) :: hr
      real(dp), allocatable :: inv_sqrt(:,:), transformed(:,:), temp_center(:), root(:,:)
      real(dp) :: center_value(size(x, 2)), norm_value
      real(dp), allocatable :: shape_value(:,:)
      logical :: do_center, do_shape
      integer :: n, p, i, status, iters

      n = size(x, 1)
      p = size(x, 2)
      result%status = icsnp_invalid_input
      allocate(result%signs(0, 0), result%center(0), result%shape(0, 0))
      if (n < 1 .or. p < 1) return
      do_center = .true.
      do_shape = .true.
      if (present(estimate_center)) do_center = estimate_center
      if (present(estimate_shape)) do_shape = estimate_shape

      if (p == 1) then
         if (present(center)) then
            if (size(center) /= 1) return
            center_value = center
         else if (do_center) then
            center_value(1) = median_value(x(:, 1))
         else
            center_value = 0.0_dp
         end if
         deallocate(result%signs, result%center, result%shape)
         allocate(result%signs(n, 1), result%center(1), result%shape(1, 1))
         result%signs(:, 1) = sign(1.0_dp, x(:, 1) - center_value(1))
         where (abs(x(:, 1) - center_value(1)) <= 16.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, abs(center_value(1)))) result%signs(:, 1) = 0.0_dp
         result%center = center_value
         result%shape(1, 1) = 1.0_dp
         result%status = icsnp_ok
         return
      end if

      if (present(center) .and. present(shape)) then
         if (size(center) /= p .or. size(shape, 1) /= p .or. size(shape, 2) /= p) return
         center_value = center
         allocate(shape_value(p, p))
         shape_value = shape
      else if (present(center)) then
         if (size(center) /= p) return
         center_value = center
         if (do_shape) then
            call tyler_shape(x, shape_value, status, location=center_value, maxiter=optional_int(maxiter, 100), &
               tolerance=optional_real(tolerance, 1.0e-6_dp))
            if (status /= icsnp_ok) then
               result%status = status
               return
            end if
         else
            shape_value = identity_matrix(p)
         end if
      else if (present(shape)) then
         if (size(shape, 1) /= p .or. size(shape, 2) /= p) return
         allocate(shape_value(p, p))
         shape_value = shape
         if (do_center) then
            call matrix_inv_sqrt(shape_value, inv_sqrt, status)
            if (status /= icsnp_ok) then
               result%status = status
               return
            end if
            transformed = matmul(x, inv_sqrt)
            call spatial_median(transformed, temp_center, status, iters, &
               maxiter=optional_int(maxiter, 500), tolerance=optional_real(tolerance, 1.0e-6_dp))
            if (status /= icsnp_ok) then
               result%status = status
               return
            end if
            call matrix_sqrt(shape_value, root, status)
            center_value = matmul(temp_center, root)
         else
            center_value = 0.0_dp
         end if
      else if (do_center .and. do_shape) then
         call HR_Mest(x, hr, maxiter=optional_int(maxiter, 100), &
            eps_scale=optional_real(tolerance, 1.0e-6_dp), &
            eps_center=optional_real(tolerance, 1.0e-6_dp))
         if (hr%status /= icsnp_ok) then
            result%status = hr%status
            return
         end if
         center_value = hr%center
         shape_value = hr%scatter
      else if (do_shape) then
         center_value = 0.0_dp
         call tyler_shape(x, shape_value, status, location=center_value, &
            maxiter=optional_int(maxiter, 100), tolerance=optional_real(tolerance, 1.0e-6_dp))
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
      else if (do_center) then
         shape_value = identity_matrix(p)
         call spatial_median(x, temp_center, status, iters, maxiter=optional_int(maxiter, 500), &
            tolerance=optional_real(tolerance, 1.0e-6_dp))
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         center_value = temp_center
      else
         center_value = 0.0_dp
         shape_value = identity_matrix(p)
      end if

      call matrix_inv_sqrt(shape_value, inv_sqrt, status)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      transformed = matmul(x - spread(center_value, 1, n), inv_sqrt)
      deallocate(result%signs, result%center, result%shape)
      allocate(result%signs(n, p), result%center(p), result%shape(p, p))
      do i = 1, n
         norm_value = sqrt(dot_product(transformed(i, :), transformed(i, :)))
         if (norm_value > sqrt(tiny(1.0_dp))) then
            result%signs(i, :) = transformed(i, :) / norm_value
         else
            result%signs(i, :) = 0.0_dp
         end if
      end do
      result%center = center_value
      result%shape = shape_value
      result%status = icsnp_ok
   end subroutine spatial_sign

   pure integer function optional_int(value, default_value) result(output)
      integer, intent(in), optional :: value
      integer, intent(in) :: default_value
      output = default_value
      if (present(value)) output = value
   end function optional_int

   pure real(dp) function optional_real(value, default_value) result(output)
      real(dp), intent(in), optional :: value
      real(dp), intent(in) :: default_value
      output = default_value
      if (present(value)) output = value
   end function optional_real

   subroutine HP1_shape(x, shape, status, location, estimate_location, maxiter, tolerance)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: shape(:,:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: location(:)
      logical, intent(in), optional :: estimate_location
      integer, intent(in), optional :: maxiter
      real(dp), intent(in), optional :: tolerance
      type(spatial_sign_result) :: signs
      real(dp), allocatable :: root(:,:), distances(:), ranks(:), weighted_outer(:,:)
      real(dp) :: score
      integer :: n, p, i

      n = size(x, 1)
      p = size(x, 2)
      allocate(shape(0, 0))
      status = icsnp_invalid_input
      if (n < p + 1 .or. p < 2) return
      if (present(location)) then
         call spatial_sign(x, signs, center=location, estimate_center=.false., estimate_shape=.true., &
            maxiter=optional_int(maxiter, 100), tolerance=optional_real(tolerance, 1.0e-6_dp))
      else if (present(estimate_location)) then
         if (estimate_location) then
            call spatial_sign(x, signs, estimate_center=.true., estimate_shape=.true., &
               maxiter=optional_int(maxiter, 100), tolerance=optional_real(tolerance, 1.0e-6_dp))
         else
            call spatial_sign(x, signs, estimate_center=.false., estimate_shape=.true., &
               maxiter=optional_int(maxiter, 100), tolerance=optional_real(tolerance, 1.0e-6_dp))
         end if
      else
         call spatial_sign(x, signs, estimate_center=.true., estimate_shape=.true., &
            maxiter=optional_int(maxiter, 100), tolerance=optional_real(tolerance, 1.0e-6_dp))
      end if
      if (signs%status /= icsnp_ok) then
         status = signs%status
         return
      end if
      call matrix_sqrt(signs%shape, root, status)
      if (status /= icsnp_ok) return
      call mahalanobis_squared(x, signs%center, signs%shape, distances, status)
      if (status /= icsnp_ok) return
      distances = sqrt(distances)
      call rank_average(distances, ranks)
      allocate(weighted_outer(p, p))
      weighted_outer = 0.0_dp
      do i = 1, n
         score = chi_square_quantile(ranks(i) / real(n + 1, dp), real(p, dp))
         weighted_outer = weighted_outer + score * spread(signs%signs(i, :), 2, p) * &
            spread(signs%signs(i, :), 1, p)
      end do
      weighted_outer = weighted_outer / real(n, dp)
      deallocate(shape)
      allocate(shape(p, p))
      shape = matmul(root, matmul(weighted_outer, root))
      call normalize_det(shape, status)
   end subroutine HP1_shape

   real(dp) function hl_loc(x, status) result(location)
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: sums(:,:), all_values(:)
      integer :: local_status, n
      n = size(x)
      if (n < 1) then
         location = 0.0_dp
         if (present(status)) status = icsnp_invalid_input
         return
      end if
      if (n == 1) then
         location = x(1)
         if (present(status)) status = icsnp_ok
         return
      end if
      call pair_sum(reshape(x, [n, 1]), sums, local_status)
      allocate(all_values(n + size(sums, 1)))
      all_values(1:n) = x
      all_values(n + 1:) = 0.5_dp * sums(:, 1)
      location = median_value(all_values)
      if (present(status)) status = icsnp_ok
   end function hl_loc

   real(dp) function vdw_loc(x, status, int_diff, maxiter) result(location)
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      integer, intent(in), optional :: int_diff, maxiter
      real(dp), allocatable :: sums(:,:), candidates(:), criteria(:)
      integer :: n, m, lower, upper, middle, width, limit, iter, i, idx1, idx2
      real(dp) :: c1, c2, x1, x2, fraction

      n = size(x)
      location = 0.0_dp
      if (n < 1) then
         if (present(status)) status = icsnp_invalid_input
         return
      end if
      if (n == 1) then
         location = x(1)
         if (present(status)) status = icsnp_ok
         return
      end if
      call pair_sum(reshape(x, [n, 1]), sums)
      allocate(candidates(n + size(sums, 1)))
      candidates(1:n) = x
      candidates(n + 1:) = 0.5_dp * sums(:, 1)
      call sort_vector(candidates)
      m = size(candidates)
      width = 10
      if (present(int_diff)) width = max(2, int_diff)
      limit = 1000
      if (present(maxiter)) limit = maxiter
      lower = 1
      upper = m
      do iter = 1, limit
         if (upper - lower <= width) then
            allocate(criteria(upper - lower + 1))
            do i = lower, upper
               criteria(i - lower + 1) = vdw_criterion(candidates(i), x)
            end do
            idx1 = minloc(abs(criteria), dim=1)
            if (criteria(idx1) < 0.0_dp) then
               idx2 = min(size(criteria), idx1 + 1)
            else
               idx2 = max(1, idx1 - 1)
            end if
            c1 = criteria(idx1)
            c2 = criteria(idx2)
            x1 = candidates(lower + idx1 - 1)
            x2 = candidates(lower + idx2 - 1)
            if (abs(x1 - x2) <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x1), abs(x2)) .or. &
                abs(c1 - c2) <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(c1), abs(c2))) then
               location = x1
            else
               fraction = abs(c1) / abs(c2 - c1)
               location = x1 + sign(1.0_dp, c2) * fraction * abs(x2 - x1)
            end if
            if (present(status)) status = icsnp_ok
            return
         end if
         middle = (lower + upper) / 2
         if (vdw_criterion(candidates(middle), x) < 0.0_dp) then
            lower = middle
         else
            upper = middle
         end if
      end do
      location = candidates((lower + upper) / 2)
      if (present(status)) status = icsnp_iteration_limit
   end function vdw_loc

   real(dp) function vdw_criterion(mu, x) result(value)
      real(dp), intent(in) :: mu, x(:)
      real(dp), allocatable :: ranks(:)
      real(dp) :: z(size(x)), probabilities(size(x))
      integer :: n, i
      n = size(x)
      z = x - mu
      call rank_average(abs(z), ranks)
      probabilities = 0.5_dp + ranks / (2.0_dp * real(n + 1, dp))
      value = 0.0_dp
      do i = 1, n
         if (z(i) > 0.0_dp) value = value - normal_quantile(probabilities(i))
         if (z(i) < 0.0_dp) value = value + normal_quantile(probabilities(i))
      end do
   end function vdw_criterion

   subroutine sort_vector(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_vector

end module icsnp_estimators
