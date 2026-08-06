! sparseIndexTracking modern Fortran translation
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-only
!
! The algorithms in this module are derived from sparseIndexTracking 0.1.1
! by Konstantinos Benidis and Daniel P. Palomar, distributed under GPL-3.

module sparse_index_tracking_core
   use sparse_index_tracking_kinds, only : dp
   use sparse_index_tracking_linalg, only : all_finite, largest_eigenvalue_psd
   use sparse_index_tracking_projection, only : bisection
   implicit none
   private

   integer, parameter, public :: sit_success = 0
   integer, parameter, public :: sit_invalid_argument = 1
   integer, parameter, public :: sit_dimension_error = 2
   integer, parameter, public :: sit_infeasible_bounds = 3
   integer, parameter, public :: sit_degenerate_data = 4
   integer, parameter, public :: sit_iteration_limit = 5
   integer, parameter, public :: sit_numerical_error = 6

   integer, parameter, public :: measure_ete = 1
   integer, parameter, public :: measure_dr = 2
   integer, parameter, public :: measure_hete = 3
   integer, parameter, public :: measure_hdr = 4

   type, public :: sparse_index_fit
      real(dp), allocatable :: weights(:)
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0
      integer :: outer_iterations = 0
      integer :: cardinality = 0
      integer :: info = sit_success
      logical :: converged = .false.
      character(len=256) :: message = 'not fitted'
   end type sparse_index_fit

   public :: fit_sparse_index_tracking
   public :: sp_index_track
   public :: tracking_objective
   public :: parse_measure

   interface spIndexTrack
      module procedure sp_index_track
   end interface spIndexTrack
   public :: spIndexTrack

contains

   subroutine fit_sparse_index_tracking(x, index_returns, lambda, fit, upper_bound, measure, &
                                        huber_parameter, initial_weights, threshold, max_iterations, &
                                        source_compatible_hdr_objective)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: index_returns(:)
      real(dp), intent(in) :: lambda
      type(sparse_index_fit), intent(out) :: fit
      real(dp), intent(in), optional :: upper_bound
      character(len=*), intent(in), optional :: measure
      real(dp), intent(in), optional :: huber_parameter
      real(dp), intent(in), optional :: initial_weights(:)
      real(dp), intent(in), optional :: threshold
      integer, intent(in), optional :: max_iterations
      logical, intent(in), optional :: source_compatible_hdr_objective

      real(dp), allocatable :: a_matrix(:, :), b_matrix(:, :), b_vector(:)
      real(dp), allocatable :: current(:), r_step(:), u_step(:), w1(:), w2(:), candidate(:)
      real(dp) :: acceleration, c1, current_objective, eigenvalue, gamma
      real(dp) :: hub, p, previous_objective, relative_change, tolerance, upper, weight_threshold
      real(dp) :: norm_r, norm_u
      integer :: eigen_info, k, maxit, measure_code, n, m, outer, projection_info
      logical :: first_in_outer, source_hdr, stopped

      call clear_fit(fit)
      m = size(x, 1)
      n = size(x, 2)
      upper = 1.0_dp
      if (present(upper_bound)) upper = upper_bound
      hub = -1.0_dp
      if (present(huber_parameter)) hub = huber_parameter
      weight_threshold = 1.0e-9_dp
      if (present(threshold)) weight_threshold = threshold
      maxit = 1000
      if (present(max_iterations)) maxit = max_iterations
      source_hdr = .false.
      if (present(source_compatible_hdr_objective)) source_hdr = source_compatible_hdr_objective

      measure_code = measure_ete
      if (present(measure)) measure_code = parse_measure(measure)

      if (m < 1 .or. n < 2) then
         call fail_fit(fit, sit_dimension_error, 'x must have at least one row and two columns', n)
         return
      end if
      if (size(index_returns) /= m) then
         call fail_fit(fit, sit_dimension_error, 'index_returns length must equal the number of rows in x', n)
         return
      end if
      if (.not. all_finite(x) .or. .not. all_finite(index_returns)) then
         call fail_fit(fit, sit_invalid_argument, 'x and index_returns must contain only finite values', n)
         return
      end if
      if (lambda <= 0.0_dp .or. .not. all_finite([lambda])) then
         call fail_fit(fit, sit_invalid_argument, 'lambda must be finite and greater than zero', n)
         return
      end if
      if (upper <= 0.0_dp .or. .not. all_finite([upper])) then
         call fail_fit(fit, sit_invalid_argument, 'upper_bound must be finite and greater than zero', n)
         return
      end if
      if (real(n, dp) * upper < 1.0_dp - 100.0_dp * epsilon(1.0_dp)) then
         call fail_fit(fit, sit_infeasible_bounds, 'n * upper_bound must be at least one', n)
         return
      end if
      if (measure_code == 0) then
         call fail_fit(fit, sit_invalid_argument, 'measure must be ete, dr, hete, or hdr', n)
         return
      end if
      if ((measure_code == measure_hete .or. measure_code == measure_hdr) .and. hub <= 0.0_dp) then
         call fail_fit(fit, sit_invalid_argument, 'huber_parameter must be greater than zero for hete and hdr', n)
         return
      end if
      if (weight_threshold < 0.0_dp .or. maxit < 1) then
         call fail_fit(fit, sit_invalid_argument, 'threshold must be nonnegative and max_iterations positive', n)
         return
      end if

      allocate(current(n), w1(n), w2(n), r_step(n), u_step(n), candidate(n), &
               b_matrix(n, n), b_vector(n))
      b_matrix = 0.0_dp
      b_vector = 0.0_dp
      if (present(initial_weights)) then
         if (size(initial_weights) /= n) then
            call fail_fit(fit, sit_dimension_error, 'initial_weights length must equal the number of columns in x', n)
            return
         end if
         if (.not. all_finite(initial_weights) .or. any(initial_weights < 0.0_dp) .or. &
             any(initial_weights > upper) .or. abs(sum(initial_weights) - 1.0_dp) > 1.0e-8_dp) then
            call fail_fit(fit, sit_invalid_argument, &
                          'initial_weights must be finite, feasible, and sum to one', n)
            return
         end if
         current = initial_weights
      else
         current = 1.0_dp / real(n, dp)
      end if

      if (measure_code == measure_ete .or. measure_code == measure_dr) then
         allocate(a_matrix(n, n))
         a_matrix = matmul(transpose(x), x) / real(m, dp)
         call largest_eigenvalue_psd(a_matrix, eigenvalue, eigen_info)
         if (eigen_info == 1 .or. eigen_info == 2) then
            call fail_fit(fit, sit_numerical_error, 'failed to compute the largest eigenvalue', n)
            return
         end if
         if (eigenvalue <= 100.0_dp * epsilon(1.0_dp)) then
            call fail_fit(fit, sit_degenerate_data, 'x has zero or numerically negligible quadratic scale', n)
            return
         end if
         b_matrix = 2.0_dp / eigenvalue * a_matrix
         do k = 1, n
            b_matrix(k, k) = b_matrix(k, k) - 2.0_dp
         end do
         b_vector = -2.0_dp / real(m, dp) * matmul(transpose(x), index_returns)
      else
         eigenvalue = 0.0_dp
         eigen_info = 0
      end if

      gamma = 7.0_dp ** (1.0_dp / 10.0_dp)
      previous_objective = huge(1.0_dp)
      stopped = .false.

      outer_loop: do outer = 0, 10
         p = 10.0_dp ** (-(gamma ** real(outer, dp)))
         c1 = log(1.0_dp + upper / p)
         tolerance = min(p / 10.0_dp, 1.0e-3_dp)
         first_in_outer = .true.

         inner_loop: do
            if (fit%iterations >= maxit) then
               fit%info = sit_iteration_limit
               fit%message = 'maximum number of MM iterations reached'
               stopped = .true.
               exit inner_loop
            end if

            fit%iterations = fit%iterations + 1

            call mm_update(current, x, index_returns, lambda, p, c1, upper, measure_code, hub, &
                           b_matrix, b_vector, eigenvalue, w1, projection_info)
            if (projection_info /= 0) then
               call fail_fit(fit, sit_numerical_error, 'failed to compute the first MM update', n)
               return
            end if
            call mm_update(w1, x, index_returns, lambda, p, c1, upper, measure_code, hub, &
                           b_matrix, b_vector, eigenvalue, w2, projection_info)
            if (projection_info /= 0) then
               call fail_fit(fit, sit_numerical_error, 'failed to compute the second MM update', n)
               return
            end if

            r_step = w1 - current
            u_step = w2 - w1 - r_step
            norm_r = norm2(r_step)
            norm_u = norm2(u_step)
            if (norm_u <= 100.0_dp * tiny(1.0_dp)) then
               if (norm_r <= 100.0_dp * tiny(1.0_dp)) then
                  acceleration = -1.0_dp
               else
                  acceleration = -300.0_dp
               end if
            else
               acceleration = max(min(-norm_r / norm_u, -1.0_dp), -300.0_dp)
            end if

            do k = 1, 200
               if (abs(acceleration + 1.0_dp) < 1.0e-6_dp) then
                  candidate = w2
               else
                  candidate = current - 2.0_dp * acceleration * r_step + &
                              acceleration * acceleration * u_step
                  call project_accelerated(candidate, upper, projection_info)
                  if (projection_info /= 0) then
                     call fail_fit(fit, sit_numerical_error, 'accelerated projection failed', n)
                     return
                  end if
               end if

               current_objective = tracking_objective(x, index_returns, candidate, lambda, p, c1, &
                                                       measure_code, hub, source_hdr)
               if (.not. all_finite([current_objective])) then
                  call fail_fit(fit, sit_numerical_error, 'objective became nonfinite', n)
                  return
               end if

               if (.not. first_in_outer .and. &
                   current_objective * (1.0_dp - sign(1.0e-9_dp, current_objective)) >= &
                   previous_objective) then
                  acceleration = 0.5_dp * (acceleration - 1.0_dp)
               else
                  exit
               end if
            end do
            if (k > 200) then
               call fail_fit(fit, sit_numerical_error, 'acceleration backtracking failed to find descent', n)
               return
            end if

            current = candidate
            if (.not. first_in_outer) then
               relative_change = abs(current_objective - previous_objective) / &
                                 max(1.0_dp, abs(previous_objective))
               if (relative_change <= tolerance) exit inner_loop
            end if
            previous_objective = current_objective
            first_in_outer = .false.
         end do inner_loop

         fit%outer_iterations = outer + 1
         if (stopped) exit outer_loop
         previous_objective = current_objective
      end do outer_loop

      where (current < weight_threshold) current = 0.0_dp
      if (sum(current) <= 0.0_dp) then
         call fail_fit(fit, sit_numerical_error, 'thresholding removed every portfolio weight', n)
         return
      end if
      current = current / sum(current)

      if (allocated(fit%weights)) deallocate(fit%weights)
      allocate(fit%weights(n))
      fit%weights = current
      fit%cardinality = count(current > 0.0_dp)
      fit%objective = tracking_objective(x, index_returns, current, lambda, p, c1, &
                                         measure_code, hub, source_hdr)
      if (fit%info == sit_success) then
         fit%converged = .true.
         fit%message = 'ok'
      else
         fit%converged = .false.
      end if
   end subroutine fit_sparse_index_tracking


   subroutine sp_index_track(x, index_returns, lambda, weights, info, message, upper_bound, measure, &
                             huber_parameter, initial_weights, threshold, max_iterations, iterations, &
                             objective, converged, source_compatible_hdr_objective)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: index_returns(:)
      real(dp), intent(in) :: lambda
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message
      real(dp), intent(in), optional :: upper_bound
      character(len=*), intent(in), optional :: measure
      real(dp), intent(in), optional :: huber_parameter
      real(dp), intent(in), optional :: initial_weights(:)
      real(dp), intent(in), optional :: threshold
      integer, intent(in), optional :: max_iterations
      integer, intent(out), optional :: iterations
      real(dp), intent(out), optional :: objective
      logical, intent(out), optional :: converged
      logical, intent(in), optional :: source_compatible_hdr_objective

      type(sparse_index_fit) :: fit

      call fit_sparse_index_tracking(x, index_returns, lambda, fit, upper_bound, measure, &
                                     huber_parameter, initial_weights, threshold, max_iterations, &
                                     source_compatible_hdr_objective)
      allocate(weights(size(fit%weights)))
      weights = fit%weights
      if (present(info)) info = fit%info
      if (present(message)) message = fit%message
      if (present(iterations)) iterations = fit%iterations
      if (present(objective)) objective = fit%objective
      if (present(converged)) converged = fit%converged
   end subroutine sp_index_track


   function tracking_objective(x, index_returns, weights, lambda, p, c1, measure_code, &
                               huber_parameter, source_compatible_hdr) result(value)
      real(dp), intent(in) :: x(:, :), index_returns(:), weights(:)
      real(dp), intent(in) :: lambda, p, c1, huber_parameter
      integer, intent(in) :: measure_code
      logical, intent(in), optional :: source_compatible_hdr
      real(dp) :: value

      real(dp), allocatable :: residual(:), loss(:)
      logical :: source_hdr

      source_hdr = .false.
      if (present(source_compatible_hdr)) source_hdr = source_compatible_hdr
      allocate(residual(size(index_returns)), loss(size(index_returns)))
      residual = index_returns - matmul(x, weights)
      loss = 0.0_dp

      select case (measure_code)
      case (measure_ete)
         loss = residual * residual
      case (measure_dr)
         loss = max(residual, 0.0_dp) ** 2
      case (measure_hete)
         where (abs(residual) <= huber_parameter)
            loss = residual * residual
         elsewhere
            loss = huber_parameter * (2.0_dp * abs(residual) - huber_parameter)
         end where
      case (measure_hdr)
         if (source_hdr) then
            if (residual(1) > 0.0_dp .and. residual(1) <= huber_parameter) loss = residual * residual
            where (residual > huber_parameter)
               loss = huber_parameter * (2.0_dp * residual - huber_parameter)
            end where
         else
            where (residual > 0.0_dp .and. residual <= huber_parameter)
               loss = residual * residual
            elsewhere (residual > huber_parameter)
               loss = huber_parameter * (2.0_dp * residual - huber_parameter)
            elsewhere
               loss = 0.0_dp
            end where
         end if
      case default
         value = huge(1.0_dp)
         return
      end select

      value = sum(loss) / lambda + real(size(index_returns), dp) / c1 * &
              sum(log(1.0_dp + weights / p))
   end function tracking_objective


   integer function parse_measure(measure) result(code)
      character(len=*), intent(in) :: measure
      character(len=:), allocatable :: lowered

      lowered = lowercase(trim(adjustl(measure)))
      select case (lowered)
      case ('ete')
         code = measure_ete
      case ('dr')
         code = measure_dr
      case ('hete')
         code = measure_hete
      case ('hdr')
         code = measure_hdr
      case default
         code = 0
      end select
   end function parse_measure


   subroutine mm_update(weights, x, index_returns, lambda, p, c1, upper_bound, measure_code, hub, &
                        b_matrix, b_vector, fixed_eigenvalue, updated, info)
      real(dp), intent(in) :: weights(:), x(:, :), index_returns(:)
      real(dp), intent(in) :: lambda, p, c1, upper_bound, hub
      integer, intent(in) :: measure_code
      real(dp), intent(in) :: b_matrix(:, :), b_vector(:), fixed_eigenvalue
      real(dp), intent(out) :: updated(:)
      integer, intent(out) :: info

      real(dp), allocatable :: alpha(:), c(:), d(:), h(:), q(:), residual(:)
      real(dp), allocatable :: q_matrix(:, :), projected(:), weighted_response(:)
      real(dp) :: eigenvalue
      integer :: eigen_info, i, m, n

      m = size(x, 1)
      n = size(x, 2)
      allocate(c(n), d(n), residual(m))
      d = lambda / ((p + abs(weights)) * c1)
      residual = index_returns - matmul(x, weights)

      select case (measure_code)
      case (measure_ete)
         c = matmul(b_matrix, weights) + (b_vector + d) / fixed_eigenvalue
      case (measure_dr)
         allocate(h(m))
         h = min(residual, 0.0_dp)
         c = matmul(b_matrix, weights) + &
             (b_vector + d + 2.0_dp / real(m, dp) * matmul(transpose(x), h)) / fixed_eigenvalue
      case (measure_hete, measure_hdr)
         allocate(alpha(m), q_matrix(n, n), weighted_response(m))
         alpha = 1.0_dp
         if (measure_code == measure_hete) then
            where (abs(residual) > hub) alpha = hub / abs(residual)
            weighted_response = alpha * index_returns
         else
            where (residual > hub) alpha = hub / residual
            where (residual < 0.0_dp) alpha = hub / (hub - 2.0_dp * residual)
            allocate(q(m))
            q = -max(-residual, 0.0_dp)
            weighted_response = alpha * (q - index_returns)
         end if

         q_matrix = 0.0_dp
         do i = 1, m
            q_matrix = q_matrix + alpha(i) * outer_product(x(i, :), x(i, :))
         end do
         q_matrix = q_matrix / real(m, dp)
         call largest_eigenvalue_psd(q_matrix, eigenvalue, eigen_info)
         if (eigen_info == 1 .or. eigen_info == 2 .or. &
             eigenvalue <= 100.0_dp * epsilon(1.0_dp)) then
            info = 4
            updated = weights
            return
         end if

         if (measure_code == measure_hete) then
            c = (2.0_dp * matmul(q_matrix, weights) - 2.0_dp * eigenvalue * weights - &
                 2.0_dp / real(m, dp) * matmul(transpose(x), weighted_response) + d) / eigenvalue
         else
            c = (2.0_dp * matmul(q_matrix, weights) - 2.0_dp * eigenvalue * weights + &
                 2.0_dp / real(m, dp) * matmul(transpose(x), weighted_response) + d) / eigenvalue
         end if
      case default
         info = 5
         updated = weights
         return
      end select

      call bisection(c, upper_bound, projected, info)
      if (info == 0) updated = projected
   end subroutine mm_update


   subroutine project_accelerated(weights, upper_bound, info)
      real(dp), intent(inout) :: weights(:)
      real(dp), intent(in) :: upper_bound
      integer, intent(out) :: info

      real(dp), allocatable :: projected(:)

      call bisection(-2.0_dp * weights, upper_bound, projected, info)
      if (info == 0) weights = projected
   end subroutine project_accelerated


   pure function outer_product(x, y) result(matrix)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: matrix(size(x), size(y))
      integer :: j

      do j = 1, size(y)
         matrix(:, j) = x * y(j)
      end do
   end function outer_product


   pure function lowercase(text) result(lowered)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lowered
      integer :: code, i

      lowered = text
      do i = 1, len(text)
         code = iachar(lowered(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lowered(i:i) = achar(code + 32)
      end do
   end function lowercase


   subroutine clear_fit(fit)
      type(sparse_index_fit), intent(out) :: fit

      allocate(fit%weights(0))
      fit%objective = huge(1.0_dp)
      fit%iterations = 0
      fit%outer_iterations = 0
      fit%cardinality = 0
      fit%info = sit_success
      fit%converged = .false.
      fit%message = 'not fitted'
   end subroutine clear_fit


   subroutine fail_fit(fit, info, message, n)
      type(sparse_index_fit), intent(inout) :: fit
      integer, intent(in) :: info, n
      character(len=*), intent(in) :: message

      if (allocated(fit%weights)) deallocate(fit%weights)
      allocate(fit%weights(max(n, 0)))
      fit%weights = 0.0_dp
      fit%info = info
      fit%message = message
      fit%converged = .false.
   end subroutine fail_fit

end module sparse_index_tracking_core
