! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_crr
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_linalg, only : inverse_matrix, solve_system
   use cmprsk_status, only : cmprsk_success, cmprsk_invalid_argument, cmprsk_singular_matrix, &
                             cmprsk_no_failure_of_interest, cmprsk_no_convergence
   use cmprsk_utils, only : stable_order_real
   use cmprsk_censoring, only : censoring_survival_left
   use cmprsk_crr_kernels, only : crr_objective_score_info, crr_objective, crr_variance_kernel, &
                                  crr_score_residuals_kernel, crr_baseline_jumps_kernel
   implicit none
   private

   type, public :: crr_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: variance(:, :)
      real(dp), allocatable :: information(:, :)
      real(dp), allocatable :: score_residuals(:, :)
      real(dp), allocatable :: failure_time(:)
      real(dp), allocatable :: baseline_jump(:)
      real(dp), allocatable :: time_functions(:, :)
      real(dp) :: loglik = 0.0_dp
      real(dp) :: loglik_null = 0.0_dp
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: n_fixed = 0
      integer :: n_time_varying = 0
   end type crr_result

   public :: fit_crr
   public :: predict_crr

contains

   pure subroutine fit_crr(ftime, fstatus, cov1, result, status, failcode, cencode, cov2, time_functions, &
                      censor_group, initial, gtol, maxiter, compute_variance)
      real(dp), intent(in) :: ftime(:) !! Failure/censoring times in arbitrary row order.
      integer, intent(in) :: fstatus(:) !! Cause code for each row; `cencode` denotes censoring and `failcode` the modeled cause.
      real(dp), intent(in) :: cov1(:, :) !! Fixed-effect covariate matrix with one row per observation; zero columns are allowed.
      type(crr_result), intent(out) :: result !! Fine-Gray fit, covariance, residuals, baseline jumps, and failure-time metadata.
      integer, intent(out) :: status !! Success or a documented `cmprsk_status` failure code.
      integer, intent(in), optional :: failcode !! Cause code to model; defaults to `1`.
      integer, intent(in), optional :: cencode !! Censoring code; defaults to `0`.
      real(dp), intent(in), optional :: cov2(:, :) !! Covariates multiplied by supplied time functions; one row per observation.
      real(dp), intent(in), optional :: time_functions(:, :) !! Time functions at distinct modeled-cause failure times.
      integer, intent(in), optional :: censor_group(:) !! Censoring groups; arbitrary positive integer labels are accepted.
      real(dp), intent(in), optional :: initial(:) !! Starting coefficient vector; defaults to all zeros.
      real(dp), intent(in), optional :: gtol !! Scaled score convergence tolerance; defaults to `1e-6`.
      integer, intent(in), optional :: maxiter !! Maximum Newton iterations; defaults to `10` and may be zero.
      logical, intent(in), optional :: compute_variance !! Calculate sandwich covariance and score residuals; defaults to true.

      integer :: cencode_value
      integer :: failcode_value
      integer :: i
      integer :: info
      integer :: iter
      integer :: maxiter_value
      integer :: n
      integer :: ncg
      integer :: ndf
      integer :: n_fixed
      integer :: n_timed
      integer :: np
      integer :: censor_status
      integer, allocatable :: event_code(:)
      integer, allocatable :: group_raw(:)
      integer, allocatable :: group_sorted(:)
      integer, allocatable :: group_work(:)
      integer, allocatable :: order(:)
      real(dp) :: convergence_scale
      real(dp) :: gtol_value
      real(dp) :: new_objective
      real(dp) :: objective
      real(dp) :: old_objective
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: candidate(:)
      real(dp), allocatable :: cov1_sorted(:, :)
      real(dp), allocatable :: cov2_sorted(:, :)
      real(dp), allocatable :: covariance_inverse(:, :)
      real(dp), allocatable :: ftime_sorted(:)
      real(dp), allocatable :: information(:, :)
      real(dp), allocatable :: meat(:, :)
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: step(:)
      real(dp), allocatable :: survival(:, :)
      real(dp), allocatable :: tf(:, :)
      real(dp), allocatable :: residual_work(:, :)
      logical :: variance_value
      logical, allocatable :: censored(:)

      n = size(ftime)
      n_fixed = size(cov1, 2)
      n_timed = 0
      if (present(cov2)) n_timed = size(cov2, 2)
      np = n_fixed + n_timed
      failcode_value = 1
      if (present(failcode)) failcode_value = failcode
      cencode_value = 0
      if (present(cencode)) cencode_value = cencode
      gtol_value = 1.0e-6_dp
      if (present(gtol)) gtol_value = gtol
      maxiter_value = 10
      if (present(maxiter)) maxiter_value = maxiter
      variance_value = .true.
      if (present(compute_variance)) variance_value = compute_variance

      if (n < 1 .or. size(fstatus) /= n .or. size(cov1, 1) /= n .or. np < 1) then
         status = cmprsk_invalid_argument
         return
      end if
      if (any(ieee_is_nan(ftime)) .or. any(ieee_is_nan(cov1)) .or. any(ftime < 0.0_dp) .or. &
          gtol_value <= 0.0_dp .or. maxiter_value < 0) then
         status = cmprsk_invalid_argument
         return
      end if
      if (present(cov2)) then
         if (size(cov2, 1) /= n .or. .not. present(time_functions) .or. any(ieee_is_nan(cov2))) then
            status = cmprsk_invalid_argument
            return
         end if
      else if (present(time_functions)) then
         status = cmprsk_invalid_argument
         return
      end if
      if (present(censor_group)) then
         if (size(censor_group) /= n .or. any(censor_group < 1)) then
            status = cmprsk_invalid_argument
            return
         end if
      end if
      if (present(initial)) then
         if (size(initial) /= np .or. any(ieee_is_nan(initial))) then
            status = cmprsk_invalid_argument
            return
         end if
      end if

      allocate(order(n), ftime_sorted(n), event_code(n), censored(n))
      allocate(cov1_sorted(n, n_fixed), cov2_sorted(n, n_timed))
      allocate(group_raw(n), group_sorted(n), group_work(n))
      call stable_order_real(ftime, order)
      do i = 1, n
         ftime_sorted(i) = ftime(order(i))
         cov1_sorted(i, :) = cov1(order(i), :)
         if (n_timed > 0) cov2_sorted(i, :) = cov2(order(i), :)
         if (present(censor_group)) then
            group_raw(i) = censor_group(order(i))
         else
            group_raw(i) = 1
         end if
         if (fstatus(order(i)) == cencode_value) then
            event_code(i) = 0
         else if (fstatus(order(i)) == failcode_value) then
            event_code(i) = 1
         else
            event_code(i) = 2
         end if
      end do
      censored = event_code == 0
      call compress_groups(group_raw, group_sorted, ncg)

      ndf = count_distinct_failures(ftime_sorted, event_code)
      if (ndf == 0) then
         status = cmprsk_no_failure_of_interest
         return
      end if
      allocate(result%failure_time(ndf), tf(ndf, n_timed))
      call collect_failure_times(ftime_sorted, event_code, result%failure_time)
      if (n_timed > 0) then
         if (size(time_functions, 1) /= ndf .or. size(time_functions, 2) /= n_timed .or. &
             any(ieee_is_nan(time_functions))) then
            status = cmprsk_invalid_argument
            return
         end if
         tf = time_functions
      end if

      allocate(survival(ncg, n))
      call censoring_survival_left(ftime_sorted, censored, group_sorted, survival, censor_status)
      if (censor_status /= cmprsk_success) then
         status = censor_status
         return
      end if

      allocate(beta(np), candidate(np), score(np), step(np), information(np, np))
      beta = 0.0_dp
      if (present(initial)) beta = initial
      call crr_objective(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                         group_sorted, beta*0.0_dp, old_objective)
      result%loglik_null = -old_objective

      call crr_objective_score_info(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                                    group_sorted, beta, objective, score, information)
      result%converged = .false.
      convergence_scale = max(abs(objective), 1.0_dp)*gtol_value
      if (maxval(abs(score)*max(abs(beta), 1.0_dp)) < convergence_scale) result%converged = .true.

      iter = 0
      do while (.not. result%converged .and. iter < maxiter_value)
         iter = iter + 1
         call solve_system(information, score, step, info)
         if (info /= 0) then
            status = cmprsk_singular_matrix
            return
         end if
         step = -step
         old_objective = objective
         candidate = beta + step
         call crr_objective(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                            group_sorted, candidate, new_objective)
         do i = 1, 20
            if (new_objective <= old_objective + 1.0e-4_dp*dot_product(step, score)) exit
            step = 0.5_dp*step
            candidate = beta + step
            call crr_objective(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                               group_sorted, candidate, new_objective)
         end do
         beta = candidate
         call crr_objective_score_info(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                                       group_sorted, beta, objective, score, information)
         convergence_scale = max(abs(objective), 1.0_dp)*gtol_value
         if (maxval(abs(score)*max(abs(beta), 1.0_dp)) < convergence_scale) result%converged = .true.
      end do

      result%iterations = iter
      result%n_fixed = n_fixed
      result%n_time_varying = n_timed
      result%loglik = -objective
      allocate(result%coefficients(np), result%information(np, np), result%time_functions(ndf, n_timed))
      result%coefficients = beta
      result%information = information
      if (n_timed > 0) result%time_functions = tf
      allocate(result%baseline_jump(ndf))
      call crr_baseline_jumps_kernel(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                                     group_sorted, beta, result%baseline_jump)

      if (variance_value) then
         allocate(meat(np, np))
         call crr_variance_kernel(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                                  group_sorted, beta, information, meat)
         call inverse_matrix(information, covariance_inverse, info)
         if (info /= 0) then
            status = cmprsk_singular_matrix
            return
         end if
         allocate(result%variance(np, np), residual_work(np, ndf), result%score_residuals(ndf, np))
         result%variance = matmul(matmul(covariance_inverse, meat), transpose(covariance_inverse))
         result%information = information
         call crr_score_residuals_kernel(ftime_sorted, event_code, cov1_sorted, cov2_sorted, tf, survival, &
                                         group_sorted, beta, residual_work)
         result%score_residuals = transpose(residual_work)
      else
         allocate(result%variance(0, 0), result%score_residuals(0, 0))
      end if

      if (result%converged) then
         status = cmprsk_success
      else
         status = cmprsk_no_convergence
      end if
   end subroutine fit_crr

   pure subroutine predict_crr(result, cov1, cif, status, cov2)
      type(crr_result), intent(in) :: result !! Fitted Fine-Gray model returned by `fit_crr`.
      real(dp), intent(in) :: cov1(:, :) !! Fixed-effect profiles, with one profile per row.
      real(dp), allocatable, intent(out) :: cif(:, :) !! Predicted cumulative incidence by failure time and profile.
      integer, intent(out) :: status !! `cmprsk_success` or `cmprsk_invalid_argument`.
      real(dp), intent(in), optional :: cov2(:, :) !! Profile covariates corresponding to fitted time-varying columns.

      integer :: i
      integer :: j
      integer :: n_profiles
      real(dp) :: cumulative_hazard
      real(dp) :: linear_predictor

      n_profiles = size(cov1, 1)
      if (size(cov1, 2) /= result%n_fixed .or. any(ieee_is_nan(cov1))) then
         status = cmprsk_invalid_argument
         return
      end if
      if (result%n_time_varying > 0) then
         if (.not. present(cov2)) then
            status = cmprsk_invalid_argument
            return
         end if
         if (size(cov2, 1) /= n_profiles .or. size(cov2, 2) /= result%n_time_varying .or. &
             any(ieee_is_nan(cov2))) then
            status = cmprsk_invalid_argument
            return
         end if
      else if (present(cov2)) then
         if (size(cov2, 2) /= 0) then
            status = cmprsk_invalid_argument
            return
         end if
      end if

      allocate(cif(size(result%failure_time), n_profiles))
      do j = 1, n_profiles
         cumulative_hazard = 0.0_dp
         do i = 1, size(result%failure_time)
            linear_predictor = 0.0_dp
            if (result%n_fixed > 0) then
               linear_predictor = dot_product(cov1(j, :), result%coefficients(1:result%n_fixed))
            end if
            if (result%n_time_varying > 0) then
               linear_predictor = linear_predictor + dot_product(cov2(j, :)*result%time_functions(i, :), &
                    result%coefficients(result%n_fixed + 1:))
            end if
            cumulative_hazard = cumulative_hazard + exp(linear_predictor)*result%baseline_jump(i)
            cif(i, j) = 1.0_dp - exp(-cumulative_hazard)
         end do
      end do
      status = cmprsk_success
   end subroutine predict_crr

   pure subroutine compress_groups(raw, compressed, number_groups)
      integer, intent(in) :: raw(:) !! Positive arbitrary group labels in sorted observation order.
      integer, intent(out) :: compressed(:) !! Consecutive one-based group codes preserving equality of raw labels.
      integer, intent(out) :: number_groups !! Number of distinct labels found in `raw`.

      integer :: i
      integer :: j
      integer, allocatable :: labels(:)
      logical :: found

      allocate(labels(size(raw)))
      number_groups = 0
      do i = 1, size(raw)
         found = .false.
         do j = 1, number_groups
            if (raw(i) == labels(j)) then
               compressed(i) = j
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            number_groups = number_groups + 1
            labels(number_groups) = raw(i)
            compressed(i) = number_groups
         end if
      end do
   end subroutine compress_groups

   pure integer function count_distinct_failures(time, event_code) result(count)
      real(dp), intent(in) :: time(:) !! Sorted follow-up times.
      integer, intent(in) :: event_code(:) !! Event codes with one denoting the modeled cause.

      integer :: i
      real(dp) :: last_time
      logical :: have_last

      count = 0
      have_last = .false.
      last_time = 0.0_dp
      do i = 1, size(time)
         if (event_code(i) /= 1) cycle
         if (.not. have_last .or. time(i) /= last_time) then
            count = count + 1
            last_time = time(i)
            have_last = .true.
         end if
      end do
   end function count_distinct_failures

   pure subroutine collect_failure_times(time, event_code, failure_time)
      real(dp), intent(in) :: time(:) !! Sorted follow-up times.
      integer, intent(in) :: event_code(:) !! Event codes with one denoting the modeled cause.
      real(dp), intent(out) :: failure_time(:) !! Ascending distinct modeled-cause failure times.

      integer :: i
      integer :: k

      k = 0
      do i = 1, size(time)
         if (event_code(i) /= 1) cycle
         if (k == 0) then
            k = 1
            failure_time(k) = time(i)
         else if (time(i) /= failure_time(k)) then
            k = k + 1
            failure_time(k) = time(i)
         end if
      end do
   end subroutine collect_failure_times

end module cmprsk_crr
