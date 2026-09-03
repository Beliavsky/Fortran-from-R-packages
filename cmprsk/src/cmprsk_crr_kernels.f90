! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_crr_kernels
   use r_kinds, only : dp
   implicit none
   private

   public :: crr_objective_score_info
   public :: crr_objective
   public :: crr_variance_kernel
   public :: crr_score_residuals_kernel
   public :: crr_baseline_jumps_kernel

contains

   pure subroutine covariate_at(row_index, failure_index, x, x2, time_functions, beta, linear_predictor, covariate)
      integer, intent(in) :: row_index !! One-based subject row whose covariates are required.
      integer, intent(in) :: failure_index !! One-based distinct target-failure time selecting the time-function row.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect covariate matrix, with subjects in rows.
      real(dp), intent(in) :: x2(:, :) !! Covariate matrix whose columns are multiplied by time functions.
      real(dp), intent(in) :: time_functions(:, :) !! Time-function values by distinct target-failure time and time-varying column.
      real(dp), intent(in) :: beta(:) !! Regression coefficients, fixed effects followed by time-varying effects.
      real(dp), intent(out) :: linear_predictor !! Linear predictor for the requested subject and failure time.
      real(dp), intent(out) :: covariate(:) !! Combined fixed and time-varying covariate vector for the requested subject/time.

      integer :: i
      integer :: n_fixed
      integer :: n_timed

      n_fixed = size(x, 2)
      n_timed = size(x2, 2)
      linear_predictor = 0.0_dp
      if (n_fixed > 0) then
         do i = 1, n_fixed
            covariate(i) = x(row_index, i)
            linear_predictor = linear_predictor + covariate(i)*beta(i)
         end do
      end if
      if (n_timed > 0) then
         do i = 1, n_timed
            covariate(n_fixed + i) = x2(row_index, i)*time_functions(failure_index, i)
            linear_predictor = linear_predictor + covariate(n_fixed + i)*beta(n_fixed + i)
         end do
      end if
   end subroutine covariate_at

   pure subroutine crr_objective_score_info(time, event_code, x, x2, time_functions, censoring_survival, &
                                             censor_group, beta, objective, score, information)
      real(dp), intent(in) :: time(:) !! Follow-up times sorted in nondecreasing order.
      integer, intent(in) :: event_code(:) !! Zero for censoring, one for target failure, and two for competing failure.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect covariates, with one row per subject.
      real(dp), intent(in) :: x2(:, :) !! Time-varying-effect covariates, with one row per subject.
      real(dp), intent(in) :: time_functions(:, :) !! Time-function values at distinct target-failure times.
      real(dp), intent(in) :: censoring_survival(:, :) !! Left-continuous censoring KM values by group and sorted subject time.
      integer, intent(in) :: censor_group(:) !! Consecutive one-based censoring-group code for each subject.
      real(dp), intent(in) :: beta(:) !! Current Fine-Gray regression coefficient vector.
      real(dp), intent(out) :: objective !! Negative Fine-Gray log pseudo-likelihood at `beta`.
      real(dp), intent(out) :: score(:) !! Gradient of the negative log pseudo-likelihood.
      real(dp), intent(out) :: information(:, :) !! Hessian/information matrix of the negative log pseudo-likelihood.

      integer :: i
      integer :: iuc
      integer :: itmp
      integer :: j
      integer :: k
      integer :: ldf
      integer :: n
      integer :: np
      real(dp) :: cft
      real(dp) :: linear_predictor
      real(dp) :: risk_sum
      real(dp) :: risk_sum_old
      real(dp) :: target_count
      real(dp) :: weight
      real(dp), allocatable :: covariate(:)
      real(dp), allocatable :: weighted_covariate(:)
      real(dp), allocatable :: weighted_cross(:, :)

      n = size(time)
      np = size(beta)
      allocate(covariate(np), weighted_covariate(np), weighted_cross(np, np))
      objective = 0.0_dp
      score = 0.0_dp
      information = 0.0_dp
      iuc = n
      ldf = size(time_functions, 1) + 1

      do
         itmp = iuc
         do i = iuc, 1, -1
            itmp = i
            if (event_code(i) == 1) then
               cft = time(i)
               exit
            end if
         end do
         if (event_code(itmp) /= 1) return
         iuc = itmp
         ldf = ldf - 1
         target_count = 0.0_dp
         do i = iuc, 1, -1
            if (time(i) < cft) exit
            itmp = i
            if (event_code(i) == 1) then
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               target_count = target_count + 1.0_dp
               objective = objective - linear_predictor
               score = score - covariate
            end if
         end do
         iuc = itmp

         risk_sum = 0.0_dp
         risk_sum_old = 0.0_dp
         weighted_covariate = 0.0_dp
         weighted_cross = 0.0_dp
         do i = 1, n
            if (time(i) < cft) then
               if (event_code(i) <= 1) cycle
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)*censoring_survival(censor_group(i), iuc) / &
                        censoring_survival(censor_group(i), i)
            else
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)
            end if
            risk_sum = risk_sum + weight
            do j = 1, np
               weighted_covariate(j) = weighted_covariate(j) + weight*covariate(j)
               covariate(j) = covariate(j) - weighted_covariate(j)/risk_sum
            end do
            if (risk_sum_old > 0.0_dp) then
               weight = risk_sum*weight/risk_sum_old
               do k = 1, np
                  do j = k, np
                     weighted_cross(k, j) = weighted_cross(k, j) + weight*covariate(k)*covariate(j)
                  end do
               end do
            end if
            risk_sum_old = risk_sum
         end do
         objective = objective + target_count*log(risk_sum)
         weight = target_count/risk_sum
         do i = 1, np
            score(i) = score(i) + weight*weighted_covariate(i)
            do j = i, np
               information(i, j) = information(i, j) + weight*weighted_cross(i, j)
               information(j, i) = information(i, j)
            end do
         end do
         iuc = iuc - 1
         if (iuc <= 0) return
      end do
   end subroutine crr_objective_score_info

   pure subroutine crr_objective(time, event_code, x, x2, time_functions, censoring_survival, &
                                 censor_group, beta, objective)
      real(dp), intent(in) :: time(:) !! Follow-up times sorted in nondecreasing order.
      integer, intent(in) :: event_code(:) !! Zero for censoring, one for target failure, and two for competing failure.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect covariates, with one row per subject.
      real(dp), intent(in) :: x2(:, :) !! Time-varying-effect covariates, with one row per subject.
      real(dp), intent(in) :: time_functions(:, :) !! Time-function values at distinct target-failure times.
      real(dp), intent(in) :: censoring_survival(:, :) !! Left-continuous censoring KM values by group and sorted subject time.
      integer, intent(in) :: censor_group(:) !! Consecutive one-based censoring-group code for each subject.
      real(dp), intent(in) :: beta(:) !! Fine-Gray regression coefficient vector.
      real(dp), intent(out) :: objective !! Negative Fine-Gray log pseudo-likelihood at `beta`.

      integer :: i
      integer :: iuc
      integer :: itmp
      integer :: ldf
      integer :: n
      real(dp) :: cft
      real(dp) :: linear_predictor
      real(dp) :: risk_sum
      real(dp) :: target_count
      real(dp) :: weight
      real(dp), allocatable :: covariate(:)

      n = size(time)
      allocate(covariate(size(beta)))
      objective = 0.0_dp
      iuc = n
      ldf = size(time_functions, 1) + 1
      do
         itmp = iuc
         do i = iuc, 1, -1
            itmp = i
            if (event_code(i) == 1) then
               cft = time(i)
               exit
            end if
         end do
         if (event_code(itmp) /= 1) return
         iuc = itmp
         ldf = ldf - 1
         target_count = 0.0_dp
         do i = iuc, 1, -1
            if (time(i) < cft) exit
            itmp = i
            if (event_code(i) == 1) then
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               target_count = target_count + 1.0_dp
               objective = objective - linear_predictor
            end if
         end do
         iuc = itmp
         risk_sum = 0.0_dp
         do i = 1, n
            if (time(i) < cft) then
               if (event_code(i) <= 1) cycle
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)*censoring_survival(censor_group(i), iuc) / &
                        censoring_survival(censor_group(i), i)
            else
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)
            end if
            risk_sum = risk_sum + weight
         end do
         objective = objective + target_count*log(risk_sum)
         iuc = iuc - 1
         if (iuc <= 0) return
      end do
   end subroutine crr_objective

   pure subroutine crr_variance_kernel(time, event_code, x, x2, time_functions, censoring_survival, &
                                       censor_group, beta, information, meat)
      real(dp), intent(in) :: time(:) !! Follow-up times sorted in nondecreasing order.
      integer, intent(in) :: event_code(:) !! Zero for censoring, one for target failure, and two for competing failure.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect covariates, with one row per subject.
      real(dp), intent(in) :: x2(:, :) !! Time-varying-effect covariates, with one row per subject.
      real(dp), intent(in) :: time_functions(:, :) !! Time-function values at distinct target-failure times.
      real(dp), intent(in) :: censoring_survival(:, :) !! Left-continuous censoring KM values by group and sorted subject time.
      integer, intent(in) :: censor_group(:) !! Consecutive one-based censoring-group code for each subject.
      real(dp), intent(in) :: beta(:) !! Fitted Fine-Gray regression coefficient vector.
      real(dp), intent(out) :: information(:, :) !! Observed information component of the sandwich covariance.
      real(dp), intent(out) :: meat(:, :) !! Empirical score covariance component of the sandwich covariance.

      integer :: i
      integer :: iflag
      integer :: j
      integer :: j1
      integer :: j2
      integer :: k
      integer :: lc
      integer :: ldf
      integer :: ldf2
      integer :: n
      integer :: ncg
      integer :: np
      real(dp) :: cft
      real(dp) :: cft2
      real(dp) :: linear_predictor
      real(dp) :: previous_time
      real(dp) :: weight
      integer, allocatable :: risk_count(:)
      real(dp), allocatable :: covariate(:)
      real(dp), allocatable :: qu(:, :)
      real(dp), allocatable :: ss2(:, :)
      real(dp), allocatable :: ss3(:, :)
      real(dp), allocatable :: ss4(:)
      real(dp), allocatable :: st(:, :)
      real(dp), allocatable :: vt(:, :)
      real(dp), allocatable :: xb(:, :)

      n = size(time)
      np = size(beta)
      ncg = size(censoring_survival, 1)
      allocate(risk_count(ncg), covariate(np), qu(np, ncg), ss2(np, ncg), ss3(np, ncg), ss4(ncg))
      allocate(st(np, 2), vt(np, np), xb(n, 0:np))
      risk_count = 0
      do i = 1, n
         risk_count(censor_group(i)) = risk_count(censor_group(i)) + 1
      end do
      information = 0.0_dp
      meat = 0.0_dp
      ss2 = 0.0_dp
      xb = 0.0_dp

      ldf = 0
      cft = min(-1.0_dp, time(1)*(1.0_dp - 1.0e-5_dp))
      do i = 1, n
         if (event_code(i) /= 1) cycle
         if (time(i) > cft) then
            cft = time(i)
            ldf = ldf + 1
         end if
         do j = 1, n
            if (time(j) < time(i)) then
               if (event_code(j) <= 1) cycle
               call covariate_at(j, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)*censoring_survival(censor_group(j), i) / &
                        censoring_survival(censor_group(j), j)
            else
               call covariate_at(j, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)
            end if
            xb(i, 0) = xb(i, 0) + weight
            do k = 1, np
               xb(i, k) = xb(i, k) + weight*covariate(k)
            end do
         end do
      end do

      lc = 1
      ldf2 = 0
      cft2 = min(-1.0_dp, time(1)*(1.0_dp - 1.0e-5_dp))
      previous_time = time(1)
      do i = 1, n
         st(:, 1) = 0.0_dp
         ldf = 0
         cft = min(-1.0_dp, time(1)*(1.0_dp - 1.0e-5_dp))
         do j = 1, n
            if (event_code(j) /= 1) cycle
            if (time(j) > cft) then
               cft = time(j)
               ldf = ldf + 1
            end if
            if (time(j) <= time(i)) then
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)
            else if (time(i) < time(j) .and. event_code(i) > 1) then
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)*censoring_survival(censor_group(i), j) / &
                        censoring_survival(censor_group(i), i)
            else
               cycle
            end if
            do k = 1, np
               st(k, 1) = st(k, 1) - (covariate(k) - xb(j, k)/xb(j, 0))*weight/xb(j, 0)
            end do
         end do

         if (event_code(i) == 1) then
            if (time(i) > cft2) then
               cft2 = time(i)
               ldf2 = ldf2 + 1
            end if
            call covariate_at(i, ldf2, x, x2, time_functions, beta, linear_predictor, covariate)
            do k = 1, np
               st(k, 1) = st(k, 1) + covariate(k) - xb(i, k)/xb(i, 0)
            end do
            vt = 0.0_dp
            do j = 1, n
               if (time(j) < time(i)) then
                  if (event_code(j) <= 1) cycle
                  call covariate_at(j, ldf2, x, x2, time_functions, beta, linear_predictor, covariate)
                  weight = exp(linear_predictor)*censoring_survival(censor_group(j), i) / &
                           censoring_survival(censor_group(j), j)
               else
                  call covariate_at(j, ldf2, x, x2, time_functions, beta, linear_predictor, covariate)
                  weight = exp(linear_predictor)
               end if
               do k = 1, np
                  covariate(k) = covariate(k) - xb(i, k)/xb(i, 0)
               end do
               do k = 1, np
                  do j2 = k, np
                     vt(k, j2) = vt(k, j2) + weight*covariate(k)*covariate(j2)
                  end do
               end do
            end do
            do j1 = 1, np
               do j2 = j1, np
                  information(j1, j2) = information(j1, j2) + vt(j1, j2)/xb(i, 0)
               end do
            end do
         end if

         iflag = 1
         if (i > 1) then
            if (time(i) <= previous_time) iflag = 0
         end if
         previous_time = time(i)

         if (iflag == 1) then
            do j = i, n
               if (time(j) > time(i)) exit
               if (event_code(j) == 0) exit
            end do
            if (j <= n) then
               if (time(j) == time(i) .and. event_code(j) == 0) then
                  ldf = ldf2
               cft = cft2
               qu = 0.0_dp
               do j1 = lc, n
                  if (event_code(j1) /= 1) cycle
                  if (time(j1) > cft) then
                     cft = time(j1)
                     ldf = ldf + 1
                  end if
                  ss4 = 0.0_dp
                  ss3 = 0.0_dp
                  do j2 = 1, n
                     if (time(j2) >= time(i)) exit
                     if (event_code(j2) <= 1) cycle
                     call covariate_at(j2, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
                     weight = exp(linear_predictor)*censoring_survival(censor_group(j2), j1) / &
                              censoring_survival(censor_group(j2), j2)
                     ss4(censor_group(j2)) = ss4(censor_group(j2)) + weight
                     do k = 1, np
                        ss3(k, censor_group(j2)) = ss3(k, censor_group(j2)) + covariate(k)*weight
                     end do
                  end do
                  do k = 1, np
                     qu(k, censor_group(j1)) = qu(k, censor_group(j1)) + &
                        (ss3(k, censor_group(j1)) - xb(j1, k)*ss4(censor_group(j1))/xb(j1, 0))/xb(j1, 0)
                  end do
               end do
               do j = i, n
                  if (time(j) > time(i)) exit
                  if (event_code(j) == 0) then
                     do k = 1, np
                        ss2(k, censor_group(j)) = ss2(k, censor_group(j)) - &
                             qu(k, censor_group(j))/real(risk_count(censor_group(j))**2, dp)
                     end do
                  end if
               end do
               end if
            end if
         end if

         st(:, 2) = ss2(:, censor_group(i))
         if (event_code(i) == 0) then
            do k = 1, np
               st(k, 2) = st(k, 2) + qu(k, censor_group(i))/real(risk_count(censor_group(i)), dp)
            end do
         end if
         do j1 = np, 1, -1
            st(j1, 1) = st(j1, 1) + st(j1, 2)
            do j2 = j1, np
               meat(j1, j2) = meat(j1, j2) + st(j1, 1)*st(j2, 1)
            end do
         end do

         if (i < n) then
            if (time(i + 1) > time(i)) then
               do j = lc, i
                  risk_count(censor_group(j)) = risk_count(censor_group(j)) - 1
               end do
               lc = i + 1
            end if
         end if
      end do

      do j1 = 1, np
         do j2 = j1 + 1, np
            information(j2, j1) = information(j1, j2)
            meat(j2, j1) = meat(j1, j2)
         end do
      end do
   end subroutine crr_variance_kernel

   pure subroutine crr_score_residuals_kernel(time, event_code, x, x2, time_functions, censoring_survival, &
                                              censor_group, beta, residuals)
      real(dp), intent(in) :: time(:) !! Follow-up times sorted in nondecreasing order.
      integer, intent(in) :: event_code(:) !! Zero for censoring, one for target failure, and two for competing failure.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect covariates, with one row per subject.
      real(dp), intent(in) :: x2(:, :) !! Time-varying-effect covariates, with one row per subject.
      real(dp), intent(in) :: time_functions(:, :) !! Time-function values at distinct target-failure times.
      real(dp), intent(in) :: censoring_survival(:, :) !! Left-continuous censoring KM values by group and sorted time.
      integer, intent(in) :: censor_group(:) !! Consecutive one-based censoring-group code for each subject.
      real(dp), intent(in) :: beta(:) !! Fitted Fine-Gray regression coefficient vector.
      real(dp), intent(out) :: residuals(:, :) !! Score contribution by coefficient and distinct target-failure time.

      integer :: i
      integer :: iuc
      integer :: itmp
      integer :: ldf
      integer :: n
      integer :: np
      real(dp) :: cft
      real(dp) :: linear_predictor
      real(dp) :: risk_sum
      real(dp) :: target_count
      real(dp) :: weight
      real(dp), allocatable :: covariate(:)
      real(dp), allocatable :: weighted_covariate(:)

      n = size(time)
      np = size(beta)
      allocate(covariate(np), weighted_covariate(np))
      residuals = 0.0_dp
      iuc = n
      ldf = size(time_functions, 1) + 1
      do
         itmp = iuc
         do i = iuc, 1, -1
            itmp = i
            if (event_code(i) == 1) then
               cft = time(i)
               exit
            end if
         end do
         if (event_code(itmp) /= 1) return
         iuc = itmp
         ldf = ldf - 1
         target_count = 0.0_dp
         do i = iuc, 1, -1
            if (time(i) < cft) exit
            itmp = i
            if (event_code(i) == 1) then
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               target_count = target_count + 1.0_dp
               residuals(:, ldf) = residuals(:, ldf) + covariate
            end if
         end do
         iuc = itmp
         risk_sum = 0.0_dp
         weighted_covariate = 0.0_dp
         do i = 1, n
            if (time(i) < cft) then
               if (event_code(i) <= 1) cycle
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)*censoring_survival(censor_group(i), iuc) / &
                        censoring_survival(censor_group(i), i)
            else
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)
            end if
            risk_sum = risk_sum + weight
            weighted_covariate = weighted_covariate + weight*covariate
         end do
         weight = -target_count/risk_sum
         residuals(:, ldf) = residuals(:, ldf) + weight*weighted_covariate
         iuc = iuc - 1
         if (iuc <= 0) return
      end do
   end subroutine crr_score_residuals_kernel

   pure subroutine crr_baseline_jumps_kernel(time, event_code, x, x2, time_functions, censoring_survival, &
                                             censor_group, beta, baseline_jump)
      real(dp), intent(in) :: time(:) !! Follow-up times sorted in nondecreasing order.
      integer, intent(in) :: event_code(:) !! Zero for censoring, one for target failure, and two for competing failure.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect covariates, with one row per subject.
      real(dp), intent(in) :: x2(:, :) !! Time-varying-effect covariates, with one row per subject.
      real(dp), intent(in) :: time_functions(:, :) !! Time-function values at distinct target-failure times.
      real(dp), intent(in) :: censoring_survival(:, :) !! Left-continuous censoring KM values by group and sorted time.
      integer, intent(in) :: censor_group(:) !! Consecutive one-based censoring-group code for each subject.
      real(dp), intent(in) :: beta(:) !! Fitted Fine-Gray regression coefficient vector.
      real(dp), intent(out) :: baseline_jump(:) !! Baseline cumulative subdistribution-hazard jumps at target failure times.

      integer :: i
      integer :: iuc
      integer :: itmp
      integer :: ldf
      integer :: n
      real(dp) :: cft
      real(dp) :: linear_predictor
      real(dp) :: risk_sum
      real(dp) :: target_count
      real(dp) :: weight
      real(dp), allocatable :: covariate(:)

      n = size(time)
      allocate(covariate(size(beta)))
      baseline_jump = 0.0_dp
      iuc = 1
      ldf = 0
      do
         itmp = iuc
         do i = iuc, n
            itmp = i
            if (event_code(i) == 1) then
               cft = time(i)
               exit
            end if
         end do
         if (event_code(itmp) /= 1) return
         iuc = itmp
         ldf = ldf + 1
         target_count = 0.0_dp
         do i = iuc, n
            if (time(i) > cft) exit
            itmp = i
            if (event_code(i) == 1) target_count = target_count + 1.0_dp
         end do
         iuc = itmp
         risk_sum = 0.0_dp
         do i = 1, n
            if (time(i) < cft) then
               if (event_code(i) <= 1) cycle
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)*censoring_survival(censor_group(i), iuc) / &
                        censoring_survival(censor_group(i), i)
            else
               call covariate_at(i, ldf, x, x2, time_functions, beta, linear_predictor, covariate)
               weight = exp(linear_predictor)
            end if
            risk_sum = risk_sum + weight
         end do
         baseline_jump(ldf) = baseline_jump(ldf) + target_count/risk_sum
         iuc = iuc + 1
         if (iuc > n) return
      end do
   end subroutine crr_baseline_jumps_kernel

end module cmprsk_crr_kernels
