! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from riskSimul 0.1.2 by Wolfgang Hormann and Ismail Basoglu.
module risksimul_simulation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ghyp_kinds, only : dp, i8
   use ghyp_rng, only : rng_state, seed_rng, uniform_rng, normal_rng, gamma_rng
   use ghyp_special, only : normal_quantile
   use risksimul_types, only : portfolio_model, estimate_row, simulation_result, &
      sis_control, importance_parameters, allocation_result, objective_mse, &
      objective_msre, objective_max_error, objective_max_relative
   use risksimul_math, only : gamma_quantile, orthogonal_completion, &
      integer_allocation, optimal_allocation_heuristic, &
      sample_variance_from_sums, sample_covariance_from_sums
   use risksimul_portfolio, only : portfolio_return_one, tail_loss_response, &
      excess_response, algorithm_3
   implicit none
   private

   real(dp), parameter :: z_975 = 1.959963984540054_dp

   public :: naive_copula, stratified_copula, search_threshold
   public :: NVTCopula, SISTCopula

   interface NVTCopula
      module procedure naive_copula
   end interface NVTCopula

contains

   pure function make_estimate(estimate, variance) result(row)
      real(dp), intent(in) :: estimate, variance
      type(estimate_row) :: row
      row%estimate = estimate
      row%variance = max(0.0_dp,variance)
      row%halfwidth = z_975*sqrt(row%variance)
      row%ci_lower = row%estimate-row%halfwidth
      row%ci_upper = row%estimate+row%halfwidth
      if (abs(row%estimate) > tiny(1.0_dp)) then
         row%relative_error_percent = 100.0_dp*row%halfwidth/abs(row%estimate)
      else
         row%relative_error_percent = huge(1.0_dp)
      end if
   end function make_estimate

   function naive_copula(n, model, thresholds, seed) result(result)
      integer, intent(in) :: n
      type(portfolio_model), intent(in) :: model
      real(dp), intent(in) :: thresholds(:)
      integer(i8), intent(in), optional :: seed
      type(simulation_result) :: result
      type(rng_state) :: rng
      real(dp), allocatable :: z(:), sum_y(:), sum_y2(:), sum_x(:), sum_x2(:), sum_xy(:)
      real(dp) :: y, portfolio_return, indicator, excess
      real(dp) :: mean_y, mean_x, var_y, var_x, cov_xy, cond_est, cond_var
      integer :: i, j, d
      integer(i8) :: local_seed

      if (.not. model%ok) then
         result%message = 'invalid portfolio model'
         return
      else if (n < 2 .or. size(thresholds) < 1 .or. .not. all(ieee_is_finite(thresholds))) then
         result%message = 'sample size must exceed one and thresholds must be finite'
         return
      end if
      d = model%dimension()
      allocate(z(d),sum_y(size(thresholds)),sum_y2(size(thresholds)))
      allocate(sum_x(size(thresholds)),sum_x2(size(thresholds)),sum_xy(size(thresholds)))
      sum_y = 0.0_dp; sum_y2 = 0.0_dp; sum_x = 0.0_dp
      sum_x2 = 0.0_dp; sum_xy = 0.0_dp
      local_seed = 123456789_i8
      if (present(seed)) local_seed = seed
      call seed_rng(rng,local_seed)

      do i = 1, n
         do j = 1, d
            z(j) = normal_rng(rng)
         end do
         y = gamma_rng(0.5_dp*model%copula_df,2.0_dp,rng)
         portfolio_return = portfolio_return_one(z,y,model)
         do j = 1, size(thresholds)
            indicator = tail_loss_response(portfolio_return,thresholds(j))
            excess = excess_response(portfolio_return,thresholds(j))
            sum_y(j) = sum_y(j)+indicator
            sum_y2(j) = sum_y2(j)+indicator*indicator
            sum_x(j) = sum_x(j)+excess
            sum_x2(j) = sum_x2(j)+excess*excess
            sum_xy(j) = sum_xy(j)+indicator*excess
         end do
      end do

      allocate(result%thresholds(size(thresholds)),result%tail_probability(size(thresholds)))
      allocate(result%unconditional_excess(size(thresholds)), &
         result%conditional_excess(size(thresholds)))
      result%thresholds = thresholds
      do j = 1, size(thresholds)
         mean_y = sum_y(j)/real(n,dp)
         mean_x = sum_x(j)/real(n,dp)
         var_y = sample_variance_from_sums(sum_y(j),sum_y2(j),n)/real(n,dp)
         var_x = sample_variance_from_sums(sum_x(j),sum_x2(j),n)/real(n,dp)
         cov_xy = sample_covariance_from_sums(sum_y(j),sum_x(j),sum_xy(j),n)/real(n,dp)
         result%tail_probability(j) = make_estimate(mean_y,var_y)
         result%unconditional_excess(j) = make_estimate(mean_x,var_x)
         if (mean_y > tiny(1.0_dp)) then
            cond_est = mean_x/mean_y-(mean_x*var_y/mean_y+cov_xy)/(mean_y*mean_y)
            cond_var = var_y*mean_x*mean_x/mean_y**4- &
               2.0_dp*mean_x*cov_xy/mean_y**3+var_x/mean_y**2
            result%conditional_excess(j) = make_estimate(cond_est,cond_var)
         else
            result%conditional_excess(j) = make_estimate(0.0_dp,0.0_dp)
         end if
      end do
      result%samples_requested = n
      result%samples_used = n
      result%ok = .true.
   end function naive_copula

   function SISTCopula(n, npilot, model, thresholds, stratasize, CEopt, beta, &
      mintype, seed) result(result)
      integer, intent(in) :: n
      integer, intent(in) :: npilot(:)
      type(portfolio_model), intent(in) :: model
      real(dp), intent(in) :: thresholds(:)
      integer, intent(in), optional :: stratasize(2)
      logical, intent(in), optional :: CEopt
      real(dp), intent(in), optional :: beta
      integer, intent(in), optional :: mintype
      integer(i8), intent(in), optional :: seed
      type(simulation_result) :: result
      type(sis_control) :: control
      integer :: remaining

      if (sum(npilot) >= n .or. any(npilot <= 0)) then
         result%message = 'pilot allocations must be positive and sum to less than n'
         return
      end if
      remaining = n-sum(npilot)
      allocate(control%allocations(size(npilot)+1))
      control%allocations(1:size(npilot)) = npilot
      control%allocations(size(npilot)+1) = remaining
      if (present(stratasize)) then
         control%normal_strata = stratasize(1)
         control%gamma_strata = stratasize(2)
      end if
      if (present(CEopt)) control%optimize_conditional_excess = CEopt
      if (present(beta)) control%intermediate_weight = beta
      if (present(mintype)) control%multi_objective = mintype
      if (present(seed)) control%seed = seed
      result = stratified_copula(model,thresholds,control)
      result%samples_requested = n
   end function SISTCopula

   function stratified_copula(model, thresholds, control) result(result)
      type(portfolio_model), intent(in) :: model
      real(dp), intent(in) :: thresholds(:)
      type(sis_control), intent(in) :: control
      type(simulation_result) :: result
      type(importance_parameters) :: importance
      type(rng_state) :: rng
      real(dp), allocatable :: basis(:,:), direction(:), z_rot(:), z(:)
      real(dp), allocatable :: sum_y(:,:,:), sum_y2(:,:,:), sum_x(:,:,:)
      real(dp), allocatable :: sum_x2(:,:,:), sum_xy(:,:,:)
      real(dp), allocatable :: mean_y(:,:,:), mean_x(:,:,:), var_y(:,:,:)
      real(dp), allocatable :: var_x(:,:,:), cov_xy(:,:,:), fractions(:)
      integer, allocatable :: counts(:,:), batch_counts(:)
      real(dp) :: target_threshold, theta, mu, u1, u2, y, weight_value
      real(dp) :: portfolio_return, indicator, excess, p_stratum
      integer :: d, jn, i1, i2, iteration, i, j, k, sample, flat, used
      logical :: basis_ok

      if (.not. model%ok) then
         result%message = 'invalid portfolio model'
         return
      else if (.not. allocated(control%allocations) .or. &
               any(control%allocations <= 0)) then
         result%message = 'SIS allocations must be allocated and positive'
         return
      else if (size(thresholds) < 1 .or. .not. all(ieee_is_finite(thresholds))) then
         result%message = 'thresholds must be finite and nonempty'
         return
      else if (control%normal_strata < 1 .or. control%gamma_strata < 1 .or. &
               control%minimum_per_stratum < 2) then
         result%message = 'strata counts must be positive and minimum count at least two'
         return
      end if

      d = model%dimension()
      jn = size(thresholds)
      i1 = control%normal_strata
      i2 = control%gamma_strata
      if (jn == 1) then
         target_threshold = thresholds(1)
      else
         target_threshold = maxval(thresholds)*control%intermediate_weight+ &
            minval(thresholds)*(1.0_dp-control%intermediate_weight)
      end if
      importance = algorithm_3(model,target_threshold,control%direction_iterations)
      if (.not. importance%ok) then
         result%message = trim(importance%message)
         return
      end if
      theta = importance%gamma_mean/(0.5_dp*model%copula_df-1.0_dp)
      mu = sqrt(dot_product(importance%shift,importance%shift))
      if (theta <= 0.0_dp .or. mu <= tiny(1.0_dp)) then
         result%message = 'invalid importance-sampling shift or gamma scale'
         return
      end if
      direction = importance%shift/mu
      call orthogonal_completion(direction,basis,basis_ok)
      if (.not. basis_ok) then
         result%message = 'could not construct the stratification basis'
         return
      end if

      allocate(z_rot(d),z(d),counts(i1,i2))
      allocate(sum_y(jn,i1,i2),sum_y2(jn,i1,i2),sum_x(jn,i1,i2))
      allocate(sum_x2(jn,i1,i2),sum_xy(jn,i1,i2))
      allocate(mean_y(jn,i1,i2),mean_x(jn,i1,i2),var_y(jn,i1,i2))
      allocate(var_x(jn,i1,i2),cov_xy(jn,i1,i2),fractions(i1*i2))
      counts = 0
      sum_y = 0.0_dp; sum_y2 = 0.0_dp; sum_x = 0.0_dp
      sum_x2 = 0.0_dp; sum_xy = 0.0_dp
      fractions = 1.0_dp/real(i1*i2,dp)
      p_stratum = 1.0_dp/real(i1*i2,dp)
      used = 0
      call seed_rng(rng,control%seed)

      do iteration = 1, size(control%allocations)
         if (iteration > 1) then
            call stratum_moments(counts,sum_y,sum_y2,sum_x,sum_x2,sum_xy, &
               mean_y,var_y,mean_x,var_x,cov_xy)
            fractions = allocation_fractions(mean_y,var_y,mean_x,var_x,cov_xy, &
               p_stratum,control)
         end if
         call integer_allocation(fractions,control%allocations(iteration), &
            control%minimum_per_stratum,batch_counts)
         flat = 0
         do i = 1, i1
            do j = 1, i2
               flat = flat+1
               do sample = 1, batch_counts(flat)
                  u1 = (real(i-1,dp)+uniform_rng(rng))/real(i1,dp)
                  z_rot(1) = mu+normal_quantile(u1)
                  do k = 2, d
                     z_rot(k) = normal_rng(rng)
                  end do
                  u2 = (real(j-1,dp)+uniform_rng(rng))/real(i2,dp)
                  y = gamma_quantile(u2,0.5_dp*model%copula_df,theta)
                  weight_value = exp((0.5_dp*mu-z_rot(1))*mu+ &
                     (2.0_dp-theta)*y/(2.0_dp*theta))* &
                     (theta/2.0_dp)**(0.5_dp*model%copula_df)
                  z = matmul(basis,z_rot)
                  portfolio_return = portfolio_return_one(z,y,model)
                  counts(i,j) = counts(i,j)+1
                  do k = 1, jn
                     indicator = weight_value*tail_loss_response(portfolio_return,thresholds(k))
                     excess = weight_value*excess_response(portfolio_return,thresholds(k))
                     sum_y(k,i,j) = sum_y(k,i,j)+indicator
                     sum_y2(k,i,j) = sum_y2(k,i,j)+indicator*indicator
                     sum_x(k,i,j) = sum_x(k,i,j)+excess
                     sum_x2(k,i,j) = sum_x2(k,i,j)+excess*excess
                     sum_xy(k,i,j) = sum_xy(k,i,j)+indicator*excess
                  end do
               end do
               used = used+batch_counts(flat)
            end do
         end do
      end do

      call stratum_moments(counts,sum_y,sum_y2,sum_x,sum_x2,sum_xy, &
         mean_y,var_y,mean_x,var_x,cov_xy)
      call finalize_stratified_result(result,thresholds,counts,mean_y,var_y, &
         mean_x,var_x,cov_xy,p_stratum)
      result%samples_requested = sum(control%allocations)
      result%samples_used = used
      result%ok = .true.
   end function stratified_copula

   subroutine stratum_moments(counts, sum_y, sum_y2, sum_x, sum_x2, sum_xy, &
      mean_y, var_y, mean_x, var_x, cov_xy)
      integer, intent(in) :: counts(:,:)
      real(dp), intent(in) :: sum_y(:,:,:), sum_y2(:,:,:), sum_x(:,:,:)
      real(dp), intent(in) :: sum_x2(:,:,:), sum_xy(:,:,:)
      real(dp), intent(out) :: mean_y(:,:,:), var_y(:,:,:), mean_x(:,:,:)
      real(dp), intent(out) :: var_x(:,:,:), cov_xy(:,:,:)
      integer :: i, j, k, n

      do k = 1, size(sum_y,1)
         do i = 1, size(counts,1)
            do j = 1, size(counts,2)
               n = counts(i,j)
               mean_y(k,i,j) = sum_y(k,i,j)/real(max(n,1),dp)
               mean_x(k,i,j) = sum_x(k,i,j)/real(max(n,1),dp)
               var_y(k,i,j) = sample_variance_from_sums(sum_y(k,i,j), &
                  sum_y2(k,i,j),n)
               var_x(k,i,j) = sample_variance_from_sums(sum_x(k,i,j), &
                  sum_x2(k,i,j),n)
               cov_xy(k,i,j) = sample_covariance_from_sums(sum_y(k,i,j), &
                  sum_x(k,i,j),sum_xy(k,i,j),n)
            end do
         end do
      end do
   end subroutine stratum_moments

   function allocation_fractions(mean_y, var_y, mean_x, var_x, cov_xy, &
      p_stratum, control) result(fractions)
      real(dp), intent(in) :: mean_y(:,:,:), var_y(:,:,:), mean_x(:,:,:)
      real(dp), intent(in) :: var_x(:,:,:), cov_xy(:,:,:), p_stratum
      type(sis_control), intent(in) :: control
      real(dp), allocatable :: fractions(:)
      real(dp), allocatable :: source_var(:,:,:), a(:,:), b(:), scores(:)
      real(dp) :: x_est, y_est, denom
      integer :: jn, i1, i2, k, i, j, flat, target
      type(allocation_result) :: allocation

      jn = size(mean_y,1); i1 = size(mean_y,2); i2 = size(mean_y,3)
      allocate(source_var(jn,i1,i2),a(jn,i1*i2),b(jn),scores(i1*i2))
      if (control%optimize_conditional_excess) then
         do k = 1, jn
            x_est = p_stratum*sum(mean_x(k,:,:))
            y_est = p_stratum*sum(mean_y(k,:,:))
            if (abs(y_est) <= tiny(1.0_dp)) then
               source_var(k,:,:) = var_x(k,:,:)
            else
               source_var(k,:,:) = x_est*x_est*var_y(k,:,:)/y_est**4- &
                  2.0_dp*x_est*cov_xy(k,:,:)/y_est**3+var_x(k,:,:)/y_est**2
               source_var(k,:,:) = max(source_var(k,:,:),0.0_dp)
            end if
         end do
      else
         source_var = max(var_y,0.0_dp)
      end if

      do k = 1, jn
         y_est = p_stratum*sum(mean_y(k,:,:))
         x_est = p_stratum*sum(mean_x(k,:,:))
         if (control%optimize_conditional_excess) then
            if (abs(x_est) > tiny(1.0_dp)) then
               b(k) = 1.0_dp/(x_est*x_est)
            else
               b(k) = 0.0_dp
            end if
         else
            if (abs(y_est) > tiny(1.0_dp)) then
               b(k) = 1.0_dp/(y_est*y_est)
            else
               b(k) = 0.0_dp
            end if
         end if
         flat = 0
         do i = 1, i1
            do j = 1, i2
               flat = flat+1
               a(k,flat) = source_var(k,i,j)*p_stratum*p_stratum
            end do
         end do
      end do

      select case(control%multi_objective)
      case(objective_mse)
         scores = sqrt(max(sum(a,dim=1),0.0_dp))
      case(objective_msre)
         scores = sqrt(max(matmul(transpose(a),b),0.0_dp))
      case(objective_max_error)
         allocation = optimal_allocation_heuristic(transpose(a), &
            control%allocation_tolerance,control%upstream_allocation_compatibility)
         if (allocation%ok) then
            fractions = allocation%fractions
            return
         end if
         scores = 1.0_dp
      case(objective_max_relative)
         do k = 1, jn
            a(k,:) = b(k)*a(k,:)
         end do
         allocation = optimal_allocation_heuristic(transpose(a), &
            control%allocation_tolerance,control%upstream_allocation_compatibility)
         if (allocation%ok) then
            fractions = allocation%fractions
            return
         end if
         scores = 1.0_dp
      case default
         target = min(jn,max(1,control%multi_objective))
         scores = sqrt(max(a(target,:),0.0_dp))
      end select
      denom = sum(scores)
      allocate(fractions(i1*i2))
      if (denom <= tiny(1.0_dp)) then
         fractions = 1.0_dp/real(i1*i2,dp)
      else
         fractions = scores/denom
      end if
   end function allocation_fractions

   subroutine finalize_stratified_result(result, thresholds, counts, mean_y, var_y, &
      mean_x, var_x, cov_xy, p_stratum)
      type(simulation_result), intent(inout) :: result
      real(dp), intent(in) :: thresholds(:)
      integer, intent(in) :: counts(:,:)
      real(dp), intent(in) :: mean_y(:,:,:), var_y(:,:,:), mean_x(:,:,:)
      real(dp), intent(in) :: var_x(:,:,:), cov_xy(:,:,:), p_stratum
      real(dp) :: tail_est, tail_var, excess_est, excess_var, cross_cov
      real(dp) :: cond_est, cond_var
      integer :: k, i, j

      allocate(result%thresholds(size(thresholds)),result%tail_probability(size(thresholds)))
      allocate(result%unconditional_excess(size(thresholds)), &
         result%conditional_excess(size(thresholds)))
      result%thresholds = thresholds
      do k = 1, size(thresholds)
         tail_est = p_stratum*sum(mean_y(k,:,:))
         excess_est = p_stratum*sum(mean_x(k,:,:))
         tail_var = 0.0_dp; excess_var = 0.0_dp; cross_cov = 0.0_dp
         do i = 1, size(counts,1)
            do j = 1, size(counts,2)
               tail_var = tail_var+p_stratum*p_stratum*var_y(k,i,j)/real(counts(i,j),dp)
               excess_var = excess_var+p_stratum*p_stratum*var_x(k,i,j)/real(counts(i,j),dp)
               cross_cov = cross_cov+p_stratum*p_stratum*cov_xy(k,i,j)/real(counts(i,j),dp)
            end do
         end do
         result%tail_probability(k) = make_estimate(tail_est,tail_var)
         result%unconditional_excess(k) = make_estimate(excess_est,excess_var)
         if (tail_est > tiny(1.0_dp)) then
            cond_est = excess_est/tail_est- &
               (excess_est*tail_var/tail_est+cross_cov)/(tail_est*tail_est)
            cond_var = tail_var*excess_est*excess_est/tail_est**4- &
               2.0_dp*excess_est*cross_cov/tail_est**3+excess_var/tail_est**2
            result%conditional_excess(k) = make_estimate(cond_est,cond_var)
         else
            result%conditional_excess(k) = make_estimate(0.0_dp,0.0_dp)
         end if
      end do
   end subroutine finalize_stratified_result

   function search_threshold(n, probability, model, lower, upper, seed, &
      relative_tolerance) result(threshold)
      integer, intent(in) :: n
      real(dp), intent(in) :: probability
      type(portfolio_model), intent(in) :: model
      real(dp), intent(in), optional :: lower, upper, relative_tolerance
      integer(i8), intent(in), optional :: seed
      real(dp) :: threshold
      real(dp) :: lo, hi, mid, tolerance, estimate
      integer :: iter
      integer(i8) :: base_seed
      type(simulation_result) :: fit

      lo = 0.85_dp; hi = 1.0_dp; tolerance = 0.05_dp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      if (present(relative_tolerance)) tolerance = relative_tolerance
      base_seed = 123456789_i8
      if (present(seed)) base_seed = seed
      do iter = 1, 80
         mid = 0.5_dp*(lo+hi)
         fit = naive_copula(n,model,[mid],base_seed+int(iter-1,i8))
         if (.not. fit%ok) exit
         estimate = fit%tail_probability(1)%estimate
         if (abs(estimate-probability)/max(probability,tiny(1.0_dp)) <= tolerance) exit
         if (estimate < probability) then
            lo = mid
         else
            hi = mid
         end if
      end do
      threshold = mid
   end function search_threshold

end module risksimul_simulation
