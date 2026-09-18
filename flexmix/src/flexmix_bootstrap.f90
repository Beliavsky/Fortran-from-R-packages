! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_bootstrap
   use flexmix_kinds, only : dp
   use flexmix_types, only : flexmix_control, flexmix_result, flexmix_boot_result
   use flexmix_types, only : flexmix_model_gaussian, flexmix_model_poisson, flexmix_model_binomial
   use flexmix_types, only : flexmix_model_gamma_regression
   use flexmix_api, only : flexmix_gaussian, flexmix_poisson, flexmix_binomial, flexmix_gamma
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: flexmix_boot_regression
   public :: flexmix_lr_test_regression

contains

   subroutine flexmix_boot_regression(template, x, y, r, k_values, simulation, boot, info, control, trials, offset, case_weights, &
                                      initialize_solution)
      type(flexmix_result), intent(in) :: template !! Fitted ordinary regression mixture supplying the parametric bootstrap law.
      real(dp), intent(in) :: x(:,:) !! Original numeric regression design matrix, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Original response vector, size `n`; binomial entries are success counts.
      integer, intent(in) :: r !! Number of bootstrap rows including the first observed-data fit; must be at least one.
      integer, intent(in) :: k_values(:) !! Component counts fitted to every bootstrap data set; all values must be positive.
      character(len=*), intent(in) :: simulation !! Bootstrap mode: `parametric`, `ordinary`, or `empirical`; the two.
      type(flexmix_boot_result), intent(out) :: boot !! Bootstrap log likelihoods, retained component counts, convergence flags.
      integer, intent(out) :: info !! Zero on success; nonzero for dimensions, unsupported models, or options.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls used for every refitted mixture.
      real(dp), intent(in), optional :: trials(:) !! Binomial trial counts, size `n`; omitted values imply Bernoulli trials.
      real(dp), intent(in), optional :: offset(:) !! Fixed regression offset, size `n`, retained or resampled with each design.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case weights, size `n`; defaults to weights.
      logical, intent(in), optional :: initialize_solution !! If true, reuse the supplied fit at `template%k`.
      type(flexmix_control) :: ctl
      type(flexmix_result) :: fit
      real(dp), allocatable :: xb(:,:), yb(:), trialsb(:), offsetb(:), weightsb(:), base_trials(:), base_offset(:), base_weights(:)
      integer, allocatable :: index(:)
      logical :: initialize
      character(len=:), allocatable :: mode
      integer :: i, j, n, p, fit_info

      info = 0
      n = size(y)
      p = size(x,2)
      mode = trim(adjustl(simulation))
      initialize = .false.
      if (present(initialize_solution)) initialize = initialize_solution
      ctl = flexmix_control()
      if (present(control)) ctl = control
      if (size(x,1) /= n .or. p /= template%p .or. r < 1 .or. size(k_values) < 1 .or. any(k_values < 1)) then
         info = 1
         return
      end if
      if (mode /= 'parametric' .and. mode /= 'ordinary' .and. mode /= 'empirical') then
         info = 2
         return
      end if
      select case (template%model_kind)
      case (flexmix_model_gaussian, flexmix_model_poisson, flexmix_model_binomial, flexmix_model_gamma_regression)
      case default
         info = 3
         return
      end select
      if (allocated(template%concomitant_coef)) then
         info = 4
         return
      end if
      if (allocated(template%group_first)) then
         if (any(.not. template%group_first)) then
            info = 5
            return
         end if
      end if
      allocate(base_trials(n), base_offset(n), base_weights(n))
      base_trials = 1.0_dp
      if (present(trials)) then
         if (size(trials) /= n .or. any(trials < 0.0_dp)) then
            info = 6
            return
         end if
         base_trials = trials
      end if
      if (template%model_kind == flexmix_model_binomial) then
         if (any(y < 0.0_dp) .or. any(y > base_trials) .or. any(abs(base_trials - real(nint(base_trials),dp)) > 1.0e-8_dp)) then
            info = 7
            return
         end if
      end if
      base_offset = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) then
            info = 8
            return
         end if
         base_offset = offset
      end if
      base_weights = 1.0_dp
      if (allocated(template%case_weights)) then
         if (size(template%case_weights) == n) base_weights = template%case_weights
      end if
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            info = 9
            return
         end if
         base_weights = case_weights
      end if
      allocate(boot%requested_k(size(k_values)), boot%loglik(r,size(k_values)), boot%fitted_k(r,size(k_values)))
      allocate(boot%converged(r,size(k_values)), boot%models(r,size(k_values)))
      allocate(xb(n,p), yb(n), trialsb(n), offsetb(n), weightsb(n), index(n))
      boot%nrep = r
      boot%requested_k = k_values
      boot%loglik = -huge(1.0_dp)
      boot%fitted_k = 0
      boot%converged = .false.
      do i = 1, r
         if (i == 1) then
            xb = x
            yb = y
            trialsb = base_trials
            offsetb = base_offset
            weightsb = base_weights
         else if (mode == 'parametric') then
            xb = x
            trialsb = base_trials
            offsetb = base_offset
            weightsb = base_weights
            call simulate_from_regression(template, xb, trialsb, offsetb, yb, fit_info)
            if (fit_info /= 0) then
               info = 20 + fit_info
               return
            end if
         else
            call sample_rows(n, index)
            do j = 1, n
               xb(j,:) = x(index(j),:)
               yb(j) = y(index(j))
               trialsb(j) = base_trials(index(j))
               offsetb(j) = base_offset(index(j))
               weightsb(j) = base_weights(index(j))
            end do
         end if
         do j = 1, size(k_values)
            if (i == 1 .and. initialize .and. k_values(j) == template%k) then
               fit = template
            else
               call fit_boot_model(template%model_kind, xb, yb, trialsb, offsetb, weightsb, k_values(j), ctl, fit, fit_info)
               if (fit_info /= 0) cycle
            end if
            boot%models(i,j) = fit
            boot%loglik(i,j) = fit%loglik
            boot%fitted_k(i,j) = fit%k
            boot%converged(i,j) = fit%converged
         end do
      end do
   end subroutine flexmix_boot_regression

   subroutine flexmix_lr_test_regression(object, x, y, r, alternative, statistic, p_value, bootstrap_statistics, boot, info, &
                                         control, trials, offset, case_weights)
      type(flexmix_result), intent(in) :: object !! Fitted null regression mixture whose adjacent component count is tested.
      real(dp), intent(in) :: x(:,:) !! Numeric regression design matrix used by the fitted null model, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Original response vector, size `n`; binomial entries are success counts.
      integer, intent(in) :: r !! Parametric bootstrap count including the observed-data row; must be at least one.
      character(len=*), intent(in) :: alternative !! Adjacent component-count alternative, `greater` or `less`.
      real(dp), intent(out) :: statistic !! Observed bootstrap likelihood-ratio statistic `2*(logLik_high-logLik_low)`.
      real(dp), intent(out) :: p_value !! Bootstrap tail probability including the observed statistic itself.
      real(dp), allocatable, intent(out) :: bootstrap_statistics(:) !! Valid likelihood-ratio statistics, with the observed.
      type(flexmix_boot_result), intent(out) :: boot !! Full two-model parametric bootstrap results used by the test.
      integer, intent(out) :: info !! Zero on success; nonzero for an invalid alternative or failed bootstrap.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls used for both fitted component counts.
      real(dp), intent(in), optional :: trials(:) !! Optional binomial trial counts, size `n`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed regression offset, size `n`.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case weights, size `n`.
      integer :: kvals(2), i, nvalid, low_col, high_col, boot_info
      logical, allocatable :: valid(:)
      real(dp), allocatable :: all_statistics(:)
      character(len=:), allocatable :: alt

      statistic = 0.0_dp
      p_value = 1.0_dp
      allocate(bootstrap_statistics(0))
      alt = trim(adjustl(alternative))
      select case (alt)
      case ('greater')
         kvals = [object%k, object%k + 1]
      case ('less')
         if (object%k <= 1) then
            info = 1
            return
         end if
         kvals = [object%k, object%k - 1]
      case default
         info = 2
         return
      end select
      call flexmix_boot_regression(object, x, y, r, kvals, 'parametric', boot, boot_info, control, trials, offset, &
                                   case_weights, initialize_solution=.true.)
      if (boot_info /= 0) then
         info = 10 + boot_info
         return
      end if
      if (kvals(1) < kvals(2)) then
         low_col = 1
         high_col = 2
      else
         low_col = 2
         high_col = 1
      end if
      allocate(valid(r), all_statistics(r))
      valid = .false.
      all_statistics = 0.0_dp
      do i = 1, r
         valid(i) = boot%fitted_k(i,1) == kvals(1) .and. boot%fitted_k(i,2) == kvals(2)
         valid(i) = valid(i) .and. boot%loglik(i,low_col) > -0.5_dp * huge(1.0_dp)
         valid(i) = valid(i) .and. boot%loglik(i,high_col) > -0.5_dp * huge(1.0_dp)
         if (valid(i)) all_statistics(i) = 2.0_dp * (boot%loglik(i,high_col) - boot%loglik(i,low_col))
      end do
      if (.not. valid(1)) then
         info = 3
         return
      end if
      nvalid = count(valid)
      deallocate(bootstrap_statistics)
      allocate(bootstrap_statistics(nvalid))
      bootstrap_statistics = pack(all_statistics, valid)
      statistic = all_statistics(1)
      p_value = real(count(bootstrap_statistics >= statistic),dp) / real(nvalid,dp)
      info = 0
   end subroutine flexmix_lr_test_regression

   subroutine fit_boot_model(model_kind, x, y, trials, offset, weights, k, control, fit, info)
      integer, intent(in) :: model_kind !! Regression-family identifier from the fitted template.
      real(dp), intent(in) :: x(:,:) !! Bootstrap regression design matrix, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Bootstrap response vector, size `n`.
      real(dp), intent(in) :: trials(:) !! Binomial trial counts or an ignored size-`n` placeholder for other families.
      real(dp), intent(in) :: offset(:) !! Fixed regression offsets, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative case weights, size `n`.
      integer, intent(in) :: k !! Requested component count for this bootstrap fit.
      type(flexmix_control), intent(in) :: control !! EM controls used for the fit.
      type(flexmix_result), intent(out) :: fit !! Fitted bootstrap mixture result.
      integer, intent(out) :: info !! Zero when the family fit returns a usable status; nonzero otherwise.
      real(dp), allocatable :: failure(:)

      select case (model_kind)
      case (flexmix_model_gaussian)
         call flexmix_gaussian(x, y, k, fit, control, case_weights=weights, offset=offset)
      case (flexmix_model_poisson)
         call flexmix_poisson(x, y, k, fit, control, case_weights=weights, offset=offset)
      case (flexmix_model_binomial)
         allocate(failure(size(y)))
         failure = trials - y
         call flexmix_binomial(x, y, failure, k, fit, control, case_weights=weights, offset=offset)
      case (flexmix_model_gamma_regression)
         call flexmix_gamma(x, y, k, fit, control, case_weights=weights, offset=offset)
      case default
         info = 1
         return
      end select
      if (fit%status /= 0) then
         info = 2
      else
         info = 0
      end if
   end subroutine fit_boot_model

   subroutine simulate_from_regression(template, x, trials, offset, y, info)
      type(flexmix_result), intent(in) :: template !! Fitted ordinary regression mixture defining priors and component response.
      real(dp), intent(in) :: x(:,:) !! Fixed design matrix at which bootstrap responses are generated, shape `(n,p)`.
      real(dp), intent(in) :: trials(:) !! Binomial trial counts, size `n`; ignored for other families.
      real(dp), intent(in) :: offset(:) !! Fixed linear-predictor offset, size `n`.
      real(dp), intent(out) :: y(:) !! Simulated response vector, size `n`.
      integer, intent(out) :: info !! Zero on success; nonzero for missing parameters or unsupported/invalid response laws.
      real(dp) :: eta, mean_value, probability, u, z
      integer :: i, j, component, ntrial, success

      if (.not. allocated(template%beta) .or. .not. allocated(template%prior)) then
         info = 1
         return
      end if
      if (size(x,2) /= size(template%beta,1) .or. size(y) /= size(x,1) .or. size(offset) /= size(y)) then
         info = 2
         return
      end if
      do i = 1, size(y)
         call random_number(u)
         component = template%k
         probability = 0.0_dp
         do j = 1, template%k
            probability = probability + template%prior(j)
            if (u <= probability) then
               component = j
               exit
            end if
         end do
         eta = dot_product(x(i,:), template%beta(:,component)) + offset(i)
         select case (template%model_kind)
         case (flexmix_model_gaussian)
            if (.not. allocated(template%sigma)) then
               info = 3
               return
            end if
            call standard_normal_random(z)
            y(i) = eta + template%sigma(component) * z
         case (flexmix_model_poisson)
            call poisson_random(exp(max(-30.0_dp,min(30.0_dp,eta))), y(i))
         case (flexmix_model_binomial)
            ntrial = nint(trials(i))
            probability = 1.0_dp / (1.0_dp + exp(-max(-30.0_dp,min(30.0_dp,eta))))
            success = 0
            do j = 1, ntrial
               call random_number(u)
               if (u < probability) success = success + 1
            end do
            y(i) = real(success,dp)
         case (flexmix_model_gamma_regression)
            if (.not. allocated(template%shape) .or. eta <= 0.0_dp) then
               info = 4
               return
            end if
            mean_value = 1.0_dp / eta
            call gamma_random(template%shape(component), mean_value / template%shape(component), y(i), success)
            if (success /= 0) then
               info = 5
               return
            end if
         case default
            info = 6
            return
         end select
      end do
      info = 0
   end subroutine simulate_from_regression

   subroutine sample_rows(n, index)
      integer, intent(in) :: n !! Positive number of source rows available for resampling.
      integer, intent(out) :: index(:) !! One-based bootstrap source-row indices, size `n`.
      real(dp) :: u
      integer :: i
      do i = 1, n
         call random_number(u)
         index(i) = min(n, 1 + int(u * real(n,dp)))
      end do
   end subroutine sample_rows

   subroutine standard_normal_random(z)
      real(dp), intent(out) :: z !! One standard-normal random variate generated by the Box-Muller transform.
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   end subroutine standard_normal_random

   subroutine poisson_random(lambda, value)
      real(dp), intent(in) :: lambda !! Nonnegative Poisson mean.
      real(dp), intent(out) :: value !! Generated nonnegative Poisson count represented in the package real kind.
      real(dp) :: limit, product_value, u, z
      integer :: count
      if (lambda <= 0.0_dp) then
         value = 0.0_dp
      else if (lambda < 30.0_dp) then
         limit = exp(-lambda)
         product_value = 1.0_dp
         count = 0
         do
            count = count + 1
            call random_number(u)
            product_value = product_value * u
            if (product_value <= limit) exit
         end do
         value = real(count - 1,dp)
      else
         call standard_normal_random(z)
         value = real(max(0,nint(lambda + sqrt(lambda) * z)),dp)
      end if
   end subroutine poisson_random

   recursive subroutine gamma_random(shape, scale, value, info)
      real(dp), intent(in) :: shape !! Positive Gamma shape parameter.
      real(dp), intent(in) :: scale !! Positive Gamma scale parameter.
      real(dp), intent(out) :: value !! Generated positive Gamma random variate.
      integer, intent(out) :: info !! Zero on success; nonzero for nonpositive parameters.
      real(dp) :: d, c, u, v, z, base
      integer :: inner_info
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         value = 0.0_dp
         info = 1
         return
      end if
      if (shape < 1.0_dp) then
         call gamma_random(shape + 1.0_dp, 1.0_dp, base, inner_info)
         if (inner_info /= 0) then
            value = 0.0_dp
            info = inner_info
            return
         end if
         call random_number(u)
         value = scale * base * max(u,tiny(1.0_dp))**(1.0_dp / shape)
         info = 0
         return
      end if
      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         call standard_normal_random(z)
         v = 1.0_dp + c * z
         if (v <= 0.0_dp) cycle
         v = v * v * v
         call random_number(u)
         if (u < 1.0_dp - 0.0331_dp * z**4) exit
         if (log(max(u,tiny(1.0_dp))) < 0.5_dp * z*z + d * (1.0_dp - v + log(v))) exit
      end do
      value = scale * d * v
      info = 0
   end subroutine gamma_random

end module flexmix_bootstrap
