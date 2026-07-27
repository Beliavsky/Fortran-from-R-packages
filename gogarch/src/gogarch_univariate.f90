! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_univariate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gogarch_kinds, only : dp
   use gogarch_optimizer, only : optimizer_result, nelder_mead
   use gogarch_distributions, only : distribution_is_valid, innovation_logpdf, random_innovation
   use gogarch_distributions, only : innovation_asym_power_moment
   use gogarch_types, only : garch11_fit, univariate_spec
   implicit none
   private

   type :: garch_context
      real(dp), allocatable :: y(:)
      type(univariate_spec) :: spec
   end type garch_context

   public :: filter_aparch, filter_garchpq, filter_garch11
   public :: fit_univariate, fit_garchpq, fit_garch11
   public :: forecast_univariate, forecast_garch11
   public :: simulate_aparch, simulate_garchpq, simulate_garch11
   public :: validate_specification

contains

   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code+32)
      end do
   end function lower_string

   pure function validate_specification(spec) result(ok)
      type(univariate_spec), intent(in) :: spec
      logical :: ok
      character(len=:), allocatable :: model
      model = trim(adjustl(lower_string(spec%model)))
      ok = (model == 'garch' .or. model == 'aparch') .and. spec%p >= 1 .and. spec%q >= 0 .and. &
         spec%o >= 0 .and. spec%o <= spec%p .and. distribution_is_valid(spec%distribution,spec%shape,spec%skew)
      if (model == 'garch') ok = ok .and. spec%o == 0
      if (model == 'aparch') ok = ok .and. spec%delta > 0.0_dp
   end function validate_specification

   subroutine filter_aparch(y, mean, omega, arch, leverage, garch, delta, distribution, shape, skew, residuals, &
      power_scale, variance, standardized, log_likelihood, valid)
      real(dp), intent(in) :: y(:), mean, omega, arch(:), leverage(:), garch(:), delta, shape, skew
      character(len=*), intent(in) :: distribution
      real(dp), intent(out) :: residuals(size(y)), power_scale(size(y)), variance(size(y)), standardized(size(y))
      real(dp), intent(out), optional :: log_likelihood
      logical, intent(out), optional :: valid
      real(dp) :: initial_variance, initial_power, ll, sample_var, persistence, shock, gamma_i
      integer :: t, i, lag, n, p, o, q
      logical :: ok
      n = size(y)
      p = size(arch)
      o = size(leverage)
      q = size(garch)
      persistence = sum(arch)+sum(garch)
      ok = n >= 2 .and. p >= 1 .and. o <= p .and. omega > 0.0_dp .and. all(arch >= 0.0_dp) .and. &
         all(garch >= 0.0_dp) .and. all(abs(leverage) < 1.0_dp) .and. delta > 0.0_dp .and. &
         persistence < 0.9995_dp .and. distribution_is_valid(distribution,shape,skew)
      if (.not. ok) then
         residuals = 0.0_dp
         power_scale = huge(1.0_dp)
         variance = huge(1.0_dp)
         standardized = 0.0_dp
         if (present(log_likelihood)) log_likelihood = -huge(1.0_dp)
         if (present(valid)) valid = .false.
         return
      end if
      residuals = y-mean
      sample_var = sum((y-sum(y)/real(n,dp))**2)/real(n,dp)
      initial_variance = max(sample_var,(omega/max(1.0_dp-persistence,1.0e-6_dp))**(2.0_dp/delta))
      initial_power = max(initial_variance**(0.5_dp*delta),1.0e-12_dp)
      power_scale(1) = initial_power
      do t = 2, n
         power_scale(t) = omega
         do i = 1, p
            lag = t-i
            gamma_i = 0.0_dp
            if (i <= o) gamma_i = leverage(i)
            if (lag >= 1) then
               shock = (abs(residuals(lag))-gamma_i*residuals(lag))**delta
            else
               shock = initial_power
            end if
            power_scale(t) = power_scale(t)+arch(i)*shock
         end do
         do i = 1, q
            lag = t-i
            if (lag >= 1) then
               power_scale(t) = power_scale(t)+garch(i)*power_scale(lag)
            else
               power_scale(t) = power_scale(t)+garch(i)*initial_power
            end if
         end do
         power_scale(t) = max(power_scale(t),1.0e-14_dp)
      end do
      variance = power_scale**(2.0_dp/delta)
      standardized = residuals/sqrt(variance)
      ll = 0.0_dp
      do t = 1, n
         ll = ll+innovation_logpdf(standardized(t),distribution,shape,skew)-0.5_dp*log(variance(t))
      end do
      ok = ieee_is_finite(ll) .and. all(ieee_is_finite(variance))
      if (present(log_likelihood)) log_likelihood = ll
      if (present(valid)) valid = ok
   end subroutine filter_aparch

   subroutine filter_garchpq(y, mean, omega, arch, garch, distribution, shape, skew, residuals, variance, &
      standardized, log_likelihood, valid)
      real(dp), intent(in) :: y(:), mean, omega, arch(:), garch(:), shape, skew
      character(len=*), intent(in) :: distribution
      real(dp), intent(out) :: residuals(size(y)), variance(size(y)), standardized(size(y))
      real(dp), intent(out), optional :: log_likelihood
      logical, intent(out), optional :: valid
      real(dp) :: power_scale(size(y))
      call filter_aparch(y,mean,omega,arch,[real(dp)::],garch,2.0_dp,distribution,shape,skew,residuals, &
         power_scale,variance,standardized,log_likelihood,valid)
   end subroutine filter_garchpq

   subroutine filter_garch11(y, mean, omega, alpha, beta, residuals, variance, standardized, log_likelihood, valid)
      real(dp), intent(in) :: y(:), mean, omega, alpha, beta
      real(dp), intent(out) :: residuals(size(y)), variance(size(y)), standardized(size(y))
      real(dp), intent(out), optional :: log_likelihood
      logical, intent(out), optional :: valid
      call filter_garchpq(y,mean,omega,[alpha],[beta],'norm',8.0_dp,1.0_dp,residuals,variance,standardized, &
         log_likelihood,valid)
   end subroutine filter_garch11

   function fit_univariate(y, specification, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      type(univariate_spec), intent(in), optional :: specification
      integer, intent(in), optional :: max_iterations
      type(garch11_fit) :: fit
      type(garch_context) :: context
      type(optimizer_result) :: opt
      type(univariate_spec) :: spec
      real(dp), allocatable :: x0(:), arch(:), leverage(:), garch(:)
      real(dp) :: mean, omega, delta, shape, skew
      integer :: maxit, n, npar
      logical :: ok
      spec = univariate_spec()
      if (present(specification)) spec = specification
      call normalize_spec(spec)
      n = size(y)
      if (.not. validate_specification(spec) .or. n < max(20,2*(spec%p+spec%q+spec%o)+5)) then
         fit%status = 8
         return
      end if
      maxit = 500
      if (present(max_iterations)) maxit = max(30,max_iterations)
      allocate(context%y(n))
      context%y = y
      context%spec = spec
      npar = parameter_count(spec)
      allocate(x0(npar))
      call initial_parameters(y,spec,x0)
      opt = nelder_mead(garch_objective,x0,context,step=0.12_dp,tolerance=1.0e-8_dp,max_iterations=maxit)
      allocate(arch(spec%p),leverage(spec%o),garch(spec%q))
      call decode_parameters(opt%x,spec,mean,omega,arch,leverage,garch,delta,shape,skew,ok)
      fit%model = spec%model
      fit%distribution = spec%distribution
      fit%p = spec%p
      fit%o = spec%o
      fit%q = spec%q
      fit%mean = mean
      fit%omega = omega
      fit%delta = delta
      fit%shape = shape
      fit%skew = skew
      fit%iterations = opt%iterations
      fit%status = opt%status
      allocate(fit%arch(spec%p),fit%leverage(spec%o),fit%garch(spec%q))
      fit%arch = arch
      fit%leverage = leverage
      fit%garch = garch
      fit%alpha = arch(1)
      if (spec%q > 0) fit%beta = garch(1)
      allocate(fit%residuals(n),fit%power_scale(n),fit%variance(n),fit%sigma(n),fit%standardized(n))
      call filter_aparch(y,mean,omega,arch,leverage,garch,delta,spec%distribution,shape,skew,fit%residuals, &
         fit%power_scale,fit%variance,fit%standardized,fit%log_likelihood,ok)
      fit%sigma = sqrt(fit%variance)
      if (.not. ok) fit%status = 2
   end function fit_univariate

   function fit_garchpq(y, p, q, distribution, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: p, q
      character(len=*), intent(in), optional :: distribution
      integer, intent(in), optional :: max_iterations
      type(garch11_fit) :: fit
      type(univariate_spec) :: spec
      spec%model = 'garch'
      spec%p = p
      spec%o = 0
      spec%q = q
      if (present(distribution)) spec%distribution = distribution
      fit = fit_univariate(y,spec,max_iterations)
   end function fit_garchpq

   function fit_garch11(y, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in), optional :: max_iterations
      type(garch11_fit) :: fit
      fit = fit_garchpq(y,1,1,'norm',max_iterations)
   end function fit_garch11

   subroutine forecast_univariate(fit, n_ahead, mean_forecast, variance_forecast, power_forecast)
      type(garch11_fit), intent(in) :: fit
      integer, intent(in) :: n_ahead
      real(dp), intent(out) :: mean_forecast(n_ahead), variance_forecast(n_ahead)
      real(dp), intent(out), optional :: power_forecast(n_ahead)
      real(dp), allocatable :: extended_power(:), shock_history(:,:)
      real(dp) :: gamma_i, moment
      integer :: k, i, index, lag, n, total
      if (n_ahead < 1 .or. .not. allocated(fit%power_scale)) return
      n = size(fit%power_scale)
      total = n+n_ahead
      allocate(extended_power(total),shock_history(total,fit%p))
      extended_power(1:n) = fit%power_scale
      shock_history = 0.0_dp
      do index = 1, n
         do i = 1, fit%p
            gamma_i = 0.0_dp
            if (i <= fit%o) gamma_i = fit%leverage(i)
            shock_history(index,i) = (abs(fit%residuals(index))-gamma_i*fit%residuals(index))**fit%delta
         end do
      end do
      do k = 1, n_ahead
         index = n+k
         extended_power(index) = fit%omega
         do i = 1, fit%p
            lag = index-i
            if (lag <= n) then
               extended_power(index) = extended_power(index)+fit%arch(i)*shock_history(lag,i)
            else
               gamma_i = 0.0_dp
               if (i <= fit%o) gamma_i = fit%leverage(i)
               moment = innovation_asym_power_moment(fit%distribution,fit%delta,gamma_i,fit%shape,fit%skew)
               extended_power(index) = extended_power(index)+fit%arch(i)*moment*extended_power(lag)
            end if
         end do
         do i = 1, fit%q
            lag = index-i
            extended_power(index) = extended_power(index)+fit%garch(i)*extended_power(lag)
         end do
         extended_power(index) = max(extended_power(index),1.0e-14_dp)
      end do
      mean_forecast = fit%mean
      variance_forecast = extended_power(n+1:total)**(2.0_dp/fit%delta)
      if (present(power_forecast)) power_forecast = extended_power(n+1:total)
   end subroutine forecast_univariate

   subroutine forecast_garch11(fit, n_ahead, mean_forecast, variance_forecast)
      type(garch11_fit), intent(in) :: fit
      integer, intent(in) :: n_ahead
      real(dp), intent(out) :: mean_forecast(n_ahead), variance_forecast(n_ahead)
      call forecast_univariate(fit,n_ahead,mean_forecast,variance_forecast)
   end subroutine forecast_garch11

   subroutine simulate_aparch(n, mean, omega, arch, leverage, garch, delta, distribution, shape, skew, y, variance, &
      burnin, power_scale)
      integer, intent(in) :: n
      real(dp), intent(in) :: mean, omega, arch(:), leverage(:), garch(:), delta, shape, skew
      character(len=*), intent(in) :: distribution
      real(dp), intent(out) :: y(n), variance(n)
      integer, intent(in), optional :: burnin
      real(dp), intent(out), optional :: power_scale(n)
      real(dp), allocatable :: yall(:), power_all(:), eps(:)
      real(dp) :: initial_power, persistence, gamma_i, moment
      integer :: i, j, lag, b, total
      logical :: ok
      persistence = sum(garch)
      do j = 1, size(arch)
         gamma_i = 0.0_dp
         if (j <= size(leverage)) gamma_i = leverage(j)
         moment = innovation_asym_power_moment(distribution,delta,gamma_i,shape,skew)
         persistence = persistence+arch(j)*moment
      end do
      ok = omega > 0.0_dp .and. all(arch >= 0.0_dp) .and. all(garch >= 0.0_dp) .and. &
         all(abs(leverage) < 1.0_dp) .and. delta > 0.0_dp .and. persistence < 0.9995_dp .and. &
         distribution_is_valid(distribution,shape,skew)
      if (.not. ok) error stop 'simulate_aparch: invalid or nonstationary parameters'
      b = 500
      if (present(burnin)) b = max(0,burnin)
      total = n+b
      allocate(yall(total),power_all(total),eps(total))
      initial_power = omega/max(1.0_dp-persistence,1.0e-6_dp)
      do i = 1, total
         power_all(i) = omega
         do j = 1, size(arch)
            lag = i-j
            gamma_i = 0.0_dp
            if (j <= size(leverage)) gamma_i = leverage(j)
            if (lag >= 1) then
               power_all(i) = power_all(i)+arch(j)*(abs(eps(lag))-gamma_i*eps(lag))**delta
            else
               power_all(i) = power_all(i)+arch(j)*moment_for_lag(distribution,delta,gamma_i,shape,skew)*initial_power
            end if
         end do
         do j = 1, size(garch)
            lag = i-j
            if (lag >= 1) then
               power_all(i) = power_all(i)+garch(j)*power_all(lag)
            else
               power_all(i) = power_all(i)+garch(j)*initial_power
            end if
         end do
         power_all(i) = max(power_all(i),1.0e-14_dp)
         eps(i) = power_all(i)**(1.0_dp/delta)*random_innovation(distribution,shape,skew)
         yall(i) = mean+eps(i)
      end do
      y = yall(b+1:total)
      variance = power_all(b+1:total)**(2.0_dp/delta)
      if (present(power_scale)) power_scale = power_all(b+1:total)
   end subroutine simulate_aparch

   subroutine simulate_garchpq(n, mean, omega, arch, garch, distribution, shape, skew, y, variance, burnin)
      integer, intent(in) :: n
      real(dp), intent(in) :: mean, omega, arch(:), garch(:), shape, skew
      character(len=*), intent(in) :: distribution
      real(dp), intent(out) :: y(n), variance(n)
      integer, intent(in), optional :: burnin
      call simulate_aparch(n,mean,omega,arch,[real(dp)::],garch,2.0_dp,distribution,shape,skew,y,variance,burnin)
   end subroutine simulate_garchpq

   subroutine simulate_garch11(n, mean, omega, alpha, beta, y, variance, burnin)
      integer, intent(in) :: n
      real(dp), intent(in) :: mean, omega, alpha, beta
      real(dp), intent(out) :: y(n), variance(n)
      integer, intent(in), optional :: burnin
      call simulate_garchpq(n,mean,omega,[alpha],[beta],'norm',8.0_dp,1.0_dp,y,variance,burnin)
   end subroutine simulate_garch11

   function garch_objective(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value, mean, omega, delta, shape, skew, ll
      real(dp), allocatable :: arch(:), leverage(:), garch(:)
      real(dp), allocatable :: residuals(:), power_scale(:), variance(:), standardized(:)
      logical :: ok
      select type (context => generic_context)
      type is (garch_context)
         allocate(arch(context%spec%p),leverage(context%spec%o),garch(context%spec%q))
         call decode_parameters(x,context%spec,mean,omega,arch,leverage,garch,delta,shape,skew,ok)
         if (.not. ok) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         allocate(residuals(size(context%y)),power_scale(size(context%y)),variance(size(context%y)), &
            standardized(size(context%y)))
         call filter_aparch(context%y,mean,omega,arch,leverage,garch,delta,context%spec%distribution,shape,skew, &
            residuals,power_scale,variance,standardized,ll,ok)
         if (ok) then
            value = -ll
         else
            value = huge(1.0_dp)/100.0_dp
         end if
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function garch_objective

   pure integer function parameter_count(spec) result(number)
      type(univariate_spec), intent(in) :: spec
      character(len=:), allocatable :: model, dist
      number = 1+spec%p+spec%q
      if (spec%include_mean) number = number+1
      model = trim(adjustl(lower_string(spec%model)))
      dist = trim(adjustl(lower_string(spec%distribution)))
      if (model == 'aparch') then
         number = number+spec%o
         if (spec%fit_delta) number = number+1
      end if
      if ((dist == 'std' .or. dist == 'sstd' .or. dist == 'ged' .or. dist == 'sged') .and. spec%fit_shape) &
         number = number+1
      if ((dist == 'snorm' .or. dist == 'sstd' .or. dist == 'sged') .and. spec%fit_skew) number = number+1
   end function parameter_count

   subroutine normalize_spec(spec)
      type(univariate_spec), intent(inout) :: spec
      spec%model = trim(adjustl(lower_string(spec%model)))
      spec%distribution = trim(adjustl(lower_string(spec%distribution)))
      if (trim(spec%model) == 'garch') then
         spec%o = 0
         spec%delta = 2.0_dp
         spec%fit_delta = .false.
      end if
      if (trim(spec%distribution) == 'norm') then
         spec%fit_shape = .false.
         spec%fit_skew = .false.
         spec%shape = 8.0_dp
         spec%skew = 1.0_dp
      else if (trim(spec%distribution) == 'snorm') then
         spec%fit_shape = .false.
         spec%shape = 8.0_dp
      else if (trim(spec%distribution) == 'std' .or. trim(spec%distribution) == 'ged') then
         spec%fit_skew = .false.
         spec%skew = 1.0_dp
      end if
   end subroutine normalize_spec

   subroutine initial_parameters(y, spec, x)
      real(dp), intent(in) :: y(:)
      type(univariate_spec), intent(in) :: spec
      real(dp), intent(out) :: x(:)
      real(dp) :: variance0, total_target, residual_weight, each_arch, each_garch
      integer :: index, i
      variance0 = max(sum((y-sum(y)/real(size(y),dp))**2)/real(size(y),dp),1.0e-8_dp)
      index = 0
      if (spec%include_mean) then
         index = index+1
         x(index) = sum(y)/real(size(y),dp)
      end if
      index = index+1
      x(index) = log(max(0.05_dp*variance0**(0.5_dp*spec%delta),1.0e-12_dp))
      total_target = 0.93_dp
      residual_weight = 0.995_dp-total_target
      each_arch = 0.08_dp/real(spec%p,dp)
      if (spec%q > 0) then
         each_garch = 0.85_dp/real(spec%q,dp)
      else
         each_arch = total_target/real(spec%p,dp)
         each_garch = 0.0_dp
      end if
      do i = 1, spec%p
         index = index+1
         x(index) = log(max(each_arch/residual_weight,1.0e-8_dp))
      end do
      do i = 1, spec%q
         index = index+1
         x(index) = log(max(each_garch/residual_weight,1.0e-8_dp))
      end do
      do i = 1, spec%o
         index = index+1
         x(index) = 0.0_dp
      end do
      if (trim(spec%model) == 'aparch' .and. spec%fit_delta) then
         index = index+1
         x(index) = inverse_logistic_bound(spec%delta,0.25_dp,4.0_dp)
      end if
      if (shape_is_used(spec%distribution) .and. spec%fit_shape) then
         index = index+1
         if (index_is_student(spec%distribution)) then
            x(index) = inverse_logistic_bound(spec%shape,2.05_dp,50.0_dp)
         else
            x(index) = inverse_logistic_bound(spec%shape,0.25_dp,6.0_dp)
         end if
      end if
      if (skew_is_used(spec%distribution) .and. spec%fit_skew) then
         index = index+1
         x(index) = inverse_logistic_bound(spec%skew,0.2_dp,5.0_dp)
      end if
   end subroutine initial_parameters

   subroutine decode_parameters(x, spec, mean, omega, arch, leverage, garch, delta, shape, skew, ok)
      real(dp), intent(in) :: x(:)
      type(univariate_spec), intent(in) :: spec
      real(dp), intent(out) :: mean, omega, arch(spec%p), leverage(spec%o), garch(spec%q), delta, shape, skew
      logical, intent(out) :: ok
      real(dp), allocatable :: e(:)
      real(dp) :: denominator
      integer :: index, i
      ok = size(x) == parameter_count(spec)
      if (.not. ok) return
      index = 0
      if (spec%include_mean) then
         index = index+1
         mean = x(index)
      else
         mean = 0.0_dp
      end if
      index = index+1
      omega = exp(max(-40.0_dp,min(20.0_dp,x(index))))
      allocate(e(spec%p+spec%q))
      do i = 1, size(e)
         index = index+1
         e(i) = exp(max(-30.0_dp,min(30.0_dp,x(index))))
      end do
      denominator = 1.0_dp+sum(e)
      arch = 0.995_dp*e(1:spec%p)/denominator
      if (spec%q > 0) garch = 0.995_dp*e(spec%p+1:spec%p+spec%q)/denominator
      do i = 1, spec%o
         index = index+1
         leverage(i) = 0.995_dp*tanh(x(index))
      end do
      delta = spec%delta
      if (trim(spec%model) == 'aparch' .and. spec%fit_delta) then
         index = index+1
         delta = logistic_bound(x(index),0.25_dp,4.0_dp)
      end if
      shape = spec%shape
      if (shape_is_used(spec%distribution) .and. spec%fit_shape) then
         index = index+1
         if (index_is_student(spec%distribution)) then
            shape = logistic_bound(x(index),2.05_dp,50.0_dp)
         else
            shape = logistic_bound(x(index),0.25_dp,6.0_dp)
         end if
      end if
      skew = spec%skew
      if (skew_is_used(spec%distribution) .and. spec%fit_skew) then
         index = index+1
         skew = logistic_bound(x(index),0.2_dp,5.0_dp)
      end if
      ok = omega > 0.0_dp .and. delta > 0.0_dp .and. distribution_is_valid(spec%distribution,shape,skew) .and. &
         sum(arch)+sum(garch) < 0.9995_dp
   end subroutine decode_parameters

   pure function logistic_bound(raw, lower, upper) result(value)
      real(dp), intent(in) :: raw, lower, upper
      real(dp) :: value, clipped
      clipped = max(-30.0_dp,min(30.0_dp,raw))
      value = lower+(upper-lower)/(1.0_dp+exp(-clipped))
   end function logistic_bound

   pure function inverse_logistic_bound(value, lower, upper) result(raw)
      real(dp), intent(in) :: value, lower, upper
      real(dp) :: raw, ratio
      ratio = (min(upper-1.0e-8_dp,max(lower+1.0e-8_dp,value))-lower)/(upper-lower)
      raw = log(ratio/(1.0_dp-ratio))
   end function inverse_logistic_bound

   pure function shape_is_used(distribution) result(used)
      character(len=*), intent(in) :: distribution
      logical :: used
      character(len=:), allocatable :: dist
      dist = trim(adjustl(lower_string(distribution)))
      used = dist == 'std' .or. dist == 'sstd' .or. dist == 'ged' .or. dist == 'sged'
   end function shape_is_used

   pure function skew_is_used(distribution) result(used)
      character(len=*), intent(in) :: distribution
      logical :: used
      character(len=:), allocatable :: dist
      dist = trim(adjustl(lower_string(distribution)))
      used = dist == 'snorm' .or. dist == 'sstd' .or. dist == 'sged'
   end function skew_is_used

   pure function index_is_student(distribution) result(used)
      character(len=*), intent(in) :: distribution
      logical :: used
      character(len=:), allocatable :: dist
      dist = trim(adjustl(lower_string(distribution)))
      used = dist == 'std' .or. dist == 'sstd'
   end function index_is_student

   function moment_for_lag(distribution, delta, gamma, shape, skew) result(moment)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in) :: delta, gamma, shape, skew
      real(dp) :: moment
      moment = innovation_asym_power_moment(distribution,delta,gamma,shape,skew)
   end function moment_for_lag

end module gogarch_univariate
