! SPDX-License-Identifier: GPL-3.0-only
module ufrisk_varcast
   use kind_mod, only : dp
   use ufrisk_types
   use ufrisk_math, only : normal_quantile, normal_pdf, student_quantile, student_pdf, &
      sample_standard_deviation, finite_vector
   use ufrisk_smoothing, only : long_memory_smooth
   use ufrisk_loggarch, only : arfilt_coefficients, log_variance_forecasts, estimate_student_df
   use smoots_estimation, only : msmooth
   use smoots_arma, only : fit_arma
   use smoots_status, only : sm_ok, sm_iteration_limit
   use fracdiff_model_api, only : fracdiff_fit
   use fracdiff_status, only : fd_ok, fd_iteration_limit
   use rugarch_mod, only : rugarch_spec_t, rugarch_filter_t, rugarch_spec, rugarch_fit, &
      rugarch_filter, rugarch_model_sgarch, rugarch_model_egarch, rugarch_model_aparch, &
      rugarch_model_figarch, rugarch_distribution_normal, rugarch_distribution_student
   implicit none
   private
   public :: varcast
contains
   function varcast(prices, options) result(out)
      real(dp), intent(in) :: prices(:)
      type(ufrisk_options), intent(in), optional :: options
      type(ufrisk_result) :: out
      type(ufrisk_options) :: controls
      type(rugarch_spec_t) :: specification
      type(rugarch_filter_t) :: combined_filter
      real(dp), allocatable :: returns(:), centered(:), zeta_in(:), zeta_out(:)
      real(dp), allocatable :: residual_in(:), residual_out(:), combined(:)
      real(dp), allocatable :: coefficients(:), log_forecast(:), conditional_in(:)
      real(dp), allocatable :: standardized(:)
      real(dp) :: mean_return, scale_forecast, mule, mulz, correction
      real(dp) :: q_es, q_var, standard_deviation, es_factor, density_es
      real(dp) :: backcast
      integer :: n_returns, n_in, n_out, status, rugarch_model, rugarch_distribution
      logical :: is_garch_model

      scale_forecast = 1.0_dp
      mule = 0.0_dp
      mulz = 0.0_dp
      correction = 1.0_dp
      controls = ufrisk_options()
      if (present(options)) controls = options
      out%model = controls%model
      out%distribution = controls%distribution
      out%smooth = controls%smooth
      out%arch_order = controls%arch_order
      out%garch_order = controls%garch_order
      out%var_tail_probability = 1.0_dp-controls%var_confidence
      out%es_tail_probability = 1.0_dp-controls%es_confidence
      out%model_name = ufrisk_model_label(controls%model,controls%smooth==ufrisk_smooth_lpr)

      if (.not.valid_controls(prices,controls)) then
         call fail(out,ufrisk_invalid_input)
         return
      end if
      n_returns = size(prices)-1
      n_out = controls%n_out
      n_in = n_returns-n_out
      allocate(returns(n_returns),centered(n_returns))
      returns = log(prices(2:)/prices(:size(prices)-1))
      mean_return = sum(returns(:n_in))/real(n_in,dp)
      centered = returns-mean_return
      out%mean_return = mean_return
      allocate(out%returns_in(n_in),out%returns_out(n_out),out%centered_in(n_in),out%centered_out(n_out))
      out%returns_in = returns(:n_in)
      out%returns_out = returns(n_in+1:)
      out%centered_in = centered(:n_in)
      out%centered_out = centered(n_in+1:)

      call estimate_scale(centered(:n_in),centered(n_in+1:),controls,out,residual_in, &
         residual_out,zeta_in,zeta_out,mule,scale_forecast,status)
      if (status /= ufrisk_ok) then
         call fail(out,status)
         return
      end if
      out%scale_forecast = scale_forecast

      is_garch_model = controls%model == ufrisk_model_sgarch .or. &
         controls%model == ufrisk_model_egarch .or. controls%model == ufrisk_model_aparch .or. &
         controls%model == ufrisk_model_figarch
      if (is_garch_model) then
         select case (controls%model)
         case (ufrisk_model_sgarch)
            rugarch_model = rugarch_model_sgarch
         case (ufrisk_model_egarch)
            rugarch_model = rugarch_model_egarch
         case (ufrisk_model_aparch)
            rugarch_model = rugarch_model_aparch
         case default
            rugarch_model = rugarch_model_figarch
         end select
         rugarch_distribution = merge(rugarch_distribution_student,rugarch_distribution_normal, &
            controls%distribution==ufrisk_distribution_student)
         specification = rugarch_spec(variance_model=rugarch_model,distribution=rugarch_distribution, &
            arch_order=controls%arch_order,garch_order=controls%garch_order,include_mean=.false., &
            truncation_lag=min(controls%truncation_lag,n_in))
         out%rugarch_fit = rugarch_fit(zeta_in,specification, &
            max_iterations=controls%max_fit_iterations,tolerance=controls%fit_tolerance)
         if (out%rugarch_fit%info /= 0) then
            call fail(out,ufrisk_model_fit_failed)
            return
         end if
         allocate(combined(n_returns)); combined = [zeta_in,zeta_out]
         backcast = sum((zeta_in-sum(zeta_in)/real(n_in,dp))**2)/real(max(1,n_in-1),dp)
         combined_filter = rugarch_filter(combined,specification,out%rugarch_fit%parameters, &
            backcast_value=backcast)
         if (combined_filter%info /= 0) then
            call fail(out,ufrisk_model_fit_failed)
            return
         end if
         allocate(out%sigma_in(n_in),out%sigma_forecast(n_out))
         out%sigma_in = combined_filter%conditional_sigma(:n_in)*out%scale
         out%sigma_forecast = combined_filter%conditional_sigma(n_in+1:)*scale_forecast
         if (controls%distribution == ufrisk_distribution_student) then
            out%degrees_freedom = max(2.0001_dp,min(100.0_dp,out%rugarch_fit%parameters%shape))
         else
            out%degrees_freedom = 0.0_dp
         end if
      else if (controls%model == ufrisk_model_lgarch) then
         if (controls%arch_order < controls%garch_order) then
            call fail(out,ufrisk_invalid_input)
            out%message = 'lGARCH requires arch_order >= garch_order in its ARMA representation'
            return
         end if
         call fit_arma(residual_in,controls%arch_order,controls%garch_order,.false., &
            out%arma_fit,status,tolerance=controls%fit_tolerance, &
            max_iterations=controls%max_fit_iterations)
         if (status /= sm_ok .and. status /= sm_iteration_limit) then
            call fail(out,ufrisk_model_fit_failed)
            return
         end if
         coefficients = arfilt_coefficients(out%arma_fit%ar,-out%arma_fit%ma,0.0_dp, &
            controls%log_filter_lag)
         mulz = -log(max(sum(exp(min(residual_in,700.0_dp)))/real(n_in,dp),tiny(1.0_dp)))
         call log_variance_forecasts(residual_in,residual_out,coefficients, &
            merge(mule,0.0_dp,controls%smooth==ufrisk_smooth_none),log_forecast)
         allocate(conditional_in(n_in),out%sigma_in(n_in),out%sigma_forecast(n_out))
         conditional_in = exp(0.5_dp*(residual_in-merge(mule,0.0_dp, &
            controls%smooth==ufrisk_smooth_none)-out%arma_fit%residuals+mule-mulz))
         out%sigma_in = conditional_in*out%scale
         out%sigma_forecast = exp(0.5_dp*(log_forecast+mule-mulz))*scale_forecast
         if (controls%distribution == ufrisk_distribution_student) then
            allocate(standardized(n_in)); standardized = centered(:n_in)/out%sigma_in
            out%degrees_freedom = estimate_student_df(standardized)
         else
            out%degrees_freedom = 0.0_dp
         end if
      else
         if (controls%arch_order < controls%garch_order) then
            call fail(out,ufrisk_invalid_input)
            out%message = 'filGARCH requires arch_order >= garch_order in its ARFIMA representation'
            return
         end if
         out%fracdiff_fit = fracdiff_fit(residual_in,nar=controls%arch_order, &
            nma=controls%garch_order,m_terms=controls%fractional_terms)
         if (out%fracdiff_fit%status /= fd_ok .and. &
            out%fracdiff_fit%status /= fd_iteration_limit) then
            call fail(out,ufrisk_model_fit_failed)
            return
         end if
         coefficients = arfilt_coefficients(out%fracdiff_fit%ar,out%fracdiff_fit%ma, &
            out%fracdiff_fit%d,controls%log_filter_lag)
         mule = sum(residual_in)/real(n_in,dp)
         call log_variance_forecasts(residual_in,residual_out,coefficients,mule,log_forecast)
         allocate(conditional_in(n_in),out%sigma_in(n_in),out%sigma_forecast(n_out))
         conditional_in = exp(0.5_dp*out%fracdiff_fit%fitted)
         correction = sample_standard_deviation(zeta_in/max(conditional_in,sqrt(tiny(1.0_dp))))
         if (correction <= sqrt(epsilon(1.0_dp))) correction = 1.0_dp
         out%sigma_in = conditional_in*out%scale*correction
         out%sigma_forecast = exp(0.5_dp*(log_forecast+mule))*scale_forecast*correction
         if (controls%distribution == ufrisk_distribution_student) then
            allocate(standardized(n_in)); standardized = centered(:n_in)/out%sigma_in
            out%degrees_freedom = estimate_student_df(standardized)
         else
            out%degrees_freedom = 0.0_dp
         end if
      end if

      if (.not.finite_vector(out%sigma_in) .or. .not.finite_vector(out%sigma_forecast) .or. &
         any(out%sigma_in <= 0.0_dp) .or. any(out%sigma_forecast <= 0.0_dp)) then
         call fail(out,ufrisk_numerical_failure)
         return
      end if
      allocate(out%var_es_level(n_out),out%var_var_level(n_out),out%expected_shortfall(n_out))
      if (controls%distribution == ufrisk_distribution_normal) then
         standard_deviation = 1.0_dp
         q_es = normal_quantile(controls%es_confidence)
         q_var = normal_quantile(controls%var_confidence)
         es_factor = normal_pdf(q_es)/(1.0_dp-controls%es_confidence)
      else
         standard_deviation = sqrt(out%degrees_freedom/(out%degrees_freedom-2.0_dp))
         q_es = student_quantile(controls%es_confidence,out%degrees_freedom)
         q_var = student_quantile(controls%var_confidence,out%degrees_freedom)
         density_es = student_pdf(q_es,out%degrees_freedom)
         es_factor = density_es/(1.0_dp-controls%es_confidence)* &
            (out%degrees_freedom+q_es*q_es)/(out%degrees_freedom-1.0_dp)/standard_deviation
      end if
      out%var_es_level = -mean_return+out%sigma_forecast*q_es/standard_deviation
      out%var_var_level = -mean_return+out%sigma_forecast*q_var/standard_deviation
      out%expected_shortfall = -mean_return+out%sigma_forecast*es_factor
      out%status = ufrisk_ok
      out%message = ufrisk_status_message(out%status)
   end function varcast

   subroutine estimate_scale(centered_in,centered_out,controls,out,residual_in,residual_out, &
      zeta_in,zeta_out,mule,scale_forecast,status)
      real(dp), intent(in) :: centered_in(:),centered_out(:)
      type(ufrisk_options), intent(in) :: controls
      type(ufrisk_result), intent(inout) :: out
      real(dp), allocatable, intent(out) :: residual_in(:),residual_out(:),zeta_in(:),zeta_out(:)
      real(dp), intent(out) :: mule,scale_forecast
      integer, intent(out) :: status
      real(dp), allocatable :: log_squared(:)
      real(dp) :: scale_correction
      integer :: n_in
      mule = 0.0_dp
      scale_forecast = 1.0_dp
      status = ufrisk_invalid_input
      n_in = size(centered_in)
      allocate(log_squared(n_in),residual_out(size(centered_out)),zeta_in(n_in),zeta_out(size(centered_out)))
      log_squared = log(max(centered_in*centered_in,tiny(1.0_dp)))
      select case (controls%smooth)
      case (ufrisk_smooth_none)
         allocate(out%scale(n_in),residual_in(n_in)); out%scale = 1.0_dp
         residual_in = log_squared
         residual_out = log(max(centered_out*centered_out,tiny(1.0_dp)))
         mule = sum(residual_in)/real(n_in,dp)
         scale_forecast = 1.0_dp
      case (ufrisk_smooth_lpr)
         if (controls%model == ufrisk_model_figarch .or. &
            controls%model == ufrisk_model_filgarch) then
            call long_memory_smooth(log_squared,out%long_memory_smooth, &
               p=controls%smoothing_order,mu=controls%smoothing_mu, &
               b_start=controls%smoothing_start,iterations=controls%smoothing_iterations, &
               p_min=controls%fractional_p_min,p_max=controls%fractional_p_max, &
               q_min=controls%fractional_q_min,q_max=controls%fractional_q_max, &
               m_terms=controls%fractional_terms)
            if (out%long_memory_smooth%status /= ufrisk_ok) then
               status = ufrisk_smoothing_failed
               return
            end if
            residual_in = out%long_memory_smooth%residuals
            allocate(out%scale(n_in))
            scale_correction = sample_standard_deviation(centered_in/ &
               exp(0.5_dp*out%long_memory_smooth%estimate))
            if (scale_correction <= sqrt(epsilon(1.0_dp))) scale_correction = 1.0_dp
            out%scale = exp(0.5_dp*out%long_memory_smooth%estimate)*scale_correction
         else
            call msmooth(log_squared,out%short_memory_smooth,p=controls%smoothing_order, &
               mu=controls%smoothing_mu,b_start=controls%smoothing_start, &
               algorithm=controls%smoothing_algorithm)
            if (out%short_memory_smooth%status /= sm_ok .and. &
               out%short_memory_smooth%status /= sm_iteration_limit) then
               status = ufrisk_smoothing_failed
               return
            end if
            residual_in = out%short_memory_smooth%residuals
            mule = -log(max(sum(exp(min(residual_in,700.0_dp)))/real(n_in,dp),tiny(1.0_dp)))
            allocate(out%scale(n_in))
            out%scale = exp(0.5_dp*(out%short_memory_smooth%estimate-mule))
         end if
         scale_forecast = out%scale(n_in)
         if (controls%model == ufrisk_model_figarch .or. &
            controls%model == ufrisk_model_filgarch) then
            residual_out = log(max(centered_out*centered_out,tiny(1.0_dp))) - &
               out%long_memory_smooth%estimate(n_in)
         else
            residual_out = log(max(centered_out*centered_out,tiny(1.0_dp))) - &
               out%short_memory_smooth%estimate(n_in)
         end if
      case default
         status = ufrisk_invalid_input
         return
      end select
      zeta_in = centered_in/out%scale
      zeta_out = centered_out/scale_forecast
      if (controls%smooth == ufrisk_smooth_lpr .and. &
         (controls%model == ufrisk_model_figarch .or. controls%model == ufrisk_model_filgarch)) &
         mule = sum(residual_in)/real(n_in,dp)
      status = ufrisk_ok
   end subroutine estimate_scale

   pure logical function valid_controls(prices,controls) result(valid)
      real(dp), intent(in) :: prices(:)
      type(ufrisk_options), intent(in) :: controls
      valid = size(prices) >= 25 .and. all(prices > 0.0_dp) .and. finite_vector(prices) .and. &
         controls%n_out >= 1 .and. controls%n_out <= size(prices)-21 .and. &
         controls%model >= ufrisk_model_sgarch .and. controls%model <= ufrisk_model_filgarch .and. &
         controls%distribution >= ufrisk_distribution_normal .and. &
         controls%distribution <= ufrisk_distribution_student .and. &
         controls%smooth >= ufrisk_smooth_none .and. controls%smooth <= ufrisk_smooth_lpr .and. &
         controls%arch_order >= 0 .and. controls%garch_order >= 0 .and. &
         controls%arch_order+controls%garch_order > 0 .and. &
         controls%var_confidence > 0.0_dp .and. controls%var_confidence < 1.0_dp .and. &
         controls%es_confidence > 0.0_dp .and. controls%es_confidence < 1.0_dp .and. &
         controls%smoothing_iterations >= 1 .and. controls%max_fit_iterations >= 1
      if (controls%model == ufrisk_model_figarch) valid = valid .and. &
         controls%arch_order >= 1 .and. controls%garch_order >= 1
   end function valid_controls

   subroutine fail(out,status)
      type(ufrisk_result), intent(inout) :: out
      integer, intent(in) :: status
      out%status = status
      out%message = ufrisk_status_message(status)
   end subroutine fail
end module ufrisk_varcast
