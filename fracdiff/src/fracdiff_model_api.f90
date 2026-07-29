! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_model_api
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fracdiff_kinds, only : dp
   use fracdiff_status, only : fd_ok, fd_invalid_input, fd_iteration_limit, &
                               fracdiff_status_message
   use fracdiff_types, only : fracdiff_model, fracdiff_summary
   use fracdiff_optimize, only : profile_evaluation, brent_profile_d
   use fracdiff_inference, only : numerical_likelihood_hessian, covariance_from_hessian, &
                                  inverse_normal_cdf
   use fracdiff_difference, only : diffseries
   use fracdiff_filter, only : conditional_arma_series_residuals
   implicit none
   private

   public :: fracdiff_fit, fracdiff_var
   public :: fracdiff_coefficients, fracdiff_confint, summarize_fracdiff
   public :: fracdiff_aic, fracdiff_bic

contains

   function fracdiff_fit(x, nar, nma, ar_initial, ma_initial, d_tolerance, d_range, &
                         h, m_terms, trace) result(model)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: nar, nma, m_terms, trace
      real(dp), intent(in), optional :: ar_initial(:), ma_initial(:)
      real(dp), intent(in), optional :: d_tolerance, d_range(2), h
      type(fracdiff_model) :: model

      real(dp), allocatable :: ar(:), ma(:), differenced(:)
      real(dp) :: dtol, drange_local(2), d_optimum, h_local
      integer :: p, q, m_value, trace_value, outer_iterations
      integer :: inference_status, residual_status
      type(profile_evaluation) :: evaluation

      p = 0
      q = 0
      if (present(nar)) p = nar
      if (present(nma)) q = nma
      m_value = 100
      if (present(m_terms)) m_value = m_terms
      trace_value = 0
      if (present(trace)) trace_value = trace
      dtol = epsilon(1.0_dp)**0.25_dp
      if (present(d_tolerance)) then
         if (d_tolerance > 0.0_dp) dtol = min(0.1_dp, max(d_tolerance, epsilon(1.0_dp)**0.5_dp))
      end if
      drange_local = [0.0_dp, 0.5_dp]
      if (present(d_range)) drange_local = d_range

      model%n = size(x)
      model%p = p
      model%q = q
      model%m_terms = m_value
      model%d_tol = dtol
      if (size(x) < max(3,p+q+3) .or. p < 0 .or. q < 0 .or. m_value < 1 .or. &
          drange_local(1) < -0.5_dp .or. drange_local(2) > 0.5_dp .or. &
          drange_local(1) >= drange_local(2)) then
         model%status = fd_invalid_input
         model%message = fracdiff_status_message(model%status)
         call allocate_empty_model(model, p, q, size(x))
         return
      end if

      allocate(ar(p), ma(q))
      ar = 0.0_dp
      ma = 0.0_dp
      if (present(ar_initial)) ar(1:min(p,size(ar_initial))) = ar_initial(1:min(p,size(ar_initial)))
      if (present(ma_initial)) ma(1:min(q,size(ma_initial))) = ma_initial(1:min(q,size(ma_initial)))

      call brent_profile_d(x, m_value, 0.0_dp, drange_local, dtol, ar, ma, d_optimum, &
                           evaluation, outer_iterations, trace_value > 0)
      model%d = d_optimum
      model%mean = evaluation%estimated_mean
      model%log_likelihood = evaluation%log_likelihood
      model%fnorm_min = evaluation%residual_norm
      model%sigma = sqrt(evaluation%white_noise_variance)
      model%outer_iterations = outer_iterations
      model%function_evaluations = evaluation%function_evaluations
      model%gradient_evaluations = evaluation%gradient_evaluations
      allocate(model%ar(p), model%ma(q))
      model%ar = evaluation%ar
      model%ma = evaluation%ma
      model%status = evaluation%status
      if (model%status == fd_iteration_limit) then
         model%message = "fit completed but an optimizer reached its iteration limit"
      else
         model%message = fracdiff_status_message(model%status)
      end if

      h_local = min(0.1_dp, sqrt(epsilon(1.0_dp)/2.0_dp)*(1.0_dp + abs(model%log_likelihood)))
      if (present(h)) then
         if (h > 0.0_dp) h_local = h
      end if
      model%h = h_local
      call allocate_inference(model)
      call numerical_likelihood_hessian(x, m_value, model%d, model%ar, model%ma, h_local, &
                                        model%hessian, inference_status)
      if (inference_status == fd_ok) then
         call covariance_from_hessian(model%hessian, model%covariance, model%std_error, &
                                      model%correlation, inference_status)
      end if
      if (inference_status /= fd_ok .and. model%status == fd_ok) then
         model%status = inference_status
         model%message = "fit succeeded; covariance calculation failed: "// &
                         fracdiff_status_message(inference_status)
      end if

      allocate(differenced(size(x)), model%residuals(size(x)), model%fitted(size(x)))
      call diffseries(x, model%d, differenced, residual_status)
      if (residual_status == fd_ok) then
         call conditional_arma_series_residuals(differenced, model%ar, model%ma, &
                                                model%residuals, model%fitted, residual_status)
         model%fitted = x - model%residuals
      else
         model%residuals = ieee_value(1.0_dp, ieee_quiet_nan)
         model%fitted = model%residuals
      end if
   end function fracdiff_fit

   subroutine fracdiff_var(x, model, h, status)
      real(dp), intent(in) :: x(:)
      type(fracdiff_model), intent(inout) :: model
      real(dp), intent(in) :: h
      integer, intent(out), optional :: status
      integer :: local_status

      if (h <= 0.0_dp .or. size(x) /= model%n) then
         local_status = fd_invalid_input
         if (present(status)) status = local_status
         return
      end if
      model%h = h
      if (.not. allocated(model%hessian)) call allocate_inference(model)
      call numerical_likelihood_hessian(x, model%m_terms, model%d, model%ar, model%ma, h, &
                                        model%hessian, local_status)
      if (local_status == fd_ok) then
         call covariance_from_hessian(model%hessian, model%covariance, model%std_error, &
                                      model%correlation, local_status)
      end if
      if (local_status /= fd_ok) then
         model%status = local_status
         model%message = "covariance calculation failed: "//fracdiff_status_message(local_status)
      end if
      if (present(status)) status = local_status
   end subroutine fracdiff_var

   function fracdiff_coefficients(model) result(coefficients)
      type(fracdiff_model), intent(in) :: model
      real(dp), allocatable :: coefficients(:)
      integer :: p, q

      p = size(model%ar)
      q = size(model%ma)
      allocate(coefficients(1+p+q))
      coefficients(1) = model%d
      if (p > 0) coefficients(2:p+1) = model%ar
      if (q > 0) coefficients(p+2:) = model%ma
   end function fracdiff_coefficients

   function fracdiff_confint(model, level) result(intervals)
      type(fracdiff_model), intent(in) :: model
      real(dp), intent(in), optional :: level
      real(dp), allocatable :: intervals(:,:)
      real(dp), allocatable :: coefficients(:)
      real(dp) :: confidence_level, alpha, z
      integer :: npar

      confidence_level = 0.95_dp
      if (present(level)) confidence_level = level
      coefficients = fracdiff_coefficients(model)
      npar = size(coefficients)
      allocate(intervals(npar,2))
      if (confidence_level <= 0.0_dp .or. confidence_level >= 1.0_dp .or. &
          .not. allocated(model%std_error)) then
         intervals = ieee_value(1.0_dp, ieee_quiet_nan)
         return
      end if
      alpha = 0.5_dp*(1.0_dp-confidence_level)
      z = inverse_normal_cdf(1.0_dp-alpha)
      intervals(:,1) = coefficients - z*model%std_error
      intervals(:,2) = coefficients + z*model%std_error
   end function fracdiff_confint

   function summarize_fracdiff(model) result(summary)
      type(fracdiff_model), intent(in) :: model
      type(fracdiff_summary) :: summary
      real(dp), allocatable :: coefficients(:)
      integer :: npar

      coefficients = fracdiff_coefficients(model)
      npar = size(coefficients)
      summary%degrees_freedom = npar + 1
      summary%aic = fracdiff_aic(model)
      summary%bic = fracdiff_bic(model)
      allocate(summary%coefficients(npar,4))
      summary%coefficients(:,1) = coefficients
      if (allocated(model%std_error)) then
         summary%coefficients(:,2) = model%std_error
         where (model%std_error > 0.0_dp)
            summary%coefficients(:,3) = coefficients/model%std_error
         elsewhere
            summary%coefficients(:,3) = ieee_value(1.0_dp, ieee_quiet_nan)
         end where
         summary%coefficients(:,4) = 2.0_dp*normal_lower_tail(-abs(summary%coefficients(:,3)))
      else
         summary%coefficients(:,2:4) = ieee_value(1.0_dp, ieee_quiet_nan)
      end if
   end function summarize_fracdiff

   pure function fracdiff_aic(model) result(value)
      type(fracdiff_model), intent(in) :: model
      real(dp) :: value
      integer :: number_parameters

      number_parameters = 1 + size(model%ar) + size(model%ma) + 1
      value = -2.0_dp*model%log_likelihood + 2.0_dp*real(number_parameters,dp)
   end function fracdiff_aic

   pure function fracdiff_bic(model) result(value)
      type(fracdiff_model), intent(in) :: model
      real(dp) :: value
      integer :: number_parameters

      number_parameters = 1 + size(model%ar) + size(model%ma) + 1
      value = -2.0_dp*model%log_likelihood + log(real(model%n,dp))*real(number_parameters,dp)
   end function fracdiff_bic

   pure elemental function normal_lower_tail(x) result(probability)
      real(dp), intent(in) :: x
      real(dp) :: probability

      probability = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_lower_tail

   subroutine allocate_inference(model)
      type(fracdiff_model), intent(inout) :: model
      integer :: npar

      npar = 1 + size(model%ar) + size(model%ma)
      if (allocated(model%covariance)) deallocate(model%covariance)
      if (allocated(model%correlation)) deallocate(model%correlation)
      if (allocated(model%std_error)) deallocate(model%std_error)
      if (allocated(model%hessian)) deallocate(model%hessian)
      allocate(model%covariance(npar,npar), model%correlation(npar,npar), &
               model%std_error(npar), model%hessian(npar,npar))
      model%covariance = 0.0_dp
      model%correlation = 0.0_dp
      model%std_error = 0.0_dp
      model%hessian = 0.0_dp
   end subroutine allocate_inference

   subroutine allocate_empty_model(model, p, q, n)
      type(fracdiff_model), intent(inout) :: model
      integer, intent(in) :: p, q, n

      if (.not. allocated(model%ar)) allocate(model%ar(max(0,p)))
      if (.not. allocated(model%ma)) allocate(model%ma(max(0,q)))
      if (.not. allocated(model%residuals)) allocate(model%residuals(max(0,n)))
      if (.not. allocated(model%fitted)) allocate(model%fitted(max(0,n)))
      model%ar = 0.0_dp
      model%ma = 0.0_dp
      model%residuals = 0.0_dp
      model%fitted = 0.0_dp
   end subroutine allocate_empty_model

end module fracdiff_model_api
