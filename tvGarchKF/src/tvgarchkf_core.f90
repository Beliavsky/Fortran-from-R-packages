! Modern Fortran translation of the computational core of tvGarchKF.
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-or-later
module tvgarchkf_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
   use fgarch_kinds, only : dp
   use fgarch_optimizer, only : optimizer_result, nelder_mead
   use fgarch_rng, only : seed_rng, random_normal
   use fgarch_fit, only : fit_garch11
   use fgarch_types, only : garch_fit_result
   use tvgarchkf_types
   use tvgarchkf_functions, only : evaluate_tv_function
   implicit none
   private

   real(dp), parameter :: penalty_value = 1.0e8_dp

   type :: fit_context
      real(dp), allocatable :: series(:)
      type(tvgarch_spec) :: template
      integer :: predict = 0
      logical :: corrected_constraints = .false.
   end type fit_context

   public :: tvgarch_kalman_filter, tvgarch_kalman_loglike
   public :: tvgarch_kalman_fit, tvgarch_kalman_print
   public :: tvgarch_simulate, tv_parameter

contains

   function tvgarch_kalman_filter(series, spec, predict, corrected_constraints) result(output)
      real(dp), intent(in) :: series(:)
      type(tvgarch_spec), intent(in) :: spec
      integer, intent(in), optional :: predict
      logical, intent(in), optional :: corrected_constraints
      type(tvgarch_filter_result) :: output
      real(dp), allocatable :: extended(:), centered_square(:)
      real(dp) :: persistence, unconditional, term
      integer :: n, n_predict, i, s1, s2, s3
      character(len=160) :: m1, m2, m3
      logical :: corrected, missing

      n_predict = 0
      if (present(predict)) n_predict = max(0,predict)
      corrected = .false.
      if (present(corrected_constraints)) corrected = corrected_constraints
      n = size(series)+n_predict
      if (size(series) < 1 .or. n < 1) then
         output%status = 10
         output%message = 'series must contain at least one observation'
         return
      end if

      allocate(extended(n),centered_square(n))
      extended(1:size(series)) = series
      if (n_predict > 0) extended(size(series)+1:n) = ieee_nan()

      output%omega = evaluate_tv_function(spec%omega,n,s1,m1)
      output%alpha = evaluate_tv_function(spec%alpha,n,s2,m2)
      output%beta = evaluate_tv_function(spec%beta,n,s3,m3)
      if (s1 /= 0 .or. s2 /= 0 .or. s3 /= 0) then
         output%status = 11
         output%message = trim(m1)//'; '//trim(m2)//'; '//trim(m3)
         return
      end if

      do i = 1, n
         persistence = output%alpha(i)+output%beta(i)
         if (corrected) then
            if (output%omega(i) <= 0.0_dp .or. output%alpha(i) < 0.0_dp .or. &
                output%beta(i) < 0.0_dp .or. persistence >= 1.0_dp) then
               output%status = 12
               output%message = 'corrected constraints require omega > 0, alpha >= 0, beta >= 0, alpha + beta < 1'
               output%criterion = penalty_value
               return
            end if
         else
            if (persistence >= 1.0_dp) then
               output%status = 13
               output%message = 'alpha + beta must be less than one'
               output%criterion = penalty_value
               return
            end if
         end if
         if (abs(1.0_dp-persistence) <= epsilon(1.0_dp)) then
            output%status = 14
            output%message = 'unconditional variance denominator is zero'
            output%criterion = penalty_value
            return
         end if
      end do

      allocate(output%state(n),output%state_variance(n),output%mse(n),output%gain(n), &
               output%conditional_variance(n),output%sigma(n))
      output%state = 0.0_dp
      output%state_variance = 0.0_dp
      output%mse = 0.0_dp
      output%gain = 0.0_dp

      persistence = output%alpha(1)+output%beta(1)
      if (1.0_dp-persistence*persistence <= 0.0_dp) then
         output%status = 15
         output%message = 'initial state-variance denominator is nonpositive'
         output%criterion = penalty_value
         return
      end if
      output%state_variance(1) = output%alpha(1)**2/(1.0_dp-persistence**2)

      do i = 1, n
         unconditional = output%omega(i)/(1.0_dp-output%alpha(i)-output%beta(i))
         missing = ieee_is_nan(extended(i))
         if (.not. missing .and. .not. ieee_is_finite(extended(i))) then
            output%status = 16
            output%message = 'series contains an infinite observation'
            output%criterion = penalty_value
            return
         end if
         if (missing) then
            centered_square(i) = 0.0_dp
         else
            centered_square(i) = extended(i)**2-unconditional
         end if

         output%mse(i) = output%state_variance(i)+1.0_dp
         if (output%mse(i) <= 0.0_dp .or. .not. ieee_is_finite(output%mse(i))) then
            output%status = 17
            output%message = 'Kalman prediction-error variance is nonpositive'
            output%criterion = penalty_value
            return
         end if
         persistence = output%alpha(i)+output%beta(i)
         output%gain(i) = (persistence*output%state_variance(i)+output%alpha(i))/output%mse(i)

         if (i < n) then
            if (missing) then
               output%state_variance(i+1) = persistence**2*output%state_variance(i)+output%alpha(i)**2
               output%state(i+1) = persistence*output%state(i)
            else
               output%state_variance(i+1) = persistence*output%state_variance(i)*(persistence-output%gain(i)) + &
                                             output%alpha(i)*(output%alpha(i)-output%gain(i))
               output%state(i+1) = (persistence-output%gain(i))*output%state(i) + &
                                    output%gain(i)*centered_square(i)
            end if
            unconditional = output%omega(i)/(1.0_dp-output%alpha(i)-output%beta(i))
            if (output%state(i+1)+unconditional < 0.0_dp) output%state(i+1) = 0.0_dp
         end if
      end do

      output%conditional_variance = output%state + output%omega/(1.0_dp-output%alpha-output%beta)
      if (any(output%conditional_variance <= 0.0_dp) .or. &
          any(.not. ieee_is_finite(output%conditional_variance))) then
         output%status = 18
         output%message = 'conditional variance is nonpositive or non-finite'
         output%criterion = penalty_value
         return
      end if
      output%sigma = sqrt(output%conditional_variance)

      output%criterion = 0.0_dp
      do i = 1, n
         if (.not. ieee_is_nan(extended(i))) then
            term = sqrt(output%mse(i))*extended(i)**2/output%conditional_variance(i) + &
                   log(output%conditional_variance(i))-0.5_dp*log(output%mse(i))
            output%criterion = output%criterion+term
         end if
      end do
      if (.not. ieee_is_finite(output%criterion)) then
         output%status = 19
         output%message = 'criterion is non-finite'
         output%criterion = penalty_value
         return
      end if
      output%status = 0
      output%message = 'ok'
   end function tvgarch_kalman_filter

   function tvgarch_kalman_loglike(parameters, series, template, predict, corrected_constraints) result(value)
      real(dp), intent(in) :: parameters(:), series(:)
      type(tvgarch_spec), intent(in) :: template
      integer, intent(in), optional :: predict
      logical, intent(in), optional :: corrected_constraints
      real(dp) :: value
      type(tvgarch_spec) :: spec
      type(tvgarch_filter_result) :: filtered
      integer :: n_predict
      logical :: corrected

      n_predict = 0
      if (present(predict)) n_predict = predict
      corrected = .false.
      if (present(corrected_constraints)) corrected = corrected_constraints
      spec = template
      if (.not. unpack_parameters(parameters,spec)) then
         value = penalty_value
         return
      end if
      filtered = tvgarch_kalman_filter(series,spec,n_predict,corrected)
      value = filtered%criterion
   end function tvgarch_kalman_loglike

   function tvgarch_kalman_fit(series, initial_spec, predict, corrected_constraints, &
                               max_iterations, tolerance) result(fit)
      real(dp), intent(in) :: series(:)
      type(tvgarch_spec), intent(in) :: initial_spec
      integer, intent(in), optional :: predict, max_iterations
      logical, intent(in), optional :: corrected_constraints
      real(dp), intent(in), optional :: tolerance
      type(tvgarch_fit_result) :: fit
      type(fit_context) :: context
      type(optimizer_result) :: opt
      real(dp), allocatable :: x0(:)
      real(dp) :: tol
      integer :: maxit, n_predict
      logical :: corrected

      n_predict = 0
      if (present(predict)) n_predict = max(0,predict)
      corrected = .false.
      if (present(corrected_constraints)) corrected = corrected_constraints
      maxit = 3000
      if (present(max_iterations)) maxit = max(1,max_iterations)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))

      x0 = pack_parameters(initial_spec)
      if (size(x0) == 0) then
         fit%status = 20
         fit%message = 'initial specification has empty coefficient arrays'
         return
      end if
      context%series = series
      context%template = initial_spec
      context%predict = n_predict
      context%corrected_constraints = corrected

      opt = nelder_mead(fit_objective,x0,context,step=0.08_dp,tolerance=tol,max_iterations=maxit)
      fit%parameters = opt%x
      allocate(fit%rounded_parameters(size(opt%x)))
      fit%rounded_parameters = anint(opt%x*1.0e5_dp)/1.0e5_dp
      fit%spec = initial_spec
      if (.not. unpack_parameters(opt%x,fit%spec)) then
         fit%status = 21
         fit%message = 'optimizer returned an invalid parameter vector'
         return
      end if
      fit%filter = tvgarch_kalman_filter(series,fit%spec,n_predict,corrected)
      fit%criterion = fit%filter%criterion
      fit%iterations = opt%iterations
      fit%evaluations = opt%evaluations
      if (fit%filter%status /= 0) then
         fit%status = fit%filter%status
         fit%message = fit%filter%message
      else if (opt%status == 0) then
         fit%status = 0
         fit%message = 'converged'
      else
         fit%status = 1
         fit%message = 'iteration limit reached; inspect estimates'
      end if
   end function tvgarch_kalman_fit

   function tvgarch_kalman_print(parameters, series, template, predict, corrected_constraints) result(output)
      real(dp), intent(in) :: parameters(:), series(:)
      type(tvgarch_spec), intent(in) :: template
      integer, intent(in), optional :: predict
      logical, intent(in), optional :: corrected_constraints
      type(tvgarch_filter_result) :: output
      type(tvgarch_spec) :: spec
      integer :: n_predict
      logical :: corrected

      n_predict = 0
      if (present(predict)) n_predict = max(0,predict)
      corrected = .false.
      if (present(corrected_constraints)) corrected = corrected_constraints
      spec = template
      if (.not. unpack_parameters(parameters,spec)) then
         output%status = 22
         output%message = 'parameter-vector length does not match the specification'
         return
      end if
      output = tvgarch_kalman_filter(series,spec,n_predict,corrected)
   end function tvgarch_kalman_print

   function tvgarch_simulate(n, spec, innovations, seed, corrected_constraints) result(simulation)
      integer, intent(in) :: n
      type(tvgarch_spec), intent(in) :: spec
      real(dp), intent(in), optional :: innovations(:)
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: corrected_constraints
      type(tvgarch_simulation_result) :: simulation
      real(dp), allocatable :: z(:), omega_values(:), alpha_values(:), beta_values(:)
      real(dp) :: persistence
      integer :: i, s1, s2, s3
      character(len=160) :: m1, m2, m3
      logical :: corrected

      corrected = .false.
      if (present(corrected_constraints)) corrected = corrected_constraints
      if (n < 1) then
         simulation%status = 30
         simulation%message = 'n must be positive'
         return
      end if
      omega_values = evaluate_tv_function(spec%omega,n,s1,m1)
      alpha_values = evaluate_tv_function(spec%alpha,n,s2,m2)
      beta_values = evaluate_tv_function(spec%beta,n,s3,m3)
      call move_alloc(omega_values,simulation%omega)
      call move_alloc(alpha_values,simulation%alpha)
      call move_alloc(beta_values,simulation%beta)
      if (s1 /= 0 .or. s2 /= 0 .or. s3 /= 0) then
         simulation%status = 31
         simulation%message = trim(m1)//'; '//trim(m2)//'; '//trim(m3)
         return
      end if
      if (corrected) then
         do i = 1, n
            persistence = simulation%alpha(i)+simulation%beta(i)
            if (simulation%omega(i) <= 0.0_dp .or. simulation%alpha(i) < 0.0_dp .or. &
                simulation%beta(i) < 0.0_dp .or. persistence >= 1.0_dp) then
               simulation%status = 32
               simulation%message = 'invalid variance parameters under corrected constraints'
               return
            end if
         end do
      end if

      allocate(z(n),simulation%returns(n),simulation%variance(n))
      if (present(innovations)) then
         if (size(innovations) < n) then
            simulation%status = 33
            simulation%message = 'innovation vector is shorter than n'
            return
         end if
         z = innovations(1:n)
      else
         if (present(seed)) call seed_rng(seed)
         do i = 1, n
            z(i) = random_normal()
         end do
      end if
      simulation%returns = 0.0_dp
      simulation%variance = 0.0_dp
      do i = 2, n
         simulation%variance(i) = simulation%omega(i)+simulation%alpha(i)*simulation%returns(i-1)**2 + &
                                  simulation%beta(i)*simulation%variance(i-1)
         if (simulation%variance(i) < 0.0_dp .or. .not. ieee_is_finite(simulation%variance(i))) then
            simulation%status = 34
            simulation%message = 'simulation produced a negative or non-finite variance'
            return
         end if
         simulation%returns(i) = z(i)*sqrt(simulation%variance(i))
      end do
      simulation%status = 0
      simulation%message = 'ok'
   end function tvgarch_simulate

   function tv_parameter(data, shift, window, max_iterations) result(result)
      real(dp), intent(in) :: data(:)
      integer, intent(in) :: shift, window
      integer, intent(in), optional :: max_iterations
      type(tv_parameter_result) :: result
      type(garch_fit_result) :: fit
      integer :: m, j, first, last, maxit

      maxit = 1200
      if (present(max_iterations)) maxit = max(1,max_iterations)
      if (shift <= 0 .or. window <= 2 .or. size(data) < window) then
         result%status = 40
         result%message = 'require shift > 0, window > 2, and data length >= window'
         return
      end if
      if (any(.not. ieee_is_finite(data))) then
         result%status = 41
         result%message = 'tv_parameter does not accept missing or infinite data'
         return
      end if
      m = int(real(size(data)-window,dp)/real(shift,dp)+1.0_dp)
      allocate(result%midpoint(m),result%omega(m),result%alpha(m),result%beta(m))

      fit = fit_garch11(data,fit_mean=.false.,max_iterations=maxit)
      if (.not. allocated(fit%spec%alpha) .or. .not. allocated(fit%spec%beta)) then
         result%status = 42
         result%message = 'global GARCH fit failed'
         return
      end if
      result%global_omega = fit%spec%omega
      result%global_alpha = fit%spec%alpha(1)
      result%global_beta = fit%spec%beta(1)

      do j = 1, m
         first = shift*(j-1)+1
         last = first+window-1
         fit = fit_garch11(data(first:last),fit_mean=.false.,max_iterations=maxit)
         if (.not. allocated(fit%spec%alpha) .or. .not. allocated(fit%spec%beta)) then
            result%status = 43
            write(result%message,'(a,i0)') 'local GARCH fit failed for block ',j
            return
         end if
         result%midpoint(j) = real(shift*(j-1),dp)+0.5_dp*real(window,dp)
         result%omega(j) = fit%spec%omega
         result%alpha(j) = fit%spec%alpha(1)
         result%beta(j) = fit%spec%beta(1)
      end do
      result%status = 0
      result%message = 'ok'
   end function tv_parameter

   function fit_objective(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value

      select type (context => generic_context)
      type is (fit_context)
         value = tvgarch_kalman_loglike(x,context%series,context%template,context%predict, &
                                        context%corrected_constraints)
      class default
         value = penalty_value
      end select
   end function fit_objective

   function pack_parameters(spec) result(parameters)
      type(tvgarch_spec), intent(in) :: spec
      real(dp), allocatable :: parameters(:)
      integer :: n1, n2, n3

      n1 = 0; n2 = 0; n3 = 0
      if (allocated(spec%omega%coefficients)) n1 = size(spec%omega%coefficients)
      if (allocated(spec%alpha%coefficients)) n2 = size(spec%alpha%coefficients)
      if (allocated(spec%beta%coefficients)) n3 = size(spec%beta%coefficients)
      allocate(parameters(n1+n2+n3))
      if (n1 > 0) parameters(1:n1) = spec%omega%coefficients
      if (n2 > 0) parameters(n1+1:n1+n2) = spec%alpha%coefficients
      if (n3 > 0) parameters(n1+n2+1:n1+n2+n3) = spec%beta%coefficients
   end function pack_parameters

   logical function unpack_parameters(parameters, spec) result(ok)
      real(dp), intent(in) :: parameters(:)
      type(tvgarch_spec), intent(inout) :: spec
      integer :: n1, n2, n3

      ok = .false.
      if (.not. allocated(spec%omega%coefficients) .or. .not. allocated(spec%alpha%coefficients) .or. &
          .not. allocated(spec%beta%coefficients)) return
      n1 = size(spec%omega%coefficients)
      n2 = size(spec%alpha%coefficients)
      n3 = size(spec%beta%coefficients)
      if (size(parameters) /= n1+n2+n3) return
      spec%omega%coefficients = parameters(1:n1)
      spec%alpha%coefficients = parameters(n1+1:n1+n2)
      spec%beta%coefficients = parameters(n1+n2+1:n1+n2+n3)
      ok = .true.
   end function unpack_parameters

   function ieee_nan() result(value)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      real(dp) :: value
      value = ieee_value(0.0_dp,ieee_quiet_nan)
   end function ieee_nan

end module tvgarchkf_core
