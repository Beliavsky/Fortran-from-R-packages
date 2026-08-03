! SPDX-License-Identifier: GPL-3.0-only
module garchito_models
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchito_kinds, only : dp
   use garchito_types, only : garchito_control, garchito_result, &
      garchito_success, garchito_max_iterations, garchito_invalid_input, &
      garchito_numerical_failure
   use garchito_utils, only : mean_value, median_value, ols_line, all_finite
   use garchito_optimizer, only : bounded_nelder_mead
   implicit none
   private

   integer, parameter :: model_unified = 1
   integer, parameter :: model_realized_nj = 2
   integer, parameter :: model_realized_jump = 3
   integer, parameter :: model_option_nj_hom = 4
   integer, parameter :: model_option_nj_het = 5
   integer, parameter :: model_option_jump_hom = 6
   integer, parameter :: model_option_jump_het = 7
   real(dp), parameter :: variance_floor = 1.0e-14_dp
   real(dp), parameter :: stationarity_margin = 1.0e-8_dp
   real(dp), parameter :: bad_objective = huge(1.0_dp) / 32.0_dp

   type :: model_context
      integer :: kind = 0
      real(dp), allocatable :: rv(:)
      real(dp), allocatable :: jv(:)
      real(dp), allocatable :: nv(:)
      real(dp), allocatable :: returns(:)
      real(dp) :: median_jv = 0.0_dp
   end type model_context

   public :: unified_est, realized_est, realized_est_option

contains

   subroutine unified_est(rv, returns, result, control)
      real(dp), intent(in) :: rv(:), returns(:)
      type(garchito_result), intent(out) :: result
      type(garchito_control), intent(in), optional :: control
      type(model_context) :: context
      type(garchito_control) :: settings
      real(dp) :: start(3), lower(3), upper(3)

      call initialize_result(result)
      if (size(rv) /= size(returns) .or. size(rv) < 2) then
         result%message = 'rv and returns must have equal length of at least two'
         return
      end if
      if (.not. all_finite(rv) .or. .not. all_finite(returns) .or. any(rv < 0.0_dp)) then
         result%message = 'rv must be finite and nonnegative; returns must be finite'
         return
      end if

      settings = garchito_control()
      if (present(control)) settings = control
      context%kind = model_unified
      context%rv = rv
      context%returns = returns
      start = [max(mean_value(rv), 1.0e-8_dp), 0.4_dp, 0.3_dp]
      lower = [1.0e-12_dp, 0.0_dp, 0.0_dp]
      upper = [1.0_dp, 1.0_dp, 1.0_dp]
      call fit_context(context, start, lower, upper, settings, result)
      call assign_names(result, ['omega_g         ', 'beta_g          ', 'gamma           '])
      if (allocated(result%coefficients)) call build_sigma(context, result)
   end subroutine unified_est

   subroutine realized_est(rv, result, jv, control)
      real(dp), intent(in) :: rv(:)
      type(garchito_result), intent(out) :: result
      real(dp), intent(in), optional :: jv(:)
      type(garchito_control), intent(in), optional :: control
      type(model_context) :: context
      type(garchito_control) :: settings
      real(dp), allocatable :: start(:), lower(:), upper(:)

      call initialize_result(result)
      if (size(rv) < 2 .or. .not. all_finite(rv) .or. any(rv < 0.0_dp)) then
         result%message = 'rv must contain at least two finite nonnegative values'
         return
      end if
      settings = garchito_control()
      if (present(control)) settings = control
      context%rv = rv

      if (present(jv)) then
         if (size(jv) /= size(rv) .or. .not. all_finite(jv) .or. any(jv < 0.0_dp)) then
            result%message = 'jv must be finite, nonnegative, and the same length as rv'
            return
         end if
         context%kind = model_realized_jump
         context%jv = jv
         context%median_jv = median_value(jv)
         allocate(start(4), lower(4), upper(4))
         start = [max(mean_value(rv), 1.0e-8_dp), 0.45_dp, 0.5_dp, 0.3_dp]
         lower = [1.0e-12_dp, 0.0_dp, 0.0_dp, 0.0_dp]
         upper = [1.0_dp, 1.0_dp, 10.0_dp, 1.0_dp]
         call fit_context(context, start, lower, upper, settings, result)
         call assign_names(result, ['omega_g         ', 'alpha_g         ', &
                                    'beta_g          ', 'gamma           '])
      else
         context%kind = model_realized_nj
         allocate(start(3), lower(3), upper(3))
         start = [max(mean_value(rv), 1.0e-8_dp), 0.45_dp, 0.3_dp]
         lower = [1.0e-12_dp, 0.0_dp, 0.0_dp]
         upper = [1.0_dp, 1.0_dp, 1.0_dp]
         call fit_context(context, start, lower, upper, settings, result)
         call assign_names(result, ['omega_g         ', 'alpha_g         ', 'gamma           '])
      end if
      if (allocated(result%coefficients)) call build_sigma(context, result)
   end subroutine realized_est

   subroutine realized_est_option(rv, nv, result, jv, homogeneous, control)
      real(dp), intent(in) :: rv(:), nv(:)
      type(garchito_result), intent(out) :: result
      real(dp), intent(in), optional :: jv(:)
      logical, intent(in), optional :: homogeneous
      type(garchito_control), intent(in), optional :: control
      type(model_context) :: context
      type(garchito_control) :: settings
      real(dp), allocatable :: start(:), lower(:), upper(:)
      real(dp) :: slope, intercept, sigma_e
      logical :: hom

      call initialize_result(result)
      if (size(rv) /= size(nv) .or. size(rv) < 3) then
         result%message = 'rv and nv must have equal length of at least three'
         return
      end if
      if (.not. all_finite(rv) .or. .not. all_finite(nv) .or. any(rv < 0.0_dp)) then
         result%message = 'rv must be finite and nonnegative; nv must be finite'
         return
      end if
      hom = .true.
      if (present(homogeneous)) hom = homogeneous
      settings = garchito_control()
      if (present(control)) settings = control
      call ols_line(rv, nv, slope, intercept, sigma_e)
      context%rv = rv
      context%nv = nv

      if (present(jv)) then
         if (size(jv) /= size(rv) .or. .not. all_finite(jv) .or. any(jv < 0.0_dp)) then
            result%message = 'jv must be finite, nonnegative, and the same length as rv'
            return
         end if
         context%jv = jv
         context%median_jv = median_value(jv)
         if (hom) then
            context%kind = model_option_jump_hom
            allocate(start(7), lower(7), upper(7))
            start = [max(mean_value(rv), 1.0e-8_dp), 0.4_dp, 0.5_dp, 0.3_dp, &
                     slope, intercept, sigma_e]
            lower = [1.0e-12_dp, 0.0_dp, 0.0_dp, 0.0_dp, -10.0_dp, -10.0_dp, 1.0e-10_dp]
            upper = [1.0_dp, 1.0_dp, 10.0_dp, 1.0_dp, 10.0_dp, 10.0_dp, 10.0_dp]
            call fit_context(context, start, lower, upper, settings, result)
            call assign_names(result, ['omega_g         ', 'alpha_g         ', &
               'beta_g          ', 'gamma           ', 'a               ', &
               'b               ', 'sigma_e         '])
         else
            context%kind = model_option_jump_het
            allocate(start(8), lower(8), upper(8))
            start = [max(mean_value(rv), 1.0e-8_dp), 0.4_dp, 0.5_dp, 0.3_dp, &
                     slope, intercept, sigma_e, sqrt(0.5_dp)]
            lower = [1.0e-12_dp, 0.0_dp, 0.0_dp, 0.0_dp, -10.0_dp, -10.0_dp, &
                     1.0e-10_dp, 0.0_dp]
            upper = [1.0_dp, 1.0_dp, 10.0_dp, 1.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp]
            call fit_context(context, start, lower, upper, settings, result)
            call assign_names(result, ['omega_g         ', 'alpha_g         ', &
               'beta_g          ', 'gamma           ', 'a               ', &
               'b               ', 'sigma_e         ', 'zeta            '])
         end if
      else
         if (hom) then
            context%kind = model_option_nj_hom
            allocate(start(6), lower(6), upper(6))
            start = [max(mean_value(rv), 1.0e-8_dp), 0.4_dp, 0.3_dp, &
                     slope, intercept, sigma_e]
            lower = [1.0e-12_dp, 0.0_dp, 0.0_dp, -10.0_dp, -10.0_dp, 1.0e-10_dp]
            upper = [1.0_dp, 1.0_dp, 1.0_dp, 10.0_dp, 10.0_dp, 10.0_dp]
            call fit_context(context, start, lower, upper, settings, result)
            call assign_names(result, ['omega_g         ', 'alpha_g         ', &
               'gamma           ', 'a               ', 'b               ', 'sigma_e         '])
         else
            context%kind = model_option_nj_het
            allocate(start(7), lower(7), upper(7))
            start = [max(mean_value(rv), 1.0e-8_dp), 0.4_dp, 0.3_dp, &
                     slope, intercept, sigma_e, sqrt(0.5_dp)]
            lower = [1.0e-12_dp, 0.0_dp, 0.0_dp, -10.0_dp, -10.0_dp, 1.0e-10_dp, 0.0_dp]
            upper = [1.0_dp, 1.0_dp, 1.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp]
            call fit_context(context, start, lower, upper, settings, result)
            call assign_names(result, ['omega_g         ', 'alpha_g         ', &
               'gamma           ', 'a               ', 'b               ', &
               'sigma_e         ', 'zeta            '])
         end if
      end if
      if (allocated(result%coefficients)) call build_sigma(context, result)
   end subroutine realized_est_option

   subroutine fit_context(context, start, lower, upper, control, result)
      type(model_context), intent(in) :: context
      real(dp), intent(in) :: start(:), lower(:), upper(:)
      type(garchito_control), intent(in) :: control
      type(garchito_result), intent(inout) :: result
      real(dp), allocatable :: candidate(:), best(:), alt(:)
      real(dp) :: value, best_value
      integer :: status, iter, eval, total_eval, attempt, best_status, best_iter

      allocate(candidate(size(start)), best(size(start)), alt(size(start)))
      best_value = bad_objective
      total_eval = 0
      best_status = garchito_numerical_failure
      best_iter = 0

      do attempt = 1, 3
         alt = start
         select case (attempt)
         case (2)
            call set_dynamic_start(context%kind, alt, 0.20_dp, 0.65_dp)
         case (3)
            call set_dynamic_start(context%kind, alt, 0.65_dp, 0.20_dp)
         end select
         call bounded_nelder_mead(model_objective, project_parameters, context, &
            alt, lower, upper, control, candidate, value, status, iter, eval)
         total_eval = total_eval + eval
         if (value < best_value) then
            best_value = value
            best = candidate
            best_status = status
            best_iter = iter
         end if
      end do

      allocate(result%coefficients(size(start)))
      result%coefficients = best
      result%objective = best_value
      result%convergence = best_status
      result%iterations = best_iter
      result%evaluations = total_eval
      if (.not. ieee_is_finite(best_value) .or. best_value >= bad_objective / 2.0_dp) then
         result%convergence = garchito_numerical_failure
         result%message = 'optimizer failed to find a finite likelihood'
      else if (best_status == garchito_success) then
         result%message = 'success'
      else if (best_status == garchito_max_iterations) then
         result%message = 'maximum iterations reached; best finite estimate returned'
      else
         result%message = 'numerical optimization failure'
      end if
   end subroutine fit_context

   subroutine set_dynamic_start(kind, x, first_weight, gamma_weight)
      integer, intent(in) :: kind
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: first_weight, gamma_weight

      x(2) = first_weight
      select case (kind)
      case (model_unified, model_realized_nj, model_option_nj_hom, model_option_nj_het)
         x(3) = gamma_weight
      case default
         x(4) = gamma_weight
      end select
   end subroutine set_dynamic_start

   subroutine project_parameters(x, data)
      real(dp), intent(inout) :: x(:)
      class(*), intent(in) :: data
      real(dp) :: total, factor
      integer :: gamma_index

      select type (context => data)
      type is (model_context)
         select case (context%kind)
         case (model_unified, model_realized_nj, model_option_nj_hom, model_option_nj_het)
            gamma_index = 3
         case default
            gamma_index = 4
         end select
         total = x(2) + x(gamma_index)
         if (total > 1.0_dp - stationarity_margin) then
            factor = (1.0_dp - stationarity_margin) / total
            x(2) = x(2) * factor
            x(gamma_index) = x(gamma_index) * factor
         end if
      class default
         error stop 'invalid GARCHIto optimizer context'
      end select
   end subroutine project_parameters

   subroutine model_objective(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      class(*), intent(in) :: data
      real(dp) :: h, log_variance, residual, denom
      integer :: i, n

      value = bad_objective
      select type (context => data)
      type is (model_context)
         n = size(context%rv)
         call initial_variance(context, x, h, denom)
         if (denom <= stationarity_margin .or. h <= variance_floor .or. &
             .not. ieee_is_finite(h)) return
         value = 0.0_dp
         do i = 1, n
            if (i > 1) h = next_variance(context, x, h, i - 1)
            if (h <= variance_floor .or. .not. ieee_is_finite(h)) then
               value = bad_objective
               return
            end if
            value = value + log(h) + context%rv(i) / h
            select case (context%kind)
            case (model_option_nj_hom)
               residual = context%nv(i) - x(5) - x(4)*h
               log_variance = 2.0_dp*log(max(x(6), 1.0e-10_dp))
            case (model_option_nj_het)
               residual = context%nv(i) - x(5) - x(4)*h
               log_variance = 2.0_dp*log(max(x(6), 1.0e-10_dp)) + x(7)*log(h)
            case (model_option_jump_hom)
               residual = context%nv(i) - x(6) - x(5)*h
               log_variance = 2.0_dp*log(max(x(7), 1.0e-10_dp))
            case (model_option_jump_het)
               residual = context%nv(i) - x(6) - x(5)*h
               log_variance = 2.0_dp*log(max(x(7), 1.0e-10_dp)) + x(8)*log(h)
            case default
               cycle
            end select
            if (log_variance > log(huge(1.0_dp)) - 4.0_dp .or. &
                log_variance < log(tiny(1.0_dp)) + 4.0_dp) then
               value = bad_objective
               return
            end if
            value = value + log_variance + residual*residual / exp(log_variance)
         end do
         if (.not. ieee_is_finite(value)) value = bad_objective
      class default
         value = bad_objective
      end select
   end subroutine model_objective

   subroutine initial_variance(context, x, h, denominator)
      type(model_context), intent(in) :: context
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: h, denominator

      select case (context%kind)
      case (model_unified, model_realized_nj, model_option_nj_hom, model_option_nj_het)
         denominator = 1.0_dp - x(2) - x(3)
         h = x(1) / denominator
      case default
         denominator = 1.0_dp - x(2) - x(4)
         h = (x(1) + x(3)*context%median_jv) / denominator
      end select
   end subroutine initial_variance

   real(dp) function next_variance(context, x, previous_h, previous_index) result(h)
      type(model_context), intent(in) :: context
      real(dp), intent(in) :: x(:), previous_h
      integer, intent(in) :: previous_index

      select case (context%kind)
      case (model_unified)
         h = x(1) + x(3)*previous_h + x(2)*context%returns(previous_index)**2
      case (model_realized_nj, model_option_nj_hom, model_option_nj_het)
         h = x(1) + x(3)*previous_h + x(2)*context%rv(previous_index)
      case default
         h = x(1) + x(4)*previous_h + x(2)*context%rv(previous_index) + &
             x(3)*context%jv(previous_index)
      end select
   end function next_variance

   subroutine build_sigma(context, result)
      type(model_context), intent(in) :: context
      type(garchito_result), intent(inout) :: result
      real(dp) :: h, denominator
      integer :: i, n

      n = size(context%rv)
      allocate(result%sigma(n))
      call initial_variance(context, result%coefficients, h, denominator)
      if (denominator <= stationarity_margin .or. h <= variance_floor) then
         result%sigma = 0.0_dp
         result%pred = 0.0_dp
         result%convergence = garchito_numerical_failure
         result%message = 'invalid fitted unconditional variance'
         return
      end if
      result%sigma(1) = h
      do i = 2, n
         h = next_variance(context, result%coefficients, h, i - 1)
         result%sigma(i) = h
      end do
      result%pred = next_variance(context, result%coefficients, h, n)
   end subroutine build_sigma

   subroutine assign_names(result, names)
      type(garchito_result), intent(inout) :: result
      character(len=16), intent(in) :: names(:)
      allocate(result%coefficient_names(size(names)))
      result%coefficient_names = names
   end subroutine assign_names

   subroutine initialize_result(result)
      type(garchito_result), intent(out) :: result
      result%pred = 0.0_dp
      result%objective = huge(1.0_dp)
      result%convergence = garchito_invalid_input
      result%iterations = 0
      result%evaluations = 0
      result%message = 'invalid input'
   end subroutine initialize_result

end module garchito_models
