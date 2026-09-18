! SPDX-License-Identifier: MIT
! SPDX-FileComment: Pure cubic smoothing-spline implementation corresponding to R stats.
module r_stats_smoothing_spline
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_linalg, only: inverse_matrix, solve_system
   use r_optional, only: optval
   use r_ordering, only: r_order
   use r_quantiles, only: r_quantile_type7
   use r_status, only: r_invalid_input, r_ok
   use r_stats_types, only: smooth_spline_fit_t
   implicit none
   private

   public :: predict_smooth_spline, smooth_spline

contains

   pure function smooth_spline(x, y, weights, spar, degrees_freedom, lambda, all_knots, &
                               n_knots, tolerance, ordinary_cv) result(fit)
      !! Fits a natural cubic smoothing spline by penalized least squares.
      real(dp), intent(in) :: x(:) !! Finite predictor values with at least four distinct entries.
      real(dp), intent(in) :: y(:) !! Finite responses conformable with `x`.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      real(dp), intent(in), optional :: spar !! R-compatible smoothing parameter.
      real(dp), intent(in), optional :: degrees_freedom !! Target effective degrees of freedom.
      real(dp), intent(in), optional :: lambda !! Nonnegative roughness multiplier.
      logical, intent(in), optional :: all_knots !! Use every distinct predictor as a knot.
      integer, intent(in), optional :: n_knots !! Explicit number of predictor knots.
      real(dp), intent(in), optional :: tolerance !! Positive predictor-grouping tolerance.
      logical, intent(in), optional :: ordinary_cv !! Select smoothing by leave-one-out CV.
      type(smooth_spline_fit_t) :: fit
      real(dp), allocatable :: basis(:, :), design_crossproduct(:, :), penalty(:, :)
      real(dp), allocatable :: group_codes(:), sorted_weights(:), sorted_x(:), sorted_y(:)
      integer, allocatable :: group_index(:), order(:)
      real(dp) :: adjusted_leverage, effective_count, lower_spar
      real(dp) :: target_df, upper_spar, within_group_ss
      integer :: group, i, selection_count

      selection_count = merge(1, 0, present(spar)) + merge(1, 0, present(degrees_freedom)) + &
                        merge(1, 0, present(lambda))
      fit%ordinary_cv = optval(ordinary_cv, .false.)
      if (size(x) /= size(y) .or. size(x) < 4 .or. selection_count > 1 .or. &
          .not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
         fit%status = r_invalid_input
         return
      end if
      call r_order(x, order)
      sorted_x = x(order)
      sorted_y = y(order)
      if (present(tolerance)) then
         fit%tolerance = tolerance
      else
         fit%tolerance = 1.0e-6_dp*(r_quantile_type7(x, 0.75_dp) - &
                                    r_quantile_type7(x, 0.25_dp))
      end if
      if (.not. ieee_is_finite(fit%tolerance) .or. fit%tolerance <= 0.0_dp) then
         fit%status = r_invalid_input
         return
      end if
      group_codes = round_to_even((sorted_x - sum(x)/real(size(x), dp))/fit%tolerance)
      if (present(weights)) then
         if (size(weights) /= size(x) .or. any(weights < 0.0_dp) .or. &
             .not. all(ieee_is_finite(weights)) .or. sum(weights) <= 0.0_dp) then
            fit%status = r_invalid_input
            return
         end if
         sorted_weights = weights(order)
         sorted_weights = sorted_weights*real(count(sorted_weights > 0.0_dp), dp)/ &
                          sum(sorted_weights)
      else
         allocate (sorted_weights(size(x)), source=1.0_dp)
      end if
      fit%n_original = size(x)
      call aggregate_duplicate_predictors(sorted_x, sorted_y, sorted_weights, group_codes, &
                                           fit%x, fit%y, fit%weights, group_index, &
                                           within_group_ss)
      if (size(fit%x) < 4) then
         fit%status = r_invalid_input
         return
      end if
      if (present(all_knots) .and. present(n_knots)) then
         fit%status = r_invalid_input
         return
      end if
      if (present(n_knots)) then
         if (n_knots < 4 .or. n_knots > size(fit%x)) then
            fit%status = r_invalid_input
            return
         end if
         fit%n_knots = n_knots
      else if (present(all_knots)) then
         if (all_knots) then
            fit%n_knots = size(fit%x)
         else
            fit%n_knots = default_smoothing_spline_knots(size(fit%x))
         end if
      else
         fit%n_knots = default_smoothing_spline_knots(size(fit%x))
      end if
      fit%minimum_x = fit%x(1)
      fit%x_range = fit%x(size(fit%x)) - fit%x(1)
      call smoothing_spline_system(fit%x, fit%weights, fit%minimum_x, fit%x_range, &
                                   fit%n_knots, fit%knots, basis, design_crossproduct, &
                                   penalty, fit%ratio)
      if (.not. ieee_is_finite(fit%ratio) .or. fit%ratio <= 0.0_dp) then
         fit%status = r_invalid_input
         return
      end if

      if (present(lambda)) then
         if (.not. ieee_is_finite(lambda) .or. lambda <= 0.0_dp) then
            fit%status = r_invalid_input
            return
         end if
         fit%lambda = lambda
         fit%spar = log(max(lambda/fit%ratio, tiny(1.0_dp)))/log(16.0_dp)/6.0_dp + &
                    1.0_dp/3.0_dp
      else if (present(spar)) then
         if (.not. ieee_is_finite(spar)) then
            fit%status = r_invalid_input
            return
         end if
         fit%spar = spar
         fit%lambda = fit%ratio*16.0_dp**(6.0_dp*fit%spar - 2.0_dp)
      else if (present(degrees_freedom)) then
         target_df = degrees_freedom
         if (.not. ieee_is_finite(target_df) .or. target_df <= 2.0_dp .or. &
             target_df >= real(size(fit%x), dp)) then
            fit%status = r_invalid_input
            return
         end if
         lower_spar = -1.5_dp
         upper_spar = 1.5_dp
         do i = 1, 80
            fit%spar = 0.5_dp*(lower_spar + upper_spar)
            fit%lambda = fit%ratio*16.0_dp**(6.0_dp*fit%spar - 2.0_dp)
            call smoothing_spline_df(basis, design_crossproduct, penalty, fit%weights, &
                                     fit%lambda, fit%degrees_freedom, fit%status)
            if (fit%status /= r_ok) return
            if (fit%degrees_freedom > target_df) then
               lower_spar = fit%spar
            else
               upper_spar = fit%spar
            end if
         end do
         fit%spar = 0.5_dp*(lower_spar + upper_spar)
         fit%lambda = fit%ratio*16.0_dp**(6.0_dp*fit%spar - 2.0_dp)
      else
         fit%automatic = .true.
         call select_smoothing_spline_criterion( &
            fit%y, fit%weights, basis, design_crossproduct, penalty, fit%ratio, &
            within_group_ss, fit%ordinary_cv, fit%spar, fit%lambda, fit%criterion, &
            fit%status)
         if (fit%status /= r_ok) return
      end if
      call solve_smoothing_spline(fit%y, fit%weights, basis, design_crossproduct, &
                                  penalty, fit%lambda, fit%coefficients, &
                                  fit%fitted_values, fit%leverage, fit%status)
      if (fit%status /= r_ok) return
      fit%residuals = fit%y - fit%fitted_values
      fit%degrees_freedom = sum(fit%leverage)
      effective_count = sum(fit%weights)
      fit%gcv = (within_group_ss + sum(fit%weights*fit%residuals**2))/effective_count/ &
                (1.0_dp - fit%degrees_freedom/effective_count)**2
      fit%cv = 0.0_dp
      do i = 1, size(sorted_y)
         group = group_index(i)
         adjusted_leverage = 0.0_dp
         if (fit%weights(group) > 0.0_dp) then
            adjusted_leverage = fit%leverage(group)*sorted_weights(i)/fit%weights(group)
         end if
         fit%cv = fit%cv + sorted_weights(i)*(sorted_y(i) - fit%fitted_values(group))**2/ &
                  (1.0_dp - adjusted_leverage)**2
      end do
      fit%cv = fit%cv/effective_count
      fit%criterion = merge(fit%cv, fit%gcv, fit%ordinary_cv)
   end function smooth_spline

   pure elemental real(dp) function round_to_even(value) result(rounded)
      !! Rounds to the nearest integer-valued real, resolving exact ties toward an even value.
      real(dp), intent(in) :: value !! Finite value to round.
      real(dp) :: fraction, lower

      lower = aint(value)
      if (value < lower) lower = lower - 1.0_dp
      fraction = value - lower
      if (fraction < 0.5_dp) then
         rounded = lower
      else if (fraction > 0.5_dp) then
         rounded = lower + 1.0_dp
      else if (modulo(lower, 2.0_dp) == 0.0_dp) then
         rounded = lower
      else
         rounded = lower + 1.0_dp
      end if
   end function round_to_even

   pure integer function default_smoothing_spline_knots(n) result(n_knots)
      !! Returns R's default reduced-knot count for a given number of distinct predictors.
      integer, intent(in) :: n !! Number of distinct predictor values.
      real(dp), parameter :: log2_50 = 5.6438561897747247_dp
      real(dp), parameter :: log2_100 = 6.6438561897747247_dp
      real(dp), parameter :: log2_140 = 7.1292830169449664_dp
      real(dp), parameter :: log2_200 = 7.6438561897747247_dp

      if (n < 50) then
         n_knots = n
      else if (n < 200) then
         n_knots = int(2.0_dp**(log2_50 + (log2_100 - log2_50)*real(n - 50, dp)/150.0_dp))
      else if (n < 800) then
         n_knots = int(2.0_dp**(log2_100 + (log2_140 - log2_100)*real(n - 200, dp)/600.0_dp))
      else if (n < 3200) then
         n_knots = int(2.0_dp**(log2_140 + (log2_200 - log2_140)*real(n - 800, dp)/2400.0_dp))
      else
         n_knots = int(200.0_dp + real(n - 3200, dp)**0.2_dp)
      end if
   end function default_smoothing_spline_knots

   pure subroutine aggregate_duplicate_predictors(x, y, weights, group_codes, unique_x, &
                                                  unique_y, unique_weights, group_index, &
                                                  within_group_ss)
      !! Combines predictors with equal rounded group codes using weighted response means.
      real(dp), intent(in) :: x(:) !! Sorted predictor values, possibly with duplicates.
      real(dp), intent(in) :: y(:) !! Responses in the same sorted order.
      real(dp), intent(in) :: weights(:) !! Normalized nonnegative observation weights.
      real(dp), intent(in) :: group_codes(:) !! Nondecreasing integer-valued grouping codes.
      real(dp), allocatable, intent(out) :: unique_x(:) !! Strictly increasing predictors.
      real(dp), allocatable, intent(out) :: unique_y(:) !! Weighted response means by predictor.
      real(dp), allocatable, intent(out) :: unique_weights(:) !! Summed weights by predictor.
      integer, allocatable, intent(out) :: group_index(:) !! Group index for each sorted input.
      real(dp), intent(out) :: within_group_ss !! Weighted variation lost during aggregation.
      real(dp) :: weighted_sum, weighted_square_sum
      integer :: first, group, i, n_unique

      n_unique = 1
      do i = 2, size(x)
         if (group_codes(i) /= group_codes(i - 1)) n_unique = n_unique + 1
      end do
      allocate (unique_x(n_unique), unique_y(n_unique), unique_weights(n_unique))
      allocate (group_index(size(x)))
      within_group_ss = 0.0_dp
      first = 1
      group = 0
      do i = 2, size(x) + 1
         if (i <= size(x)) then
            if (group_codes(i) == group_codes(i - 1)) cycle
         end if
         group = group + 1
         group_index(first:i - 1) = group
         unique_x(group) = x(first)
         unique_weights(group) = sum(weights(first:i - 1))
         weighted_sum = sum(weights(first:i - 1)*y(first:i - 1))
         weighted_square_sum = sum(weights(first:i - 1)*y(first:i - 1)**2)
         if (unique_weights(group) > 0.0_dp) then
            unique_y(group) = weighted_sum/unique_weights(group)
            within_group_ss = within_group_ss + weighted_square_sum - &
                              unique_weights(group)*unique_y(group)**2
         else
            unique_y(group) = 0.0_dp
         end if
         first = i
      end do
   end subroutine aggregate_duplicate_predictors

   pure subroutine select_smoothing_spline_criterion( &
      y, weights, basis, design_crossproduct, penalty, ratio, within_group_ss, ordinary_cv, &
      spar, lambda, criterion, status)
      !! Selects `spar` by Brent minimization of generalized or ordinary CV.
      real(dp), intent(in) :: y(:) !! Aggregated responses at distinct predictors.
      real(dp), intent(in) :: weights(:) !! Aggregated normalized observation weights.
      real(dp), intent(in) :: basis(:, :) !! Training B-spline design matrix.
      real(dp), intent(in) :: design_crossproduct(:, :) !! Weighted design product.
      real(dp), intent(in) :: penalty(:, :) !! Roughness Gram matrix.
      real(dp), intent(in) :: ratio !! Scale mapping `spar` to `lambda`.
      real(dp), intent(in) :: within_group_ss !! Weighted variation within duplicate groups.
      logical, intent(in) :: ordinary_cv !! Use leave-one-out CV instead of GCV when true.
      real(dp), intent(out) :: spar !! GCV-selected smoothing parameter.
      real(dp), intent(out) :: lambda !! Corresponding roughness multiplier.
      real(dp), intent(out) :: criterion !! Minimized cross-validation criterion.
      integer, intent(out) :: status !! Zero on success or a numerical failure status.
      real(dp), parameter :: golden_step = 0.3819660112501051518_dp
      real(dp), parameter :: relative_tolerance = 2.0e-8_dp
      real(dp), parameter :: tolerance = 1.0e-4_dp
      real(dp) :: a, b, d, e, fx, fu, fv, fw, p, q, r, tolerance1, tolerance2
      real(dp) :: u, v, w, x, xmid
      integer :: iteration

      a = -1.5_dp
      b = 1.5_dp
      x = a + golden_step*(b - a)
      v = x
      w = x
      e = 0.0_dp
      call smoothing_spline_criterion(y, weights, basis, design_crossproduct, penalty, ratio, &
                                      within_group_ss, ordinary_cv, x, fx)
      fv = fx
      fw = fx
      do iteration = 1, 200
         xmid = 0.5_dp*(a + b)
         tolerance1 = relative_tolerance*abs(x) + tolerance/3.0_dp
         tolerance2 = 2.0_dp*tolerance1
         if (abs(x - xmid) <= tolerance2 - 0.5_dp*(b - a)) exit
         if (abs(e) > tolerance1) then
            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            q = 2.0_dp*(q - r)
            if (q > 0.0_dp) p = -p
            q = abs(q)
            r = e
            e = d
            if (abs(p) >= abs(0.5_dp*q*r) .or. p <= q*(a - x) .or. p >= q*(b - x)) then
               e = merge(b - x, a - x, x < xmid)
               d = golden_step*e
            else
               d = p/q
               u = x + d
               if (u - a < tolerance2 .or. b - u < tolerance2) then
                  d = sign(tolerance1, xmid - x)
               end if
            end if
         else
            e = merge(b - x, a - x, x < xmid)
            d = golden_step*e
         end if
         if (abs(d) >= tolerance1) then
            u = x + d
         else
            u = x + sign(tolerance1, d)
         end if
         call smoothing_spline_criterion(y, weights, basis, design_crossproduct, penalty, &
                                         ratio, within_group_ss, ordinary_cv, u, fu)
         if (fu <= fx) then
            if (u < x) then
               b = x
            else
               a = x
            end if
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
         else
            if (u < x) then
               a = u
            else
               b = u
            end if
            if (fu <= fw .or. w == x) then
               v = w
               fv = fw
               w = u
               fw = fu
            else if (fu <= fv .or. v == x .or. v == w) then
               v = u
               fv = fu
            end if
         end if
      end do
      spar = x
      criterion = fx
      lambda = ratio*16.0_dp**(6.0_dp*spar - 2.0_dp)
      status = merge(r_ok, 1, ieee_is_finite(criterion) .and. criterion < huge(1.0_dp))
   end subroutine select_smoothing_spline_criterion

   pure subroutine smoothing_spline_criterion( &
      y, weights, basis, design_crossproduct, penalty, ratio, within_group_ss, ordinary_cv, &
      spar, criterion)
      !! Evaluates generalized or ordinary cross-validation at one `spar` value.
      real(dp), intent(in) :: y(:) !! Aggregated responses at distinct predictors.
      real(dp), intent(in) :: weights(:) !! Aggregated normalized observation weights.
      real(dp), intent(in) :: basis(:, :) !! Training B-spline design matrix.
      real(dp), intent(in) :: design_crossproduct(:, :) !! Weighted design product.
      real(dp), intent(in) :: penalty(:, :) !! Roughness Gram matrix.
      real(dp), intent(in) :: ratio !! Scale mapping `spar` to `lambda`.
      real(dp), intent(in) :: within_group_ss !! Weighted variation within duplicate groups.
      logical, intent(in) :: ordinary_cv !! Use leave-one-out CV instead of GCV when true.
      real(dp), intent(in) :: spar !! Candidate smoothing parameter.
      real(dp), intent(out) :: criterion !! Candidate cross-validation criterion.
      real(dp), allocatable :: coefficients(:), fitted(:), leverage(:)
      real(dp) :: degrees_freedom, effective_count, lambda
      integer :: status

      lambda = ratio*16.0_dp**(6.0_dp*spar - 2.0_dp)
      call solve_smoothing_spline(y, weights, basis, design_crossproduct, penalty, lambda, &
                                  coefficients, fitted, leverage, status)
      if (status /= r_ok) then
         criterion = huge(1.0_dp)
         return
      end if
      effective_count = sum(weights)
      degrees_freedom = sum(leverage)
      if (ordinary_cv) then
         criterion = sum(weights*((y - fitted)/(1.0_dp - leverage))**2)/real(size(y), dp)
      else
         criterion = (within_group_ss + sum(weights*(y - fitted)**2))/effective_count/ &
                     (1.0_dp - degrees_freedom/effective_count)**2
      end if
   end subroutine smoothing_spline_criterion

   pure function predict_smooth_spline(fit, x_new, derivative) result(prediction)
      !! Evaluates a smoothing spline or its first two derivatives.
      type(smooth_spline_fit_t), intent(in) :: fit !! Previously fitted smoothing spline.
      real(dp), intent(in) :: x_new(:) !! Finite prediction locations.
      integer, intent(in), optional :: derivative !! Derivative order from zero through two.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: basis(:), first_derivative(:), second_derivative(:)
      real(dp) :: endpoint, endpoint_slope, normalized_x
      integer :: derivative_order, i

      allocate (prediction(size(x_new)))
      derivative_order = 0
      if (present(derivative)) derivative_order = derivative
      if (fit%status /= r_ok .or. .not. allocated(fit%knots) .or. &
          .not. allocated(fit%coefficients) .or. fit%x_range <= 0.0_dp .or. &
          .not. all(ieee_is_finite(x_new)) .or. derivative_order < 0 .or. &
          derivative_order > 2) then
         prediction = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      do i = 1, size(x_new)
         normalized_x = (x_new(i) - fit%minimum_x)/fit%x_range
         if (normalized_x < 0.0_dp) then
            endpoint = 0.0_dp
         else if (normalized_x > 1.0_dp) then
            endpoint = 1.0_dp
         else
            select case (derivative_order)
            case (0)
               if (.not. allocated(basis)) allocate (basis(size(fit%coefficients)))
               call bspline_basis(fit%knots, 3, normalized_x, basis)
               prediction(i) = dot_product(fit%coefficients, basis)
            case (1)
               call bspline_basis_and_derivatives(fit%knots, normalized_x, basis, &
                                                  first_derivative)
               prediction(i) = dot_product(fit%coefficients, first_derivative)/fit%x_range
            case (2)
               call bspline_second_derivative(fit%knots, normalized_x, second_derivative)
               prediction(i) = dot_product(fit%coefficients, second_derivative)/fit%x_range**2
            end select
            cycle
         end if
         select case (derivative_order)
         case (0, 1)
            call bspline_basis_and_derivatives(fit%knots, endpoint, basis, first_derivative)
            endpoint_slope = dot_product(fit%coefficients, first_derivative)/fit%x_range
            if (derivative_order == 0) then
               prediction(i) = dot_product(fit%coefficients, basis) + &
                               endpoint_slope*(x_new(i) - fit%minimum_x - &
                               endpoint*fit%x_range)
            else
               prediction(i) = endpoint_slope
            end if
         case (2)
            prediction(i) = 0.0_dp
         end select
      end do
   end function predict_smooth_spline

   pure subroutine smoothing_spline_system(x, weights, minimum_x, x_range, n_knots, knots, &
                                           basis, design_crossproduct, penalty, ratio)
      !! Constructs reduced- or all-knots cubic B-spline design and roughness matrices.
      real(dp), intent(in) :: x(:) !! Sorted distinct training predictors.
      real(dp), intent(in) :: weights(:) !! Normalized observation weights.
      real(dp), intent(in) :: minimum_x !! Lower normalization endpoint.
      real(dp), intent(in) :: x_range !! Positive normalization range.
      integer, intent(in) :: n_knots !! Number of predictor knots in the spline basis.
      real(dp), allocatable, intent(out) :: knots(:) !! Clamped normalized knot sequence.
      real(dp), allocatable, intent(out) :: basis(:, :) !! Training B-spline design matrix.
      real(dp), allocatable, intent(out) :: design_crossproduct(:, :) !! Weighted design product.
      real(dp), allocatable, intent(out) :: penalty(:, :) !! Roughness Gram matrix.
      real(dp), intent(out) :: ratio !! Interior diagonal trace ratio.
      real(dp), allocatable :: delta(:), left_second(:), normalized_x(:), right_second(:)
      real(dp) :: interval_width, trace_design, trace_penalty
      integer :: i, j, knot_index, n, n_basis

      n = size(x)
      n_basis = n_knots + 2
      normalized_x = (x - minimum_x)/x_range
      allocate (knots(n_knots + 6))
      knots(:3) = 0.0_dp
      do i = 1, n_knots
         knot_index = 1 + int(real((i - 1)*(n - 1), dp)/real(n_knots - 1, dp))
         knots(i + 3) = normalized_x(knot_index)
      end do
      knots(n_knots + 4:) = 1.0_dp
      allocate (basis(n, n_basis))
      do i = 1, n
         call bspline_basis(knots, 3, normalized_x(i), basis(i, :))
      end do
      allocate (design_crossproduct(n_basis, n_basis), penalty(n_basis, n_basis), &
                source=0.0_dp)
      do j = 1, n_basis
         do i = 1, n_basis
            design_crossproduct(i, j) = dot_product(weights*basis(:, i), basis(:, j))
         end do
      end do
      do i = 1, n_knots - 1
         interval_width = knots(i + 4) - knots(i + 3)
         call bspline_second_derivative(knots, knots(i + 3), left_second)
         call bspline_second_derivative(knots, knots(i + 4), right_second)
         delta = right_second - left_second
         penalty = penalty + interval_width*( &
            spread(left_second, 2, n_basis)*spread(left_second, 1, n_basis) + &
            0.5_dp*spread(delta, 2, n_basis)*spread(left_second, 1, n_basis) + &
            0.5_dp*spread(left_second, 2, n_basis)*spread(delta, 1, n_basis) + &
            0.3330_dp*spread(delta, 2, n_basis)*spread(delta, 1, n_basis))
      end do
      trace_design = 0.0_dp
      trace_penalty = 0.0_dp
      do i = 3, n_basis - 3
         trace_design = trace_design + design_crossproduct(i, i)
         trace_penalty = trace_penalty + penalty(i, i)
      end do
      ratio = trace_design/trace_penalty
   end subroutine smoothing_spline_system

   pure subroutine solve_smoothing_spline(y, weights, basis, design_crossproduct, penalty, &
                                          lambda, coefficients, fitted, leverage, status)
      !! Solves one penalized spline system and returns fitted values and leverages.
      real(dp), intent(in) :: y(:) !! Responses in sorted predictor order.
      real(dp), intent(in) :: weights(:) !! Normalized observation weights.
      real(dp), intent(in) :: basis(:, :) !! Training B-spline design matrix.
      real(dp), intent(in) :: design_crossproduct(:, :) !! Weighted design product.
      real(dp), intent(in) :: penalty(:, :) !! Roughness Gram matrix.
      real(dp), intent(in) :: lambda !! Nonnegative penalty multiplier.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Fitted B-spline coefficients.
      real(dp), allocatable, intent(out) :: fitted(:) !! Smoothed training responses.
      real(dp), allocatable, intent(out) :: leverage(:) !! Smoothing-matrix diagonal.
      integer, intent(out) :: status !! Zero on success or one for a singular system.
      real(dp), allocatable :: inverse(:, :), rhs(:), system(:, :)
      integer :: i, info

      system = design_crossproduct + lambda*penalty
      rhs = matmul(transpose(basis), weights*y)
      allocate (coefficients(size(rhs)))
      call solve_system(system, rhs, coefficients, info)
      if (info /= 0) then
         status = 1
         return
      end if
      call inverse_matrix(system, inverse, info)
      if (info /= 0) then
         status = 1
         return
      end if
      fitted = matmul(basis, coefficients)
      allocate (leverage(size(y)))
      do i = 1, size(y)
         leverage(i) = weights(i)*dot_product(basis(i, :), matmul(inverse, basis(i, :)))
      end do
      status = r_ok
   end subroutine solve_smoothing_spline

   pure subroutine smoothing_spline_df(basis, design_crossproduct, penalty, weights, &
                                       lambda, degrees_freedom, status)
      !! Computes the smoothing-matrix trace for a candidate penalty multiplier.
      real(dp), intent(in) :: basis(:, :) !! Training B-spline design matrix.
      real(dp), intent(in) :: design_crossproduct(:, :) !! Weighted design product.
      real(dp), intent(in) :: penalty(:, :) !! Roughness Gram matrix.
      real(dp), intent(in) :: weights(:) !! Normalized observation weights.
      real(dp), intent(in) :: lambda !! Candidate penalty multiplier.
      real(dp), intent(out) :: degrees_freedom !! Effective degrees of freedom.
      integer, intent(out) :: status !! Zero on success or one for a singular system.
      real(dp), allocatable :: inverse(:, :), system(:, :)
      integer :: i, info

      system = design_crossproduct + lambda*penalty
      call inverse_matrix(system, inverse, info)
      if (info /= 0) then
         status = 1
         return
      end if
      degrees_freedom = 0.0_dp
      do i = 1, size(weights)
         degrees_freedom = degrees_freedom + &
            weights(i)*dot_product(basis(i, :), matmul(inverse, basis(i, :)))
      end do
      status = r_ok
   end subroutine smoothing_spline_df

   pure subroutine bspline_basis(knots, degree, x, values)
      !! Evaluates every B-spline of a selected degree at one normalized location.
      real(dp), intent(in) :: knots(:) !! Nondecreasing knot sequence.
      integer, intent(in) :: degree !! Nonnegative B-spline degree.
      real(dp), intent(in) :: x !! Evaluation location within the knot range.
      real(dp), intent(out) :: values(:) !! Basis values in coefficient order.
      real(dp), allocatable :: current(:), next(:)
      real(dp) :: left_denominator, right_denominator
      integer :: basis_count, i, pass

      if (x >= knots(size(knots))) then
         values = 0.0_dp
         values(size(values)) = 1.0_dp
         return
      end if
      allocate (current(size(knots) - 1), source=0.0_dp)
      do i = 1, size(current)
         if (knots(i) <= x .and. x < knots(i + 1)) current(i) = 1.0_dp
      end do
      do pass = 1, degree
         basis_count = size(knots) - pass - 1
         allocate (next(basis_count), source=0.0_dp)
         do i = 1, basis_count
            left_denominator = knots(i + pass) - knots(i)
            right_denominator = knots(i + pass + 1) - knots(i + 1)
            if (left_denominator > 0.0_dp) then
               next(i) = next(i) + (x - knots(i))*current(i)/left_denominator
            end if
            if (right_denominator > 0.0_dp) then
               next(i) = next(i) + (knots(i + pass + 1) - x)*current(i + 1)/ &
                         right_denominator
            end if
         end do
         call move_alloc(next, current)
      end do
      values = current
   end subroutine bspline_basis

   pure subroutine bspline_basis_and_derivatives(knots, x, basis, derivative)
      !! Evaluates cubic B-splines and their first derivatives.
      real(dp), intent(in) :: knots(:) !! Clamped cubic knot sequence.
      real(dp), intent(in) :: x !! Normalized evaluation location.
      real(dp), allocatable, intent(out) :: basis(:) !! Cubic basis values.
      real(dp), allocatable, intent(out) :: derivative(:) !! First derivatives.
      real(dp), allocatable :: quadratic(:)
      real(dp) :: denominator, derivative_x
      integer :: i

      allocate (basis(size(knots) - 4), derivative(size(knots) - 4), &
                quadratic(size(knots) - 3))
      call bspline_basis(knots, 3, x, basis)
      derivative_x = x
      if (x >= knots(size(knots))) derivative_x = nearest(knots(size(knots)), -1.0_dp)
      call bspline_basis(knots, 2, derivative_x, quadratic)
      derivative = 0.0_dp
      do i = 1, size(basis)
         denominator = knots(i + 3) - knots(i)
         if (denominator > 0.0_dp) derivative(i) = 3.0_dp*quadratic(i)/denominator
         denominator = knots(i + 4) - knots(i + 1)
         if (denominator > 0.0_dp) then
            derivative(i) = derivative(i) - 3.0_dp*quadratic(i + 1)/denominator
         end if
      end do
   end subroutine bspline_basis_and_derivatives

   pure subroutine bspline_second_derivative(knots, x, second_derivative)
      !! Evaluates second derivatives of every cubic B-spline.
      real(dp), intent(in) :: knots(:) !! Clamped cubic knot sequence.
      real(dp), intent(in) :: x !! Normalized evaluation location.
      real(dp), allocatable, intent(out) :: second_derivative(:) !! Cubic second derivatives.
      real(dp), allocatable :: linear(:), quadratic_derivative(:)
      real(dp) :: denominator, derivative_x
      integer :: i, quadratic_count

      allocate (linear(size(knots) - 2))
      derivative_x = x
      if (x >= knots(size(knots))) derivative_x = nearest(knots(size(knots)), -1.0_dp)
      call bspline_basis(knots, 1, derivative_x, linear)
      quadratic_count = size(knots) - 3
      allocate (quadratic_derivative(quadratic_count), source=0.0_dp)
      do i = 1, quadratic_count
         denominator = knots(i + 2) - knots(i)
         if (denominator > 0.0_dp) quadratic_derivative(i) = 2.0_dp*linear(i)/denominator
         denominator = knots(i + 3) - knots(i + 1)
         if (denominator > 0.0_dp) then
            quadratic_derivative(i) = quadratic_derivative(i) - &
                                      2.0_dp*linear(i + 1)/denominator
         end if
      end do
      allocate (second_derivative(size(knots) - 4), source=0.0_dp)
      do i = 1, size(second_derivative)
         denominator = knots(i + 3) - knots(i)
         if (denominator > 0.0_dp) then
            second_derivative(i) = 3.0_dp*quadratic_derivative(i)/denominator
         end if
         denominator = knots(i + 4) - knots(i + 1)
         if (denominator > 0.0_dp) then
            second_derivative(i) = second_derivative(i) - &
                                   3.0_dp*quadratic_derivative(i + 1)/denominator
         end if
      end do
   end subroutine bspline_second_derivative

end module r_stats_smoothing_spline
