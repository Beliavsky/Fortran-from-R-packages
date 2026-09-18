! SPDX-License-Identifier: MIT
! SPDX-FileComment: Pure scalar numerical methods corresponding to R stats interfaces.
module r_stats_numerical
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_stats_types, only: integrate_result_t, optimize_result_t, uniroot_result_t
   implicit none
   private

   integer, parameter, public :: numerical_invalid_input = 1, numerical_iteration_limit = 2
   integer, parameter, public :: root_extend_none = 0, root_extend_both = 1
   integer, parameter, public :: root_extend_lower = 2, root_extend_upper = 3

   public :: integrate, optimize, scalar_function_interface, uniroot

   abstract interface
      pure function scalar_function_interface(x) result(value)
         !! Evaluates a real scalar objective or integrand.
         import dp
         real(dp), intent(in) :: x !! Point at which to evaluate the function.
         real(dp) :: value !! Function value at `x`.
      end function scalar_function_interface
   end interface

contains

   pure function uniroot(function, lower, upper, tolerance, max_iterations, extend_interval, &
                         max_extensions) result(search)
      !! Finds a zero in a sign-changing finite interval by safeguarded bisection.
      procedure(scalar_function_interface) :: function !! Pure scalar objective function.
      real(dp), intent(in) :: lower !! Finite lower endpoint of the search interval.
      real(dp), intent(in) :: upper !! Finite upper endpoint of the search interval.
      real(dp), intent(in), optional :: tolerance !! Positive absolute root tolerance.
      integer, intent(in), optional :: max_iterations !! Positive refinement limit.
      integer, intent(in), optional :: extend_interval !! Interval-extension rule; defaults to none.
      integer, intent(in), optional :: max_extensions !! Maximum number of bracket expansions.
      type(uniroot_result_t) :: search !! Root estimate and convergence diagnostics.
      real(dp) :: a, b, fa, fb, fm, midpoint, requested_tolerance, step
      integer :: extension, extension_limit, extension_rule, iteration, iteration_limit

      requested_tolerance = optval(tolerance, epsilon(1.0_dp)**0.25_dp)
      iteration_limit = optval(max_iterations, 1000)
      extension_rule = optval(extend_interval, root_extend_none)
      extension_limit = optval(max_extensions, 100)
      if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. &
          lower >= upper .or. requested_tolerance <= 0.0_dp .or. iteration_limit <= 0 .or. &
          extension_rule < root_extend_none .or. extension_rule > root_extend_upper .or. &
          extension_limit < 0) then
         search%status = numerical_invalid_input
         return
      end if
      a = lower
      b = upper
      fa = function(a)
      fb = function(b)
      if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb) .or. &
          ieee_is_nan(fa) .or. ieee_is_nan(fb)) then
         search%status = numerical_invalid_input
         return
      end if
      if (fa == 0.0_dp) then
         call set_root_result(search, a, fa, 0.0_dp, 0)
         return
      else if (fb == 0.0_dp) then
         call set_root_result(search, b, fb, 0.0_dp, 0)
         return
      end if
      if (same_nonzero_sign(fa, fb)) then
         if (extension_rule == root_extend_none) then
            search%status = numerical_invalid_input
            return
         end if
         step = b - a
         do extension = 1, extension_limit
            select case (extension_rule)
            case (root_extend_both)
               a = a - step
               b = b + step
               fa = function(a)
               fb = function(b)
            case (root_extend_lower)
               a = a - step
               fa = function(a)
            case (root_extend_upper)
               b = b + step
               fb = function(b)
            end select
            if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) then
               search%status = numerical_invalid_input
               return
            end if
            if (fa == 0.0_dp) then
               call set_root_result(search, a, fa, 0.0_dp, 0)
               return
            else if (fb == 0.0_dp) then
               call set_root_result(search, b, fb, 0.0_dp, 0)
               return
            else if (.not. same_nonzero_sign(fa, fb)) then
               exit
            end if
            step = 2.0_dp*step
         end do
         if (same_nonzero_sign(fa, fb)) then
            search%status = numerical_iteration_limit
            return
         end if
      end if
      do iteration = 1, iteration_limit
         midpoint = a + 0.5_dp*(b - a)
         fm = function(midpoint)
         if (.not. ieee_is_finite(fm)) then
            search%status = numerical_invalid_input
            return
         end if
         if (fm == 0.0_dp .or. 0.5_dp*(b - a) <= requested_tolerance) then
            call set_root_result(search, midpoint, fm, 0.5_dp*(b - a), iteration)
            return
         end if
         if ((fa < 0.0_dp .and. fm > 0.0_dp) .or. &
             (fa > 0.0_dp .and. fm < 0.0_dp)) then
            b = midpoint
            fb = fm
         else
            a = midpoint
            fa = fm
         end if
      end do
      search%root = a + 0.5_dp*(b - a)
      search%function_value = function(search%root)
      search%estimated_precision = 0.5_dp*(b - a)
      search%iterations = iteration_limit
      search%status = numerical_iteration_limit
   end function uniroot

   pure function same_nonzero_sign(a, b) result(same)
      !! Reports whether two nonzero finite values have the same sign without multiplying them.
      real(dp), intent(in) :: a !! First finite nonzero value.
      real(dp), intent(in) :: b !! Second finite nonzero value.
      logical :: same !! True when both values are positive or both are negative.

      same = (a > 0.0_dp .and. b > 0.0_dp) .or. (a < 0.0_dp .and. b < 0.0_dp)
   end function same_nonzero_sign

   pure subroutine set_root_result(search, root, function_value, precision, iterations)
      !! Records a successfully converged root search.
      type(uniroot_result_t), intent(out) :: search !! Root result to initialize.
      real(dp), intent(in) :: root !! Converged root estimate.
      real(dp), intent(in) :: function_value !! Objective value at the root estimate.
      real(dp), intent(in) :: precision !! Final half-width of the bracket.
      integer, intent(in) :: iterations !! Number of refinements performed.

      search%root = root
      search%function_value = function_value
      search%estimated_precision = precision
      search%iterations = iterations
      search%status = 0
      search%converged = .true.
   end subroutine set_root_result

   pure function optimize(function, lower, upper, maximum, tolerance, max_iterations) result(search)
      !! Optimizes a scalar function over a finite closed interval by safeguarded Brent search.
      procedure(scalar_function_interface) :: function !! Pure scalar objective function.
      real(dp), intent(in) :: lower !! Finite lower endpoint of the search interval.
      real(dp), intent(in) :: upper !! Finite upper endpoint of the search interval.
      logical, intent(in), optional :: maximum !! Maximizes instead of minimizing when true.
      real(dp), intent(in), optional :: tolerance !! Positive absolute location tolerance.
      integer, intent(in), optional :: max_iterations !! Positive refinement limit.
      type(optimize_result_t) :: search !! Optimum estimate and convergence diagnostics.
      real(dp), parameter :: golden_section = 0.3819660112501051518_dp
      real(dp) :: a, b, d, e, e_previous, fu, fv, fw, fx, midpoint
      real(dp) :: objective_sign, p, q, r, requested_tolerance, tolerance1, tolerance2
      real(dp) :: u, v, w, x
      integer :: iteration, iteration_limit

      requested_tolerance = optval(tolerance, epsilon(1.0_dp)**0.25_dp)
      iteration_limit = optval(max_iterations, 200)
      if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. &
          lower >= upper .or. requested_tolerance <= 0.0_dp .or. iteration_limit <= 0) then
         search%status = numerical_invalid_input
         return
      end if
      objective_sign = merge(-1.0_dp, 1.0_dp, optval(maximum, .false.))
      a = lower
      b = upper
      x = a + golden_section*(b - a)
      w = x
      v = x
      fx = objective_sign*function(x)
      fw = fx
      fv = fx
      d = 0.0_dp
      e = 0.0_dp
      search%evaluations = 1
      if (.not. ieee_is_finite(fx)) then
         search%status = numerical_invalid_input
         return
      end if
      do iteration = 1, iteration_limit
         midpoint = 0.5_dp*(a + b)
         tolerance1 = sqrt(epsilon(1.0_dp))*abs(x) + requested_tolerance/3.0_dp
         tolerance2 = 2.0_dp*tolerance1
         if (abs(x - midpoint) <= tolerance2 - 0.5_dp*(b - a)) then
            search%location = x
            search%objective = objective_sign*fx
            search%status = 0
            search%converged = .true.
            return
         end if
         if (abs(e) > tolerance1) then
            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            q = 2.0_dp*(q - r)
            if (q > 0.0_dp) p = -p
            q = abs(q)
            e_previous = e
            e = d
            if (abs(p) >= abs(0.5_dp*q*e_previous) .or. &
                p <= q*(a - x) .or. p >= q*(b - x)) then
               e = merge(b - x, a - x, x < midpoint)
               d = golden_section*e
            else
               d = p/q
               u = x + d
               if (u - a < tolerance2 .or. b - u < tolerance2) then
                  d = sign(tolerance1, midpoint - x)
               end if
            end if
         else
            e = merge(b - x, a - x, x < midpoint)
            d = golden_section*e
         end if
         u = x + sign(max(abs(d), tolerance1), d)
         fu = objective_sign*function(u)
         search%evaluations = search%evaluations + 1
         if (.not. ieee_is_finite(fu)) then
            search%status = numerical_invalid_input
            return
         end if
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
      search%location = x
      search%objective = objective_sign*fx
      search%status = numerical_iteration_limit
   end function optimize

   pure function integrate(function, lower, upper, abs_tolerance, rel_tolerance, &
                           subdivisions) result(quadrature)
      !! Integrates a pure scalar function over finite or infinite bounds.
      procedure(scalar_function_interface) :: function !! Pure scalar integrand.
      real(dp), intent(in) :: lower !! Lower integration bound, possibly negative infinity.
      real(dp), intent(in) :: upper !! Upper integration bound, possibly positive infinity.
      real(dp), intent(in), optional :: abs_tolerance !! Positive absolute error target.
      real(dp), intent(in), optional :: rel_tolerance !! Positive relative error target.
      integer, intent(in), optional :: subdivisions !! Maximum number of accepted leaf intervals.
      type(integrate_result_t) :: quadrature !! Integral estimate and error diagnostics.

      if (ieee_is_nan(lower) .or. ieee_is_nan(upper) .or. lower >= upper) then
         quadrature%status = numerical_invalid_input
      else if (ieee_is_finite(lower) .and. ieee_is_finite(upper)) then
         quadrature = integrate_finite(function, lower, upper, abs_tolerance, &
                                       rel_tolerance, subdivisions)
      else if (ieee_is_finite(lower) .and. upper > 0.0_dp) then
         quadrature = integrate_finite(right_infinite_transform, 0.0_dp, 1.0_dp, &
                                       abs_tolerance, rel_tolerance, subdivisions)
      else if (lower < 0.0_dp .and. ieee_is_finite(upper)) then
         quadrature = integrate_finite(left_infinite_transform, 0.0_dp, 1.0_dp, &
                                       abs_tolerance, rel_tolerance, subdivisions)
      else if (lower < 0.0_dp .and. upper > 0.0_dp) then
         quadrature = integrate_finite(two_sided_infinite_transform, 0.0_dp, 1.0_dp, &
                                       abs_tolerance, rel_tolerance, subdivisions)
      else
         quadrature%status = numerical_invalid_input
      end if

   contains

      pure function right_infinite_transform(t) result(value)
         !! Maps a semi-infinite interval to the unit interval.
         real(dp), intent(in) :: t !! Unit-interval coordinate.
         real(dp) :: value !! Jacobian-weighted integrand.
         real(dp) :: denominator

         if (t >= 1.0_dp) then
            value = 0.0_dp
         else
            denominator = 1.0_dp - t
            value = function(lower + t/denominator)/denominator**2
         end if
      end function right_infinite_transform

      pure function left_infinite_transform(t) result(value)
         !! Maps a negative semi-infinite interval to the unit interval.
         real(dp), intent(in) :: t !! Unit-interval coordinate.
         real(dp) :: value !! Jacobian-weighted integrand.
         real(dp) :: denominator

         if (t >= 1.0_dp) then
            value = 0.0_dp
         else
            denominator = 1.0_dp - t
            value = function(upper - t/denominator)/denominator**2
         end if
      end function left_infinite_transform

      pure function two_sided_infinite_transform(t) result(value)
         !! Maps the whole real line to the unit interval.
         real(dp), intent(in) :: t !! Unit-interval coordinate.
         real(dp) :: value !! Jacobian-weighted integrand.
         real(dp) :: angle, cosine

         if (t <= 0.0_dp .or. t >= 1.0_dp) then
            value = 0.0_dp
         else
            angle = acos(-1.0_dp)*(t - 0.5_dp)
            cosine = cos(angle)
            value = acos(-1.0_dp)*function(tan(angle))/cosine**2
         end if
      end function two_sided_infinite_transform

   end function integrate

   pure function integrate_finite(function, lower, upper, abs_tolerance, rel_tolerance, &
                                  subdivisions) result(quadrature)
      !! Integrates a pure scalar function over a finite interval with adaptive Simpson quadrature.
      procedure(scalar_function_interface) :: function !! Pure scalar integrand.
      real(dp), intent(in) :: lower !! Finite lower integration bound.
      real(dp), intent(in) :: upper !! Finite upper integration bound.
      real(dp), intent(in), optional :: abs_tolerance !! Positive absolute error target.
      real(dp), intent(in), optional :: rel_tolerance !! Positive relative error target.
      integer, intent(in), optional :: subdivisions !! Maximum number of accepted leaf intervals.
      type(integrate_result_t) :: quadrature !! Integral estimate and error diagnostics.
      real(dp) :: absolute_target, fa, fb, fm, midpoint, relative_target, whole
      integer :: maximum_subdivisions
      logical :: all_converged, valid

      absolute_target = optval(abs_tolerance, epsilon(1.0_dp)**0.25_dp)
      relative_target = optval(rel_tolerance, epsilon(1.0_dp)**0.25_dp)
      maximum_subdivisions = optval(subdivisions, 1000)
      if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. &
          absolute_target <= 0.0_dp .or. relative_target <= 0.0_dp .or. &
          maximum_subdivisions <= 0) then
         quadrature%status = numerical_invalid_input
         return
      end if
      if (lower == upper) then
         quadrature%converged = .true.
         return
      end if
      midpoint = lower + 0.5_dp*(upper - lower)
      fa = function(lower)
      fm = function(midpoint)
      fb = function(upper)
      quadrature%evaluations = 3
      if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fm) .or. &
          .not. ieee_is_finite(fb)) then
         quadrature%status = numerical_invalid_input
         return
      end if
      whole = (upper - lower)*(fa + 4.0_dp*fm + fb)/6.0_dp
      quadrature%subdivisions = 1
      all_converged = .true.
      valid = .true.
      call refine_simpson(function, lower, upper, fa, fm, fb, whole, &
                          max(absolute_target, relative_target*abs(whole)), &
                          maximum_subdivisions, quadrature, all_converged, valid)
      quadrature%converged = all_converged .and. valid
      if (.not. valid) then
         quadrature%status = numerical_invalid_input
      else
         quadrature%status = merge(0, numerical_iteration_limit, all_converged)
      end if
   end function integrate_finite

   pure recursive subroutine refine_simpson(function, lower, upper, fa, fm, fb, whole, &
                                            tolerance, maximum_subdivisions, quadrature, converged, valid)
      !! Recursively bisects one Simpson panel until its local error target is met.
      procedure(scalar_function_interface) :: function !! Pure scalar integrand.
      real(dp), intent(in) :: lower !! Lower panel endpoint.
      real(dp), intent(in) :: upper !! Upper panel endpoint.
      real(dp), intent(in) :: fa !! Integrand at `lower`.
      real(dp), intent(in) :: fm !! Integrand at the panel midpoint.
      real(dp), intent(in) :: fb !! Integrand at `upper`.
      real(dp), intent(in) :: whole !! Simpson estimate on the whole panel.
      real(dp), intent(in) :: tolerance !! Local absolute error target.
      integer, intent(in) :: maximum_subdivisions !! Maximum leaf-panel count.
      type(integrate_result_t), intent(inout) :: quadrature !! Accumulated integral and diagnostics.
      logical, intent(inout) :: converged !! False if refinement exhausts the panel budget.
      logical, intent(inout) :: valid !! False if the integrand returns a nonfinite value.
      real(dp) :: delta, left_estimate, left_midpoint, midpoint, right_estimate, right_midpoint
      real(dp) :: f_left_midpoint, f_right_midpoint, refined

      midpoint = lower + 0.5_dp*(upper - lower)
      left_midpoint = lower + 0.5_dp*(midpoint - lower)
      right_midpoint = midpoint + 0.5_dp*(upper - midpoint)
      f_left_midpoint = function(left_midpoint)
      f_right_midpoint = function(right_midpoint)
      quadrature%evaluations = quadrature%evaluations + 2
      if (.not. ieee_is_finite(f_left_midpoint) .or. .not. ieee_is_finite(f_right_midpoint)) then
         converged = .false.
         valid = .false.
         return
      end if
      left_estimate = (midpoint - lower)*(fa + 4.0_dp*f_left_midpoint + fm)/6.0_dp
      right_estimate = (upper - midpoint)*(fm + 4.0_dp*f_right_midpoint + fb)/6.0_dp
      refined = left_estimate + right_estimate
      delta = refined - whole
      if (abs(delta) <= 15.0_dp*tolerance) then
         quadrature%value = quadrature%value + refined + delta/15.0_dp
         quadrature%absolute_error = quadrature%absolute_error + abs(delta)/15.0_dp
      else if (quadrature%subdivisions >= maximum_subdivisions) then
         quadrature%value = quadrature%value + refined + delta/15.0_dp
         quadrature%absolute_error = quadrature%absolute_error + abs(delta)/15.0_dp
         converged = .false.
      else
         quadrature%subdivisions = quadrature%subdivisions + 1
         call refine_simpson(function, lower, midpoint, fa, f_left_midpoint, fm, &
                             left_estimate, 0.5_dp*tolerance, maximum_subdivisions, &
                             quadrature, converged, valid)
         call refine_simpson(function, midpoint, upper, fm, f_right_midpoint, fb, &
                             right_estimate, 0.5_dp*tolerance, maximum_subdivisions, &
                             quadrature, converged, valid)
      end if
   end subroutine refine_simpson

end module r_stats_numerical
