! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of binom 1.1-2 R/binom.optim.R; see NOTICE.md.
module binom_optimization
   use r_kinds, only: dp
   use r_special, only: r_log_beta
   use r_ordering, only: r_sort_values_in_place
   use binom_intervals, only: binom_interval, binom_confint
   use binom_numerics, only: beta_probability, beta_quantile, binom_nan, logistic, logit
   use lbfgsb3_mod, only: lbfgsb_control_t, lbfgsb_result_t, lbfgsb_minimize
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private
   public :: binom_optim, binom_optim_result, binom_optim_control
   public :: binom_optim_objective, integrate_binomial_error, binom_optim_max_n, binom_optim_max_memory
   integer, parameter :: binom_optim_max_n = 200 !! Resource limit for the piecewise optimization objective.
   integer, parameter :: binom_optim_max_memory = 1000 !! Limit preventing excessive or overflowing solver workspace.

   type :: binom_optim_control
      !! R optim L-BFGS-B defaults, plus a strict total objective-evaluation cap.
      integer :: memory = 5 !! L-BFGS-B history length, restricted to 1:binom_optim_max_memory.
      integer :: max_iterations = 100
      integer :: max_evaluations = 10000
      real(dp) :: factr = 1.0e7_dp
      real(dp) :: pgtol = 0.0_dp
      real(dp) :: finite_difference_step = 1.0e-3_dp
   end type binom_optim_control

   type :: binom_optim_result
      !! Arrays use R order, with the full interval arrays indexed from zero.
      real(dp), allocatable :: par(:)
      real(dp), allocatable :: conf_level(:)
      type(binom_interval), allocatable :: confint(:)
      real(dp) :: initial_value = 0.0_dp
      real(dp) :: value = 0.0_dp
      integer :: convergence = 52
      integer :: iterations = 0
      integer :: function_evaluations = 0
      integer :: gradient_evaluations = 0
      integer :: objective_evaluations = 0
      character(len=:), allocatable :: message
   end type binom_optim_result

   type :: optimization_context
      integer :: n = 0
      integer :: evaluations = 0
      integer :: evaluation_cap = 10000
      real(dp) :: coverage = 0.95_dp
      real(dp) :: p0 = 0.0_dp
      real(dp) :: step = 1.0e-3_dp
      logical :: transform = .true.
      logical :: capped = .false.
      character(len=:), allocatable :: method
      real(dp), allocatable :: last_complete_x(:)
      real(dp) :: last_complete_value = 0.0_dp
   end type optimization_context
contains
   pure real(dp) function integrate_binomial_error(x, n, lower, upper, conf_level) result(value)
      !! Integral of (sum_{j in x} dbinom(j,n,p) - conf_level)**2 over the bounds.
      integer, intent(in) :: x(:) !! Nonempty success counts, each in zero through n; duplicates are retained.
      integer, intent(in) :: n !! Nonnegative number of binomial trials.
      real(dp), intent(in), optional :: lower !! First probability bound, default zero; the two bounds are sorted.
      real(dp), intent(in), optional :: upper !! Second probability bound, default one; both bounds must lie in [0,1].
      real(dp), intent(in), optional :: conf_level !! Target coverage probability, default 0.95.
      real(dp) :: lo, hi, coverage, a, b, ri, rj, weight, log_np1
      real(dp) :: squared, linear, cross
      integer :: i, j
      value = binom_nan()
      lo = 0.0_dp
      hi = 1.0_dp
      coverage = 0.95_dp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      if (present(conf_level)) coverage = conf_level
      if (.not. all(ieee_is_finite([lo, hi, coverage]))) return
      if (n < 0 .or. size(x) < 1) return
      if (any(x < 0) .or. any(x > n)) return
      if (min(lo, hi) < 0.0_dp .or. max(lo, hi) > 1.0_dp) return
      if (lo > hi) then
         weight = lo
         lo = hi
         hi = weight
      end if
      log_np1 = log(real(n, dp) + 1.0_dp)
      squared = 0.0_dp
      linear = 0.0_dp
      cross = 0.0_dp
      do i = 1, size(x)
         ri = r_log_beta(real(x(i), dp) + 1.0_dp, real(n - x(i), dp) + 1.0_dp)
         a = 2.0_dp * real(x(i), dp) + 1.0_dp
         b = 2.0_dp * real(n - x(i), dp) + 1.0_dp
         weight = exp(r_log_beta(a, b) - 2.0_dp * (ri + log_np1))
         squared = squared + weight * beta_interval(lo, hi, a, b)
         linear = linear - 2.0_dp * coverage / (real(n, dp) + 1.0_dp) * &
            beta_interval(lo, hi, real(x(i), dp) + 1.0_dp, real(n - x(i), dp) + 1.0_dp)
         do j = i + 1, size(x)
            rj = r_log_beta(real(x(j), dp) + 1.0_dp, real(n - x(j), dp) + 1.0_dp)
            a = real(x(i), dp) + real(x(j), dp) + 1.0_dp
            b = 2.0_dp * real(n, dp) - real(x(i), dp) - real(x(j), dp) + 1.0_dp
            weight = exp(r_log_beta(a, b) - ri - rj - 2.0_dp * log_np1)
            cross = cross + 2.0_dp * weight * beta_interval(lo, hi, a, b)
         end do
      end do
      value = squared + linear + cross + (hi - lo) * coverage ** 2
   end function integrate_binomial_error

   pure elemental real(dp) function beta_interval(lo, hi, a, b) result(value)
      real(dp), intent(in) :: lo !! Lower integration bound in the unit interval.
      real(dp), intent(in) :: hi !! Upper integration bound, at least lo and at most one.
      real(dp), intent(in) :: a !! Positive first beta shape.
      real(dp), intent(in) :: b !! Positive second beta shape.
      if (lo > a / (a + b)) then
         value = beta_probability(lo, a, b, upper_tail=.true.) - &
            beta_probability(hi, a, b, upper_tail=.true.)
      else
         value = beta_probability(hi, a, b) - beta_probability(lo, a, b)
      end if
   end function beta_interval

   pure real(dp) function binom_optim_objective(alpha, n, conf_level, method, p0, transform, half) result(value)
      !! Piecewise integrated squared coverage error; divergent trials or wholly invalid intervals return NaN.
      !! Raw closed-form trials retain R algebra, dropped NaN bounds and three-valued interval membership.
      !! Profile intervals use scalar per-count defaults, avoiding R vector lazy-default/recycling behavior.
      real(dp), intent(in) :: alpha(:) !! Length k <= n/2+1; log odds of error probabilities, or raw errors if transform is false.
      integer, intent(in) :: n !! Number of binomial trials in 1:binom_optim_max_n (200).
      real(dp), intent(in), optional :: conf_level !! Target coverage and fixed confidence levels; default 0.95.
      character(len=*), intent(in), optional :: method !! Confidence-interval method accepted by binom_confint; default lrt.
      real(dp), intent(in), optional :: p0 !! Symmetric integration cutoff in [0,0.5]; default zero.
      logical, intent(in), optional :: transform !! Apply the R binomial logit inverse to alpha; default true.
      logical, intent(in), optional :: half !! Integrate below one half and double, as upstream; default true.
      type(binom_interval), allocatable :: ci(:)
      real(dp), allocatable :: levels(:), bounds(:)
      integer, allocatable :: covered(:)
      real(dp) :: coverage, cutoff, mid, lo, hi, part, error_probability
      integer :: i, j, included, count_bounds, valid_intervals, k
      logical :: linked, half_domain, missing_count, raw_algebra
      character(len=:), allocatable :: ci_method
      value = binom_nan()
      coverage = 0.95_dp
      cutoff = 0.0_dp
      linked = .true.
      half_domain = .true.
      ci_method = 'lrt'
      if (present(conf_level)) coverage = conf_level
      if (present(p0)) cutoff = p0
      if (present(transform)) linked = transform
      if (present(half)) half_domain = half
      if (present(method)) ci_method = method
      raw_algebra = .not. linked .and. raw_algebraic_method(ci_method)
      k = size(alpha)
      if (n < 1 .or. n > binom_optim_max_n) return
      if (k < 1 .or. k > n / 2 + 1) return
      if (.not. all(ieee_is_finite([coverage, cutoff]))) return
      if (coverage <= 0.0_dp .or. coverage >= 1.0_dp) return
      if (cutoff < 0.0_dp .or. cutoff > 0.5_dp) return
      if (any(ieee_is_nan(alpha))) return
      if (.not. all(ieee_is_finite(alpha))) return
      allocate(levels(0:n), ci(0:n), covered(n + 1), bounds(2 * (n + 1) + 1))
      count_bounds = 0
      valid_intervals = 0
      do i = 0, n
         j = min(i, n - i) + 1
         levels(i) = coverage
         if (j <= k) then
            error_probability = alpha(j)
            if (linked) error_probability = r_binomial_inverse(alpha(j))
            levels(i) = 1.0_dp - error_probability
         end if
         ci(i) = optimization_interval(real(i, dp), real(n, dp), ci_method, levels(i))
         if (.not. raw_algebra) then
            if (ci(i)%status /= 0) return
            if (.not. all(ieee_is_finite([ci(i)%lower, ci(i)%upper]))) return
         end if
         if (all(ieee_is_finite([ci(i)%lower, ci(i)%upper]))) valid_intervals = valid_intervals + 1
         ! R sort() drops NA endpoints for raw algebraic interval trials.
         if (.not. ieee_is_nan(ci(i)%lower)) then
            if (.not. ieee_is_finite(ci(i)%lower)) return
            count_bounds = count_bounds + 1
            bounds(count_bounds) = ci(i)%lower
         end if
         if (.not. ieee_is_nan(ci(i)%upper)) then
            if (.not. ieee_is_finite(ci(i)%upper)) return
            count_bounds = count_bounds + 1
            bounds(count_bounds) = ci(i)%upper
         end if
      end do
      ! Reject the R artifact in which no finite interval remains and the empty objective becomes zero.
      if (valid_intervals == 0) return
      if (count_bounds > 0) call sort_unique(bounds, count_bounds)
      if (half_domain) then
         included = count(bounds(:count_bounds) < 0.5_dp)
         count_bounds = included + 1
         bounds(count_bounds) = 0.5_dp
      end if
      value = 0.0_dp
      do i = 1, count_bounds - 1
         mid = bounds(i) + 0.5_dp * (bounds(i + 1) - bounds(i))
         included = 0
         missing_count = .false.
         do j = 0, n
            ! R's three-valued comparisons retain an NA count unless a bound excludes the midpoint.
            if (ieee_is_nan(ci(j)%lower)) then
               if (ieee_is_nan(ci(j)%upper)) then
                  missing_count = .true.
               else if (mid < ci(j)%upper) then
                  missing_count = .true.
               end if
            else if (ieee_is_nan(ci(j)%upper)) then
               if (mid > ci(j)%lower) missing_count = .true.
            else if (mid > ci(j)%lower .and. mid < ci(j)%upper) then
               included = included + 1
               covered(included) = j
            end if
         end do
         if (included == 0 .and. .not. missing_count) then
            value = binom_nan()
            return
         end if
         lo = max(cutoff, bounds(i))
         hi = min(1.0_dp - cutoff, bounds(i + 1))
         if (hi < lo) cycle
         ! R's beta integration errors if sorting removes every selected NA count.
         if (included == 0) then
            value = binom_nan()
            return
         end if
         part = integrate_binomial_error(covered(:included), n, lo, hi, coverage)
         if (ieee_is_nan(part)) then
            value = part
            return
         end if
         value = value + part
      end do
      if (half_domain) value = 2.0_dp * value
   end function binom_optim_objective

   pure elemental logical function raw_algebraic_method(method) result(supported)
      character(len=*), intent(in) :: method !! Method whose upstream closed-form algebra permits raw out-of-domain trials.
      select case (method)
      case ('exact', 'wilson', 'agresti-coull', 'ac', 'asymp', 'asymptotic', 'cloglog', 'logit', 'probit')
         supported = .true.
      case default
         supported = .false.
      end select
   end function raw_algebraic_method

   pure elemental function optimization_interval(x, n, method, level) result(ci)
      !! R raw optimization evaluates algebraic formulas outside the public confidence domain.
      real(dp), intent(in) :: x !! Success count between zero and n, provided by the optimization count grid.
      real(dp), intent(in) :: n !! Positive trial count from the validated optimization request.
      character(len=*), intent(in) :: method !! Interval method requested for the coverage objective.
      real(dp), intent(in) :: level !! Trial confidence; raw algebra may yield reversed intervals or NaN bounds.
      type(binom_interval) :: ci
      real(dp) :: tail, temporary, endpoint
      ci = binom_confint(x, n, method, level)
      if (.not. raw_algebraic_method(method)) return
      if (.not. ieee_is_finite(level)) return
      if (level > 0.0_dp .and. level < 1.0_dp) return
      tail = 0.5_dp * (1.0_dp - level)
      ci%mean = x / n
      ci%lower = binom_nan()
      ci%upper = binom_nan()
      if (method == 'exact') then
         ci%lower = 0.0_dp
         ci%upper = 1.0_dp
         ! Keep R's complement expression: at zero confidence adjacent exact bounds must tie.
         if (x > 0.0_dp) ci%lower = 1.0_dp - beta_quantile(1.0_dp - tail, n - x + 1.0_dp, x)
         if (x < n) ci%upper = 1.0_dp - beta_quantile(tail, n - x, x + 1.0_dp)
      else if ((method == 'cloglog' .or. method == 'logit' .or. method == 'probit') .and. &
         (x == 0.0_dp .or. x == n)) then
         ! R uses the exact boundary formula even when its confidence level is nonstandard.
         endpoint = binom_nan()
         if (tail >= 0.0_dp) endpoint = tail ** (1.0_dp / n)
         if (n == 1.0_dp) endpoint = tail
         ci%lower = 0.0_dp
         ci%upper = 1.0_dp
         if (x == 0.0_dp) ci%upper = 1.0_dp - endpoint
         if (x == n) ci%lower = endpoint
      else if ((method == 'cloglog' .or. method == 'probit') .and. abs(level) == 1.0_dp) then
         ci%lower = 0.5_dp * (1.0_dp - level)
         ci%upper = 0.5_dp * (1.0_dp + level)
      else if (method == 'logit' .and. level == 1.0_dp) then
         ! The R odds formula retains the NaN infinity/infinity upper endpoint.
         ci%lower = 0.0_dp
      else if (method == 'logit' .and. level == -1.0_dp) then
         ci%upper = 0.0_dp
      else if (level == 0.0_dp) then
         ci%lower = x / n
         ci%upper = x / n
      else if (level > -1.0_dp .and. level < 0.0_dp) then
         ! Symmetry of qnorm reverses the ordinary normal-link interval at negative confidence.
         ci = binom_confint(x, n, method, -level)
         temporary = ci%lower
         ci%lower = ci%upper
         ci%upper = temporary
      end if
      if (all(ieee_is_finite([ci%lower, ci%upper]))) ci%status = 0
   end function optimization_interval

   pure subroutine sort_unique(values, count_values)
      real(dp), intent(inout) :: values(:) !! Work array with finite values in its first count_values positions.
      integer, intent(inout) :: count_values !! Number of input values, replaced with the number of unique sorted values.
      integer :: i, unique_count
      call r_sort_values_in_place(values(:count_values))
      unique_count = 1
      do i = 2, count_values
         if (values(i) /= values(unique_count)) then
            unique_count = unique_count + 1
            values(unique_count) = values(i)
         end if
      end do
      count_values = unique_count
   end subroutine sort_unique

   pure elemental real(dp) function r_binomial_inverse(x) result(value)
      real(dp), intent(in) :: x !! Log odds, with R binomial-link tail saturation beyond +/-30.
      value = logistic(x)
      if (x < -30.0_dp) value = epsilon(1.0_dp) / (1.0_dp + epsilon(1.0_dp))
      if (x > 30.0_dp) value = 1.0_dp / (1.0_dp + epsilon(1.0_dp))
   end function r_binomial_inverse

   subroutine binom_optim(n, fit, conf_level, method, k, p0, transform, start, control)
      !! L-BFGS-B optimization using the repository sibling, without the unused R lower/upper bounds.
      integer, intent(in) :: n !! Number of binomial trials in 1:binom_optim_max_n (200).
      type(binom_optim_result), intent(out) :: fit !! Parameters, all n+1 intervals, objective, counts and convergence status.
      real(dp), intent(in), optional :: conf_level !! Target coverage and fixed confidence levels, default 0.95.
      character(len=*), intent(in), optional :: method !! Confidence-interval method accepted by binom_confint, default lrt.
      integer, intent(in), optional :: k !! Number of free symmetric confidence levels, from one through n/2+1; default maximum.
      real(dp), intent(in), optional :: p0 !! Integration cutoff in [0,0.5], default zero.
      logical, intent(in), optional :: transform !! Optimize log odds of errors if true, raw errors otherwise; default true.
      real(dp), intent(in), optional :: start(:) !! Nonempty initial confidence probabilities in (0,1), recycled to length k.
      type(binom_optim_control), intent(in), optional :: control !! L-BFGS-B and fixed central finite-difference controls.
      type(binom_optim_control) :: settings
      type(lbfgsb_control_t) :: solver_control
      type(lbfgsb_result_t) :: solver_result
      type(optimization_context) :: context
      real(dp), allocatable :: alpha(:)
      integer :: free_count, i, j
      fit%value = binom_nan()
      fit%initial_value = binom_nan()
      fit%message = 'invalid optimization input'
      if (n < 1 .or. n > binom_optim_max_n) return
      context%n = n
      context%method = 'lrt'
      if (present(conf_level)) context%coverage = conf_level
      if (present(method)) context%method = method
      if (present(p0)) context%p0 = p0
      if (present(transform)) context%transform = transform
      free_count = n / 2 + 1
      if (present(k)) free_count = k
      if (free_count < 1 .or. free_count > n / 2 + 1) return
      if (.not. all(ieee_is_finite([context%coverage, context%p0]))) return
      if (context%coverage <= 0.0_dp .or. context%coverage >= 1.0_dp) return
      if (context%p0 < 0.0_dp .or. context%p0 > 0.5_dp) return
      settings = binom_optim_control()
      if (present(control)) settings = control
      if (settings%memory <= 0 .or. settings%memory > binom_optim_max_memory) return
      if (settings%max_iterations < 0 .or. settings%max_evaluations < 1) return
      if (.not. all(ieee_is_finite([settings%factr, settings%pgtol, settings%finite_difference_step]))) return
      if (settings%factr < 0.0_dp .or. settings%pgtol < 0.0_dp .or. settings%finite_difference_step <= 0.0_dp) return
      if (settings%finite_difference_step > 0.5_dp * huge(1.0_dp)) return
      context%step = settings%finite_difference_step
      context%evaluation_cap = settings%max_evaluations
      allocate(alpha(free_count))
      alpha = 1.0_dp - context%coverage
      if (present(start)) then
         if (size(start) == 0) return
         if (.not. all(ieee_is_finite(start))) return
         if (any(start <= 0.0_dp) .or. any(start >= 1.0_dp)) return
         do i = 1, free_count
            alpha(i) = 1.0_dp - start(mod(i - 1, size(start)) + 1)
         end do
      end if
      if (context%transform) alpha = logit(alpha)
      fit%initial_value = context_objective(alpha, context)
      if (.not. ieee_is_finite(fit%initial_value)) then
         fit%convergence = 53
         fit%message = 'initial objective is invalid or intervals leave uncovered segments'
         fit%objective_evaluations = context%evaluations
         return
      end if
      context%last_complete_x = alpha
      context%last_complete_value = fit%initial_value
      if (settings%max_iterations == 0) then
         fit%value = fit%initial_value
         fit%convergence = 1
         fit%message = 'maximum iteration count is zero'
      else
         solver_control%memory = settings%memory
         solver_control%factr = settings%factr
         solver_control%pgtol = settings%pgtol
         solver_control%abstol = 0.0_dp
         solver_control%reltol = 0.0_dp
         solver_control%max_iterations = settings%max_iterations
         solver_control%max_evaluations = huge(0)
         call lbfgsb_minimize(alpha, optimizer_gradient, solver_result, control=solver_control, user_data=context)
         fit%value = solver_result%value
         fit%convergence = solver_result%convergence
         fit%message = solver_result%message
         fit%iterations = solver_result%iterations
         fit%function_evaluations = solver_result%function_evaluations
         fit%gradient_evaluations = solver_result%gradient_evaluations
         if (context%capped) then
            alpha = context%last_complete_x
            fit%value = context%last_complete_value
            fit%convergence = 1
            fit%message = 'maximum objective-evaluation count reached'
         end if
      end if
      fit%objective_evaluations = context%evaluations
      allocate(fit%par(free_count), fit%conf_level(0:n), fit%confint(0:n))
      fit%par = 1.0_dp - alpha
      if (context%transform) fit%par = 1.0_dp - r_binomial_inverse(alpha)
      do i = 0, n
         j = min(i, n - i) + 1
         fit%conf_level(i) = context%coverage
         if (j <= free_count) fit%conf_level(i) = fit%par(j)
         fit%confint(i) = optimization_interval(real(i, dp), real(n, dp), context%method, fit%conf_level(i))
      end do
      if (fit%convergence == 0 .and. any(fit%confint%status /= 0)) then
         fit%convergence = 53
         fit%message = 'solver stopped with nonfinite final interval bounds'
      end if
   end subroutine binom_optim

   real(dp) function context_objective(alpha, context) result(value)
      real(dp), intent(in) :: alpha(:) !! Free error probabilities or their log odds.
      type(optimization_context), intent(inout) :: context !! Model settings and actual objective-evaluation counter.
      value = binom_optim_objective(alpha, context%n, context%coverage, context%method, context%p0, context%transform)
      context%evaluations = context%evaluations + 1
   end function context_objective

   subroutine optimizer_gradient(x, f, g, user_data)
      real(dp), intent(in) :: x(:) !! Solver coordinates: errors or log odds of errors.
      real(dp), intent(out) :: f !! Objective at x, or NaN for a divergent evaluation or reached evaluation cap.
      real(dp), intent(out) :: g(:) !! Fixed-step central finite-difference gradient, matching R optim default ndeps.
      class(*), intent(inout), optional :: user_data !! Required optimization_context supplied by binom_optim.
      real(dp) :: plus(size(x)), minus(size(x)), fp, fm
      integer :: j
      f = binom_nan()
      g = binom_nan()
      if (.not. present(user_data)) return
      select type (context => user_data)
      type is (optimization_context)
         if (context%evaluation_cap - context%evaluations < 1 + 2 * size(x)) then
            context%capped = .true.
            return
         end if
         f = context_objective(x, context)
         if (.not. ieee_is_finite(f)) return
         do j = 1, size(x)
            plus = x
            minus = x
            plus(j) = x(j) + context%step
            minus(j) = x(j) - context%step
            fp = context_objective(plus, context)
            fm = context_objective(minus, context)
            if (.not. all(ieee_is_finite([fp, fm]))) return
            g(j) = (fp - fm) / (2.0_dp * context%step)
         end do
         context%last_complete_x = x
         context%last_complete_value = f
      end select
   end subroutine optimizer_gradient
end module binom_optimization
