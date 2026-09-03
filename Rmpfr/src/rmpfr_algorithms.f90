module rmpfr_algorithms
  use rmpfr_kinds, only: i64
  use rmpfr_types
  use rmpfr_probability, only: mpfr_pnorm
  implicit none
  private

  integer, parameter, public :: optimize_brent = 1
  integer, parameter, public :: optimize_golden = 2
  integer, parameter, public :: extend_none = 0
  integer, parameter, public :: extend_both = 1
  integer, parameter, public :: extend_down = 2
  integer, parameter, public :: extend_up = 3

  type, public :: mpfr_integral_result
    type(mpfr_real) :: value
    type(mpfr_real) :: abs_error
    integer :: subdivisions = 0
    integer :: order = 0
    logical :: converged = .false.
  end type mpfr_integral_result

  type, public :: mpfr_optimize_result
    type(mpfr_real) :: location
    type(mpfr_real) :: objective
    type(mpfr_real) :: estimated_precision
    integer :: iterations = 0
    logical :: converged = .false.
  end type mpfr_optimize_result

  type, public :: mpfr_root_result
    type(mpfr_real) :: root
    type(mpfr_real) :: f_root
    type(mpfr_real) :: estimated_precision
    integer :: iterations = 0
    logical :: converged = .false.
  end type mpfr_root_result

  type, public :: mpfr_hjk_result
    type(mpfr_real), allocatable :: par(:)
    type(mpfr_real) :: value
    integer :: evaluations = 0
    integer :: iterations = 0
    logical :: converged = .false.
  end type mpfr_hjk_result

  interface assignment(=)
    module procedure assign_integral_result
    module procedure assign_optimize_result
    module procedure assign_root_result
    module procedure assign_hjk_result
  end interface assignment(=)

  public :: assignment(=)

  abstract interface
    function mpfr_scalar_function(x) result(y)
      import :: mpfr_real
      type(mpfr_real), intent(in) :: x !! Arbitrary-precision scalar argument supplied by a numerical algorithm.
      type(mpfr_real) :: y
    end function mpfr_scalar_function

    function mpfr_vector_function(x) result(y)
      import :: mpfr_real
      type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision parameter vector supplied by an optimizer.
      type(mpfr_real) :: y
    end function mpfr_vector_function
  end interface

  public :: mpfr_scalar_function, mpfr_vector_function
  public :: mpfr_integrate_romberg, mpfr_optimize, mpfr_uniroot, mpfr_qnorm
  public :: mpfr_hjk

contains

  subroutine assign_integral_result(lhs, rhs)
    type(mpfr_integral_result), intent(out) :: lhs !! Destination integral-result object receiving a deep copy.
    type(mpfr_integral_result), intent(in) :: rhs !! Source integral-result object to copy.

    lhs%value = rhs%value
    lhs%abs_error = rhs%abs_error
    lhs%subdivisions = rhs%subdivisions
    lhs%order = rhs%order
    lhs%converged = rhs%converged
  end subroutine assign_integral_result

  subroutine assign_optimize_result(lhs, rhs)
    type(mpfr_optimize_result), intent(out) :: lhs !! Destination optimization-result object receiving a deep copy.
    type(mpfr_optimize_result), intent(in) :: rhs !! Source optimization-result object to copy.

    lhs%location = rhs%location
    lhs%objective = rhs%objective
    lhs%estimated_precision = rhs%estimated_precision
    lhs%iterations = rhs%iterations
    lhs%converged = rhs%converged
  end subroutine assign_optimize_result

  subroutine assign_root_result(lhs, rhs)
    type(mpfr_root_result), intent(out) :: lhs !! Destination root-result object receiving a deep copy.
    type(mpfr_root_result), intent(in) :: rhs !! Source root-result object to copy.

    lhs%root = rhs%root
    lhs%f_root = rhs%f_root
    lhs%estimated_precision = rhs%estimated_precision
    lhs%iterations = rhs%iterations
    lhs%converged = rhs%converged
  end subroutine assign_root_result

  subroutine assign_hjk_result(lhs, rhs)
    type(mpfr_hjk_result), intent(out) :: lhs !! Destination Hooke-Jeeves result receiving a deep copy.
    type(mpfr_hjk_result), intent(in) :: rhs !! Source Hooke-Jeeves result to copy.
    integer :: i

    if (allocated(rhs%par)) then
      allocate(lhs%par(size(rhs%par)))
      do i = 1, size(rhs%par)
        lhs%par(i) = rhs%par(i)
      end do
    end if
    lhs%value = rhs%value
    lhs%evaluations = rhs%evaluations
    lhs%iterations = rhs%iterations
    lhs%converged = rhs%converged
  end subroutine assign_hjk_result

  function mpfr_integrate_romberg(f, lower, upper, order, rel_tol, abs_tol, max_order) result(out)
    procedure(mpfr_scalar_function) :: f !! Scalar integrand evaluated at arbitrary-precision abscissae.
    type(mpfr_real), intent(in) :: lower !! Finite lower integration limit.
    type(mpfr_real), intent(in) :: upper !! Finite upper integration limit greater than the lower limit.
    integer, intent(in), optional :: order !! Requested Romberg order; convergence may stop earlier when tolerances are given.
    type(mpfr_real), intent(in), optional :: rel_tol !! Positive relative convergence tolerance.
    type(mpfr_real), intent(in), optional :: abs_tol !! Positive absolute convergence tolerance.
    integer, intent(in), optional :: max_order !! Maximum Romberg order; defaults to 19 as in upstream integrateR.
    type(mpfr_integral_result) :: out
    type(mpfr_real), allocatable :: table(:)
    type(mpfr_real) :: width, spacing, endpoint_average, one, four_power
    type(mpfr_real) :: sum_new, previous, current, error, rtol, atol, x
    integer :: h, j, n, k, p, target_order, max_ord
    logical :: check_convergence, use_order

    if (lower >= upper) error stop "Rmpfr: integration requires lower < upper"
    if (.not. mpfr_is_finite(lower) .or. .not. mpfr_is_finite(upper)) then
      error stop "Rmpfr: integration limits must be finite"
    end if

    p = max(mpfr_precision(lower), mpfr_precision(upper))
    max_ord = 19
    if (present(max_order)) max_ord = max_order
    if (max_ord < 0) error stop "Rmpfr: max_order must be nonnegative"

    use_order = present(order)
    target_order = 13
    if (present(order)) target_order = order
    if (target_order < 0) error stop "Rmpfr: integration order must be nonnegative"
    target_order = min(target_order, max_ord)

    check_convergence = present(rel_tol) .or. present(abs_tol) .or. .not. use_order
    rtol = mpfr_from_string("0.0001220703125", p)
    if (present(rel_tol)) rtol = mpfr_copy(rel_tol, p)
    atol = rtol
    if (present(abs_tol)) atol = mpfr_copy(abs_tol, p)
    if (mpfr_sign(rtol) < 0 .or. mpfr_sign(atol) <= 0) then
      error stop "Rmpfr: invalid integration tolerances"
    end if

    allocate(table(0:max_ord))
    width = mpfr_copy(upper, p) - mpfr_copy(lower, p)
    endpoint_average = (f(lower) + f(upper)) / mpfr_from_integer(2_i64, p)
    table(0) = endpoint_average
    previous = table(0) * width
    current = previous
    error = mpfr_inf(1, p)
    one = mpfr_from_integer(1_i64, p)
    n = 1

    if (target_order == 0 .and. use_order) then
      out%value = current
      out%abs_error = mpfr_zero(1, p)
      out%subdivisions = 2
      out%order = 0
      out%converged = .true.
      return
    end if

    do h = 1, max_ord
      spacing = width / mpfr_from_integer(int(2 * n, i64), p)
      sum_new = mpfr_zero(1, p)
      do k = 1, n
        x = mpfr_copy(lower, p) + mpfr_from_integer(int(2 * k - 1, i64), p) * spacing
        sum_new = sum_new + f(x)
      end do
      table(h) = (sum_new / mpfr_from_integer(int(n, i64), p) + table(h - 1)) / &
                 mpfr_from_integer(2_i64, p)

      four_power = one
      do j = h - 1, 0, -1
        four_power = four_power * mpfr_from_integer(4_i64, p)
        table(j) = table(j + 1) + (table(j + 1) - table(j)) / (four_power - one)
      end do

      current = table(0) * width
      error = mpfr_abs(current - previous)
      out%order = h
      out%subdivisions = 2 * n + 1
      if (check_convergence) then
        if (error < mpfr_minimum(mpfr_abs(current) * rtol, atol)) then
          out%converged = .true.
          exit
        end if
      end if
      if (use_order .and. h >= target_order) then
        out%converged = .true.
        exit
      end if
      previous = current
      n = 2 * n
    end do

    out%value = current
    out%abs_error = error
  end function mpfr_integrate_romberg

  function mpfr_optimize(f, lower, upper, tol, method, maximum, maxiter) result(out)
    procedure(mpfr_scalar_function) :: f !! Scalar objective function evaluated at arbitrary-precision arguments.
    type(mpfr_real), intent(in) :: lower !! Lower endpoint of the closed search interval.
    type(mpfr_real), intent(in) :: upper !! Upper endpoint of the closed search interval.
    type(mpfr_real), intent(in), optional :: tol !! Positive requested absolute optimizer tolerance; defaults to 1e-20.
    integer, intent(in), optional :: method !! Method code optimize_brent or optimize_golden; defaults to Brent.
    logical, intent(in), optional :: maximum !! If true maximize instead of minimize.
    integer, intent(in), optional :: maxiter !! Maximum iteration count; defaults to 1000.
    type(mpfr_optimize_result) :: out
    integer :: chosen_method

    chosen_method = optimize_brent
    if (present(method)) chosen_method = method
    select case (chosen_method)
    case (optimize_brent)
      out = optimize_brent_impl(f, lower, upper, tol, maximum, maxiter)
    case (optimize_golden)
      out = optimize_golden_impl(f, lower, upper, tol, maximum, maxiter)
    case default
      error stop "Rmpfr: unknown optimization method"
    end select
  end function mpfr_optimize

  function optimize_golden_impl(f, lower, upper, tol, maximum, maxiter) result(out)
    procedure(mpfr_scalar_function) :: f !! Scalar objective function evaluated by golden-section search.
    type(mpfr_real), intent(in) :: lower !! Lower search endpoint.
    type(mpfr_real), intent(in) :: upper !! Upper search endpoint.
    type(mpfr_real), intent(in), optional :: tol !! Positive requested absolute search tolerance.
    logical, intent(in), optional :: maximum !! If true maximize instead of minimize.
    integer, intent(in), optional :: maxiter !! Maximum number of interval reductions.
    type(mpfr_optimize_result) :: out
    type(mpfr_real) :: a, b, x2, x3, y2, y3, phi, one, two, five, tolerance, midpoint
    type(mpfr_real) :: objective_value
    integer :: iterations, limit, p
    logical :: maximize

    if (lower > upper) error stop "Rmpfr: optimization requires lower <= upper"
    p = max(mpfr_precision(lower), mpfr_precision(upper))
    tolerance = mpfr_from_string("1e-20", p)
    if (present(tol)) tolerance = mpfr_copy(tol, p)
    if (mpfr_sign(tolerance) <= 0) error stop "Rmpfr: optimization tolerance must be positive"
    maximize = .false.
    if (present(maximum)) maximize = maximum
    limit = 1000
    if (present(maxiter)) limit = maxiter
    if (limit < 1) error stop "Rmpfr: maxiter must be positive"

    a = mpfr_copy(lower, p)
    b = mpfr_copy(upper, p)
    one = mpfr_from_integer(1_i64, p)
    two = mpfr_from_integer(2_i64, p)
    five = mpfr_from_integer(5_i64, p)
    phi = one - (mpfr_sqrt(five) - one) / two
    x2 = a + phi * (b - a)
    x3 = b - phi * (b - a)
    y2 = signed_objective(f, x2, maximize)
    y3 = signed_objective(f, x3, maximize)

    iterations = 0
    do while (x3 - x2 > tolerance .and. iterations < limit)
      iterations = iterations + 1
      if (y3 > y2) then
        b = x3
        x3 = x2
        y3 = y2
        x2 = a + phi * (x3 - a)
        y2 = signed_objective(f, x2, maximize)
      else
        a = x2
        x2 = x3
        y2 = y3
        x3 = b - phi * (b - x2)
        y3 = signed_objective(f, x3, maximize)
      end if
    end do

    midpoint = (x2 + x3) / two
    objective_value = f(midpoint)
    out%location = midpoint
    out%objective = objective_value
    out%estimated_precision = mpfr_abs(x3 - x2)
    out%iterations = iterations
    out%converged = x3 - x2 <= tolerance
  end function optimize_golden_impl

  function optimize_brent_impl(f, lower, upper, tol, maximum, maxiter) result(out)
    procedure(mpfr_scalar_function) :: f !! Scalar objective function evaluated by Brent minimization.
    type(mpfr_real), intent(in) :: lower !! Lower search endpoint.
    type(mpfr_real), intent(in) :: upper !! Upper search endpoint.
    type(mpfr_real), intent(in), optional :: tol !! Positive requested absolute search tolerance.
    logical, intent(in), optional :: maximum !! If true maximize instead of minimize.
    integer, intent(in), optional :: maxiter !! Maximum number of Brent iterations.
    type(mpfr_optimize_result) :: out
    type(mpfr_real) :: a, b, cgold, eps, tolerance, tol1, tol2, tol3, xm
    type(mpfr_real) :: x, w, v, fx, fw, fv, d, e, par_p, par_q, par_r, u, fu, old_e
    type(mpfr_real) :: zero, two, three, five, estimated
    type(mpfr_real) :: abs_par_p, abs_par_step, par_lower, par_upper, gap_left, gap_right
    integer :: iterations, limit, prec_bits
    logical :: maximize, use_parabola

    if (lower > upper) error stop "Rmpfr: optimization requires lower <= upper"
    prec_bits = max(mpfr_precision(lower), mpfr_precision(upper))
    tolerance = mpfr_from_string("1e-20", prec_bits)
    if (present(tol)) tolerance = mpfr_copy(tol, prec_bits)
    if (mpfr_sign(tolerance) <= 0) error stop "Rmpfr: optimization tolerance must be positive"
    maximize = .false.
    if (present(maximum)) maximize = maximum
    limit = 1000
    if (present(maxiter)) limit = maxiter
    if (limit < 1) error stop "Rmpfr: maxiter must be positive"

    zero = mpfr_zero(1, prec_bits)
    two = mpfr_from_integer(2_i64, prec_bits)
    three = mpfr_from_integer(3_i64, prec_bits)
    five = mpfr_from_integer(5_i64, prec_bits)
    a = mpfr_copy(lower, prec_bits)
    b = mpfr_copy(upper, prec_bits)
    cgold = (three - mpfr_sqrt(five)) / two
    eps = mpfr_sqrt(mpfr_from_integer(2_i64, prec_bits) ** (-prec_bits))
    tol3 = tolerance / three

    x = a + cgold * (b - a)
    w = x
    v = x
    fx = signed_objective(f, x, maximize)
    fw = fx
    fv = fx
    d = zero
    e = zero
    estimated = b - a
    iterations = 0
    out%converged = .false.

    do while (iterations < limit)
      iterations = iterations + 1
      xm = (a + b) / two
      tol1 = eps * mpfr_abs(x) + tol3
      tol2 = two * tol1
      estimated = (b - a) / two
      if (mpfr_abs(x - xm) <= tol2 - estimated) then
        out%converged = .true.
        exit
      end if

      par_p = zero
      par_q = zero
      par_r = zero
      use_parabola = .false.
      if (mpfr_abs(e) > tol1) then
        par_r = (x - w) * (fx - fv)
        par_q = (x - v) * (fx - fw)
        par_p = (x - v) * par_q - (x - w) * par_r
        par_q = two * (par_q - par_r)
        if (par_q > zero) then
          par_p = -par_p
        else
          par_q = -par_q
        end if
        old_e = e
        e = d
        abs_par_p = mpfr_abs(par_p)
        abs_par_step = mpfr_abs(par_q * old_e / two)
        par_lower = par_q * (a - x)
        par_upper = par_q * (b - x)
        use_parabola = abs_par_p < abs_par_step
        if (use_parabola) use_parabola = par_p > par_lower
        if (use_parabola) use_parabola = par_p < par_upper
      end if

      if (use_parabola) then
        d = par_p / par_q
        u = x + d
        gap_left = u - a
        gap_right = b - u
        if (gap_left < tol2 .or. gap_right < tol2) then
          d = tol1
          if (x >= xm) d = -d
        end if
      else
        if (x < xm) then
          e = b - x
        else
          e = a - x
        end if
        d = cgold * e
      end if

      if (mpfr_abs(d) >= tol1) then
        u = x + d
      else if (d > zero) then
        u = x + tol1
      else
        u = x - tol1
      end if
      fu = signed_objective(f, u, maximize)

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

    out%location = x
    out%objective = f(x)
    out%estimated_precision = mpfr_abs(estimated)
    out%iterations = iterations
  end function optimize_brent_impl

  function mpfr_uniroot(f, lower, upper, tol, maxiter, extend_mode) result(out)
    procedure(mpfr_scalar_function) :: f !! Scalar function whose sign-changing root is sought.
    type(mpfr_real), intent(in) :: lower !! Initial lower endpoint of the root-search interval.
    type(mpfr_real), intent(in) :: upper !! Initial upper endpoint of the root-search interval.
    type(mpfr_real), intent(in), optional :: tol !! Positive absolute stopping tolerance.
    integer, intent(in), optional :: maxiter !! Maximum root iterations and interval-extension steps.
    integer, intent(in), optional :: extend_mode !! Interval-extension code extend_none, extend_both, extend_down, or extend_up.
    type(mpfr_root_result) :: out
    type(mpfr_real) :: a, b, c, fa, fb, fc, tolerance, epsc, tol2, dprev, dnew
    type(mpfr_real) :: par_p, par_q, t1, t2, cb, delta, zero, one, two
    type(mpfr_real) :: hundred, tiny_start, abs_dprev, abs_dnew, abs_fa, abs_fb, interp_limit, prev_limit
    integer :: p, limit, mode, iterations, ext_count

    if (lower >= upper) error stop "Rmpfr: root interval requires lower < upper"
    p = max(mpfr_precision(lower), mpfr_precision(upper))
    zero = mpfr_zero(1, p)
    one = mpfr_from_integer(1_i64, p)
    two = mpfr_from_integer(2_i64, p)
    hundred = mpfr_from_integer(100_i64, p)
    tiny_start = mpfr_from_string("1e-4", p)
    tolerance = mpfr_from_string("0.0001220703125", p)
    if (present(tol)) tolerance = mpfr_copy(tol, p)
    if (mpfr_sign(tolerance) <= 0) error stop "Rmpfr: root tolerance must be positive"
    limit = 1000
    if (present(maxiter)) limit = maxiter
    if (limit < 1) error stop "Rmpfr: maxiter must be positive"
    mode = extend_none
    if (present(extend_mode)) mode = extend_mode

    a = mpfr_copy(lower, p)
    b = mpfr_copy(upper, p)
    fa = f(a)
    fb = f(b)
    if (mpfr_is_nan(fa) .or. mpfr_is_nan(fb)) error stop "Rmpfr: root endpoint evaluates to NaN"

    ext_count = 0
    select case (mode)
    case (extend_none)
      continue
    case (extend_both)
      do while (same_nonzero_sign(fa, fb))
        if (ext_count >= limit) error stop "Rmpfr: no sign change found while extending root interval"
        delta = mpfr_abs(a)
        if (delta < tiny_start) delta = tiny_start
        delta = delta / hundred
        a = a - delta
        fa = f(a)
        delta = mpfr_abs(b)
        if (delta < tiny_start) delta = tiny_start
        delta = delta / hundred
        b = b + delta
        fb = f(b)
        ext_count = ext_count + 1
      end do
    case (extend_down)
      delta = mpfr_abs(a)
      if (delta < tiny_start) delta = tiny_start
      delta = delta / hundred
      do while (mpfr_sign(fa) > 0)
        if (ext_count >= limit) error stop "Rmpfr: no sign change found while extending root interval downward"
        a = a - delta
        fa = f(a)
        delta = two * delta
        ext_count = ext_count + 1
      end do
    case (extend_up)
      delta = mpfr_abs(b)
      if (delta < tiny_start) delta = tiny_start
      delta = delta / hundred
      do while (mpfr_sign(fb) < 0)
        if (ext_count >= limit) error stop "Rmpfr: no sign change found while extending root interval upward"
        b = b + delta
        fb = f(b)
        delta = two * delta
        ext_count = ext_count + 1
      end do
    case default
      error stop "Rmpfr: invalid root interval-extension mode"
    end select

    if (same_nonzero_sign(fa, fb)) error stop "Rmpfr: root endpoints do not bracket a sign change"
    if (mpfr_is_zero(fa)) then
      out%root = a
      out%f_root = fa
      out%estimated_precision = zero
      out%converged = .true.
      return
    end if
    if (mpfr_is_zero(fb)) then
      out%root = b
      out%f_root = fb
      out%estimated_precision = zero
      out%converged = .true.
      return
    end if

    epsc = mpfr_from_integer(2_i64, p) ** (-min(mpfr_precision(fa), mpfr_precision(fb)))
    c = a
    fc = fa
    iterations = 0
    out%converged = .false.

    do while (iterations < limit)
      dprev = b - a
      if (mpfr_abs(fc) < mpfr_abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      end if
      tol2 = two * epsc * mpfr_abs(b) + tolerance / two
      dnew = (c - b) / two
      abs_dnew = mpfr_abs(dnew)
      if (abs_dnew <= tol2 .or. mpfr_is_zero(fb)) then
        out%converged = .true.
        exit
      end if

      abs_dprev = mpfr_abs(dprev)
      abs_fa = mpfr_abs(fa)
      abs_fb = mpfr_abs(fb)
      if (abs_dprev >= tol2 .and. abs_fa > abs_fb) then
        cb = c - b
        if (a == c) then
          t1 = fb / fa
          par_p = cb * t1
          par_q = one - t1
        else
          par_q = fa / fc
          t1 = fb / fc
          t2 = fb / fa
          par_p = t2 * (cb * par_q * (par_q - t1) - (b - a) * (t1 - one))
          par_q = (par_q - one) * (t1 - one) * (t2 - one)
        end if
        if (par_p > zero) then
          par_q = -par_q
        else
          par_p = -par_p
        end if
        interp_limit = mpfr_from_string("0.75", p) * cb * par_q - mpfr_abs(tol2 * par_q) / two
        prev_limit = mpfr_abs(dprev * par_q / two)
        if (par_p < interp_limit .and. par_p < prev_limit) then
          dnew = par_p / par_q
        end if
      end if

      if (mpfr_abs(dnew) < tol2) then
        if (dnew > zero) then
          dnew = tol2
        else
          dnew = -tol2
        end if
      end if
      a = b
      fa = fb
      b = b + dnew
      fb = f(b)
      iterations = iterations + 1
      if (same_nonzero_sign(fb, fc)) then
        c = a
        fc = fa
      end if
    end do

    out%root = b
    out%f_root = f(b)
    out%estimated_precision = mpfr_abs(c - b)
    out%iterations = iterations
  end function mpfr_uniroot

  function mpfr_qnorm(probability, mean, sd, lower_tail, log_p, tol, maxiter) result(out)
    type(mpfr_real), intent(in) :: probability !! Requested normal probability, or its logarithm when log_p is true.
    type(mpfr_real), intent(in), optional :: mean !! Normal mean; defaults to zero.
    type(mpfr_real), intent(in), optional :: sd !! Positive normal standard deviation; defaults to one.
    logical, intent(in), optional :: lower_tail !! If true invert the lower tail; defaults to true.
    logical, intent(in), optional :: log_p !! If true probability is supplied on the logarithmic scale.
    type(mpfr_real), intent(in), optional :: tol !! Positive root-finding tolerance.
    integer, intent(in), optional :: maxiter !! Maximum bracketing/root iterations.
    type(mpfr_root_result) :: out
    type(mpfr_real) :: mu, sigma, target, pzero, pone, lo, hi, flo, fhi, step
    integer :: p, limit, iter
    logical :: lower, logarithmic

    p = mpfr_precision(probability)
    if (present(mean)) p = max(p, mpfr_precision(mean))
    if (present(sd)) p = max(p, mpfr_precision(sd))
    mu = mpfr_zero(1, p)
    if (present(mean)) mu = mpfr_copy(mean, p)
    sigma = mpfr_from_integer(1_i64, p)
    if (present(sd)) sigma = mpfr_copy(sd, p)
    if (mpfr_sign(sigma) <= 0) then
      out%root = mpfr_nan(p)
      out%f_root = mpfr_nan(p)
      out%estimated_precision = mpfr_nan(p)
      return
    end if
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    logarithmic = .false.
    if (present(log_p)) logarithmic = log_p
    target = mpfr_copy(probability, p)
    pzero = mpfr_zero(1, p)
    pone = mpfr_from_integer(1_i64, p)

    if (logarithmic) then
      if (target > pzero) then
        out%root = mpfr_nan(p)
        return
      end if
      if (mpfr_is_infinite(target) .and. mpfr_sign(target) < 0) then
        if (lower) then
          out%root = mpfr_inf(-1, p)
        else
          out%root = mpfr_inf(1, p)
        end if
        out%f_root = pzero
        out%estimated_precision = pzero
        out%converged = .true.
        return
      end if
      if (mpfr_is_zero(target)) then
        if (lower) then
          out%root = mpfr_inf(1, p)
        else
          out%root = mpfr_inf(-1, p)
        end if
        out%f_root = pzero
        out%estimated_precision = pzero
        out%converged = .true.
        return
      end if
    else
      if (target < pzero .or. target > pone) then
        out%root = mpfr_nan(p)
        return
      end if
      if (mpfr_is_zero(target) .or. target == pone) then
        if ((mpfr_is_zero(target) .and. lower) .or. (target == pone .and. .not. lower)) then
          out%root = mpfr_inf(-1, p)
        else
          out%root = mpfr_inf(1, p)
        end if
        out%f_root = pzero
        out%estimated_precision = pzero
        out%converged = .true.
        return
      end if
    end if

    lo = mu - sigma
    hi = mu + sigma
    step = sigma
    limit = 1000
    if (present(maxiter)) limit = maxiter
    do iter = 1, limit
      flo = qnorm_residual(lo, target, mu, sigma, lower, logarithmic)
      fhi = qnorm_residual(hi, target, mu, sigma, lower, logarithmic)
      if (.not. same_nonzero_sign(flo, fhi)) exit
      lo = lo - step
      hi = hi + step
      step = step * mpfr_from_integer(2_i64, p)
    end do
    if (same_nonzero_sign(flo, fhi)) error stop "Rmpfr: qnorm could not bracket requested probability"
    out = mpfr_uniroot_qnorm(target, mu, sigma, lower, logarithmic, lo, hi, tol, maxiter)
  end function mpfr_qnorm

  function mpfr_uniroot_qnorm(target, mean, sd, lower_tail, log_p, lower, upper, tol, maxiter) result(out)
    type(mpfr_real), intent(in) :: target !! Target normal probability on the requested scale.
    type(mpfr_real), intent(in) :: mean !! Normal mean used by the inversion callback.
    type(mpfr_real), intent(in) :: sd !! Positive normal standard deviation used by the inversion callback.
    logical, intent(in) :: lower_tail !! Tail convention used by the inversion callback.
    logical, intent(in) :: log_p !! Probability-scale convention used by the inversion callback.
    type(mpfr_real), intent(in) :: lower !! Bracketing lower quantile.
    type(mpfr_real), intent(in) :: upper !! Bracketing upper quantile.
    type(mpfr_real), intent(in), optional :: tol !! Positive root-finding tolerance.
    integer, intent(in), optional :: maxiter !! Maximum Brent root iterations.
    type(mpfr_root_result) :: out
    type(mpfr_real) :: a, b, c, fa, fb, fc, tolerance, epsc, tol2, dprev, dnew
    type(mpfr_real) :: par_p, par_q, t1, t2, cb, zero, one, two
    type(mpfr_real) :: abs_dprev, abs_dnew, abs_fa, abs_fb, interp_limit, prev_limit
    integer :: p, limit, iterations

    p = max(mpfr_precision(lower), mpfr_precision(upper))
    zero = mpfr_zero(1, p)
    one = mpfr_from_integer(1_i64, p)
    two = mpfr_from_integer(2_i64, p)
    tolerance = two ** (-(p + 2))
    if (present(tol)) tolerance = mpfr_copy(tol, p)
    limit = 1000
    if (present(maxiter)) limit = maxiter
    a = mpfr_copy(lower, p)
    b = mpfr_copy(upper, p)
    fa = qnorm_residual(a, target, mean, sd, lower_tail, log_p)
    fb = qnorm_residual(b, target, mean, sd, lower_tail, log_p)
    c = a
    fc = fa
    epsc = two ** (-min(mpfr_precision(fa), mpfr_precision(fb)))
    iterations = 0
    out%converged = .false.

    do while (iterations < limit)
      dprev = b - a
      if (mpfr_abs(fc) < mpfr_abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      end if
      tol2 = two * epsc * mpfr_abs(b) + tolerance / two
      dnew = (c - b) / two
      abs_dnew = mpfr_abs(dnew)
      if (abs_dnew <= tol2 .or. mpfr_is_zero(fb)) then
        out%converged = .true.
        exit
      end if
      abs_dprev = mpfr_abs(dprev)
      abs_fa = mpfr_abs(fa)
      abs_fb = mpfr_abs(fb)
      if (abs_dprev >= tol2 .and. abs_fa > abs_fb) then
        cb = c - b
        if (a == c) then
          t1 = fb / fa
          par_p = cb * t1
          par_q = one - t1
        else
          par_q = fa / fc
          t1 = fb / fc
          t2 = fb / fa
          par_p = t2 * (cb * par_q * (par_q - t1) - (b - a) * (t1 - one))
          par_q = (par_q - one) * (t1 - one) * (t2 - one)
        end if
        if (par_p > zero) then
          par_q = -par_q
        else
          par_p = -par_p
        end if
        interp_limit = mpfr_from_string("0.75", p) * cb * par_q - mpfr_abs(tol2 * par_q) / two
        prev_limit = mpfr_abs(dprev * par_q / two)
        if (par_p < interp_limit .and. par_p < prev_limit) dnew = par_p / par_q
      end if
      if (mpfr_abs(dnew) < tol2) then
        if (dnew > zero) then
          dnew = tol2
        else
          dnew = -tol2
        end if
      end if
      a = b
      fa = fb
      b = b + dnew
      fb = qnorm_residual(b, target, mean, sd, lower_tail, log_p)
      iterations = iterations + 1
      if (same_nonzero_sign(fb, fc)) then
        c = a
        fc = fa
      end if
    end do

    out%root = b
    out%f_root = qnorm_residual(b, target, mean, sd, lower_tail, log_p)
    out%estimated_precision = mpfr_abs(c - b)
    out%iterations = iterations
  end function mpfr_uniroot_qnorm

  function qnorm_residual(x, target, mean, sd, lower_tail, log_p) result(r)
    type(mpfr_real), intent(in) :: x !! Trial quantile used while inverting the normal CDF.
    type(mpfr_real), intent(in) :: target !! Target normal probability on the requested scale.
    type(mpfr_real), intent(in) :: mean !! Normal mean.
    type(mpfr_real), intent(in) :: sd !! Positive normal standard deviation.
    logical, intent(in) :: lower_tail !! Tail convention for the normal CDF.
    logical, intent(in) :: log_p !! Probability-scale convention for the normal CDF.
    type(mpfr_real) :: r

    r = mpfr_pnorm(x, mean=mean, sd=sd, lower_tail=lower_tail, log_p=log_p) - target
  end function qnorm_residual

  subroutine mpfr_hjk(par, f, out, tol, maxfeval, maximize, target)
    type(mpfr_real), intent(in) :: par(:) !! Initial parameter vector; at least two parameters are required.
    procedure(mpfr_vector_function) :: f !! Multivariate objective function evaluated by Hooke-Jeeves search.
    type(mpfr_hjk_result), intent(out) :: out !! Returned optimum, objective value, counts, and convergence flag.
    type(mpfr_real), intent(in), optional :: tol !! Positive final exploratory step size; defaults to 1e-6.
    integer, intent(in), optional :: maxfeval !! Maximum objective evaluations; defaults to a large finite count.
    logical, intent(in), optional :: maximize !! If true maximize instead of minimize.
    type(mpfr_real), intent(in), optional :: target !! Absolute scaled-objective limit used as an early-stop guard.
    type(mpfr_real), allocatable :: x(:), xb(:), xc(:), direction(:)
    type(mpfr_real) :: fx, fb, step, tolerance, target_limit, scale, two, abs_fx
    integer :: p, n, i, evals, limit, iter
    logical :: found, maximize_flag

    n = size(par)
    if (n < 2) error stop "Rmpfr: Hooke-Jeeves requires at least two parameters"
    p = vector_precision(par)
    tolerance = mpfr_from_string("1e-6", p)
    if (present(tol)) tolerance = mpfr_copy(tol, p)
    if (mpfr_sign(tolerance) <= 0) error stop "Rmpfr: Hooke-Jeeves tolerance must be positive"
    limit = 500000000
    if (present(maxfeval)) limit = maxfeval
    if (limit < 1) error stop "Rmpfr: maxfeval must be positive"
    maximize_flag = .false.
    if (present(maximize)) maximize_flag = maximize
    scale = mpfr_from_integer(1_i64, p)
    if (maximize_flag) scale = -scale
    target_limit = mpfr_inf(1, p)
    if (present(target)) target_limit = mpfr_copy(target, p)
    two = mpfr_from_integer(2_i64, p)
    allocate(x(n), xb(n), xc(n), direction(n))
    call copy_mpfr_vector(x, par)
    fx = scale * f(x)
    evals = 1
    iter = 0
    step = mpfr_from_integer(1_i64, p)

    do
      if (step < tolerance) exit
      if (evals >= limit) exit
      abs_fx = mpfr_abs(fx)
      if (abs_fx >= target_limit) exit
      iter = iter + 1
      call copy_mpfr_vector(xb, x)
      call copy_mpfr_vector(xc, x)
      call hj_explore(xb, xc, step, f, scale, limit, evals, target_limit, x, fx, found)
      do
        if (.not. found) exit
        if (evals >= limit) exit
        abs_fx = mpfr_abs(fx)
        if (abs_fx >= target_limit) exit
        do i = 1, n
          direction(i) = x(i) - xb(i)
        end do
        call copy_mpfr_vector(xb, x)
        do i = 1, n
          xc(i) = x(i) + direction(i)
        end do
        fb = fx
        call hj_explore(xb, xc, step, f, scale, limit, evals, target_limit, x, fx, found, fb)
        if (.not. found .and. evals < limit) then
          call hj_explore(xb, xb, step, f, scale, limit, evals, target_limit, x, fx, found, fb)
        end if
      end do
      step = step / two
    end do

    allocate(out%par(n))
    call copy_mpfr_vector(out%par, x)
    out%value = fx / scale
    out%evaluations = evals
    out%iterations = iter
    abs_fx = mpfr_abs(fx)
    out%converged = evals <= limit .and. abs_fx <= target_limit
  end subroutine mpfr_hjk

  subroutine hj_explore(xb, xc, h, f, scale, maxfeval, evals, target, x, fx, found, known_fb)
    type(mpfr_real), intent(in) :: xb(:) !! Base point for the exploratory Hooke-Jeeves move.
    type(mpfr_real), intent(in) :: xc(:) !! Candidate center from which coordinate directions are explored.
    type(mpfr_real), intent(in) :: h !! Positive coordinate step size.
    procedure(mpfr_vector_function) :: f !! Multivariate objective function evaluated during exploration.
    type(mpfr_real), intent(in) :: scale !! +1 for minimization or -1 for maximization.
    integer, intent(in) :: maxfeval !! Maximum total number of objective evaluations permitted.
    integer, intent(inout) :: evals !! Running objective-evaluation count, updated in place.
    type(mpfr_real), intent(in) :: target !! Absolute scaled-objective early-stop threshold.
    type(mpfr_real), intent(out) :: x(:) !! Best point found by the exploratory move.
    type(mpfr_real), intent(out) :: fx !! Scaled objective at the returned point.
    logical, intent(out) :: found !! True when the exploratory move improves the base objective.
    type(mpfr_real), intent(in), optional :: known_fb !! Optional already-computed scaled objective at the base point.
    type(mpfr_real), allocatable :: xt(:), candidate(:)
    type(mpfr_real) :: fb, fp, abs_fb
    integer :: k, n

    n = size(xb)
    if (size(xc) /= n .or. size(x) /= n) error stop "Rmpfr: Hooke-Jeeves vector shape mismatch"
    allocate(xt(n), candidate(n))
    if (present(known_fb)) then
      fb = known_fb
    else
      fb = scale * f(xb)
      evals = evals + 1
    end if
    call copy_mpfr_vector(xt, xc)
    found = .false.

    do k = 1, n
      if (evals >= maxfeval) exit
      abs_fb = mpfr_abs(fb)
      if (abs_fb >= target) exit
      call copy_mpfr_vector(candidate, xt)
      candidate(k) = candidate(k) + h
      fp = scale * f(candidate)
      evals = evals + 1
      if (fp >= fb .and. evals < maxfeval) then
        call copy_mpfr_vector(candidate, xt)
        candidate(k) = candidate(k) - h
        fp = scale * f(candidate)
        evals = evals + 1
      end if
      if (fp < fb) then
        found = .true.
        call copy_mpfr_vector(xt, candidate)
        fb = fp
      end if
    end do

    if (found) then
      call copy_mpfr_vector(x, xt)
      fx = fb
    else
      call copy_mpfr_vector(x, xb)
      fx = fb
    end if
  end subroutine hj_explore

  function signed_objective(f, x, maximize) result(y)
    procedure(mpfr_scalar_function) :: f !! Scalar objective evaluated at x.
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision objective argument.
    logical, intent(in) :: maximize !! If true negate the objective for minimization machinery.
    type(mpfr_real) :: y

    y = f(x)
    if (maximize) y = -y
  end function signed_objective

  function mpfr_minimum(a, b) result(r)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision value.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision value.
    type(mpfr_real) :: r

    if (a <= b) then
      r = a
    else
      r = b
    end if
  end function mpfr_minimum

  logical function same_nonzero_sign(a, b) result(same)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision value whose sign is compared.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision value whose sign is compared.
    integer :: sa, sb

    sa = mpfr_sign(a)
    sb = mpfr_sign(b)
    same = sa /= 0 .and. sb /= 0 .and. sa == sb
  end function same_nonzero_sign

  subroutine copy_mpfr_vector(destination, source)
    type(mpfr_real), intent(out) :: destination(:) !! Allocated destination vector receiving independent MPFR values.
    type(mpfr_real), intent(in) :: source(:) !! Source vector whose elements are deep-copied.
    integer :: i

    if (size(destination) /= size(source)) error stop "Rmpfr: vector copy shape mismatch"
    do i = 1, size(source)
      destination(i) = source(i)
    end do
  end subroutine copy_mpfr_vector

  integer function vector_precision(x) result(p)
    type(mpfr_real), intent(in) :: x(:) !! Arbitrary-precision vector inspected for maximum element precision.
    integer :: i

    p = 2
    do i = 1, size(x)
      p = max(p, mpfr_precision(x(i)))
    end do
  end function vector_precision

end module rmpfr_algorithms
