module rcppnumerical_integration_1d
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rcppnumerical_kinds, only : dp
  use rcppnumerical_callbacks, only : scalar_function_interface
  use rcppnumerical_gk_data, only : get_gk_rule, valid_gk_rule
  implicit none
  private

  integer, parameter, public :: gk15 = 1
  integer, parameter, public :: gk21 = 2
  integer, parameter, public :: gk31 = 3
  integer, parameter, public :: gk41 = 4
  integer, parameter, public :: gk51 = 5
  integer, parameter, public :: gk61 = 6
  integer, parameter, public :: gk71 = 7
  integer, parameter, public :: gk81 = 8
  integer, parameter, public :: gk91 = 9
  integer, parameter, public :: gk101 = 10
  integer, parameter, public :: gk121 = 11
  integer, parameter, public :: gk201 = 12

  type, public :: integration_result_t
    real(dp) :: value = 0.0_dp
    real(dp) :: error_estimate = huge(1.0_dp)
    integer :: error_code = 6
    integer :: evaluations = 0
    integer :: subintervals = 0
  contains
    procedure :: successful => integration_successful
  end type integration_result_t

  public :: integrate_1d

contains

  pure logical function integration_successful(self)
    class(integration_result_t), intent(in) :: self
    integration_successful = self%error_code == 0
  end function integration_successful

  subroutine integrate_1d(f, lower, upper, result, subdiv, eps_abs, eps_rel, rule, user_data)
    procedure(scalar_function_interface) :: f
    real(dp), intent(in) :: lower, upper
    type(integration_result_t), intent(out) :: result
    integer, intent(in), optional :: subdiv, rule
    real(dp), intent(in), optional :: eps_abs, eps_rel
    class(*), intent(inout), optional :: user_data

    integer :: max_sub, local_rule
    real(dp) :: abs_tol, rel_tol, sign, lo, hi
    logical :: lo_finite, hi_finite

    max_sub = 100
    if (present(subdiv)) max_sub = subdiv
    abs_tol = 1.0e-8_dp
    if (present(eps_abs)) abs_tol = eps_abs
    rel_tol = 1.0e-6_dp
    if (present(eps_rel)) rel_tol = eps_rel
    local_rule = gk41
    if (present(rule)) local_rule = rule

    result = integration_result_t()
    if (lower == upper) then
      result%error_estimate = 0.0_dp
      result%error_code = 0
      return
    end if
    if (max_sub < 1 .or. (.not. valid_gk_rule(local_rule)) .or. &
        (abs_tol <= 0.0_dp .and. rel_tol < 50.0_dp*epsilon(1.0_dp))) then
      return
    end if

    sign = 1.0_dp
    lo = lower
    hi = upper
    if (hi < lo) then
      lo = upper
      hi = lower
      sign = -1.0_dp
    end if
    lo_finite = ieee_is_finite(lo)
    hi_finite = ieee_is_finite(hi)

    if (lo_finite .and. hi_finite) then
      call adaptive_integrate(finite_eval, lo, hi, result, max_sub, abs_tol, rel_tol, local_rule)
    else
      call adaptive_integrate(infinite_eval, 0.0_dp, 1.0_dp, result, max_sub, &
                              abs_tol, rel_tol, local_rule)
    end if
    result%value = sign*result%value

  contains

    function finite_eval(x, ignored) result(y)
      real(dp), intent(in) :: x
      class(*), intent(inout), optional :: ignored
      real(dp) :: y
      if (present(user_data)) then
        y = f(x, user_data)
      else
        y = f(x)
      end if
    end function finite_eval

    function infinite_eval(t, ignored) result(y)
      real(dp), intent(in) :: t
      class(*), intent(inout), optional :: ignored
      real(dp) :: y, x, fp, fm
      x = (1.0_dp - t)/t
      if (lo_finite) then
        if (present(user_data)) then
          y = f(lo + x, user_data)/(t*t)
        else
          y = f(lo + x)/(t*t)
        end if
      else if (hi_finite) then
        if (present(user_data)) then
          y = f(hi - x, user_data)/(t*t)
        else
          y = f(hi - x)/(t*t)
        end if
      else
        if (present(user_data)) then
          fp = f(x, user_data)
          fm = f(-x, user_data)
        else
          fp = f(x)
          fm = f(-x)
        end if
        y = (fp + fm)/(t*t)
      end if
    end function infinite_eval

  end subroutine integrate_1d

  subroutine adaptive_integrate(f, lower, upper, result, max_sub, eps_abs, eps_rel, rule)
    procedure(scalar_function_interface) :: f
    real(dp), intent(in) :: lower, upper, eps_abs, eps_rel
    type(integration_result_t), intent(out) :: result
    integer, intent(in) :: max_sub, rule

    real(dp), allocatable :: lowers(:), uppers(:), areas(:), errors(:)
    real(dp) :: area, error_sum, error_bound, error_max
    real(dp) :: lower1, upper1, lower2, upper2
    real(dp) :: area1, area2, error1, error2, absres, absdiff1, absdiff2
    real(dp) :: area12, error12
    integer :: nsub, imax, roundoff1, roundoff2

    allocate(lowers(max_sub), uppers(max_sub), areas(max_sub), errors(max_sub))
    lowers = 0.0_dp
    uppers = 0.0_dp
    areas = 0.0_dp
    errors = 0.0_dp

    call kronrod_rule(f, lower, upper, rule, area, error_sum, absres, absdiff1, &
                      result%evaluations)
    nsub = 1
    lowers(1) = lower
    uppers(1) = upper
    areas(1) = area
    errors(1) = error_sum
    result%value = area
    result%error_estimate = error_sum
    result%subintervals = 1
    result%error_code = 0

    error_bound = max(eps_abs, eps_rel*abs(area))
    if (max_sub == 1 .and. error_sum > error_bound) then
      result%error_code = 1
      return
    end if
    if (error_sum <= 50.0_dp*epsilon(1.0_dp)*absdiff1 .and. error_sum > error_bound) then
      result%error_code = 2
      return
    end if
    if ((error_sum <= error_bound .and. error_sum /= absres) .or. error_sum == 0.0_dp) return

    roundoff1 = 0
    roundoff2 = 0
    do while (nsub < max_sub)
      imax = maxloc(errors(1:nsub), dim=1)
      error_max = errors(imax)
      lower1 = lowers(imax)
      upper2 = uppers(imax)
      upper1 = 0.5_dp*(lower1 + upper2)
      lower2 = upper1

      call kronrod_rule(f, lower1, upper1, rule, area1, error1, absres, absdiff1, &
                        result%evaluations)
      call kronrod_rule(f, lower2, upper2, rule, area2, error2, absres, absdiff2, &
                        result%evaluations)
      area12 = area1 + area2
      error12 = error1 + error2
      error_sum = error_sum + error12 - error_max
      area = area + area12 - areas(imax)

      if (absdiff1 /= error1 .and. absdiff2 /= error2) then
        if (abs(areas(imax) - area12) <= abs(area12)*1.0e-5_dp .and. &
            error12 >= 0.99_dp*error_max) roundoff1 = roundoff1 + 1
        if (nsub + 1 > 10 .and. error12 > error_max) roundoff2 = roundoff2 + 1
      end if

      areas(imax) = area1
      errors(imax) = error1
      uppers(imax) = upper1
      nsub = nsub + 1
      lowers(nsub) = lower2
      uppers(nsub) = upper2
      areas(nsub) = area2
      errors(nsub) = error2

      result%value = area
      result%error_estimate = error_sum
      result%subintervals = nsub
      error_bound = max(eps_abs, eps_rel*abs(area))
      if (error_sum <= error_bound) return
      if (roundoff1 >= 6 .or. roundoff2 >= 20) then
        result%error_code = 2
        return
      end if
      if (max(abs(lower1), abs(upper2)) <= &
          (100.0_dp*epsilon(1.0_dp) + 1.0_dp)* &
          (abs(lower2) + 1000.0_dp*tiny(1.0_dp))) then
        result%error_code = 3
        return
      end if
    end do
    result%error_code = 1
  end subroutine adaptive_integrate

  subroutine kronrod_rule(f, lower, upper, rule, value, error, abs_integral, &
                          abs_diff_integral, evaluations)
    procedure(scalar_function_interface) :: f
    real(dp), intent(in) :: lower, upper
    integer, intent(in) :: rule
    real(dp), intent(out) :: value, error, abs_integral, abs_diff_integral
    integer, intent(inout) :: evaluations

    real(dp), allocatable :: xgk(:), wgk(:), wg(:), fminus(:), fplus(:)
    real(dp) :: half_length, center, fcenter, result_gauss, result_kronrod
    real(dp) :: result_mean, abscissa, tmp
    integer :: i, nk

    call get_gk_rule(rule, xgk, wgk, wg)
    nk = size(xgk)
    allocate(fminus(nk - 1), fplus(nk - 1))
    half_length = 0.5_dp*(upper - lower)
    center = 0.5_dp*(lower + upper)
    fcenter = f(center)
    do i = 1, nk - 1
      abscissa = half_length*xgk(i)
      fminus(i) = f(center - abscissa)
      fplus(i) = f(center + abscissa)
    end do
    evaluations = evaluations + 2*nk - 1

    if (any(rule == [gk15, gk31, gk51, gk71, gk91])) then
      result_gauss = wg(size(wg))*fcenter
    else
      result_gauss = 0.0_dp
    end if
    result_kronrod = wgk(nk)*fcenter
    abs_integral = abs(result_kronrod)
    result_kronrod = result_kronrod + sum((fminus + fplus)*wgk(1:nk - 1))
    result_mean = 0.5_dp*result_kronrod
    abs_diff_integral = wgk(nk)*abs(fcenter - result_mean)

    do i = 1, nk - 1
      if (mod(i - 1, 2) == 1) then
        result_gauss = result_gauss + wg(ishft(i, -1))*(fminus(i) + fplus(i))
      end if
      abs_integral = abs_integral + wgk(i)*(abs(fminus(i)) + abs(fplus(i)))
      abs_diff_integral = abs_diff_integral + wgk(i)* &
        (abs(fminus(i) - result_mean) + abs(fplus(i) - result_mean))
    end do

    value = result_kronrod*half_length
    abs_integral = abs_integral*abs(half_length)
    abs_diff_integral = abs_diff_integral*abs(half_length)
    error = abs((result_kronrod - result_gauss)*half_length)
    if (abs_diff_integral /= 0.0_dp .and. error /= 0.0_dp) then
      tmp = 200.0_dp*error/abs_diff_integral
      error = abs_diff_integral*min(1.0_dp, tmp*sqrt(tmp))
    end if
    if (abs_integral > tiny(1.0_dp)/(50.0_dp*epsilon(1.0_dp))) then
      error = max(50.0_dp*epsilon(1.0_dp)*abs_integral, error)
    end if
  end subroutine kronrod_rule

end module rcppnumerical_integration_1d
