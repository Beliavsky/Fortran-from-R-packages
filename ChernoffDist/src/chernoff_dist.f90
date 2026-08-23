module chernoff_dist
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
  use chernoff_kinds, only: dp
  use chernoff_data, only: n_small, n_airy, acoef, bcoef, airy_ai_zeros, airy_ai_deriv_at_zeros
  use chernoff_quad, only: integrate_adaptive
  implicit none
  private

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: two_third = 2.0_dp / 3.0_dp
  real(dp), parameter :: cbrt2 = 2.0_dp ** (1.0_dp / 3.0_dp)
  real(dp), parameter :: cbrt4 = 4.0_dp ** (1.0_dp / 3.0_dp)
  real(dp), parameter :: inner_upper = 12.0_dp

  public :: dchern, pchern, qchern
  public :: dchern_vector, pchern_vector, qchern_vector

contains

  function dchern(x, log_value) result(value)
    real(dp), intent(in) :: x
    logical, intent(in), optional :: log_value
    real(dp) :: value
    real(dp) :: gx, gmx
    logical :: want_log

    want_log = .false.
    if (present(log_value)) want_log = log_value

    gx = gfun(x)
    gmx = gfun(-x)
    value = 0.5_dp * gx * gmx
    if (value < 0.0_dp .and. abs(value) < 1.0e-14_dp) value = 0.0_dp

    if (want_log) then
      if (value > 0.0_dp) then
        value = log(value)
      else
        value = ieee_value(value, ieee_negative_inf)
      end if
    end if
  end function dchern

  function pchern(q, lower_tail, log_p) result(value)
    real(dp), intent(in) :: q
    logical, intent(in), optional :: lower_tail, log_p
    real(dp) :: value
    real(dp) :: area
    logical :: lower, want_log

    lower = .true.
    want_log = .false.
    if (present(lower_tail)) lower = lower_tail
    if (present(log_p)) want_log = log_p

    if (abs(q) <= tiny(1.0_dp)) then
      value = 0.5_dp
    else if (abs(q) >= 4.0_dp) then
      if (q > 0.0_dp) then
        value = 1.0_dp
      else
        value = 0.0_dp
      end if
    else
      if (abs(q) <= 1.0_dp) then
        area = integrate_adaptive(density_integrand, 0.0_dp, abs(q), 0.0_dp, 2.0e-10_dp, 2.0e-9_dp)
      else
        area = integrate_adaptive(density_integrand, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, 2.0e-9_dp)
        area = area + integrate_adaptive(density_integrand, 1.0_dp, abs(q), 0.0_dp, 1.0e-10_dp, 2.0e-9_dp)
      end if
      if (q > 0.0_dp) then
        value = min(1.0_dp, 0.5_dp + area)
      else
        value = max(0.0_dp, 0.5_dp - area)
      end if
    end if

    if (.not. lower) value = 1.0_dp - value
    if (want_log) then
      if (value > 0.0_dp) then
        value = log(value)
      else
        value = ieee_value(value, ieee_negative_inf)
      end if
    end if
  end function pchern

  function qchern(p, lower_tail, log_p) result(value)
    real(dp), intent(in) :: p
    logical, intent(in), optional :: lower_tail, log_p
    real(dp) :: value
    real(dp) :: prob, target, lo, hi, x, fx, den, trial
    logical :: lower, input_log
    integer :: iter

    lower = .true.
    input_log = .false.
    if (present(lower_tail)) lower = lower_tail
    if (present(log_p)) input_log = log_p

    if (input_log) then
      prob = exp(p)
    else
      prob = p
    end if
    if (.not. lower) prob = 1.0_dp - prob

    if (prob < 0.0_dp .or. prob > 1.0_dp) then
      value = ieee_value(prob, ieee_quiet_nan)
      return
    end if
    if (prob <= 0.0_dp) then
      value = ieee_value(prob, ieee_negative_inf)
      return
    end if
    if (prob >= 1.0_dp) then
      value = ieee_value(prob, ieee_positive_inf)
      return
    end if
    if (abs(prob - 0.5_dp) <= epsilon(prob)) then
      value = 0.0_dp
      return
    end if

    if (prob < 0.5_dp) then
      target = 1.0_dp - prob
    else
      target = prob
    end if

    lo = 0.0_dp
    hi = 2.0_dp
    do while (pchern(hi) < target .and. hi < 4.0_dp)
      hi = min(4.0_dp, 1.5_dp * hi)
    end do
    x = 0.5_dp * (lo + hi)

    do iter = 1, 30
      fx = pchern(x) - target
      if (abs(fx) < 2.0e-10_dp) exit
      if (fx > 0.0_dp) then
        hi = x
      else
        lo = x
      end if
      den = dchern(x)
      if (den > 1.0e-12_dp) then
        trial = x - fx / den
        if (trial > lo .and. trial < hi) then
          x = trial
        else
          x = 0.5_dp * (lo + hi)
        end if
      else
        x = 0.5_dp * (lo + hi)
      end if
    end do

    if (prob < 0.5_dp) then
      value = -x
    else
      value = x
    end if
  end function qchern

  subroutine dchern_vector(x, value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value(size(x))
    integer :: i

    do i = 1, size(x)
      value(i) = dchern(x(i))
    end do
  end subroutine dchern_vector

  subroutine pchern_vector(q, value)
    real(dp), intent(in) :: q(:)
    real(dp), intent(out) :: value(size(q))
    integer :: i

    do i = 1, size(q)
      value(i) = pchern(q(i))
    end do
  end subroutine pchern_vector

  subroutine qchern_vector(p, value)
    real(dp), intent(in) :: p(:)
    real(dp), intent(out) :: value(size(p))
    integer :: i

    do i = 1, size(p)
      value(i) = qchern(p(i))
    end do
  end subroutine qchern_vector

  function density_integrand(x, context) result(value)
    real(dp), intent(in) :: x, context
    real(dp) :: value

    value = dchern(x) + 0.0_dp * context
  end function density_integrand

  function gfun(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    if (x > -1.0_dp) then
      value = g2fun(x)
    else
      value = g3fun(x)
    end if
  end function gfun

  function g2fun(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    value = 2.0_dp * x - g1fun(x) + hfun(x)
  end function g2fun

  function g3fun(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    real(dp) :: s
    integer :: k

    s = 0.0_dp
    do k = 1, n_airy
      s = s + exp(-cbrt2 * airy_ai_zeros(k) * x) / airy_ai_deriv_at_zeros(k)
    end do
    value = cbrt4 * exp(two_third * x ** 3) * s
  end function g3fun

  function hfun(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    value = 2.0_dp * sqrt(2.0_dp / pi) * &
      integrate_adaptive(h_integrand, 0.0_dp, inner_upper, x, 2.0e-11_dp, 2.0e-10_dp)
  end function hfun

  function h_integrand(v, x) result(value)
    real(dp), intent(in) :: v, x
    real(dp) :: value
    real(dp) :: z

    z = 2.0_dp * x + v * v
    value = (z * v * v + 0.5_dp * z * z) * exp(-0.5_dp * v * v * z * z)
  end function h_integrand

  function g1fun(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    value = integrate_adaptive(g1_integrand, 0.0_dp, inner_upper, x, 2.0e-11_dp, 2.0e-10_dp) / sqrt(2.0_dp * pi)
  end function g1fun

  function g1_integrand(t, x) result(value)
    real(dp), intent(in) :: t, x
    real(dp) :: value

    value = pfun(t) * exp(-0.5_dp * t * (2.0_dp * x + t) ** 2)
  end function g1_integrand

  function pfun(t) result(value)
    real(dp), intent(in) :: t
    real(dp) :: value
    real(dp) :: sum_a, sum_b, p2
    integer :: k

    if (t <= 1.0_dp) then
      sum_a = 0.0_dp
      sum_b = 0.0_dp
      do k = 1, n_small
        sum_a = sum_a + acoef(k) * t ** (3 * (k - 1))
        sum_b = sum_b + bcoef(k) * t ** (3.0_dp * (real(k, dp) - 0.5_dp))
      end do
      value = -sqrt(pi / 2.0_dp) * sum_a + sum_b
    else
      p2 = 0.0_dp
      do k = 1, n_small
        p2 = p2 + exp(cbrt2 * t * airy_ai_zeros(k))
      end do
      value = -t ** (-1.5_dp) + 2.0_dp * sqrt(2.0_dp * pi) * exp(-t ** 3 / 6.0_dp) * p2
    end if
  end function pfun

end module chernoff_dist
