! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_mortality
  use r4gpf_kinds, only: dp
  use r4gpf_status, only: r4gpf_success, r4gpf_invalid_argument, r4gpf_numerical_error
  use r4gpf_optimization, only: golden_section_minimize, bisection_root, nelder_mead_minimize
  implicit none
  private
  public :: gompertz_fit, gompertz_survival_probability, gompertz_pdf, life_expectancy
  public :: gompertz_mode_from_life_expectancy, fit_gompertz_mortality
  public :: fit_joint_gompertz, upper_incomplete_gamma, retirement_ruin_probability
  public :: gompertz_annuity_factor, regularized_gamma_p

  type :: gompertz_fit
    real(dp) :: mode = 0.0_dp
    real(dp) :: dispersion = 0.0_dp
    real(dp) :: current_age = 0.0_dp
    real(dp) :: max_age = -1.0_dp
    real(dp) :: objective = huge(1.0_dp)
    integer :: status = r4gpf_success
    real(dp), allocatable :: ages(:)
    real(dp), allocatable :: survival(:)
    real(dp), allocatable :: fitted_survival(:)
    real(dp), allocatable :: probability_of_death(:)
  end type gompertz_fit
contains

  elemental real(dp) function gompertz_survival_probability(current_age, target_age, mode, dispersion, max_age) result(probability)
    real(dp), intent(in) :: current_age, target_age, mode, dispersion
    real(dp), intent(in), optional :: max_age
    real(dp) :: s_target, s_max, denominator

    if (dispersion <= 0.0_dp .or. target_age < current_age) then
      probability = 0.0_dp
      return
    end if
    s_target = exp(exp((current_age - mode) / dispersion) * &
      (1.0_dp - exp((target_age - current_age) / dispersion)))
    if (.not. present(max_age)) then
      probability = min(1.0_dp, max(0.0_dp, s_target))
      return
    end if
    if (target_age > max_age) then
      probability = 0.0_dp
      return
    end if
    s_max = exp(exp((current_age - mode) / dispersion) * &
      (1.0_dp - exp((max_age - current_age) / dispersion)))
    denominator = 1.0_dp - s_max
    if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp)) then
      probability = 0.0_dp
    else
      probability = (s_target - s_max) / denominator
      if (.not. (probability >= 0.0_dp .and. probability <= 1.0_dp)) probability = 0.0_dp
    end if
  end function gompertz_survival_probability

  elemental real(dp) function gompertz_pdf(current_age, target_age, mode, dispersion) result(value)
    real(dp), intent(in) :: current_age, target_age, mode, dispersion
    real(dp) :: hazard
    if (dispersion <= 0.0_dp .or. target_age < current_age) then
      value = 0.0_dp
      return
    end if
    hazard = exp((target_age - mode) / dispersion) / dispersion
    value = hazard * gompertz_survival_probability(current_age, target_age, mode, dispersion)
  end function gompertz_pdf

  real(dp) function life_expectancy(current_age, mode, dispersion, max_age) result(value)
    real(dp), intent(in) :: current_age, mode, dispersion
    real(dp), intent(in), optional :: max_age
    integer :: age, maxage_i
    real(dp) :: m_age

    m_age = 120.0_dp
    if (present(max_age)) m_age = max_age
    maxage_i = floor(m_age)
    value = current_age + 0.5_dp
    do age = floor(current_age) + 1, maxage_i
      value = value + gompertz_survival_probability(current_age, real(age, dp), mode, dispersion, m_age)
    end do
  end function life_expectancy

  subroutine gompertz_mode_from_life_expectancy(desired_life_expectancy, current_age, dispersion, mode, status, max_age)
    real(dp), intent(in) :: desired_life_expectancy, current_age, dispersion
    real(dp), intent(out) :: mode
    integer, intent(out) :: status
    real(dp), intent(in), optional :: max_age
    real(dp) :: m_age

    m_age = 120.0_dp
    if (present(max_age)) m_age = max_age
    call bisection_root(objective, current_age + 1.0_dp, 150.0_dp, mode, status, 1.0e-8_dp, 500)
  contains
    function objective(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v
      v = life_expectancy(current_age, x, dispersion, m_age) - desired_life_expectancy
    end function objective
  end subroutine gompertz_mode_from_life_expectancy

  subroutine fit_gompertz_mortality(ages, mortality_rates, current_age, fit, estimate_max_age)
    real(dp), intent(in) :: ages(:), mortality_rates(:), current_age
    type(gompertz_fit), intent(out) :: fit
    logical, intent(in), optional :: estimate_max_age
    real(dp), allocatable :: a(:), q(:), xbest(:)
    real(dp) :: fbest, dispersion
    integer :: i, n, first, st
    logical :: estimate_max

    fit%status = r4gpf_invalid_argument
    if (size(ages) /= size(mortality_rates) .or. size(ages) < 2) return
    first = 0
    do i = 1, size(ages)
      if (ages(i) >= current_age) then
        first = i
        exit
      end if
    end do
    if (first == 0) return
    n = size(ages) - first + 1
    allocate(a(n), q(n), fit%ages(n), fit%survival(n), fit%probability_of_death(n), fit%fitted_survival(n))
    a = ages(first:)
    q = min(1.0_dp, max(0.0_dp, mortality_rates(first:)))
    fit%ages = a
    fit%survival(1) = 1.0_dp
    do i = 2, n
      fit%survival(i) = fit%survival(i - 1) * (1.0_dp - q(i))
    end do
    fit%probability_of_death(1) = 0.0_dp
    do i = 2, n
      fit%probability_of_death(i) = fit%survival(i - 1) - fit%survival(i)
    end do
    fit%mode = a(maxloc(fit%probability_of_death, dim=1))
    fit%current_age = current_age
    estimate_max = .false.
    if (present(estimate_max_age)) estimate_max = estimate_max_age
    if (.not. estimate_max) then
      call golden_section_minimize(dispersion_objective, 0.05_dp, 100.0_dp, dispersion, fbest, st, 1.0e-9_dp, 1000)
      fit%dispersion = dispersion
      fit%max_age = -1.0_dp
    else
      call nelder_mead_minimize(full_objective, [10.0_dp, max(100.0_dp, maxval(a))], xbest, fbest, st, &
        step=[1.0_dp, 2.0_dp], lower=[0.05_dp, max(current_age + 1.0_dp, maxval(a) - 20.0_dp)], &
        upper=[100.0_dp, 160.0_dp], tolerance=1.0e-9_dp, max_iterations=5000)
      fit%dispersion = xbest(1)
      fit%max_age = xbest(2)
    end if
    fit%objective = fbest
    fit%status = st
    do i = 1, n
      if (fit%max_age > 0.0_dp) then
        fit%fitted_survival(i) = gompertz_survival_probability(current_age, a(i), fit%mode, fit%dispersion, fit%max_age)
      else
        fit%fitted_survival(i) = gompertz_survival_probability(current_age, a(i), fit%mode, fit%dispersion)
      end if
    end do
  contains
    function dispersion_objective(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v
      integer :: k
      v = 0.0_dp
      do k = 1, n
        v = v + (gompertz_survival_probability(current_age, a(k), fit%mode, x) - fit%survival(k))**2
      end do
    end function dispersion_objective

    function full_objective(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v
      integer :: k
      v = 0.0_dp
      do k = 1, n
        v = v + (gompertz_survival_probability(current_age, a(k), fit%mode, x(1), x(2)) - fit%survival(k))**2
      end do
    end function full_objective
  end subroutine fit_gompertz_mortality

  subroutine fit_joint_gompertz(age, mode, dispersion, max_age, fit, member_survival, joint_survival)
    real(dp), intent(in) :: age(:), mode(:), dispersion(:), max_age
    type(gompertz_fit), intent(out) :: fit
    real(dp), allocatable, intent(out), optional :: member_survival(:, :), joint_survival(:)
    real(dp), allocatable :: surv(:, :), joint(:), xbest(:)
    real(dp) :: fbest, min_age
    integer :: i, j, n_years, st

    if (size(age) /= size(mode) .or. size(age) /= size(dispersion) .or. size(age) < 1) then
      fit%status = r4gpf_invalid_argument
      return
    end if
    min_age = minval(age)
    n_years = max(1, floor(max_age - min_age) + 1)
    allocate(surv(size(age), n_years), joint(n_years))
    do j = 1, n_years
      do i = 1, size(age)
        surv(i, j) = gompertz_survival_probability(age(i), age(i) + real(j - 1, dp), mode(i), dispersion(i))
      end do
      joint(j) = 1.0_dp - product(1.0_dp - surv(:, j))
    end do
    call nelder_mead_minimize(objective, [sum(mode)/real(size(mode), dp), sum(dispersion)/real(size(dispersion), dp)], &
      xbest, fbest, st, step=[1.0_dp, 0.5_dp], lower=[min_age, 0.05_dp], upper=[160.0_dp, 100.0_dp], &
      tolerance=1.0e-10_dp, max_iterations=5000)
    fit%mode = xbest(1)
    fit%dispersion = xbest(2)
    fit%current_age = min_age
    fit%max_age = max_age
    fit%objective = fbest
    fit%status = st
    allocate(fit%ages(n_years), fit%survival(n_years), fit%fitted_survival(n_years), fit%probability_of_death(n_years))
    do j = 1, n_years
      fit%ages(j) = min_age + real(j - 1, dp)
      fit%survival(j) = joint(j)
      fit%fitted_survival(j) = gompertz_survival_probability(min_age, fit%ages(j), fit%mode, fit%dispersion)
    end do
    fit%probability_of_death(1) = 0.0_dp
    do j = 2, n_years
      fit%probability_of_death(j) = joint(j - 1) - joint(j)
    end do
    if (present(member_survival)) then
      allocate(member_survival(size(age), n_years))
      member_survival = surv
    end if
    if (present(joint_survival)) then
      allocate(joint_survival(n_years))
      joint_survival = joint
    end if
  contains
    function objective(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v
      integer :: k
      v = 0.0_dp
      do k = 1, n_years
        v = v + (gompertz_survival_probability(min_age, min_age + real(k - 1, dp), x(1), x(2)) - joint(k))**2
      end do
    end function objective
  end subroutine fit_joint_gompertz

  real(dp) function upper_incomplete_gamma(a, c, status) result(value)
    real(dp), intent(in) :: a, c
    integer, intent(out), optional :: status
    real(dp) :: upper, fa, fb, fm, whole

    if (c < 0.0_dp .or. (c <= tiny(1.0_dp) .and. a <= 0.0_dp)) then
      value = huge(1.0_dp)
      if (present(status)) status = r4gpf_invalid_argument
      return
    end if
    if (c <= tiny(1.0_dp)) then
      value = gamma(a)
      if (present(status)) status = r4gpf_success
      return
    end if
    upper = 1.0_dp - 1.0e-10_dp
    fa = transformed_integrand(0.0_dp)
    fb = transformed_integrand(upper)
    fm = transformed_integrand(0.5_dp * upper)
    whole = upper * (fa + 4.0_dp * fm + fb) / 6.0_dp
    value = adaptive_simpson(0.0_dp, upper, fa, fm, fb, whole, 1.0e-11_dp, 24)
    if (present(status)) then
      if (value >= 0.0_dp .and. value < huge(1.0_dp)) then
        status = r4gpf_success
      else
        status = r4gpf_numerical_error
      end if
    end if
  contains
    function transformed_integrand(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y, t, one_minus
      one_minus = 1.0_dp - x
      if (one_minus <= 0.0_dp) then
        y = 0.0_dp
        return
      end if
      t = c + x / one_minus
      if (t > 745.0_dp) then
        y = 0.0_dp
      else
        y = exp((a - 1.0_dp) * log(t) - t) / (one_minus * one_minus)
      end if
    end function transformed_integrand

    recursive function adaptive_simpson(left, right, fleft, fmid, fright, total, tol, depth) result(ans)
      real(dp), intent(in) :: left, right, fleft, fmid, fright, total, tol
      integer, intent(in) :: depth
      real(dp) :: ans, mid, lmid, rmid, flmid, frmid, left_s, right_s, delta
      mid = 0.5_dp * (left + right)
      lmid = 0.5_dp * (left + mid)
      rmid = 0.5_dp * (mid + right)
      flmid = transformed_integrand(lmid)
      frmid = transformed_integrand(rmid)
      left_s = (mid - left) * (fleft + 4.0_dp * flmid + fmid) / 6.0_dp
      right_s = (right - mid) * (fmid + 4.0_dp * frmid + fright) / 6.0_dp
      delta = left_s + right_s - total
      if (depth <= 0 .or. abs(delta) <= 15.0_dp * tol) then
        ans = left_s + right_s + delta / 15.0_dp
      else
        ans = adaptive_simpson(left, mid, fleft, flmid, fmid, left_s, 0.5_dp * tol, depth - 1) + &
          adaptive_simpson(mid, right, fmid, frmid, fright, right_s, 0.5_dp * tol, depth - 1)
      end if
    end function adaptive_simpson
  end function upper_incomplete_gamma

  real(dp) function gompertz_annuity_factor(v, age, mode, dispersion) result(value)
    real(dp), intent(in) :: v, age, mode, dispersion
    value = dispersion * exp(exp((age - mode) / dispersion) + (age - mode) * v) * &
      upper_incomplete_gamma(-dispersion * v, exp((age - mode) / dispersion))
  end function gompertz_annuity_factor

  real(dp) function regularized_gamma_p(shape, x) result(value)
    real(dp), intent(in) :: shape, x
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    real(dp) :: ap, del, sumv, b, c, d, h, an, gln
    integer :: n

    if (shape <= 0.0_dp .or. x < 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (x <= tiny(1.0_dp)) then
      value = 0.0_dp
      return
    end if
    gln = log_gamma(shape)
    if (x < shape + 1.0_dp) then
      ap = shape
      sumv = 1.0_dp / shape
      del = sumv
      do n = 1, 10000
        ap = ap + 1.0_dp
        del = del * x / ap
        sumv = sumv + del
        if (abs(del) < abs(sumv) * eps) exit
      end do
      value = sumv * exp(-x + shape * log(x) - gln)
    else
      b = x + 1.0_dp - shape
      c = 1.0_dp / fpmin
      d = 1.0_dp / b
      h = d
      do n = 1, 10000
        an = -real(n, dp) * (real(n, dp) - shape)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an / c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp / d
        del = d * c
        h = h * del
        if (abs(del - 1.0_dp) < eps) exit
      end do
      value = 1.0_dp - exp(-x + shape * log(x) - gln) * h
    end if
    value = min(1.0_dp, max(0.0_dp, value))
  end function regularized_gamma_p

  real(dp) function retirement_ruin_probability(portfolio_return_mean, portfolio_return_sd, age, mode, dispersion, &
      spending_rate) result(probability)
    real(dp), intent(in) :: portfolio_return_mean, portfolio_return_sd, age, mode, dispersion, spending_rate
    real(dp) :: nu, sigma, mu, m1, m2, alpha, beta, denom

    nu = portfolio_return_mean
    sigma = portfolio_return_sd
    mu = nu + 0.5_dp * sigma**2
    m1 = gompertz_annuity_factor(mu - sigma**2, age, mode, dispersion)
    denom = mu / 2.0_dp - sigma**2
    if (abs(denom) <= 100.0_dp * epsilon(1.0_dp)) then
      probability = 0.0_dp
      return
    end if
    m2 = (m1 - gompertz_annuity_factor(2.0_dp * mu - 3.0_dp * sigma**2, age, mode, dispersion)) / denom
    alpha = (2.0_dp * m2 - m1**2) / (m2 - m1**2)
    beta = (m2 - m1**2) / (m2 * m1)
    probability = regularized_gamma_p(alpha, spending_rate / beta)
  end function retirement_ruin_probability

end module r4gpf_mortality
