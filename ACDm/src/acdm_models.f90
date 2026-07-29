! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_models
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use acdm_kinds, only : dp, tiny_pos, huge_penalty, ACDM_SUCCESS, &
                         ACDM_BAD_INPUT, ACDM_BAD_PARAMETER, &
                         ACDM_NUMERIC_FAILURE
  use acdm_math, only : rng_state
  use acdm_distributions, only : distribution_logpdf, sample_distribution
  implicit none
  private

  integer, parameter, public :: MODEL_ACD = 1
  integer, parameter, public :: MODEL_LACD1 = 2
  integer, parameter, public :: MODEL_LACD2 = 3
  integer, parameter, public :: MODEL_EXACD = 4
  integer, parameter, public :: MODEL_AMACD = 5
  integer, parameter, public :: MODEL_ABACD = 6
  integer, parameter, public :: MODEL_AACD = 7
  integer, parameter, public :: MODEL_TACD = 8
  integer, parameter, public :: MODEL_BACD = 9
  integer, parameter, public :: MODEL_BCACD = 10
  integer, parameter, public :: MODEL_SNIACD = 11
  integer, parameter, public :: MODEL_LSNIACD = 12
  integer, parameter, public :: MODEL_TAMACD = 13

  type, public :: acd_order
    integer :: p = 1
    integer :: r = 0
    integer :: q = 1
  end type acd_order

  public :: model_code, model_name, model_parameter_count
  public :: default_model_parameters, filter_acd, acd_loglik
  public :: simulate_acd

contains

  pure integer function model_code(name) result(code)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: s

    s = uppercase(trim(adjustl(name)))
    select case (s)
    case ('ACD')
      code = MODEL_ACD
    case ('LACD1')
      code = MODEL_LACD1
    case ('LACD2')
      code = MODEL_LACD2
    case ('EXACD')
      code = MODEL_EXACD
    case ('AMACD')
      code = MODEL_AMACD
    case ('ABACD')
      code = MODEL_ABACD
    case ('AACD')
      code = MODEL_AACD
    case ('TACD')
      code = MODEL_TACD
    case ('BACD')
      code = MODEL_BACD
    case ('BCACD')
      code = MODEL_BCACD
    case ('SNIACD')
      code = MODEL_SNIACD
    case ('LSNIACD')
      code = MODEL_LSNIACD
    case ('TAMACD')
      code = MODEL_TAMACD
    case default
      code = 0
    end select
  end function model_code

  pure function model_name(code) result(name)
    integer, intent(in) :: code
    character(len=8) :: name

    select case (code)
    case (MODEL_ACD)
      name = 'ACD'
    case (MODEL_LACD1)
      name = 'LACD1'
    case (MODEL_LACD2)
      name = 'LACD2'
    case (MODEL_EXACD)
      name = 'EXACD'
    case (MODEL_AMACD)
      name = 'AMACD'
    case (MODEL_ABACD)
      name = 'ABACD'
    case (MODEL_AACD)
      name = 'AACD'
    case (MODEL_TACD)
      name = 'TACD'
    case (MODEL_BACD)
      name = 'BACD'
    case (MODEL_BCACD)
      name = 'BCACD'
    case (MODEL_SNIACD)
      name = 'SNIACD'
    case (MODEL_LSNIACD)
      name = 'LSNIACD'
    case (MODEL_TAMACD)
      name = 'TAMACD'
    case default
      name = 'UNKNOWN'
    end select
  end function model_name

  pure function uppercase(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i, c

    out = s
    do i = 1, len(s)
      c = iachar(out(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) then
        out(i:i) = achar(c - 32)
      end if
    end do
  end function uppercase

  pure integer function model_parameter_count(model, order, nbreak) result(n)
    integer, intent(in) :: model
    type(acd_order), intent(in) :: order
    integer, intent(in), optional :: nbreak
    integer :: j

    j = 1
    if (present(nbreak)) j = nbreak + 1
    select case (model)
    case (MODEL_ACD, MODEL_LACD1, MODEL_LACD2)
      n = 1 + order%p + order%q
    case (MODEL_EXACD)
      n = 1 + 2 * order%p + order%q
    case (MODEL_AMACD)
      n = 1 + order%p + order%r + order%q
    case (MODEL_ABACD, MODEL_AACD)
      n = 1 + order%p + order%q + 4
    case (MODEL_BACD)
      n = 1 + order%p + order%q + 2
    case (MODEL_BCACD)
      n = 1 + order%p + order%q + 1
    case (MODEL_TACD)
      n = j * (1 + order%p + order%q)
    case (MODEL_TAMACD)
      n = j * (1 + order%p + order%r + order%q)
    case (MODEL_SNIACD, MODEL_LSNIACD)
      n = 1 + j + max(0, order%p - 1) + order%q
    case default
      n = -1
    end select
  end function model_parameter_count

  function default_model_parameters(model, order, mean_x, nbreak) result(par)
    integer, intent(in) :: model
    type(acd_order), intent(in) :: order
    real(dp), intent(in) :: mean_x
    integer, intent(in), optional :: nbreak
    real(dp), allocatable :: par(:)
    integer :: n, j, p, r, q, pos, reg

    p = order%p
    r = order%r
    q = order%q
    j = 1
    if (present(nbreak)) j = nbreak + 1
    n = model_parameter_count(model, order, j - 1)
    allocate(par(n))
    par = 0.0_dp

    select case (model)
    case (MODEL_ACD)
      par(1) = mean_x / 10.0_dp
      if (p > 0) par(2:1 + p) = 0.15_dp / real(p, dp)
      if (q > 0) par(2 + p:1 + p + q) = 0.8_dp / real(q, dp)
    case (MODEL_LACD1)
      par(1) = 0.03_dp
      if (p > 0) par(2:1 + p) = 0.03_dp / real(p, dp)
      if (q > 0) par(2 + p:1 + p + q) = 0.95_dp / real(q, dp)
    case (MODEL_LACD2)
      par(1) = 0.0_dp
      if (p > 0) par(2:1 + p) = 0.03_dp / real(p, dp)
      if (q > 0) par(2 + p:1 + p + q) = 0.95_dp / real(q, dp)
    case (MODEL_EXACD)
      par(1) = 0.0_dp
      if (p > 0) then
        par(2:1 + p) = 0.05_dp / real(p, dp)
        par(2 + p:1 + 2 * p) = 0.01_dp
      end if
      if (q > 0) par(2 + 2 * p:1 + 2 * p + q) = 0.9_dp / real(q, dp)
    case (MODEL_AMACD)
      par(1) = mean_x / 20.0_dp
      if (p > 0) par(2:1 + p) = 0.10_dp / real(p, dp)
      if (r > 0) par(2 + p:1 + p + r) = 0.02_dp / real(r, dp)
      if (q > 0) par(2 + p + r:1 + p + r + q) = 0.8_dp / real(q, dp)
    case (MODEL_ABACD, MODEL_AACD)
      par(1) = mean_x / 20.0_dp
      if (p > 0) par(2:1 + p) = 0.10_dp / real(p, dp)
      if (q > 0) par(2 + p:1 + p + q) = 0.8_dp / real(q, dp)
      pos = 2 + p + q
      par(pos) = 0.8_dp
      par(pos + 1) = 0.1_dp
      par(pos + 2) = 1.1_dp
      par(pos + 3) = 1.1_dp
    case (MODEL_BACD)
      par(1) = mean_x / 20.0_dp
      if (p > 0) par(2:1 + p) = 0.10_dp / real(p, dp)
      if (q > 0) par(2 + p:1 + p + q) = 0.8_dp / real(q, dp)
      par(2 + p + q) = 1.0_dp
      par(3 + p + q) = 1.0_dp
    case (MODEL_BCACD)
      par(1) = 0.0_dp
      if (p > 0) par(2:1 + p) = 0.10_dp / real(p, dp)
      if (q > 0) par(2 + p:1 + p + q) = 0.8_dp / real(q, dp)
      par(2 + p + q) = 1.0_dp
    case (MODEL_TACD)
      pos = 1
      do reg = 1, j
        par(pos) = mean_x / 10.0_dp
        pos = pos + 1
      end do
      do reg = 1, j
        if (p > 0) par(pos:pos + p - 1) = 0.15_dp / real(p, dp)
        pos = pos + p
      end do
      do reg = 1, j
        if (q > 0) par(pos:pos + q - 1) = 0.8_dp / real(q, dp)
        pos = pos + q
      end do
    case (MODEL_TAMACD)
      pos = 1
      do reg = 1, j
        par(pos) = mean_x / 10.0_dp
        pos = pos + 1
      end do
      do reg = 1, j
        if (p > 0) par(pos:pos + p - 1) = 0.10_dp / real(p, dp)
        pos = pos + p
      end do
      do reg = 1, j
        if (r > 0) par(pos:pos + r - 1) = 0.02_dp / real(r, dp)
        pos = pos + r
      end do
      do reg = 1, j
        if (q > 0) par(pos:pos + q - 1) = 0.8_dp / real(q, dp)
        pos = pos + q
      end do
    case (MODEL_SNIACD, MODEL_LSNIACD)
      if (model == MODEL_SNIACD) then
        par(1) = mean_x / 10.0_dp
      else
        par(1) = 0.03_dp
      end if
      par(2) = 0.15_dp
      if (j > 1) par(3:1 + j) = 0.02_dp
      pos = 2 + j
      if (p > 1) then
        par(pos:pos + p - 2) = 0.1_dp
        pos = pos + p - 1
      end if
      if (q > 0) par(pos:pos + q - 1) = 0.8_dp / real(q, dp)
    end select
  end function default_model_parameters

  subroutine filter_acd(x, model, order, par, mu, residual, status, &
                        breakpoints, exogenous, new_day, threshold, &
                        threshold_on_mu)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: model
    type(acd_order), intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp), intent(out) :: mu(:), residual(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: breakpoints(:)
    real(dp), intent(in), optional :: exogenous(:, :)
    integer, intent(in), optional :: new_day(:)
    real(dp), intent(in), optional :: threshold(:)
    logical, intent(in), optional :: threshold_on_mu

    integer :: n, base_n, exo_n, maxlag, i, j, segment_start
    integer :: nd_pos, next_day, reg, jreg, pos, m
    real(dp) :: mean_x, value, z, c, v, d1, d2, holder, tvar
    real(dp), allocatable :: logmu(:), powered(:)
    logical :: use_mu_threshold

    status = ACDM_SUCCESS
    n = size(x)
    if (size(mu) /= n .or. size(residual) /= n .or. n < 1 .or. &
        any(x <= 0.0_dp)) then
      status = ACDM_BAD_INPUT
      return
    end if
    m = 0
    if (present(breakpoints)) m = size(breakpoints)
    base_n = model_parameter_count(model, order, m)
    if (base_n < 1 .or. size(par) < base_n) then
      status = ACDM_BAD_PARAMETER
      return
    end if
    exo_n = size(par) - base_n
    if (exo_n > 0) then
      if (.not. present(exogenous)) then
        status = ACDM_BAD_INPUT
        return
      end if
      if (size(exogenous, 1) /= n .or. size(exogenous, 2) /= exo_n) then
        status = ACDM_BAD_INPUT
        return
      end if
    end if
    if (present(threshold)) then
      if (size(threshold) /= n) then
        status = ACDM_BAD_INPUT
        return
      end if
    end if

    maxlag = max(order%p, order%q, order%r)
    if (maxlag < 1) maxlag = 1
    mean_x = sum(x) / real(n, dp)
    allocate(logmu(n), powered(n))
    logmu = log(mean_x)
    powered = mean_x
    use_mu_threshold = .false.
    if (present(threshold_on_mu)) use_mu_threshold = threshold_on_mu

    segment_start = 1
    nd_pos = 1
    next_day = n + 1
    if (present(new_day)) then
      if (size(new_day) > 0) next_day = new_day(1)
    end if

    do i = 1, n
      if (i == next_day) then
        segment_start = i
        nd_pos = nd_pos + 1
        next_day = n + 1
        if (present(new_day)) then
          if (nd_pos <= size(new_day)) next_day = new_day(nd_pos)
        end if
      end if

      if (i - segment_start < maxlag) then
        mu(i) = mean_x
        logmu(i) = log(mean_x)
        residual(i) = x(i) / mu(i)
        powered(i) = mean_x
        cycle
      end if

      select case (model)
      case (MODEL_ACD)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * x(i - j)
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * mu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        mu(i) = value

      case (MODEL_LACD1)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * log(max(tiny_pos, residual(i - j)))
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)

      case (MODEL_LACD2)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * residual(i - j)
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)

      case (MODEL_EXACD)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * residual(i - j) + &
                  par(1 + order%p + j) * abs(residual(i - j) - 1.0_dp)
        end do
        do j = 1, order%q
          value = value + par(1 + 2 * order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)

      case (MODEL_AMACD)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * x(i - j)
        end do
        do j = 1, order%r
          value = value + par(1 + order%p + j) * residual(i - j)
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + order%r + j) * mu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        mu(i) = value

      case (MODEL_ABACD, MODEL_AACD)
        c = par(2 + order%p + order%q)
        v = par(3 + order%p + order%q)
        d1 = par(4 + order%p + order%q)
        d2 = par(5 + order%p + order%q)
        value = par(1)
        do j = 1, order%p
          z = residual(i - j) - v
          holder = abs(z) + c * z
          if (holder < 0.0_dp) then
            status = ACDM_BAD_PARAMETER
            return
          end if
          if (model == MODEL_ABACD) then
            value = value + par(1 + j) * holder**d2
          else
            value = value + par(1 + j) * powered(i - j) * holder**d2
          end if
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * powered(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        if (value <= 0.0_dp .or. d1 <= 0.0_dp) then
          status = ACDM_BAD_PARAMETER
          return
        end if
        powered(i) = value
        mu(i) = value**(1.0_dp / d1)

      case (MODEL_BACD)
        d1 = par(2 + order%p + order%q)
        d2 = par(3 + order%p + order%q)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * residual(i - j)**d2
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * powered(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        if (value <= 0.0_dp .or. d1 <= 0.0_dp) then
          status = ACDM_BAD_PARAMETER
          return
        end if
        powered(i) = value
        mu(i) = value**(1.0_dp / d1)

      case (MODEL_BCACD)
        d1 = par(2 + order%p + order%q)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * residual(i - j)**d1
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)

      case (MODEL_TACD, MODEL_TAMACD)
        if (.not. present(breakpoints)) then
          status = ACDM_BAD_INPUT
          return
        end if
        if (use_mu_threshold) then
          tvar = mu(max(1, i - 1))
        else if (present(threshold)) then
          tvar = threshold(max(1, i - 1))
        else
          tvar = x(max(1, i - 1))
        end if
        reg = threshold_regime(tvar, breakpoints)
        jreg = m + 1
        value = par(reg)
        pos = jreg + (reg - 1) * order%p
        do j = 1, order%p
          value = value + par(pos + j) * x(i - j)
        end do
        pos = jreg * (1 + order%p)
        if (model == MODEL_TAMACD) then
          pos = pos + (reg - 1) * order%r
          do j = 1, order%r
            value = value + par(pos + j) * residual(i - j)
          end do
          pos = jreg * (1 + order%p + order%r)
        end if
        pos = pos + (reg - 1) * order%q
        do j = 1, order%q
          value = value + par(pos + j) * mu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        mu(i) = value

      case (MODEL_SNIACD, MODEL_LSNIACD)
        if (.not. present(breakpoints)) then
          status = ACDM_BAD_INPUT
          return
        end if
        value = par(1)
        do j = 1, order%p
          holder = spline_lag_value(residual(i - j), par(2:2 + m), &
                                    breakpoints)
          if (j > 1) then
            pos = 2 + (m + 1) + (j - 2)
            holder = holder * par(pos)
          end if
          value = value + holder
        end do
        pos = 2 + (m + 1) + max(0, order%p - 1) - 1
        do j = 1, order%q
          if (model == MODEL_SNIACD) then
            value = value + par(pos + j) * mu(i - j)
          else
            value = value + par(pos + j) * logmu(i - j)
          end if
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        if (model == MODEL_SNIACD) then
          mu(i) = value
        else
          logmu(i) = value
          mu(i) = exp(value)
        end if

      case default
        status = ACDM_BAD_INPUT
        return
      end select

      if (.not. ieee_is_finite(mu(i)) .or. mu(i) <= tiny_pos .or. &
          mu(i) > 1.0e100_dp) then
        status = ACDM_NUMERIC_FAILURE
        return
      end if
      residual(i) = x(i) / mu(i)
      if (.not. ieee_is_finite(residual(i)) .or. residual(i) <= 0.0_dp) then
        status = ACDM_NUMERIC_FAILURE
        return
      end if
      if (model /= MODEL_LACD1 .and. model /= MODEL_LACD2 .and. &
          model /= MODEL_EXACD .and. model /= MODEL_BCACD .and. &
          model /= MODEL_LSNIACD) logmu(i) = log(mu(i))
      if (model /= MODEL_ABACD .and. model /= MODEL_AACD .and. &
          model /= MODEL_BACD) powered(i) = mu(i)
    end do
  end subroutine filter_acd

  subroutine add_exogenous(value, par, base_n, exogenous, i)
    real(dp), intent(inout) :: value
    real(dp), intent(in) :: par(:)
    integer, intent(in) :: base_n, i
    real(dp), intent(in), optional :: exogenous(:, :)
    integer :: k

    if (.not. present(exogenous)) return
    do k = 1, size(exogenous, 2)
      value = value + par(base_n + k) * exogenous(i, k)
    end do
  end subroutine add_exogenous

  pure integer function threshold_regime(value, bp) result(reg)
    real(dp), intent(in) :: value
    real(dp), intent(in) :: bp(:)
    integer :: k

    reg = 1
    do k = 1, size(bp)
      if (value > bp(k)) reg = k + 1
    end do
  end function threshold_regime

  pure function spline_lag_value(e, c, bp) result(value)
    real(dp), intent(in) :: e, c(:), bp(:)
    real(dp) :: value
    integer :: k

    value = c(1) * e
    do k = 1, size(bp)
      if (e > bp(k)) value = value + c(k + 1) * (e - bp(k))
    end do
  end function spline_lag_value

  function acd_loglik(x, model, order, par, dist_code, dist_para, &
                      force_mean, mu, residual, status, breakpoints, &
                      exogenous, new_day, threshold, threshold_on_mu) result(ll)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: model, dist_code
    type(acd_order), intent(in) :: order
    real(dp), intent(in) :: par(:), dist_para(:)
    logical, intent(in) :: force_mean
    real(dp), intent(out) :: mu(:), residual(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: breakpoints(:)
    real(dp), intent(in), optional :: exogenous(:, :)
    integer, intent(in), optional :: new_day(:)
    real(dp), intent(in), optional :: threshold(:)
    logical, intent(in), optional :: threshold_on_mu
    real(dp) :: ll, logf
    integer :: i

    call filter_acd(x, model, order, par, mu, residual, status, &
                    breakpoints, exogenous, new_day, threshold, &
                    threshold_on_mu)
    if (status /= ACDM_SUCCESS) then
      ll = -huge_penalty
      return
    end if
    ll = 0.0_dp
    do i = 1, size(x)
      logf = distribution_logpdf(residual(i), dist_code, dist_para, force_mean)
      if (.not. ieee_is_finite(logf) .or. logf < -1.0e90_dp) then
        status = ACDM_BAD_PARAMETER
        ll = -huge_penalty
        return
      end if
      ll = ll + logf - log(mu(i))
    end do
  end function acd_loglik

  subroutine simulate_acd(n, model, order, par, dist_code, dist_para, &
                          durations, status, rng, burn, errors, force_mean, &
                          start_x, start_mu, breakpoints, exogenous)
    integer, intent(in) :: n, model, dist_code
    type(acd_order), intent(in) :: order
    real(dp), intent(in) :: par(:), dist_para(:)
    real(dp), intent(out) :: durations(:)
    integer, intent(out) :: status
    type(rng_state), intent(inout) :: rng
    integer, intent(in), optional :: burn
    real(dp), intent(in), optional :: errors(:)
    logical, intent(in), optional :: force_mean
    real(dp), intent(in), optional :: start_x(:), start_mu(:)
    real(dp), intent(in), optional :: breakpoints(:)
    real(dp), intent(in), optional :: exogenous(:, :)

    integer :: nb, nt, maxlag, i, j, base_n, m, reg, jreg, pos
    real(dp), allocatable :: x(:), mu(:), e(:), logmu(:), powered(:)
    real(dp) :: value, z, c, v, d1, d2, holder
    logical :: fm

    status = ACDM_SUCCESS
    nb = 50
    if (present(burn)) nb = max(0, burn)
    nt = n + nb
    if (size(durations) /= n .or. n < 1) then
      status = ACDM_BAD_INPUT
      return
    end if
    m = 0
    if (present(breakpoints)) m = size(breakpoints)
    base_n = model_parameter_count(model, order, m)
    if (size(par) < base_n) then
      status = ACDM_BAD_PARAMETER
      return
    end if
    if (present(exogenous)) then
      if (size(exogenous, 1) /= nt .or. &
          size(exogenous, 2) /= size(par) - base_n) then
        status = ACDM_BAD_INPUT
        return
      end if
    end if

    allocate(x(nt), mu(nt), e(nt), logmu(nt), powered(nt))
    fm = .true.
    if (present(force_mean)) fm = force_mean
    if (present(errors)) then
      if (size(errors) /= nt) then
        status = ACDM_BAD_INPUT
        return
      end if
      e = errors
    else
      do i = 1, nt
        e(i) = sample_distribution(rng, dist_code, dist_para, fm)
      end do
    end if
    if (any(e <= 0.0_dp)) then
      status = ACDM_BAD_PARAMETER
      return
    end if

    maxlag = max(order%p, order%q, order%r)
    if (maxlag < 1) maxlag = 1
    do i = 1, min(maxlag, nt)
      if (present(start_mu)) then
        mu(i) = start_mu(1 + modulo(i - 1, size(start_mu)))
      else
        mu(i) = 1.0_dp
      end if
      if (present(start_x)) then
        x(i) = start_x(1 + modulo(i - 1, size(start_x)))
        e(i) = x(i) / mu(i)
      else
        x(i) = mu(i) * e(i)
      end if
      logmu(i) = log(mu(i))
      powered(i) = mu(i)
    end do

    do i = maxlag + 1, nt
      select case (model)
      case (MODEL_ACD)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * x(i - j)
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * mu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        mu(i) = value
      case (MODEL_LACD1)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * log(e(i - j))
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)
      case (MODEL_LACD2)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * e(i - j)
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)
      case (MODEL_EXACD)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * e(i - j) + &
                  par(1 + order%p + j) * abs(e(i - j) - 1.0_dp)
        end do
        do j = 1, order%q
          value = value + par(1 + 2 * order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)
      case (MODEL_AMACD)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * x(i - j)
        end do
        do j = 1, order%r
          value = value + par(1 + order%p + j) * e(i - j)
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + order%r + j) * mu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        mu(i) = value
      case (MODEL_ABACD, MODEL_AACD)
        c = par(2 + order%p + order%q)
        v = par(3 + order%p + order%q)
        d1 = par(4 + order%p + order%q)
        d2 = par(5 + order%p + order%q)
        value = par(1)
        do j = 1, order%p
          z = e(i - j) - v
          holder = abs(z) + c * z
          if (model == MODEL_ABACD) then
            value = value + par(1 + j) * holder**d2
          else
            value = value + par(1 + j) * powered(i - j) * holder**d2
          end if
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * powered(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        powered(i) = value
        mu(i) = value**(1.0_dp / d1)
      case (MODEL_BACD)
        d1 = par(2 + order%p + order%q)
        d2 = par(3 + order%p + order%q)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * e(i - j)**d2
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * powered(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        powered(i) = value
        mu(i) = value**(1.0_dp / d1)
      case (MODEL_BCACD)
        d1 = par(2 + order%p + order%q)
        value = par(1)
        do j = 1, order%p
          value = value + par(1 + j) * e(i - j)**d1
        end do
        do j = 1, order%q
          value = value + par(1 + order%p + j) * logmu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        logmu(i) = value
        mu(i) = exp(value)
      case (MODEL_TACD, MODEL_TAMACD)
        if (.not. present(breakpoints)) then
          status = ACDM_BAD_INPUT
          return
        end if
        reg = threshold_regime(x(i - 1), breakpoints)
        jreg = m + 1
        value = par(reg)
        pos = jreg + (reg - 1) * order%p
        do j = 1, order%p
          value = value + par(pos + j) * x(i - j)
        end do
        pos = jreg * (1 + order%p)
        if (model == MODEL_TAMACD) then
          pos = pos + (reg - 1) * order%r
          do j = 1, order%r
            value = value + par(pos + j) * e(i - j)
          end do
          pos = jreg * (1 + order%p + order%r)
        end if
        pos = pos + (reg - 1) * order%q
        do j = 1, order%q
          value = value + par(pos + j) * mu(i - j)
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        mu(i) = value
      case (MODEL_SNIACD, MODEL_LSNIACD)
        if (.not. present(breakpoints)) then
          status = ACDM_BAD_INPUT
          return
        end if
        value = par(1)
        do j = 1, order%p
          holder = spline_lag_value(e(i - j), par(2:2 + m), breakpoints)
          if (j > 1) then
            pos = 2 + (m + 1) + (j - 2)
            holder = holder * par(pos)
          end if
          value = value + holder
        end do
        pos = 2 + (m + 1) + max(0, order%p - 1) - 1
        do j = 1, order%q
          if (model == MODEL_SNIACD) then
            value = value + par(pos + j) * mu(i - j)
          else
            value = value + par(pos + j) * logmu(i - j)
          end if
        end do
        call add_exogenous(value, par, base_n, exogenous, i)
        if (model == MODEL_SNIACD) then
          mu(i) = value
        else
          logmu(i) = value
          mu(i) = exp(value)
        end if
      case default
        status = ACDM_BAD_INPUT
        return
      end select

      if (.not. ieee_is_finite(mu(i)) .or. mu(i) <= tiny_pos) then
        status = ACDM_NUMERIC_FAILURE
        return
      end if
      x(i) = mu(i) * e(i)
      if (model /= MODEL_LACD1 .and. model /= MODEL_LACD2 .and. &
          model /= MODEL_EXACD .and. model /= MODEL_BCACD .and. &
          model /= MODEL_LSNIACD) logmu(i) = log(mu(i))
      if (model /= MODEL_ABACD .and. model /= MODEL_AACD .and. &
          model /= MODEL_BACD) powered(i) = mu(i)
    end do
    durations = x(nb + 1:nt)
  end subroutine simulate_acd

end module acdm_models
