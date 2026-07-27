! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
! Numerical translation of strand 0.2.3.
module strand_stats
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use strand_kinds, only : dp
  use strand_linalg, only : least_squares, weighted_least_squares
  use strand_types, only : performance_stats, exposure_result
  implicit none
  private
  public :: rank_normal, normalize, normalize_grouped, adjust_numeric
  public :: maximum_drawdown, factor_exposure, category_exposure, calculate_exposures
  public :: summarize_performance

contains

  function rank_normal(x, out_sd) result(z)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: out_sd
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: work(:), q(:), ranks(:)
    integer, allocatable :: original(:), order(:)
    real(dp) :: mu, sd, target_sd, qmu, qsd, rounded
    integer :: i, j, k, nvalid, first, last

    target_sd = 1.0_dp
    if (present(out_sd)) target_sd = out_sd
    allocate(z(size(x)))
    z = ieee_value(0.0_dp, ieee_quiet_nan)
    nvalid = count(ieee_is_finite(x))
    if (nvalid == 0) return
    allocate(work(nvalid), q(nvalid), ranks(nvalid), original(nvalid), order(nvalid))
    j = 0
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        j = j + 1
        work(j) = x(i)
        original(j) = i
        order(j) = j
      end if
    end do
    if (nvalid == 1) then
      z(original(1)) = 0.0_dp
      return
    end if

    mu = sum(work) / real(nvalid, dp)
    sd = sample_sd(work)
    if (sd <= 100.0_dp * epsilon(1.0_dp)) then
      do i = 1, nvalid
        z(original(i)) = 0.0_dp
      end do
      return
    end if
    do i = 1, nvalid
      rounded = anint(((work(i) - mu) / sd) * 1.0e11_dp) / 1.0e11_dp
      work(i) = rounded
    end do
    call sort_indices(work, order)
    first = 1
    do while (first <= nvalid)
      last = first
      do while (last < nvalid)
        if (abs(work(order(last + 1)) - work(order(first))) > 0.0_dp) exit
        last = last + 1
      end do
      do k = first, last
        ranks(order(k)) = 0.5_dp * real(first + last, dp)
      end do
      first = last + 1
    end do
    do i = 1, nvalid
      q(i) = inverse_normal_cdf(ranks(i) / real(nvalid + 1, dp))
    end do
    qmu = sum(q) / real(nvalid, dp)
    qsd = sample_sd(q)
    if (qsd <= tiny(1.0_dp)) then
      q = 0.0_dp
    else
      q = (q - qmu) / qsd * target_sd
    end if
    do i = 1, nvalid
      z(original(i)) = q(i)
    end do
  end function rank_normal

  function normalize(x) result(z)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: z(:)
    z = rank_normal(x)
  end function normalize

  function normalize_grouped(x, within_groups, loops, scale_group, scale_level, &
    scale_value, adjustment) result(z)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: within_groups(:, :)
    integer, intent(in), optional :: loops
    integer, intent(in), optional :: scale_group(:), scale_level(:)
    real(dp), intent(in), optional :: scale_value(:)
    real(dp), intent(in), optional :: adjustment(:, :)
    real(dp), allocatable :: z(:), subset(:), adjusted(:)
    logical, allocatable :: mask(:)
    integer :: nloop, iter, g, level, i
    real(dp) :: current_sd, factor

    z = rank_normal(x)
    nloop = 3
    if (present(loops)) nloop = max(1, loops)
    allocate(mask(size(x)))
    do iter = 1, nloop
      if (present(within_groups)) then
        if (size(within_groups, 1) /= size(x)) error stop 'normalize_grouped: row mismatch'
        do g = 1, size(within_groups, 2)
          do level = minval(within_groups(:, g)), maxval(within_groups(:, g))
            mask = within_groups(:, g) == level
            if (count(mask .and. ieee_is_finite(z)) == 0) cycle
            allocate(subset(count(mask)))
            subset = pack(z, mask)
            subset = rank_normal(subset)
            z = unpack(subset, mask, z)
            deallocate(subset)
          end do
          call apply_requested_scales(z, within_groups, g, scale_group, scale_level, scale_value)
        end do
      end if
      z = rank_normal(z)

      if (present(scale_group) .and. present(scale_level) .and. present(scale_value)) then
        if (.not. present(within_groups)) error stop 'normalize_grouped: scales require groups'
        do i = 1, size(scale_group)
          g = scale_group(i)
          level = scale_level(i)
          if (g < 1 .or. g > size(within_groups, 2)) cycle
          mask = within_groups(:, g) == level .and. ieee_is_finite(z)
          if (count(mask) < 2) cycle
          current_sd = sample_sd(pack(z, mask))
          if (current_sd > tiny(1.0_dp)) then
            factor = scale_value(i) / current_sd
            where (mask) z = z * factor
          end if
        end do
        z = rank_normal(z)
      end if

      if (present(adjustment)) then
        adjusted = adjust_numeric(z, adjustment, 5)
        z = rank_normal(adjusted)
      end if
    end do
  end function normalize_grouped

  subroutine apply_requested_scales(z, groups, group_index, scale_group, scale_level, scale_value)
    real(dp), intent(inout) :: z(:)
    integer, intent(in) :: groups(:, :), group_index
    integer, intent(in), optional :: scale_group(:), scale_level(:)
    real(dp), intent(in), optional :: scale_value(:)
    integer :: i

    if (.not. present(scale_group) .or. .not. present(scale_level) .or. .not. present(scale_value)) return
    do i = 1, size(scale_group)
      if (scale_group(i) == group_index) then
        where (groups(:, group_index) == scale_level(i)) z = z * scale_value(i)
      end if
    end do
  end subroutine apply_requested_scales

  function adjust_numeric(y, x, loops) result(adjusted)
    real(dp), intent(in) :: y(:), x(:, :)
    integer, intent(in), optional :: loops
    real(dp), allocatable :: adjusted(:), yn(:), xn(:, :), design(:, :), beta(:), residual(:), weights(:)
    logical, allocatable :: valid(:)
    logical :: ok
    integer :: nloop, iter, i, j, nvalid, p

    if (size(x, 1) /= size(y)) error stop 'adjust_numeric: row mismatch'
    allocate(adjusted(size(y)), valid(size(y)))
    adjusted = ieee_value(0.0_dp, ieee_quiet_nan)
    valid = ieee_is_finite(y)
    do j = 1, size(x, 2)
      valid = valid .and. ieee_is_finite(x(:, j))
    end do
    nvalid = count(valid)
    p = size(x, 2)
    if (nvalid < 2) return
    allocate(yn(nvalid), xn(nvalid, p), design(nvalid, p + 1), beta(p + 1), residual(nvalid), weights(nvalid))
    yn = pack(y, valid)
    yn = rank_normal(yn)
    do j = 1, p
      xn(:, j) = rank_normal(pack(x(:, j), valid))
    end do
    design(:, 1) = 1.0_dp
    if (p > 0) design(:, 2:p + 1) = xn
    nloop = 0
    if (present(loops)) nloop = max(0, loops)
    do iter = 1, nloop
      do i = 1, nvalid
        weights(i) = 1.0_dp - 1.0_dp / (1.0_dp + exp(4.8_dp * (abs(yn(i)) - 1.5_dp)))
      end do
      call weighted_least_squares(design, yn, weights, beta, ok)
      if (ok) then
        residual = yn - matmul(design, beta)
        yn = rank_normal(residual)
      end if
      call least_squares(design, yn, beta, ok)
      if (ok) then
        residual = yn - matmul(design, beta)
        yn = rank_normal(residual)
      end if
    end do
    call least_squares(design, yn, beta, ok)
    if (ok) then
      residual = yn - matmul(design, beta)
      yn = rank_normal(residual)
    end if
    adjusted = unpack(yn, valid, adjusted)
  end function adjust_numeric

  function maximum_drawdown(returns) result(dd)
    real(dp), intent(in) :: returns(:)
    real(dp) :: dd, cumulative, peak
    integer :: i

    cumulative = 0.0_dp
    peak = 0.0_dp
    dd = 0.0_dp
    do i = 1, size(returns)
      cumulative = cumulative + returns(i)
      peak = max(peak, cumulative)
      dd = min(dd, cumulative - peak)
    end do
  end function maximum_drawdown

  function factor_exposure(nmv, factor, divisor) result(exposure)
    real(dp), intent(in) :: nmv(:), factor(:), divisor
    real(dp) :: exposure
    if (size(nmv) /= size(factor) .or. abs(divisor) <= tiny(1.0_dp)) error stop 'factor_exposure: invalid input'
    exposure = dot_product(nmv, factor) / divisor
  end function factor_exposure

  function category_exposure(nmv, category, level, divisor) result(exposure)
    real(dp), intent(in) :: nmv(:), divisor
    integer, intent(in) :: category(:), level
    real(dp) :: exposure
    if (size(nmv) /= size(category) .or. abs(divisor) <= tiny(1.0_dp)) error stop 'category_exposure: invalid input'
    exposure = sum(nmv, mask=category == level) / divisor
  end function category_exposure


  function calculate_exposures(nmv, capital, factors, categories) result(exposures)
    real(dp), intent(in) :: nmv(:, :), capital(:)
    real(dp), intent(in), optional :: factors(:, :)
    integer, intent(in), optional :: categories(:, :)
    type(exposure_result) :: exposures
    integer :: n, s, nf, nc, max_count, i, j, k, count_k, level
    integer, allocatable :: level_work(:, :)

    n = size(nmv, 1)
    s = size(nmv, 2)
    if (size(capital) /= s .or. any(abs(capital) <= tiny(1.0_dp))) then
      error stop 'calculate_exposures: invalid capital'
    end if
    nf = 0
    if (present(factors)) then
      if (size(factors, 1) /= n) error stop 'calculate_exposures: factor row mismatch'
      nf = size(factors, 2)
    end if
    allocate(exposures%factor(nf, s))
    exposures%factor = 0.0_dp
    if (present(factors)) then
      do j = 1, s
        exposures%factor(:, j) = matmul(transpose(factors), nmv(:, j)) / capital(j)
      end do
    end if

    nc = 0
    max_count = 0
    if (present(categories)) then
      if (size(categories, 1) /= n) error stop 'calculate_exposures: category row mismatch'
      nc = size(categories, 2)
      allocate(exposures%category_count(nc), level_work(n, nc))
      exposures%category_count = 0
      level_work = 0
      do k = 1, nc
        count_k = 0
        do i = 1, n
          if (count_k == 0 .or. .not. any(level_work(1:count_k, k) == categories(i, k))) then
            count_k = count_k + 1
            level_work(count_k, k) = categories(i, k)
          end if
        end do
        exposures%category_count(k) = count_k
        max_count = max(max_count, count_k)
      end do
      allocate(exposures%category_level(max_count, nc))
      exposures%category_level = 0
      do k = 1, nc
        exposures%category_level(1:exposures%category_count(k), k) = &
          level_work(1:exposures%category_count(k), k)
      end do
    else
      allocate(exposures%category_count(0), exposures%category_level(0, 0))
    end if
    allocate(exposures%category(max_count, nc, s))
    exposures%category = 0.0_dp
    if (present(categories)) then
      do j = 1, s
        do k = 1, nc
          do i = 1, exposures%category_count(k)
            level = exposures%category_level(i, k)
            exposures%category(i, k, j) = sum(nmv(:, j), mask=categories(:, k) == level) / capital(j)
          end do
        end do
      end do
    end if
  end function calculate_exposures

  function summarize_performance(pnl, gmv, nmv, turnover, periods_per_year) result(stats)
    real(dp), intent(in) :: pnl(:), gmv(:), nmv(:), turnover(:)
    real(dp), intent(in), optional :: periods_per_year
    type(performance_stats) :: stats
    real(dp), allocatable :: ret(:)
    real(dp) :: annual, vol
    logical, allocatable :: valid(:)

    if (size(pnl) /= size(gmv) .or. size(nmv) /= size(gmv) .or. size(turnover) /= size(gmv)) then
      error stop 'summarize_performance: dimension mismatch'
    end if
    annual = 252.0_dp
    if (present(periods_per_year)) annual = periods_per_year
    allocate(valid(size(pnl)), ret(size(pnl)))
    valid = gmv > tiny(1.0_dp) .and. ieee_is_finite(pnl) .and. ieee_is_finite(gmv)
    ret = 0.0_dp
    where (valid) ret = pnl / gmv
    stats%total_pnl = sum(pnl)
    if (count(valid) > 0) then
      stats%total_return = sum(ret, mask=valid)
      stats%annualized_return = sum(ret, mask=valid) / real(count(valid), dp) * annual
      if (count(valid) > 1) then
        vol = sample_sd(pack(ret, valid)) * sqrt(annual)
      else
        vol = 0.0_dp
      end if
      stats%annualized_volatility = vol
      if (vol > tiny(1.0_dp)) stats%annualized_sharpe = stats%annualized_return / vol
      stats%max_drawdown = maximum_drawdown(pack(ret, valid))
    end if
    if (size(gmv) > 0) then
      stats%average_gmv = sum(gmv) / real(size(gmv), dp)
      stats%average_nmv = sum(nmv) / real(size(nmv), dp)
      stats%average_turnover = sum(turnover) / real(size(turnover), dp)
      if (stats%average_turnover > tiny(1.0_dp)) then
        stats%holding_period_months = 12.0_dp / (stats%average_turnover / max(stats%average_gmv, tiny(1.0_dp)) * annual / 2.0_dp)
      end if
    end if
  end function summarize_performance

  real(dp) function sample_sd(x) result(sd)
    real(dp), intent(in) :: x(:)
    real(dp) :: mu
    if (size(x) < 2) then
      sd = 0.0_dp
      return
    end if
    mu = sum(x) / real(size(x), dp)
    sd = sqrt(sum((x - mu)**2) / real(size(x) - 1, dp))
  end function sample_sd

  subroutine sort_indices(values, order)
    real(dp), intent(in) :: values(:)
    integer, intent(inout) :: order(:)
    integer :: i, j, key
    do i = 2, size(order)
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (values(order(j)) <= values(key)) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = key
    end do
  end subroutine sort_indices

  real(dp) function inverse_normal_cdf(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a1=-3.969683028665376e1_dp, a2=2.209460984245205e2_dp
    real(dp), parameter :: a3=-2.759285104469687e2_dp, a4=1.383577518672690e2_dp
    real(dp), parameter :: a5=-3.066479806614716e1_dp, a6=2.506628277459239_dp
    real(dp), parameter :: b1=-5.447609879822406e1_dp, b2=1.615858368580409e2_dp
    real(dp), parameter :: b3=-1.556989798598866e2_dp, b4=6.680131188771972e1_dp
    real(dp), parameter :: b5=-1.328068155288572e1_dp
    real(dp), parameter :: c1=-7.784894002430293e-3_dp, c2=-3.223964580411365e-1_dp
    real(dp), parameter :: c3=-2.400758277161838_dp, c4=-2.549732539343734_dp
    real(dp), parameter :: c5=4.374664141464968_dp, c6=2.938163982698783_dp
    real(dp), parameter :: d1=7.784695709041462e-3_dp, d2=3.224671290700398e-1_dp
    real(dp), parameter :: d3=2.445134137142996_dp, d4=3.754408661907416_dp
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
  end function inverse_normal_cdf

end module strand_stats
