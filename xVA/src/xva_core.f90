module xva_core
  use trading, only : dp, trade_t, csa_t, normal_cdf
  use xva_types, only : regulatory_data_t
  use xva_math, only : inverse_normal_cdf, same_text
  implicit none
  private

  public :: calc_ngr
  public :: calc_pd
  public :: calc_va
  public :: generate_time_grid
  public :: is_eligible_currency
  public :: is_investment_grade
  public :: calc_effective_maturity
  public :: calc_default_capital
  public :: calc_kva
  public :: aggregate_geometric_pd
  public :: cem_addon_factor

contains

  real(dp) function calc_ngr(mtm_vector) result(ngr)
    real(dp), intent(in) :: mtm_vector(:)
    real(dp) :: gross_positive
    real(dp) :: net_value

    if (size(mtm_vector) == 0) then
      ngr = 1.0_dp
      return
    end if
    if (maxval(abs(mtm_vector)) <= tiny(1.0_dp)) then
      ngr = 1.0_dp
      return
    end if
    if (all(mtm_vector < 0.0_dp)) then
      ngr = 0.0_dp
      return
    end if

    net_value = max(sum(mtm_vector), 0.0_dp)
    gross_positive = sum(max(mtm_vector, 0.0_dp))
    if (gross_positive <= 0.0_dp) then
      ngr = 0.0_dp
    else
      ngr = min(max(net_value / gross_positive, 0.0_dp), 1.0_dp)
    end if
  end function calc_ngr

  subroutine calc_pd(spread_basis_points, lgd, time_points, pd)
    real(dp), intent(in) :: spread_basis_points(:)
    real(dp), intent(in) :: lgd
    real(dp), intent(in) :: time_points(:)
    real(dp), allocatable, intent(out) :: pd(:)
    real(dp), allocatable :: spread(:)
    integer :: n

    n = size(time_points)
    if (n < 2) error stop "calc_pd: at least two time points are required"
    if (size(spread_basis_points) /= n) then
      error stop "calc_pd: spread and time_points must have the same size"
    end if
    if (lgd <= 0.0_dp) error stop "calc_pd: lgd must be positive"
    if (any(time_points(2:) < time_points(:n - 1))) then
      error stop "calc_pd: time_points must be nondecreasing"
    end if

    allocate(spread(n), pd(n - 1))
    spread = spread_basis_points / 10000.0_dp
    pd = exp(-spread(:n - 1) * time_points(:n - 1) / lgd) - &
      exp(-spread(2:) * time_points(2:) / lgd)
    where (pd < 0.0_dp) pd = 0.0_dp
  end subroutine calc_pd

  real(dp) function calc_va(exposure, discount_factors, pd, lgd) result(value)
    real(dp), intent(in) :: exposure(:)
    real(dp), intent(in) :: discount_factors(:)
    real(dp), intent(in) :: pd(:)
    real(dp), intent(in), optional :: lgd
    real(dp) :: loss_given_default
    integer :: n

    n = size(discount_factors)
    if (n < 2) error stop "calc_va: at least two discount factors are required"
    if (size(pd) /= n - 1) then
      error stop "calc_va: pd must have one fewer element than discount_factors"
    end if
    if (size(exposure) < n - 1) then
      error stop "calc_va: exposure vector is too short"
    end if

    loss_given_default = 1.0_dp
    if (present(lgd)) loss_given_default = lgd
    if (loss_given_default < 0.0_dp) error stop "calc_va: lgd must be nonnegative"

    value = -loss_given_default * &
      sum(exposure(:n - 1) * discount_factors(:n - 1) * pd)
  end function calc_va

  subroutine generate_time_grid(csa, maturity, time_points)
    type(csa_t), intent(in) :: csa
    real(dp), intent(in) :: maturity
    real(dp), allocatable, intent(out) :: time_points(:)
    real(dp), allocatable :: regular(:)
    real(dp), allocatable :: lookback(:)
    real(dp) :: margin_period_years
    real(dp) :: remargin_years
    integer :: i
    integer :: j

    if (maturity < 0.0_dp) error stop "generate_time_grid: maturity must be nonnegative"
    remargin_years = csa%remargin_frequency / 360.0_dp
    margin_period_years = csa%mpor_days / 360.0_dp
    if (remargin_years <= 0.0_dp) then
      error stop "generate_time_grid: remargin_frequency must be positive"
    end if

    call arithmetic_sequence(0.0_dp, maturity, remargin_years, regular)
    call arithmetic_sequence(remargin_years - margin_period_years, &
      maturity - margin_period_years, remargin_years, lookback)

    allocate(time_points(size(regular) + size(lookback)))
    time_points(:size(regular)) = regular
    time_points(size(regular) + 1:) = lookback

    do i = 2, size(time_points)
      j = i
      do while (j > 1 .and. time_points(j) < time_points(j - 1))
        call swap_real(time_points(j), time_points(j - 1))
        j = j - 1
      end do
    end do
    call remove_duplicate_points(time_points)
  end subroutine generate_time_grid

  pure logical function is_eligible_currency(currency) result(eligible)
    character(len=*), intent(in) :: currency
    character(len=16), parameter :: eligible_currencies(8) = [ &
      character(len=16) :: "USD", "EUR", "GBP", "AUD", "CAD", "SEK", &
      "JPY", "REPORTING" ]
    integer :: i

    eligible = .false.
    do i = 1, size(eligible_currencies)
      if (same_text(currency, eligible_currencies(i))) then
        eligible = .true.
        return
      end if
    end do
  end function is_eligible_currency

  pure logical function is_investment_grade(rating) result(investment_grade)
    character(len=*), intent(in) :: rating
    character(len=8), parameter :: ratings(10) = [ &
      character(len=8) :: "AAA", "AA", "A", "BBB", "AA+", "A+", "BBB+", &
      "AA-", "A-", "BBB-" ]
    integer :: i

    investment_grade = .false.
    do i = 1, size(ratings)
      if (same_text(rating, ratings(i))) then
        investment_grade = .true.
        return
      end if
    end do
  end function is_investment_grade

  real(dp) function calc_effective_maturity(trades, time_points, framework, &
      simulated_exposure) result(effective_maturity)
    type(trade_t), intent(in) :: trades(:)
    real(dp), intent(in) :: time_points(:)
    character(len=*), intent(in) :: framework
    real(dp), intent(in), optional :: simulated_exposure(:)
    real(dp) :: denominator
    real(dp) :: numerator
    real(dp) :: total_notional
    integer :: i

    if (size(trades) == 0) error stop "calc_effective_maturity: trades are required"

    if (same_text(framework, "IMM")) then
      if (.not. present(simulated_exposure)) then
        error stop "calc_effective_maturity: simulated_exposure is required for IMM"
      end if
      if (size(simulated_exposure) /= size(time_points)) then
        error stop "calc_effective_maturity: exposure and time grid sizes differ"
      end if
      if (size(time_points) < 2) then
        error stop "calc_effective_maturity: at least two time points are required"
      end if

      numerator = 0.0_dp
      denominator = 0.0_dp
      do i = 1, size(time_points) - 1
        if (time_points(i) < 1.0_dp .and. time_points(i + 1) <= 1.0_dp) then
          denominator = denominator + simulated_exposure(i) * &
            (time_points(i + 1) - time_points(i))
        else if (time_points(i) >= 1.0_dp) then
          numerator = numerator + simulated_exposure(i + 1) * &
            (time_points(i + 1) - time_points(i))
        end if
      end do
      if (denominator <= tiny(1.0_dp)) then
        effective_maturity = 1.0_dp
      else
        effective_maturity = max(1.0_dp, min(numerator / denominator, 5.0_dp))
      end if
    else
      total_notional = sum(abs(trades%notional))
      if (total_notional <= tiny(1.0_dp)) then
        effective_maturity = max(1.0_dp, maxval(trades%ei))
      else
        effective_maturity = max(1.0_dp, &
          sum(abs(trades%notional) * trades%ei) / total_notional)
      end if
    end if
  end function calc_effective_maturity

  real(dp) function calc_default_capital(ead, reg_data, effective_maturity) &
      result(default_capital_charge)
    real(dp), intent(in) :: ead
    type(regulatory_data_t), intent(in) :: reg_data
    real(dp), intent(in) :: effective_maturity
    real(dp) :: b
    real(dp) :: capital_factor
    real(dp) :: correlation_basel_ii
    real(dp) :: correlation_stressed
    real(dp) :: pd_value

    pd_value = reg_data%pd_counterparty
    if (pd_value <= 0.0_dp .or. pd_value >= 1.0_dp) then
      error stop "calc_default_capital: pd_counterparty must be in (0,1)"
    end if
    if (reg_data%lgd < 0.0_dp) error stop "calc_default_capital: lgd must be nonnegative"

    correlation_basel_ii = 0.12_dp * (1.0_dp - exp(-50.0_dp * pd_value)) / &
      (1.0_dp - exp(-50.0_dp)) + 0.24_dp * &
      (1.0_dp - (1.0_dp - exp(-50.0_dp * pd_value)) / &
      (1.0_dp - exp(-50.0_dp)))
    correlation_stressed = min(1.25_dp * correlation_basel_ii, 0.999999_dp)
    b = (0.11852_dp - 0.05478_dp * log(pd_value))**2
    capital_factor = normal_cdf(inverse_normal_cdf(pd_value) / &
      sqrt(1.0_dp - correlation_stressed) + &
      sqrt(correlation_stressed / (1.0_dp - correlation_stressed)) * &
      inverse_normal_cdf(0.999_dp)) - pd_value
    capital_factor = capital_factor * &
      (1.0_dp + (effective_maturity - 2.5_dp) * b) / (1.0_dp - 1.5_dp * b)

    default_capital_charge = 0.08_dp * reg_data%lgd * capital_factor * ead
  end function calc_default_capital

  real(dp) function calc_kva(ead, reg_data, effective_maturity, &
      cva_capital_charge) result(kva)
    real(dp), intent(in) :: ead
    type(regulatory_data_t), intent(in) :: reg_data
    real(dp), intent(in) :: effective_maturity
    real(dp), intent(in) :: cva_capital_charge
    real(dp) :: default_capital_charge

    if (reg_data%ignore_default_charge) then
      default_capital_charge = 0.0_dp
    else
      default_capital_charge = calc_default_capital(ead, reg_data, effective_maturity)
    end if
    kva = -(default_capital_charge + cva_capital_charge) * 0.5_dp * &
      sqrt(max(effective_maturity, 0.0_dp)) * reg_data%return_on_capital
  end function calc_kva

  pure real(dp) function aggregate_geometric_pd(pd_value, years, survival_pd) &
      result(value)
    real(dp), intent(in) :: pd_value
    integer, intent(in) :: years
    real(dp), intent(in), optional :: survival_pd
    real(dp) :: survival_probability
    integer :: i

    survival_probability = pd_value
    if (present(survival_pd)) survival_probability = survival_pd
    value = 0.0_dp
    do i = 1, max(years, 0)
      value = value + pd_value * (1.0_dp - survival_probability)**(i - 1)
    end do
  end function aggregate_geometric_pd

  pure real(dp) function cem_addon_factor(maturity) result(value)
    real(dp), intent(in) :: maturity

    if (maturity <= 1.0_dp) then
      value = 0.0_dp
    else if (maturity <= 5.0_dp) then
      value = 0.005_dp
    else
      value = 0.015_dp
    end if
  end function cem_addon_factor

  subroutine arithmetic_sequence(first, last, increment, values)
    real(dp), intent(in) :: first
    real(dp), intent(in) :: last
    real(dp), intent(in) :: increment
    real(dp), allocatable, intent(out) :: values(:)
    real(dp) :: tolerance
    integer :: count
    integer :: i

    if (increment <= 0.0_dp) error stop "arithmetic_sequence: increment must be positive"
    tolerance = 100.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(first), abs(last), abs(increment))
    if (first > last + tolerance) then
      allocate(values(0))
      return
    end if
    count = floor((last - first + tolerance) / increment) + 1
    allocate(values(count))
    do i = 1, count
      values(i) = first + real(i - 1, dp) * increment
    end do
  end subroutine arithmetic_sequence


  subroutine remove_duplicate_points(values)
    real(dp), allocatable, intent(inout) :: values(:)
    real(dp), allocatable :: work(:)
    real(dp) :: tolerance
    integer :: i
    integer :: n_unique

    if (size(values) <= 1) return
    allocate(work(size(values)))
    work(1) = values(1)
    n_unique = 1
    do i = 2, size(values)
      tolerance = 100.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, abs(values(i)), abs(work(n_unique)))
      if (abs(values(i) - work(n_unique)) > tolerance) then
        n_unique = n_unique + 1
        work(n_unique) = values(i)
      end if
    end do
    deallocate(values)
    allocate(values(n_unique))
    values = work(:n_unique)
  end subroutine remove_duplicate_points

  pure subroutine swap_real(left, right)
    real(dp), intent(inout) :: left
    real(dp), intent(inout) :: right
    real(dp) :: temp

    temp = left
    left = right
    right = temp
  end subroutine swap_real

end module xva_core
