! SPDX-License-Identifier: GPL-3.0-or-later
module frapo_risk
  use frapo_kinds, only : dp
  use frapo_types, only : frapo_ok, frapo_invalid_input
  use frapo_statistics, only : is_symmetric, average_ranks
  implicit none
  private

  integer, parameter, public :: tdc_empirical = 1
  integer, parameter, public :: tdc_evt = 2

  public :: marginal_risk_contribution
  public :: diversification_ratio, concentration_ratio
  public :: volatility_weighted_correlation, tail_dependence_coefficient
  public :: mrc, dr, cr, rhow, tdc

  interface mrc
    module procedure marginal_risk_contribution
  end interface
  interface dr
    module procedure diversification_ratio
  end interface
  interface cr
    module procedure concentration_ratio
  end interface
  interface rhow
    module procedure volatility_weighted_correlation
  end interface
  interface tdc
    module procedure tail_dependence_coefficient
  end interface

contains

  function marginal_risk_contribution(weights, covariance, percentage, status) result(contribution)
    real(dp), intent(in) :: weights(:), covariance(:, :)
    logical, intent(in), optional :: percentage
    integer, intent(out), optional :: status
    real(dp), allocatable :: contribution(:)
    real(dp), allocatable :: sigma_w(:)
    real(dp) :: portfolio_sigma
    logical :: pct

    allocate(contribution(size(weights)))
    pct = .true.
    if (present(percentage)) pct = percentage
    if (.not. is_symmetric(covariance) .or. size(covariance, 1) /= size(weights)) then
      contribution = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    sigma_w = matmul(covariance, weights)
    portfolio_sigma = sqrt(dot_product(weights, sigma_w))
    if (portfolio_sigma <= 0.0_dp) then
      contribution = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    contribution = weights * sigma_w / portfolio_sigma
    if (pct) contribution = 100.0_dp * contribution / sum(contribution)
    if (present(status)) status = frapo_ok
  end function marginal_risk_contribution

  function diversification_ratio(weights, covariance, status) result(ratio)
    real(dp), intent(in) :: weights(:), covariance(:, :)
    integer, intent(out), optional :: status
    real(dp) :: ratio, numerator, denominator
    real(dp), allocatable :: sd(:)
    integer :: i, n

    n = size(weights)
    if (.not. is_symmetric(covariance) .or. size(covariance, 1) /= n) then
      ratio = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    allocate(sd(n))
    do i = 1, n
      sd(i) = sqrt(covariance(i, i))
    end do
    numerator = dot_product(weights, sd)
    denominator = sqrt(dot_product(weights, matmul(covariance, weights)))
    if (denominator <= 0.0_dp) then
      ratio = 0.0_dp
      if (present(status)) status = frapo_invalid_input
    else
      ratio = numerator / denominator
      if (present(status)) status = frapo_ok
    end if
  end function diversification_ratio

  function concentration_ratio(weights, covariance, status) result(ratio)
    real(dp), intent(in) :: weights(:), covariance(:, :)
    integer, intent(out), optional :: status
    real(dp) :: ratio, denominator
    real(dp), allocatable :: products(:)
    integer :: i, n

    n = size(weights)
    if (.not. is_symmetric(covariance) .or. size(covariance, 1) /= n) then
      ratio = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    allocate(products(n))
    do i = 1, n
      products(i) = weights(i) * sqrt(covariance(i, i))
    end do
    denominator = sum(products)**2
    if (denominator <= 0.0_dp) then
      ratio = 0.0_dp
      if (present(status)) status = frapo_invalid_input
    else
      ratio = sum(products**2) / denominator
      if (present(status)) status = frapo_ok
    end if
  end function concentration_ratio

  function volatility_weighted_correlation(weights, covariance, status) result(ratio)
    real(dp), intent(in) :: weights(:), covariance(:, :)
    integer, intent(out), optional :: status
    real(dp) :: ratio, numerator, denominator, product, corr
    real(dp), allocatable :: sd(:)
    integer :: i, j, n

    n = size(weights)
    if (.not. is_symmetric(covariance) .or. size(covariance, 1) /= n) then
      ratio = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    allocate(sd(n))
    do i = 1, n
      sd(i) = sqrt(covariance(i, i))
    end do
    numerator = 0.0_dp
    denominator = 0.0_dp
    do j = 2, n
      do i = 1, j - 1
        product = weights(i) * sd(i) * weights(j) * sd(j)
        corr = covariance(i, j) / (sd(i) * sd(j))
        numerator = numerator + product * corr
        denominator = denominator + product
      end do
    end do
    if (abs(denominator) <= tiny(1.0_dp)) then
      ratio = 0.0_dp
      if (present(status)) status = frapo_invalid_input
    else
      ratio = numerator / denominator
      if (present(status)) status = frapo_ok
    end if
  end function volatility_weighted_correlation

  function tail_dependence_coefficient(x, method, lower_tail, k, status) result(matrix)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in), optional :: method, k
    logical, intent(in), optional :: lower_tail
    integer, intent(out), optional :: status
    real(dp), allocatable :: matrix(:, :), ranks(:, :), column_ranks(:)
    integer :: meth, kk, m, n, i, j, t, count_value
    logical :: lower

    m = size(x, 1)
    n = size(x, 2)
    meth = tdc_empirical
    if (present(method)) meth = method
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    kk = int(sqrt(real(m, dp)))
    if (present(k)) kk = k
    allocate(matrix(n, n), ranks(m, n))
    matrix = 0.0_dp
    if (m < 1 .or. n < 1 .or. kk < 1 .or. kk > m .or. &
        (meth /= tdc_empirical .and. meth /= tdc_evt)) then
      if (present(status)) status = frapo_invalid_input
      return
    end if
    do j = 1, n
      call average_ranks(x(:, j), column_ranks)
      ranks(:, j) = column_ranks
    end do
    do i = 1, n
      matrix(i, i) = 1.0_dp
    end do
    do j = 2, n
      do i = 1, j - 1
        count_value = 0
        do t = 1, m
          if (meth == tdc_empirical) then
            if (lower) then
              if (ranks(t, i) <= real(kk, dp) .and. ranks(t, j) <= real(kk, dp)) count_value = count_value + 1
            else
              if (ranks(t, i) > real(m - kk, dp) .and. ranks(t, j) > real(m - kk, dp)) count_value = count_value + 1
            end if
          else
            if (lower) then
              if (ranks(t, i) <= real(kk, dp) .or. ranks(t, j) <= real(kk, dp)) count_value = count_value + 1
            else
              if (ranks(t, i) > real(m - kk, dp) .or. ranks(t, j) > real(m - kk, dp)) count_value = count_value + 1
            end if
          end if
        end do
        if (meth == tdc_empirical) then
          matrix(i, j) = real(count_value, dp) / real(kk, dp)
        else
          matrix(i, j) = 2.0_dp - real(count_value, dp) / real(kk, dp)
        end if
        matrix(j, i) = matrix(i, j)
      end do
    end do
    if (present(status)) status = frapo_ok
  end function tail_dependence_coefficient
end module frapo_risk
