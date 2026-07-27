! SPDX-License-Identifier: GPL-3.0-only
! Derived from BLModel 1.0.2, Copyright (C) 2017 Andrzej Palczewski and Jan Palczewski.
module blmodel_equilibrium
  use blmodel_kinds, only : dp
  use blmodel_linalg, only : sort_indices_ascending, symmetric_matrix
  use blmodel_types, only : moment_result, equilibrium_result
  implicit none
  private

  public :: discrete_variance, equilibrium_mean, equilibrium_mean_elliptic

contains

  function discrete_variance(returns, probabilities, returns_coef) result(res)
    real(dp), intent(in) :: returns(:,:), probabilities(:), returns_coef
    type(moment_result) :: res
    real(dp), allocatable :: weights(:), centered(:,:)
    real(dp) :: total_probability
    integer :: i, n, k

    n = size(returns, 1)
    k = size(returns, 2)
    if (n == 0 .or. k == 0 .or. size(probabilities) /= n) then
      call fail_moment(res, 'returns and probabilities have incompatible dimensions')
      return
    end if
    if (returns_coef <= 0.0_dp .or. any(probabilities < 0.0_dp)) then
      call fail_moment(res, 'returns_coef must be positive and probabilities nonnegative')
      return
    end if

    total_probability = sum(probabilities)
    if (total_probability <= 0.0_dp) then
      call fail_moment(res, 'probabilities must have positive sum')
      return
    end if

    allocate(weights(n), res%mean(k), res%covariance(k, k), centered(n, k))
    weights = probabilities / total_probability
    res%mean = matmul(weights, returns) * returns_coef
    centered = returns
    do i = 1, n
      centered(i, :) = centered(i, :) - res%mean / returns_coef
    end do
    res%covariance = 0.0_dp
    do i = 1, n
      res%covariance = res%covariance + weights(i) * outer_product(centered(i, :), centered(i, :))
    end do
    res%covariance = symmetric_matrix(res%covariance * returns_coef)
    res%ok = .true.
    res%message = ''
  end function discrete_variance

  function equilibrium_mean(returns, probabilities, weights_in, market_return, risk, alpha) result(res)
    real(dp), intent(in) :: returns(:,:), probabilities(:), weights_in(:), market_return
    character(len=*), intent(in) :: risk
    real(dp), intent(in), optional :: alpha
    type(equilibrium_result) :: res
    real(dp), allocatable :: weights(:), losses(:,:), centered(:,:), portfolio_losses(:)
    real(dp), allocatable :: sorted_weights(:), tail_vector(:)
    integer, allocatable :: order(:)
    real(dp) :: total_probability, alpha_value, cumulative, scale
    integer :: i, j, n, k, boundary
    character(len=:), allocatable :: risk_name

    n = size(returns, 1)
    k = size(returns, 2)
    if (n == 0 .or. k == 0 .or. size(probabilities) /= n .or. size(weights_in) /= k) then
      call fail_equilibrium(res, 'incompatible returns, probabilities, or portfolio dimensions')
      return
    end if
    if (market_return <= 0.0_dp .or. any(probabilities < 0.0_dp)) then
      call fail_equilibrium(res, 'market_return must be positive and probabilities nonnegative')
      return
    end if

    total_probability = sum(probabilities)
    if (total_probability <= 0.0_dp) then
      call fail_equilibrium(res, 'probabilities must have positive sum')
      return
    end if

    allocate(weights(n), losses(n, k), centered(n, k), portfolio_losses(n), tail_vector(k))
    weights = probabilities / total_probability
    losses = -returns
    centered = losses
    do i = 1, n
      centered(i, :) = centered(i, :) - matmul(weights, losses)
    end do
    portfolio_losses = matmul(centered, weights_in)
    call sort_indices_ascending(portfolio_losses, order)
    allocate(sorted_weights(n))
    sorted_weights = weights(order)
    tail_vector = 0.0_dp
    risk_name = uppercase(trim(risk))

    select case (risk_name)
    case ('CVAR', 'DCVAR')
      alpha_value = 0.95_dp
      if (present(alpha)) alpha_value = alpha
      if (alpha_value <= 0.0_dp .or. alpha_value >= 1.0_dp) then
        call fail_equilibrium(res, 'alpha must lie strictly between zero and one')
        return
      end if
      cumulative = 0.0_dp
      boundary = n
      do i = 1, n
        cumulative = cumulative + sorted_weights(i)
        if (cumulative >= alpha_value) then
          boundary = i
          exit
        end if
      end do
      if (boundary < n) then
        do j = boundary + 1, n
          tail_vector = tail_vector + sorted_weights(j) * centered(order(j), :)
        end do
      end if
      tail_vector = tail_vector + (cumulative - alpha_value) * centered(order(boundary), :)
    case ('LSAD', 'MAD')
      do j = 1, n
        if (portfolio_losses(order(j)) <= 0.0_dp) then
          tail_vector = tail_vector + sorted_weights(j) * centered(order(j), :)
        end if
      end do
    case default
      call fail_equilibrium(res, 'risk must be CVAR, DCVAR, LSAD, or MAD')
      return
    end select

    scale = dot_product(tail_vector, weights_in) / market_return
    if (abs(scale) <= 100.0_dp * epsilon(1.0_dp)) then
      call fail_equilibrium(res, 'inverse optimization is singular for these inputs')
      return
    end if

    allocate(res%market_returns(k), res%portfolio(k))
    res%market_returns = tail_vector / scale
    res%portfolio = weights_in
    res%ok = .true.
    res%message = ''
  end function equilibrium_mean

  function equilibrium_mean_elliptic(covariance, market_portfolio, market_price_of_risk) result(res)
    real(dp), intent(in) :: covariance(:,:), market_portfolio(:), market_price_of_risk
    type(equilibrium_result) :: res
    integer :: k

    k = size(market_portfolio)
    if (k == 0 .or. size(covariance, 1) /= k .or. size(covariance, 2) /= k) then
      call fail_equilibrium(res, 'covariance and portfolio dimensions are incompatible')
      return
    end if
    if (market_price_of_risk < 0.0_dp) then
      call fail_equilibrium(res, 'market price of risk must be nonnegative')
      return
    end if

    allocate(res%market_returns(k), res%portfolio(k))
    res%portfolio = market_portfolio
    res%market_returns = market_price_of_risk * matmul(covariance, market_portfolio)
    res%ok = .true.
    res%message = ''
  end function equilibrium_mean_elliptic

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: i

    do i = 1, size(x)
      a(i, :) = x(i) * y
    end do
  end function outer_product

  pure function uppercase(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code

    out = text
    do i = 1, len(text)
      code = iachar(out(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) out(i:i) = achar(code - 32)
    end do
  end function uppercase

  subroutine fail_moment(res, message)
    type(moment_result), intent(out) :: res
    character(len=*), intent(in) :: message

    res%ok = .false.
    res%message = message
  end subroutine fail_moment

  subroutine fail_equilibrium(res, message)
    type(equilibrium_result), intent(out) :: res
    character(len=*), intent(in) :: message

    res%ok = .false.
    res%message = message
  end subroutine fail_equilibrium

end module blmodel_equilibrium
