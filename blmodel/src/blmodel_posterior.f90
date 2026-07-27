! SPDX-License-Identifier: GPL-3.0-only
! Derived from BLModel 1.0.2, Copyright (C) 2017 Andrzej Palczewski and Jan Palczewski.
module blmodel_posterior
  use blmodel_kinds, only : dp
  use blmodel_distributions, only : view_density_interface
  use blmodel_linalg, only : symmetric_matrix
  use blmodel_utils, only : diag_of, make_diag
  use blmodel_equilibrium, only : discrete_variance, equilibrium_mean, equilibrium_mean_elliptic
  use blmodel_types, only : moment_result, equilibrium_result, posterior_result
  implicit none
  private

  public :: post_distribution, bl_post_distribution

contains

  function post_distribution(returns, probabilities, equilibrium_returns, q, pick, covariance, tau, &
      view_density, view_covariance_type, view_params) result(res)
    real(dp), intent(in) :: returns(:,:), probabilities(:), equilibrium_returns(:), q(:)
    real(dp), intent(in) :: pick(:,:), covariance(:,:), tau
    procedure(view_density_interface) :: view_density
    character(len=*), intent(in) :: view_covariance_type
    real(dp), intent(in), optional :: view_params(:)
    type(posterior_result) :: res
    real(dp), allocatable :: weights(:), prior_mean(:), shifted(:,:), auxiliary(:,:), points(:,:), density(:)
    real(dp) :: total_probability, posterior_total
    integer :: i, info, n, k, m
    character(len=:), allocatable :: covariance_type

    n = size(returns, 1)
    k = size(returns, 2)
    m = size(pick, 1)
    if (n == 0 .or. k == 0 .or. size(probabilities) /= n .or. size(equilibrium_returns) /= k) then
      call fail_posterior(res, 'incompatible returns, probabilities, or equilibrium mean dimensions')
      return
    end if
    if (size(pick, 2) /= k .or. size(q) /= m) then
      call fail_posterior(res, 'pick matrix and view vector dimensions are incompatible')
      return
    end if
    if (size(covariance, 1) /= k .or. size(covariance, 2) /= k .or. tau <= 0.0_dp) then
      call fail_posterior(res, 'covariance dimensions are invalid or tau is not positive')
      return
    end if
    if (any(probabilities < 0.0_dp)) then
      call fail_posterior(res, 'probabilities must be nonnegative')
      return
    end if

    total_probability = sum(probabilities)
    if (total_probability <= 0.0_dp) then
      call fail_posterior(res, 'probabilities must have positive sum')
      return
    end if

    allocate(weights(n), prior_mean(k), shifted(n, k), auxiliary(m, m), points(m, n))
    weights = probabilities / total_probability
    prior_mean = matmul(weights, returns)
    shifted = returns
    do i = 1, n
      shifted(i, :) = shifted(i, :) - prior_mean + equilibrium_returns
    end do

    auxiliary = symmetric_matrix(matmul(pick, matmul(covariance, transpose(pick))) / tau)
    covariance_type = lowercase(trim(view_covariance_type))
    allocate(res%view_covariance(m, m), source=0.0_dp)
    select case (covariance_type)
    case ('diag', 'diagonal')
      res%view_covariance = make_diag(diag_of(auxiliary))
    case ('full')
      res%view_covariance = auxiliary
    case default
      call fail_posterior(res, 'view covariance type must be diag or full')
      return
    end select

    points = matmul(pick, transpose(shifted))
    if (present(view_params)) then
      call view_density(points, q, res%view_covariance, view_params, density, info)
    else
      call view_density(points, q, res%view_covariance, density=density, info=info)
    end if
    if (info /= 0 .or. size(density) /= n .or. any(density < 0.0_dp)) then
      call fail_posterior(res, 'view-density callback failed or returned invalid values')
      return
    end if

    allocate(res%returns(n, k), res%probabilities(n), res%equilibrium_returns(k))
    res%returns = shifted
    res%equilibrium_returns = equilibrium_returns
    res%probabilities = weights * density
    posterior_total = sum(res%probabilities)
    if (posterior_total <= tiny(1.0_dp)) then
      call fail_posterior(res, 'posterior density weights sum to zero')
      return
    end if
    res%probabilities = res%probabilities / posterior_total
    res%ok = .true.
    res%message = ''
  end function post_distribution

  function bl_post_distribution(returns, probabilities, returns_freq, prior_type, market_portfolio, &
      sharpe_ratio, pick, q, tau, risk, alpha, view_density, view_covariance_type, covariance, view_params) result(res)
    real(dp), intent(in) :: returns(:,:), probabilities(:), returns_freq
    character(len=*), intent(in) :: prior_type, risk, view_covariance_type
    real(dp), intent(in) :: market_portfolio(:), sharpe_ratio, pick(:,:), q(:), tau
    real(dp), intent(in), optional :: alpha, covariance(:,:), view_params(:)
    procedure(view_density_interface) :: view_density
    type(posterior_result) :: res
    type(moment_result) :: moments
    type(equilibrium_result) :: equilibrium
    real(dp), allocatable :: sample_covariance(:,:), normalized_portfolio(:), q_period(:)
    real(dp) :: portfolio_sum, portfolio_variance, market_return, alpha_value, market_price_of_risk
    character(len=:), allocatable :: prior_name
    integer :: k

    k = size(returns, 2)
    if (returns_freq <= 0.0_dp .or. sharpe_ratio < 0.0_dp .or. size(market_portfolio) /= k) then
      call fail_posterior(res, 'invalid frequency, Sharpe ratio, or portfolio dimension')
      return
    end if

    portfolio_sum = sum(market_portfolio)
    if (abs(portfolio_sum) <= 100.0_dp * epsilon(1.0_dp)) then
      call fail_posterior(res, 'market portfolio weights must have nonzero sum')
      return
    end if
    allocate(normalized_portfolio(k), sample_covariance(k, k), q_period(size(q)))
    normalized_portfolio = market_portfolio / portfolio_sum

    if (present(covariance)) then
      if (size(covariance, 1) /= k .or. size(covariance, 2) /= k) then
        call fail_posterior(res, 'supplied covariance has incompatible dimensions')
        return
      end if
      sample_covariance = covariance
    else
      moments = discrete_variance(returns, probabilities, returns_freq)
      if (.not. moments%ok) then
        call fail_posterior(res, moments%message)
        return
      end if
      sample_covariance = moments%covariance
    end if

    portfolio_variance = dot_product(normalized_portfolio, matmul(sample_covariance, normalized_portfolio))
    if (portfolio_variance <= 0.0_dp) then
      call fail_posterior(res, 'market portfolio variance must be positive')
      return
    end if

    prior_name = lowercase(trim(prior_type))
    select case (prior_name)
    case ('general', 'none', 'nonelliptic', 'non-elliptic')
      market_return = sharpe_ratio * sqrt(portfolio_variance) / returns_freq
      alpha_value = 0.95_dp
      if (present(alpha)) alpha_value = alpha
      equilibrium = equilibrium_mean(returns, probabilities, normalized_portfolio, market_return, risk, alpha_value)
    case ('elliptic', 'elliptical')
      market_price_of_risk = sharpe_ratio / sqrt(portfolio_variance)
      equilibrium = equilibrium_mean_elliptic(sample_covariance / returns_freq, normalized_portfolio, &
        market_price_of_risk)
    case default
      call fail_posterior(res, 'prior_type must be elliptic or general')
      return
    end select
    if (.not. equilibrium%ok) then
      call fail_posterior(res, equilibrium%message)
      return
    end if

    q_period = q / returns_freq
    if (present(view_params)) then
      res = post_distribution(returns, probabilities, equilibrium%market_returns, q_period, pick, &
        sample_covariance / returns_freq, tau, view_density, view_covariance_type, view_params)
    else
      res = post_distribution(returns, probabilities, equilibrium%market_returns, q_period, pick, &
        sample_covariance / returns_freq, tau, view_density, view_covariance_type)
    end if
  end function bl_post_distribution

  pure function lowercase(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code

    out = text
    do i = 1, len(text)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
    end do
  end function lowercase

  subroutine fail_posterior(res, message)
    type(posterior_result), intent(out) :: res
    character(len=*), intent(in) :: message

    res%ok = .false.
    res%message = message
  end subroutine fail_posterior

end module blmodel_posterior
