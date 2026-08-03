! SPDX-License-Identifier: GPL-3.0-only
module pa_statistics
  use pa_kinds, only: dp
  use pa_linalg, only: sort_real, normal_quantile, normal_pdf
  implicit none
  private
  public :: sample_moments, weighted_moments, portfolio_returns
  public :: portfolio_mean, portfolio_variance, portfolio_stddev
  public :: historical_var, historical_es, gaussian_var, gaussian_es
  public :: semideviation, maximum_drawdown, turnover, hhi, diversification
  public :: risk_contributions, component_expected_shortfall
  public :: sample_coskewness, sample_cokurtosis, portfolio_skewness
  public :: portfolio_kurtosis, transaction_cost_value
  public :: conditional_second_moment, expected_quadratic_shortfall

contains

  subroutine sample_moments(r, mu, sigma)
    real(dp), intent(in) :: r(:,:)
    real(dp), intent(out) :: mu(:), sigma(:,:)
    real(dp), allocatable :: xc(:,:)
    integer :: nobs, nassets, i
    nobs = size(r,1)
    nassets = size(r,2)
    if (size(mu) /= nassets .or. size(sigma,1) /= nassets .or. size(sigma,2) /= nassets) return
    mu = sum(r, dim=1) / real(nobs, dp)
    allocate(xc(nobs,nassets))
    do i = 1, nobs
      xc(i,:) = r(i,:) - mu
    end do
    if (nobs > 1) then
      sigma = matmul(transpose(xc), xc) / real(nobs - 1, dp)
    else
      sigma = 0.0_dp
    end if
    sigma = 0.5_dp * (sigma + transpose(sigma))
  end subroutine sample_moments

  subroutine weighted_moments(r, probabilities, mu, sigma)
    real(dp), intent(in) :: r(:,:), probabilities(:)
    real(dp), intent(out) :: mu(:), sigma(:,:)
    real(dp), allocatable :: p(:), xc(:,:)
    real(dp) :: psum
    integer :: nobs, nassets, i, j
    nobs = size(r,1)
    nassets = size(r,2)
    if (size(probabilities) /= nobs) return
    allocate(p(nobs), xc(nobs,nassets))
    psum = sum(probabilities)
    if (psum <= 0.0_dp) return
    p = probabilities / psum
    mu = matmul(transpose(r), p)
    do i = 1, nobs
      xc(i,:) = r(i,:) - mu
    end do
    sigma = 0.0_dp
    do i = 1, nobs
      do j = 1, nassets
        sigma(j,:) = sigma(j,:) + p(i) * xc(i,j) * xc(i,:)
      end do
    end do
    sigma = 0.5_dp * (sigma + transpose(sigma))
  end subroutine weighted_moments

  subroutine portfolio_returns(r, weights, rp)
    real(dp), intent(in) :: r(:,:), weights(:)
    real(dp), intent(out) :: rp(:)
    if (size(r,2) /= size(weights) .or. size(rp) /= size(r,1)) return
    rp = matmul(r, weights)
  end subroutine portfolio_returns

  pure real(dp) function portfolio_mean(weights, mu) result(value)
    real(dp), intent(in) :: weights(:), mu(:)
    value = dot_product(weights, mu)
  end function portfolio_mean

  pure real(dp) function portfolio_variance(weights, sigma) result(value)
    real(dp), intent(in) :: weights(:), sigma(:,:)
    value = dot_product(weights, matmul(sigma, weights))
    value = max(value, 0.0_dp)
  end function portfolio_variance

  pure real(dp) function portfolio_stddev(weights, sigma) result(value)
    real(dp), intent(in) :: weights(:), sigma(:,:)
    value = sqrt(portfolio_variance(weights, sigma))
  end function portfolio_stddev

  real(dp) function historical_var(returns, alpha) result(value)
    real(dp), intent(in) :: returns(:), alpha
    real(dp), allocatable :: x(:)
    real(dp) :: pos, frac, q
    integer :: n, k
    n = size(returns)
    if (n == 0 .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    allocate(x(n))
    x = returns
    call sort_real(x)
    pos = 1.0_dp + alpha * real(n - 1, dp)
    k = floor(pos)
    frac = pos - real(k, dp)
    if (k >= n) then
      q = x(n)
    else
      q = (1.0_dp-frac)*x(k) + frac*x(k+1)
    end if
    value = -q
  end function historical_var

  real(dp) function historical_es(returns, alpha) result(value)
    real(dp), intent(in) :: returns(:), alpha
    real(dp), allocatable :: x(:)
    real(dp) :: count, tailmass, remaining, take
    integer :: n, i
    n = size(returns)
    if (n == 0 .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    allocate(x(n))
    x = returns
    call sort_real(x)
    tailmass = alpha * real(n, dp)
    count = 0.0_dp
    value = 0.0_dp
    remaining = tailmass
    do i = 1, n
      if (remaining <= 0.0_dp) exit
      take = min(1.0_dp, remaining)
      value = value - take*x(i)
      count = count + take
      remaining = remaining - take
    end do
    if (count > 0.0_dp) value = value / count
  end function historical_es

  pure real(dp) function gaussian_var(mean_return, sd_return, alpha) result(value)
    real(dp), intent(in) :: mean_return, sd_return, alpha
    value = -(mean_return + sd_return*normal_quantile(alpha))
  end function gaussian_var

  pure real(dp) function gaussian_es(mean_return, sd_return, alpha) result(value)
    real(dp), intent(in) :: mean_return, sd_return, alpha
    real(dp) :: z
    z = normal_quantile(alpha)
    value = -mean_return + sd_return*normal_pdf(z)/alpha
  end function gaussian_es

  real(dp) function semideviation(returns, target) result(value)
    real(dp), intent(in) :: returns(:), target
    real(dp) :: d
    integer :: i, n
    value = 0.0_dp
    n = size(returns)
    if (n == 0) return
    do i = 1, n
      d = min(returns(i)-target, 0.0_dp)
      value = value + d*d
    end do
    value = sqrt(value / real(n,dp))
  end function semideviation

  real(dp) function maximum_drawdown(returns) result(value)
    real(dp), intent(in) :: returns(:)
    real(dp) :: wealth, peak, dd
    integer :: i
    wealth = 1.0_dp
    peak = 1.0_dp
    value = 0.0_dp
    do i = 1, size(returns)
      wealth = wealth * (1.0_dp + returns(i))
      peak = max(peak, wealth)
      if (peak > 0.0_dp) then
        dd = 1.0_dp - wealth/peak
        value = max(value, dd)
      end if
    end do
  end function maximum_drawdown

  pure real(dp) function turnover(weights, initial_weights) result(value)
    real(dp), intent(in) :: weights(:), initial_weights(:)
    if (size(weights) /= size(initial_weights) .or. size(weights) == 0) then
      value = huge(1.0_dp)
    else
      value = sum(abs(weights-initial_weights)) / real(size(weights),dp)
    end if
  end function turnover

  pure real(dp) function hhi(weights) result(value)
    real(dp), intent(in) :: weights(:)
    value = sum(weights*weights)
  end function hhi

  pure real(dp) function diversification(weights) result(value)
    real(dp), intent(in) :: weights(:)
    value = 1.0_dp - hhi(weights)
  end function diversification

  subroutine risk_contributions(weights, sigma, contributions)
    real(dp), intent(in) :: weights(:), sigma(:,:)
    real(dp), intent(out) :: contributions(:)
    real(dp), allocatable :: marginal(:)
    real(dp) :: sd
    allocate(marginal(size(weights)))
    marginal = matmul(sigma, weights)
    sd = sqrt(max(dot_product(weights,marginal), 0.0_dp))
    if (sd > epsilon(1.0_dp)) then
      contributions = weights * marginal / sd
    else
      contributions = 0.0_dp
    end if
  end subroutine risk_contributions

  subroutine component_expected_shortfall(r, weights, alpha, contributions)
    real(dp), intent(in) :: r(:,:), weights(:), alpha
    real(dp), intent(out) :: contributions(:)
    real(dp), allocatable :: rp(:), sorted(:)
    real(dp) :: threshold, total, scale
    integer :: i, n, count
    n = size(r,1)
    allocate(rp(n), sorted(n))
    call portfolio_returns(r, weights, rp)
    sorted = rp
    call sort_real(sorted)
    threshold = sorted(max(1, ceiling(alpha*real(n,dp))))
    contributions = 0.0_dp
    count = 0
    do i = 1, n
      if (rp(i) <= threshold) then
        contributions = contributions - weights*r(i,:)
        count = count + 1
      end if
    end do
    if (count > 0) contributions = contributions / real(count,dp)
    total = sum(contributions)
    scale = historical_es(rp, alpha)
    if (abs(total) > epsilon(1.0_dp)) contributions = contributions * scale / total
  end subroutine component_expected_shortfall

  subroutine sample_coskewness(r, m3)
    real(dp), intent(in) :: r(:,:)
    real(dp), intent(out) :: m3(:,:)
    real(dp), allocatable :: mu(:), sigma(:,:), x(:,:)
    integer :: nobs, n, t, i, j, k, col
    nobs = size(r,1)
    n = size(r,2)
    allocate(mu(n),sigma(n,n),x(nobs,n))
    call sample_moments(r,mu,sigma)
    do t = 1, nobs
      x(t,:) = r(t,:) - mu
    end do
    m3 = 0.0_dp
    do i = 1, n
      do j = 1, n
        do k = 1, n
          col = (j-1)*n + k
          m3(i,col) = sum(x(:,i)*x(:,j)*x(:,k)) / real(nobs,dp)
        end do
      end do
    end do
  end subroutine sample_coskewness

  subroutine sample_cokurtosis(r, m4)
    real(dp), intent(in) :: r(:,:)
    real(dp), intent(out) :: m4(:,:)
    real(dp), allocatable :: mu(:), sigma(:,:), x(:,:)
    integer :: nobs, n, t, i, j, k, l, col
    nobs = size(r,1)
    n = size(r,2)
    allocate(mu(n),sigma(n,n),x(nobs,n))
    call sample_moments(r,mu,sigma)
    do t = 1, nobs
      x(t,:) = r(t,:) - mu
    end do
    m4 = 0.0_dp
    do i = 1, n
      do j = 1, n
        do k = 1, n
          do l = 1, n
            col = ((j-1)*n + (k-1))*n + l
            m4(i,col) = sum(x(:,i)*x(:,j)*x(:,k)*x(:,l)) / real(nobs,dp)
          end do
        end do
      end do
    end do
  end subroutine sample_cokurtosis

  pure real(dp) function portfolio_skewness(weights, m3) result(value)
    real(dp), intent(in) :: weights(:), m3(:,:)
    integer :: n, i, j, k, col
    n = size(weights)
    value = 0.0_dp
    do i = 1, n
      do j = 1, n
        do k = 1, n
          col = (j-1)*n+k
          value = value + weights(i)*weights(j)*weights(k)*m3(i,col)
        end do
      end do
    end do
  end function portfolio_skewness

  pure real(dp) function portfolio_kurtosis(weights, m4) result(value)
    real(dp), intent(in) :: weights(:), m4(:,:)
    integer :: n, i, j, k, l, col
    n = size(weights)
    value = 0.0_dp
    do i = 1, n
      do j = 1, n
        do k = 1, n
          do l = 1, n
            col = ((j-1)*n+(k-1))*n+l
            value = value + weights(i)*weights(j)*weights(k)*weights(l)*m4(i,col)
          end do
        end do
      end do
    end do
  end function portfolio_kurtosis

  pure real(dp) function transaction_cost_value(weights, initial_weights, rates) result(value)
    real(dp), intent(in) :: weights(:), initial_weights(:), rates(:)
    if (size(weights) /= size(initial_weights) .or. size(weights) /= size(rates)) then
      value = huge(1.0_dp)
    else
      value = sum(rates * abs(weights-initial_weights))
    end if
  end function transaction_cost_value

  real(dp) function conditional_second_moment(returns, alpha) result(value)
    real(dp), intent(in) :: returns(:), alpha
    real(dp) :: lo, hi, sd, phi, x1, x2, f1, f2
    integer :: iter
    if (size(returns)==0 .or. alpha<=0.0_dp .or. alpha>=1.0_dp) then
      value=huge(1.0_dp)
      return
    end if
    sd=sqrt(sum((returns-sum(returns)/real(size(returns),dp))**2)/max(1.0_dp,real(size(returns)-1,dp)))
    lo=minval(-returns)-10.0_dp*max(sd,1.0e-6_dp)
    hi=maxval(-returns)+max(sd,1.0e-6_dp)
    phi=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
    x1=hi-phi*(hi-lo)
    x2=lo+phi*(hi-lo)
    f1=csm_at(x1,returns,alpha)
    f2=csm_at(x2,returns,alpha)
    do iter=1,160
      if (f1>f2) then
        lo=x1; x1=x2; f1=f2; x2=lo+phi*(hi-lo); f2=csm_at(x2,returns,alpha)
      else
        hi=x2; x2=x1; f2=f1; x1=hi-phi*(hi-lo); f1=csm_at(x1,returns,alpha)
      end if
      if (abs(hi-lo)<1.0e-12_dp*max(1.0_dp,abs(lo)+abs(hi))) exit
    end do
    value=min(f1,f2)
  end function conditional_second_moment

  pure real(dp) function csm_at(zeta,returns,alpha) result(value)
    real(dp),intent(in)::zeta,returns(:),alpha
    real(dp)::z(size(returns))
    z=max(-returns-zeta,0.0_dp)
    value=zeta+sqrt(sum(z*z))/alpha
  end function csm_at

  real(dp) function expected_quadratic_shortfall(returns, alpha) result(value)
    real(dp), intent(in) :: returns(:), alpha
    real(dp), allocatable :: losses(:)
    integer :: n
    n=size(returns)
    if (n==0 .or. alpha<=0.0_dp .or. alpha>=1.0_dp) then
      value=huge(1.0_dp)
      return
    end if
    allocate(losses(n))
    losses=max(returns,0.0_dp)**2
    value=historical_es(-losses,alpha)
  end function expected_quadratic_shortfall

end module pa_statistics
