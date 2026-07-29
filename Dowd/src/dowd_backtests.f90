! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_backtests
  use dowd_kinds, only: dp
  use dowd_math, only: binomial_cdf, chi_square_cdf, normal_cdf, random_normal, &
       sorted_copy, quantile_linear
  implicit none
  private

  public :: binomial_backtest
  public :: christoffersen_unconditional_coverage
  public :: christoffersen_independence
  public :: christoffersen_conditional_coverage
  public :: lopez_backtest, blanco_ihle_backtest, jarque_bera_backtest
  public :: ks_statistic_normal, kuiper_statistic_normal, anderson_darling_statistic_normal
  public :: simulate_normal_gof_interval

contains

  pure real(dp) function xlogy(x, y) result(value)
    real(dp), intent(in) :: x, y
    if (abs(x) <= tiny(1.0_dp)) then
      value = 0.0_dp
    else if (y <= 0.0_dp) then
      value = -huge(1.0_dp)
    else
      value = x*log(y)
    end if
  end function xlogy

  real(dp) function binomial_backtest(exceptions, observations, cl) result(p_value)
    integer, intent(in) :: exceptions, observations
    real(dp), intent(in) :: cl
    real(dp) :: p
    if (exceptions < 0 .or. exceptions > observations .or. observations <= 0 .or. &
        cl <= 0.0_dp .or. cl >= 1.0_dp) error stop "binomial_backtest: invalid parameters"
    p = 1.0_dp-cl
    if (real(exceptions,dp) >= real(observations,dp)*p) then
      p_value = 1.0_dp-binomial_cdf(exceptions-1,observations,p)
    else
      p_value = binomial_cdf(exceptions,observations,p)
    end if
  end function binomial_backtest

  real(dp) function christoffersen_unconditional_coverage(profit_loss, var_values, cl, lr_stat) result(p_value)
    real(dp), intent(in) :: profit_loss(:), var_values(:), cl
    real(dp), intent(out), optional :: lr_stat
    integer :: n, x
    real(dp) :: p, phat, lr, log_null, log_alt
    n = size(profit_loss)
    if (size(var_values) /= n .or. n == 0) error stop "christoffersen_unconditional_coverage: size mismatch"
    p = 1.0_dp-cl
    x = count(-profit_loss > var_values)
    phat = real(x,dp)/real(n,dp)
    log_null = xlogy(real(x,dp),p)+xlogy(real(n-x,dp),1.0_dp-p)
    log_alt = xlogy(real(x,dp),phat)+xlogy(real(n-x,dp),1.0_dp-phat)
    lr = max(0.0_dp,-2.0_dp*(log_null-log_alt))
    if (present(lr_stat)) lr_stat = lr
    p_value = 1.0_dp-chi_square_cdf(lr,1.0_dp)
  end function christoffersen_unconditional_coverage

  real(dp) function christoffersen_independence(profit_loss, var_values, lr_stat) result(p_value)
    real(dp), intent(in) :: profit_loss(:), var_values(:)
    real(dp), intent(out), optional :: lr_stat
    integer :: n, i, t00, t01, t10, t11
    logical, allocatable :: hit(:)
    real(dp) :: pi0, pi1, pi_all, log_null, log_alt, lr
    n = size(profit_loss)
    if (size(var_values) /= n .or. n < 2) error stop "christoffersen_independence: size mismatch"
    allocate(hit(n))
    hit = -profit_loss > var_values
    t00=0; t01=0; t10=0; t11=0
    do i = 2, n
      if (.not.hit(i-1) .and. .not.hit(i)) t00=t00+1
      if (.not.hit(i-1) .and. hit(i)) t01=t01+1
      if (hit(i-1) .and. .not.hit(i)) t10=t10+1
      if (hit(i-1) .and. hit(i)) t11=t11+1
    end do
    pi0 = real(t01,dp)/real(max(1,t00+t01),dp)
    pi1 = real(t11,dp)/real(max(1,t10+t11),dp)
    pi_all = real(t01+t11,dp)/real(max(1,t00+t01+t10+t11),dp)
    log_null = xlogy(real(t00+t10,dp),1.0_dp-pi_all)+xlogy(real(t01+t11,dp),pi_all)
    log_alt = xlogy(real(t00,dp),1.0_dp-pi0)+xlogy(real(t01,dp),pi0)+ &
              xlogy(real(t10,dp),1.0_dp-pi1)+xlogy(real(t11,dp),pi1)
    lr = max(0.0_dp,-2.0_dp*(log_null-log_alt))
    if (present(lr_stat)) lr_stat = lr
    p_value = 1.0_dp-chi_square_cdf(lr,1.0_dp)
  end function christoffersen_independence

  real(dp) function christoffersen_conditional_coverage(profit_loss, var_values, cl, lr_stat) result(p_value)
    real(dp), intent(in) :: profit_loss(:), var_values(:), cl
    real(dp), intent(out), optional :: lr_stat
    real(dp) :: lr_uc, lr_ind, ignored, lr
    ignored = christoffersen_unconditional_coverage(profit_loss,var_values,cl,lr_uc)
    ignored = christoffersen_independence(profit_loss,var_values,lr_ind)
    lr = lr_uc+lr_ind
    if (present(lr_stat)) lr_stat = lr
    p_value = 1.0_dp-chi_square_cdf(lr,2.0_dp)
  end function christoffersen_conditional_coverage

  pure real(dp) function lopez_backtest(profit_loss, var_values, cl) result(score)
    real(dp), intent(in) :: profit_loss(:), var_values(:), cl
    integer :: n, exceptions
    real(dp) :: expected
    n = size(profit_loss)
    if (size(var_values) /= n .or. n == 0) error stop "lopez_backtest: size mismatch"
    exceptions = count(-profit_loss > var_values)
    expected = real(n,dp)*(1.0_dp-cl)
    score = 2.0_dp*(real(exceptions,dp)-expected)**2/real(n,dp)
  end function lopez_backtest

  pure real(dp) function blanco_ihle_backtest(profit_loss, var_values, es_values, cl) result(score)
    real(dp), intent(in) :: profit_loss(:), var_values(:), es_values(:), cl
    real(dp) :: s
    integer :: i, n
    n = size(profit_loss)
    if (size(var_values) /= n .or. size(es_values) /= n .or. n == 0) &
      error stop "blanco_ihle_backtest: size mismatch"
    if (cl <= 0.0_dp .or. cl >= 1.0_dp) error stop "blanco_ihle_backtest: invalid confidence level"
    s = 0.0_dp
    do i = 1, n
      if (-profit_loss(i) > var_values(i) .and. abs(var_values(i)) > tiny(1.0_dp)) then
        s = s+((-profit_loss(i)-var_values(i))/var_values(i) - &
               (es_values(i)-var_values(i))/var_values(i))
      end if
    end do
    score = 2.0_dp*s*s/real(n,dp)
  end function blanco_ihle_backtest

  real(dp) function jarque_bera_backtest(sample_skewness, sample_kurtosis, n, statistic) result(p_value)
    real(dp), intent(in) :: sample_skewness, sample_kurtosis
    integer, intent(in) :: n
    real(dp), intent(out), optional :: statistic
    real(dp) :: jb
    if (n <= 0) error stop "jarque_bera_backtest: n must be positive"
    jb = real(n,dp)/6.0_dp*(sample_skewness**2+(sample_kurtosis-3.0_dp)**2/4.0_dp)
    if (present(statistic)) statistic = jb
    p_value = 1.0_dp-chi_square_cdf(jb,2.0_dp)
  end function jarque_bera_backtest

  real(dp) function ks_statistic_normal(sample) result(statistic)
    real(dp), intent(in) :: sample(:)
    real(dp), allocatable :: x(:)
    real(dp) :: dplus, dminus
    integer :: i, n
    n = size(sample)
    if (n == 0) error stop "ks_statistic_normal: empty sample"
    x = sorted_copy(sample)
    dplus = 0.0_dp
    dminus = 0.0_dp
    do i = 1, n
      dplus = max(dplus,real(i,dp)/real(n,dp)-normal_cdf(x(i)))
      dminus = max(dminus,normal_cdf(x(i))-real(i-1,dp)/real(n,dp))
    end do
    statistic = max(dplus,dminus)
  end function ks_statistic_normal

  real(dp) function kuiper_statistic_normal(sample) result(statistic)
    real(dp), intent(in) :: sample(:)
    real(dp), allocatable :: x(:)
    real(dp) :: dplus, dminus
    integer :: i, n
    n = size(sample)
    if (n == 0) error stop "kuiper_statistic_normal: empty sample"
    x = sorted_copy(sample)
    dplus = 0.0_dp
    dminus = 0.0_dp
    do i = 1, n
      dplus = max(dplus,real(i,dp)/real(n,dp)-normal_cdf(x(i)))
      dminus = max(dminus,normal_cdf(x(i))-real(i-1,dp)/real(n,dp))
    end do
    statistic = dplus+dminus
  end function kuiper_statistic_normal

  real(dp) function anderson_darling_statistic_normal(sample) result(statistic)
    real(dp), intent(in) :: sample(:)
    real(dp), allocatable :: x(:)
    real(dp) :: p1, p2, total
    integer :: i, n
    n = size(sample)
    if (n == 0) error stop "anderson_darling_statistic_normal: empty sample"
    x = sorted_copy(sample)
    total = 0.0_dp
    do i = 1, n
      p1 = min(max(normal_cdf(x(i)),tiny(1.0_dp)),1.0_dp-epsilon(1.0_dp))
      p2 = min(max(normal_cdf(x(n+1-i)),tiny(1.0_dp)),1.0_dp-epsilon(1.0_dp))
      total = total+real(2*i-1,dp)*(log(p1)+log(1.0_dp-p2))
    end do
    statistic = -real(n,dp)-total/real(n,dp)
  end function anderson_darling_statistic_normal

  subroutine simulate_normal_gof_interval(number_trials, sample_size, confidence_level, &
      ks_interval, kuiper_interval, ad_interval)
    integer, intent(in) :: number_trials, sample_size
    real(dp), intent(in) :: confidence_level
    real(dp), intent(out) :: ks_interval(2), kuiper_interval(2), ad_interval(2)
    real(dp), allocatable :: sample(:), ks(:), kuiper(:), ad(:)
    real(dp) :: lower_p, upper_p
    integer :: i, j
    if (number_trials <= 0 .or. sample_size <= 0 .or. confidence_level <= 0.0_dp .or. &
        confidence_level >= 1.0_dp) error stop "simulate_normal_gof_interval: invalid input"
    allocate(sample(sample_size),ks(number_trials),kuiper(number_trials),ad(number_trials))
    do i = 1, number_trials
      do j = 1, sample_size
        sample(j) = random_normal()
      end do
      ks(i) = ks_statistic_normal(sample)
      kuiper(i) = kuiper_statistic_normal(sample)
      ad(i) = anderson_darling_statistic_normal(sample)
    end do
    lower_p = 0.5_dp*(1.0_dp-confidence_level)
    upper_p = 1.0_dp-lower_p
    ks_interval = [quantile_linear(ks,lower_p),quantile_linear(ks,upper_p)]
    kuiper_interval = [quantile_linear(kuiper,lower_p),quantile_linear(kuiper,upper_p)]
    ad_interval = [quantile_linear(ad,lower_p),quantile_linear(ad,upper_p)]
  end subroutine simulate_normal_gof_interval

end module dowd_backtests
