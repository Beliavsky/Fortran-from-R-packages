! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_dowd
  use dowd
  implicit none

  integer :: failures
  real(dp) :: value, expected, var_value, es_value, p_value, lr
  real(dp) :: cov(2,2), mu_vec(2), positions(2), inc(2)
  real(dp) :: eig(2), evec(2,2), explained(2)
  real(dp) :: pl(8), vars(8), ess(8)
  real(dp) :: ci_var(2), ci_es(2)
  real(dp) :: data_evt(12)

  failures = 0
  call set_random_seed_scalar(20260727)

  call assert_close("normal quantile", normal_quantile(0.975_dp), 1.959963984540054_dp, 2.0e-9_dp)
  call assert_close("student t quantile", student_t_quantile(0.975_dp,6.0_dp), 2.446911851144969_dp, 2.0e-8_dp)

  value = normal_var(0.001_dp,0.02_dp,0.99_dp,10.0_dp)
  call assert_close("normal VaR",value,0.13713115823719108_dp,2.0e-10_dp)
  value = normal_es(0.001_dp,0.02_dp,0.99_dp,10.0_dp)
  call assert_close("normal ES",value,0.15856294777125254_dp,2.0e-10_dp)
  value = student_t_var(0.001_dp,0.02_dp,6.0_dp,0.99_dp,10.0_dp)
  call assert_close("student t VaR",value,0.1522866985146443_dp,2.0e-8_dp)
  value = student_t_es(0.001_dp,0.02_dp,6.0_dp,0.99_dp,10.0_dp)
  call assert_close("student t ES",value,0.19823883394483235_dp,3.0e-8_dp)
  value = lognormal_var(100.0_dp,0.001_dp,0.02_dp,0.99_dp,10.0_dp)
  call assert_close("lognormal VaR",value,12.814413244593448_dp,2.0e-8_dp)
  value = lognormal_es(100.0_dp,0.001_dp,0.02_dp,0.99_dp,10.0_dp)
  call assert_close("lognormal ES",value,14.64671516452266_dp,3.0e-8_dp)

  pl = [-1.0_dp,2.0_dp,-3.0_dp,4.0_dp,-5.0_dp,0.5_dp,-0.25_dp,1.5_dp]
  call historical_var_es(pl,0.75_dp,var_value,es_value)
  call assert_true("historical ES >= VaR",es_value >= var_value)
  call assert_close("historical VaR",var_value,1.5_dp,1.0e-12_dp)

  value = cornish_fisher_var(0.0_dp,1.0_dp,0.0_dp,3.0_dp,0.99_dp)
  call assert_close("Cornish-Fisher normal limit",value,-normal_quantile(0.01_dp),2.0e-10_dp)
  call assert_true("Cornish-Fisher ES ordering",cornish_fisher_es(0.0_dp,1.0_dp,0.2_dp,3.5_dp,0.99_dp) > &
                   cornish_fisher_var(0.0_dp,1.0_dp,0.2_dp,3.5_dp,0.99_dp))

  call assert_close("Gumbel VaR",gumbel_var(0.0_dp,1.0_dp,100.0_dp,0.99_dp), &
                    -log(-100.0_dp*log(0.99_dp)),1.0e-12_dp)
  call assert_true("Frechet ES ordering",frechet_es(0.0_dp,1.0_dp,0.2_dp,100.0_dp,0.99_dp) > &
                   frechet_var(0.0_dp,1.0_dp,0.2_dp,100.0_dp,0.99_dp))
  data_evt = [1.0_dp,1.2_dp,1.4_dp,1.8_dp,2.0_dp,2.4_dp,3.0_dp,4.0_dp,5.0_dp,7.0_dp,10.0_dp,15.0_dp]
  call assert_true("Hill positive",hill_estimator(data_evt,4) > 0.0_dp)
  call assert_true("GPD ES ordering",generalized_pareto_es(data_evt,1.0_dp,0.2_dp,0.25_dp,0.99_dp) > &
                   generalized_pareto_var(data_evt,1.0_dp,0.2_dp,0.25_dp,0.99_dp))

  value = kernel_var(-pl,0.75_dp,kernel_gaussian)
  call assert_true("kernel VaR finite",abs(value) < huge(1.0_dp))
  call assert_true("kernel ES ordering",kernel_es(-pl,0.75_dp,kernel_gaussian) >= value)
  call bootstrap_var_es(pl,100,0.75_dp,var_value,es_value)
  call assert_true("bootstrap finite",abs(var_value) < huge(1.0_dp) .and. es_value >= var_value)
  call bootstrap_confidence_interval(pl,100,0.75_dp,0.05_dp,0.95_dp,ci_var,ci_es)
  call assert_true("bootstrap intervals ordered",ci_var(1) <= ci_var(2) .and. ci_es(1) <= ci_es(2))
  value = boxcox_var(pl,0.75_dp)
  call assert_true("Box-Cox finite",abs(value) < huge(1.0_dp))

  cov = reshape([0.04_dp,0.01_dp,0.01_dp,0.09_dp],[2,2])
  mu_vec = [0.001_dp,0.002_dp]
  positions = [1.0_dp,2.0_dp]
  expected = -0.005_dp-normal_quantile(0.05_dp)*sqrt(0.44_dp)
  value = variance_covariance_var(cov,mu_vec,positions,0.95_dp,1.0_dp)
  call assert_close("portfolio VaR",value,expected,3.0e-10_dp)
  call assert_true("portfolio ES ordering",variance_covariance_es(cov,mu_vec,positions,0.95_dp,1.0_dp) > value)
  call normal_var_hotspots(cov,mu_vec,positions,0.95_dp,1.0_dp,inc)
  call assert_true("hotspots finite",all(abs(inc) < huge(1.0_dp)))
  call pca_prelim(cov,eig,evec,explained)
  call assert_close("PCA explained sum",sum(explained),1.0_dp,1.0e-12_dp)
  call assert_close("PCA full VaR",pca_var(cov,mu_vec,positions,2,0.95_dp,1.0_dp),value,2.0e-10_dp)

  vars = 2.0_dp
  ess = 2.5_dp
  p_value = christoffersen_unconditional_coverage(pl,vars,0.75_dp,lr)
  call assert_true("UC p-value",p_value >= 0.0_dp .and. p_value <= 1.0_dp .and. lr >= 0.0_dp)
  p_value = christoffersen_independence(pl,vars,lr)
  call assert_true("independence p-value",p_value >= 0.0_dp .and. p_value <= 1.0_dp)
  call assert_true("binomial p-value",binomial_backtest(1,100,0.99_dp) > 0.0_dp)
  p_value = jarque_bera_backtest(0.0_dp,3.0_dp,100,lr)
  call assert_close("JB normal moments",p_value,1.0_dp,1.0e-12_dp)
  call assert_true("GOF stats positive",ks_statistic_normal([-1.0_dp,0.0_dp,1.0_dp]) > 0.0_dp .and. &
                   kuiper_statistic_normal([-1.0_dp,0.0_dp,1.0_dp]) > 0.0_dp .and. &
                   anderson_darling_statistic_normal([-1.0_dp,0.0_dp,1.0_dp]) > 0.0_dp)

  call assert_close("Black-Scholes call",black_scholes_call_price(100.0_dp,105.0_dp,0.03_dp,0.2_dp,0.5_dp), &
                    4.178299715513496_dp,3.0e-10_dp)
  call assert_close("Black-Scholes put",black_scholes_put_price(100.0_dp,105.0_dp,0.03_dp,0.2_dp,0.5_dp), &
                    7.615053373835067_dp,3.0e-10_dp)
  value = american_put_price_binomial(100.0_dp,105.0_dp,0.03_dp,0.2_dp,180.0_dp,500)
  call assert_true("American put >= European",value >= black_scholes_put_price(100.0_dp,105.0_dp,0.03_dp,0.2_dp,0.5_dp))
  call assert_true("American put ES >= VaR",american_put_es_binomial(10000.0_dp,100.0_dp,105.0_dp,0.03_dp,0.2_dp, &
                   180.0_dp,100,0.95_dp,10.0_dp) >= american_put_var_binomial(10000.0_dp,100.0_dp,105.0_dp,0.03_dp, &
                   0.2_dp,180.0_dp,100,0.95_dp,10.0_dp))

  call assert_close("product copula",product_copula(0.2_dp,0.3_dp),0.06_dp,1.0e-15_dp)
  call assert_close("Gumbel beta 1",gumbel_copula(0.2_dp,0.3_dp,1.0_dp),0.06_dp,1.0e-14_dp)
  call assert_close("product sum CDF",cdf_sum_product_copula(0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp),0.5_dp,1.0e-14_dp)
  call assert_close("Gaussian copula VaR",gaussian_copula_var(0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp,0.99_dp), &
                    -sqrt(3.0_dp)*normal_quantile(0.01_dp),2.0e-10_dp)
  call assert_close("Gumbel beta1 sum CDF",cdf_sum_gumbel_copula(0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,800), &
                    0.5_dp,2.0e-3_dp)

  if (failures > 0) then
    write(*,'(a,i0)') "FAILURES: ",failures
    error stop 1
  end if
  write(*,'(a)') "All Dowd tests passed."

contains

  subroutine assert_close(name, actual, target, tolerance)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: actual, target, tolerance
    if (abs(actual-target) > tolerance*max(1.0_dp,abs(target))) then
      failures = failures+1
      write(*,'(a,2(1x,es24.16))') "FAIL "//trim(name)//":",actual,target
    end if
  end subroutine assert_close

  subroutine assert_true(name, condition)
    character(len=*), intent(in) :: name
    logical, intent(in) :: condition
    if (.not.condition) then
      failures = failures+1
      write(*,'(a)') "FAIL "//trim(name)
    end if
  end subroutine assert_true

end program test_dowd
