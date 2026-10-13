! SPDX-License-Identifier: GPL-2.0-or-later
program test_properties
   use esg
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   integer, parameter :: n = 40000
   type(randompack_rng_type) :: rng1, rng2, other
   type(esg_parameters) :: p, q
   type(esg_scenario) :: s, deterministic, missing
   real(dp), allocatable :: e(:,:,:), copy(:,:,:), paths(:,:), forward(:), value(:), other_value(:), third(:)
   real(dp), allocatable :: mean_path(:), residual(:), rate_change(:)
   real(dp) :: zero(61), short_e(2,2,5), tiny_e(2,0), small(2,2), check_draw(4)
   real(dp) :: mean, variance, target, covariance, deviation, nan, inf, discount, price, time
   integer :: status, i, j
   character(18), parameter :: kinds(9) = [character(18) :: 'Zero-Coupon', 'Bond', 'CBond', 'ConvBond', &
      'EuroCall_Stock', 'EuroPut_Stock', 'EuroCall_ZC', 'EuroPut_ZC', 'CDS']
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   zero = 0.03_dp
   p = esg_parameters()
   p%rho = -0.45_dp
   allocate(e(n,2,5), copy(n,2,5))
   rng1 = randompack_rng('pcg64', seed=314159)
   rng2 = randompack_rng('pcg64', seed=314159)
   other = randompack_rng('pcg64', seed=777)
   call esg_innovations(rng1, e, status)
   call require(status == 0, 'local RNG succeeds')
   call other%normal(check_draw, info=status)
   call esg_innovations(rng2, copy, status)
   call require(all(e == copy), 'local streams reproducible despite another stream being used')
   call esg_innovations(rng1, copy, status)
   call require(any(e /= copy), 'local stream advances')
   mean = sum(e(:,:,1)) / real(2 * n, dp)
   variance = sum((e(:,:,1) - mean) ** 2) / real(2 * n, dp)
   call require(abs(mean) < 0.015_dp .and. abs(variance - 1.0_dp) < 0.025_dp, 'normal moments')
   call esg_generate(p, zero, e, s, status, 'stock')
   call require(status == 0, 'stock generation succeeds')
   mean = sum(s%rate(:,2)) / real(n, dp)
   variance = sum((s%rate(:,2) - mean) ** 2) / real(n, dp)
   target = p%rate_vol ** 2 * (1.0_dp - exp(-2.0_dp * p%k)) / (2.0_dp * p%k)
   call require(abs(mean - 0.03_dp) < 0.0002_dp, 'one-step rate mean')
   call require(abs(variance / target - 1.0_dp) < 0.025_dp, 'one-step rate variance')
   residual = (log(s%stock(:,2) / p%stock0) - s%rate(:,2) + p%stock_vol ** 2 / 2.0_dp) / p%stock_vol
   rate_change = e(:,1,1)
   covariance = sum((residual - sum(residual) / real(n, dp)) &
      * (rate_change - sum(rate_change) / real(n, dp))) / real(n, dp)
   call require(abs(covariance - p%rho) < 0.02_dp, 'requested stock-rate innovation covariance')
   residual = log(s%stock(:,2) / p%stock0)
   mean = sum(residual) / real(n, dp)
   variance = sum((residual - mean) ** 2) / real(n, dp)
   target = target + p%stock_vol ** 2 + 2.0_dp * sqrt(target) * p%stock_vol * p%rho
   call require(abs(mean - (0.03_dp - p%stock_vol ** 2 / 2.0_dp)) < 0.004_dp, 'log-stock mean')
   call require(abs(variance / target - 1.0_dp) < 0.025_dp, 'log-stock variance including contemporaneous rate')
   call esg_martingale(s, mean_path, deviation, status)
   call require(status == 0 .and. deviation < 0.012_dp, 'discounted stock martingale moments')
   call esg_liquidity(p, e(:,:,4), paths, status)
   mean = sum(paths(:,3)) / real(n, dp)
   variance = sum((paths(:,3) - mean) ** 2) / real(n, dp)
   call require(abs(mean - p%liquidity0) < 0.00003_dp, 'liquidity random-walk mean')
   call require(abs(variance / (2.0_dp * p%eta ** 2) - 1.0_dp) < 0.03_dp, 'liquidity random-walk variance')

   call esg_default(p, e(:,:,5), paths, status)
   call require(status == 0, 'default one-step moments sample remains positive')
   mean = sum(paths(:,2)) / real(n, dp)
   variance = sum((paths(:,2) - mean) ** 2) / real(n, dp)
   target = p%default0 + p%alpha * (p%beta - p%default0) - p%default_vol / 4.0_dp + p%default_vol ** 2 / 4.0_dp
   call require(abs(mean - target) < 0.00003_dp, 'default recurrence analytic mean retains unsquared drift')
   target = p%default_vol ** 2 * p%default0 + p%default_vol ** 4 / 8.0_dp
   call require(abs(variance / target - 1.0_dp) < 0.03_dp, 'default recurrence analytic variance')

   short_e = 0.0_dp
   q = p
   q%rate_vol = 0.0_dp
   q%stock_vol = 0.0_dp
   q%estate_vol = 0.0_dp
   q%eta = 0.0_dp
   q%default_vol = 0.0_dp
   call esg_generate(q, zero, short_e, deterministic, status)
   call require(status == 0, 'zero vol deterministic generation')
   do i = 1, 3
      call require(all(abs(deterministic%rate(:,i) - 0.03_dp) < 1.0e-15_dp), 'constant rates')
      call require(all(abs(deterministic%stock(:,i) - p%stock0 * exp(0.03_dp * real(i - 1, dp))) &
         < 1.0e-12_dp), 'deterministic stock compound growth')
   end do
   call esg_zc_bond(deterministic, 1, 2.0_dp, value, status)
   discount = exp(-0.03_dp)
   call require(status == 0 .and. maxval(abs(value - discount)) < 1.0e-14_dp, 'analytic constant-rate discount')
   call esg_bond(deterministic, 1, 2.0_dp, 2, 0.025_dp, value, status)
   price = discount + 0.025_dp * (discount + exp(-0.015_dp) + 1.0_dp) &
      + 0.025_dp * (0.5_dp + 3.0_dp / 360.0_dp)
   call require(maxval(abs(value - price)) < 1.0e-14_dp, 'independent coupon schedule including accrued amount')
   call esg_stock_call(deterministic, 1, 2.0_dp, 100.0_dp, value, status)
   call require(status == 0 .and. maxval(abs(value - 100.0_dp * (1.0_dp - discount))) < 1.0e-12_dp, &
      'zero-variance stock call limit with source lag')
   call esg_zc_call(deterministic, 0, 1.0_dp, 1.0_dp, 0.9_dp, value, status)
   call require(status == 0 .and. abs(value(1) - 0.1_dp * exp(-0.03_dp)) < 1.0e-14_dp, 'zero H bond-option limit')

   call esg_generate(p, zero, short_e, s, status)
   call require(status == 0, 'all risks deterministic innovations')
   call esg_stock_call(s, 1, 2.0_dp, 95.0_dp, value, status)
   call esg_stock_put(s, 1, 2.0_dp, 95.0_dp, other_value, status)
   call esg_zc_bond(s, 1, 2.0_dp, third, status)
   call require(maxval(abs(value - other_value - (s%stock(:,1) - 95.0_dp * third))) < 1.0e-12_dp, &
      'stock call-put parity under source conventions')
   call esg_zc_call(s, 0, 1.0_dp, 2.0_dp, 0.95_dp, value, status)
   call esg_zc_put(s, 0, 1.0_dp, 2.0_dp, 0.95_dp, other_value, status)
   target = exp(-0.06_dp) - 0.95_dp * exp(-0.03_dp)
   call require(abs(value(1) - other_value(1) - target) < 1.0e-14_dp, 'bond option call-put parity')
   ! End-to-end generated interface and dispatcher must agree for every supported asset type.
   do j = 1, 9
      call esg_asset_price(s, trim(kinds(j)), 0, 1.0_dp, value, status, &
         coupons=2, coupon_rate=0.02_dp, recovery=0.4_dp, underlying=2.0_dp, strike=0.95_dp)
      call require(status == 0, 'asset dispatch status')
      call esg_simulate_asset(p, zero, short_e, trim(kinds(j)), 0, 1.0_dp, other_value, status, &
         coupons=2, coupon_rate=0.02_dp, recovery=0.4_dp, underlying=2.0_dp, strike=0.95_dp)
      call require(status == 0 .and. maxval(abs(value - other_value)) < 1.0e-13_dp, 'generate-and-price equivalence')
   end do
   call esg_asset_price(s, 'EuroCall_S', 0, 1.0_dp, value, status, strike=100.0_dp)
   call esg_stock_call(s, 0, 1.0_dp, 100.0_dp, other_value, status)
   call require(maxval(abs(value - other_value)) < 1.0e-13_dp, 'repaired upstream call alias')

   missing = s
   missing%p%alpha = 0.03_dp
   missing%p%beta = 0.08_dp
   missing%p%default_vol = 1.0e-10_dp
   call esg_corporate_bond(missing, 0, 2.0_dp, 2, 0.0_dp, 1.0_dp, value, status)
   target = exp(0.03_dp * ((1.0_dp - exp(-0.08_dp * 2.0_dp)) / 0.08_dp ** 2 - 2.0_dp / 0.08_dp))
   call require(status == 0 .and. abs(value(1) - target) < 1.0e-12_dp, 'tiny-vol LMN analytic principal limit')
   missing%p%default_vol = 1.0e-200_dp
   call esg_corporate_bond(missing, 0, 2.0_dp, 2, 0.0_dp, 1.0_dp, value, status)
   call require(status == 0 .and. abs(value(1) - target) < 1.0e-12_dp, 'underflow-vol LMN analytic limit')
   missing%p%beta = 0.0_dp
   call esg_corporate_bond(missing, 0, 2.0_dp, 2, 0.0_dp, 1.0_dp, value, status)
   call require(status == 0 .and. abs(value(1) - exp(-0.03_dp * 2.0_dp)) < 1.0e-12_dp, &
      'zero-beta underflow-vol LMN limit')

   missing%p%default_vol = 1.0e-10_dp
   call esg_corporate_bond(missing, 0, 2.0_dp, 2, 0.0_dp, 1.0_dp, value, status)
   call require(status == 0 .and. abs(value(1) - exp(-0.03_dp * 2.0_dp)) < 1.0e-12_dp, &
      'zero-beta tiny-vol LMN cancellation avoided')

   ! Explicit invalid, singular, nonfinite, shape, tail and overflow behavior.
   call esg_forward_rates(zero(:12), 2, forward, status)
   call require(status == 2, 'insufficient forward curve')
   call esg_zero_rates(zero, -1, forward, status)
   call require(status == 1, 'negative horizon')
   call esg_forward_rates([nan, zero], 0, forward, status)
   call require(status == 1, 'NaN curve rejected')
   q = p
   q%k = 0.0_dp
   call esg_generate(q, zero, short_e, missing, status)
   call require(status == 1, 'singular k rejected')
   q = p
   q%rho = 1.1_dp
   call esg_generate(q, zero, short_e, missing, status, 'stock')
   call require(status == 1, 'correlation outside interval rejected')
   q%rho = 1.0_dp
   call esg_generate(q, zero, short_e, missing, status, 'stock')
   call require(status == 0, 'perfect correlation allowed')
   q%rho = -1.0_dp
   call esg_generate(q, zero, short_e, missing, status, 'stock')
   call require(status == 0, 'perfect anticorrelation allowed')
   q = p
   q%stock0 = -1.0_dp
   call esg_generate(q, zero, short_e, missing, status, 'stock')
   call require(status == 1, 'negative stock rejected')
   q = p
   q%eta = inf
   call esg_liquidity(q, short_e(:,:,4), paths, status)
   call require(status == 1, 'infinite volatility rejected')
   call esg_generate(p, zero, short_e(:,:,:4), missing, status)
   call require(status == 2, 'risk dimension rejected')
   call esg_generate(p, zero, short_e, missing, status, 'unrecognized')
   call require(status == 1, 'unrecognized risk rejected')
   call esg_short_rate(p, [0.03_dp], tiny_e, paths, status)
   call require(status == 0 .and. size(paths, 2) == 1, 'zero horizon returns initial column')
   q = p
   q%default0 = 0.0_dp
   q%alpha = 0.0_dp
   q%default_vol = 0.5_dp
   call esg_default(q, short_e(:,:,5), paths, status)
   call require(status == 3 .and. all(paths(:,2) < 0.0_dp) .and. all(ieee_is_nan(paths(:,3))), &
      'negative default recurrence reported without clipping or sqrt of negative')
   small = 1.0e6_dp
   call esg_stock(p, s%rate, short_e(:,:,1), small, paths, status)
   call require(status == 4 .and. any(ieee_is_nan(paths)), 'stock overflow reported')
   call esg_asset_price(s, 'unknown', 0, 1.0_dp, value, status)
   call require(status == 1, 'unknown asset rejected')
   call esg_asset_price(s, 'Bond', 0, 1.0_dp, value, status)
   call require(status == 1, 'missing coupon arguments rejected')
   call esg_stock_put(s, 0, 1.0_dp, -1.0_dp, value, status)
   call require(status == 1, 'negative strike rejected')
   call esg_zc_bond(s, -1, 1.0_dp, value, status)
   call require(status == 1, 'negative valuation date rejected')
   call esg_zc_bond(s, 0, nan, value, status)
   call require(status == 1, 'NaN maturity rejected')
   call esg_zc_bond(s, 0, 5.0_dp, value, status)
   call require(status == 2, 'maturity past horizon rejected')
   call esg_bond(s, 0, 1.0_dp, 0, 0.02_dp, value, status)
   call require(status == 1, 'zero coupon frequency rejected')
   call esg_cds(s, 0, 1.0_dp, 1.1_dp, value, status)
   call require(status == 1, 'out of range omega rejected')
   call esg_cds(deterministic, 0, 1.0_dp, 0.4_dp, value, status)
   call require(status == 1, 'singular LMN zero volatility rejected')
   call esg_stock_call(s, 0, 1.0_dp, 1.0e100_dp, value, status)
   call require(status == 0 .and. all(abs(value) < 1.0e-100_dp), 'far out of money tail')
   time = esg_k(1.0_dp, 1.0e-12_dp)
   call require(abs(time - 1.0_dp) < 1.0e-12_dp, 'small mean-reversion stable kernel')
   print *, 'PASS: RNG, analytic moments, covariance, discounting, pricing, edge and invalid cases'
contains
   subroutine require(condition, label)
      logical, intent(in) :: condition !! Boolean invariant or expected status that must be true.
      character(*), intent(in) :: label !! Human-readable independent property or invalid-input case.
      if (.not. condition) then
         print *, 'FAIL ', label
         error stop 1
      end if
   end subroutine
end program
