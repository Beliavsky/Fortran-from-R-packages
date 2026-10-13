! SPDX-License-Identifier: LGPL-2.0-or-later
! Regression tests against the unchanged installed R fTrading 3042.79 package.
program test_r_fixtures
   use r_kinds, only: dp
   use ftrading_base
   use ftrading_indicators
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   type :: reference_record
      character(len=80) :: name = ''
      real(dp), allocatable :: values(:)
      logical :: used = .false.
   end type
   type(reference_record), allocatable :: records(:)
   real(dp), allocatable :: open(:), high(:), low(:), close(:), volume(:), levels(:), rates(:)
   real(dp), allocatable :: missing(:), positions(:), returns(:), positions_missing(:), returns_missing(:)
   real(dp), allocatable :: ohlc(:, :), signals(:), times(:)
   integer, allocatable :: indices(:)
   integer :: n, comparisons, k
   procedure(rolling_operation), pointer :: reducer

   reducer => square_sum
   comparisons = 0
   call load_references()
   open = fixture('input_open')
   high = fixture('input_high')
   low = fixture('input_low')
   close = fixture('input_close')
   volume = fixture('input_volume')
   levels = fixture('input_levels')
   rates = fixture('input_rates')
   missing = fixture('input_missing')
   positions = fixture('input_positions')
   returns = fixture('input_returns')
   positions_missing = fixture('input_positions_missing')
   returns_missing = fixture('input_returns_missing')
   n = size(close)
   allocate(ohlc(n, 4))
   ohlc(:, 1) = open
   ohlc(:, 2) = high
   ohlc(:, 3) = low
   ohlc(:, 4) = close
   call check('sharpeRatio_default', [sharpe_ratio(levels)])
   call check('sharpeRatio_scalar', &
      [sharpe_ratio(levels, risk_free = 0.003_dp, scale = 7.0_dp)])
   call check('sharpeRatio_rates', &
      [sharpe_ratio(levels, scale = 3.0_dp, rates = rates)])
   call check('sharpeRatio_short', [sharpe_ratio(levels(:2))])
   call check('sterlingRatio_default', [sterling_ratio(levels)])
   call check('sterlingRatio_rising', [sterling_ratio([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp])])
   call check_drawdown('maxDrawDown_levels', levels)
   call check_drawdown('maxDrawDown_ties', &
      [10.0_dp, 12.0_dp, 9.0_dp, 12.0_dp, 9.0_dp, 14.0_dp, 11.0_dp, 14.0_dp, 11.0_dp])
   call check_drawdown('maxDrawDown_rising', [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp])
   call check('rollFun_mean', roll_fun(close, 7, 'mean'))
   call check('rollFun_sum', roll_fun(close, 7, 'sum'))
   call check('rollFun_min', roll_fun(close, 7, 'min'))
   call check('rollFun_max', roll_fun(close, 7, 'max'))
   call check('rollFun_var', roll_fun(close, 7, 'var'))
   call check('rollFun_sd', roll_fun(close, 7, 'sd'))
   call check('rollFun_padded', roll_fun(close, 4, 'mean', trim = .false.))
   call check('rollFun_na_windows', &
      roll_fun(missing, 3, 'mean', trim = .false.))
   call check('rollFun_na_removed', &
      roll_fun(missing, 3, 'mean', trim = .false., na_rm = .true.))
   call check('rollFun_callback', roll_fun(close, 6, reducer))
   call check('rollFun_unit', roll_fun(close, 1, 'sum'))
   call check('rollVar_default', roll_var(close))
   call check('rollVar_biased', &
      roll_var(close, 7, trim = .false., unbiased = .false.))
   call check('rollVar_na_removed', &
      roll_var(missing, 3, na_rm = .true., unbiased = .false.))
   call check('rollVar_unit', roll_var(close, 1))
   call check('SMA_default', sma(close))
   call check('SMA_custom', sma(close, 13))
   call check('emaTA_default', ema_ta(close))
   call check('emaTA_zero', ema_ta(close, lambda = 0.0_dp))
   call check('emaTA_period', ema_ta(close, lambda = 7.0_dp))
   call check('emaTA_one', ema_ta(close, lambda = 1.0_dp))
   call check('emaTA_decay', ema_ta(close, lambda = 0.35_dp))
   call check('emaTA_startup', ema_ta(close, lambda = 0.3_dp, startup = 4))
   call check('emaTA_oversized', ema_ta(close, lambda = 0.3_dp, startup = 97))
   call check('emaTA_short', ema_ta(close(:6)))
   call check('EWMA_default_startup', ewma(close, 0.4_dp))
   call check('EWMA_custom_startup', ewma(close, 11.0_dp, 6))
   call check('biasTA_default', bias_ta(close))
   call check('biasTA_custom', bias_ta(close, 3))
   call check('rocTA_default', roc_ta(close))
   call check('rocTA_custom', roc_ta(close, 3))
   call check('momTA_default', mom_ta(close))
   call check('momTA_custom', mom_ta(close, 3))
   call check('medpriceTA_default', medprice_ta(high, low))
   call check('vohlTA_default', vohl_ta(high, low))
   call check('vorTA_default', vor_ta(high, low))
   call check('typicalpriceTA_default', typicalprice_ta(high, low, close))
   call check('wcloseTA_default', wclose_ta(high, low, close))
   call check('oscTA_default', osc_ta(close))
   call check('oscTA_custom', osc_ta(close, 4, 9))
   call check('macdTA_default', macd_ta(close))
   call check('macdTA_custom', macd_ta(close, 4, 9))
   call check('cdsTA_default', cds_ta(close))
   call check('cdsTA_custom', cds_ta(close, 4, 9, 3))
   call check('cdoTA_default', cdo_ta(close))
   call check('cdoTA_custom', cdo_ta(close, 4, 9, 3))
   call check('fpkTA_default', fpk_ta(close, high, low))
   call check('fpkTA_custom', fpk_ta(close, high, low, 3))
   call check('wprTA_default', wpr_ta(close, high, low))
   call check('wprTA_custom', wpr_ta(close, high, low, 3))
   call check('fpdTA_default', fpd_ta(close, high, low))
   call check('fpdTA_custom', fpd_ta(close, high, low, 3, 4))
   call check('spdTA_default', spd_ta(close, high, low))
   call check('spdTA_custom', spd_ta(close, high, low, 3, 4, 6))
   call check('apdTA_default', apd_ta(close, high, low))
   call check('apdTA_custom', apd_ta(close, high, low, 3, 4, 6, 7))
   call check('stochasticTA_fast', reshape(stochastic_ta(close, high, low), [2 * n]))
   call check('stochasticTA_slow', &
      reshape(stochastic_ta(close, high, low, 3, 4, 6, slow = .true.), [2 * n]))
   call check('rsiTA_default', rsi_ta(close))
   call check('rsiTA_custom', rsi_ta(close, 4))
   call check('accelTA_default', accel_ta(close))
   call check('accelTA_custom', accel_ta(close, 4))
   call check('adiTA_default', adi_ta(high, low, close, volume))
   call check('adoscillatorTA_default', adoscillator_ta(open, high, low, close))
   call check('bollingerTA_default', reshape(bollinger_ta(close), [3 * n]))
   call check('bollingerTA_custom', &
      reshape(bollinger_ta(close, 7, 1.5_dp), [3 * n]))
   call check('chaikinoTA_default', chaikino_ta(high, low, close, volume))
   call check('chaikinoTA_custom', chaikino_ta(high, low, close, volume, 6, 4))
   call check('chaikinvTA_default_ts', chaikinv_ta(high, low))
   call check('chaikinvTA_custom_ts', &
      chaikinv_ta(high, low, 6, 4))
   call check('garmanklassTA_default', garmanklass_ta(open, high, low, close))
   call check('nviTA_default', nvi_ta(close, volume))
   call check('obvTA_default', obv_ta(close, volume))
   call check('pviTA_default', pvi_ta(close, volume))
   call check('pvtrendTA_default', pvtrend_ta(close, volume))
   call check('williamsadTA_ts', williamsad_ta(high, low, close))
   call check('williamsrTA_default_ts', &
      williamsr_ta(high, low, close))
   call check('williamsrTA_custom_ts', &
      williamsr_ta(high, low, close, 7))
   call check('dailyTA_ema', daily_ta(ohlc, 'ema', [6]))
   call check('dailyTA_bias', daily_ta(ohlc, 'bias', [6]))
   call check('dailyTA_medprice', daily_ta(ohlc, 'medprice', [6]))
   call check('dailyTA_typicalprice', daily_ta(ohlc, 'typicalprice', [6]))
   call check('dailyTA_wclose', daily_ta(ohlc, 'wclose', [6]))
   call check('dailyTA_roc', daily_ta(ohlc, 'roc', [6]))
   call check('dailyTA_mom', daily_ta(ohlc, 'mom', [6]))
   call check('dailyTA_osc', daily_ta(ohlc, 'osc', [4, 9]))
   call check('dailyTA_macd', daily_ta(ohlc, 'macd', [4, 9]))
   call check('dailyTA_cds', daily_ta(ohlc, 'cds', [4, 9, 3]))
   call check('dailyTA_cdo', daily_ta(ohlc, 'cdo', [4, 9, 3]))
   call check('dailyTA_vohl', daily_ta(ohlc, 'vohl', [6]))
   call check('dailyTA_vor', daily_ta(ohlc, 'vor', [6]))
   call check('dailyTA_rsi', daily_ta(ohlc, 'rsi', [6]))
   call check('dailyTA_fpk', daily_ta(ohlc, 'fpk', [6]))
   call check('dailyTA_fpd', daily_ta(ohlc, 'fpd', [4, 3]))
   call check('dailyTA_spd', daily_ta(ohlc, 'spd', [4, 3, 6]))
   call check('dailyTA_apd', daily_ta(ohlc, 'apd', [4, 3, 6, 7]))
   call check('dailyTA_wpr', daily_ta(ohlc, 'wpr', [6]))
   call check('dailyTA_select', &
      daily_ta(ohlc, 'ema', [6], select_column = 1))
   call check('hitRate_default', [hit_rate(returns, positions)])
   call check('hitRate_missing', &
      [hit_rate(returns_missing, positions_missing)])
   call trade_signals(positions, signals, indices)
   call check('tradeSignals_values', signals)
   call check('tradeSignals_indices', real(indices, dp))
   times = fixture('input_signal_times')
   call check('tradeLengths', trade_lengths(times))

   do k = 1, size(records)
      if (.not. records(k)%used) then
         write(*, '(a)') 'Unused R fixture: ' // trim(records(k)%name)
         error stop 1
      end if
   end do
   write(*, '(a,i0,a)') 'Passed ', comparisons, ' genuine R fixture comparisons.'
contains
   subroutine load_references()
      integer :: unit, status, count_records, i, j, count_values
      open(newunit = unit, file = 'test/reference/fixtures.txt', action = 'read', status = 'old', iostat = status)
      if (status /= 0) error stop 'Run fpm test from the fTrading package root.'
      read(unit, *, iostat = status) count_records
      if (status /= 0 .or. count_records < 1) error stop 'Invalid R fixture header.'
      allocate(records(count_records))
      do i = 1, count_records
         read(unit, '(a)', iostat = status) records(i)%name
         if (status /= 0) error stop 'Missing R fixture name.'
         read(unit, *, iostat = status) count_values
         if (status /= 0 .or. count_values < 0) error stop 'Invalid R fixture length.'
         allocate(records(i)%values(count_values))
         do j = 1, count_values
            read(unit, *, iostat = status) records(i)%values(j)
            if (status /= 0) error stop 'Invalid R fixture scalar.'
         end do
      end do
      close(unit)
   end subroutine load_references

   function fixture(name) result(values)
      character(len=*), intent(in) :: name !! Unique record label in the checked-in R oracle data.
      real(dp), allocatable :: values(:)
      integer :: i
      do i = 1, size(records)
         if (trim(records(i)%name) == name) then
            records(i)%used = .true.
            values = records(i)%values
            return
         end if
      end do
      write(*, '(a)') 'R fixture was not found: ' // name
      error stop 1
   end function fixture

   subroutine check(name, actual)
      character(len=*), intent(in) :: name !! Label selecting the unchanged R calculation to compare.
      real(dp), intent(in) :: actual(:) !! Fortran result, flattened column-major for matrix outputs.
      real(dp), allocatable :: expected(:)
      real(dp), parameter :: tolerance = 3.0e-11_dp
      logical :: matches
      integer :: i
      expected = fixture(name)
      if (size(actual) /= size(expected)) then
         write(*, '(a,2(1x,i0))') 'Shape mismatch: ' // name, size(actual), size(expected)
         error stop 1
      end if
      do i = 1, size(expected)
         if (ieee_is_nan(expected(i))) then
            matches = ieee_is_nan(actual(i))
         else if (.not. ieee_is_finite(expected(i))) then
            matches = .not. ieee_is_nan(actual(i)) .and. .not. ieee_is_finite(actual(i))
            if (matches) matches = sign(1.0_dp, actual(i)) == sign(1.0_dp, expected(i))
         else
            matches = ieee_is_finite(actual(i))
            if (matches) matches = abs(actual(i) - expected(i)) <= tolerance * max(1.0_dp, abs(expected(i)))
         end if
         if (.not. matches) then
            write(*, '(a,1x,i0,2(1x,es24.16))') 'R mismatch: ' // name, i, actual(i), expected(i)
            error stop 1
         end if
      end do
      comparisons = comparisons + 1
   end subroutine check

   subroutine check_drawdown(name, values)
      character(len=*), intent(in) :: name !! Prefix of the drawdown magnitude and tied-index R fixtures.
      real(dp), intent(in) :: values(:) !! Chronological finite levels passed unchanged to max_draw_down.
      real(dp) :: drawdown
      integer, allocatable :: from(:), to(:)
      call max_draw_down(values, drawdown, from, to)
      call check(name // '_value', [drawdown])
      call check(name // '_from', real(from, dp))
      call check(name // '_to', real(to, dp))
   end subroutine check_drawdown

   pure function square_sum(window) result(value)
      real(dp), intent(in) :: window(:) !! Complete rolling sample whose squared observations are summed.
      real(dp) :: value
      value = sum(window ** 2)
   end function square_sum
end program test_r_fixtures
