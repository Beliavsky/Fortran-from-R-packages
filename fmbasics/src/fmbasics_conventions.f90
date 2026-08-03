! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_conventions
   use fmbasics_kinds, only : FM_OK, FM_INVALID_ARGUMENT
   use fmbasics_dates
   implicit none
   private

   type, public :: currency_t
      character(len=3) :: iso = '   '
      type(calendar_t) :: calendar
   end type currency_t

   type, public :: currency_pair_t
      type(currency_t) :: base_ccy
      type(currency_t) :: quote_ccy
      type(calendar_t) :: calendar
   end type currency_pair_t

   type, public :: index_t
      character(len=16) :: name = ''
      type(currency_t) :: currency
      type(period_t) :: tenor
      type(period_t) :: spot_lag
      type(calendar_t) :: calendar
      type(calendar_t) :: pfc_calendar
      character(len=12) :: day_basis = 'act/365'
      character(len=4) :: day_convention = 'f'
      logical :: is_eom = .false.
      logical :: is_cash = .false.
   end type index_t

   public :: currency, currency_pair, pair_iso, invert
   public :: aud, eur, gbp, jpy, nzd, usd, chf, hkd, nok
   public :: audusd, eurusd, nzdusd, gbpusd, usdjpy, gbpjpy
   public :: eurgbp, audnzd, eurchf, usdchf, usdhkd, eurnok, usdnok
   public :: is_t1, to_spot, to_spot_next, to_today, to_tomorrow
   public :: to_forward, to_fx_value
   public :: ibor_index, cash_index, to_reset, to_value, to_maturity
   public :: audbbsw, audbbsw1b, euribor, gbplibor, jpylibor, jpytibor
   public :: nzdbkbm, usdlibor, chflibor, hkdhibor, noknibor
   public :: aonia, eonia, sonia, tonar, nziona, fedfunds, chftois, honix

   interface to_spot
      module procedure to_spot_scalar
      module procedure to_spot_vector
   end interface to_spot

   interface to_spot_next
      module procedure to_spot_next_scalar
      module procedure to_spot_next_vector
   end interface to_spot_next

   interface to_today
      module procedure to_today_scalar
      module procedure to_today_vector
   end interface to_today

   interface to_tomorrow
      module procedure to_tomorrow_scalar
      module procedure to_tomorrow_vector
   end interface to_tomorrow

   interface to_forward
      module procedure to_forward_scalar
      module procedure to_forward_vector
   end interface to_forward

   interface to_fx_value
      module procedure to_fx_value_named_scalar
      module procedure to_fx_value_named_vector
      module procedure to_fx_value_period_scalar
      module procedure to_fx_value_period_vector
   end interface to_fx_value

   interface to_reset
      module procedure to_reset_scalar
      module procedure to_reset_vector
   end interface to_reset

   interface to_value
      module procedure to_value_scalar
      module procedure to_value_vector
   end interface to_value

   interface to_maturity
      module procedure to_maturity_scalar
      module procedure to_maturity_vector
   end interface to_maturity

contains

   function currency(iso, cal) result(value)
      character(len=*), intent(in) :: iso
      type(calendar_t), intent(in) :: cal
      type(currency_t) :: value
      value%iso = upper3(iso)
      value%calendar = cal
   end function currency

   function aud() result(value)
      type(currency_t) :: value
      value = currency('AUD', calendar('AUSY'))
   end function aud

   function eur() result(value)
      type(currency_t) :: value
      value = currency('EUR', calendar('EUTA'))
   end function eur

   function gbp() result(value)
      type(currency_t) :: value
      value = currency('GBP', calendar('GBLO'))
   end function gbp

   function jpy() result(value)
      type(currency_t) :: value
      value = currency('JPY', calendar('JPTO'))
   end function jpy

   function nzd() result(value)
      type(currency_t) :: value
      value = currency('NZD', calendar([character(len=8) :: 'NZAU', 'NZWE']))
   end function nzd

   function usd() result(value)
      type(currency_t) :: value
      value = currency('USD', calendar('USNY'))
   end function usd

   function chf() result(value)
      type(currency_t) :: value
      value = currency('CHF', calendar('CHZH'))
   end function chf

   function hkd() result(value)
      type(currency_t) :: value
      value = currency('HKD', calendar('HKHK'))
   end function hkd

   function nok() result(value)
      type(currency_t) :: value
      value = currency('NOK', calendar('NOOS'))
   end function nok

   function currency_pair(base_ccy, quote_ccy, cal) result(value)
      type(currency_t), intent(in) :: base_ccy, quote_ccy
      type(calendar_t), intent(in), optional :: cal
      type(currency_pair_t) :: value
      type(calendar_t) :: combined
      value%base_ccy = base_ccy
      value%quote_ccy = quote_ccy
      if (present(cal)) then
         combined = cal
      else
         combined = joint_calendar(base_ccy%calendar, quote_ccy%calendar)
      end if
      value%calendar = remove_calendar(combined, 'USNY')
   end function currency_pair

   pure function pair_iso(pair) result(value)
      type(currency_pair_t), intent(in) :: pair
      character(len=6) :: value
      value = pair%base_ccy%iso // pair%quote_ccy%iso
   end function pair_iso

   function invert(pair) result(value)
      type(currency_pair_t), intent(in) :: pair
      type(currency_pair_t) :: value
      value = currency_pair(pair%quote_ccy, pair%base_ccy, pair%calendar)
   end function invert

   function audusd() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(aud(), usd())
   end function audusd
   function eurusd() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(eur(), usd())
   end function eurusd
   function nzdusd() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(nzd(), usd())
   end function nzdusd
   function gbpusd() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(gbp(), usd())
   end function gbpusd
   function usdjpy() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(usd(), jpy())
   end function usdjpy
   function gbpjpy() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(gbp(), jpy())
   end function gbpjpy
   function eurgbp() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(eur(), gbp())
   end function eurgbp
   function audnzd() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(aud(), nzd())
   end function audnzd
   function eurchf() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(eur(), chf())
   end function eurchf
   function usdchf() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(usd(), chf())
   end function usdchf
   function usdhkd() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(usd(), hkd())
   end function usdhkd
   function eurnok() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(eur(), nok())
   end function eurnok
   function usdnok() result(value)
      type(currency_pair_t) :: value
      value = currency_pair(usd(), nok())
   end function usdnok

   pure logical function is_t1(pair) result(value)
      type(currency_pair_t), intent(in) :: pair
      character(len=6), parameter :: pairs(12) = [character(len=6) :: &
         'USDCAD', 'USDTRY', 'USDPHP', 'USDRUB', 'USDKZT', 'USDPKR', &
         'CADUSD', 'TRYUSD', 'PHPUSD', 'RUBUSD', 'KZTUSD', 'PKRUSD']
      value = any(pair_iso(pair) == pairs)
   end function is_t1

   integer function to_spot_scalar(date, pair) result(value)
      integer, intent(in) :: date
      type(currency_pair_t), intent(in) :: pair
      type(period_t) :: lag
      type(calendar_t) :: with_usny
      lag = days_period(merge(1, 2, is_t1(pair)))
      value = shift_date(date, lag, 'f', pair%calendar, .false.)
      with_usny = add_calendar(pair%calendar, 'USNY')
      value = shift_date(value, days_period(0), 'f', with_usny, .false.)
   end function to_spot_scalar

   function to_spot_vector(date, pair) result(value)
      integer, intent(in) :: date(:)
      type(currency_pair_t), intent(in) :: pair
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_spot_scalar(date(i), pair)
      end do
   end function to_spot_vector

   integer function to_spot_next_scalar(date, pair) result(value)
      integer, intent(in) :: date
      type(currency_pair_t), intent(in) :: pair
      type(calendar_t) :: with_usny
      with_usny = add_calendar(pair%calendar, 'USNY')
      value = shift_date(to_spot_scalar(date, pair), days_period(1), 'f', with_usny, .false.)
   end function to_spot_next_scalar

   function to_spot_next_vector(date, pair) result(value)
      integer, intent(in) :: date(:)
      type(currency_pair_t), intent(in) :: pair
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_spot_next_scalar(date(i), pair)
      end do
   end function to_spot_next_vector

   integer function to_forward_scalar(date, tenor, pair) result(value)
      integer, intent(in) :: date
      type(period_t), intent(in) :: tenor
      type(currency_pair_t), intent(in) :: pair
      type(calendar_t) :: with_usny
      with_usny = add_calendar(pair%calendar, 'USNY')
      value = shift_date(to_spot_scalar(date, pair), tenor, 'f', with_usny, .true.)
   end function to_forward_scalar

   function to_forward_vector(date, tenor, pair) result(value)
      integer, intent(in) :: date(:)
      type(period_t), intent(in) :: tenor
      type(currency_pair_t), intent(in) :: pair
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_forward_scalar(date(i), tenor, pair)
      end do
   end function to_forward_vector

   integer function to_today_scalar(date, pair) result(value)
      integer, intent(in) :: date
      type(currency_pair_t), intent(in) :: pair
      type(calendar_t) :: with_usny
      with_usny = add_calendar(pair%calendar, 'USNY')
      if (is_good_day(date, with_usny)) then
         value = date
      else
         value = DATE_NA
      end if
   end function to_today_scalar

   function to_today_vector(date, pair) result(value)
      integer, intent(in) :: date(:)
      type(currency_pair_t), intent(in) :: pair
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_today_scalar(date(i), pair)
      end do
   end function to_today_vector

   integer function to_tomorrow_scalar(date, pair) result(value)
      integer, intent(in) :: date
      type(currency_pair_t), intent(in) :: pair
      type(calendar_t) :: with_usny
      with_usny = add_calendar(pair%calendar, 'USNY')
      value = shift_date(date, days_period(1), 'f', with_usny, .false.)
      if (value >= to_spot_scalar(date, pair)) value = DATE_NA
   end function to_tomorrow_scalar

   function to_tomorrow_vector(date, pair) result(value)
      integer, intent(in) :: date(:)
      type(currency_pair_t), intent(in) :: pair
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_tomorrow_scalar(date(i), pair)
      end do
   end function to_tomorrow_vector

   integer function to_fx_value_named_scalar(date, tenor, pair, status) result(value)
      integer, intent(in) :: date
      character(len=*), intent(in) :: tenor
      type(currency_pair_t), intent(in) :: pair
      integer, intent(out), optional :: status
      select case (trim(lower(tenor)))
      case ('spot')
         value = to_spot_scalar(date, pair)
      case ('spot_next')
         value = to_spot_next_scalar(date, pair)
      case ('today')
         value = to_today_scalar(date, pair)
      case ('tomorrow')
         value = to_tomorrow_scalar(date, pair)
      case default
         value = DATE_NA
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end select
      if (present(status)) status = FM_OK
   end function to_fx_value_named_scalar

   function to_fx_value_named_vector(date, tenor, pair, status) result(value)
      integer, intent(in) :: date(:)
      character(len=*), intent(in) :: tenor
      type(currency_pair_t), intent(in) :: pair
      integer, intent(out), optional :: status
      integer, allocatable :: value(:)
      integer :: i, stat_i
      allocate(value(size(date)))
      stat_i = FM_OK
      do i = 1, size(date)
         value(i) = to_fx_value_named_scalar(date(i), tenor, pair, stat_i)
      end do
      if (present(status)) status = stat_i
   end function to_fx_value_named_vector

   integer function to_fx_value_period_scalar(date, tenor, pair) result(value)
      integer, intent(in) :: date
      type(period_t), intent(in) :: tenor
      type(currency_pair_t), intent(in) :: pair
      value = to_forward_scalar(date, tenor, pair)
   end function to_fx_value_period_scalar

   function to_fx_value_period_vector(date, tenor, pair) result(value)
      integer, intent(in) :: date(:)
      type(period_t), intent(in) :: tenor
      type(currency_pair_t), intent(in) :: pair
      integer, allocatable :: value(:)
      value = to_forward_vector(date, tenor, pair)
   end function to_fx_value_period_vector

   function ibor_index(name, ccy, tenor, spot_lag, cal, day_basis, convention, is_eom) result(value)
      character(len=*), intent(in) :: name, day_basis, convention
      type(currency_t), intent(in) :: ccy
      type(period_t), intent(in) :: tenor, spot_lag
      type(calendar_t), intent(in) :: cal
      logical, intent(in) :: is_eom
      type(index_t) :: value
      value%name = name
      value%currency = ccy
      value%tenor = tenor
      value%spot_lag = spot_lag
      value%calendar = cal
      value%pfc_calendar = ccy%calendar
      value%day_basis = day_basis
      value%day_convention = convention
      value%is_eom = is_eom
      value%is_cash = .false.
   end function ibor_index

   function cash_index(name, ccy, spot_lag, cal, day_basis, convention) result(value)
      character(len=*), intent(in) :: name, day_basis, convention
      type(currency_t), intent(in) :: ccy
      type(period_t), intent(in) :: spot_lag
      type(calendar_t), intent(in) :: cal
      type(index_t) :: value
      value = ibor_index(name, ccy, days_period(1), spot_lag, cal, day_basis, convention, .false.)
      value%is_cash = .true.
   end function cash_index

   pure integer function to_reset_scalar(date, index) result(value)
      integer, intent(in) :: date
      type(index_t), intent(in) :: index
      type(period_t) :: lag
      lag = index%spot_lag
      lag%days = -lag%days
      lag%months = -lag%months
      value = shift_date(date, lag, index%day_convention, index%calendar, index%is_eom)
   end function to_reset_scalar

   function to_reset_vector(date, index) result(value)
      integer, intent(in) :: date(:)
      type(index_t), intent(in) :: index
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_reset_scalar(date(i), index)
      end do
   end function to_reset_vector

   pure integer function to_value_scalar(date, index) result(value)
      integer, intent(in) :: date
      type(index_t), intent(in) :: index
      value = shift_date(date, index%spot_lag, index%day_convention, index%calendar, index%is_eom)
   end function to_value_scalar

   function to_value_vector(date, index) result(value)
      integer, intent(in) :: date(:)
      type(index_t), intent(in) :: index
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_value_scalar(date(i), index)
      end do
   end function to_value_vector

   integer function to_maturity_scalar(date, index) result(value)
      integer, intent(in) :: date
      type(index_t), intent(in) :: index
      type(calendar_t) :: combined
      combined = joint_calendar(index%pfc_calendar, index%calendar)
      value = shift_date(date, index%tenor, index%day_convention, combined, index%is_eom)
   end function to_maturity_scalar

   function to_maturity_vector(date, index) result(value)
      integer, intent(in) :: date(:)
      type(index_t), intent(in) :: index
      integer, allocatable :: value(:)
      integer :: i
      allocate(value(size(date)))
      do i = 1, size(date)
         value(i) = to_maturity_scalar(date(i), index)
      end do
   end function to_maturity_vector

   function audbbsw(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      value = ibor_index('BBSW', aud(), tenor, days_period(0), calendar('AUSY'), &
                         'act/365', 'ms', .false.)
   end function audbbsw

   function audbbsw1b(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      value = ibor_index('BBSW1b', aud(), tenor, days_period(1), calendar('AUSY'), &
                         'act/365', 'ms', .false.)
   end function audbbsw1b

   function euribor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      value = ibor_index('EURIBOR', eur(), tenor, days_period(2), calendar('EUTA'), &
                         'act/360', 'mf', .true.)
   end function euribor

   function gbplibor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      character(len=4) :: convention
      convention = merge('f   ', 'mf  ', tenor%months == 0)
      value = ibor_index('LIBOR', gbp(), tenor, days_period(0), calendar('GBLO'), &
                         'act/365', convention, .true.)
   end function gbplibor

   function jpylibor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      integer :: lag
      character(len=4) :: convention
      lag = 2
      convention = 'mf'
      if (tenor%months == 0) then
         convention = 'f'
         if (tenor%days == 1) lag = 0
      end if
      value = ibor_index('LIBOR', jpy(), tenor, days_period(lag), calendar('GBLO'), &
                         'act/360', convention, .true.)
   end function jpylibor

   function jpytibor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      integer :: lag
      character(len=4) :: convention
      lag = 2
      convention = 'mf'
      if (tenor%months == 0) then
         convention = 'f'
         if (tenor%days == 1) lag = 0
      end if
      value = ibor_index('TIBOR', jpy(), tenor, days_period(lag), calendar('JPTO'), &
                         'act/365', convention, .false.)
   end function jpytibor

   function nzdbkbm(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      value = ibor_index('BKBM', nzd(), tenor, days_period(0), &
                         calendar([character(len=8) :: 'NZAU', 'NZWE']), &
                         'act/365', 'mf', .false.)
   end function nzdbkbm

   function usdlibor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      type(calendar_t) :: cal
      integer :: lag
      character(len=4) :: convention
      lag = 2
      if (tenor%months == 0) then
         convention = 'f'
         if (tenor%days == 1) then
            lag = 0
            cal = calendar([character(len=8) :: 'USNY', 'GBLO'])
         else
            cal = calendar('GBLO')
         end if
      else
         convention = 'mf'
         cal = calendar('GBLO')
      end if
      value = ibor_index('LIBOR', usd(), tenor, days_period(lag), cal, &
                         'act/360', convention, .true.)
   end function usdlibor

   function chflibor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      integer :: lag
      character(len=4) :: convention
      if (tenor%months == 0) then
         lag = merge(2, 0, tenor%days == 1)
         convention = 'f'
      else
         lag = 2
         convention = 'mf'
      end if
      value = ibor_index('LIBOR', chf(), tenor, days_period(lag), calendar('GBLO'), &
                         'act/360', convention, .true.)
   end function chflibor

   function hkdhibor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      character(len=4) :: convention
      convention = merge('f   ', 'mf  ', tenor%months == 0)
      value = ibor_index('HIBOR', hkd(), tenor, days_period(0), calendar('HKHK'), &
                         'act/365', convention, .false.)
   end function hkdhibor

   function noknibor(tenor) result(value)
      type(period_t), intent(in) :: tenor
      type(index_t) :: value
      integer :: lag
      character(len=4) :: convention
      if (tenor%months == 0) then
         lag = merge(1, 2, tenor%days == 1)
         convention = 'f'
      else
         lag = 2
         convention = 'mf'
      end if
      value = ibor_index('NIBOR', nok(), tenor, days_period(lag), calendar('NOOS'), &
                         'act/360', convention, .false.)
   end function noknibor

   function aonia() result(value)
      type(index_t) :: value
      value = cash_index('AONIA', aud(), days_period(0), calendar('AUSY'), 'act/365', 'f')
   end function aonia
   function eonia() result(value)
      type(index_t) :: value
      value = cash_index('EONIA', eur(), days_period(0), calendar('EUTA'), 'act/360', 'f')
   end function eonia
   function sonia() result(value)
      type(index_t) :: value
      value = cash_index('SONIA', gbp(), days_period(0), calendar('GBLO'), 'act/365', 'f')
   end function sonia
   function tonar() result(value)
      type(index_t) :: value
      value = cash_index('TONAR', jpy(), days_period(0), calendar('JPTO'), 'act/365', 'f')
   end function tonar
   function nziona() result(value)
      type(index_t) :: value
      value = cash_index('NZIONA', nzd(), days_period(0), &
                         calendar([character(len=8) :: 'NZAU', 'NZWE']), 'act/365', 'f')
   end function nziona
   function fedfunds() result(value)
      type(index_t) :: value
      value = cash_index('FedFunds', usd(), days_period(0), calendar('USNY'), 'act/360', 'f')
   end function fedfunds
   function chftois() result(value)
      type(index_t) :: value
      value = cash_index('CHFTOIS', chf(), days_period(1), calendar('CHZH'), 'act/360', 'f')
   end function chftois
   function honix() result(value)
      type(index_t) :: value
      value = cash_index('HONIX', hkd(), days_period(0), calendar('HKHK'), 'act/365', 'f')
   end function honix

   function add_calendar(cal, code) result(value)
      type(calendar_t), intent(in) :: cal
      character(len=*), intent(in) :: code
      type(calendar_t) :: value
      if (calendar_contains(cal, code)) then
         value = cal
      else
         value = joint_calendar(calendar(code), cal)
      end if
   end function add_calendar

   function remove_calendar(cal, code) result(value)
      type(calendar_t), intent(in) :: cal
      character(len=*), intent(in) :: code
      type(calendar_t) :: value
      integer :: i, n
      character(len=8), allocatable :: tmp(:)
      allocate(tmp(cal%size()))
      n = 0
      do i = 1, cal%size()
         if (trim(cal%code(i)) /= trim(upper8(code))) then
            n = n + 1
            tmp(n) = cal%code(i)
         end if
      end do
      allocate(value%code(n))
      if (n > 0) value%code = tmp(:n)
   end function remove_calendar

   pure function upper3(text) result(out)
      character(len=*), intent(in) :: text
      character(len=3) :: out
      integer :: i, k
      out = '   '
      do i = 1, min(3, len_trim(text))
         k = iachar(text(i:i))
         if (k >= iachar('a') .and. k <= iachar('z')) then
            out(i:i) = achar(k - 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function upper3

   pure function upper8(text) result(out)
      character(len=*), intent(in) :: text
      character(len=8) :: out
      integer :: i, k
      out = '        '
      do i = 1, min(8, len_trim(text))
         k = iachar(text(i:i))
         if (k >= iachar('a') .and. k <= iachar('z')) then
            out(i:i) = achar(k - 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function upper8

   pure function lower(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, k
      do i = 1, len(text)
         k = iachar(text(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) then
            out(i:i) = achar(k + 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lower

end module fmbasics_conventions
