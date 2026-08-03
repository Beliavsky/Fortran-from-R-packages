! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_cashflows
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument, qbc_size_mismatch
   use qbc_dates, only : qbc_date, add_months, add_years, add_days, adjust_business_day, year_fraction, &
      days_between, operator(>), operator(<=), operator(==)
   use qbc_types, only : qbc_coupon_schedule, qbc_schedule_short_first, qbc_schedule_long_first, &
      qbc_schedule_short_last, qbc_schedule_long_last, qbc_asset_tes
   implicit none
   private
   public :: coupon_dates, coupon_dates_from_years, coupon_cashflows, accrued_interest

contains

   subroutine coupon_dates(maturity, analysis_date, frequency, schedule, convention, result, &
                           trade_date, settlement_lag, status)
      type(qbc_date), intent(in) :: maturity, analysis_date
      integer, intent(in) :: frequency
      integer, intent(in), optional :: schedule
      character(len=*), intent(in), optional :: convention
      type(qbc_coupon_schedule), intent(out) :: result
      type(qbc_date), intent(in), optional :: trade_date
      integer, intent(in), optional :: settlement_lag
      integer, intent(out), optional :: status
      type(qbc_date), allocatable :: work(:), tmp(:)
      type(qbc_date) :: start_date, settlement, d
      integer :: months_per_coupon, nmax, n, i, sched, lag, total_months, remainder, first_month
      character(len=2) :: conv
      integer :: st

      st = qbc_success
      sched = qbc_schedule_short_first
      if (present(schedule)) sched = schedule
      conv = 'F '
      if (present(convention)) conv = adjustl(convention)
      lag = 0
      if (present(settlement_lag)) lag = settlement_lag

      if (allocated(result%dates)) deallocate(result%dates)
      if (allocated(result%effective_dates)) deallocate(result%effective_dates)
      if (frequency < 0) then
         st = qbc_invalid_argument
         allocate(result%dates(0), result%effective_dates(0))
         if (present(status)) status = st
         return
      end if
      if (frequency > 0) then
         if (mod(12, frequency) /= 0) then
            st = qbc_invalid_argument
            allocate(result%dates(0), result%effective_dates(0))
            if (present(status)) status = st
            return
         end if
      end if
      if (frequency == 0) then
         if (maturity > analysis_date) then
            allocate(result%dates(1), result%effective_dates(1))
            result%dates(1) = maturity
            result%effective_dates(1) = adjust_business_day(maturity, conv)
         else
            allocate(result%dates(0), result%effective_dates(0))
         end if
         if (present(status)) status = st
         return
      end if

      months_per_coupon = 12 / frequency
      if (present(trade_date)) then
         start_date = trade_date
         settlement = trade_date
         do i = 1, lag
            settlement = add_days(settlement, 1)
            do while (.not. (adjust_business_day(settlement, 'F') == settlement))
               settlement = add_days(settlement, 1)
            end do
         end do
         start_date = settlement
         total_months = (maturity%year - start_date%year) * 12 + maturity%month - start_date%month
         if (maturity%day < start_date%day) total_months = total_months - 1
         total_months = max(total_months, 0)
         remainder = mod(total_months, months_per_coupon)
         if (remainder == 0) sched = qbc_schedule_long_first
         select case (sched)
         case (qbc_schedule_short_first)
            first_month = remainder
            if (first_month == 0) first_month = months_per_coupon
            nmax = max(1, (total_months - first_month) / months_per_coupon + 1)
            allocate(work(nmax + 1))
            n = 0
            i = first_month
            do while (i <= total_months)
               d = add_months(start_date, i)
               if (d > analysis_date .and. d > start_date .and. d <= maturity) then
                  n = n + 1; work(n) = d
               end if
               i = i + months_per_coupon
            end do
         case (qbc_schedule_long_first)
            first_month = remainder + months_per_coupon
            if (remainder == 0) first_month = months_per_coupon
            nmax = max(1, total_months / months_per_coupon + 1)
            allocate(work(nmax + 1))
            n = 0
            i = first_month
            do while (i <= total_months)
               d = add_months(start_date, i)
               if (d > analysis_date .and. d > start_date .and. d <= maturity) then
                  n = n + 1; work(n) = d
               end if
               i = i + months_per_coupon
            end do
         case (qbc_schedule_short_last)
            nmax = max(1, total_months / months_per_coupon + 2)
            allocate(work(nmax))
            n = 0
            i = months_per_coupon
            do while (i < total_months)
               d = add_months(start_date, i)
               if (d > analysis_date .and. d > start_date) then
                  n = n + 1; work(n) = d
               end if
               i = i + months_per_coupon
            end do
            if (maturity > analysis_date) then
               n = n + 1; work(n) = maturity
            end if
         case (qbc_schedule_long_last)
            nmax = max(1, total_months / months_per_coupon + 1)
            allocate(work(nmax))
            n = 0
            i = months_per_coupon
            do while (i <= total_months - remainder - months_per_coupon)
               d = add_months(start_date, i)
               if (d > analysis_date .and. d > start_date) then
                  n = n + 1; work(n) = d
               end if
               i = i + months_per_coupon
            end do
            if (maturity > analysis_date) then
               n = n + 1; work(n) = maturity
            end if
         case default
            st = qbc_invalid_argument
            allocate(work(0)); n = 0
         end select
      else
         nmax = max(1, ((maturity%year - analysis_date%year + 2) * 12) / months_per_coupon + 2)
         allocate(work(nmax))
         n = 0
         d = maturity
         do while (d > analysis_date)
            n = n + 1
            if (n > size(work)) then
               allocate(tmp(2 * size(work)))
               tmp(1:n-1) = work(1:n-1)
               call move_alloc(tmp, work)
            end if
            work(n) = d
            d = add_months(d, -months_per_coupon)
         end do
         if (n > 1) work(1:n) = work(n:1:-1)
      end if

      allocate(result%dates(n), result%effective_dates(n))
      if (n > 0) result%dates = work(1:n)
      do i = 1, n
         result%effective_dates(i) = adjust_business_day(result%dates(i), conv)
      end do
      if (allocated(work)) deallocate(work)
      if (present(status)) status = st
   end subroutine coupon_dates

   subroutine coupon_dates_from_years(years_to_maturity, analysis_date, frequency, result, &
                                      schedule, convention, trade_date, status)
      real(dp), intent(in) :: years_to_maturity
      type(qbc_date), intent(in) :: analysis_date
      integer, intent(in) :: frequency
      type(qbc_coupon_schedule), intent(out) :: result
      integer, intent(in), optional :: schedule
      character(len=*), intent(in), optional :: convention
      type(qbc_date), intent(in), optional :: trade_date
      integer, intent(out), optional :: status
      type(qbc_date) :: base, maturity
      integer :: st
      base = analysis_date
      if (present(trade_date)) base = trade_date
      maturity = add_months(base, nint(12.0_dp * years_to_maturity))
      call coupon_dates(maturity, analysis_date, frequency, schedule, convention, result, &
                        trade_date=base, status=st)
      if (present(status)) status = st
   end subroutine coupon_dates_from_years

   subroutine coupon_cashflows(dates, coupon_rates, principal, daycount, cashflows, &
                               previous_date, tes_style, status)
      type(qbc_date), intent(in) :: dates(:)
      real(dp), intent(in) :: coupon_rates(:)
      real(dp), intent(in) :: principal
      character(len=*), intent(in) :: daycount
      real(dp), allocatable, intent(out) :: cashflows(:)
      type(qbc_date), intent(in), optional :: previous_date
      logical, intent(in), optional :: tes_style
      integer, intent(out), optional :: status
      type(qbc_date) :: prev
      integer :: i, st
      logical :: tes
      st = qbc_success
      tes = .false.
      if (present(tes_style)) tes = tes_style
      if (size(coupon_rates) /= 1 .and. size(coupon_rates) /= size(dates)) then
         allocate(cashflows(0)); st = qbc_size_mismatch
         if (present(status)) status = st
         return
      end if
      allocate(cashflows(size(dates)))
      cashflows = 0.0_dp
      if (size(dates) == 0) then
         if (present(status)) status = st
         return
      end if
      if (tes) then
         cashflows = principal * coupon_rates(1)
      else
         if (present(previous_date)) then
            prev = previous_date
         else if (size(dates) >= 2) then
            prev = add_months(dates(1), -((dates(2)%year - dates(1)%year) * 12 + dates(2)%month - dates(1)%month))
         else
            prev = add_years(dates(1), -1)
         end if
         do i = 1, size(dates)
            if (size(coupon_rates) == 1) then
               cashflows(i) = coupon_rates(1) * principal * year_fraction(prev, dates(i), daycount)
            else
               cashflows(i) = coupon_rates(i) * principal * year_fraction(prev, dates(i), daycount)
            end if
            prev = dates(i)
         end do
      end if
      cashflows(size(cashflows)) = cashflows(size(cashflows)) + principal
      if (present(status)) status = st
   end subroutine coupon_cashflows

   real(dp) function accrued_interest(previous_coupon, next_coupon, analysis_date, &
                                           coupon_rate, principal, daycount) result(value)
      type(qbc_date), intent(in) :: previous_coupon, next_coupon, analysis_date
      real(dp), intent(in) :: coupon_rate, principal
      character(len=*), intent(in) :: daycount
      real(dp) :: elapsed, full
      elapsed = year_fraction(previous_coupon, analysis_date, daycount)
      full = year_fraction(previous_coupon, next_coupon, daycount)
      if (full <= 0.0_dp) then
         value = 0.0_dp
      else
         value = coupon_rate * principal * elapsed
      end if
   end function accrued_interest

end module qbc_cashflows
