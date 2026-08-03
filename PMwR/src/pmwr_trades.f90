module pmwr_trades
   use pmwr_kinds, only : dp
   use pmwr_types, only : journal_type, pl_path_result, pl_summary_result
   implicit none
   private
   public :: make_journal, sort_journal, append_journal
   public :: positions_at, journal_cashflows
   public :: pl_summary, pl_path
   public :: split_trade_runs, scale_trades_to_unity, close_on_first, limit_position
   public :: time_weighted_exposure

contains

   subroutine make_journal(timestamp, amount, price, instrument, journal, account)
      real(dp), intent(in) :: timestamp(:), amount(:), price(:)
      integer, intent(in) :: instrument(:)
      type(journal_type), intent(out) :: journal
      integer, intent(in), optional :: account(:)
      integer :: n
      n = size(amount)
      if (size(timestamp) /= n .or. size(price) /= n .or. size(instrument) /= n) &
         error stop "make_journal: size mismatch"
      journal%n = n
      allocate(journal%timestamp(n), journal%amount(n), journal%price(n), journal%instrument(n), journal%account(n))
      journal%timestamp = timestamp; journal%amount = amount; journal%price = price
      journal%instrument = instrument; journal%account = 1
      if (present(account)) then
         if (size(account) /= n) error stop "make_journal: account size mismatch"
         journal%account = account
      end if
   end subroutine make_journal

   subroutine sort_journal(journal)
      type(journal_type), intent(inout) :: journal
      call stable_sort_full(journal)
   end subroutine sort_journal

   subroutine stable_sort_full(journal)
      type(journal_type), intent(inout) :: journal
      integer :: i, j, inst_tmp, acc_tmp
      real(dp) :: time_tmp, amount_tmp, price_tmp
      type(journal_type) :: src
      if (journal%n < 2) return
      call copy_journal(journal, src)
      ! Rebuild from the original copy with stable insertion ordering.
      journal%timestamp = src%timestamp
      journal%amount = src%amount
      journal%price = src%price
      journal%instrument = src%instrument
      journal%account = src%account
      do i = 2, journal%n
         time_tmp = journal%timestamp(i); amount_tmp = journal%amount(i); price_tmp = journal%price(i)
         inst_tmp = journal%instrument(i); acc_tmp = journal%account(i)
         j = i - 1
         do while (j >= 1 .and. journal%timestamp(j) > time_tmp)
            journal%timestamp(j + 1) = journal%timestamp(j)
            journal%amount(j + 1) = journal%amount(j)
            journal%price(j + 1) = journal%price(j)
            journal%instrument(j + 1) = journal%instrument(j)
            journal%account(j + 1) = journal%account(j)
            j = j - 1
         end do
         journal%timestamp(j + 1) = time_tmp
         journal%amount(j + 1) = amount_tmp
         journal%price(j + 1) = price_tmp
         journal%instrument(j + 1) = inst_tmp
         journal%account(j + 1) = acc_tmp
      end do
   end subroutine stable_sort_full

   subroutine copy_journal(a, b)
      type(journal_type), intent(in) :: a
      type(journal_type), intent(out) :: b
      b%n = a%n
      allocate(b%timestamp(a%n), b%amount(a%n), b%price(a%n), b%instrument(a%n), b%account(a%n))
      b%timestamp = a%timestamp; b%amount = a%amount; b%price = a%price
      b%instrument = a%instrument; b%account = a%account
   end subroutine copy_journal

   subroutine append_journal(a, b, output)
      type(journal_type), intent(in) :: a, b
      type(journal_type), intent(out) :: output
      integer :: n
      n = a%n + b%n
      output%n = n
      allocate(output%timestamp(n), output%amount(n), output%price(n), output%instrument(n), output%account(n))
      if (a%n > 0) then
         output%timestamp(1:a%n) = a%timestamp; output%amount(1:a%n) = a%amount
         output%price(1:a%n) = a%price; output%instrument(1:a%n) = a%instrument
         output%account(1:a%n) = a%account
      end if
      if (b%n > 0) then
         output%timestamp(a%n + 1:n) = b%timestamp; output%amount(a%n + 1:n) = b%amount
         output%price(a%n + 1:n) = b%price; output%instrument(a%n + 1:n) = b%instrument
         output%account(a%n + 1:n) = b%account
      end if
      call stable_sort_full(output)
   end subroutine append_journal

   subroutine positions_at(journal, query_time, n_instruments, position, account_filter)
      type(journal_type), intent(in) :: journal
      real(dp), intent(in) :: query_time(:)
      integer, intent(in) :: n_instruments
      real(dp), allocatable, intent(out) :: position(:,:)
      integer, intent(in), optional :: account_filter
      integer :: i, j, k
      allocate(position(size(query_time), n_instruments)); position = 0.0_dp
      do i = 1, size(query_time)
         do j = 1, journal%n
            if (journal%timestamp(j) <= query_time(i)) then
               if (present(account_filter)) then
                  if (journal%account(j) /= account_filter) cycle
               end if
               k = journal%instrument(j)
               if (k >= 1 .and. k <= n_instruments) position(i, k) = position(i, k) + journal%amount(j)
            end if
         end do
      end do
   end subroutine positions_at

   subroutine journal_cashflows(journal, multiplier, cash_journal, cash_instrument, flip_sign)
      type(journal_type), intent(in) :: journal
      real(dp), intent(in) :: multiplier(:)
      type(journal_type), intent(out) :: cash_journal
      integer, intent(in), optional :: cash_instrument
      logical, intent(in), optional :: flip_sign
      integer :: i, ci
      logical :: flip
      ci = 1; if (present(cash_instrument)) ci = cash_instrument
      flip = .true.; if (present(flip_sign)) flip = flip_sign
      cash_journal%n = journal%n
      allocate(cash_journal%timestamp(journal%n), cash_journal%amount(journal%n), cash_journal%price(journal%n), &
               cash_journal%instrument(journal%n), cash_journal%account(journal%n))
      cash_journal%timestamp = journal%timestamp
      cash_journal%price = 1.0_dp
      cash_journal%instrument = ci
      cash_journal%account = journal%account
      do i = 1, journal%n
         if (journal%instrument(i) < 1 .or. journal%instrument(i) > size(multiplier)) &
            error stop "journal_cashflows: invalid instrument"
         cash_journal%amount(i) = journal%amount(i) * journal%price(i) * multiplier(journal%instrument(i))
         if (flip) cash_journal%amount(i) = -cash_journal%amount(i)
      end do
   end subroutine journal_cashflows

   subroutine pl_summary(amount, price, result, tolerance)
      real(dp), intent(in) :: amount(:), price(:)
      type(pl_summary_result), intent(out) :: result
      real(dp), intent(in), optional :: tolerance
      real(dp) :: tol, buys, sells
      logical, allocatable :: buy(:)
      if (size(amount) /= size(price)) error stop "pl_summary: size mismatch"
      tol = 1.0e-10_dp; if (present(tolerance)) tol = tolerance
      allocate(buy(size(amount))); buy = amount > 0.0_dp
      result%closed = abs(sum(amount)) <= tol
      result%volume = sum(abs(amount))
      if (result%closed) result%total_pl = -sum(amount * price)
      buys = sum(amount, mask=buy)
      sells = sum(amount, mask=.not. buy)
      if (abs(buys) > tol) result%average_buy = sum(price * amount, mask=buy) / buys
      if (abs(sells) > tol) result%average_sell = sum(price * amount, mask=.not. buy) / sells
   end subroutine pl_summary

   subroutine pl_path(amount, price, result)
      real(dp), intent(in) :: amount(:), price(:)
      type(pl_path_result), intent(out) :: result
      real(dp) :: prev_pos, pos
      integer :: n, i
      n = size(amount)
      if (size(price) /= n) error stop "pl_path: size mismatch"
      allocate(result%cumulative_position(n), result%average_price(n), result%realized(n), &
               result%unrealized(n), result%total(n), result%volume(n))
      if (n == 0) return
      result%cumulative_position = 0.0_dp; result%average_price = 0.0_dp
      result%realized = 0.0_dp; result%unrealized = 0.0_dp; result%total = 0.0_dp
      result%volume = 0.0_dp
      pos = amount(1)
      result%cumulative_position(1) = pos
      result%average_price(1) = price(1)
      result%volume(1) = abs(amount(1))
      result%unrealized(1) = pos * (price(1) - result%average_price(1))
      result%total(1) = result%realized(1) + result%unrealized(1)
      do i = 2, n
         prev_pos = pos
         pos = prev_pos + amount(i)
         result%cumulative_position(i) = pos
         result%realized(i) = result%realized(i - 1)
         if (pos * prev_pos < 0.0_dp .and. abs(prev_pos) > epsilon(1.0_dp)) then
            result%realized(i) = result%realized(i) + (price(i) - result%average_price(i - 1)) * prev_pos
            result%average_price(i) = price(i)
         else if (abs(pos) > abs(prev_pos)) then
            if (abs(pos) > epsilon(1.0_dp)) then
               result%average_price(i) = (result%average_price(i - 1) * prev_pos + price(i) * amount(i)) / pos
            else
               result%average_price(i) = price(i)
            end if
         else
            result%average_price(i) = result%average_price(i - 1)
            result%realized(i) = result%realized(i) + amount(i) * (result%average_price(i) - price(i))
         end if
         result%unrealized(i) = pos * (price(i) - result%average_price(i))
         result%total(i) = result%realized(i) + result%unrealized(i)
         result%volume(i) = result%volume(i - 1) + abs(amount(i))
      end do
   end subroutine pl_path

   subroutine split_trade_runs(amount, price, timestamp, out_amount, out_price, out_timestamp, first, last)
      real(dp), intent(in) :: amount(:), price(:), timestamp(:)
      real(dp), allocatable, intent(out) :: out_amount(:), out_price(:), out_timestamp(:)
      integer, allocatable, intent(out) :: first(:), last(:)
      real(dp), allocatable :: a(:), p(:), t(:), cum(:)
      integer :: n, i, m, nr

      n = size(amount)
      if (size(price) /= n .or. size(timestamp) /= n) error stop "split_trade_runs: size mismatch"
      allocate(a(2*n), p(2*n), t(2*n), cum(n))
      cum = 0.0_dp
      if (n > 0) then
         cum(1) = amount(1)
         do i = 2, n
            cum(i) = cum(i - 1) + amount(i)
         end do
      end if
      m = 0
      do i = 1, n
         if (i > 1) then
            if (cum(i - 1) * cum(i) < 0.0_dp) then
               m = m + 1; a(m) = -cum(i - 1); p(m) = price(i); t(m) = timestamp(i)
               m = m + 1; a(m) = cum(i); p(m) = price(i); t(m) = timestamp(i)
            else
               m = m + 1; a(m) = amount(i); p(m) = price(i); t(m) = timestamp(i)
            end if
         else
            m = m + 1; a(m) = amount(i); p(m) = price(i); t(m) = timestamp(i)
         end if
      end do
      allocate(out_amount(m), out_price(m), out_timestamp(m))
      out_amount = a(1:m); out_price = p(1:m); out_timestamp = t(1:m)
      nr = 0
      do i = 1, m
         if (abs(sum(out_amount(1:i))) <= 1.0e-12_dp) nr = nr + 1
      end do
      if (m > 0 .and. abs(sum(out_amount)) > 1.0e-12_dp) nr = nr + 1
      allocate(first(nr), last(nr))
      if (nr == 0) return
      nr = 0; first(1) = 1
      do i = 1, m
         if (abs(sum(out_amount(first(nr + 1):i))) <= 1.0e-12_dp) then
            nr = nr + 1; last(nr) = i
            if (i < m) first(nr + 1) = i + 1
         end if
      end do
      if (nr < size(first)) then
         nr = nr + 1; last(nr) = m
      end if
   end subroutine split_trade_runs

   subroutine scale_trades_to_unity(amount, scaled)
      real(dp), intent(in) :: amount(:)
      real(dp), allocatable, intent(out) :: scaled(:)
      real(dp) :: c, max_abs
      integer :: i
      allocate(scaled(size(amount)))
      c = 0.0_dp; max_abs = 0.0_dp
      do i = 1, size(amount)
         c = c + amount(i); max_abs = max(max_abs, abs(c))
      end do
      if (max_abs > 0.0_dp) then
         scaled = amount / max_abs
      else
         scaled = amount
      end if
   end subroutine scale_trades_to_unity

   subroutine close_on_first(amount, closed_amount)
      real(dp), intent(in) :: amount(:)
      real(dp), allocatable, intent(out) :: closed_amount(:)
      real(dp) :: c
      integer :: i
      allocate(closed_amount(size(amount))); closed_amount = amount
      if (size(amount) == 0) return
      c = amount(1)
      do i = 2, size(amount)
         if (amount(1) * amount(i) < 0.0_dp) then
            closed_amount(i) = -c
            if (i < size(amount)) closed_amount(i + 1:) = 0.0_dp
            exit
         end if
         c = c + amount(i)
      end do
   end subroutine close_on_first

   subroutine limit_position(amount, limit_value, limited)
      real(dp), intent(in) :: amount(:), limit_value
      real(dp), allocatable, intent(out) :: limited(:)
      real(dp) :: c, capped, previous
      integer :: i
      allocate(limited(size(amount)))
      c = 0.0_dp; previous = 0.0_dp
      do i = 1, size(amount)
         c = c + amount(i)
         if (amount(1) >= 0.0_dp) then
            capped = min(c, limit_value)
         else
            capped = max(c, -limit_value)
         end if
         limited(i) = capped - previous
         previous = capped
      end do
   end subroutine limit_position

   real(dp) function time_weighted_exposure(amount, timestamp, start_time, end_time, absolute_value) result(ans)
      real(dp), intent(in) :: amount(:), timestamp(:)
      real(dp), intent(in), optional :: start_time, end_time
      logical, intent(in), optional :: absolute_value
      real(dp), allocatable :: t(:), a(:)
      real(dp) :: start_, end_, pos
      logical :: absval
      integer :: n, i

      if (size(amount) /= size(timestamp)) error stop "time_weighted_exposure: size mismatch"
      start_ = minval(timestamp); if (present(start_time)) start_ = start_time
      end_ = maxval(timestamp); if (present(end_time)) end_ = end_time
      absval = .true.; if (present(absolute_value)) absval = absolute_value
      n = size(amount) + 2
      allocate(t(n), a(n)); t(1) = start_; a(1) = 0.0_dp
      t(2:n - 1) = timestamp; a(2:n - 1) = amount
      t(n) = end_; a(n) = 0.0_dp
      ans = 0.0_dp; pos = 0.0_dp
      do i = 1, n - 1
         pos = pos + a(i)
         if (absval) then
            ans = ans + abs(pos) * (t(i + 1) - t(i))
         else
            ans = ans + pos * (t(i + 1) - t(i))
         end if
      end do
      if (end_ > start_) ans = ans / (end_ - start_)
   end function time_weighted_exposure

end module pmwr_trades
