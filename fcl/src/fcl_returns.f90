! MIT License. Copyright (c) 2024 fcl authors.
module fcl_returns
   use fcl_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   implicit none
   private
   public :: return_series_type

   type :: return_series_type
      integer, allocatable :: dates(:)
      real(dp), allocatable :: market_values(:)
      real(dp), allocatable :: pnl(:)
      integer :: status = 0
   contains
      procedure :: initialize => returns_initialize
      procedure :: twrr_daily
      procedure :: twrr_cumulative
      procedure :: cumulative_pnl
      procedure :: dietz_average_capital
      procedure :: dietz
      procedure, private :: range_indices
      procedure, private :: cashflow
   end type return_series_type

contains

   subroutine returns_initialize(self, dates, market_values, pnl)
      class(return_series_type), intent(out) :: self
      integer, intent(in) :: dates(:)
      real(dp), intent(in) :: market_values(:), pnl(:)
      integer, allocatable :: order(:)
      integer :: i, j, n, dmin, dmax, pos

      self%status = 0
      n = size(dates)
      if (n == 0 .or. size(market_values) /= n .or. size(pnl) /= n) then
         self%status = 1
         return
      end if
      allocate(order(n))
      order = [(i, i = 1, n)]
      do i = 2, n
         pos = order(i)
         j = i - 1
         do while (j >= 1 .and. dates(order(j)) > dates(pos))
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = pos
      end do
      do i = 2, n
         if (dates(order(i)) == dates(order(i - 1))) then
            self%status = 2
            return
         end if
      end do
      dmin = dates(order(1))
      dmax = dates(order(n))
      allocate(self%dates(dmax - dmin + 1), self%market_values(dmax - dmin + 1), &
         self%pnl(dmax - dmin + 1))
      self%dates = [(i, i = dmin, dmax)]
      j = 1
      do i = 1, size(self%dates)
         if (j <= n .and. self%dates(i) == dates(order(j))) then
            self%market_values(i) = market_values(order(j))
            self%pnl(i) = pnl(order(j))
            j = j + 1
         else
            self%market_values(i) = self%market_values(i - 1)
            self%pnl(i) = 0.0_dp
         end if
      end do
   end subroutine returns_initialize

   subroutine range_indices(self, from_date, to_date, first, last, status)
      class(return_series_type), intent(in) :: self
      integer, intent(in) :: from_date, to_date
      integer, intent(out) :: first, last, status
      first = 0
      last = 0
      status = 0
      if (from_date > to_date .or. .not. allocated(self%dates)) then
         status = 1
         return
      end if
      first = from_date - self%dates(1) + 1
      last = to_date - self%dates(1) + 1
      if (first < 1 .or. last > size(self%dates)) status = 2
   end subroutine range_indices

   pure function cashflow(self, i) result(cf)
      class(return_series_type), intent(in) :: self
      integer, intent(in) :: i
      real(dp) :: cf
      if (i <= 1 .or. i > size(self%dates)) then
         cf = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         cf = self%market_values(i) - self%market_values(i - 1) - self%pnl(i)
      end if
   end function cashflow

   function twrr_daily(self, from_date, to_date, status) result(out)
      class(return_series_type), intent(in) :: self
      integer, intent(in) :: from_date, to_date
      integer, intent(out), optional :: status
      real(dp), allocatable :: out(:)
      integer :: first, last, stat, i, k
      real(dp) :: cf, denominator, nan
      call self%range_indices(from_date, to_date, first, last, stat)
      if (stat /= 0) then
         allocate(out(0))
         if (present(status)) status = stat
         return
      end if
      allocate(out(last - first + 1))
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = first, last
         k = i - first + 1
         cf = self%cashflow(i)
         if (.not. ieee_is_finite(cf)) then
            out(k) = nan
         else
            denominator = self%market_values(i - 1) + max(cf, 0.0_dp)
            if (abs(denominator) <= tiny(1.0_dp)) then
               out(k) = nan
            else
               out(k) = self%pnl(i) / denominator
            end if
         end if
      end do
      if (present(status)) status = 0
   end function twrr_daily

   function twrr_cumulative(self, from_date, to_date, status) result(out)
      class(return_series_type), intent(in) :: self
      integer, intent(in) :: from_date, to_date
      integer, intent(out), optional :: status
      real(dp), allocatable :: out(:), daily(:)
      integer :: i, stat
      daily = self%twrr_daily(from_date, to_date, stat)
      allocate(out(size(daily)))
      if (size(daily) > 0) then
         out(1) = daily(1)
         do i = 2, size(daily)
            if (ieee_is_finite(daily(i)) .and. ieee_is_finite(out(i - 1))) then
               out(i) = (1.0_dp + out(i - 1)) * (1.0_dp + daily(i)) - 1.0_dp
            else
               out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
         end do
      end if
      if (present(status)) status = stat
   end function twrr_cumulative

   function cumulative_pnl(self, from_date, to_date, status) result(out)
      class(return_series_type), intent(in) :: self
      integer, intent(in) :: from_date, to_date
      integer, intent(out), optional :: status
      real(dp), allocatable :: out(:)
      integer :: first, last, stat, i, k
      call self%range_indices(from_date, to_date, first, last, stat)
      if (stat /= 0) then
         allocate(out(0))
         if (present(status)) status = stat
         return
      end if
      allocate(out(last - first + 1))
      out(1) = self%pnl(first)
      do i = first + 1, last
         k = i - first + 1
         out(k) = out(k - 1) + self%pnl(i)
      end do
      if (present(status)) status = 0
   end function cumulative_pnl

   function dietz_average_capital(self, from_date, to_date, status) result(out)
      class(return_series_type), intent(in) :: self
      integer, intent(in) :: from_date, to_date
      integer, intent(out), optional :: status
      real(dp), allocatable :: out(:)
      integer :: first, last, stat, i, j, k, total_days
      real(dp) :: weighted, cf, weight
      call self%range_indices(from_date, to_date, first, last, stat)
      if (stat /= 0 .or. first <= 1) then
         allocate(out(0))
         if (present(status)) status = merge(stat, 3, stat /= 0)
         return
      end if
      allocate(out(last - first + 1))
      do i = first, last
         k = i - first + 1
         total_days = i - first + 1
         weighted = 0.0_dp
         do j = first, i
            cf = self%cashflow(j)
            weight = real(i - j + merge(1, 0, cf > 0.0_dp), dp) / real(total_days, dp)
            weighted = weighted + cf * weight
         end do
         out(k) = self%market_values(first - 1) + weighted
      end do
      if (present(status)) status = 0
   end function dietz_average_capital

   function dietz(self, from_date, to_date, status) result(out)
      class(return_series_type), intent(in) :: self
      integer, intent(in) :: from_date, to_date
      integer, intent(out), optional :: status
      real(dp), allocatable :: out(:), cpnl(:), avc(:)
      integer :: i, stat
      cpnl = self%cumulative_pnl(from_date, to_date, stat)
      avc = self%dietz_average_capital(from_date, to_date, stat)
      allocate(out(size(cpnl)))
      do i = 1, size(out)
         if (abs(avc(i)) <= tiny(1.0_dp)) then
            out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            out(i) = cpnl(i) / avc(i)
         end if
      end do
      if (present(status)) status = stat
   end function dietz

end module fcl_returns
