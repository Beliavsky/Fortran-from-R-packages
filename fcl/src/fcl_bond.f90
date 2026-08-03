! MIT License. Copyright (c) 2024 fcl authors.
module fcl_bond
   use fcl_kinds, only : dp
   use fcl_dates, only : date_type, add_months, year_frac, days_between
   use fcl_xirr, only : xirr, xnpv
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: fixed_bond_type, bond_value_type, cashflow_type

   type :: cashflow_type
      type(date_type), allocatable :: dates(:)
      real(dp), allocatable :: values(:)
   end type cashflow_type

   type :: bond_value_type
      real(dp) :: ytm = 0.0_dp
      real(dp) :: macaulay_duration = 0.0_dp
      real(dp) :: modified_duration = 0.0_dp
      integer :: status = 0
   end type bond_value_type

   type :: fixed_bond_type
      type(date_type) :: value_date
      type(date_type) :: maturity_date
      real(dp) :: redemption_value = 100.0_dp
      real(dp) :: coupon_rate = 0.0_dp
      integer :: coupon_frequency = 0
   contains
      procedure :: valid => bond_valid
      procedure :: coupon_cashflow
      procedure :: redemption_cashflow
      procedure :: all_cashflow
      procedure :: value => bond_value
      procedure, private :: next_coupon_date
      procedure, private :: coupon_value
      procedure, private :: accrued
   end type fixed_bond_type

contains

   pure logical function bond_valid(self)
      class(fixed_bond_type), intent(in) :: self
      bond_valid = any(self%coupon_frequency == [0, 1, 2, 4, 6, 12]) .and. &
         self%maturity_date%serial() > self%value_date%serial()
   end function bond_valid

   function next_coupon_date(self, ref_date, adjust, found) result(date)
      class(fixed_bond_type), intent(in) :: self
      type(date_type), intent(in) :: ref_date
      logical, intent(in) :: adjust
      logical, intent(out) :: found
      type(date_type) :: date

      if (ref_date%serial() >= self%maturity_date%serial()) then
         found = .false.
         date = self%maturity_date
         return
      end if
      found = .true.
      if (self%coupon_frequency == 0) then
         date = self%maturity_date
      else
         date = add_months(ref_date, 12 / self%coupon_frequency)
         if (adjust .and. date%serial() > self%maturity_date%serial()) date = self%maturity_date
      end if
   end function next_coupon_date

   pure real(dp) function coupon_value(self) result(value)
      class(fixed_bond_type), intent(in) :: self
      if (self%coupon_frequency == 0) then
         value = self%redemption_value * self%coupon_rate * &
            year_frac(self%maturity_date, self%value_date)
      else
         value = self%redemption_value * self%coupon_rate / real(self%coupon_frequency, dp)
      end if
   end function coupon_value

   function accrued(self, ref_date, eod) result(value)
      class(fixed_bond_type), intent(in) :: self
      type(date_type), intent(in) :: ref_date
      logical, intent(in) :: eod
      real(dp) :: value
      type(date_type) :: previous, following
      logical :: found
      integer :: coupon_days, elapsed_days

      value = 0.0_dp
      if (ref_date%serial() > self%maturity_date%serial() .or. &
          ref_date%serial() <= self%value_date%serial()) return
      if (eod .and. ref_date%serial() == self%maturity_date%serial()) return

      previous = self%value_date
      do
         following = self%next_coupon_date(previous, .false., found)
         if (.not. found) return
         if (ref_date%serial() <= following%serial()) exit
         previous = following
      end do
      if (eod .and. ref_date%serial() == following%serial()) return
      coupon_days = days_between(following, previous)
      elapsed_days = days_between(ref_date, previous)
      value = self%coupon_value() * real(elapsed_days, dp) / real(coupon_days, dp)
   end function accrued

   function make_cashflow(self, kind) result(cf)
      class(fixed_bond_type), intent(in) :: self
      integer, intent(in) :: kind
      type(cashflow_type) :: cf
      type(date_type), allocatable :: temp_dates(:)
      real(dp), allocatable :: temp_values(:)
      type(date_type) :: date
      logical :: found
      integer :: n, capacity
      real(dp) :: coupon, redemption

      capacity = 16
      allocate(temp_dates(capacity), temp_values(capacity))
      n = 0
      date = self%value_date
      do
         date = self%next_coupon_date(date, .true., found)
         if (.not. found) exit
         n = n + 1
         if (n > capacity) call grow_arrays(temp_dates, temp_values, capacity)
         coupon = self%accrued(date, .false.)
         redemption = merge(self%redemption_value, 0.0_dp, &
            date%serial() == self%maturity_date%serial())
         temp_dates(n) = date
         select case (kind)
         case (1)
            temp_values(n) = coupon
         case (2)
            temp_values(n) = redemption
         case default
            temp_values(n) = coupon + redemption
         end select
      end do
      allocate(cf%dates(n), cf%values(n))
      if (n > 0) then
         cf%dates = temp_dates(:n)
         cf%values = temp_values(:n)
      end if
   end function make_cashflow

   subroutine grow_arrays(dates, values, capacity)
      type(date_type), allocatable, intent(inout) :: dates(:)
      real(dp), allocatable, intent(inout) :: values(:)
      integer, intent(inout) :: capacity
      type(date_type), allocatable :: dtmp(:)
      real(dp), allocatable :: vtmp(:)
      integer :: old
      old = capacity
      capacity = 2 * capacity
      allocate(dtmp(capacity), vtmp(capacity))
      dtmp(:old) = dates
      vtmp(:old) = values
      call move_alloc(dtmp, dates)
      call move_alloc(vtmp, values)
   end subroutine grow_arrays

   function coupon_cashflow(self) result(cf)
      class(fixed_bond_type), intent(in) :: self
      type(cashflow_type) :: cf
      cf = make_cashflow(self, 1)
   end function coupon_cashflow

   function redemption_cashflow(self) result(cf)
      class(fixed_bond_type), intent(in) :: self
      type(cashflow_type) :: cf
      cf = make_cashflow(self, 2)
   end function redemption_cashflow

   function all_cashflow(self) result(cf)
      class(fixed_bond_type), intent(in) :: self
      type(cashflow_type) :: cf
      cf = make_cashflow(self, 3)
   end function all_cashflow

   function bond_value(self, ref_date, clean_price) result(out)
      class(fixed_bond_type), intent(in) :: self
      type(date_type), intent(in) :: ref_date
      real(dp), intent(in) :: clean_price
      type(bond_value_type) :: out
      type(cashflow_type) :: full
      type(date_type), allocatable :: dates(:)
      real(dp), allocatable :: values(:)
      real(dp) :: dirty_price, shift, npv_plus, npv_minus
      integer :: i, n, j, stat

      if (.not. self%valid() .or. ref_date%serial() >= self%maturity_date%serial()) then
         call set_invalid(out, 1)
         return
      end if
      full = self%all_cashflow()
      n = count([(full%dates(i)%serial() > ref_date%serial(), i = 1, size(full%dates))])
      if (n == 0) then
         call set_invalid(out, 2)
         return
      end if
      allocate(dates(n + 1), values(n + 1))
      dirty_price = clean_price + self%accrued(ref_date, .true.)
      dates(1) = ref_date
      values(1) = -dirty_price
      j = 1
      do i = 1, size(full%dates)
         if (full%dates(i)%serial() > ref_date%serial()) then
            j = j + 1
            dates(j) = full%dates(i)
            values(j) = full%values(i)
         end if
      end do
      out%ytm = xirr(values, dates, status=stat)
      if (stat /= 0) then
         call set_invalid(out, 3)
         return
      end if
      shift = 1.0e-6_dp
      npv_plus = xnpv(out%ytm + shift, values, dates)
      npv_minus = xnpv(out%ytm - shift, values, dates)
      out%modified_duration = -(npv_plus - npv_minus) / (2.0_dp * shift * dirty_price)
      out%macaulay_duration = 0.0_dp
      do i = 2, size(values)
         out%macaulay_duration = out%macaulay_duration + values(i) * &
            year_frac(dates(i), ref_date) * &
            (1.0_dp + out%ytm)**(-year_frac(dates(i), ref_date))
      end do
      out%macaulay_duration = out%macaulay_duration / dirty_price
      out%status = 0
   end function bond_value

   subroutine set_invalid(out, status)
      type(bond_value_type), intent(out) :: out
      integer, intent(in) :: status
      real(dp) :: nan
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      out%ytm = nan
      out%macaulay_duration = nan
      out%modified_duration = nan
      out%status = status
   end subroutine set_invalid

end module fcl_bond
