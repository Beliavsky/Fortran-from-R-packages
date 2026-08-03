! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_money
   use fmbasics_kinds, only : dp, FM_OK, FM_INVALID_ARGUMENT
   use fmbasics_conventions, only : currency_t
   implicit none
   private

   type, public :: single_currency_money_t
      real(dp), allocatable :: value(:)
      type(currency_t) :: currency
   contains
      procedure :: size => single_money_size
   end type single_currency_money_t

   type, public :: multi_currency_money_t
      real(dp), allocatable :: value(:)
      type(currency_t), allocatable :: currency(:)
   contains
      procedure :: size => multi_money_size
   end type multi_currency_money_t

   type, public :: cash_flow_t
      integer, allocatable :: date(:)
      type(multi_currency_money_t) :: money
   contains
      procedure :: size => cash_flow_size
   end type cash_flow_t

   public :: single_currency_money, multi_currency_money, cash_flow
   public :: combine_money, aggregate_by_currency

   interface single_currency_money
      module procedure single_money_scalar
      module procedure single_money_vector
   end interface single_currency_money

   interface cash_flow
      module procedure cash_flow_scalar_date
      module procedure cash_flow_vector_date
   end interface cash_flow

contains

   function single_money_scalar(value, currency) result(money)
      real(dp), intent(in) :: value
      type(currency_t), intent(in) :: currency
      type(single_currency_money_t) :: money
      allocate(money%value(1))
      money%value = value
      money%currency = currency
   end function single_money_scalar

   function single_money_vector(value, currency) result(money)
      real(dp), intent(in) :: value(:)
      type(currency_t), intent(in) :: currency
      type(single_currency_money_t) :: money
      allocate(money%value(size(value)))
      money%value = value
      money%currency = currency
   end function single_money_vector

   function multi_currency_money(monies, status) result(multi)
      type(single_currency_money_t), intent(in) :: monies(:)
      integer, intent(out), optional :: status
      type(multi_currency_money_t) :: multi
      integer :: i
      allocate(multi%value(size(monies)), multi%currency(size(monies)))
      do i = 1, size(monies)
         if (monies(i)%size() /= 1) then
            allocate(multi%value(0), multi%currency(0))
            if (present(status)) status = FM_INVALID_ARGUMENT
            return
         end if
         multi%value(i) = monies(i)%value(1)
         multi%currency(i) = monies(i)%currency
      end do
      if (present(status)) status = FM_OK
   end function multi_currency_money

   function combine_money(a, b, status) result(multi)
      type(single_currency_money_t), intent(in) :: a, b
      integer, intent(out), optional :: status
      type(multi_currency_money_t) :: multi
      multi = multi_currency_money([a, b], status)
   end function combine_money

   function cash_flow_scalar_date(date, money, status) result(flow)
      integer, intent(in) :: date
      type(multi_currency_money_t), intent(in) :: money
      integer, intent(out), optional :: status
      type(cash_flow_t) :: flow
      allocate(flow%date(money%size()))
      flow%date = date
      flow%money = money
      if (present(status)) status = FM_OK
   end function cash_flow_scalar_date

   function cash_flow_vector_date(date, money, status) result(flow)
      integer, intent(in) :: date(:)
      type(multi_currency_money_t), intent(in) :: money
      integer, intent(out), optional :: status
      type(cash_flow_t) :: flow
      if (size(date) /= money%size()) then
         allocate(flow%date(0), flow%money%value(0), flow%money%currency(0))
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      allocate(flow%date(size(date)))
      flow%date = date
      flow%money = money
      if (present(status)) status = FM_OK
   end function cash_flow_vector_date

   function aggregate_by_currency(money) result(aggregated)
      type(multi_currency_money_t), intent(in) :: money
      type(multi_currency_money_t) :: aggregated
      real(dp), allocatable :: values(:)
      type(currency_t), allocatable :: currencies(:)
      integer :: i, j, n
      allocate(values(money%size()), currencies(money%size()))
      values = 0.0_dp
      n = 0
      do i = 1, money%size()
         j = 0
         if (n > 0) then
            do j = 1, n
               if (currencies(j)%iso == money%currency(i)%iso) exit
            end do
            if (j > n) j = 0
         end if
         if (j == 0) then
            n = n + 1
            currencies(n) = money%currency(i)
            values(n) = money%value(i)
         else
            values(j) = values(j) + money%value(i)
         end if
      end do
      allocate(aggregated%value(n), aggregated%currency(n))
      aggregated%value = values(:n)
      aggregated%currency = currencies(:n)
   end function aggregate_by_currency

   pure integer function single_money_size(self) result(value)
      class(single_currency_money_t), intent(in) :: self
      if (allocated(self%value)) then
         value = size(self%value)
      else
         value = 0
      end if
   end function single_money_size

   pure integer function multi_money_size(self) result(value)
      class(multi_currency_money_t), intent(in) :: self
      if (allocated(self%value)) then
         value = size(self%value)
      else
         value = 0
      end if
   end function multi_money_size

   pure integer function cash_flow_size(self) result(value)
      class(cash_flow_t), intent(in) :: self
      if (allocated(self%date)) then
         value = size(self%date)
      else
         value = 0
      end if
   end function cash_flow_size

end module fmbasics_money
