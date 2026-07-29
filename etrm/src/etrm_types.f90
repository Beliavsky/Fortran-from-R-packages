! SPDX-License-Identifier: MIT
! Derived from etrm 1.0.2, Copyright (c) 2021 etrm authors.
module etrm_types
   use etrm_kinds, only : dp
   implicit none
   private

   type, public :: strategy_result
      character(len=4) :: name = "    "
      real(dp) :: volume = 0.0_dp
      real(dp) :: transaction_cost = 0.0_dp
      logical :: integer_trades = .true.
      real(dp) :: strike_price = 0.0_dp
      real(dp) :: annual_volatility = 0.0_dp
      real(dp) :: interest_rate = 0.0_dp
      integer :: trading_days = 0
      real(dp), allocatable :: market(:)
      real(dp), allocatable :: trade(:)
      real(dp), allocatable :: exposed(:)
      real(dp), allocatable :: position(:)
      real(dp), allocatable :: hedge(:)
      real(dp), allocatable :: target(:)
      real(dp), allocatable :: portfolio(:)
      real(dp), allocatable :: risk_factor(:)
   contains
      procedure :: churn_rate
   end type strategy_result

   type, public :: strategy_summary
      real(dp) :: volume = 0.0_dp
      real(dp) :: churn = 0.0_dp
      real(dp) :: first(7) = 0.0_dp
      real(dp) :: maximum(7) = 0.0_dp
      real(dp) :: minimum(7) = 0.0_dp
      real(dp) :: last(7) = 0.0_dp
   end type strategy_summary

   type, public :: msfc_result
      integer :: n_days = 0
      integer :: n_contracts = 0
      integer :: n_polynomials = 0
      integer, allocatable :: day(:)
      integer, allocatable :: original_index(:)
      integer, allocatable :: start_day(:)
      integer, allocatable :: end_day(:)
      character(len=:), allocatable :: contract(:)
      real(dp), allocatable :: market_price(:)
      real(dp), allocatable :: computed_price(:)
      real(dp), allocatable :: prior(:)
      real(dp), allocatable :: curve(:)
      real(dp), allocatable :: knots(:)
      real(dp), allocatable :: coefficients(:, :)
   end type msfc_result

   public :: summarize_strategy

contains

   pure real(dp) function churn_rate(self) result(churn)
      class(strategy_result), intent(in) :: self
      if (.not. allocated(self%trade) .or. abs(self%volume) <= tiny(1.0_dp)) then
         churn = 0.0_dp
      else
         churn = sum(abs(self%trade)) / abs(self%volume)
      end if
   end function churn_rate

   subroutine summarize_strategy(result, summary)
      type(strategy_result), intent(in) :: result
      type(strategy_summary), intent(out) :: summary
      real(dp), allocatable :: values(:, :)
      integer :: n

      summary%volume = result%volume
      summary%churn = result%churn_rate()
      if (.not. allocated(result%market)) return
      n = size(result%market)
      allocate(values(7, n))
      values(1, :) = result%market
      values(2, :) = result%trade
      values(3, :) = result%exposed
      values(4, :) = result%position
      values(5, :) = result%hedge
      values(6, :) = result%target
      values(7, :) = result%portfolio
      summary%first = values(:, 1)
      summary%maximum = maxval(values, dim=2)
      summary%minimum = minval(values, dim=2)
      summary%last = values(:, n)
   end subroutine summarize_strategy

end module etrm_types
