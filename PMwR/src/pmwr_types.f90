module pmwr_types
   use pmwr_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: state_unknown = 0
   integer, parameter, public :: state_up = 1
   integer, parameter, public :: state_down = -1

   type, public :: journal_type
      integer :: n = 0
      real(dp), allocatable :: timestamp(:)
      real(dp), allocatable :: amount(:)
      real(dp), allocatable :: price(:)
      integer, allocatable :: instrument(:)
      integer, allocatable :: account(:)
   end type journal_type

   type, public :: rebalance_result
      real(dp), allocatable :: current(:)
      real(dp), allocatable :: target(:)
      real(dp), allocatable :: difference(:)
      real(dp), allocatable :: current_value(:)
      real(dp), allocatable :: target_value(:)
      real(dp) :: notional = 0.0_dp
      real(dp) :: turnover = 0.0_dp
      real(dp) :: target_net_value = 0.0_dp
   end type rebalance_result

   type, public :: pl_path_result
      real(dp), allocatable :: cumulative_position(:)
      real(dp), allocatable :: average_price(:)
      real(dp), allocatable :: realized(:)
      real(dp), allocatable :: unrealized(:)
      real(dp), allocatable :: total(:)
      real(dp), allocatable :: volume(:)
   end type pl_path_result

   type, public :: pl_summary_result
      real(dp) :: total_pl = 0.0_dp
      real(dp) :: volume = 0.0_dp
      real(dp) :: average_buy = 0.0_dp
      real(dp) :: average_sell = 0.0_dp
      logical :: closed = .false.
   end type pl_summary_result

   type, public :: drawdown_table
      integer, allocatable :: peak(:)
      integer, allocatable :: trough(:)
      integer, allocatable :: recover(:)
      real(dp), allocatable :: depth(:)
   end type drawdown_table

   type, public :: streak_table
      integer, allocatable :: first(:)
      integer, allocatable :: last(:)
      integer, allocatable :: state(:)
      real(dp), allocatable :: change(:)
   end type streak_table

   type, public :: unit_price_result
      real(dp), allocatable :: price(:)
      real(dp), allocatable :: units(:)
      real(dp), allocatable :: issued_units(:)
   end type unit_price_result

   type, public :: nav_summary_result
      real(dp) :: initial_value = 0.0_dp
      real(dp) :: final_value = 0.0_dp
      real(dp) :: total_return = 0.0_dp
      real(dp) :: annualized_return = 0.0_dp
      real(dp) :: annualized_volatility = 0.0_dp
      real(dp) :: sharpe = 0.0_dp
      real(dp) :: max_drawdown = 0.0_dp
      integer :: max_drawdown_peak = 0
      integer :: max_drawdown_trough = 0
      integer :: max_drawdown_recover = 0
   end type nav_summary_result

   type, public :: backtest_result
      real(dp), allocatable :: position(:,:)
      real(dp), allocatable :: suggested_position(:,:)
      real(dp), allocatable :: cash(:)
      real(dp), allocatable :: wealth(:)
      real(dp), allocatable :: cumulative_cost(:)
      type(journal_type) :: journal
      real(dp) :: initial_wealth = 0.0_dp
   end type backtest_result

   type, public :: attribution_result
      real(dp), allocatable :: allocation(:,:)
      real(dp), allocatable :: selection(:,:)
      real(dp), allocatable :: interaction(:,:)
      real(dp), allocatable :: linked_total(:)
   end type attribution_result

end module pmwr_types
