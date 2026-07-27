! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_types
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use ob_kinds, only : dp, i8
   implicit none
   private

   integer, parameter, public :: action_created = 1
   integer, parameter, public :: action_changed = 2
   integer, parameter, public :: action_deleted = 3

   integer, parameter, public :: side_bid = 1
   integer, parameter, public :: side_ask = 2

   integer, parameter, public :: type_unknown = 0
   integer, parameter, public :: type_flashed_limit = 1
   integer, parameter, public :: type_resting_limit = 2
   integer, parameter, public :: type_market_limit = 3
   integer, parameter, public :: type_pacman = 4
   integer, parameter, public :: type_market = 5

   integer, parameter, public :: trade_sell = -1
   integer, parameter, public :: trade_buy = 1

   type, public :: event_t
      integer :: event_id = 0
      integer(i8) :: id = 0_i8
      integer(i8) :: timestamp_ms = 0_i8
      integer(i8) :: exchange_timestamp_ms = 0_i8
      real(dp) :: price = 0.0_dp
      real(dp) :: volume = 0.0_dp
      integer :: action = action_created
      integer :: side = side_bid
      real(dp) :: fill = 0.0_dp
      integer :: matching_event = 0
      integer :: order_type = type_unknown
      real(dp) :: aggressiveness_bps = 0.0_dp
      logical :: has_aggressiveness = .false.
   end type event_t

   type, public :: trade_t
      integer(i8) :: timestamp_ms = 0_i8
      real(dp) :: price = 0.0_dp
      real(dp) :: volume = 0.0_dp
      integer :: direction = trade_buy
      integer :: maker_event_id = 0
      integer :: taker_event_id = 0
      integer(i8) :: maker_id = 0_i8
      integer(i8) :: taker_id = 0_i8
   end type trade_t

   type, public :: depth_update_t
      integer(i8) :: timestamp_ms = 0_i8
      real(dp) :: price = 0.0_dp
      real(dp) :: volume = 0.0_dp
      integer :: side = side_bid
   end type depth_update_t

   type, public :: depth_summary_t
      integer :: bins = 0
      integer :: bps = 0
      integer(i8), allocatable :: timestamp_ms(:)
      real(dp), allocatable :: best_bid_price(:)
      real(dp), allocatable :: best_bid_volume(:)
      real(dp), allocatable :: bid_volume(:,:)
      real(dp), allocatable :: best_ask_price(:)
      real(dp), allocatable :: best_ask_volume(:)
      real(dp), allocatable :: ask_volume(:,:)
   end type depth_summary_t

   type, public :: spread_t
      integer(i8), allocatable :: timestamp_ms(:)
      real(dp), allocatable :: bid_price(:)
      real(dp), allocatable :: bid_volume(:)
      real(dp), allocatable :: ask_price(:)
      real(dp), allocatable :: ask_volume(:)
   end type spread_t

   type, public :: order_level_t
      integer(i8) :: id = 0_i8
      integer(i8) :: timestamp_ms = 0_i8
      integer(i8) :: exchange_timestamp_ms = 0_i8
      real(dp) :: price = 0.0_dp
      real(dp) :: volume = 0.0_dp
      real(dp) :: liquidity = 0.0_dp
      real(dp) :: bps = 0.0_dp
   end type order_level_t

   type, public :: order_book_t
      integer(i8) :: timestamp_ms = 0_i8
      type(order_level_t), allocatable :: asks(:)
      type(order_level_t), allocatable :: bids(:)
   end type order_book_t

   type, public :: impact_t
      integer(i8) :: id = 0_i8
      real(dp) :: min_price = 0.0_dp
      real(dp) :: max_price = 0.0_dp
      real(dp) :: vwap = 0.0_dp
      integer :: hits = 0
      real(dp) :: volume = 0.0_dp
      integer(i8) :: start_time_ms = 0_i8
      integer(i8) :: end_time_ms = 0_i8
      integer :: direction = trade_buy
   end type impact_t

   type, public :: processing_result_t
      type(event_t), allocatable :: events(:)
      type(trade_t), allocatable :: trades(:)
      type(depth_update_t), allocatable :: depth(:)
      type(depth_summary_t) :: depth_summary
   end type processing_result_t

   public :: action_name, side_name, order_type_name, trade_direction_name
   public :: nan_dp

contains

   pure function nan_dp() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   pure function action_name(action) result(name)
      integer, intent(in) :: action
      character(len=:), allocatable :: name
      select case (action)
      case (action_created); name = 'created'
      case (action_changed); name = 'changed'
      case (action_deleted); name = 'deleted'
      case default; name = 'unknown'
      end select
   end function action_name

   pure function side_name(side) result(name)
      integer, intent(in) :: side
      character(len=:), allocatable :: name
      select case (side)
      case (side_bid); name = 'bid'
      case (side_ask); name = 'ask'
      case default; name = 'unknown'
      end select
   end function side_name

   pure function order_type_name(order_type) result(name)
      integer, intent(in) :: order_type
      character(len=:), allocatable :: name
      select case (order_type)
      case (type_flashed_limit); name = 'flashed-limit'
      case (type_resting_limit); name = 'resting-limit'
      case (type_market_limit); name = 'market-limit'
      case (type_pacman); name = 'pacman'
      case (type_market); name = 'market'
      case default; name = 'unknown'
      end select
   end function order_type_name

   pure function trade_direction_name(direction) result(name)
      integer, intent(in) :: direction
      character(len=:), allocatable :: name
      if (direction == trade_sell) then
         name = 'sell'
      else
         name = 'buy'
      end if
   end function trade_direction_name

end module ob_types
