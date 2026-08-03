! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_types
   use ragtop_kinds, only : dp
   use ragtop_constants, only : instrument_european_option
   implicit none
   private

   type, public :: option_value
      real(dp) :: price = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: vega = 0.0_dp
      integer :: status = 0
   end type option_value

   type, public :: dividend_schedule
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: fixed(:)
      real(dp), allocatable :: proportional(:)
   contains
      procedure :: size => dividend_count
   end type dividend_schedule

   type, public :: cashflow_schedule
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: amount(:)
   contains
      procedure :: size => cashflow_count
   end type cashflow_schedule

   type, public :: exercise_schedule
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: value(:)
   contains
      procedure :: size => exercise_count
   end type exercise_schedule

   type, public :: discount_curve
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: rate(:)
      real(dp), allocatable :: df(:)
      real(dp), allocatable :: forward_rate(:)
   end type discount_curve

   type, public :: volatility_curve
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: volatility(:)
      real(dp), allocatable :: cumulative_variance(:)
      real(dp), allocatable :: forward_volatility(:)
   end type volatility_curve

   type, public :: hazard_spec
      real(dp) :: base_intensity = 0.0_dp
      real(dp) :: constant_fraction = 1.0_dp
      real(dp) :: power = 0.0_dp
      real(dp) :: reference_spot = 1.0_dp
   end type hazard_spec

   type, public :: market_spec
      real(dp) :: short_rate = 0.0_dp
      real(dp) :: volatility = 0.5_dp
      real(dp) :: default_intensity = 0.0_dp
      real(dp) :: dividend_rate = 0.0_dp
      real(dp) :: borrow_cost = 0.0_dp
      type(discount_curve) :: rates
      type(volatility_curve) :: vols
      type(hazard_spec) :: hazard
      logical :: use_rate_curve = .false.
      logical :: use_vol_curve = .false.
      logical :: use_hazard_link = .false.
   end type market_spec

   type, public :: instrument_spec
      integer :: kind = instrument_european_option
      character(len=64) :: name = ''
      real(dp) :: maturity = 0.0_dp
      integer :: callput = 1
      real(dp) :: strike = 0.0_dp
      real(dp) :: notional = 0.0_dp
      real(dp) :: recovery_rate = 0.0_dp
      real(dp) :: conversion_ratio = 0.0_dp
      real(dp) :: dividend_ceiling = huge(1.0_dp)
      logical :: accelerate_future_coupons = .false.
      real(dp) :: acceleration_time = huge(1.0_dp)
      type(cashflow_schedule) :: coupons
      type(exercise_schedule) :: calls
      type(exercise_schedule) :: puts
   end type instrument_spec

   type, public :: grid_spec
      integer :: n_time_steps = 0
      integer :: n_space = 0
      real(dp) :: t_max = 0.0_dp
      real(dp) :: dt = 0.0_dp
      real(dp) :: dz = 0.0_dp
      real(dp) :: z0 = 0.0_dp
      real(dp) :: z_width = 0.0_dp
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: time(:)
   end type grid_spec

   type, public :: price_grid
      real(dp), allocatable :: stock(:)
      real(dp), allocatable :: value(:)
      integer :: status = 0
   end type price_grid

   public :: make_dividend_schedule, make_cashflow_schedule, make_exercise_schedule

   type, public :: greek_result
      real(dp) :: price = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      real(dp) :: vega = 0.0_dp
      real(dp) :: rho = 0.0_dp
      real(dp) :: hazard_sensitivity = 0.0_dp
      real(dp) :: theta = 0.0_dp
      integer :: status = 0
   end type greek_result

contains

   pure integer function dividend_count(this) result(n)
      class(dividend_schedule), intent(in) :: this
      if (allocated(this%time)) then
         n = size(this%time)
      else
         n = 0
      end if
   end function dividend_count

   pure integer function cashflow_count(this) result(n)
      class(cashflow_schedule), intent(in) :: this
      if (allocated(this%time)) then
         n = size(this%time)
      else
         n = 0
      end if
   end function cashflow_count

   pure integer function exercise_count(this) result(n)
      class(exercise_schedule), intent(in) :: this
      if (allocated(this%time)) then
         n = size(this%time)
      else
         n = 0
      end if
   end function exercise_count

   function make_dividend_schedule(time, fixed, proportional) result(schedule)
      real(dp), intent(in) :: time(:), fixed(:), proportional(:)
      type(dividend_schedule) :: schedule
      integer :: n
      n = min(size(time),min(size(fixed),size(proportional)))
      allocate(schedule%time(n),schedule%fixed(n),schedule%proportional(n))
      schedule%time = time(1:n)
      schedule%fixed = fixed(1:n)
      schedule%proportional = proportional(1:n)
   end function make_dividend_schedule

   function make_cashflow_schedule(time, amount) result(schedule)
      real(dp), intent(in) :: time(:), amount(:)
      type(cashflow_schedule) :: schedule
      integer :: n
      n = min(size(time),size(amount))
      allocate(schedule%time(n),schedule%amount(n))
      schedule%time = time(1:n)
      schedule%amount = amount(1:n)
   end function make_cashflow_schedule

   function make_exercise_schedule(time, value) result(schedule)
      real(dp), intent(in) :: time(:), value(:)
      type(exercise_schedule) :: schedule
      integer :: n
      n = min(size(time),size(value))
      allocate(schedule%time(n),schedule%value(n))
      schedule%time = time(1:n)
      schedule%value = value(1:n)
   end function make_exercise_schedule

end module ragtop_types
