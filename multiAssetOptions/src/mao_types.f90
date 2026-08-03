! This file is part of multiAssetOptions-fortran.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module mao_types
   use mao_kinds, only: dp
   use mao_status, only: status_type, clear_status, set_status, &
      mao_invalid_argument, mao_allocation_error
   implicit none
   private

   integer, parameter, public :: payoff_digital = 0
   integer, parameter, public :: payoff_best_of = 1
   integer, parameter, public :: payoff_worst_of = 2
   integer, parameter, public :: exercise_european = 0
   integer, parameter, public :: exercise_american = 1
   integer, parameter, public :: option_call = 0
   integer, parameter, public :: option_put = 1
   integer, parameter, public :: timestep_constant = 0
   integer, parameter, public :: timestep_adaptive = 1

   type, public :: option_spec
      integer :: n_asset = 0
      integer :: pay_type = payoff_best_of
      integer :: exercise_type = exercise_european
      integer, allocatable :: pc_flag(:)
      real(dp) :: ttm = 1.0_dp
      real(dp), allocatable :: strike(:)
      real(dp) :: rf = 0.0_dp
      real(dp), allocatable :: q(:)
      real(dp), allocatable :: vol(:)
      real(dp), allocatable :: rho(:,:)
   end type option_spec

   type, public :: fd_spec
      integer, allocatable :: m(:)
      real(dp), allocatable :: left_bound(:)
      real(dp), allocatable :: k_mult(:)
      real(dp), allocatable :: density(:)
      integer, allocatable :: k_shift(:)
      real(dp) :: theta = 0.5_dp
      integer :: max_smooth = 2
      real(dp) :: tol = 1.0e-7_dp
      integer :: max_iter = 10
   end type fd_spec

   type, public :: time_spec
      integer :: ts_type = timestep_constant
      integer :: n_steps = 100
      real(dp) :: dt_init = 0.01_dp
      real(dp) :: d_norm = 1.0_dp
      real(dp) :: scale_d = 0.05_dp
   end type time_spec

   type, public :: pricing_config
      type(option_spec) :: opt
      type(fd_spec) :: fd
      type(time_spec) :: time
   end type pricing_config

   type, public :: asset_grid
      real(dp), allocatable :: x(:)
   end type asset_grid

   type, public :: grid_set
      type(asset_grid), allocatable :: asset(:)
      integer, allocatable :: dims(:)
      integer, allocatable :: strides(:)
      integer :: n_nodes = 0
   end type grid_set

   type, public :: pricing_result
      type(grid_set) :: grid
      real(dp), allocatable :: value(:,:)
      real(dp), allocatable :: time(:)
      integer, allocatable :: linear_iterations(:)
      integer, allocatable :: penalty_iterations(:)
   end type pricing_result

   public :: initialize_config, validate_config

contains

   subroutine initialize_config(config, n_asset, status)
      type(pricing_config), intent(out) :: config
      integer, intent(in) :: n_asset
      type(status_type), intent(out) :: status
      integer :: i, stat

      call clear_status(status)
      if (n_asset < 1) then
         call set_status(status, mao_invalid_argument, &
            'n_asset must be at least one')
         return
      end if

      config%opt%n_asset = n_asset
      allocate(config%opt%pc_flag(n_asset), config%opt%strike(n_asset), &
         config%opt%q(n_asset), config%opt%vol(n_asset), &
         config%opt%rho(n_asset,n_asset), config%fd%m(n_asset), &
         config%fd%left_bound(n_asset), config%fd%k_mult(n_asset), &
         config%fd%density(n_asset), config%fd%k_shift(n_asset), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate pricing configuration')
         return
      end if

      config%opt%pc_flag = option_call
      config%opt%strike = 100.0_dp
      config%opt%q = 0.0_dp
      config%opt%vol = 0.2_dp
      config%opt%rho = 0.0_dp
      do i = 1, n_asset
         config%opt%rho(i,i) = 1.0_dp
      end do
      config%fd%m = 40
      config%fd%left_bound = 0.0_dp
      config%fd%k_mult = 4.0_dp
      config%fd%density = 5.0_dp
      config%fd%k_shift = 1
   end subroutine initialize_config

   subroutine validate_config(config, status)
      type(pricing_config), intent(in) :: config
      type(status_type), intent(out) :: status
      integer :: n, i, j
      real(dp), parameter :: corr_tol = 1.0e-10_dp

      call clear_status(status)
      n = config%opt%n_asset
      if (n < 1) then
         call set_status(status, mao_invalid_argument, &
            'n_asset must be at least one')
         return
      end if
      if (.not. allocated(config%opt%pc_flag) .or. &
          .not. allocated(config%opt%strike) .or. &
          .not. allocated(config%opt%q) .or. &
          .not. allocated(config%opt%vol) .or. &
          .not. allocated(config%opt%rho) .or. &
          .not. allocated(config%fd%m) .or. &
          .not. allocated(config%fd%left_bound) .or. &
          .not. allocated(config%fd%k_mult) .or. &
          .not. allocated(config%fd%density) .or. &
          .not. allocated(config%fd%k_shift)) then
         call set_status(status, mao_invalid_argument, &
            'all per-asset configuration arrays must be allocated')
         return
      end if
      if (size(config%opt%pc_flag) /= n .or. &
          size(config%opt%strike) /= n .or. size(config%opt%q) /= n .or. &
          size(config%opt%vol) /= n .or. size(config%opt%rho,1) /= n .or. &
          size(config%opt%rho,2) /= n .or. size(config%fd%m) /= n .or. &
          size(config%fd%left_bound) /= n .or. &
          size(config%fd%k_mult) /= n .or. &
          size(config%fd%density) /= n .or. &
          size(config%fd%k_shift) /= n) then
         call set_status(status, mao_invalid_argument, &
            'all per-asset arrays must have length n_asset')
         return
      end if
      if (config%opt%pay_type < payoff_digital .or. &
          config%opt%pay_type > payoff_worst_of) then
         call set_status(status, mao_invalid_argument, 'invalid payoff type')
         return
      end if
      if (config%opt%exercise_type /= exercise_european .and. &
          config%opt%exercise_type /= exercise_american) then
         call set_status(status, mao_invalid_argument, 'invalid exercise type')
         return
      end if
      if (config%opt%ttm <= 0.0_dp) then
         call set_status(status, mao_invalid_argument, &
            'time to maturity must be positive')
         return
      end if
      if (config%fd%theta < 0.0_dp .or. config%fd%theta > 1.0_dp) then
         call set_status(status, mao_invalid_argument, &
            'theta must lie between zero and one')
         return
      end if
      if (config%fd%tol <= 0.0_dp .or. config%fd%max_iter < 1 .or. &
          config%fd%max_smooth < 0) then
         call set_status(status, mao_invalid_argument, &
            'invalid finite-difference iteration controls')
         return
      end if
      if (config%time%ts_type == timestep_constant) then
         if (config%time%n_steps < 1) then
            call set_status(status, mao_invalid_argument, &
               'constant stepping requires n_steps >= 1')
            return
         end if
      else if (config%time%ts_type == timestep_adaptive) then
         if (config%time%dt_init <= 0.0_dp .or. &
             config%time%d_norm <= 0.0_dp .or. &
             config%time%scale_d <= 0.0_dp) then
            call set_status(status, mao_invalid_argument, &
               'invalid adaptive timestep controls')
            return
         end if
      else
         call set_status(status, mao_invalid_argument, 'invalid timestep type')
         return
      end if

      do i = 1, n
         if (config%opt%pc_flag(i) /= option_call .and. &
             config%opt%pc_flag(i) /= option_put) then
            call set_status(status, mao_invalid_argument, &
               'pc_flag entries must be zero or one')
            return
         end if
         if (config%opt%strike(i) <= 0.0_dp .or. &
             config%opt%vol(i) < 0.0_dp .or. config%fd%m(i) < 2 .or. &
             config%fd%left_bound(i) < 0.0_dp .or. &
             config%fd%left_bound(i) >= config%opt%strike(i) .or. &
             config%fd%k_mult(i) < 0.0_dp .or. &
             config%fd%density(i) < 0.0_dp .or. &
             config%fd%k_shift(i) < 0 .or. config%fd%k_shift(i) > 2) then
            call set_status(status, mao_invalid_argument, &
               'invalid per-asset option or grid input')
            return
         end if
         if (abs(config%opt%rho(i,i) - 1.0_dp) > corr_tol) then
            call set_status(status, mao_invalid_argument, &
               'correlation matrix diagonal must equal one')
            return
         end if
         do j = 1, n
            if (abs(config%opt%rho(i,j)) > 1.0_dp + corr_tol .or. &
                abs(config%opt%rho(i,j) - config%opt%rho(j,i)) > corr_tol) then
               call set_status(status, mao_invalid_argument, &
                  'rho must be a symmetric correlation matrix')
               return
            end if
         end do
      end do
   end subroutine validate_config

end module mao_types
