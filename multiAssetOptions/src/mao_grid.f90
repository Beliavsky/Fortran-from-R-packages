! This file is part of multiAssetOptions-fortran.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module mao_grid
   use mao_kinds, only: dp
   use mao_status, only: status_type, clear_status, set_status, &
      mao_invalid_argument, mao_allocation_error
   use mao_types, only: pricing_config, grid_set
   implicit none
   private

   public :: node_spacer, build_grid, linear_index, decode_index
   public :: interpolate_value

contains

   subroutine node_spacer(strike, left_bound, right_bound, nodes, density, &
      k_shift, x, status)
      real(dp), intent(in) :: strike, left_bound, right_bound, density
      integer, intent(in) :: nodes, k_shift
      real(dp), allocatable, intent(out) :: x(:)
      type(status_type), intent(out) :: status
      real(dp) :: c, dxi, xi0, mid, factor
      integer :: i, stat

      call clear_status(status)
      if (.not. (left_bound < strike .and. strike < right_bound) .or. &
          nodes <= 2 .or. density < 0.0_dp .or. &
          k_shift < 0 .or. k_shift > 2) then
         call set_status(status, mao_invalid_argument, &
            'invalid node_spacer arguments')
         return
      end if

      allocate(x(nodes), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate node vector')
         return
      end if

      if (density <= 0.0_dp) then
         do i = 1, nodes
            x(i) = left_bound + real(i-1,dp) * &
               (right_bound-left_bound) / real(nodes-1,dp)
         end do
      else
         c = (strike-left_bound) / density
         dxi = (asinh((right_bound-strike)/c) - &
            asinh(-(strike-left_bound)/c)) / real(nodes-1,dp)
         xi0 = asinh(-(strike-left_bound)/c)
         do i = 1, nodes
            x(i) = strike + c * sinh(xi0 + real(i-1,dp)*dxi)
         end do
      end if

      if (k_shift /= 0) then
         do i = 1, nodes-1
            if (x(i) < strike .and. x(i+1) >= strike) then
               if (k_shift == 1) then
                  mid = 0.5_dp * (x(i)+x(i+1))
                  factor = 1.0_dp + (strike-mid)/mid
               else
                  if (abs(x(i)) > 100.0_dp*tiny(1.0_dp)) then
                     factor = strike/x(i)
                  else
                     factor = strike/x(i+1)
                  end if
               end if
               x = x * factor
               exit
            end if
         end do
      end if
   end subroutine node_spacer

   subroutine build_grid(config, grid, status)
      type(pricing_config), intent(in) :: config
      type(grid_set), intent(out) :: grid
      type(status_type), intent(out) :: status
      type(status_type) :: local_status
      real(dp) :: right_bound
      integer :: i, n, stat

      call clear_status(status)
      n = config%opt%n_asset
      allocate(grid%asset(n), grid%dims(n), grid%strides(n), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate grid metadata')
         return
      end if

      grid%n_nodes = 1
      do i = 1, n
         if (config%fd%k_mult(i) <= 0.0_dp) then
            right_bound = config%opt%strike(i) * exp(sqrt(2.0_dp * &
               config%opt%vol(i)**2 * config%opt%ttm * abs(log(0.001_dp))))
         else
            right_bound = config%fd%k_mult(i) * config%opt%strike(i)
         end if
         if (right_bound <= config%opt%strike(i)) then
            call set_status(status, mao_invalid_argument, &
               'right grid boundary must exceed strike')
            return
         end if
         call node_spacer(config%opt%strike(i), config%fd%left_bound(i), &
            right_bound, config%fd%m(i)+1, config%fd%density(i), &
            config%fd%k_shift(i), grid%asset(i)%x, local_status)
         if (.not. local_status%ok()) then
            status = local_status
            return
         end if
         grid%dims(i) = size(grid%asset(i)%x)
         if (i == 1) then
            grid%strides(i) = 1
         else
            grid%strides(i) = grid%strides(i-1) * grid%dims(i-1)
         end if
         if (grid%n_nodes > huge(grid%n_nodes) / grid%dims(i)) then
            call set_status(status, mao_invalid_argument, &
               'grid contains too many nodes for default integer indexing')
            return
         end if
         grid%n_nodes = grid%n_nodes * grid%dims(i)
      end do
   end subroutine build_grid

   pure integer function linear_index(indices, strides) result(index_value)
      integer, intent(in) :: indices(:), strides(:)
      integer :: i

      index_value = 1
      do i = 1, size(indices)
         index_value = index_value + (indices(i)-1) * strides(i)
      end do
   end function linear_index

   pure subroutine decode_index(index_value, dims, indices)
      integer, intent(in) :: index_value
      integer, intent(in) :: dims(:)
      integer, intent(out) :: indices(:)
      integer :: remaining, i

      remaining = index_value - 1
      do i = 1, size(dims)
         indices(i) = mod(remaining,dims(i)) + 1
         remaining = remaining / dims(i)
      end do
   end subroutine decode_index

   subroutine interpolate_value(grid, values, spot, result_value, status)
      type(grid_set), intent(in) :: grid
      real(dp), intent(in) :: values(:)
      real(dp), intent(in) :: spot(:)
      real(dp), intent(out) :: result_value
      type(status_type), intent(out) :: status
      integer, allocatable :: lower(:), upper(:), indices(:)
      real(dp), allocatable :: weight(:)
      integer :: n, i, j, corner, ncorner, stat, idx
      real(dp) :: corner_weight

      call clear_status(status)
      n = size(grid%asset)
      if (size(values) /= grid%n_nodes .or. size(spot) /= n) then
         call set_status(status, mao_invalid_argument, &
            'invalid interpolation dimensions')
         return
      end if
      if (n >= bit_size(ncorner)-1) then
         call set_status(status, mao_invalid_argument, &
            'too many dimensions for multilinear interpolation')
         return
      end if
      allocate(lower(n), upper(n), indices(n), weight(n), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate interpolation workspace')
         return
      end if

      do i = 1, n
         if (spot(i) <= grid%asset(i)%x(1)) then
            lower(i) = 1
            upper(i) = 1
            weight(i) = 0.0_dp
         else if (spot(i) >= grid%asset(i)%x(grid%dims(i))) then
            lower(i) = grid%dims(i)
            upper(i) = grid%dims(i)
            weight(i) = 0.0_dp
         else
            do j = 1, grid%dims(i)-1
               if (grid%asset(i)%x(j) <= spot(i) .and. &
                   spot(i) <= grid%asset(i)%x(j+1)) then
                  lower(i) = j
                  upper(i) = j+1
                  weight(i) = (spot(i)-grid%asset(i)%x(j)) / &
                     (grid%asset(i)%x(j+1)-grid%asset(i)%x(j))
                  exit
               end if
            end do
         end if
      end do

      result_value = 0.0_dp
      ncorner = 2**n
      do corner = 0, ncorner-1
         corner_weight = 1.0_dp
         do i = 1, n
            if (btest(corner,i-1) .and. upper(i) /= lower(i)) then
               indices(i) = upper(i)
               corner_weight = corner_weight * weight(i)
            else
               indices(i) = lower(i)
               if (upper(i) /= lower(i)) &
                  corner_weight = corner_weight * (1.0_dp-weight(i))
            end if
         end do
         idx = linear_index(indices,grid%strides)
         result_value = result_value + corner_weight * values(idx)
      end do
   end subroutine interpolate_value

end module mao_grid
