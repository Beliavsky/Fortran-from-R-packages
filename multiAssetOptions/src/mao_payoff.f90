! This file is part of multiAssetOptions-fortran.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module mao_payoff
   use mao_kinds, only: dp
   use mao_status, only: status_type, clear_status, set_status, &
      mao_invalid_argument, mao_allocation_error
   use mao_types, only: grid_set, payoff_digital, payoff_best_of, &
      payoff_worst_of, option_call, option_put
   use mao_grid, only: decode_index
   implicit none
   private

   public :: payoff_values

contains

   subroutine payoff_values(pay_type, pc_flag, strike, grid, values, status)
      integer, intent(in) :: pay_type
      integer, intent(in) :: pc_flag(:)
      real(dp), intent(in) :: strike(:)
      type(grid_set), intent(in) :: grid
      real(dp), allocatable, intent(out) :: values(:)
      type(status_type), intent(out) :: status
      integer, allocatable :: indices(:)
      real(dp), allocatable :: moneyness(:)
      integer :: row, i, n, stat
      logical :: digital_in_money

      call clear_status(status)
      n = size(grid%asset)
      if (size(pc_flag) /= n .or. size(strike) /= n .or. &
          pay_type < payoff_digital .or. pay_type > payoff_worst_of) then
         call set_status(status, mao_invalid_argument, &
            'invalid payoff dimensions or type')
         return
      end if
      allocate(values(grid%n_nodes), indices(n), moneyness(n), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate payoff workspace')
         return
      end if

      do row = 1, grid%n_nodes
         call decode_index(row,grid%dims,indices)
         digital_in_money = .true.
         do i = 1, n
            select case (pc_flag(i))
            case (option_call)
               moneyness(i) = max(grid%asset(i)%x(indices(i))-strike(i),0.0_dp)
               digital_in_money = digital_in_money .and. &
                  grid%asset(i)%x(indices(i)) >= strike(i)
            case (option_put)
               moneyness(i) = max(strike(i)-grid%asset(i)%x(indices(i)),0.0_dp)
               digital_in_money = digital_in_money .and. &
                  grid%asset(i)%x(indices(i)) <= strike(i)
            case default
               call set_status(status, mao_invalid_argument, &
                  'pc_flag entries must be zero or one')
               return
            end select
         end do

         select case (pay_type)
         case (payoff_digital)
            if (digital_in_money) then
               values(row) = 1.0_dp
            else
               values(row) = 0.0_dp
            end if
         case (payoff_best_of)
            values(row) = maxval(moneyness)
         case (payoff_worst_of)
            values(row) = minval(moneyness)
         end select
      end do
   end subroutine payoff_values

end module mao_payoff
