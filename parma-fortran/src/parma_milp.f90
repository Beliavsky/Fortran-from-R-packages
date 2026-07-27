! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Exact binary optimization for small models, replacing common parma MILP uses.
module parma_milp
   use parma_kinds, only: dp
   use parma_types, only: milp_result
   implicit none
   private
   public :: milp_binary_solve

contains

   subroutine milp_binary_solve(c,result,a,lower,upper,maximize,max_variables)
      real(dp), intent(in) :: c(:)
      type(milp_result), intent(out) :: result
      real(dp), intent(in), optional :: a(:,:),lower(:),upper(:)
      logical, intent(in), optional :: maximize
      integer, intent(in), optional :: max_variables
      integer :: n,limit,i
      integer(kind=8) :: states,state
      integer, allocatable :: x(:),bestx(:)
      real(dp) :: value,best
      logical :: do_max,feasible
      real(dp), allocatable :: ax(:)

      n = size(c)
      limit = 24
      if (present(max_variables)) limit = max_variables
      allocate(result%x(n),x(n),bestx(n))
      result%x = 0
      if (n > limit .or. n >= 62) then
         result%status = 2
         result%message = 'too many binary variables for exact enumeration'
         return
      end if
      do_max = .false.
      if (present(maximize)) do_max = maximize
      if (do_max) then
         best = -huge(1.0_dp)
      else
         best = huge(1.0_dp)
      end if
      states = shiftl(1_8,n)
      do state = 0_8,states-1_8
         do i = 1,n
            if (btest(state,i-1)) then
               x(i) = 1
            else
               x(i) = 0
            end if
         end do
         feasible = .true.
         if (present(a)) then
            ax = matmul(a,real(x,dp))
            if (present(lower)) feasible = feasible .and. all(ax >= lower-1.0e-12_dp)
            if (present(upper)) feasible = feasible .and. all(ax <= upper+1.0e-12_dp)
         end if
         if (.not. feasible) cycle
         value = dot_product(c,real(x,dp))
         result%evaluations = result%evaluations+1
         if ((do_max .and. value > best) .or. (.not. do_max .and. value < best)) then
            best = value
            bestx = x
         end if
      end do
      if ((do_max .and. best <= -huge(1.0_dp)/2.0_dp) .or. &
          (.not. do_max .and. best >= huge(1.0_dp)/2.0_dp)) then
         result%status = 1
         result%message = 'no feasible binary point'
      else
         result%status = 0
         result%message = 'optimal by exact enumeration'
         result%objective = best
         result%x = bestx
      end if
   end subroutine milp_binary_solve

end module parma_milp
