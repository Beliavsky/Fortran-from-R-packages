! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_knapsack
   use nilde_kinds, only : i8, dp
   use nilde_types, only : knapsack_result_t
   use nilde_collect, only : append_i8_column
   implicit none
   private
   public :: get_knapsack
contains

   function get_knapsack(objective, a, n, problem, bounds, legacy_unbounded_all) result(res)
      real(dp), intent(in) :: objective(:)
      integer(i8), intent(in) :: a(:), n
      character(len=*), intent(in), optional :: problem
      integer(i8), intent(in), optional :: bounds(:)
      logical, intent(in), optional :: legacy_unbounded_all
      type(knapsack_result_t) :: res
      character(len=:), allocatable :: p
      integer(i8), allocatable :: b(:), x(:), store(:,:)
      real(dp), allocatable :: vals(:)
      integer :: i, ns, nbest
      real(dp) :: best, tol
      logical :: legacy

      if (size(a) < 2 .or. size(objective) /= size(a)) error stop "get_knapsack: incompatible inputs"
      if (any(a <= 0_i8) .or. n <= 0_i8) error stop "get_knapsack: positive weights/capacity required"
      p='uknap'; if (present(problem)) p=trim(problem)
      legacy=.true.; if (present(legacy_unbounded_all)) legacy=legacy_unbounded_all
      allocate(b(size(a)), x(size(a))); x=0_i8
      select case(p)
      case('uknap')
         b = n/minval(a)
      case('knap01')
         b = 1_i8
      case('bknap')
         if (.not. present(bounds)) error stop "get_knapsack: bounds required for bknap"
         if (size(bounds) /= size(a)) error stop "get_knapsack: bounds mismatch"
         b = bounds
      case default
         error stop "get_knapsack: unknown problem"
      end select
      ns=0
      call rec(size(a), n)
      if (ns == 0) return
      allocate(vals(ns))
      do i=1,ns
         vals(i)=sum(objective*real(store(:,i),dp))
      end do
      best=maxval(vals); res%objective=best
      tol=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(best))
      if (p == 'uknap' .and. legacy) then
         res%nsol=ns; res%legacy_unbounded_all=.true.
         allocate(res%x(size(a),ns)); res%x=store(:,1:ns)
      else
         nbest=count(abs(vals-best)<=tol)
         res%nsol=nbest; allocate(res%x(size(a),nbest))
         nbest=0
         do i=1,ns
            if (abs(vals(i)-best)<=tol) then
               nbest=nbest+1; res%x(:,nbest)=store(:,i)
            end if
         end do
      end if

   contains
      recursive subroutine rec(idx, rem)
         integer, intent(in) :: idx
         integer(i8), intent(in) :: rem
         integer(i8) :: v, vmax
         if (idx == 0) then
            call append_i8_column(store, ns, x)
            return
         end if
         vmax=min(b(idx),rem/a(idx))
         do v=0_i8,vmax
            x(idx)=v
            call rec(idx-1,rem-v*a(idx))
         end do
         x(idx)=0_i8
      end subroutine rec
   end function get_knapsack

end module nilde_knapsack
