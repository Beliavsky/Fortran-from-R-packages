! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_binpacking
   use nilde_kinds, only : i8
   use nilde_types, only : bin_packing_result_t
   use nilde_collect, only : append_int_column
   implicit none
   private
   public :: bin_packing
contains

   function bin_packing(weights, capacity) result(res)
      integer(i8), intent(in) :: weights(:), capacity
      type(bin_packing_result_t) :: res
      integer, allocatable :: assign(:), store(:,:)
      integer(i8), allocatable :: loads(:), suffix(:)
      integer :: n, best, ns, j, b

      n=size(weights)
      if (n == 0 .or. capacity <= 0_i8) error stop "bin_packing: invalid input"
      if (any(weights <= 0_i8) .or. any(weights > capacity)) return
      allocate(assign(n), loads(n), suffix(n+1)); assign=0; loads=0_i8
      suffix(n+1)=0_i8
      do j=n,1,-1; suffix(j)=suffix(j+1)+weights(j); end do
      best=n+1; ns=0
      call rec(1,0)
      if (ns == 0) return
      res%min_bins=best; res%nsol=ns
      allocate(res%assignment(n,ns)); res%assignment=store(:,1:ns)
      allocate(res%bin_ineff(best,ns), res%total_ineff(ns))
      do j=1,ns
         do b=1,best
            res%bin_ineff(b,j)=capacity-sum(weights, mask=res%assignment(:,j)==b)
         end do
         res%total_ineff(j)=sum(res%bin_ineff(:,j))
      end do

   contains
      recursive subroutine rec(item, nbins)
         integer, intent(in) :: item, nbins
         integer :: k
         integer(i8) :: free_now, remaining, extra
         if (nbins > best) return
         if (item > n) then
            if (nbins < best) then
               best=nbins; ns=0
               if (allocated(store)) deallocate(store)
            end if
            if (nbins == best) call append_int_column(store,ns,assign)
            return
         end if

         free_now=0_i8
         do k=1,nbins; free_now=free_now+(capacity-loads(k)); end do
         remaining=suffix(item)
         if (remaining > free_now) then
            extra=(remaining-free_now+capacity-1_i8)/capacity
         else
            extra=0_i8
         end if
         if (nbins+int(extra) > best) return

         do k=1,nbins
            if (loads(k)+weights(item) <= capacity) then
               assign(item)=k; loads(k)=loads(k)+weights(item)
               call rec(item+1,nbins)
               loads(k)=loads(k)-weights(item); assign(item)=0
            end if
         end do
         if (nbins+1 <= best) then
            assign(item)=nbins+1; loads(nbins+1)=weights(item)
            call rec(item+1,nbins+1)
            loads(nbins+1)=0_i8; assign(item)=0
         end if
      end subroutine rec
   end function bin_packing

end module nilde_binpacking
