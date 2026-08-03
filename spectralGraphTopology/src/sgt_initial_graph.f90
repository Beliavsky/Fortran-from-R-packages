! SPDX-License-Identifier: GPL-3.0-only
module sgt_initial_graph
   use sgt_kinds, only : dp
   use sgt_status, only : sgt_ok, sgt_invalid_input
   use sgt_utils, only : pairwise_matrix_rownorm2
   use sgt_linalg, only : sort_indices_ascending
   implicit none
   private
   public :: build_initial_graph, laplacian_from_directed_affinity
contains
   subroutine build_initial_graph(y,m,affinity,status)
      real(dp), intent(in) :: y(:,:)
      integer, intent(in) :: m
      real(dp), allocatable, intent(out) :: affinity(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: distances(:,:)
      integer, allocatable :: order(:)
      real(dp) :: den,ei
      integer :: n,i,j,t
      n=size(y,1)
      if (n<3 .or. size(y,2)<1 .or. m<1 .or. m>n-2) then
         allocate(affinity(0,0)); if (present(status)) status=sgt_invalid_input; return
      end if
      distances=pairwise_matrix_rownorm2(y)
      allocate(affinity(n,n),order(n)); affinity=0.0_dp
      do i=1,n
         call sort_indices_ascending(distances(i,:),order)
         ei=distances(i,order(m+2))
         den=real(m,dp)*ei
         do t=2,m+1
            den=den-distances(i,order(t))
         end do
         if (den<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,ei)) then
            do t=2,m+1
               affinity(i,order(t))=1.0_dp/real(m,dp)
            end do
         else
            do t=2,m+1
               j=order(t)
               affinity(i,j)=(ei-distances(i,j))/den
            end do
         end if
      end do
      if (present(status)) status=sgt_ok
   end subroutine build_initial_graph

   pure function laplacian_from_directed_affinity(affinity) result(laplacian)
      real(dp), intent(in) :: affinity(:,:)
      real(dp), allocatable :: laplacian(:,:)
      real(dp), allocatable :: sym(:,:)
      integer :: n,i
      n=size(affinity,1)
      if (n<1 .or. size(affinity,2)/=n) then
         allocate(laplacian(0,0)); return
      end if
      allocate(sym(n,n),laplacian(n,n))
      sym=0.5_dp*(affinity+transpose(affinity))
      laplacian=-sym
      do i=1,n
         laplacian(i,i)=sum(sym(:,i))-sym(i,i)
      end do
   end function laplacian_from_directed_affinity
end module sgt_initial_graph
