! SPDX-License-Identifier: GPL-3.0-only
module fingraph_operators
   use fingraph_kinds, only : dp
   implicit none
   private
   public :: L, A, D, Lstar, Astar, Dstar, Linv, Ainv
   public :: Mmat, Pmat, Dmat, vec, vecLmat
contains
   pure integer function node_count_from_edge_count(m) result(n)
      integer, intent(in) :: m
      n=nint(0.5_dp*(1.0_dp+sqrt(1.0_dp+8.0_dp*real(m,dp))))
      if (n*(n-1)/2/=m) n=0
   end function node_count_from_edge_count

   pure function L(w) result(lw)
      real(dp), intent(in) :: w(:)
      real(dp), allocatable :: lw(:,:)
      integer :: n,i,j,k
      n=node_count_from_edge_count(size(w))
      if (n==0) then
         allocate(lw(0,0)); return
      end if
      allocate(lw(n,n)); lw=0.0_dp; k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1
            lw(i,j)=-w(k); lw(j,i)=-w(k)
            lw(i,i)=lw(i,i)+w(k); lw(j,j)=lw(j,j)+w(k)
         end do
      end do
   end function L

   pure function A(w) result(aw)
      real(dp), intent(in) :: w(:)
      real(dp), allocatable :: aw(:,:)
      integer :: n,i,j,k
      n=node_count_from_edge_count(size(w))
      if (n==0) then
         allocate(aw(0,0)); return
      end if
      allocate(aw(n,n)); aw=0.0_dp; k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1; aw(i,j)=w(k); aw(j,i)=w(k)
         end do
      end do
   end function A

   pure function Lstar(m) result(w)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: w(:)
      integer :: n,i,j,k
      n=size(m,1)
      if (n<1 .or. size(m,2)/=n) then
         allocate(w(0)); return
      end if
      allocate(w(n*(n-1)/2)); k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1
            w(k)=m(i,i)+m(j,j)-m(j,i)-m(i,j)
         end do
      end do
   end function Lstar

   pure function Astar(m) result(w)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: w(:)
      integer :: n,i,j,k
      n=size(m,1)
      if (n<1 .or. size(m,2)/=n) then
         allocate(w(0)); return
      end if
      allocate(w(n*(n-1)/2)); k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1; w(k)=m(j,i)+m(i,j)
         end do
      end do
   end function Astar

   pure function Linv(m) result(w)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: w(:)
      integer :: n,i,j,k
      n=size(m,1)
      if (n<1 .or. size(m,2)/=n) then
         allocate(w(0)); return
      end if
      allocate(w(n*(n-1)/2)); k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1; w(k)=-m(i,j)
         end do
      end do
   end function Linv

   pure function Ainv(m) result(w)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: w(:)
      integer :: n,i,j,k
      n=size(m,1)
      if (n<1 .or. size(m,2)/=n) then
         allocate(w(0)); return
      end if
      allocate(w(n*(n-1)/2)); k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1; w(k)=m(i,j)
         end do
      end do
   end function Ainv

   pure function D(w) result(dw)
      real(dp), intent(in) :: w(:)
      real(dp), allocatable :: dw(:)
      real(dp), allocatable :: aw(:,:)
      aw=A(w)
      allocate(dw(size(aw,1)))
      if (size(aw,1)>0) dw=sum(aw,dim=1)
   end function D

   pure function Dstar(v) result(w)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: w(:)
      integer :: n,i,j,k
      n=size(v); allocate(w(n*(n-1)/2)); k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1; w(k)=v(i)+v(j)
         end do
      end do
   end function Dstar

   pure function vec(m) result(v)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: v(:)
      integer :: i,j,k
      allocate(v(size(m)))
      k=0
      do j=1,size(m,2)
         do i=1,size(m,1)
            k=k+1; v(k)=m(i,j)
         end do
      end do
   end function vec

   function Mmat(n) result(m)
      integer, intent(in) :: n
      real(dp), allocatable :: m(:,:)
      real(dp), allocatable :: e(:)
      integer :: j
      if (n<1) then
         allocate(m(0,0)); return
      end if
      allocate(m(n,n),e(n)); m=0.0_dp
      do j=1,n
         e=0.0_dp; e(j)=1.0_dp; m(:,j)=Lstar(L(e))
      end do
   end function Mmat

   function Pmat(n) result(m)
      integer, intent(in) :: n
      real(dp), allocatable :: m(:,:)
      real(dp), allocatable :: e(:)
      integer :: j
      if (n<1) then
         allocate(m(0,0)); return
      end if
      allocate(m(n,n),e(n)); m=0.0_dp
      do j=1,n
         e=0.0_dp; e(j)=1.0_dp; m(:,j)=Astar(A(e))
      end do
   end function Pmat

   function Dmat(n) result(m)
      integer, intent(in) :: n
      real(dp), allocatable :: m(:,:)
      real(dp), allocatable :: e(:)
      integer :: j
      if (n<1) then
         allocate(m(0,0)); return
      end if
      allocate(m(n,n),e(n)); m=0.0_dp
      do j=1,n
         e=0.0_dp; e(j)=1.0_dp; m(:,j)=Dstar(D(e))
      end do
   end function Dmat

   function vecLmat(n) result(r)
      integer, intent(in) :: n
      real(dp), allocatable :: r(:,:)
      real(dp), allocatable :: e(:)
      integer :: m,j
      m=n*(n-1)/2
      if (n<1) then
         allocate(r(0,0)); return
      end if
      allocate(r(n*n,m),e(m)); r=0.0_dp
      do j=1,m
         e=0.0_dp; e(j)=1.0_dp; r(:,j)=vec(L(e))
      end do
   end function vecLmat
end module fingraph_operators
