! SPDX-License-Identifier: GPL-3.0-only
! Anderson acceleration translated from upstream SCS aa.c (MIT license).
module scs_acceleration
   use scs_kinds, only : dp
   use scs_linalg, only : norm_2
   implicit none
   private
   public :: aa_workspace

   type :: aa_workspace
      logical :: type1 = .true.
      integer :: mem = 0
      integer :: dim = 0
      integer :: iter = 0
      logical :: success = .false.
      real(dp) :: regularization = 1.0e-6_dp
      real(dp) :: safeguard_factor = 1.0_dp
      real(dp) :: max_weight_norm = 1.0e10_dp
      real(dp) :: norm_g = 0.0_dp
      real(dp), allocatable :: x(:), f(:), g(:), g_prev(:)
      real(dp), allocatable :: y(:), s(:), d(:)
      real(dp), allocatable :: Ymat(:,:), Smat(:,:), Dmat(:,:)
   contains
      procedure :: init => aa_init
      procedure :: reset => aa_reset
      procedure :: apply => aa_apply
      procedure :: safeguard => aa_safeguard
   end type aa_workspace
contains
   subroutine aa_init(a,dim,mem,type1)
      class(aa_workspace),intent(inout)::a
      integer,intent(in)::dim,mem
      logical,intent(in)::type1
      a%dim=dim;a%mem=min(abs(mem),dim);a%type1=type1;a%iter=0;a%success=.false.
      a%regularization=merge(1.0e-6_dp,1.0e-10_dp,type1)
      if(a%mem<=0)return
      allocate(a%x(dim),a%f(dim),a%g(dim),a%g_prev(dim),a%y(dim),a%s(dim),a%d(dim))
      allocate(a%Ymat(dim,a%mem),a%Smat(dim,a%mem),a%Dmat(dim,a%mem))
      a%x=0.0_dp;a%f=0.0_dp;a%g=0.0_dp;a%g_prev=0.0_dp
      a%Ymat=0.0_dp;a%Smat=0.0_dp;a%Dmat=0.0_dp
   end subroutine aa_init

   subroutine aa_reset(a)
      class(aa_workspace),intent(inout)::a
      a%iter=0;a%success=.false.
   end subroutine aa_reset

   real(dp) function aa_apply(a,f,xin) result(aa_norm)
      class(aa_workspace),intent(inout)::a
      real(dp),intent(inout)::f(:)
      real(dp),intent(in)::xin(:)
      integer::len,idx
      real(dp),allocatable::M(:,:),rhs(:),gamma(:)
      logical::ok
      aa_norm=0.0_dp;a%success=.false.
      if(a%mem<=0)return
      len=min(a%iter,a%mem)
      if(a%iter==0)then
         a%x=xin;a%f=f;a%g_prev=xin-f;a%iter=1;return
      end if
      idx=mod(a%iter-1,a%mem)+1
      a%g=xin-f;a%s=xin-a%x;a%d=f-a%f;a%y=a%g-a%g_prev
      a%Ymat(:,idx)=a%y;a%Smat(:,idx)=a%s;a%Dmat(:,idx)=a%d
      a%x=xin;a%f=f;a%g_prev=a%g;a%norm_g=norm_2(a%g)
      if(a%iter>=a%mem)then
         allocate(M(len,len),rhs(len),gamma(len))
         if(a%type1)then
            M=matmul(transpose(a%Smat(:,1:len)),a%Ymat(:,1:len))
            rhs=matmul(transpose(a%Smat(:,1:len)),a%g)
         else
            M=matmul(transpose(a%Ymat(:,1:len)),a%Ymat(:,1:len))
            rhs=matmul(transpose(a%Ymat(:,1:len)),a%g)
         end if
         call regularize(M,a%regularization)
         gamma=rhs
         call solve_general(M,gamma,ok)
         aa_norm=norm_2(gamma)
         if(ok .and. aa_norm<a%max_weight_norm)then
            f=f-matmul(a%Dmat(:,1:len),gamma)
            a%success=.true.
         else
            aa_norm=-aa_norm
            call a%reset()
            return
         end if
      end if
      a%iter=a%iter+1
   end function aa_apply

   integer function aa_safeguard(a,f_new,x_new) result(status)
      class(aa_workspace),intent(inout)::a
      real(dp),intent(inout)::f_new(:),x_new(:)
      real(dp)::nd
      status=0
      if(.not.a%success)return
      a%success=.false.;nd=norm_2(x_new-f_new)
      if(nd>a%safeguard_factor*a%norm_g)then
         f_new=a%f;x_new=a%x;call a%reset();status=-1
      end if
   end function aa_safeguard

   subroutine regularize(M,reg)
      real(dp),intent(inout)::M(:,:)
      real(dp),intent(in)::reg
      real(dp)::r
      integer::i
      if(reg<=0.0_dp)return
      r=reg*sqrt(sum(M*M))
      do i=1,size(M,1);M(i,i)=M(i,i)+r;end do
   end subroutine regularize

   subroutine solve_general(A,b,ok)
      real(dp),intent(inout)::A(:,:)
      real(dp),intent(inout)::b(:)
      logical,intent(out)::ok
      integer::n,k,i,p
      real(dp)::mx,f,tmp
      real(dp),allocatable::row(:)
      n=size(b);ok=.true.;allocate(row(n))
      do k=1,n-1
         p=k;mx=abs(A(k,k))
         do i=k+1,n
            if(abs(A(i,k))>mx)then;p=i;mx=abs(A(i,k));end if
         end do
         if(mx<1.0e-18_dp)then;ok=.false.;return;end if
         if(p/=k)then
            row=A(k,:);A(k,:)=A(p,:);A(p,:)=row;tmp=b(k);b(k)=b(p);b(p)=tmp
         end if
         do i=k+1,n
            f=A(i,k)/A(k,k);A(i,k)=0.0_dp;A(i,k+1:n)=A(i,k+1:n)-f*A(k,k+1:n);b(i)=b(i)-f*b(k)
         end do
      end do
      if(abs(A(n,n))<1.0e-18_dp)then;ok=.false.;return;end if
      do i=n,1,-1
         if(i<n)b(i)=b(i)-dot_product(A(i,i+1:n),b(i+1:n))
         b(i)=b(i)/A(i,i)
      end do
   end subroutine solve_general
end module scs_acceleration
