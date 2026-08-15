! Additive and tensor-product P-spline design front ends.
! These are matrix-first counterparts of common additive smoother composition.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_additive_v03
   use gamlss_kinds, only : dp
   use gamlss_smoothers, only : p_spline_spec_t, fit_p_spline_basis, predict_p_spline_basis
   implicit none
   private
   public :: additive_pspline_spec_t, tensor_pspline_spec_t
   public :: build_additive_p_splines, predict_additive_p_splines
   public :: tensor_p_spline_2d, predict_tensor_p_spline_2d

   type, public :: additive_pspline_spec_t
      type(p_spline_spec_t), allocatable :: term(:)
      integer, allocatable :: ncoef(:)
      real(dp), allocatable :: lambda(:)
      logical :: intercept = .true.
   end type additive_pspline_spec_t

   type, public :: tensor_pspline_spec_t
      type(p_spline_spec_t) :: x_spec,y_spec
      integer :: nx = 0, ny = 0
      real(dp) :: lambda_x = 1.0_dp, lambda_y = 1.0_dp
   end type tensor_pspline_spec_t
contains

   subroutine build_additive_p_splines(x,basis,penalty,spec,df,degree,order,lambda,intercept,status)
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable,intent(out)::basis(:,:),penalty(:,:)
      type(additive_pspline_spec_t),intent(out)::spec
      integer,intent(in),optional::df(:),degree(:),order(:)
      real(dp),intent(in),optional::lambda(:)
      logical,intent(in),optional::intercept
      integer,intent(out),optional::status
      real(dp),allocatable::b(:,:),bc(:,:),t(:,:),pc(:,:)
      integer::m,j,k,total,pos,istat,dj,gj,oj
      logical::inc

      m=size(x,2);inc=.true.;if(present(intercept))inc=intercept
      if(size(x,1)<1.or.m<1)then
         allocate(basis(0,0),penalty(0,0));if(present(status))status=1;return
      end if
      if(present(df))then;if(size(df)/=m)then;if(present(status))status=2;return;end if;end if
      if(present(degree))then;if(size(degree)/=m)then;if(present(status))status=3;return;end if;end if
      if(present(order))then;if(size(order)/=m)then;if(present(status))status=4;return;end if;end if
      if(present(lambda))then;if(size(lambda)/=m)then;if(present(status))status=5;return;end if;end if
      allocate(spec%term(m),spec%ncoef(m),spec%lambda(m));spec%intercept=inc
      spec%lambda=1.0_dp;if(present(lambda))spec%lambda=lambda
      total=merge(1,0,inc)
      ! First pass obtains persistent knot specs and coefficient counts.
      do j=1,m
         dj=10;if(present(df))dj=df(j);gj=3;if(present(degree))gj=degree(j);oj=2;if(present(order))oj=order(j)
         call fit_p_spline_basis(x(:,j),spec%term(j),b,df=dj,degree=gj,order=oj,status=istat)
         if(istat/=0.or.size(b,2)<2)then
            allocate(basis(0,0),penalty(0,0));if(present(status))status=10+j;return
         end if
         spec%ncoef(j)=size(b,2)-1;total=total+spec%ncoef(j)
      end do
      allocate(basis(size(x,1),total),penalty(total,total));basis=0.0_dp;penalty=0.0_dp
      pos=0
      if(inc)then;pos=1;basis(:,1)=1.0_dp;end if
      do j=1,m
         call predict_p_spline_basis(x(:,j),spec%term(j),b,istat)
         k=size(b,2);call sum_zero_transform(k,t)
         bc=matmul(b,t);pc=matmul(transpose(t),matmul(spec%term(j)%penalty,t))
         basis(:,pos+1:pos+k-1)=bc
         penalty(pos+1:pos+k-1,pos+1:pos+k-1)=spec%lambda(j)*pc
         pos=pos+k-1
      end do
      if(present(status))status=0
   end subroutine build_additive_p_splines

   subroutine predict_additive_p_splines(xnew,spec,basis,status)
      real(dp),intent(in)::xnew(:,:)
      type(additive_pspline_spec_t),intent(in)::spec
      real(dp),allocatable,intent(out)::basis(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::b(:,:),t(:,:)
      integer::m,j,total,pos,k,istat
      m=size(spec%term)
      if(size(xnew,2)/=m)then;allocate(basis(0,0));if(present(status))status=1;return;end if
      total=merge(1,0,spec%intercept)+sum(spec%ncoef)
      allocate(basis(size(xnew,1),total));basis=0.0_dp;pos=0
      if(spec%intercept)then;pos=1;basis(:,1)=1.0_dp;end if
      do j=1,m
         call predict_p_spline_basis(xnew(:,j),spec%term(j),b,istat)
         if(istat/=0)then;deallocate(basis);allocate(basis(0,0));if(present(status))status=10+j;return;end if
         k=size(b,2);call sum_zero_transform(k,t)
         basis(:,pos+1:pos+k-1)=matmul(b,t);pos=pos+k-1
      end do
      if(present(status))status=0
   end subroutine predict_additive_p_splines

   subroutine tensor_p_spline_2d(x,y,basis,penalty,spec,df_x,df_y,degree_x,degree_y, &
      order_x,order_y,lambda_x,lambda_y,status)
      real(dp),intent(in)::x(:),y(:)
      real(dp),allocatable,intent(out)::basis(:,:),penalty(:,:)
      type(tensor_pspline_spec_t),intent(out)::spec
      integer,intent(in),optional::df_x,df_y,degree_x,degree_y,order_x,order_y
      real(dp),intent(in),optional::lambda_x,lambda_y
      integer,intent(out),optional::status
      real(dp),allocatable::bx(:,:),by(:,:),ix(:,:),iy(:,:),p1(:,:),p2(:,:)
      integer::i,j,k,idx,nx,ny,istat,dx,dy,gx,gy,ox,oy
      if(size(y)/=size(x).or.size(x)<1)then
         allocate(basis(0,0),penalty(0,0));if(present(status))status=1;return
      end if
      dx=8;if(present(df_x))dx=df_x;dy=8;if(present(df_y))dy=df_y
      gx=3;if(present(degree_x))gx=degree_x;gy=3;if(present(degree_y))gy=degree_y
      ox=2;if(present(order_x))ox=order_x;oy=2;if(present(order_y))oy=order_y
      spec%lambda_x=1.0_dp;if(present(lambda_x))spec%lambda_x=max(0.0_dp,lambda_x)
      spec%lambda_y=1.0_dp;if(present(lambda_y))spec%lambda_y=max(0.0_dp,lambda_y)
      call fit_p_spline_basis(x,spec%x_spec,bx,df=dx,degree=gx,order=ox,status=istat)
      if(istat/=0)then;allocate(basis(0,0),penalty(0,0));if(present(status))status=2;return;end if
      call fit_p_spline_basis(y,spec%y_spec,by,df=dy,degree=gy,order=oy,status=istat)
      if(istat/=0)then;allocate(basis(0,0),penalty(0,0));if(present(status))status=3;return;end if
      nx=size(bx,2);ny=size(by,2);spec%nx=nx;spec%ny=ny
      allocate(basis(size(x),nx*ny));basis=0.0_dp
      do i=1,size(x)
         idx=0
         do j=1,nx;do k=1,ny;idx=idx+1;basis(i,idx)=bx(i,j)*by(i,k);end do;end do
      end do
      call identity_matrix(nx,ix);call identity_matrix(ny,iy)
      call kronecker(spec%x_spec%penalty,iy,p1);call kronecker(ix,spec%y_spec%penalty,p2)
      penalty=spec%lambda_x*p1+spec%lambda_y*p2
      if(present(status))status=0
   end subroutine tensor_p_spline_2d

   subroutine predict_tensor_p_spline_2d(x,y,spec,basis,status)
      real(dp),intent(in)::x(:),y(:)
      type(tensor_pspline_spec_t),intent(in)::spec
      real(dp),allocatable,intent(out)::basis(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::bx(:,:),by(:,:)
      integer::i,j,k,idx,istat
      if(size(y)/=size(x))then;allocate(basis(0,0));if(present(status))status=1;return;end if
      call predict_p_spline_basis(x,spec%x_spec,bx,istat)
      if(istat/=0)then;allocate(basis(0,0));if(present(status))status=2;return;end if
      call predict_p_spline_basis(y,spec%y_spec,by,istat)
      if(istat/=0)then;allocate(basis(0,0));if(present(status))status=3;return;end if
      allocate(basis(size(x),size(bx,2)*size(by,2)));basis=0.0_dp
      do i=1,size(x)
         idx=0
         do j=1,size(bx,2);do k=1,size(by,2);idx=idx+1;basis(i,idx)=bx(i,j)*by(i,k);end do;end do
      end do
      if(present(status))status=0
   end subroutine predict_tensor_p_spline_2d

   subroutine sum_zero_transform(k,t)
      integer,intent(in)::k
      real(dp),allocatable,intent(out)::t(:,:)
      integer::j
      allocate(t(k,k-1));t=0.0_dp
      do j=1,k-1;t(j,j)=1.0_dp;t(k,j)=-1.0_dp;end do
   end subroutine sum_zero_transform

   subroutine identity_matrix(n,a)
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::a(:,:)
      integer::i
      allocate(a(n,n));a=0.0_dp;do i=1,n;a(i,i)=1.0_dp;end do
   end subroutine identity_matrix

   subroutine kronecker(a,b,c)
      real(dp),intent(in)::a(:,:),b(:,:)
      real(dp),allocatable,intent(out)::c(:,:)
      integer::i,j,r1,r2,c1,c2,rs,re,cs,ce
      r1=size(a,1);c1=size(a,2);r2=size(b,1);c2=size(b,2)
      allocate(c(r1*r2,c1*c2));c=0.0_dp
      do i=1,r1;do j=1,c1
         rs=(i-1)*r2+1;re=i*r2;cs=(j-1)*c2+1;ce=j*c2
         c(rs:re,cs:ce)=a(i,j)*b
      end do;end do
   end subroutine kronecker
end module gamlss_additive_v03
