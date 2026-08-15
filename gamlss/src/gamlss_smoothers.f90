! Penalized-smoother construction used by the GAMLSS fitting engine.
! Algorithms correspond to the computational parts of pb(), ps(), ridge(),
! random(), ri(), and cyclic difference penalties in upstream gamlss.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_smoothers
   use gamlss_kinds, only : dp
   use splines, only : bs_basis, natural_spline_basis
   implicit none
   private
   public :: difference_matrix, difference_penalty, cyclic_difference_penalty
   public :: p_spline_design, natural_spline_design, ridge_penalty
   public :: random_intercept_design, combine_design_penalty, all_pair_difference_matrix
   public :: p_spline_spec_t, fit_p_spline_basis, predict_p_spline_basis

   type, public :: p_spline_spec_t
      integer :: degree = 3
      integer :: order = 2
      real(dp) :: boundary(2) = 0.0_dp
      real(dp), allocatable :: knots(:)
      real(dp), allocatable :: penalty(:,:)
   end type p_spline_spec_t

contains

   subroutine difference_matrix(ncoef, order, d)
      integer, intent(in) :: ncoef, order
      real(dp), allocatable, intent(out) :: d(:,:)
      real(dp), allocatable :: cur(:,:), nxt(:,:)
      integer :: r, i, nr

      if (ncoef <= 0 .or. order < 0 .or. order >= ncoef) then
         allocate(d(0,0))
         return
      end if
      allocate(cur(ncoef,ncoef)); cur = 0.0_dp
      do i=1,ncoef
         cur(i,i)=1.0_dp
      end do
      if (order == 0) then
         d = cur
         return
      end if
      do r=1,order
         nr = size(cur,1)-1
         allocate(nxt(nr,ncoef)); nxt = 0.0_dp
         do i=1,nr
            nxt(i,:) = cur(i+1,:) - cur(i,:)
         end do
         call move_alloc(nxt,cur)
      end do
      call move_alloc(cur,d)
   end subroutine difference_matrix

   subroutine difference_penalty(ncoef, order, penalty)
      integer, intent(in) :: ncoef, order
      real(dp), allocatable, intent(out) :: penalty(:,:)
      real(dp), allocatable :: d(:,:)
      call difference_matrix(ncoef,order,d)
      if (size(d,1)==0) then
         allocate(penalty(max(0,ncoef),max(0,ncoef))); penalty=0.0_dp
      else
         allocate(penalty(ncoef,ncoef))
         penalty = matmul(transpose(d),d)
      end if
   end subroutine difference_penalty

   subroutine cyclic_difference_penalty(ncoef, order, penalty)
      integer, intent(in) :: ncoef, order
      real(dp), allocatable, intent(out) :: penalty(:,:)
      real(dp), allocatable :: d(:,:), next(:,:)
      integer :: i,r,j
      if(ncoef<=1 .or. order<1)then
         allocate(penalty(max(0,ncoef),max(0,ncoef))); penalty=0.0_dp; return
      end if
      allocate(d(ncoef,ncoef)); d=0.0_dp
      do i=1,ncoef
         d(i,i)=-1.0_dp
         j=mod(i,ncoef)+1
         d(i,j)=1.0_dp
      end do
      do r=2,order
         allocate(next(ncoef,ncoef)); next=matmul(d,first_difference(ncoef))
         call move_alloc(next,d)
      end do
      allocate(penalty(ncoef,ncoef)); penalty=matmul(transpose(d),d)
   contains
      function first_difference(n) result(a)
         integer,intent(in)::n
         real(dp)::a(n,n)
         integer::k,l
         a=0.0_dp
         do k=1,n
            a(k,k)=-1.0_dp; l=mod(k,n)+1; a(k,l)=1.0_dp
         end do
      end function first_difference
   end subroutine cyclic_difference_penalty

   subroutine p_spline_design(x, basis, penalty, df, degree, order, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: basis(:,:), penalty(:,:)
      integer, intent(in), optional :: df, degree, order
      integer, intent(out), optional :: status
      integer :: ndf, deg, ord, ierr
      ndf=10; if(present(df))ndf=df
      deg=3; if(present(degree))deg=degree
      ord=2; if(present(order))ord=order
      call bs_basis(x,basis,degree=deg,df=ndf,intercept=.true.,status=ierr)
      if(ierr/=0)then
         if(present(status))status=ierr
         allocate(penalty(0,0)); return
      end if
      call difference_penalty(size(basis,2),ord,penalty)
      if(present(status))status=0
   end subroutine p_spline_design

   subroutine natural_spline_design(x, basis, df, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: basis(:,:)
      integer, intent(in), optional :: df
      integer, intent(out), optional :: status
      integer :: ndf,ierr
      ndf=5; if(present(df))ndf=df
      call natural_spline_basis(x,basis,df=ndf,intercept=.true.,status=ierr)
      if(present(status))status=ierr
   end subroutine natural_spline_design

   subroutine ridge_penalty(ncoef, penalty, penalize_intercept)
      integer,intent(in)::ncoef
      real(dp),allocatable,intent(out)::penalty(:,:)
      logical,intent(in),optional::penalize_intercept
      logical::pi
      integer::i
      pi=.false.; if(present(penalize_intercept))pi=penalize_intercept
      allocate(penalty(ncoef,ncoef)); penalty=0.0_dp
      do i=1,ncoef
         if(i>1 .or. pi) penalty(i,i)=1.0_dp
      end do
   end subroutine ridge_penalty

   subroutine random_intercept_design(group,z,penalty,nlevels)
      integer,intent(in)::group(:)
      real(dp),allocatable,intent(out)::z(:,:),penalty(:,:)
      integer,intent(out),optional::nlevels
      integer::gmax,i
      if(size(group)==0 .or. minval(group)<1)then
         allocate(z(0,0),penalty(0,0)); if(present(nlevels))nlevels=0; return
      end if
      gmax=maxval(group)
      allocate(z(size(group),gmax)); z=0.0_dp
      do i=1,size(group); z(i,group(i))=1.0_dp; end do
      allocate(penalty(gmax,gmax)); penalty=0.0_dp
      do i=1,gmax; penalty(i,i)=1.0_dp; end do
      if(present(nlevels))nlevels=gmax
   end subroutine random_intercept_design

   subroutine combine_design_penalty(x1,p1,x2,p2,x,p)
      real(dp),intent(in)::x1(:,:),p1(:,:),x2(:,:),p2(:,:)
      real(dp),allocatable,intent(out)::x(:,:),p(:,:)
      integer::n,a,b
      n=size(x1,1); a=size(x1,2); b=size(x2,2)
      if(size(x2,1)/=n)then; allocate(x(0,0),p(0,0)); return; end if
      allocate(x(n,a+b)); x(:,1:a)=x1; x(:,a+1:a+b)=x2
      allocate(p(a+b,a+b)); p=0.0_dp
      if(size(p1,1)==a .and. size(p1,2)==a)p(1:a,1:a)=p1
      if(size(p2,1)==b .and. size(p2,2)==b)p(a+1:a+b,a+1:a+b)=p2
   end subroutine combine_design_penalty

   subroutine all_pair_difference_matrix(ncat,d)
      integer,intent(in)::ncat
      real(dp),allocatable,intent(out)::d(:,:)
      integer::i,j,k,nr
      if(ncat<2)then;allocate(d(0,max(0,ncat)));return;end if
      nr=ncat*(ncat-1)/2;allocate(d(nr,ncat));d=0.0_dp;k=0
      do i=2,ncat
         do j=1,i-1
            k=k+1;d(k,i)=1.0_dp;d(k,j)=-1.0_dp
         end do
      end do
   end subroutine all_pair_difference_matrix

   subroutine fit_p_spline_basis(x,spec,basis,df,degree,order,status)
      real(dp),intent(in)::x(:)
      type(p_spline_spec_t),intent(out)::spec
      real(dp),allocatable,intent(out)::basis(:,:)
      integer,intent(in),optional::df,degree,order
      integer,intent(out),optional::status
      integer::ndf,deg,ord,ierr
      real(dp),allocatable::ik(:)
      ndf=10;if(present(df))ndf=df
      deg=3;if(present(degree))deg=degree
      ord=2;if(present(order))ord=order
      spec%degree=deg;spec%order=ord;spec%boundary=[minval(x),maxval(x)]
      call bs_basis(x,basis,degree=deg,df=ndf,intercept=.true., &
         boundary_knots=spec%boundary,interior_knots=ik,status=ierr)
      if(ierr==0)then
         spec%knots=ik
         call difference_penalty(size(basis,2),ord,spec%penalty)
      else
         allocate(spec%knots(0),spec%penalty(0,0))
      end if
      if(present(status))status=ierr
   end subroutine fit_p_spline_basis

   subroutine predict_p_spline_basis(x,spec,basis,status)
      real(dp),intent(in)::x(:)
      type(p_spline_spec_t),intent(in)::spec
      real(dp),allocatable,intent(out)::basis(:,:)
      integer,intent(out),optional::status
      integer::ierr
      call bs_basis(x,basis,degree=spec%degree,knots=spec%knots,intercept=.true., &
         boundary_knots=spec%boundary,status=ierr)
      if(present(status))status=ierr
   end subroutine predict_p_spline_basis

end module gamlss_smoothers
