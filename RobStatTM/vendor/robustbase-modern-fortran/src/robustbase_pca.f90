! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_pca
   use robustbase_kinds, only: dp
   implicit none
   private
   public :: pca_result, rank_mm, classical_pca

   interface
      subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
         import dp
         character(len=1),intent(in)::jobu,jobvt
         integer,intent(in)::m,n,lda,ldu,ldvt,lwork
         real(dp),intent(inout)::a(lda,*),work(*)
         real(dp),intent(out)::s(*),u(ldu,*),vt(ldvt,*)
         integer,intent(out)::info
      end subroutine dgesvd
   end interface

   type :: pca_result
      integer :: rank = 0
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: loadings(:,:)
      real(dp), allocatable :: scores(:,:)
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: scale(:)
   end type pca_result
contains
   integer function rank_mm(a, tol) result(r)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: singular(:), vt(:,:)
      real(dp) :: threshold
      integer :: n,p,info
      n=size(a,1);p=size(a,2)
      if(min(n,p)==0) then
         r=0
         return
      end if
      call compute_svd(a,singular,vt,info)
      if(info/=0) then
         r=0
         return
      end if
      if(present(tol)) then
         threshold=max(0.0_dp,tol)
      else
         threshold=real(max(n,p),dp)*epsilon(1.0_dp)*singular(1)
      end if
      r=count(singular>=threshold)
   end function rank_mm

   subroutine classical_pca(x, result, center_data, scale_data, signflip, return_scores)
      real(dp), intent(in) :: x(:,:)
      type(pca_result), intent(out) :: result
      logical, intent(in), optional :: center_data, scale_data, signflip, return_scores
      real(dp), allocatable :: z(:,:), singular(:), vt(:,:)
      real(dp) :: threshold
      logical :: do_center,do_scale,do_flip,do_scores
      integer :: n,p,info,j,k,r
      n=size(x,1);p=size(x,2)
      if(n<=1 .or. p<1) error stop "classical_pca: invalid dimensions"
      do_center=.true.;do_scale=.false.;do_flip=.true.;do_scores=.false.
      if(present(center_data)) do_center=center_data
      if(present(scale_data)) do_scale=scale_data
      if(present(signflip)) do_flip=signflip
      if(present(return_scores)) do_scores=return_scores
      allocate(result%center(p),result%scale(p),z(n,p))
      if(do_center) then
         do j=1,p
            result%center(j)=sum(x(:,j))/real(n,dp)
         end do
      else
         result%center=0.0_dp
      end if
      do j=1,p
         z(:,j)=x(:,j)-result%center(j)
      end do
      result%scale=1.0_dp
      if(do_scale) then
         do j=1,p
            result%scale(j)=sqrt(sum(z(:,j)**2)/real(n-1,dp))
            if(result%scale(j)<=sqrt(epsilon(1.0_dp))) result%scale(j)=1.0_dp
            z(:,j)=z(:,j)/result%scale(j)
         end do
      end if
      call compute_svd(z,singular,vt,info)
      if(info/=0) error stop "classical_pca: SVD failed"
      threshold=real(max(n,p),dp)*epsilon(1.0_dp)*singular(1)
      r=count(singular>=threshold)
      result%rank=r
      allocate(result%eigenvalues(r),result%loadings(p,r))
      do k=1,r
         result%eigenvalues(k)=singular(k)**2/real(n-1,dp)
         result%loadings(:,k)=vt(k,:)
         if(do_flip) then
            j=maxloc(abs(result%loadings(:,k)),dim=1)
            if(result%loadings(j,k)<0.0_dp) result%loadings(:,k)=-result%loadings(:,k)
         end if
      end do
      if(do_scores) then
         allocate(result%scores(n,r))
         if(r>0) result%scores=matmul(z,result%loadings)
      else
         allocate(result%scores(0,0))
      end if
   end subroutine classical_pca

   subroutine compute_svd(a,singular,vt,info)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::singular(:),vt(:,:)
      integer,intent(out)::info
      real(dp),allocatable::aa(:,:),work(:),u(:,:)
      real(dp)::query(1)
      integer::m,n,k,lwork
      m=size(a,1);n=size(a,2);k=min(m,n)
      allocate(aa(m,n),singular(k),u(1,1),vt(k,n))
      aa=a
      call dgesvd('N','S',m,n,aa,m,singular,u,1,vt,k,query,-1,info)
      if(info/=0) return
      lwork=max(1,int(query(1)))
      allocate(work(lwork))
      aa=a
      call dgesvd('N','S',m,n,aa,m,singular,u,1,vt,k,work,lwork,info)
   end subroutine compute_svd
end module robustbase_pca
