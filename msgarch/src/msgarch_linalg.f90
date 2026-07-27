! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_linalg
   use msgarch_kinds, only : dp
   implicit none
   private
   abstract interface
      function scalar_function(x) result(f)
         import dp
         real(dp),intent(in)::x(:)
         real(dp)::f
      end function scalar_function
   end interface
   public :: numerical_hessian, invert_matrix, sample_mean_sd
contains
   function numerical_hessian(fun,x,relative_step) result(hessian)
      procedure(scalar_function)::fun
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::relative_step
      real(dp),allocatable::hessian(:,:)
      real(dp),allocatable::xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:),step(:)
      real(dp)::f0
      integer::i,j,n
      n=size(x);allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),step(n))
      step=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x))
      if(present(relative_step))step=max(relative_step*max(1.0_dp,abs(x)),1.0e-7_dp)
      f0=fun(x);hessian=0.0_dp
      do i=1,n
         xp=x;xm=x;xp(i)=xp(i)+step(i);xm(i)=xm(i)-step(i)
         hessian(i,i)=(fun(xp)-2.0_dp*f0+fun(xm))/(step(i)**2)
         do j=i+1,n
            xpp=x;xpm=x;xmp=x;xmm=x
            xpp(i)=xpp(i)+step(i);xpp(j)=xpp(j)+step(j)
            xpm(i)=xpm(i)+step(i);xpm(j)=xpm(j)-step(j)
            xmp(i)=xmp(i)-step(i);xmp(j)=xmp(j)+step(j)
            xmm(i)=xmm(i)-step(i);xmm(j)=xmm(j)-step(j)
            hessian(i,j)=(fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*step(i)*step(j))
            hessian(j,i)=hessian(i,j)
         end do
      end do
   end function numerical_hessian

   subroutine invert_matrix(a,inverse,success)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::inverse(:,:)
      logical,intent(out)::success
      real(dp),allocatable::aug(:,:),temp(:)
      real(dp)::pivot
      integer::n,i,j,k,imax
      n=size(a,1);success=.false.
      if(size(a,2)/=n)return
      allocate(aug(n,2*n),inverse(n,n),temp(2*n));aug=0.0_dp;aug(:,1:n)=a
      do i=1,n;aug(i,n+i)=1.0_dp;end do
      do i=1,n
         imax=i
         do k=i+1,n;if(abs(aug(k,i))>abs(aug(imax,i)))imax=k;end do
         if(abs(aug(imax,i))<1.0e-12_dp)return
         if(imax/=i)then;temp=aug(i,:);aug(i,:)=aug(imax,:);aug(imax,:)=temp;end if
         pivot=aug(i,i);aug(i,:)=aug(i,:)/pivot
         do j=1,n
            if(j/=i)aug(j,:)=aug(j,:)-aug(j,i)*aug(i,:)
         end do
      end do
      inverse=aug(:,n+1:2*n);success=.true.
   end subroutine invert_matrix

   subroutine sample_mean_sd(draws,mean,sd)
      real(dp),intent(in)::draws(:,:)
      real(dp),allocatable,intent(out)::mean(:),sd(:)
      integer::j,n
      n=size(draws,1);allocate(mean(size(draws,2)),sd(size(draws,2)))
      mean=sum(draws,dim=1)/real(n,dp)
      do j=1,size(draws,2)
         if(n>1)then;sd(j)=sqrt(sum((draws(:,j)-mean(j))**2)/real(n-1,dp));else;sd(j)=0.0_dp;end if
      end do
   end subroutine sample_mean_sd
end module msgarch_linalg
