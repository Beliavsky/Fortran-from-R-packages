! Translation of distr/src/ks.c, which was taken from R Core stats/src/ks.c.
! Original copyright (C) 1999-2012 The R Core Team.
! SPDX-License-Identifier: GPL-2.0-or-later
!
! This module intentionally remains separately licensed from the LGPL-3.0-only
! distr-derived modules. See LICENSES.md.
module distr_ks
   use distr_kinds, only : dp, pi
   implicit none
   private
   public :: p_ks2_asymptotic, p_smirnov2x, p_kolmogorov2x

contains

   elemental real(dp) function p_ks2_asymptotic(x,tol) result(p)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: tol
      real(dp) :: epsv,z,w,s,old,newv,sgn
      integer :: k,kmax
      epsv=1.0e-12_dp; if (present(tol)) epsv=tol
      if (x<=0.0_dp) then; p=0.0_dp; return; end if
      kmax=int(sqrt(2.0_dp-log(epsv)))
      if (x<1.0_dp) then
         z=-(0.5_dp*pi)*(0.25_dp*pi)/(x*x)
         w=log(x); s=0.0_dp
         do k=1,kmax-1,2; s=s+exp(real(k*k,dp)*z-w); end do
         p=s*sqrt(2.0_dp*pi)
      else
         z=-2.0_dp*x*x; sgn=-1.0_dp; k=1; old=0.0_dp; newv=1.0_dp
         do while(abs(old-newv)>epsv)
            old=newv; newv=newv+2.0_dp*sgn*exp(z*real(k*k,dp)); sgn=-sgn; k=k+1
            if (k>1000000) exit
         end do
         p=newv
      end if
      p=min(1.0_dp,max(0.0_dp,p))
   end function p_ks2_asymptotic

   real(dp) function p_smirnov2x(x,m,n) result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: m,n
      integer :: mm,nn,i,j,it
      real(dp) :: md,nd,q,w
      real(dp), allocatable :: u(:)
      mm=m; nn=n
      if (mm<=0.or.nn<=0) then; p=0.0_dp; return; end if
      if (mm>nn) then; it=nn; nn=mm; mm=it; end if
      md=real(mm,dp); nd=real(nn,dp)
      q=(0.5_dp+floor(x*md*nd-1.0e-7_dp))/(md*nd)
      allocate(u(0:nn))
      do j=0,nn; u(j)=merge(0.0_dp,1.0_dp,real(j,dp)/nd>q); end do
      do i=1,mm
         w=real(i,dp)/real(i+nn,dp)
         if (real(i,dp)/md>q) then; u(0)=0.0_dp; else; u(0)=w*u(0); end if
         do j=1,nn
            if (abs(real(i,dp)/md-real(j,dp)/nd)>q) then
               u(j)=0.0_dp
            else
               u(j)=w*u(j)+u(j-1)
            end if
         end do
      end do
      p=u(nn)
   end function p_smirnov2x

   real(dp) function p_kolmogorov2x(d,n) result(p)
      real(dp), intent(in) :: d
      integer, intent(in) :: n
      integer :: k,m,i,j,g,e_h,e_q
      real(dp) :: h,s
      real(dp), allocatable :: hh(:,:),qq(:,:)
      if (d<=0.0_dp) then; p=0.0_dp; return; end if
      if (d>=1.0_dp) then; p=1.0_dp; return; end if
      if (n<=0) then; p=0.0_dp; return; end if
      s=d*d*real(n,dp)
      if (s>7.24_dp .or. (s>3.76_dp .and. n>99)) then
         p=1.0_dp-2.0_dp*exp(-(2.000071_dp+0.331_dp/sqrt(real(n,dp))+1.409_dp/real(n,dp))*s)
         p=min(1.0_dp,max(0.0_dp,p)); return
      end if
      k=int(real(n,dp)*d)+1; m=2*k-1; h=real(k,dp)-real(n,dp)*d
      allocate(hh(m,m),qq(m,m)); hh=0.0_dp; qq=0.0_dp
      do i=1,m
         do j=1,m
            if ((i-1)-(j-1)+1>=0) hh(i,j)=1.0_dp
         end do
      end do
      do i=1,m
         hh(i,1)=hh(i,1)-h**i
         hh(m,i)=hh(m,i)-h**(m-i+1)
      end do
      if (2.0_dp*h-1.0_dp>0.0_dp) hh(m,1)=hh(m,1)+(2.0_dp*h-1.0_dp)**m
      do i=1,m
         do j=1,m
            if (i-j+1>0) then
               do g=1,i-j+1; hh(i,j)=hh(i,j)/real(g,dp); end do
            end if
         end do
      end do
      e_h=0
      call matrix_power_scaled(hh,e_h,qq,e_q,n)
      s=qq(k,k)
      do i=1,n
         s=s*real(i,dp)/real(n,dp)
         if (s<1.0e-140_dp) then; s=s*1.0e140_dp; e_q=e_q-140; end if
      end do
      p=s*10.0_dp**e_q
      p=min(1.0_dp,max(0.0_dp,p))
   end function p_kolmogorov2x

   subroutine matrix_multiply(a,b,c)
      real(dp), intent(in) :: a(:,:),b(:,:)
      real(dp), intent(out) :: c(size(a,1),size(b,2))
      integer :: i,j,k,n
      n=size(a,2); c=0.0_dp
      do i=1,size(a,1); do j=1,size(b,2); do k=1,n; c(i,j)=c(i,j)+a(i,k)*b(k,j); end do; end do; end do
   end subroutine matrix_multiply

   recursive subroutine matrix_power_scaled(a,e_a,v,e_v,n)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: e_a,n
      real(dp), intent(out) :: v(size(a,1),size(a,2))
      integer, intent(out) :: e_v
      real(dp), allocatable :: b(:,:)
      integer :: e_b
      if (n==1) then; v=a; e_v=e_a; return; end if
      call matrix_power_scaled(a,e_a,v,e_v,n/2)
      allocate(b(size(a,1),size(a,2))); call matrix_multiply(v,v,b); e_b=2*e_v
      if (mod(n,2)==0) then; v=b; e_v=e_b
      else; call matrix_multiply(a,b,v); e_v=e_a+e_b; end if
      if (v(size(v,1)/2+1,size(v,2)/2+1)>1.0e140_dp) then; v=v*1.0e-140_dp; e_v=e_v+140; end if
   end subroutine matrix_power_scaled

end module distr_ks
