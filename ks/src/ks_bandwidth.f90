! SPDX-License-Identifier: GPL-2.0-only
module ks_bandwidth
   use ks_kinds, only: dp, pi
   use ks_linalg, only: covariance_matrix, spd_inverse, determinant_spd, matrix_sqrt, trace_matrix
   use ks_normal, only: mvn_pdf, normal_derivative, psins_1d, mvn_derivative_tensor
   use ks_utils, only: vech, invvech
   use ks_optimize, only: golden_minimize, nelder_mead
   implicit none
   private
   public :: hns_matrix, hns_diag, hns_1d, lscv_value, hlscv_1d, hlscv_matrix
   public :: kfe_1d, hpi, hscv, hpi_matrix, hpi_diag
contains
   subroutine hns_matrix(x,H,deriv_order)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::H(size(x,2),size(x,2))
      integer,intent(in),optional::deriv_order
      real(dp)::s(size(x,2),size(x,2)),fac
      integer::n,d,r
      n=size(x,1);d=size(x,2);r=0;if(present(deriv_order))r=deriv_order
      call covariance_matrix(x,s)
      fac=(4.0_dp/(real(n,dp)*real(d+2*r+2,dp)))**(2.0_dp/real(d+2*r+4,dp))
      H=fac*s
   end subroutine hns_matrix

   function hns_1d(x,deriv_order) result(h)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::deriv_order
      real(dp)::h,mu,sd
      integer::n,r
      n=size(x);r=0;if(present(deriv_order))r=deriv_order
      mu=sum(x)/real(n,dp)
      if(n>1)then;sd=sqrt(sum((x-mu)**2)/real(n-1,dp));else;sd=0.0_dp;end if
      h=(4.0_dp/(real(n,dp)*real(2*r+3,dp)))**(1.0_dp/real(2*r+5,dp))*sd
   end function hns_1d

   subroutine hns_diag(x,H)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::H(size(x,2),size(x,2))
      real(dp)::s(size(x,2),size(x,2)),dmat(size(x,2),size(x,2)),dinv(size(x,2),size(x,2))
      real(dp)::delta(size(x,2),size(x,2)),deltainv(size(x,2),size(x,2)),detd,fac
      integer::d,n,i,info
      n=size(x,1);d=size(x,2);call covariance_matrix(x,s);dmat=0.0_dp
      do i=1,d;dmat(i,i)=max(s(i,i),tiny(1.0_dp));end do
      call spd_inverse(dmat,dinv,info)
      delta=matmul(dinv,s);call spd_inverse(0.5_dp*(delta+transpose(delta)),deltainv,info)
      detd=determinant_spd(0.5_dp*(delta+transpose(delta)),info)
      fac=((4.0_dp*real(d,dp)*sqrt(max(detd,tiny(1.0_dp))))/ &
         (2.0_dp*trace_matrix(matmul(deltainv,deltainv))+trace_matrix(deltainv)**2))**(2.0_dp/real(d+4,dp)) &
         *real(n,dp)**(-2.0_dp/real(d+4,dp))
      H=fac*dmat
   end subroutine hns_diag

   function lscv_value(x,H) result(v)
      real(dp),intent(in)::x(:,:),H(:,:)
      real(dp)::v,h2(size(H,1),size(H,2)),zero(size(x,2)),term1,term2
      integer::i,j,n
      n=size(x,1);zero=0.0_dp;h2=2.0_dp*H;term1=0.0_dp;term2=0.0_dp
      do i=1,n
         do j=1,n
            term1=term1+mvn_pdf(x(i,:)-x(j,:),zero,h2)
            if(i/=j)term2=term2+mvn_pdf(x(i,:)-x(j,:),zero,H)
         end do
      end do
      v=term1/real(n*n,dp)-2.0_dp*term2/real(n*(n-1),dp)
   end function lscv_value

   function kfe_1d(x,g,deriv_order,inc) result(v)
      real(dp),intent(in)::x(:),g
      integer,intent(in)::deriv_order
      logical,intent(in),optional::inc
      real(dp)::v
      integer::i,j,n,cnt
      logical::include
      include=.true.;if(present(inc))include=inc
      n=size(x);v=0.0_dp;cnt=0
      do i=1,n
         do j=1,n
            if(.not.include .and. i==j)cycle
            v=v+normal_derivative(x(i)-x(j),0.0_dp,g,deriv_order);cnt=cnt+1
         end do
      end do
      if(cnt>0)v=v/real(cnt,dp)
   end function kfe_1d

   function lscv_1d_value(x,h,deriv_order) result(v)
      real(dp),intent(in)::x(:),h
      integer,intent(in)::deriv_order
      real(dp)::v
      integer::r
      r=deriv_order
      v=(-1.0_dp)**r*(kfe_1d(x,sqrt(2.0_dp)*h,2*r,.true.)-2.0_dp*kfe_1d(x,h,2*r,.false.))
   end function lscv_1d_value

   function min_positive_distance(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v,d
      integer::i,j
      v=huge(1.0_dp)
      do i=1,size(x)-1;do j=i+1,size(x);d=abs(x(i)-x(j));if(d>0.0_dp)v=min(v,d);end do;end do
      if(v>=0.5_dp*huge(1.0_dp))v=max(hns_1d(x)*1e-3_dp,tiny(1.0_dp))
   end function min_positive_distance

   function hlscv_1d(x,deriv_order) result(h)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::deriv_order
      real(dp)::h,fmin,lo,hi,hn
      integer::r
      r=0;if(present(deriv_order))r=deriv_order
      hn=hns_1d(x,r);lo=min_positive_distance(x);hi=max(2.0_dp*hn,lo*1.01_dp)
      call golden_minimize(obj,lo,hi,h,fmin,1e-9_dp,500)
   contains
      function obj(z) result(v)
         real(dp),intent(in)::z
         real(dp)::v
         if(z<=0.0_dp)then;v=huge(1.0_dp);else;v=lscv_1d_value(x,z,r);end if
      end function obj
   end function hlscv_1d

   subroutine hlscv_matrix(x,H,Hstart,maxiter)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::H(size(x,2),size(x,2))
      real(dp),intent(in),optional::Hstart(:,:)
      integer,intent(in),optional::maxiter
      real(dp)::h0(size(x,2),size(x,2)),root(size(x,2),size(x,2)),fbest
      real(dp),allocatable::p0(:),pbest(:),b(:,:)
      integer::info,nit
      if(present(Hstart))then;h0=Hstart;else;call hns_matrix(x,h0);end if
      call matrix_sqrt(h0,root,info);p0=vech(root);allocate(pbest(size(p0)))
      nit=1500;if(present(maxiter))nit=maxiter
      call nelder_mead(obj,p0,pbest,fbest,1e-8_dp,nit,0.08_dp)
      call invvech(pbest,b,info);H=matmul(b,b)
   contains
      function obj(p) result(v)
         real(dp),intent(in)::p(:)
         real(dp)::v,det
         real(dp),allocatable::bb(:,:)
         integer::ii
         call invvech(p,bb,ii)
         if(ii/=0)then;v=huge(1.0_dp);return;end if
         det=determinant_spd(matmul(bb,bb),ii)
         if(ii/=0 .or. det<=sqrt(tiny(1.0_dp)))then;v=huge(1.0_dp);return;end if
         v=lscv_value(x,matmul(bb,bb))
      end function obj
   end subroutine hlscv_matrix

   function hpi(x,nstage,deriv_order) result(h)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::nstage,deriv_order
      real(dp)::h,sd,mu,k2r4,k2r6,psi2r8,psi2r6,psi2r4,g1,g2,mr
      integer::n,ns,r
      n=size(x);ns=2;if(present(nstage))ns=nstage;r=0;if(present(deriv_order))r=deriv_order
      mu=sum(x)/real(n,dp);sd=sqrt(sum((x-mu)**2)/real(max(n-1,1),dp))
      k2r4=normal_derivative(0.0_dp,0.0_dp,1.0_dp,2*r+4)
      k2r6=normal_derivative(0.0_dp,0.0_dp,1.0_dp,2*r+6)
      if(ns>=2)then
         psi2r8=psins_1d(2*r+8,sd)
         g1=(2.0_dp*k2r6/(-psi2r8*real(n,dp)))**(1.0_dp/real(2*r+9,dp))
         psi2r6=kfe_1d(x,g1,2*r+6,.true.)
      else
         psi2r6=psins_1d(2*r+6,sd)
      end if
      g2=(2.0_dp*k2r4/(-psi2r6*real(n,dp)))**(1.0_dp/real(2*r+7,dp))
      psi2r4=kfe_1d(x,g2,2*r+4,.true.)
      mr=psins_1d(2*r,1.0_dp)
      h=(real(2*r+1,dp)*mr/(psi2r4*real(n,dp)))**(1.0_dp/real(2*r+5,dp))
      if(.not.(h>0.0_dp))h=hns_1d(x,r)
   end function hpi

   function hscv(x,nstage) result(h)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::nstage
      real(dp)::h,sd,mu,hn,hmin,hmax,psi6,psi10,g1,g2,g3,g4,psi4,psi8,c,fmin
      integer::n,ns
      n=size(x);ns=2;if(present(nstage))ns=nstage
      mu=sum(x)/real(n,dp);sd=sqrt(sum((x-mu)**2)/real(max(n-1,1),dp));hn=hns_1d(x)
      hmin=0.1_dp*hn;hmax=2.0_dp*hn
      if(ns==1)then
         psi6=psins_1d(6,sd);psi10=psins_1d(10,sd)
      else
         g1=(2.0_dp/(7.0_dp*real(n,dp)))**(1.0_dp/9.0_dp)*sqrt(2.0_dp)*sd
         g2=(2.0_dp/(11.0_dp*real(n,dp)))**(1.0_dp/13.0_dp)*sqrt(2.0_dp)*sd
         psi6=kfe_1d(x,g1,6,.true.);psi10=kfe_1d(x,g2,10,.true.)
      end if
      g3=(-6.0_dp/(sqrt(2.0_dp*pi)*psi6*real(n,dp)))**(1.0_dp/7.0_dp)
      g4=(-210.0_dp/(sqrt(2.0_dp*pi)*psi10*real(n,dp)))**(1.0_dp/11.0_dp)
      psi4=kfe_1d(x,g3,4,.true.);psi8=kfe_1d(x,g4,8,.true.)
      c=(441.0_dp/(64.0_dp*pi))**(1.0_dp/18.0_dp)*(4.0_dp*pi)**(-1.0_dp/5.0_dp) &
         *abs(psi4)**(-2.0_dp/5.0_dp)*abs(psi8)**(-1.0_dp/9.0_dp)
      call golden_minimize(obj,hmin,hmax,h,fmin,1e-8_dp,500)
   contains
      function obj(z) result(v)
         real(dp),intent(in)::z
         real(dp)::v,g,bias
         g=c*real(n,dp)**(-23.0_dp/45.0_dp)*z**(-2.0_dp)
         bias=kfe_1d(x,sqrt(2*z*z+2*g*g),0,.true.)-2*kfe_1d(x,sqrt(z*z+2*g*g),0,.true.) &
              +kfe_1d(x,sqrt(2*g*g),0,.true.)
         bias=max(bias,0.0_dp)
         v=1.0_dp/(real(n,dp)*z*sqrt(4.0_dp*pi))+bias
      end function obj
   end function hscv

   subroutine hpi_matrix(x,H,Hstart,maxiter)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::H(size(x,2),size(x,2))
      real(dp),intent(in),optional::Hstart(:,:)
      integer,intent(in),optional::maxiter
      integer::n,d,i,j,a,b,c,e,idx,row,col,info,nit
      real(dp)::G(size(x,2),size(x,2)),h0(size(x,2),size(x,2)),root(size(x,2),size(x,2)), &
                zero(size(x,2)),delta(size(x,2)),rk,fbest
      real(dp),allocatable::psi(:,:),tmp(:),p0(:),pbest(:),bb(:,:)
      n=size(x,1);d=size(x,2);zero=0.0_dp
      call hns_matrix(x,G,2)
      allocate(psi(d*d,d*d));psi=0.0_dp
      do i=1,n
         do j=1,n
            delta=x(i,:)-x(j,:)
            call mvn_derivative_tensor(delta,zero,G,4,tmp,info)
            if(info/=0) error stop 'hpi_matrix: pilot derivative'
            do b=1,d; do a=1,d
               row=a+(b-1)*d
               do e=1,d; do c=1,d
                  col=c+(e-1)*d
                  idx=1+(a-1)+d*(b-1)+d*d*(c-1)+d*d*d*(e-1)
                  psi(row,col)=psi(row,col)+tmp(idx)
               end do; end do
            end do; end do
         end do
      end do
      psi=psi/real(n*n,dp)
      if(present(Hstart))then;h0=Hstart;else;call hns_matrix(x,h0);end if
      call matrix_sqrt(h0,root,info);p0=vech(root);allocate(pbest(size(p0)))
      nit=1800;if(present(maxiter))nit=maxiter
      rk=(4.0_dp*pi)**(-0.5_dp*real(d,dp))
      call nelder_mead(obj,p0,pbest,fbest,1.0e-8_dp,nit,0.08_dp)
      call invvech(pbest,bb,info);H=matmul(bb,bb)
   contains
      function obj(p) result(v)
         real(dp),intent(in)::p(:)
         real(dp)::v,det
         real(dp),allocatable::Bmat(:,:),HH(:,:)
         real(dp)::vhloc(d*d)
         integer::ii,jj,kk
         call invvech(p,Bmat,ii)
         if(ii/=0)then;v=huge(1.0_dp);return;end if
         HH=matmul(Bmat,Bmat);det=determinant_spd(HH,ii)
         if(ii/=0.or.det<=sqrt(tiny(1.0_dp)))then;v=huge(1.0_dp);return;end if
         kk=0
         do jj=1,d;do ii=1,d;kk=kk+1;vhloc(kk)=HH(ii,jj);end do;end do
         v=rk/(real(n,dp)*sqrt(det))+0.25_dp*dot_product(vhloc,matmul(psi,vhloc))
      end function obj
   end subroutine hpi_matrix

   subroutine hpi_diag(x,H,maxiter)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::H(size(x,2),size(x,2))
      integer,intent(in),optional::maxiter
      real(dp)::Hfull(size(x,2),size(x,2))
      integer::i
      call hpi_matrix(x,Hfull,maxiter=maxiter)
      H=0.0_dp
      do i=1,size(x,2);H(i,i)=max(Hfull(i,i),tiny(1.0_dp));end do
   end subroutine hpi_diag

end module ks_bandwidth
