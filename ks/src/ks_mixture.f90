! SPDX-License-Identifier: GPL-2.0-only
module ks_mixture
   use ks_kinds, only: dp, pi
   use ks_normal, only: mvn_pdf, mvn_derivative_tensor
   use ks_linalg, only: determinant_spd, covariance_matrix, matrix_sqrt
   use ks_utils, only: vech, invvech
   use ks_optimize, only: nelder_mead, golden_minimize
   implicit none
   private
   public :: mise_normal_mixture, amise_normal_mixture, ise_normal_mixture
   public :: hmise_normal_mixture, hamise_normal_mixture
   public :: mise_normal_mixture_1d, amise_normal_mixture_1d, normal_mixture_modes, normal_mixture_moments
contains
   function weighted_omega(mus,sigmas,props,H,a) result(v)
      real(dp),intent(in)::mus(:,:),sigmas(:,:,:),props(:),H(:,:),a
      real(dp)::v,cov(size(H,1),size(H,2)),zero(size(H,1))
      integer::i,j
      zero=0.0_dp;v=0.0_dp
      do i=1,size(props)
         do j=1,size(props)
            cov=sigmas(:,:,i)+sigmas(:,:,j)+a*H
            v=v+props(i)*props(j)*mvn_pdf(mus(i,:)-mus(j,:),zero,cov)
         end do
      end do
   end function weighted_omega

   function mise_normal_mixture(H,mus,sigmas,props,samp) result(v)
      real(dp),intent(in)::H(:,:),mus(:,:),sigmas(:,:,:),props(:)
      integer,intent(in)::samp
      real(dp)::v,detH
      integer::info,d
      d=size(H,1);detH=determinant_spd(H,info)
      if(info/=0 .or. detH<=0.0_dp .or. samp<=0)then;v=huge(1.0_dp);return;end if
      v=1.0_dp/(real(samp,dp)*(4.0_dp*pi)**(0.5_dp*real(d,dp))*sqrt(detH)) &
        +(1.0_dp-1.0_dp/real(samp,dp))*weighted_omega(mus,sigmas,props,H,2.0_dp) &
        -2.0_dp*weighted_omega(mus,sigmas,props,H,1.0_dp) &
        +weighted_omega(mus,sigmas,props,H,0.0_dp)
   end function mise_normal_mixture

   function amise_bias_pair(mu,cov,H) result(v)
      real(dp),intent(in)::mu(:),cov(:,:),H(:,:)
      real(dp)::v,zero(size(mu))
      real(dp),allocatable::d4(:)
      integer::d,a,b,c,e,code
      zero=0.0_dp;d=size(mu)
      call mvn_derivative_tensor(mu,zero,cov,4,d4)
      v=0.0_dp;code=0
      do e=1,d
         do c=1,d
            do b=1,d
               do a=1,d
                  code=code+1
                  v=v+H(a,b)*H(c,e)*d4(code)
               end do
            end do
         end do
      end do
   end function amise_bias_pair

   function amise_normal_mixture(H,mus,sigmas,props,samp) result(v)
      real(dp),intent(in)::H(:,:),mus(:,:),sigmas(:,:,:),props(:)
      integer,intent(in)::samp
      real(dp)::v,detH,bias
      integer::info,d,i,j
      d=size(H,1);detH=determinant_spd(H,info)
      if(info/=0 .or. detH<=0.0_dp .or. samp<=0)then;v=huge(1.0_dp);return;end if
      bias=0.0_dp
      do i=1,size(props)
         do j=1,size(props)
            bias=bias+props(i)*props(j)*amise_bias_pair(mus(i,:)-mus(j,:),sigmas(:,:,i)+sigmas(:,:,j),H)
         end do
      end do
      v=1.0_dp/(real(samp,dp)*(4.0_dp*pi)**(0.5_dp*real(d,dp))*sqrt(detH))+0.25_dp*bias
   end function amise_normal_mixture

   function ise_normal_mixture(x,H,mus,sigmas,props) result(v)
      real(dp),intent(in)::x(:,:),H(:,:),mus(:,:),sigmas(:,:,:),props(:)
      real(dp)::v,t1,t2,t3,cov(size(H,1),size(H,2)),zero(size(H,1))
      integer::i,j,k,n
      n=size(x,1);zero=0.0_dp;t1=0.0_dp;t2=0.0_dp;t3=0.0_dp
      do i=1,n;do j=1,n;t1=t1+mvn_pdf(x(i,:)-x(j,:),zero,2.0_dp*H);end do;end do
      t1=t1/real(n*n,dp)
      do i=1,n
         do j=1,size(props)
            cov=H+sigmas(:,:,j);t2=t2+props(j)*mvn_pdf(x(i,:)-mus(j,:),zero,cov)
         end do
      end do
      t2=t2/real(n,dp)
      do j=1,size(props);do k=1,size(props)
         cov=sigmas(:,:,j)+sigmas(:,:,k)
         t3=t3+props(j)*props(k)*mvn_pdf(mus(j,:)-mus(k,:),zero,cov)
      end do;end do
      v=t1-2.0_dp*t2+t3
   end function ise_normal_mixture

   subroutine default_start(mus,sigmas,props,samp,H)
      real(dp),intent(in)::mus(:,:),sigmas(:,:,:),props(:)
      integer,intent(in)::samp
      real(dp),intent(out)::H(size(mus,2),size(mus,2))
      real(dp)::mean(size(mus,2)),cov(size(mus,2),size(mus,2)),dx(size(mus,2)),fac
      integer::k,d
      d=size(mus,2);mean=0.0_dp;cov=0.0_dp
      do k=1,size(props);mean=mean+props(k)*mus(k,:);end do
      do k=1,size(props)
         dx=mus(k,:)-mean
         cov=cov+props(k)*(sigmas(:,:,k)+spread(dx,2,d)*spread(dx,1,d))
      end do
      fac=(4.0_dp/(real(samp,dp)*real(d+2,dp)))**(2.0_dp/real(d+4,dp));H=fac*cov
   end subroutine default_start

   subroutine hmise_normal_mixture(mus,sigmas,props,samp,H,Hstart,maxiter)
      real(dp),intent(in)::mus(:,:),sigmas(:,:,:),props(:)
      integer,intent(in)::samp
      real(dp),intent(out)::H(size(mus,2),size(mus,2))
      real(dp),intent(in),optional::Hstart(:,:)
      integer,intent(in),optional::maxiter
      real(dp)::h0(size(mus,2),size(mus,2)),root(size(mus,2),size(mus,2)),fb
      real(dp),allocatable::p0(:),pb(:),b(:,:)
      integer::info,nit
      if(present(Hstart))then;h0=Hstart;else;call default_start(mus,sigmas,props,samp,h0);end if
      call matrix_sqrt(h0,root,info);p0=vech(root);allocate(pb(size(p0)));nit=1500;if(present(maxiter))nit=maxiter
      call nelder_mead(obj,p0,pb,fb,1e-8_dp,nit,0.08_dp);call invvech(pb,b,info);H=matmul(b,b)
   contains
      function obj(p) result(f)
         real(dp),intent(in)::p(:);real(dp)::f
         real(dp),allocatable::bb(:,:);integer::ii
         call invvech(p,bb,ii);if(ii/=0)then;f=huge(1.0_dp);return;end if
         f=mise_normal_mixture(matmul(bb,bb),mus,sigmas,props,samp)
      end function obj
   end subroutine hmise_normal_mixture

   subroutine hamise_normal_mixture(mus,sigmas,props,samp,H,Hstart,maxiter)
      real(dp),intent(in)::mus(:,:),sigmas(:,:,:),props(:)
      integer,intent(in)::samp
      real(dp),intent(out)::H(size(mus,2),size(mus,2))
      real(dp),intent(in),optional::Hstart(:,:)
      integer,intent(in),optional::maxiter
      real(dp)::h0(size(mus,2),size(mus,2)),root(size(mus,2),size(mus,2)),fb
      real(dp),allocatable::p0(:),pb(:),b(:,:)
      integer::info,nit
      if(size(props)==1)then
         H=(4.0_dp/(real(samp,dp)*real(size(mus,2)+2,dp)))**(2.0_dp/real(size(mus,2)+4,dp))*sigmas(:,:,1);return
      end if
      if(present(Hstart))then;h0=Hstart;else;call default_start(mus,sigmas,props,samp,h0);end if
      call matrix_sqrt(h0,root,info);p0=vech(root);allocate(pb(size(p0)));nit=1500;if(present(maxiter))nit=maxiter
      call nelder_mead(obj,p0,pb,fb,1e-8_dp,nit,0.08_dp);call invvech(pb,b,info);H=matmul(b,b)
   contains
      function obj(p) result(f)
         real(dp),intent(in)::p(:);real(dp)::f
         real(dp),allocatable::bb(:,:);integer::ii
         call invvech(p,bb,ii);if(ii/=0)then;f=huge(1.0_dp);return;end if
         f=amise_normal_mixture(matmul(bb,bb),mus,sigmas,props,samp)
      end function obj
   end subroutine hamise_normal_mixture

   function mise_normal_mixture_1d(h,mus,sigmas,props,samp) result(v)
      real(dp),intent(in)::h,mus(:),sigmas(:),props(:)
      integer,intent(in)::samp
      real(dp)::v
      real(dp)::mm(size(mus),1),ss(1,1,size(mus)),hh(1,1)
      integer::k
      mm(:,1)=mus;hh(1,1)=h*h
      do k=1,size(mus);ss(1,1,k)=sigmas(k)**2;end do
      v=mise_normal_mixture(hh,mm,ss,props,samp)
   end function mise_normal_mixture_1d

   function amise_normal_mixture_1d(h,mus,sigmas,props,samp) result(v)
      real(dp),intent(in)::h,mus(:),sigmas(:),props(:)
      integer,intent(in)::samp
      real(dp)::v
      real(dp)::mm(size(mus),1),ss(1,1,size(mus)),hh(1,1)
      integer::k
      mm(:,1)=mus;hh(1,1)=h*h
      do k=1,size(mus);ss(1,1,k)=sigmas(k)**2;end do
      v=amise_normal_mixture(hh,mm,ss,props,samp)
   end function amise_normal_mixture_1d

   subroutine normal_mixture_modes(mus,sigmas,props,modes,maxiter)
      real(dp),intent(in)::mus(:,:),sigmas(:,:,:),props(:)
      real(dp),intent(out)::modes(size(mus,1),size(mus,2))
      integer,intent(in),optional::maxiter
      real(dp)::fbest,tmpmode(size(mus,2))
      integer::i,nit
      nit=1000;if(present(maxiter))nit=maxiter
      if(size(props)/=size(mus,1).or.size(sigmas,3)/=size(props)) error stop 'normal_mixture_modes: shape'
      do i=1,size(props)
         call nelder_mead(obj,mus(i,:),tmpmode,fbest,1.0e-9_dp,nit,0.05_dp)
         modes(i,:)=tmpmode
      end do
   contains
      function obj(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v,zero(size(z))
         integer::j
         zero=0.0_dp;v=0.0_dp
         do j=1,size(props);v=v-props(j)*mvn_pdf(z,mus(j,:),sigmas(:,:,j));end do
      end function obj
   end subroutine normal_mixture_modes

   subroutine normal_mixture_moments(mus,sigmas,props,mean,covariance)
      real(dp),intent(in)::mus(:,:),sigmas(:,:,:),props(:)
      real(dp),intent(out)::mean(size(mus,2)),covariance(size(mus,2),size(mus,2))
      integer::i
      real(dp)::sw
      if(size(props)/=size(mus,1).or.size(sigmas,3)/=size(props)) error stop 'normal_mixture_moments: shape'
      sw=sum(props);if(sw<=0.0_dp)error stop 'normal_mixture_moments: proportions'
      mean=0.0_dp;covariance=0.0_dp
      do i=1,size(props);mean=mean+(props(i)/sw)*mus(i,:);end do
      do i=1,size(props)
         covariance=covariance+(props(i)/sw)*(sigmas(:,:,i)+spread(mus(i,:)-mean,2,size(mean))*spread(mus(i,:)-mean,1,size(mean)))
      end do
   end subroutine normal_mixture_moments

end module ks_mixture
