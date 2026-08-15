! Multivariate-normal density and rectangle probabilities for copula likelihoods.
! Rectangle probabilities use the Genz sequential conditioning transform with
! deterministic Halton points and antithetic pairing.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_mvn_v07
   use gamlss_kinds, only : dp, pi, log2pi
   use gamlss_special, only : normal_cdf, normal_quantile
   use gamlss_linalg, only : cholesky_factor, invert_matrix
   implicit none
   private
   public :: mvn_logpdf, mvn_rectangle_probability, mvn_conditional
contains

   real(dp) function mvn_logpdf(x,mean,cov,status) result(lp)
      real(dp),intent(in) :: x(:),mean(:),cov(:,:)
      integer,intent(out) :: status
      real(dp),allocatable :: inv(:,:),d(:),l(:,:)
      real(dp) :: logdet
      integer :: i,n,istat
      n=size(x);status=0;lp=-huge(1.0_dp)
      if(size(mean)/=n.or.any(shape(cov)/=[n,n]))then;status=1;return;end if
      call stable_cholesky(cov,l,istat)
      if(istat/=0)then;status=2;return;end if
      logdet=0.0_dp
      do i=1,n;logdet=logdet+2.0_dp*log(l(i,i));end do
      call invert_matrix(cov,inv,istat)
      if(istat/=0)then;status=3;return;end if
      d=x-mean
      lp=-0.5_dp*(real(n,dp)*log2pi+logdet+dot_product(d,matmul(inv,d)))
   end function mvn_logpdf

   subroutine mvn_conditional(cov,continuous,discrete,z_cont,mean_d,cov_d,status)
      real(dp),intent(in) :: cov(:,:),z_cont(:)
      integer,intent(in) :: continuous(:),discrete(:)
      real(dp),allocatable,intent(out) :: mean_d(:),cov_d(:,:)
      integer,intent(out) :: status
      real(dp),allocatable :: rcc(:,:),rdc(:,:),rdd(:,:),invcc(:,:)
      integer :: nc,nd,i,j,istat
      status=0;nc=size(continuous);nd=size(discrete)
      if(size(cov,1)/=size(cov,2).or.nc/=size(z_cont))then;status=1;return;end if
      allocate(mean_d(nd),cov_d(nd,nd));mean_d=0.0_dp;cov_d=0.0_dp
      if(nd==0)return
      allocate(rdd(nd,nd))
      do i=1,nd;do j=1,nd;rdd(i,j)=cov(discrete(i),discrete(j));end do;end do
      if(nc==0)then;cov_d=rdd;return;end if
      allocate(rcc(nc,nc),rdc(nd,nc))
      do i=1,nc;do j=1,nc;rcc(i,j)=cov(continuous(i),continuous(j));end do;end do
      do i=1,nd;do j=1,nc;rdc(i,j)=cov(discrete(i),continuous(j));end do;end do
      call invert_matrix(rcc,invcc,istat)
      if(istat/=0)then;status=2;return;end if
      mean_d=matmul(rdc,matmul(invcc,z_cont))
      cov_d=rdd-matmul(matmul(rdc,invcc),transpose(rdc))
      call symmetrize_and_jitter(cov_d)
   end subroutine mvn_conditional

   real(dp) function mvn_rectangle_probability(lower,upper,mean,cov,status,n_qmc,error_estimate) result(prob)
      real(dp),intent(in) :: lower(:),upper(:),mean(:),cov(:,:)
      integer,intent(out) :: status
      integer,intent(in),optional :: n_qmc
      real(dp),intent(out),optional :: error_estimate
      real(dp),allocatable :: l(:,:),a(:),b(:),z(:)
      real(dp) :: s1,s2,val,val2
      integer :: n,m,k,istat
      n=size(lower);status=0;prob=0.0_dp
      if(present(error_estimate))error_estimate=huge(1.0_dp)
      if(size(upper)/=n.or.size(mean)/=n.or.any(shape(cov)/=[n,n]))then;status=1;return;end if
      if(any(upper<=lower))then;prob=0.0_dp;if(present(error_estimate))error_estimate=0.0_dp;return;end if
      if(n==0)then;prob=1.0_dp;if(present(error_estimate))error_estimate=0.0_dp;return;end if
      call stable_cholesky(cov,l,istat)
      if(istat/=0)then;status=2;return;end if
      allocate(a(n),b(n),z(n));a=lower-mean;b=upper-mean
      if(n==1)then
         prob=normal_cdf(b(1)/l(1,1))-normal_cdf(a(1)/l(1,1))
         prob=max(0.0_dp,min(1.0_dp,prob))
         if(present(error_estimate))error_estimate=0.0_dp
         return
      end if
      m=2048;if(present(n_qmc))m=max(64,n_qmc)
      if(mod(m,2)/=0)m=m+1
      s1=0.0_dp;s2=0.0_dp
      do k=1,m/2
         call genz_integrand(k,m/2,a,b,l,z,.false.,val)
         call genz_integrand(k,m/2,a,b,l,z,.true.,val2)
         s1=s1+val;s2=s2+val2
      end do
      prob=(s1+s2)/real(m,dp)
      prob=max(0.0_dp,min(1.0_dp,prob))
      if(present(error_estimate))error_estimate=abs(s1-s2)/real(m,dp)
   end function mvn_rectangle_probability

   subroutine genz_integrand(k,m,a,b,l,z,antithetic,val)
      integer,intent(in) :: k,m
      real(dp),intent(in) :: a(:),b(:),l(:,:)
      real(dp),intent(inout) :: z(:)
      logical,intent(in) :: antithetic
      real(dp),intent(out) :: val
      real(dp) :: lo,hi,w,u,cond
      integer :: i,n
      n=size(a);z=0.0_dp;val=1.0_dp
      do i=1,n
         cond=0.0_dp
         if(i>1)cond=dot_product(l(i,1:i-1),z(1:i-1))
         lo=normal_cdf((a(i)-cond)/l(i,i));hi=normal_cdf((b(i)-cond)/l(i,i))
         w=max(0.0_dp,hi-lo);val=val*w
         if(w<=0.0_dp)return
         if(i<n)then
            u=halton_value(k,i,m)
            if(antithetic)u=1.0_dp-u
            u=min(1.0_dp-1.0e-13_dp,max(1.0e-13_dp,lo+u*w))
            z(i)=normal_quantile(u)
         end if
      end do
   end subroutine genz_integrand

   real(dp) function halton_value(k,dim,m) result(u)
      integer,intent(in) :: k,dim,m
      integer,parameter :: primes(24)=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89]
      integer :: kk,base,digit
      real(dp) :: f
      base=primes(1+mod(dim-1,size(primes)))
      ! Cranley-like deterministic shift depends on dimension and run length.
      kk=k;u=0.0_dp;f=1.0_dp/real(base,dp)
      do while(kk>0)
         digit=mod(kk,base);u=u+f*real(digit,dp);kk=kk/base;f=f/real(base,dp)
      end do
      u=modulo(u+modulo(0.6180339887498949_dp*real(dim+m,dp),1.0_dp),1.0_dp)
      u=min(1.0_dp-1.0e-13_dp,max(1.0e-13_dp,u))
   end function halton_value

   subroutine stable_cholesky(a,l,status)
      real(dp),intent(in) :: a(:,:)
      real(dp),allocatable,intent(out) :: l(:,:)
      integer,intent(out) :: status
      real(dp),allocatable :: work(:,:)
      real(dp) :: scale,jitter
      integer :: k,istat,n
      n=size(a,1);status=0
      if(size(a,2)/=n)then;allocate(l(0,0));status=1;return;end if
      work=0.5_dp*(a+transpose(a));scale=max(1.0_dp,maxval(abs(work)))
      do k=0,8
         if(k==0)then;jitter=0.0_dp;else;jitter=scale*10.0_dp**(-13+k);end if
         if(k>0)then;work=0.5_dp*(a+transpose(a));work=work+jitter*identity_matrix(n);end if
         call cholesky_factor(work,l,istat)
         if(istat==0)then;status=0;return;end if
         if(allocated(l))deallocate(l)
      end do
      allocate(l(0,0));status=2
   end subroutine stable_cholesky

   subroutine symmetrize_and_jitter(a)
      real(dp),intent(inout) :: a(:,:)
      real(dp) :: eps
      integer :: i
      a=0.5_dp*(a+transpose(a));eps=max(1.0e-12_dp,1.0e-10_dp*max(1.0_dp,maxval(abs(a))))
      do i=1,size(a,1);a(i,i)=max(a(i,i),eps);end do
   end subroutine symmetrize_and_jitter

   pure function identity_matrix(n) result(a)
      integer,intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a=0.0_dp;do i=1,n;a(i,i)=1.0_dp;end do
   end function identity_matrix
end module gamlss_mvn_v07
