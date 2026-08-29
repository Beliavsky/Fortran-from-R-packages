module spatialextremes_mvprob
   use spatialextremes_base, only: dp, chol_upper
   use r_compat, only: runif1, normal_cdf, qnorm, qchisq
   implicit none
   private
   public :: mvnorm_cdf_qmc, mvstudent_cdf_qmc
   integer,parameter :: primes(100)=[ &
      2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71, &
      73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173, &
      179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281, &
      283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409, &
      419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541]
contains
   function mvnorm_cdf_qmc(bounds,cov,mean,n_per_dim,info) result(prob)
      ! Randomized lattice/antithetic Genz-style integration corresponding
      ! to upstream pmvnorm().  cov may be a covariance or correlation matrix.
      real(dp),intent(in)::bounds(:),cov(:,:)
      real(dp),intent(in),optional::mean(:)
      integer,intent(in),optional::n_per_dim
      integer,intent(out),optional::info
      real(dp)::prob
      integer::n,nmc,j,k,l,istat
      integer,allocatable::idx(:)
      real(dp),allocatable::b(:),a(:,:),r(:,:),delta(:),q(:),w(:),wa(:),e(:),ea(:),y(:),ya(:)
      real(dp)::f,fa,dummy,dummya,u,epsu
      n=size(bounds)
      if(n==0)then
      prob=1.0_dp
      if(present(info))info=0
      return
      end if
      if(size(cov,1)/=n .or. size(cov,2)/=n .or. n>size(primes))then
         prob=0.0_dp
         if(present(info))info=1
         return
      end if
      nmc=n*1000
      if(present(n_per_dim))nmc=max(1,n*max(1,n_per_dim))
      allocate(idx(n),b(n),a(n,n),r(n,n),delta(n),q(n),w(n),wa(n),e(n),ea(n))
      if(n>1)allocate(y(n-1),ya(n-1))
      idx=[(k,k=1,n)]
      b=bounds
      if(present(mean))then
         if(size(mean)/=n)then
         prob=0.0_dp
         if(present(info))info=2
         return
         end if
         b=b-mean
      end if
      call sort_bounds_index(b,idx)
      do j=1,n
      do k=1,n
      a(j,k)=cov(idx(j),idx(k))
      end do
      end do
      call chol_upper(a,r,istat)
      if(istat/=0)then
      prob=0.0_dp
      if(present(info))info=3
      return
      end if
      r=transpose(r) ! lower Cholesky
      do k=1,n
      q(k)=sqrt(real(primes(k),dp))
      delta(k)=runif1()
      end do
      prob=0.0_dp
      epsu=32.0_dp*epsilon(1.0_dp)
      do j=0,nmc-1
         do k=1,n
            u=real(j,dp)*q(k)+delta(k)
            u=u-floor(u)
            w(k)=abs(2.0_dp*u-1.0_dp)
            wa(k)=1.0_dp-w(k)
            w(k)=min(1.0_dp-epsu,max(epsu,w(k)))
            wa(k)=min(1.0_dp-epsu,max(epsu,wa(k)))
         end do
         e(1)=normal_cdf(b(1)/r(1,1))
         ea(1)=e(1)
         f=e(1)
         fa=ea(1)
         do k=2,n
            y(k-1)=qnorm(min(1.0_dp-epsu,max(epsu,w(k-1)*e(k-1))))
            ya(k-1)=qnorm(min(1.0_dp-epsu,max(epsu,wa(k-1)*ea(k-1))))
            dummy=0.0_dp
            dummya=0.0_dp
            do l=1,k-1
               dummy=dummy+r(k,l)*y(l)
               dummya=dummya+r(k,l)*ya(l)
            end do
            e(k)=normal_cdf((b(k)-dummy)/r(k,k))
            ea(k)=normal_cdf((b(k)-dummya)/r(k,k))
            f=f*e(k)
            fa=fa*ea(k)
         end do
         prob=prob+(0.5_dp*(f+fa)-prob)/real(j+1,dp)
      end do
      prob=min(1.0_dp,max(0.0_dp,prob))
      if(present(info))info=0
   end function mvnorm_cdf_qmc

   function mvstudent_cdf_qmc(bounds,df,mean,scale,n_per_dim,info) result(prob)
      ! Randomized lattice/antithetic integration corresponding to upstream
      ! pmvt().  scale is the Student scale matrix, not its covariance.
      real(dp),intent(in)::bounds(:),df,mean(:),scale(:,:)
      integer,intent(in),optional::n_per_dim
      integer,intent(out),optional::info
      real(dp)::prob
      integer::n,nmc,j,k,l,istat
      real(dp),allocatable::b0(:),r(:,:),delta(:),q(:),w(:),wa(:),e(:),ea(:),y(:),ya(:),b(:),ba(:)
      real(dp)::f,fa,dummy,dummya,u,s,sa,epsu
      n=size(bounds)
      if(n==0)then
      prob=1.0_dp
      if(present(info))info=0
      return
      end if
      if(df<=0.0_dp .or. size(mean)/=n .or. size(scale,1)/=n .or. size(scale,2)/=n .or. n>size(primes))then
         prob=0.0_dp
         if(present(info))info=1
         return
      end if
      nmc=n*1000
      if(present(n_per_dim))nmc=max(1,n*max(1,n_per_dim))
      allocate(b0(n),r(n,n),delta(n),q(n),w(n),wa(n),e(n),ea(n),b(n),ba(n))
      if(n>1)allocate(y(n-1),ya(n-1))
      b0=bounds-mean
      call chol_upper(scale,r,istat)
      if(istat/=0)then
      prob=0.0_dp
      if(present(info))info=2
      return
      end if
      r=transpose(r)
      do k=1,n
      q(k)=sqrt(real(primes(k),dp))
      delta(k)=runif1()
      end do
      prob=0.0_dp
      epsu=32.0_dp*epsilon(1.0_dp)
      do j=0,nmc-1
         do k=1,n
            u=real(j,dp)*q(k)+delta(k)
            u=u-floor(u)
            w(k)=abs(2.0_dp*u-1.0_dp)
            wa(k)=1.0_dp-w(k)
            w(k)=min(1.0_dp-epsu,max(epsu,w(k)))
            wa(k)=min(1.0_dp-epsu,max(epsu,wa(k)))
         end do
         s=sqrt(qchisq(w(n),df)/df)
         sa=sqrt(qchisq(wa(n),df)/df)
         b=b0*s
         ba=b0*sa
         e(1)=normal_cdf(b(1)/r(1,1))
         ea(1)=normal_cdf(ba(1)/r(1,1))
         f=e(1)
         fa=ea(1)
         do k=2,n
            y(k-1)=qnorm(min(1.0_dp-epsu,max(epsu,w(k-1)*e(k-1))))
            ya(k-1)=qnorm(min(1.0_dp-epsu,max(epsu,wa(k-1)*ea(k-1))))
            dummy=0.0_dp
            dummya=0.0_dp
            do l=1,k-1
               dummy=dummy+r(k,l)*y(l)
               dummya=dummya+r(k,l)*ya(l)
            end do
            e(k)=normal_cdf((b(k)-dummy)/r(k,k))
            ea(k)=normal_cdf((ba(k)-dummya)/r(k,k))
            f=f*e(k)
            fa=fa*ea(k)
         end do
         prob=prob+(0.5_dp*(f+fa)-prob)/real(j+1,dp)
      end do
      prob=min(1.0_dp,max(0.0_dp,prob))
      if(present(info))info=0
   end function mvstudent_cdf_qmc

   subroutine sort_bounds_index(x,idx)
      real(dp),intent(inout)::x(:)
      integer,intent(inout)::idx(:)
      integer::i,j,it
      real(dp)::xt
      do i=2,size(x)
         xt=x(i)
         it=idx(i)
         j=i-1
         do while(j>=1)
            if(x(j)<=xt)exit
            x(j+1)=x(j)
            idx(j+1)=idx(j)
            j=j-1
         end do
         x(j+1)=xt
         idx(j+1)=it
      end do
   end subroutine sort_bounds_index
end module spatialextremes_mvprob
