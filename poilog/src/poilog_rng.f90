! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog_rng
   use poilog_kinds, only : dp
   use poilog_math, only : pi_dp, safe_exp
   implicit none
   private
   public :: poilog_seed, rpoilog, rbipoilog, randn, rand_poisson

contains

   subroutine poilog_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729*i + 8191*i*i, huge(1)-1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine poilog_seed

   real(dp) function randn() result(z)
      real(dp) :: u1,u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1,tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi_dp*u2)
   end function randn

   integer function rand_poisson(lambda) result(k)
      real(dp), intent(in) :: lambda
      real(dp) :: l,p,u,b,a,inv_alpha,v_r,v,us,lhs,rhs
      integer :: kk

      if (lambda <= 0.0_dp) then
         k = 0
         return
      end if
      if (lambda < 30.0_dp) then
         l = exp(-lambda)
         p = 1.0_dp
         k = 0
         do
            k = k+1
            call random_number(u)
            p = p*u
            if (p <= l) exit
         end do
         k = k-1
         return
      end if

      b = 0.931_dp + 2.53_dp*sqrt(lambda)
      a = -0.059_dp + 0.02483_dp*b
      inv_alpha = 1.1239_dp + 1.1328_dp/(b-3.4_dp)
      v_r = 0.9277_dp - 3.6224_dp/(b-2.0_dp)
      do
         call random_number(u)
         call random_number(v)
         u = u-0.5_dp
         us = 0.5_dp-abs(u)
         if (us <= 0.0_dp) cycle
         kk = floor((2.0_dp*a/us+b)*u+lambda+0.43_dp)
         if (us >= 0.07_dp .and. v <= v_r .and. kk >= 0) then
            k = kk
            return
         end if
         if (kk < 0) cycle
         if (us < 0.013_dp .and. v > us) cycle
         lhs = log(v*inv_alpha/(a/(us*us)+b))
         rhs = -lambda + real(kk,dp)*log(lambda) - log_gamma(real(kk+1,dp))
         if (lhs <= rhs) then
            k = kk
            return
         end if
      end do
   end function rand_poisson

   function rpoilog(s,mu,sig,nu,cond_s,keep0) result(x)
      integer, intent(in) :: s
      real(dp), intent(in) :: mu,sig
      real(dp), intent(in), optional :: nu
      logical, intent(in), optional :: cond_s,keep0
      integer, allocatable :: x(:)
      real(dp) :: nuv,lambda
      logical :: cond,keep
      integer :: i,k,nkeep
      integer, allocatable :: tmp(:)

      if (s < 1 .or. sig < 0.0_dp) then
         allocate(x(0)); return
      end if
      nuv = 1.0_dp; if (present(nu)) nuv = nu
      cond = .false.; if (present(cond_s)) cond = cond_s
      keep = .false.; if (present(keep0)) keep = keep0
      if (nuv < 0.0_dp) then
         allocate(x(0)); return
      end if

      if (cond .and. .not. keep) then
         allocate(x(s))
         nkeep = 0
         do while (nkeep < s)
            lambda = safe_exp(mu + sig*randn())*nuv
            k = rand_poisson(lambda)
            if (k > 0) then
               nkeep = nkeep+1
               x(nkeep) = k
            end if
         end do
      else
         allocate(tmp(s))
         nkeep = 0
         do i = 1, s
            lambda = safe_exp(mu + sig*randn())*nuv
            k = rand_poisson(lambda)
            if (keep .or. k > 0) then
               nkeep = nkeep+1
               tmp(nkeep) = k
            end if
         end do
         allocate(x(nkeep))
         if (nkeep > 0) x = tmp(:nkeep)
      end if
   end function rpoilog

   function rbipoilog(s,mu1,mu2,sig1,sig2,rho,nu1,nu2,cond_s,keep0) result(xy)
      integer, intent(in) :: s
      real(dp), intent(in) :: mu1,mu2,sig1,sig2,rho
      real(dp), intent(in), optional :: nu1,nu2
      logical, intent(in), optional :: cond_s,keep0
      integer, allocatable :: xy(:,:)
      integer, allocatable :: tmp(:,:)
      real(dp) :: n1v,n2v,z1,z2,lambda1,lambda2
      logical :: cond,keep
      integer :: i,k1,k2,nkeep

      if (s < 1 .or. sig1 < 0.0_dp .or. sig2 < 0.0_dp .or. abs(rho)>1.0_dp) then
         allocate(xy(0,2)); return
      end if
      n1v=1.0_dp; if (present(nu1)) n1v=nu1
      n2v=1.0_dp; if (present(nu2)) n2v=nu2
      cond=.false.; if (present(cond_s)) cond=cond_s
      keep=.false.; if (present(keep0)) keep=keep0
      if (n1v < 0.0_dp .or. n2v < 0.0_dp) then
         allocate(xy(0,2)); return
      end if

      if (cond .and. .not. keep) then
         allocate(xy(s,2))
         nkeep=0
         do while (nkeep<s)
            z1=randn()
            z2=rho*z1+sqrt(max(0.0_dp,1.0_dp-rho*rho))*randn()
            lambda1=safe_exp(mu1+sig1*z1)*n1v
            lambda2=safe_exp(mu2+sig2*z2)*n2v
            k1=rand_poisson(lambda1); k2=rand_poisson(lambda2)
            if (k1+k2>0) then
               nkeep=nkeep+1
               xy(nkeep,:)=[k1,k2]
            end if
         end do
      else
         allocate(tmp(s,2))
         nkeep=0
         do i=1,s
            z1=randn()
            z2=rho*z1+sqrt(max(0.0_dp,1.0_dp-rho*rho))*randn()
            k1=rand_poisson(safe_exp(mu1+sig1*z1)*n1v)
            k2=rand_poisson(safe_exp(mu2+sig2*z2)*n2v)
            if (keep .or. k1+k2>0) then
               nkeep=nkeep+1
               tmp(nkeep,:)=[k1,k2]
            end if
         end do
         allocate(xy(nkeep,2))
         if (nkeep>0) xy=tmp(:nkeep,:)
      end if
   end function rbipoilog

end module poilog_rng
