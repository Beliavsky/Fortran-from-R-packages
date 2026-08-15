! distr-fortran -- computational translation of the R package distr.
! Copyright (C) 2005-2025 distr authors.
! SPDX-License-Identifier: LGPL-3.0-only
module distr_rng
   use distr_kinds, only : dp, pi, nan_dp
   implicit none
   private
   public :: seed_rng, rand_uniform, rand_normal, rand_exponential, rand_gamma
   public :: rand_poisson, rand_binomial, rand_negative_binomial, rand_geometric
   public :: rand_chisq, rand_noncentral_chisq, rand_beta, rand_f, rand_t
   public :: rand_cauchy, rand_logistic, rand_weibull, rand_hypergeometric

contains

   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729*i + 7919*i*i, huge(1)-1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   real(dp) function rand_uniform() result(u)
      call random_number(u)
      do while (u <= 0.0_dp .or. u >= 1.0_dp)
         call random_number(u)
      end do
   end function rand_uniform

   real(dp) function rand_normal() result(z)
      real(dp), save :: spare = 0.0_dp
      logical, save :: has_spare = .false.
      real(dp) :: u, v, s, fac
      if (has_spare) then
         z = spare
         has_spare = .false.
         return
      end if
      do
         u = 2.0_dp*rand_uniform()-1.0_dp
         v = 2.0_dp*rand_uniform()-1.0_dp
         s = u*u+v*v
         if (s > 0.0_dp .and. s < 1.0_dp) exit
      end do
      fac = sqrt(-2.0_dp*log(s)/s)
      z = u*fac
      spare = v*fac
      has_spare = .true.
   end function rand_normal

   real(dp) function rand_exponential(rate) result(x)
      real(dp), intent(in) :: rate
      if (rate <= 0.0_dp) then
         x = nan_dp()
      else
         x = -log(rand_uniform())/rate
      end if
   end function rand_exponential

   recursive real(dp) function rand_gamma(shape, scale) result(x)
      real(dp), intent(in) :: shape, scale
      real(dp) :: d, c, z, v, u
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x = nan_dp()
         return
      end if
      if (shape < 1.0_dp) then
         x = rand_gamma(shape+1.0_dp,scale)*rand_uniform()**(1.0_dp/shape)
         return
      end if
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z = rand_normal()
            v = 1.0_dp + c*z
            if (v > 0.0_dp) exit
         end do
         v = v*v*v
         u = rand_uniform()
         if (u < 1.0_dp - 0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
      end do
      x = scale*d*v
   end function rand_gamma

   integer function rand_poisson(lambda) result(k)
      real(dp), intent(in) :: lambda
      real(dp) :: l, p, u, v, us, b, a, inv_alpha, vr, loglam, lhs, rhs, slam
      if (lambda < 0.0_dp) then
         k = -huge(1)
         return
      else if (lambda == 0.0_dp) then
         k = 0
         return
      end if
      if (lambda < 30.0_dp) then
         l = exp(-lambda)
         k = 0
         p = 1.0_dp
         do
            k = k + 1
            p = p*rand_uniform()
            if (p <= l) exit
         end do
         k = k - 1
         return
      end if
      ! PTRS transformed rejection, W. Hoermann (1993).
      slam = sqrt(lambda)
      loglam = log(lambda)
      b = 0.931_dp + 2.53_dp*slam
      a = -0.059_dp + 0.02483_dp*b
      inv_alpha = 1.1239_dp + 1.1328_dp/(b-3.4_dp)
      vr = 0.9277_dp - 3.6224_dp/(b-2.0_dp)
      do
         u = rand_uniform()-0.5_dp
         v = rand_uniform()
         us = 0.5_dp-abs(u)
         k = int(floor((2.0_dp*a/us+b)*u+lambda+0.43_dp))
         if (us >= 0.07_dp .and. v <= vr) return
         if (k < 0 .or. (us < 0.013_dp .and. v > us)) cycle
         lhs = log(v*inv_alpha/(a/(us*us)+b))
         rhs = -lambda + real(k,dp)*loglam - log_gamma(real(k+1,dp))
         if (lhs <= rhs) return
      end do
   end function rand_poisson

   integer function rand_binomial(n, prob) result(k)
      integer, intent(in) :: n
      real(dp), intent(in) :: prob
      integer :: i
      if (n < 0 .or. prob < 0.0_dp .or. prob > 1.0_dp) then
         k = -huge(1)
         return
      end if
      if (prob == 0.0_dp) then
         k = 0
      else if (prob == 1.0_dp) then
         k = n
      else
         k = 0
         do i = 1, n
            if (rand_uniform() < prob) k = k+1
         end do
      end if
   end function rand_binomial

   integer function rand_negative_binomial(size, prob) result(k)
      real(dp), intent(in) :: size, prob
      real(dp) :: lambda
      if (size <= 0.0_dp .or. prob <= 0.0_dp .or. prob > 1.0_dp) then
         k = -huge(1)
         return
      end if
      if (prob == 1.0_dp) then
         k = 0
      else
         lambda = rand_gamma(size,(1.0_dp-prob)/prob)
         k = rand_poisson(lambda)
      end if
   end function rand_negative_binomial

   integer function rand_geometric(prob) result(k)
      real(dp), intent(in) :: prob
      if (prob <= 0.0_dp .or. prob > 1.0_dp) then
         k = -huge(1)
      else if (prob == 1.0_dp) then
         k = 0
      else
         k = int(floor(log(rand_uniform())/log(1.0_dp-prob)))
      end if
   end function rand_geometric

   real(dp) function rand_chisq(df) result(x)
      real(dp), intent(in) :: df
      x = rand_gamma(0.5_dp*df,2.0_dp)
   end function rand_chisq

   real(dp) function rand_noncentral_chisq(df,ncp) result(x)
      real(dp), intent(in) :: df, ncp
      integer :: k
      if (df <= 0.0_dp .or. ncp < 0.0_dp) then
         x = nan_dp()
         return
      end if
      if (ncp == 0.0_dp) then
         x = rand_chisq(df)
      else
         k = rand_poisson(0.5_dp*ncp)
         x = rand_chisq(df+2.0_dp*real(k,dp))
      end if
   end function rand_noncentral_chisq

   real(dp) function rand_beta(a,b,ncp) result(x)
      real(dp), intent(in) :: a,b,ncp
      real(dp) :: u,v
      if (ncp == 0.0_dp) then
         u = rand_gamma(a,1.0_dp)
         v = rand_gamma(b,1.0_dp)
      else
         u = rand_noncentral_chisq(2.0_dp*a,ncp)
         v = rand_chisq(2.0_dp*b)
      end if
      x = u/(u+v)
   end function rand_beta

   real(dp) function rand_f(df1,df2,ncp) result(x)
      real(dp), intent(in) :: df1,df2,ncp
      x = (rand_noncentral_chisq(df1,ncp)/df1)/(rand_chisq(df2)/df2)
   end function rand_f

   real(dp) function rand_t(df,ncp) result(x)
      real(dp), intent(in) :: df,ncp
      x = (rand_normal()+ncp)/sqrt(rand_chisq(df)/df)
   end function rand_t

   real(dp) function rand_cauchy(location,scale) result(x)
      real(dp), intent(in) :: location,scale
      x = location + scale*tan(pi*(rand_uniform()-0.5_dp))
   end function rand_cauchy

   real(dp) function rand_logistic(location,scale) result(x)
      real(dp), intent(in) :: location,scale
      real(dp) :: u
      u = rand_uniform()
      x = location + scale*log(u/(1.0_dp-u))
   end function rand_logistic

   real(dp) function rand_weibull(shape,scale) result(x)
      real(dp), intent(in) :: shape,scale
      x = scale*(-log(rand_uniform()))**(1.0_dp/shape)
   end function rand_weibull

   integer function rand_hypergeometric(m,n,kdraw) result(x)
      integer, intent(in) :: m,n,kdraw
      integer :: good, bad, i
      real(dp) :: prob
      if (m < 0 .or. n < 0 .or. kdraw < 0 .or. kdraw > m+n) then
         x = -huge(1)
         return
      end if
      good = m
      bad = n
      x = 0
      do i = 1, kdraw
         if (good+bad <= 0) exit
         prob = real(good,dp)/real(good+bad,dp)
         if (rand_uniform() < prob) then
            x = x+1
            good = good-1
         else
            bad = bad-1
         end if
      end do
   end function rand_hypergeometric

end module distr_rng
