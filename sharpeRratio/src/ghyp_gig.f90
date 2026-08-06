! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_gig
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use ghyp_kinds, only : dp, i8
   use ghyp_special, only : log_bessel_k, digamma_fn, gauss_legendre_rule
   use ghyp_rng, only : rng_state, seed_rng, uniform_rng, gamma_rng
   implicit none
   private

   public :: gig_valid, dgig, log_dgig, pgig, qgig, rgig
   public :: gig_raw_moment, gig_mean, gig_variance, gig_mean_log, gig_mean_inverse
   public :: esgig

contains

   pure function gig_valid(lambda, chi, psi) result(ok)
      real(dp), intent(in) :: lambda, chi, psi
      logical :: ok
      ok = chi >= 0.0_dp .and. psi >= 0.0_dp
      if (chi <= tiny(1.0_dp)) ok = ok .and. lambda > 0.0_dp .and. psi > 0.0_dp
      if (psi <= tiny(1.0_dp)) ok = ok .and. lambda < 0.0_dp .and. chi > 0.0_dp
      if (chi > 0.0_dp .and. psi > 0.0_dp) ok = ok
      if (chi <= tiny(1.0_dp) .and. psi <= tiny(1.0_dp)) ok = .false.
   end function gig_valid

   function log_dgig(x, lambda, chi, psi) result(value)
      real(dp), intent(in) :: x, lambda, chi, psi
      real(dp) :: value, alpha, beta, z
      if (.not. gig_valid(lambda,chi,psi) .or. x <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      if (psi <= tiny(1.0_dp)) then
         alpha = -lambda
         beta = 0.5_dp*chi
         value = alpha*log(beta)-log_gamma(alpha)-(alpha+1.0_dp)*log(x)-beta/x
      else if (chi <= tiny(1.0_dp)) then
         alpha = lambda
         beta = 0.5_dp*psi
         value = alpha*log(beta)-log_gamma(alpha)+(alpha-1.0_dp)*log(x)-beta*x
      else
         z = sqrt(chi*psi)
         value = 0.5_dp*lambda*log(psi/chi)-log(2.0_dp)- &
            log_bessel_k(lambda,z)+(lambda-1.0_dp)*log(x)- &
            0.5_dp*(chi/x+psi*x)
      end if
   end function log_dgig

   function dgig(x, lambda, chi, psi) result(value)
      real(dp), intent(in) :: x, lambda, chi, psi
      real(dp) :: value, lv
      lv = log_dgig(x,lambda,chi,psi)
      if (lv <= log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(lv)
      end if
   end function dgig

   function gig_raw_moment(order, lambda, chi, psi) result(value)
      real(dp), intent(in) :: order, lambda, chi, psi
      real(dp) :: value, alpha, beta, z, lv
      if (.not. gig_valid(lambda,chi,psi)) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      if (psi <= tiny(1.0_dp)) then
         alpha = -lambda
         beta = 0.5_dp*chi
         if (alpha <= order) then
            value = huge(1.0_dp)
         else
            value = exp(order*log(beta)+log_gamma(alpha-order)-log_gamma(alpha))
         end if
      else if (chi <= tiny(1.0_dp)) then
         alpha = lambda
         beta = 0.5_dp*psi
         if (alpha+order <= 0.0_dp) then
            value = huge(1.0_dp)
         else
            value = exp(log_gamma(alpha+order)-log_gamma(alpha)-order*log(beta))
         end if
      else
         z = sqrt(chi*psi)
         lv = 0.5_dp*order*log(chi/psi)+log_bessel_k(lambda+order,z)- &
            log_bessel_k(lambda,z)
         if (lv > log(huge(1.0_dp))) then
            value = huge(1.0_dp)
         else
            value = exp(lv)
         end if
      end if
   end function gig_raw_moment

   function gig_mean(lambda, chi, psi) result(value)
      real(dp), intent(in) :: lambda, chi, psi
      real(dp) :: value
      value = gig_raw_moment(1.0_dp,lambda,chi,psi)
   end function gig_mean

   function gig_variance(lambda, chi, psi) result(value)
      real(dp), intent(in) :: lambda, chi, psi
      real(dp) :: value, m1, m2
      m1 = gig_raw_moment(1.0_dp,lambda,chi,psi)
      m2 = gig_raw_moment(2.0_dp,lambda,chi,psi)
      value = max(0.0_dp,m2-m1*m1)
   end function gig_variance

   function gig_mean_inverse(lambda, chi, psi) result(value)
      real(dp), intent(in) :: lambda, chi, psi
      real(dp) :: value
      value = gig_raw_moment(-1.0_dp,lambda,chi,psi)
   end function gig_mean_inverse

   function gig_mean_log(lambda, chi, psi) result(value)
      real(dp), intent(in) :: lambda, chi, psi
      real(dp) :: value, h, z
      if (.not. gig_valid(lambda,chi,psi)) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
      else if (psi <= tiny(1.0_dp)) then
         value = log(0.5_dp*chi)-digamma_fn(-lambda)
      else if (chi <= tiny(1.0_dp)) then
         value = digamma_fn(lambda)-log(0.5_dp*psi)
      else
         z = sqrt(chi*psi)
         h = max(1.0e-5_dp,abs(lambda)*1.0e-5_dp)
         value = 0.5_dp*log(chi/psi)+(log_bessel_k(lambda+h,z)- &
            log_bessel_k(lambda-h,z))/(2.0_dp*h)
      end if
   end function gig_mean_log

   function pgig(q, lambda, chi, psi) result(value)
      real(dp), intent(in) :: q, lambda, chi, psi
      real(dp) :: value, lo, hi, half, mid, y, lv
      real(dp), allocatable :: nodes(:), weights(:)
      integer :: i
      if (q <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (.not. gig_valid(lambda,chi,psi)) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      lo = min(-40.0_dp,log(q)-30.0_dp)
      hi = log(q)
      call gauss_legendre_rule(192,nodes,weights)
      half=0.5_dp*(hi-lo);mid=0.5_dp*(hi+lo);value=0.0_dp
      do i=1,size(nodes)
         y=mid+half*nodes(i)
         lv=log_dgig(exp(y),lambda,chi,psi)+y
         if(lv>log(tiny(1.0_dp)).and.lv<log(huge(1.0_dp))) &
            value=value+weights(i)*exp(lv)
      end do
      value=half*value
      value = min(1.0_dp,max(0.0_dp,value))
   end function pgig

   function qgig(p, lambda, chi, psi) result(value)
      real(dp), intent(in) :: p, lambda, chi, psi
      real(dp) :: value, lo, hi, mid, m, v
      integer :: iter
      if (p <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      m = gig_mean(lambda,chi,psi)
      v = gig_variance(lambda,chi,psi)
      if (.not. ieee_is_finite(m)) m = 1.0_dp
      if (.not. ieee_is_finite(v)) v = m*m
      lo = max(tiny(1.0_dp),m-10.0_dp*sqrt(max(v,0.0_dp)))
      hi = max(1.0_dp,m+10.0_dp*sqrt(max(v,0.0_dp)))
      do while (pgig(hi,lambda,chi,psi) < p .and. hi < huge(1.0_dp)/4.0_dp)
         hi = 2.0_dp*hi
      end do
      do iter = 1, 100
         mid = sqrt(lo*hi)
         if (pgig(mid,lambda,chi,psi) < p) then
            lo = mid
         else
            hi = mid
         end if
         if (abs(hi-lo) <= 1.0e-10_dp*max(1.0_dp,mid)) exit
      end do
      value = 0.5_dp*(lo+hi)
   end function qgig

   function rgig(lambda, chi, psi, rng, seed) result(value)
      real(dp), intent(in) :: lambda, chi, psi
      type(rng_state), intent(inout), optional :: rng
      integer(i8), intent(in), optional :: seed
      real(dp) :: value
      type(rng_state) :: local_rng
      real(dp) :: u
      if (present(seed)) call seed_rng(local_rng,seed)
      if (present(rng)) then
         if (psi <= tiny(1.0_dp)) then
            value = 1.0_dp/gamma_rng(-lambda,2.0_dp/chi,rng)
         else if (chi <= tiny(1.0_dp)) then
            value = gamma_rng(lambda,2.0_dp/psi,rng)
         else
            u = uniform_rng(rng)
            value = qgig(u,lambda,chi,psi)
         end if
      else
         if (.not. present(seed)) call seed_rng(local_rng,918273645_i8)
         if (psi <= tiny(1.0_dp)) then
            value = 1.0_dp/gamma_rng(-lambda,2.0_dp/chi,local_rng)
         else if (chi <= tiny(1.0_dp)) then
            value = gamma_rng(lambda,2.0_dp/psi,local_rng)
         else
            u = uniform_rng(local_rng)
            value = qgig(u,lambda,chi,psi)
         end if
      end if
   end function rgig

   function esgig(alpha, lambda, chi, psi, loss) result(value)
      real(dp), intent(in) :: alpha, lambda, chi, psi
      logical, intent(in), optional :: loss
      real(dp) :: value, q, lo, hi, denom, half, mid, y, lv
      real(dp), allocatable :: nodes(:), weights(:)
      logical :: upper
      integer :: i
      upper = .false.
      if (present(loss)) upper = loss
      q = qgig(alpha,lambda,chi,psi)
      if (upper) then
         lo = log(max(q,tiny(1.0_dp)))
         hi = max(lo+30.0_dp,40.0_dp)
         denom = 1.0_dp-alpha
      else
         lo = min(-40.0_dp,log(max(q,tiny(1.0_dp)))-30.0_dp)
         hi = log(max(q,tiny(1.0_dp)))
         denom = alpha
      end if
      call gauss_legendre_rule(192,nodes,weights)
      half=0.5_dp*(hi-lo);mid=0.5_dp*(hi+lo);value=0.0_dp
      do i=1,size(nodes)
         y=mid+half*nodes(i)
         lv=log_dgig(exp(y),lambda,chi,psi)+2.0_dp*y
         if(lv>log(tiny(1.0_dp)).and.lv<log(huge(1.0_dp))) &
            value=value+weights(i)*exp(lv)
      end do
      value=half*value/max(denom,tiny(1.0_dp))
   end function esgig

end module ghyp_gig
