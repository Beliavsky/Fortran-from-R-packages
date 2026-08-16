module polya_aeppli_numerics
   use polya_aeppli_kinds, only : dp
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: log1p_stable, logaddexp, log1mexp
   public :: rand_uniform, rand_normal, rand_poisson, rand_geometric_failures
   public :: set_polya_aeppli_seed

contains

   pure elemental real(dp) function log1p_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: x2
      if (abs(x) < 1.0e-5_dp) then
         x2 = x*x
         y = x - 0.5_dp*x2 + x*x2/3.0_dp - x2*x2/4.0_dp + &
             x2*x2*x/5.0_dp - x2*x2*x2/6.0_dp
      else
         y = log(1.0_dp+x)
      end if
   end function log1p_stable

   pure elemental real(dp) function logaddexp(a,b) result(y)
      real(dp), intent(in) :: a,b
      real(dp) :: m
      if (a <= -0.5_dp*huge(1.0_dp)) then
         y = b
      else if (b <= -0.5_dp*huge(1.0_dp)) then
         y = a
      else
         m = max(a,b)
         y = m + log1p_stable(exp(min(a,b)-m))
      end if
   end function logaddexp

   pure elemental real(dp) function log1mexp(x) result(y)
      real(dp), intent(in) :: x
      real(dp), parameter :: log2 = log(2.0_dp)
      if (x >= 0.0_dp) then
         y = -huge(1.0_dp)
      else if (x < -log2) then
         y = log1p_stable(-exp(x))
      else
         y = log(-expm1_series(x))
      end if
   end function log1mexp

   pure elemental real(dp) function expm1_series(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: x2
      if (abs(x) < 1.0e-5_dp) then
         x2 = x*x
         y = x + 0.5_dp*x2 + x*x2/6.0_dp + x2*x2/24.0_dp + &
             x2*x2*x/120.0_dp
      else
         y = exp(x)-1.0_dp
      end if
   end function expm1_series

   real(dp) function rand_uniform() result(u)
      call random_number(u)
      if (u <= 0.0_dp) u = tiny(1.0_dp)
      if (u >= 1.0_dp) u = 1.0_dp-spacing(1.0_dp)
   end function rand_uniform

   real(dp) function rand_normal() result(z)
      real(dp) :: u1,u2
      u1 = rand_uniform()
      u2 = rand_uniform()
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function rand_normal

   integer function rand_poisson(lambda) result(k)
      real(dp), intent(in) :: lambda
      real(dp) :: l,p,u,a,b,inv_alpha,vr,us,v,x,logv
      integer :: n
      if (lambda < 0.0_dp) error stop "rand_poisson: lambda must be nonnegative"
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
            p = p*rand_uniform()
            if (p <= l) exit
         end do
         k = k-1
         return
      end if

      ! PTRS transformed rejection, Hoermann (1993).
      b = 0.931_dp + 2.53_dp*sqrt(lambda)
      a = -0.059_dp + 0.02483_dp*b
      inv_alpha = 1.1239_dp + 1.1328_dp/(b-3.4_dp)
      vr = 0.9277_dp - 3.6224_dp/(b-2.0_dp)
      do
         u = rand_uniform()-0.5_dp
         v = rand_uniform()
         us = 0.5_dp-abs(u)
         if (us <= 0.0_dp) cycle
         x = (2.0_dp*a/us+b)*u+lambda+0.43_dp
         n = floor(x)
         if (n < 0) cycle
         if (us >= 0.07_dp .and. v <= vr) then
            k = n
            return
         end if
         if (us < 0.013_dp .and. v > us) cycle
         logv = log(v*inv_alpha/(a/(us*us)+b))
         if (logv <= -lambda + real(n,dp)*log(lambda) - log_gamma(real(n+1,dp))) then
            k = n
            return
         end if
      end do
   end function rand_poisson

   integer function rand_geometric_failures(success_prob) result(k)
      real(dp), intent(in) :: success_prob
      real(dp) :: u
      if (success_prob <= 0.0_dp .or. success_prob > 1.0_dp) then
         error stop "rand_geometric_failures: success probability must be in (0,1]"
      end if
      if (success_prob >= 1.0_dp) then
         k = 0
         return
      end if
      u = rand_uniform()
      k = floor(log(u)/log(1.0_dp-success_prob))
   end function rand_geometric_failures

   subroutine set_polya_aeppli_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: state(:)
      call random_seed(size=n)
      allocate(state(n))
      do i = 1, n
         state(i) = mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=state)
   end subroutine set_polya_aeppli_seed

end module polya_aeppli_numerics
