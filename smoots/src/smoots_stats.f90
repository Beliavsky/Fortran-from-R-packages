! SPDX-License-Identifier: GPL-3.0-only
module smoots_stats
   use smoots_kinds, only : dp
   implicit none
   private
   public :: mean_value, variance_value, empirical_quantile, normal_quantile
   public :: seed_rng, random_uniform, random_normal, sample_with_replacement
   public :: factorial_real
   integer(kind=8), save :: rng_state = 88172645463393265_8
contains
   pure function mean_value(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x)/real(size(x),dp)
      end if
   end function mean_value

   pure function variance_value(x, sample) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: sample
      real(dp) :: value, mu, den
      logical :: use_sample
      use_sample = .true.
      if (present(sample)) use_sample = sample
      if (size(x) < merge(2,1,use_sample)) then
         value = 0.0_dp
         return
      end if
      mu = mean_value(x)
      den = real(size(x) - merge(1,0,use_sample), dp)
      value = sum((x-mu)**2)/den
   end function variance_value

   function empirical_quantile(x, probability) result(value)
      real(dp), intent(in) :: x(:), probability
      real(dp) :: value, h, frac
      real(dp), allocatable :: work(:)
      integer :: n, lo
      n = size(x)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      allocate(work(n)); work = x
      call quicksort(work, 1, n)
      if (probability <= 0.0_dp) then
         value = work(1)
      else if (probability >= 1.0_dp) then
         value = work(n)
      else
         h = 1.0_dp + real(n-1,dp)*probability
         lo = int(floor(h))
         frac = h-real(lo,dp)
         if (lo >= n) then
            value = work(n)
         else
            value = (1.0_dp-frac)*work(lo) + frac*work(lo+1)
         end if
      end if
   end function empirical_quantile

   recursive subroutine quicksort(a, left, right)
      real(dp), intent(inout) :: a(:)
      integer, intent(in) :: left, right
      integer :: i, j
      real(dp) :: pivot, temp
      if (left >= right) return
      i = left; j = right; pivot = a((left+right)/2)
      do
         do while (a(i) < pivot); i = i + 1; end do
         do while (a(j) > pivot); j = j - 1; end do
         if (i <= j) then
            temp = a(i); a(i) = a(j); a(j) = temp
            i = i + 1; j = j - 1
         end if
         if (i > j) exit
      end do
      if (left < j) call quicksort(a,left,j)
      if (i < right) call quicksort(a,i,right)
   end subroutine quicksort

   pure function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x, q, r
      real(dp), parameter :: a(6) = [ -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ 7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
         2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p > phigh) then
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else
         q = p-0.5_dp; r=q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      end if
   end function normal_quantile

   subroutine seed_rng(seed)
      integer(kind=8), intent(in) :: seed
      if (seed == 0_8) then
         rng_state = 88172645463393265_8
      else
         rng_state = seed
      end if
   end subroutine seed_rng

   function random_uniform() result(u)
      real(dp) :: u
      integer(kind=8) :: x
      x = rng_state
      x = ieor(x, shiftl(x,13))
      x = ieor(x, shiftr(x,7))
      x = ieor(x, shiftl(x,17))
      rng_state = x
      u = real(iand(x, int(z'001FFFFFFFFFFFFF',kind=8)),dp)/real(int(z'0020000000000000',kind=8),dp)
      if (u <= 0.0_dp) u = epsilon(1.0_dp)
   end function random_uniform

   function random_normal() result(z)
      real(dp) :: z
      real(dp), save :: spare = 0.0_dp
      logical, save :: has_spare = .false.
      real(dp) :: u1,u2
      if (has_spare) then
         z = spare; has_spare = .false.; return
      end if
      u1 = random_uniform(); u2 = random_uniform()
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
      spare = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*acos(-1.0_dp)*u2)
      has_spare = .true.
   end function random_normal

   subroutine sample_with_replacement(source, output)
      real(dp), intent(in) :: source(:)
      real(dp), intent(out) :: output(:)
      integer :: i, j, n
      n = size(source)
      if (n == 0) then
         output = 0.0_dp
         return
      end if
      do i = 1, size(output)
         j = 1 + int(random_uniform()*real(n,dp))
         j = min(n,max(1,j))
         output(i) = source(j)
      end do
   end subroutine sample_with_replacement

   pure function factorial_real(n) result(value)
      integer, intent(in) :: n
      real(dp) :: value
      integer :: i
      value = 1.0_dp
      do i = 2, n
         value = value*real(i,dp)
      end do
   end function factorial_real
end module smoots_stats
