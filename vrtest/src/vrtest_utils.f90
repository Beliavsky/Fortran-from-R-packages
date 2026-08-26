! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
!
! Numerical support for the Fortran translation of vrtest 1.2.
module vrtest_utils
   use vrtest_kinds, only : dp, pi
   use r_descriptive, only : r_mean, r_variance
   use r_distributions, only : r_pnorm
   use r_sorting, only : r_average_ranks, r_quantile_type7
   use r_stability, only : r_core_log_mean_exp => r_log_mean_exp
   use r_status, only : r_core_ok => r_ok
   use r_time_series, only : core_autocorrelation => r_autocorrelation
   implicit none
   private

   public :: mean_value, variance_value, autocorrelation
   public :: normal_cdf, normal_quantile, chi_square_cdf, chi_square_quantile
   public :: sample_quantile, average_ranks, solve_linear, inverse_matrix
   public :: seed_random, random_normal_vector, wild_weights, random_permutation
   public :: log_mean_exp, simpson_integral

contains

   pure real(dp) function mean_value(x) result(ans)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         ans = 0.0_dp
      else
         ans = r_mean(x)
      end if
   end function mean_value

   pure real(dp) function variance_value(x, sample) result(ans)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: sample
      logical :: use_sample

      use_sample = .true.
      if (present(sample)) use_sample = sample
      if (size(x) < 2) then
         ans = 0.0_dp
         return
      end if
      if (use_sample) then
         ans = r_variance(x, ddof=1)
      else
         ans = r_variance(x, ddof=0)
      end if
   end function variance_value

   pure real(dp) function autocorrelation(x, lag) result(ans)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lag
      integer :: n
      integer :: core_status
      real(dp), allocatable :: values(:)

      n = size(x)
      if (lag < 0 .or. lag >= n .or. n < 2) then
         ans = 0.0_dp
         return
      end if
      call core_autocorrelation(x, values, lag_max=lag, status=core_status)
      if (core_status /= r_core_ok) then
         ans = 0.0_dp
      else
         ans = values(lag)
      end if
   end function autocorrelation

   pure real(dp) function normal_cdf(x) result(ans)
      real(dp), intent(in) :: x
      ans = r_pnorm(x)
   end function normal_cdf

   pure real(dp) function normal_quantile(p) result(x)
      ! Peter J. Acklam's inverse-normal approximation with one Newton step.
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q, r, e

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if

      if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      e = normal_cdf(x) - p
      x = x - e * sqrt(2.0_dp*pi) * exp(0.5_dp*x*x)
   end function normal_quantile

   pure real(dp) function regularized_gamma_p(a, x) result(ans)
      real(dp), intent(in) :: a, x
      integer, parameter :: max_iter = 10000
      real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: i
      real(dp) :: ap, del, summ, b, c, d, h, an

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         ans = 0.0_dp
         return
      else if (x <= tiny(1.0_dp)) then
         ans = 0.0_dp
         return
      end if

      if (x < a + 1.0_dp) then
         ap = a
         summ = 1.0_dp/a
         del = summ
         do i = 1, max_iter
            ap = ap + 1.0_dp
            del = del*x/ap
            summ = summ + del
            if (abs(del) <= abs(summ)*eps) exit
         end do
         ans = summ*exp(-x + a*log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do i = 1, max_iter
            an = -real(i,dp)*(real(i,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         ans = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
      end if
      ans = max(0.0_dp, min(1.0_dp, ans))
   end function regularized_gamma_p

   pure real(dp) function chi_square_cdf(x, df) result(ans)
      real(dp), intent(in) :: x, df
      if (x <= 0.0_dp) then
         ans = 0.0_dp
      else
         ans = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
      end if
   end function chi_square_cdf

   pure real(dp) function chi_square_quantile(p, df) result(ans)
      real(dp), intent(in) :: p, df
      integer :: iter
      real(dp) :: lo, hi, mid

      if (p <= 0.0_dp) then
         ans = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         ans = huge(1.0_dp)
         return
      end if
      lo = 0.0_dp
      hi = max(df, 1.0_dp)
      do while (chi_square_cdf(hi, df) < p)
         hi = 2.0_dp*hi
         if (hi > 1.0e12_dp) exit
      end do
      do iter = 1, 160
         mid = 0.5_dp*(lo+hi)
         if (chi_square_cdf(mid, df) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      ans = 0.5_dp*(lo+hi)
   end function chi_square_quantile

   pure real(dp) function sample_quantile(x, p) result(ans)
      ! R's default type-7 sample quantile.
      real(dp), intent(in) :: x(:), p
      real(dp) :: pp

      if (size(x) == 0) then
         ans = 0.0_dp
         return
      end if
      pp = max(0.0_dp, min(1.0_dp, p))
      ans = r_quantile_type7(x, pp)
   end function sample_quantile

   pure subroutine average_ranks(x, ranks)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: ranks(size(x))
      real(dp), allocatable :: shared_ranks(:)

      call r_average_ranks(x, shared_ranks)
      ranks = shared_ranks
   end subroutine average_ranks

   pure subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer, intent(out), optional :: info
      real(dp), allocatable :: aa(:,:), bb(:)
      real(dp) :: factor, pivot, scale
      integer :: i, j, k, n, p

      n = size(b)
      if (present(info)) info = 0
      if (size(a,1) /= n .or. size(a,2) /= n) then
         x = 0.0_dp
         if (present(info)) info = -1
         return
      end if
      allocate(aa(n,n), bb(n))
      aa = a
      bb = b
      do k = 1, n-1
         p = k - 1 + maxloc(abs(aa(k:n,k)), dim=1)
         scale = max(1.0_dp, maxval(abs(aa)))
         if (abs(aa(p,k)) <= epsilon(1.0_dp)*scale) then
            x = 0.0_dp
            if (present(info)) info = k
            return
         end if
         if (p /= k) then
            do j = 1, n
               pivot = aa(k,j)
               aa(k,j) = aa(p,j)
               aa(p,j) = pivot
            end do
            pivot = bb(k)
            bb(k) = bb(p)
            bb(p) = pivot
         end if
         do i = k+1, n
            factor = aa(i,k)/aa(k,k)
            aa(i,k:n) = aa(i,k:n) - factor*aa(k,k:n)
            bb(i) = bb(i) - factor*bb(k)
         end do
      end do
      scale = max(1.0_dp, maxval(abs(aa)))
      if (abs(aa(n,n)) <= epsilon(1.0_dp)*scale) then
         x = 0.0_dp
         if (present(info)) info = n
         return
      end if
      x(n) = bb(n)/aa(n,n)
      do i = n-1, 1, -1
         x(i) = (bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
      end do
   end subroutine solve_linear

   pure subroutine inverse_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(size(a,1),size(a,2))
      integer, intent(out), optional :: info
      real(dp), allocatable :: e(:), col(:)
      integer :: i, ierr, n

      n = size(a,1)
      if (present(info)) info = 0
      if (size(a,2) /= n) then
         ainv = 0.0_dp
         if (present(info)) info = -1
         return
      end if
      allocate(e(n),col(n))
      do i = 1, n
         e = 0.0_dp
         e(i) = 1.0_dp
         call solve_linear(a,e,col,ierr)
         if (ierr /= 0) then
            ainv = 0.0_dp
            if (present(info)) info = ierr
            return
         end if
         ainv(:,i) = col
      end do
   end subroutine inverse_matrix

   subroutine seed_random(seed)
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer :: i, n
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_random

   subroutine random_normal_vector(z)
      real(dp), intent(out) :: z(:)
      integer :: i, n
      real(dp) :: u1, u2, radius, angle
      n = size(z)
      i = 1
      do while (i <= n)
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1,tiny(1.0_dp))
         radius = sqrt(-2.0_dp*log(u1))
         angle = 2.0_dp*pi*u2
         z(i) = radius*cos(angle)
         if (i+1 <= n) z(i+1) = radius*sin(angle)
         i = i + 2
      end do
   end subroutine random_normal_vector

   subroutine wild_weights(kind, z)
      character(len=*), intent(in) :: kind
      real(dp), intent(out) :: z(:)
      real(dp), allocatable :: u(:)
      real(dp) :: p, sqrt5

      select case (trim(adjustl(kind)))
      case ('Normal','normal','NORMAL')
         call random_normal_vector(z)
      case ('Mammen','mammen','MAMMEN')
         sqrt5 = sqrt(5.0_dp)
         p = (sqrt5+1.0_dp)/(2.0_dp*sqrt5)
         allocate(u(size(z)))
         call random_number(u)
         z = -(sqrt5-1.0_dp)/2.0_dp
         where (u > p) z = (sqrt5+1.0_dp)/2.0_dp
      case ('Rademacher','rademacher','RADEMACHER')
         allocate(u(size(z)))
         call random_number(u)
         z = 1.0_dp
         where (u > 0.5_dp) z = -1.0_dp
      case default
         error stop 'wild_weights: unknown bootstrap distribution'
      end select
   end subroutine wild_weights

   subroutine random_permutation(n, p)
      integer, intent(in) :: n
      integer, intent(out) :: p(n)
      integer :: i, j, tmp
      real(dp) :: u
      p = [(i,i=1,n)]
      do i = n, 2, -1
         call random_number(u)
         j = 1 + int(u*real(i,dp))
         j = min(j,i)
         tmp = p(i)
         p(i) = p(j)
         p(j) = tmp
      end do
   end subroutine random_permutation

   pure real(dp) function log_mean_exp(x) result(ans)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         ans = -huge(1.0_dp)
      else
         ans = r_core_log_mean_exp(x)
      end if
   end function log_mean_exp

   pure real(dp) function simpson_integral(a, b, f) result(ans)
      real(dp), intent(in) :: a, b, f(:)
      integer :: i, ninterval
      real(dp) :: h, total

      ninterval = size(f)-1
      if (ninterval < 2 .or. modulo(ninterval,2) /= 0) then
         ans = 0.0_dp
         return
      end if
      total = f(1)+f(size(f))
      do i = 2, size(f)-1
         if (modulo(i,2) == 0) then
            total = total + 4.0_dp*f(i)
         else
            total = total + 2.0_dp*f(i)
         end if
      end do
      h = (b-a)/real(ninterval,dp)
      ans = h*total/3.0_dp
   end function simpson_integral

end module vrtest_utils
