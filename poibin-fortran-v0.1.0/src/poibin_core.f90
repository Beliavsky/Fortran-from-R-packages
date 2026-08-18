! SPDX-License-Identifier: GPL-2.0-only
module poibin_core
   use poibin_kinds, only : dp
   use poibin_special, only : normal_cdf, normal_pdf, poisson_cdf
   implicit none
   private

   character(len=*), parameter, public :: poibin_method_dft_cf = 'DFT-CF'
   character(len=*), parameter, public :: poibin_method_rf = 'RF'
   character(len=*), parameter, public :: poibin_method_rna = 'RNA'
   character(len=*), parameter, public :: poibin_method_na = 'NA'
   character(len=*), parameter, public :: poibin_method_pa = 'PA'

   public :: dpoibin, ppoibin, qpoibin, rpoibin
   public :: dpoibin_vec, ppoibin_vec, qpoibin_vec, rpoibin_sample
   public :: poibin_pmf_dft_cf, poibin_pmf_rf, poibin_cdf_all
   public :: poibin_mean, poibin_variance, poibin_skew_numerator
   public :: poibin_support_max, poibin_valid_parameters, poibin_seed

contains

   pure logical function poibin_valid_parameters(pp, wts) result(ok)
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)

      ok = size(pp) > 0
      if (.not. ok) return
      ok = all(pp >= 0.0_dp .and. pp <= 1.0_dp)
      if (.not. ok) return
      if (present(wts)) then
         ok = size(wts) == size(pp)
         if (.not. ok) return
         ok = all(wts >= 0)
      end if
   end function poibin_valid_parameters

   pure integer function poibin_support_max(pp, wts) result(n)
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)
      if (present(wts)) then
         n = sum(wts)
      else
         n = size(pp)
      end if
   end function poibin_support_max

   pure real(dp) function poibin_mean(pp, wts) result(mu)
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)
      integer :: j

      mu = 0.0_dp
      if (present(wts)) then
         do j = 1, size(pp)
            mu = mu + real(wts(j), dp)*pp(j)
         end do
      else
         mu = sum(pp)
      end if
   end function poibin_mean

   pure real(dp) function poibin_variance(pp, wts) result(v)
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)
      integer :: j

      v = 0.0_dp
      if (present(wts)) then
         do j = 1, size(pp)
            v = v + real(wts(j), dp)*pp(j)*(1.0_dp - pp(j))
         end do
      else
         v = sum(pp*(1.0_dp - pp))
      end if
   end function poibin_variance

   pure real(dp) function poibin_skew_numerator(pp, wts) result(gamma3)
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)
      integer :: j

      gamma3 = 0.0_dp
      if (present(wts)) then
         do j = 1, size(pp)
            gamma3 = gamma3 + real(wts(j), dp)*pp(j)*(1.0_dp - pp(j))*(1.0_dp - 2.0_dp*pp(j))
         end do
      else
         gamma3 = sum(pp*(1.0_dp - pp)*(1.0_dp - 2.0_dp*pp))
      end if
   end function poibin_skew_numerator

   subroutine poibin_pmf_dft_cf(pp, pmf, wts, status)
      real(dp), intent(in) :: pp(:)
      real(dp), allocatable, intent(out) :: pmf(:)
      integer, intent(in), optional :: wts(:)
      integer, intent(out), optional :: status

      complex(dp), allocatable :: phi(:)
      complex(dp) :: z
      real(dp) :: theta, logmag, phase, mag, total
      real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)
      integer :: n, m, j, h, k, wt
      logical :: zero_phi

      if (present(status)) status = 0
      if (.not. poibin_valid_parameters(pp, wts)) then
         allocate(pmf(0:0)); pmf = 0.0_dp
         if (present(status)) status = 1
         return
      end if

      n = poibin_support_max(pp, wts)
      m = n + 1
      allocate(pmf(0:n), phi(0:n))

      do j = 0, n
         theta = twopi*real(j, dp)/real(m, dp)
         logmag = 0.0_dp
         phase = 0.0_dp
         zero_phi = .false.
         do h = 1, size(pp)
            if (present(wts)) then
               wt = wts(h)
            else
               wt = 1
            end if
            if (wt == 0) cycle
            z = cmplx(1.0_dp - pp(h) + pp(h)*cos(theta), pp(h)*sin(theta), kind=dp)
            mag = abs(z)
            if (mag <= tiny(1.0_dp)) then
               zero_phi = .true.
               exit
            end if
            logmag = logmag + real(wt, dp)*log(mag)
            phase = phase + real(wt, dp)*atan2(aimag(z), real(z, dp))
         end do
         if (zero_phi) then
            phi(j) = cmplx(0.0_dp, 0.0_dp, kind=dp)
         else
            phi(j) = exp(logmag)*cmplx(cos(phase), sin(phase), kind=dp)
         end if
      end do

      do k = 0, n
         z = cmplx(0.0_dp, 0.0_dp, kind=dp)
         do j = 0, n
            theta = -twopi*real(j*k, dp)/real(m, dp)
            z = z + phi(j)*cmplx(cos(theta), sin(theta), kind=dp)
         end do
         pmf(k) = max(0.0_dp, real(z, dp)/real(m, dp))
      end do

      total = sum(pmf)
      if (total > 0.0_dp) pmf = pmf/total
   end subroutine poibin_pmf_dft_cf

   subroutine poibin_pmf_rf(pp, pmf, wts, status)
      real(dp), intent(in) :: pp(:)
      real(dp), allocatable, intent(out) :: pmf(:)
      integer, intent(in), optional :: wts(:)
      integer, intent(out), optional :: status

      integer :: n, j, r, k, used, wt
      real(dp) :: p

      if (present(status)) status = 0
      if (.not. poibin_valid_parameters(pp, wts)) then
         allocate(pmf(0:0)); pmf = 0.0_dp
         if (present(status)) status = 1
         return
      end if

      n = poibin_support_max(pp, wts)
      allocate(pmf(0:n))
      pmf = 0.0_dp
      pmf(0) = 1.0_dp
      used = 0

      do j = 1, size(pp)
         p = pp(j)
         if (present(wts)) then
            wt = wts(j)
         else
            wt = 1
         end if
         do r = 1, wt
            used = used + 1
            do k = used, 1, -1
               pmf(k) = (1.0_dp - p)*pmf(k) + p*pmf(k - 1)
            end do
            pmf(0) = (1.0_dp - p)*pmf(0)
         end do
      end do
   end subroutine poibin_pmf_rf

   subroutine poibin_cdf_all(pp, cdf, wts, method, status)
      real(dp), intent(in) :: pp(:)
      real(dp), allocatable, intent(out) :: cdf(:)
      integer, intent(in), optional :: wts(:)
      character(len=*), intent(in), optional :: method
      integer, intent(out), optional :: status

      real(dp), allocatable :: pmf(:)
      character(len=:), allocatable :: meth
      integer :: k, n, istat

      if (present(status)) status = 0
      if (.not. poibin_valid_parameters(pp, wts)) then
         allocate(cdf(0:0)); cdf = 0.0_dp
         if (present(status)) status = 1
         return
      end if

      if (present(method)) then
         meth = upper_ascii(trim(method))
      else
         meth = poibin_method_dft_cf
      end if

      n = poibin_support_max(pp, wts)
      allocate(cdf(0:n))

      select case (meth)
      case ('DFT-CF', 'DFTCF')
         call poibin_pmf_dft_cf(pp, pmf, wts, istat)
      case ('RF')
         call poibin_pmf_rf(pp, pmf, wts, istat)
      case default
         cdf = 0.0_dp
         if (present(status)) status = 2
         return
      end select

      cdf(0) = pmf(0)
      do k = 1, n
         cdf(k) = cdf(k - 1) + pmf(k)
      end do
      cdf = max(0.0_dp, min(1.0_dp, cdf))
      cdf(n) = 1.0_dp
   end subroutine poibin_cdf_all

   real(dp) function dpoibin(k, pp, wts) result(prob)
      integer, intent(in) :: k
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)
      real(dp), allocatable :: pmf(:)
      integer :: n, status

      if (.not. poibin_valid_parameters(pp, wts)) then
         prob = 0.0_dp
         return
      end if
      n = poibin_support_max(pp, wts)
      if (k < 0 .or. k > n) then
         prob = 0.0_dp
         return
      end if
      call poibin_pmf_dft_cf(pp, pmf, wts, status)
      prob = pmf(k)
   end function dpoibin

   subroutine dpoibin_vec(k, pp, prob, wts)
      integer, intent(in) :: k(:)
      real(dp), intent(in) :: pp(:)
      real(dp), intent(out) :: prob(:)
      integer, intent(in), optional :: wts(:)
      real(dp), allocatable :: pmf(:)
      integer :: i, n, status

      if (size(prob) /= size(k)) error stop 'dpoibin_vec: output size mismatch'
      if (.not. poibin_valid_parameters(pp, wts)) then
         prob = 0.0_dp
         return
      end if
      n = poibin_support_max(pp, wts)
      call poibin_pmf_dft_cf(pp, pmf, wts, status)
      do i = 1, size(k)
         if (k(i) < 0 .or. k(i) > n) then
            prob(i) = 0.0_dp
         else
            prob(i) = pmf(k(i))
         end if
      end do
   end subroutine dpoibin_vec

   real(dp) function ppoibin(k, pp, method, wts) result(prob)
      integer, intent(in) :: k
      real(dp), intent(in) :: pp(:)
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: wts(:)
      character(len=:), allocatable :: meth
      real(dp), allocatable :: cdf(:)
      real(dp) :: mu, var, sigma, gamma3, z
      integer :: n, status

      if (.not. poibin_valid_parameters(pp, wts)) then
         prob = 0.0_dp
         return
      end if
      n = poibin_support_max(pp, wts)
      if (k < 0) then
         prob = 0.0_dp
         return
      else if (k >= n) then
         prob = 1.0_dp
         return
      end if

      if (present(method)) then
         meth = upper_ascii(trim(method))
      else
         meth = poibin_method_dft_cf
      end if

      select case (meth)
      case ('DFT-CF', 'DFTCF', 'RF')
         call poibin_cdf_all(pp, cdf, wts, meth, status)
         prob = cdf(k)
      case ('RNA', 'NA')
         mu = poibin_mean(pp, wts)
         var = poibin_variance(pp, wts)
         if (var <= 0.0_dp) then
            if (real(k, dp) >= mu) then
               prob = 1.0_dp
            else
               prob = 0.0_dp
            end if
            return
         end if
         sigma = sqrt(var)
         z = (real(k, dp) + 0.5_dp - mu)/sigma
         prob = normal_cdf(z)
         if (meth == 'RNA') then
            gamma3 = poibin_skew_numerator(pp, wts)
            prob = prob + gamma3/(6.0_dp*sigma**3)*(1.0_dp - z*z)*normal_pdf(z)
         end if
         prob = max(0.0_dp, min(1.0_dp, prob))
      case ('PA')
         mu = poibin_mean(pp, wts)
         prob = poisson_cdf(k, mu)
      case default
         prob = 0.0_dp
      end select
   end function ppoibin

   subroutine ppoibin_vec(k, pp, prob, method, wts)
      integer, intent(in) :: k(:)
      real(dp), intent(in) :: pp(:)
      real(dp), intent(out) :: prob(:)
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: wts(:)
      integer :: i

      if (size(prob) /= size(k)) error stop 'ppoibin_vec: output size mismatch'
      do i = 1, size(k)
         prob(i) = ppoibin(k(i), pp, method, wts)
      end do
   end subroutine ppoibin_vec

   integer function qpoibin(q, pp, wts) result(kq)
      real(dp), intent(in) :: q
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)
      real(dp), allocatable :: cdf(:)
      integer :: lo, hi, mid, status

      if (.not. poibin_valid_parameters(pp, wts) .or. q < 0.0_dp .or. q > 1.0_dp) then
         kq = -1
         return
      end if
      call poibin_cdf_all(pp, cdf, wts, poibin_method_dft_cf, status)
      if (q <= cdf(0)) then
         kq = 0
         return
      end if
      lo = 0
      hi = ubound(cdf, 1)
      do while (lo + 1 < hi)
         mid = (lo + hi)/2
         if (cdf(mid) >= q) then
            hi = mid
         else
            lo = mid
         end if
      end do
      kq = hi
   end function qpoibin

   subroutine qpoibin_vec(q, pp, kq, wts)
      real(dp), intent(in) :: q(:)
      real(dp), intent(in) :: pp(:)
      integer, intent(out) :: kq(:)
      integer, intent(in), optional :: wts(:)
      real(dp), allocatable :: cdf(:)
      integer :: i, lo, hi, mid, status

      if (size(kq) /= size(q)) error stop 'qpoibin_vec: output size mismatch'
      if (.not. poibin_valid_parameters(pp, wts)) then
         kq = -1
         return
      end if
      call poibin_cdf_all(pp, cdf, wts, poibin_method_dft_cf, status)
      do i = 1, size(q)
         if (q(i) < 0.0_dp .or. q(i) > 1.0_dp) then
            kq(i) = -1
         else if (q(i) <= cdf(0)) then
            kq(i) = 0
         else
            lo = 0
            hi = ubound(cdf, 1)
            do while (lo + 1 < hi)
               mid = (lo + hi)/2
               if (cdf(mid) >= q(i)) then
                  hi = mid
               else
                  lo = mid
               end if
            end do
            kq(i) = hi
         end if
      end do
   end subroutine qpoibin_vec

   integer function rpoibin(pp, wts) result(x)
      real(dp), intent(in) :: pp(:)
      integer, intent(in), optional :: wts(:)
      real(dp) :: u
      call random_number(u)
      x = qpoibin(u, pp, wts)
   end function rpoibin

   subroutine rpoibin_sample(m, pp, x, wts)
      integer, intent(in) :: m
      real(dp), intent(in) :: pp(:)
      integer, intent(out) :: x(:)
      integer, intent(in), optional :: wts(:)
      real(dp), allocatable :: cdf(:), u(:)
      integer :: i, lo, hi, mid, status

      if (size(x) /= m) error stop 'rpoibin_sample: output size mismatch'
      if (.not. poibin_valid_parameters(pp, wts)) then
         x = -1
         return
      end if
      call poibin_cdf_all(pp, cdf, wts, poibin_method_dft_cf, status)
      allocate(u(m))
      call random_number(u)
      do i = 1, m
         if (u(i) <= cdf(0)) then
            x(i) = 0
         else
            lo = 0
            hi = ubound(cdf, 1)
            do while (lo + 1 < hi)
               mid = (lo + hi)/2
               if (cdf(mid) >= u(i)) then
                  hi = mid
               else
                  lo = mid
               end if
            end do
            x(i) = hi
         end if
      end do
   end subroutine rpoibin_sample

   subroutine poibin_seed(seed)
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer :: i, n
      integer(kind=8) :: x

      call random_seed(size=n)
      allocate(put(n))
      x = int(seed, kind=8)
      if (x == 0_8) x = 104729_8
      do i = 1, n
         x = modulo(1664525_8*x + 1013904223_8, 2147483647_8)
         put(i) = int(max(1_8, x))
      end do
      call random_seed(put=put)
   end subroutine poibin_seed

   pure function upper_ascii(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c

      out = s
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) out(i:i) = achar(c - 32)
      end do
   end function upper_ascii

end module poibin_core
