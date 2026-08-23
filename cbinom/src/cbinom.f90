module cbinom
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_negative_inf
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter :: eps = epsilon(1.0_dp)
   real(dp), parameter :: fd_h = 1.0e-6_dp

   public :: dcbinom, pcbinom, qcbinom, rcbinom
   public :: regularized_beta

   interface dcbinom
      module procedure dcbinom_scalar
      module procedure dcbinom_vector
   end interface

   interface pcbinom
      module procedure pcbinom_scalar
      module procedure pcbinom_vector
   end interface

   interface qcbinom
      module procedure qcbinom_scalar
      module procedure qcbinom_vector
   end interface

contains

   pure function regularized_beta(x, a, b) result(p)
      real(dp), intent(in) :: x, a, b
      real(dp) :: p
      real(dp) :: bt

      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x >= 1.0_dp) then
         p = 1.0_dp
         return
      end if

      bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
               a*log(x) + b*log1p_safe(-x))

      if (x < (a + 1.0_dp)/(a + b + 2.0_dp)) then
         p = bt*beta_cf(a, b, x)/a
      else
         p = 1.0_dp - bt*beta_cf(b, a, 1.0_dp - x)/b
      end if
      p = max(0.0_dp, min(1.0_dp, p))
   end function regularized_beta

   pure function beta_cf(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      real(dp) :: cf
      integer, parameter :: maxit = 400
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      real(dp), parameter :: tol = 8.0_dp*eps
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap

      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d

      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/ &
              ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c

         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del - 1.0_dp) <= tol) exit
      end do
      cf = h
   end function beta_cf

   pure function log1p_safe(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (abs(x) > 1.0e-4_dp) then
         y = log(1.0_dp + x)
      else
         y = x*(1.0_dp + x*(-0.5_dp + x*(1.0_dp/3.0_dp + &
             x*(-0.25_dp + 0.2_dp*x))))
      end if
   end function log1p_safe

   pure function cdf_core(x, size, prob) result(f)
      real(dp), intent(in) :: x, size, prob
      real(dp) :: f
      if (x < eps) then
         f = 0.0_dp
      else if (x > size + 1.0_dp - eps) then
         f = 1.0_dp
      else
         f = 1.0_dp - regularized_beta(prob, x, size - x + 1.0_dp)
      end if
   end function cdf_core

   pure function valid_params(nsize, prob) result(ok)
      real(dp), intent(in) :: nsize, prob
      logical :: ok
      ok = (nsize >= 0.0_dp .and. nsize < huge(1.0_dp) .and. &
            prob >= 0.0_dp .and. prob <= 1.0_dp)
   end function valid_params

   function pcbinom_scalar(q, size, prob, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: q, size, prob
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: ans, f
      logical :: lt, lp

      lt = .true.; if (present(lower_tail)) lt = lower_tail
      lp = .false.; if (present(log_p)) lp = log_p
      if (.not. valid_params(size, prob)) then
         ans = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      f = cdf_core(q, size, prob)
      if (.not. lt) f = 1.0_dp - f
      if (lp) then
         if (f <= 0.0_dp) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
         else
            ans = log(f)
         end if
      else
         ans = f
      end if
   end function pcbinom_scalar

   function pcbinom_vector(q, nsize, prob, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: q(:)
      real(dp), intent(in) :: nsize, prob
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: ans(size(q))
      integer :: i
      do i = 1, size(q)
         ans(i) = pcbinom_scalar(q(i), nsize, prob, lower_tail, log_p)
      end do
   end function pcbinom_vector

   function dcbinom_scalar(x, size, prob, log_p) result(ans)
      real(dp), intent(in) :: x, size, prob
      logical, intent(in), optional :: log_p
      real(dp) :: ans
      real(dp) :: fhi, flo, den, xx
      logical :: lp

      lp = .false.; if (present(log_p)) lp = log_p
      if (.not. valid_params(size, prob)) then
         ans = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (x < 0.0_dp .or. x > size + 1.0_dp) then
         if (lp) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
         else
            ans = 0.0_dp
         end if
         return
      end if

      if (x <= eps) then
         xx = eps*1.0000001_dp
         fhi = cdf_core(xx + fd_h, size, prob)
         flo = cdf_core(xx, size, prob)
         den = fd_h
      else if (x >= fd_h .and. x <= size + 1.0_dp - fd_h) then
         fhi = cdf_core(x + fd_h, size, prob)
         flo = cdf_core(x - fd_h, size, prob)
         den = 2.0_dp*fd_h
      else if (x <= fd_h) then
         fhi = cdf_core(x + fd_h, size, prob)
         flo = cdf_core(x, size, prob)
         den = fd_h
      else
         fhi = cdf_core(x, size, prob)
         flo = cdf_core(x - fd_h, size, prob)
         den = fd_h
      end if
      ans = max(0.0_dp, (fhi - flo)/den)
      if (lp) then
         if (ans <= 0.0_dp) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
         else
            ans = log(ans)
         end if
      end if
   end function dcbinom_scalar

   function dcbinom_vector(x, nsize, prob, log_p) result(ans)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: nsize, prob
      logical, intent(in), optional :: log_p
      real(dp) :: ans(size(x))
      integer :: i
      do i = 1, size(x)
         ans(i) = dcbinom_scalar(x(i), nsize, prob, log_p)
      end do
   end function dcbinom_vector

   function qcbinom_scalar(p, size, prob, lower_tail, log_p) result(x)
      real(dp), intent(in) :: p, size, prob
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: x
      real(dp) :: target, lo, hi, mid, fm
      logical :: lt, lp
      integer :: iter

      lt = .true.; if (present(lower_tail)) lt = lower_tail
      lp = .false.; if (present(log_p)) lp = log_p
      if (.not. valid_params(size, prob)) then
         x = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      if (lp) then
         if (p > 0.0_dp) then
            x = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         if (p < log(tiny(1.0_dp))) then
            target = 0.0_dp
         else
            target = exp(p)
         end if
      else
         if (p < 0.0_dp .or. p > 1.0_dp) then
            x = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         target = p
      end if
      if (.not. lt) target = 1.0_dp - target

      if (prob >= 1.0_dp - 2.0_dp*eps) then
         x = size + 1.0_dp
         return
      else if (prob < 2.0_dp*eps) then
         x = 0.0_dp
         return
      end if
      if (target <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (target >= 1.0_dp) then
         x = size + 1.0_dp
         return
      end if

      lo = 0.0_dp
      hi = size + 1.0_dp
      do iter = 1, 200
         mid = 0.5_dp*(lo + hi)
         fm = cdf_core(mid, size, prob)
         if (fm < target) then
            lo = mid
         else
            hi = mid
         end if
         if (hi - lo <= sqrt(eps)*max(1.0_dp, abs(mid))) exit
      end do
      x = 0.5_dp*(lo + hi)
   end function qcbinom_scalar

   function qcbinom_vector(p, nsize, prob, lower_tail, log_p) result(x)
      real(dp), intent(in) :: p(:)
      real(dp), intent(in) :: nsize, prob
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: x(size(p))
      integer :: i
      do i = 1, size(p)
         x(i) = qcbinom_scalar(p(i), nsize, prob, lower_tail, log_p)
      end do
   end function qcbinom_vector

   subroutine rcbinom(x, nsize, prob)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: nsize, prob
      real(dp) :: u(size(x))
      integer :: i
      if (.not. valid_params(nsize, prob)) then
         x = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      call random_number(u)
      do i = 1, size(x)
         x(i) = qcbinom_scalar(u(i), nsize, prob)
      end do
   end subroutine rcbinom

end module cbinom
