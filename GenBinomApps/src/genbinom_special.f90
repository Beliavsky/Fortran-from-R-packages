module genbinom_special
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use genbinom_kinds, only : dp
   implicit none
   private

   public :: beta_cdf, beta_quantile
   public :: bisect_probability_root

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
   end interface

contains

   pure elemental real(dp) function log1p_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: x2
      if (abs(x) < 1.0e-5_dp) then
         x2 = x*x
         y = x-0.5_dp*x2+x*x2/3.0_dp-x2*x2/4.0_dp+x2*x2*x/5.0_dp
      else
         y = log(1.0_dp+x)
      end if
   end function log1p_stable

   real(dp) function beta_cf(a,b,x,status) result(h)
      real(dp), intent(in) :: a,b,x
      integer, intent(out), optional :: status
      integer, parameter :: maxit = 400
      real(dp), parameter :: eps = 4.0_dp*epsilon(1.0_dp)
      real(dp), parameter :: fpmin = tiny(1.0_dp)/epsilon(1.0_dp)
      real(dp) :: qab,qap,qam,c,d,aa,del
      integer :: m,m2

      if (present(status)) status = 0
      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d

      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/ &
              ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c

         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) return
      end do
      if (present(status)) status = 1
   end function beta_cf

   real(dp) function beta_cdf(x,a,b,status) result(p)
      real(dp), intent(in) :: x,a,b
      integer, intent(out), optional :: status
      real(dp) :: logbt,bt,cf
      integer :: istat

      if (present(status)) status = 0
      if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(a) .or. &
          .not. ieee_is_finite(b) .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         p = 0.0_dp
         if (present(status)) status = 2
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

      ! The upstream countermeasure code uses a=1e-100 solely as a numerical
      ! proxy for the limiting Beta(0,b) CDF. For every x>0 that limit is 1.
      if (a < 1.0e-50_dp) then
         p = 1.0_dp
         return
      end if

      logbt = log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
              a*log(x)+b*log1p_stable(-x)
      if (logbt < log(tiny(1.0_dp))) then
         bt = 0.0_dp
      else
         bt = exp(logbt)
      end if

      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         cf = beta_cf(a,b,x,istat)
         p = bt*cf/a
      else
         cf = beta_cf(b,a,1.0_dp-x,istat)
         p = 1.0_dp-bt*cf/b
      end if

      p = min(1.0_dp,max(0.0_dp,p))
      if (present(status)) status = istat
   end function beta_cdf

   real(dp) function beta_quantile(prob,a,b,tol,max_iter,status) result(x)
      real(dp), intent(in) :: prob,a,b
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iter
      integer, intent(out), optional :: status
      real(dp) :: lo,hi,mid,fm,eps
      integer :: it,imax,istat

      if (present(status)) status = 0
      if (prob < 0.0_dp .or. prob > 1.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         x = 0.0_dp
         if (present(status)) status = 2
         return
      end if
      if (prob <= 0.0_dp) then
         x = 0.0_dp
         return
      end if
      if (prob >= 1.0_dp) then
         x = 1.0_dp
         return
      end if

      eps = 1.0e-13_dp
      imax = 200
      if (present(tol)) eps = tol
      if (present(max_iter)) imax = max_iter
      lo = 0.0_dp
      hi = 1.0_dp

      do it = 1, imax
         mid = 0.5_dp*(lo+hi)
         fm = beta_cdf(mid,a,b,istat)
         if (fm >= prob) then
            hi = mid
         else
            lo = mid
         end if
         if (hi-lo <= eps*max(1.0_dp,abs(mid))) exit
      end do
      x = 0.5_dp*(lo+hi)
      if (it > imax .and. present(status)) status = 1
   end function beta_quantile

   real(dp) function bisect_probability_root(fn,lower,upper,tol,max_iter,status) result(root)
      procedure(scalar_function) :: fn
      real(dp), intent(in) :: lower,upper
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iter
      integer, intent(out), optional :: status
      real(dp) :: lo,hi,mid,flo,fhi,fmid,eps
      integer :: it,imax

      eps = 1.0e-12_dp
      imax = 200
      if (present(tol)) eps = tol
      if (present(max_iter)) imax = max_iter
      if (present(status)) status = 0

      lo = lower
      hi = upper
      flo = fn(lo)
      fhi = fn(hi)
      if (abs(flo) <= tiny(1.0_dp)) then
         root = lo
         return
      end if
      if (abs(fhi) <= tiny(1.0_dp)) then
         root = hi
         return
      end if
      if (flo*fhi > 0.0_dp) then
         root = 0.5_dp*(lo+hi)
         if (present(status)) status = 2
         return
      end if

      do it = 1, imax
         mid = 0.5_dp*(lo+hi)
         fmid = fn(mid)
         if (abs(fmid) <= eps .or. hi-lo <= eps*max(1.0_dp,abs(mid))) exit
         if (flo*fmid <= 0.0_dp) then
            hi = mid
            fhi = fmid
         else
            lo = mid
            flo = fmid
         end if
      end do
      root = 0.5_dp*(lo+hi)
      if (it > imax .and. present(status)) status = 1
   end function bisect_probability_root

end module genbinom_special
