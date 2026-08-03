! MIT License. Copyright (c) 2024 fcl authors.
module fcl_xirr
   use fcl_kinds, only : dp
   use fcl_dates, only : date_type, year_frac
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: xnpv, xirr

   real(dp), parameter :: tol = 1.0e-7_dp

contains

   function xnpv(rate, values, dates, status) result(value)
      real(dp), intent(in) :: rate
      real(dp), intent(in) :: values(:)
      type(date_type), intent(in) :: dates(:)
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: i

      value = 0.0_dp
      if (size(values) /= size(dates) .or. size(values) == 0 .or. rate <= -1.0_dp) then
         if (present(status)) status = 1
         return
      end if
      if (rate == 0.0_dp) then
         value = sum(values)
      else
         do i = 1, size(values)
            value = value + values(i) / (1.0_dp + rate)**year_frac(dates(i), dates(1))
         end do
      end if
      if (present(status)) status = 0
   end function xnpv

   function xirr(values, dates, guess, status) result(root)
      real(dp), intent(in) :: values(:)
      type(date_type), intent(in) :: dates(:)
      real(dp), intent(in), optional :: guess
      integer, intent(out), optional :: status
      real(dp) :: root
      real(dp) :: x, fx, dfx, xn, lo, hi, flo, fhi, mid, fmid, step
      integer :: iter, stat

      root = 0.0_dp
      if (size(values) /= size(dates) .or. size(values) < 2 .or. &
          .not. any(values > 0.0_dp) .or. .not. any(values < 0.0_dp)) then
         if (present(status)) status = 1
         return
      end if

      x = 0.0_dp
      if (present(guess)) x = guess
      do iter = 1, 30
         fx = xnpv(x, values, dates, stat)
         if (stat /= 0 .or. .not. ieee_is_finite(fx)) exit
         dfx = (xnpv(x + tol, values, dates) - xnpv(x - tol, values, dates)) / (2.0_dp * tol)
         if (abs(dfx) <= tiny(1.0_dp)) exit
         xn = x - fx / dfx
         if (xn <= -1.0_dp .or. .not. ieee_is_finite(xn)) exit
         if (abs(xn - x) <= tol .or. abs(fx) <= tol) then
            root = xn
            if (present(status)) status = 0
            return
         end if
         x = xn
      end do

      lo = -0.9999999_dp
      flo = xnpv(lo, values, dates)
      hi = 0.01_dp
      fhi = xnpv(hi, values, dates)
      step = 0.02_dp
      do iter = 1, 80
         if (ieee_is_finite(flo) .and. ieee_is_finite(fhi)) then
            if (flo * fhi <= 0.0_dp) exit
         end if
         hi = hi + step
         step = step * 1.6_dp
         fhi = xnpv(hi, values, dates)
      end do
      if (.not. ieee_is_finite(flo) .or. .not. ieee_is_finite(fhi) .or. flo * fhi > 0.0_dp) then
         if (present(status)) status = 2
         return
      end if

      do iter = 1, 2000
         mid = lo + 0.5_dp * (hi - lo)
         fmid = xnpv(mid, values, dates)
         if (abs(fmid) <= tol .or. abs(hi - lo) <= tol) then
            root = mid
            if (present(status)) status = 0
            return
         end if
         if (flo * fmid <= 0.0_dp) then
            hi = mid
            fhi = fmid
         else
            lo = mid
            flo = fmid
         end if
      end do
      if (present(status)) status = 3
   end function xirr

end module fcl_xirr
