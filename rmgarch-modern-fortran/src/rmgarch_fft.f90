! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_fft
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : pi, normal_pdf
   use rmgarch_types, only : grid_distribution
   implicit none
   private

   public :: make_grid_distribution, normal_grid_distribution
   public :: convolve_grid_distributions, grid_density, grid_cdf, grid_quantile
   public :: grid_moments, fft_transform

contains

   function make_grid_distribution(x, density) result(dist)
      real(dp), intent(in) :: x(:), density(:)
      type(grid_distribution) :: dist
      real(dp) :: dx, total, tolerance
      integer :: i, n

      n = size(x)
      if (n < 2 .or. size(density) /= n) then
         dist%status = 2
         return
      end if
      dx = x(2)-x(1)
      tolerance = 1.0e-8_dp*max(1.0_dp,abs(dx))
      if (dx <= 0.0_dp) then
         dist%status = 2
         return
      end if
      do i = 2, n-1
         if (abs((x(i+1)-x(i))-dx) > tolerance) then
            dist%status = 3
            return
         end if
      end do
      allocate(dist%x(n),dist%density(n),dist%cdf(n))
      dist%x = x
      dist%density = max(density,0.0_dp)
      total = sum(dist%density)*dx
      if (total <= tiny(1.0_dp)) then
         dist%density = 0.0_dp
         dist%cdf = 0.0_dp
         dist%status = 4
         return
      end if
      dist%density = dist%density/total
      call build_cdf(dist%density,dx,dist%cdf)
      dist%status = 0
   end function make_grid_distribution

   function normal_grid_distribution(mean, sigma, lower, upper, npoints) result(dist)
      real(dp), intent(in) :: mean, sigma, lower, upper
      integer, intent(in) :: npoints
      type(grid_distribution) :: dist
      real(dp), allocatable :: x(:), d(:)
      integer :: i

      if (sigma <= 0.0_dp .or. upper <= lower .or. npoints < 2) then
         dist%status = 2
         return
      end if
      allocate(x(npoints),d(npoints))
      do i = 1, npoints
         x(i) = lower+(upper-lower)*real(i-1,dp)/real(npoints-1,dp)
      end do
      d = normal_pdf((x-mean)/sigma)/sigma
      dist = make_grid_distribution(x,d)
   end function normal_grid_distribution

   function convolve_grid_distributions(left, right) result(output)
      type(grid_distribution), intent(in) :: left, right
      type(grid_distribution) :: output
      complex(dp), allocatable :: a(:), b(:)
      real(dp), allocatable :: x(:), density(:)
      real(dp) :: dx_left, dx_right, dx, tolerance
      integer :: nleft, nright, nout, nfft, i

      if (left%status /= 0 .or. right%status /= 0) then
         output%status = 2
         return
      end if
      nleft = size(left%x)
      nright = size(right%x)
      if (nleft < 2 .or. nright < 2) then
         output%status = 2
         return
      end if
      dx_left = left%x(2)-left%x(1)
      dx_right = right%x(2)-right%x(1)
      tolerance = 1.0e-8_dp*max(1.0_dp,abs(dx_left),abs(dx_right))
      if (abs(dx_left-dx_right) > tolerance) then
         output%status = 3
         return
      end if
      dx = 0.5_dp*(dx_left+dx_right)
      nout = nleft+nright-1
      nfft = next_power_of_two(nout)
      allocate(a(nfft),b(nfft),x(nout),density(nout))
      a = cmplx(0.0_dp,0.0_dp,dp)
      b = cmplx(0.0_dp,0.0_dp,dp)
      a(1:nleft) = cmplx(left%density,0.0_dp,dp)
      b(1:nright) = cmplx(right%density,0.0_dp,dp)
      call fft_transform(a,.false.)
      call fft_transform(b,.false.)
      a = a*b
      call fft_transform(a,.true.)
      density = max(real(a(1:nout),dp)*dx,0.0_dp)
      do i = 1, nout
         x(i) = left%x(1)+right%x(1)+real(i-1,dp)*dx
      end do
      output = make_grid_distribution(x,density)
   end function convolve_grid_distributions

   function grid_density(dist, x) result(value)
      type(grid_distribution), intent(in) :: dist
      real(dp), intent(in) :: x
      real(dp) :: value
      if (dist%status /= 0 .or. .not. allocated(dist%x)) then
         value = 0.0_dp
      else
         value = linear_interpolate(dist%x,dist%density,x,0.0_dp,0.0_dp)
      end if
   end function grid_density

   function grid_cdf(dist, x) result(value)
      type(grid_distribution), intent(in) :: dist
      real(dp), intent(in) :: x
      real(dp) :: value
      if (dist%status /= 0 .or. .not. allocated(dist%x)) then
         value = 0.0_dp
      else
         value = linear_interpolate(dist%x,dist%cdf,x,0.0_dp,1.0_dp)
         value = min(1.0_dp,max(0.0_dp,value))
      end if
   end function grid_cdf

   function grid_quantile(dist, probability) result(value)
      type(grid_distribution), intent(in) :: dist
      real(dp), intent(in) :: probability
      real(dp) :: value, p, fraction
      integer :: lo, hi, mid, n

      if (dist%status /= 0 .or. .not. allocated(dist%x)) then
         value = 0.0_dp
         return
      end if
      p = min(1.0_dp,max(0.0_dp,probability))
      n = size(dist%x)
      if (p <= dist%cdf(1)) then
         value = dist%x(1)
         return
      else if (p >= dist%cdf(n)) then
         value = dist%x(n)
         return
      end if
      lo = 1
      hi = n
      do while (hi-lo > 1)
         mid = (lo+hi)/2
         if (dist%cdf(mid) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      if (dist%cdf(hi) <= dist%cdf(lo)+tiny(1.0_dp)) then
         value = dist%x(hi)
      else
         fraction = (p-dist%cdf(lo))/(dist%cdf(hi)-dist%cdf(lo))
         value = dist%x(lo)+fraction*(dist%x(hi)-dist%x(lo))
      end if
   end function grid_quantile

   subroutine grid_moments(dist, mean, variance, skewness, kurtosis)
      type(grid_distribution), intent(in) :: dist
      real(dp), intent(out) :: mean, variance, skewness, kurtosis
      real(dp) :: dx, third, fourth

      mean = 0.0_dp
      variance = 0.0_dp
      skewness = 0.0_dp
      kurtosis = 0.0_dp
      if (dist%status /= 0 .or. size(dist%x) < 2) return
      dx = dist%x(2)-dist%x(1)
      mean = sum(dist%x*dist%density)*dx
      variance = sum((dist%x-mean)**2*dist%density)*dx
      if (variance > tiny(1.0_dp)) then
         third = sum((dist%x-mean)**3*dist%density)*dx
         fourth = sum((dist%x-mean)**4*dist%density)*dx
         skewness = third/variance**1.5_dp
         kurtosis = fourth/(variance*variance)
      end if
   end subroutine grid_moments

   subroutine fft_transform(values, inverse)
      complex(dp), intent(inout) :: values(:)
      logical, intent(in) :: inverse
      complex(dp) :: temp, w, wm
      real(dp) :: angle
      integer :: n, i, j, bit, len, half, k

      n = size(values)
      if (n < 1 .or. iand(n,n-1) /= 0) return
      j = 1
      do i = 1, n
         if (i < j) then
            temp = values(i)
            values(i) = values(j)
            values(j) = temp
         end if
         bit = n/2
         do while (bit >= 1 .and. j > bit)
            j = j-bit
            bit = bit/2
         end do
         j = j+bit
      end do

      len = 2
      do while (len <= n)
         half = len/2
         angle = merge(2.0_dp*pi/real(len,dp),-2.0_dp*pi/real(len,dp),inverse)
         wm = cmplx(cos(angle),sin(angle),dp)
         do i = 1, n, len
            w = cmplx(1.0_dp,0.0_dp,dp)
            do k = 0, half-1
               temp = w*values(i+k+half)
               values(i+k+half) = values(i+k)-temp
               values(i+k) = values(i+k)+temp
               w = w*wm
            end do
         end do
         len = 2*len
      end do
      if (inverse) values = values/real(n,dp)
   end subroutine fft_transform

   subroutine build_cdf(density, dx, cdf)
      real(dp), intent(in) :: density(:), dx
      real(dp), intent(out) :: cdf(size(density))
      integer :: i
      cdf(1) = density(1)*dx
      do i = 2, size(density)
         cdf(i) = cdf(i-1)+density(i)*dx
      end do
      if (cdf(size(cdf)) > tiny(1.0_dp)) cdf = cdf/cdf(size(cdf))
      cdf(size(cdf)) = 1.0_dp
   end subroutine build_cdf

   function linear_interpolate(x, y, target, left_value, right_value) result(value)
      real(dp), intent(in) :: x(:), y(:), target, left_value, right_value
      real(dp) :: value, fraction
      integer :: lo, hi, mid, n
      n = size(x)
      if (target < x(1)) then
         value = left_value
         return
      else if (target > x(n)) then
         value = right_value
         return
      end if
      lo = 1
      hi = n
      do while (hi-lo > 1)
         mid = (lo+hi)/2
         if (x(mid) <= target) then
            lo = mid
         else
            hi = mid
         end if
      end do
      if (x(hi) <= x(lo)+tiny(1.0_dp)) then
         value = y(lo)
      else
         fraction = (target-x(lo))/(x(hi)-x(lo))
         value = y(lo)+fraction*(y(hi)-y(lo))
      end if
   end function linear_interpolate

   pure integer function next_power_of_two(n) result(value)
      integer, intent(in) :: n
      value = 1
      do while (value < n)
         value = 2*value
      end do
   end function next_power_of_two

end module rmgarch_fft
