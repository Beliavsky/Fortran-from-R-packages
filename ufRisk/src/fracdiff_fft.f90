! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_fft
   use fracdiff_kinds, only : dp, pi_dp
   implicit none
   private

   public :: fft_inplace, next_power_of_two

contains

   pure function next_power_of_two(n) result(value)
      integer, intent(in) :: n
      integer :: value

      value = 1
      do while (value < max(1, n))
         value = value*2
      end do
   end function next_power_of_two

   subroutine fft_inplace(z, inverse)
      complex(dp), intent(inout) :: z(:)
      logical, intent(in) :: inverse

      complex(dp) :: temp, w, wm
      real(dp) :: angle
      integer :: n, i, j, m, half, k

      n = size(z)
      if (n <= 1) return
      if (iand(n, n - 1) /= 0) error stop "fft_inplace: size must be a power of two"

      j = 1
      do i = 1, n
         if (i < j) then
            temp = z(i)
            z(i) = z(j)
            z(j) = temp
         end if
         m = n/2
         do while (m >= 1 .and. j > m)
            j = j - m
            m = m/2
         end do
         j = j + m
      end do

      m = 2
      do while (m <= n)
         half = m/2
         if (inverse) then
            angle = 2.0_dp*pi_dp/real(m, dp)
         else
            angle = -2.0_dp*pi_dp/real(m, dp)
         end if
         wm = cmplx(cos(angle), sin(angle), kind=dp)
         w = cmplx(1.0_dp, 0.0_dp, kind=dp)
         do j = 1, half
            do k = j, n, m
               temp = w*z(k + half)
               z(k + half) = z(k) - temp
               z(k) = z(k) + temp
            end do
            w = w*wm
         end do
         m = m*2
      end do

      if (inverse) z = z/real(n, dp)
   end subroutine fft_inplace

end module fracdiff_fft
