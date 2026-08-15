! distr-fortran -- internal radix-2 FFT utilities.
! SPDX-License-Identifier: LGPL-3.0-only
module distr_fft
   use distr_kinds, only : dp, pi
   implicit none
   private
   public :: next_power_of_two, real_convolution_fft, real_convolution_power_fft
contains

   integer function next_power_of_two(n) result(m)
      integer, intent(in) :: n
      if (n < 1) error stop 'next_power_of_two requires n >= 1'
      m = 1
      do while (m < n)
         if (m > huge(m)/2) error stop 'FFT size overflow'
         m = 2*m
      end do
   end function next_power_of_two

   subroutine fft_inplace(z, inverse)
      complex(dp), intent(inout) :: z(:)
      logical, intent(in) :: inverse
      integer :: n, i, j, m, len, half, k
      complex(dp) :: tmp, w, wlen
      real(dp) :: angle

      n = size(z)
      if (n < 1) return
      if (iand(n,n-1) /= 0) error stop 'FFT length must be a power of two'

      j = 1
      do i = 2, n
         m = n/2
         do while (j > m .and. m >= 1)
            j = j - m
            m = m/2
         end do
         j = j + m
         if (i < j) then
            tmp = z(i); z(i) = z(j); z(j) = tmp
         end if
      end do

      len = 2
      do while (len <= n)
         angle = 2.0_dp*pi/real(len,dp)
         if (.not. inverse) angle = -angle
         wlen = cmplx(cos(angle), sin(angle), kind=dp)
         half = len/2
         do i = 1, n, len
            w = cmplx(1.0_dp,0.0_dp,kind=dp)
            do k = 0, half-1
               tmp = w*z(i+k+half)
               z(i+k+half) = z(i+k)-tmp
               z(i+k) = z(i+k)+tmp
               w = w*wlen
            end do
         end do
         if (len > n/2) exit
         len = 2*len
      end do

      if (inverse) z = z/real(n,dp)
   end subroutine fft_inplace

   subroutine real_convolution_fft(a, b, c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp), allocatable, intent(out) :: c(:)
      complex(dp), allocatable :: za(:), zb(:)
      integer :: na, nb, nc, nfft

      na = size(a); nb = size(b)
      if (na < 1 .or. nb < 1) then
         allocate(c(0)); return
      end if
      nc = na + nb - 1
      nfft = next_power_of_two(nc)
      allocate(za(nfft),zb(nfft)); za = cmplx(0.0_dp,0.0_dp,kind=dp); zb = za
      za(:na) = cmplx(a,0.0_dp,kind=dp)
      zb(:nb) = cmplx(b,0.0_dp,kind=dp)
      call fft_inplace(za,.false.)
      call fft_inplace(zb,.false.)
      za = za*zb
      call fft_inplace(za,.true.)
      allocate(c(nc)); c = real(za(:nc),dp)
      where (c < 0.0_dp .and. c > -256.0_dp*epsilon(1.0_dp)) c = 0.0_dp
   end subroutine real_convolution_fft

   subroutine real_convolution_power_fft(a, exponent, c)
      real(dp), intent(in) :: a(:)
      integer, intent(in) :: exponent
      real(dp), allocatable, intent(out) :: c(:)
      complex(dp), allocatable :: z(:)
      integer :: na, nc, nfft

      if (exponent < 0) error stop 'FFT convolution power exponent must be nonnegative'
      if (exponent == 0) then
         allocate(c(1)); c = 1.0_dp; return
      end if
      na = size(a)
      if (na < 1) then
         allocate(c(0)); return
      end if
      nc = exponent*(na-1)+1
      nfft = next_power_of_two(nc)
      allocate(z(nfft)); z = cmplx(0.0_dp,0.0_dp,kind=dp)
      z(:na) = cmplx(a,0.0_dp,kind=dp)
      call fft_inplace(z,.false.)
      z = z**exponent
      call fft_inplace(z,.true.)
      allocate(c(nc)); c = real(z(:nc),dp)
      where (c < 0.0_dp .and. c > -512.0_dp*epsilon(1.0_dp)) c = 0.0_dp
   end subroutine real_convolution_power_fft

end module distr_fft
