! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_difference
   use fracdiff_kinds, only : dp
   use fracdiff_fft, only : fft_inplace, next_power_of_two
   use fracdiff_status, only : fd_ok, fd_invalid_input
   use r_descriptive, only : r_mean
   implicit none
   private

   public :: diffseries, diffseries_direct, fractional_weights

contains

   subroutine fractional_weights(n, d, weights, status)
      integer, intent(in) :: n
      real(dp), intent(in) :: d
      real(dp), intent(out) :: weights(:)
      integer, intent(out), optional :: status
      integer :: k

      if (present(status)) status = fd_ok
      if (n < 1 .or. size(weights) < n) then
         if (present(status)) status = fd_invalid_input
         return
      end if

      weights(1) = 1.0_dp
      do k = 2, n
         weights(k) = weights(k - 1)*(real(k - 2, dp) - d)/real(k - 1, dp)
      end do
   end subroutine fractional_weights

   subroutine diffseries_direct(x, d, dx, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: d
      real(dp), intent(out) :: dx(:)
      integer, intent(out), optional :: status

      real(dp), allocatable :: centered(:), weights(:)
      real(dp) :: xmean
      integer :: n, i, j, local_status

      n = size(x)
      if (present(status)) status = fd_ok
      if (n < 2 .or. size(dx) /= n) then
         if (present(status)) status = fd_invalid_input
         return
      end if

      allocate(centered(n), weights(n))
      xmean = r_mean(x)
      centered = x - xmean
      call fractional_weights(n, d, weights, local_status)
      if (local_status /= fd_ok) then
         if (present(status)) status = local_status
         return
      end if

      do i = 1, n
         dx(i) = 0.0_dp
         do j = 1, i
            dx(i) = dx(i) + weights(j)*centered(i - j + 1)
         end do
      end do
   end subroutine diffseries_direct

   subroutine diffseries(x, d, dx, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: d
      real(dp), intent(out) :: dx(:)
      integer, intent(out), optional :: status

      complex(dp), allocatable :: packed(:), product(:)
      complex(dp) :: zk, zmirror, x_spectrum, w_spectrum
      real(dp) :: xmean
      real(dp) :: weight
      integer :: n, nfft, k, mirror

      n = size(x)
      if (present(status)) status = fd_ok
      if (n < 2 .or. size(dx) /= n) then
         if (present(status)) status = fd_invalid_input
         return
      end if

      nfft = next_power_of_two(2*n - 1)
      allocate(packed(nfft), product(nfft))

      xmean = r_mean(x)
      packed = cmplx(0.0_dp, 0.0_dp, kind=dp)
      packed(1:n) = cmplx(x-xmean,0.0_dp,kind=dp)
      weight=1.0_dp
      packed(1)=cmplx(real(packed(1),dp),weight,kind=dp)
      do k=2,n
         weight=weight*(real(k-2,dp)-d)/real(k-1,dp)
         packed(k)=cmplx(real(packed(k),dp),weight,kind=dp)
      end do

      ! Two real transforms can be recovered from a single complex transform:
      ! X(k)=(Z(k)+conjg(Z(-k)))/2 and W(k)=(Z(k)-conjg(Z(-k)))/(2i).
      ! This reduces the convolution from three FFTs to two.
      call fft_inplace(packed,.false.)
      do k=1,nfft
         mirror=modulo(nfft-(k-1),nfft)+1
         zk=packed(k)
         zmirror=conjg(packed(mirror))
         x_spectrum=0.5_dp*(zk+zmirror)
         w_spectrum=cmplx(0.0_dp,-0.5_dp,kind=dp)*(zk-zmirror)
         product(k)=x_spectrum*w_spectrum
      end do
      call fft_inplace(product,.true.)
      dx=real(product(1:n),dp)
   end subroutine diffseries

end module fracdiff_difference
