! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_difference
   use fracdiff_kinds, only : dp
   use fracdiff_fft, only : fft_inplace, next_power_of_two
   use fracdiff_status, only : fd_ok, fd_invalid_input
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
      xmean = sum(x)/real(n, dp)
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

      complex(dp), allocatable :: zx(:), zw(:)
      real(dp), allocatable :: weights(:)
      real(dp) :: xmean
      integer :: n, nfft, local_status

      n = size(x)
      if (present(status)) status = fd_ok
      if (n < 2 .or. size(dx) /= n) then
         if (present(status)) status = fd_invalid_input
         return
      end if

      nfft = next_power_of_two(2*n - 1)
      allocate(zx(nfft), zw(nfft), weights(n))
      call fractional_weights(n, d, weights, local_status)
      if (local_status /= fd_ok) then
         if (present(status)) status = local_status
         return
      end if

      xmean = sum(x)/real(n, dp)
      zx = cmplx(0.0_dp, 0.0_dp, kind=dp)
      zw = cmplx(0.0_dp, 0.0_dp, kind=dp)
      zx(1:n) = cmplx(x - xmean, 0.0_dp, kind=dp)
      zw(1:n) = cmplx(weights, 0.0_dp, kind=dp)

      call fft_inplace(zx, .false.)
      call fft_inplace(zw, .false.)
      zx = zx*zw
      call fft_inplace(zx, .true.)
      dx = real(zx(1:n), dp)
   end subroutine diffseries

end module fracdiff_difference
