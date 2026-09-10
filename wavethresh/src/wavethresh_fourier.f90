! SPDX-License-Identifier: GPL-2.0-or-later
! Direct real Fourier-series helpers corresponding to wavethresh rfft routines.
module wavethresh_fourier
   use wavethresh_types, only : dp
   implicit none
   private
   public :: rfft, rfftinv, rfftwt
contains

   pure function rfft(x) result(coefficients)
      real(dp), intent(in) :: x(:) !! Periodic samples on the regular grid [0,1), in wavethresh rfft convention.
      real(dp), allocatable :: coefficients(:)
      real(dp) :: angle
      real(dp) :: pi
      integer :: n
      integer :: k
      integer :: j
      integer :: position
      n = size(x)
      allocate(coefficients(n))
      if (n == 0) return
      pi = acos(-1.0_dp)
      coefficients = 0.0_dp
      coefficients(1) = sum(x) / real(n, dp)
      position = 2
      do k = 1, n / 2
         if (position > n) exit
         do j = 0, n - 1
            angle = 2.0_dp * pi * real(k * j, dp) / real(n, dp)
            coefficients(position) = coefficients(position) + sqrt(2.0_dp) * x(j + 1) * cos(angle) / real(n, dp)
         end do
         position = position + 1
         if (position > n) exit
         if (mod(n, 2) == 0 .and. k == n / 2) exit
         do j = 0, n - 1
            angle = 2.0_dp * pi * real(k * j, dp) / real(n, dp)
            coefficients(position) = coefficients(position) + sqrt(2.0_dp) * x(j + 1) * sin(angle) / real(n, dp)
         end do
         position = position + 1
      end do
   end function rfft

   pure function rfftinv(coefficients, n_points) result(x)
      real(dp), intent(in) :: coefficients(:) !! Packed real Fourier coefficients returned by rfft.
      integer, intent(in), optional :: n_points !! Number of regular output grid points; default is coefficient count.
      real(dp), allocatable :: x(:)
      real(dp) :: angle
      real(dp) :: pi
      integer :: n
      integer :: m
      integer :: k
      integer :: j
      integer :: position
      n = size(coefficients)
      if (present(n_points)) n = n_points
      if (n < size(coefficients) .or. n < 1) then
         allocate(x(0))
         return
      end if
      allocate(x(n))
      m = size(coefficients)
      pi = acos(-1.0_dp)
      x = coefficients(1)
      position = 2
      k = 1
      do while (position <= m)
         do j = 0, n - 1
            angle = 2.0_dp * pi * real(k * j, dp) / real(n, dp)
            if (n == m .and. mod(m, 2) == 0 .and. k == m / 2) then
               x(j + 1) = x(j + 1) + coefficients(position) * cos(angle) / sqrt(2.0_dp)
            else
               x(j + 1) = x(j + 1) + sqrt(2.0_dp) * coefficients(position) * cos(angle)
            end if
         end do
         position = position + 1
         if (position > m) exit
         if (n == m .and. mod(m, 2) == 0 .and. k == m / 2) exit
         do j = 0, n - 1
            angle = 2.0_dp * pi * real(k * j, dp) / real(n, dp)
            x(j + 1) = x(j + 1) + sqrt(2.0_dp) * coefficients(position) * sin(angle)
         end do
         position = position + 1
         k = k + 1
      end do
   end function rfftinv

   pure function rfftwt(coefficients, weights) result(weighted)
      real(dp), intent(in) :: coefficients(:) !! Even-length packed rfft coefficient vector.
      real(dp), intent(in) :: weights(:) !! Frequency weights, normally one per positive Fourier frequency.
      real(dp), allocatable :: weighted(:)
      integer :: position
      integer :: k
      weighted = coefficients
      if (size(weighted) == 0) return
      position = 2
      do k = 1, size(weights)
         if (position > size(weighted)) exit
         weighted(position) = weighted(position) * weights(k)
         position = position + 1
         if (position > size(weighted)) exit
         if (k == size(weights) .and. mod(size(weighted), 2) == 0) exit
         weighted(position) = weighted(position) * weights(k)
         position = position + 1
      end do
   end function rfftwt

end module wavethresh_fourier
