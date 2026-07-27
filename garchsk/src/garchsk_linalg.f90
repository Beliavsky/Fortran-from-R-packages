! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from GARCHSK 0.1.0, Copyright (C) 2021 Kei Nakagawa.
module garchsk_linalg
   use garchsk_kinds, only : dp
   implicit none
   private
   public :: invert_matrix, symmetrize_matrix

contains

   subroutine invert_matrix(a, ainv, success)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      logical, intent(out) :: success
      real(dp), allocatable :: aug(:, :), temp(:)
      real(dp) :: pivot
      integer :: n, i, k, p

      n = size(a, 1)
      success = .false.
      allocate(ainv(n, n))
      ainv = 0.0_dp
      if (n == 0 .or. size(a, 2) /= n) return
      allocate(aug(n, 2*n), temp(2*n))
      aug(:, :n) = a
      aug(:, n+1:) = 0.0_dp
      do i = 1, n
         aug(i, n+i) = 1.0_dp
      end do

      do k = 1, n
         p = k
         do i = k + 1, n
            if (abs(aug(i, k)) > abs(aug(p, k))) p = i
         end do
         if (abs(aug(p, k)) <= 100.0_dp * epsilon(1.0_dp)) return
         if (p /= k) then
            temp = aug(k, :)
            aug(k, :) = aug(p, :)
            aug(p, :) = temp
         end if
         pivot = aug(k, k)
         aug(k, :) = aug(k, :) / pivot
         do i = 1, n
            if (i == k) cycle
            pivot = aug(i, k)
            if (abs(pivot) > 0.0_dp) aug(i, :) = aug(i, :) - pivot * aug(k, :)
         end do
      end do
      ainv = aug(:, n+1:)
      success = .true.
   end subroutine invert_matrix

   pure subroutine symmetrize_matrix(a)
      real(dp), intent(inout) :: a(:, :)
      real(dp) :: v
      integer :: i, j
      do j = 1, size(a, 2)
         do i = j + 1, size(a, 1)
            v = 0.5_dp * (a(i, j) + a(j, i))
            a(i, j) = v
            a(j, i) = v
         end do
      end do
   end subroutine symmetrize_matrix

end module garchsk_linalg
