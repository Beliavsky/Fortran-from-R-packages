! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_orthogonal
   use gogarch_kinds, only : dp
   use gogarch_linalg, only : identity_matrix, determinant_matrix
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: rd2, uprod_r, umatch, unvech, vech, angle_dimension
contains

   pure function rd2(theta) result(r)
      real(dp), intent(in) :: theta
      real(dp) :: r(2,2)
      r(1,1) = cos(theta)
      r(1,2) = -sin(theta)
      r(2,1) = sin(theta)
      r(2,2) = cos(theta)
   end function rd2

   pure integer function angle_dimension(number_of_angles) result(d)
      integer, intent(in) :: number_of_angles
      real(dp) :: candidate
      candidate = 0.5_dp+sqrt(0.25_dp+2.0_dp*real(number_of_angles,dp))
      d = nint(candidate)
      if (d*(d-1)/2 /= number_of_angles) d = -1
   end function angle_dimension

   function uprod_r(theta, ok) result(u)
      real(dp), intent(in) :: theta(:)
      logical, intent(out), optional :: ok
      integer :: d, i, j, k
      real(dp) :: u(max(1,angle_dimension(size(theta))),max(1,angle_dimension(size(theta))))
      real(dp), allocatable :: rotation(:,:)
      d = angle_dimension(size(theta))
      if (d < 1) then
         u = 0.0_dp
         if (present(ok)) ok = .false.
         return
      end if
      u = identity_matrix(d)
      allocate(rotation(d,d))
      k = 0
      do i = 1, d-1
         do j = i+1, d
            k = k+1
            rotation = identity_matrix(d)
            rotation(i,i) = cos(theta(k))
            rotation(i,j) = -sin(theta(k))
            rotation(j,i) = sin(theta(k))
            rotation(j,j) = cos(theta(k))
            u = matmul(u,rotation)
         end do
      end do
      if (present(ok)) ok = .true.
   end function uprod_r

   function umatch(from, to, ok) result(matched)
      real(dp), intent(in) :: from(:,:), to(:,:)
      logical, intent(out), optional :: ok
      real(dp) :: matched(size(from,1),size(from,2))
      logical :: used(size(from,2)), det_ok
      real(dp) :: score, best_score, sign_value, det
      integer :: i, j, best_j, weakest
      real(dp) :: alignment(size(from,2))
      if (size(from,1) /= size(from,2) .or. any(shape(to) /= shape(from))) then
         matched = 0.0_dp
         if (present(ok)) ok = .false.
         return
      end if
      used = .false.
      matched = 0.0_dp
      alignment = 0.0_dp
      do i = 1, size(from,2)
         best_score = -1.0_dp
         best_j = 0
         do j = 1, size(to,2)
            if (.not. used(j)) then
               score = abs(dot_product(from(:,i),to(:,j)))
               if (score > best_score) then
                  best_score = score
                  best_j = j
               end if
            end if
         end do
         if (best_j == 0) then
            matched = 0.0_dp
            if (present(ok)) ok = .false.
            return
         end if
         sign_value = sign(1.0_dp,dot_product(from(:,i),to(:,best_j)))
         matched(:,i) = sign_value*to(:,best_j)
         used(best_j) = .true.
         alignment(i) = abs(dot_product(from(:,i),matched(:,i)))
      end do
      det = determinant_matrix(matched,det_ok)
      if (det_ok .and. det < 0.0_dp) then
         weakest = minloc(alignment,dim=1)
         matched(:,weakest) = -matched(:,weakest)
      end if
      if (present(ok)) ok = det_ok
   end function umatch

   function unvech(v, ok) result(a)
      real(dp), intent(in) :: v(:)
      logical, intent(out), optional :: ok
      integer :: n, i, j, k
      real(dp) :: a(max(1,nint((-1.0_dp+sqrt(1.0_dp+8.0_dp*real(size(v),dp)))/2.0_dp)), &
                   max(1,nint((-1.0_dp+sqrt(1.0_dp+8.0_dp*real(size(v),dp)))/2.0_dp)))
      n = nint((-1.0_dp+sqrt(1.0_dp+8.0_dp*real(size(v),dp)))/2.0_dp)
      if (n*(n+1)/2 /= size(v)) then
         a = 0.0_dp
         if (present(ok)) ok = .false.
         return
      end if
      a = 0.0_dp
      k = 0
      do j = 1, n
         do i = j, n
            k = k+1
            a(i,j) = v(k)
            a(j,i) = v(k)
         end do
      end do
      if (present(ok)) ok = .true.
   end function unvech

   function vech(a, ok) result(v)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out), optional :: ok
      real(dp) :: v(size(a,1)*(size(a,1)+1)/2)
      integer :: i, j, k, n
      n = size(a,1)
      if (size(a,2) /= n) then
         v = 0.0_dp
         if (present(ok)) ok = .false.
         return
      end if
      k = 0
      do j = 1, n
         do i = j, n
            k = k+1
            v(k) = a(i,j)
         end do
      end do
      if (present(ok)) ok = .true.
   end function vech

end module gogarch_orthogonal
