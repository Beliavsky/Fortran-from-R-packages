! SPDX-License-Identifier: GPL-2.0-or-later
module coneproj_shape
   use coneproj_kinds, only : dp
   use coneproj_types, only : cone_result, coneproj_success, coneproj_invalid_input
   use coneproj_linalg, only : solve_spd, matrix_rank
   use coneproj_core, only : cone_b
   implicit none
   private
   public :: make_delta, check_irreducible
   integer, parameter, public :: shape_increasing = 1
   integer, parameter, public :: shape_decreasing = 2
   integer, parameter, public :: shape_convex = 3
   integer, parameter, public :: shape_concave = 4
   integer, parameter, public :: shape_increasing_convex = 5
   integer, parameter, public :: shape_decreasing_convex = 6
   integer, parameter, public :: shape_increasing_concave = 7
   integer, parameter, public :: shape_decreasing_concave = 8

contains

   subroutine make_delta(x, shape, delta, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: shape
      real(dp), allocatable, intent(out) :: delta(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: xu(:), amat(:,:), wmat(:,:), atil(:,:), gram(:,:), rhs(:), sol(:)
      real(dp), allocatable :: pr(:,:), pg(:,:), tmpv(:)
      integer, allocatable :: first_idx(:)
      integer :: n, n1, i, j, c1, c2, c3, info, m
      real(dp), parameter :: sm = 1.0e-8_dp
      real(dp) :: nrm

      n = size(x)
      status = coneproj_success
      if (n < 2 .or. shape < 1 .or. shape > 8) then
         status = coneproj_invalid_input
         allocate(delta(0,n))
         return
      end if
      call sorted_unique(x, xu, first_idx)
      n1 = size(xu)
      if ((shape <= 2 .and. n1 < 2) .or. (shape >= 3 .and. n1 < 3)) then
         status = coneproj_invalid_input
         allocate(delta(0,n))
         return
      end if

      if (shape <= 2) then
         m = n1 - 1
      else if (shape <= 4) then
         m = n1 - 2
      else
         m = n1 - 1
      end if
      allocate(amat(m,n))
      amat = 0.0_dp
      if (shape <= 2) then
         do i = 1, n1 - 1
            c1 = first_idx(i)
            c2 = first_idx(i+1)
            amat(i,c1) = -1.0_dp
            amat(i,c2) = 1.0_dp
         end do
         if (shape == shape_decreasing) amat = -amat
      else if (shape <= 4) then
         do i = 1, n1 - 2
            c1 = first_idx(i)
            c2 = first_idx(i+1)
            c3 = first_idx(i+2)
            amat(i,c1) = xu(i+2) - xu(i+1)
            amat(i,c2) = xu(i) - xu(i+2)
            amat(i,c3) = xu(i+1) - xu(i)
         end do
         if (shape == shape_concave) amat = -amat
      else
         do i = 1, n1 - 2
            c1 = first_idx(i)
            c2 = first_idx(i+1)
            c3 = first_idx(i+2)
            amat(i,c1) = xu(i+2) - xu(i+1)
            amat(i,c2) = xu(i) - xu(i+2)
            amat(i,c3) = xu(i+1) - xu(i)
         end do
         select case (shape)
         case (shape_increasing_convex)
            c1 = first_idx(1); c2 = first_idx(2)
            amat(m,c1) = -1.0_dp; amat(m,c2) = 1.0_dp
         case (shape_decreasing_convex)
            c1 = first_idx(n1); c2 = first_idx(n1-1)
            amat(m,c1) = -1.0_dp; amat(m,c2) = 1.0_dp
         case (shape_increasing_concave)
            amat = -amat
            c1 = first_idx(n1); c2 = first_idx(n1-1)
            amat(m,c1) = 1.0_dp; amat(m,c2) = -1.0_dp
         case (shape_decreasing_concave)
            amat = -amat
            c1 = first_idx(1); c2 = first_idx(2)
            amat(m,c1) = 1.0_dp; amat(m,c2) = -1.0_dp
         end select
      end if

      if (n1 < n) then
         allocate(wmat(n,n1))
         wmat = 0.0_dp
         do i = 1, n1
            do j = 1, n
               if (abs(x(j)-xu(i)) < sm) wmat(j,i) = 1.0_dp
            end do
         end do
         atil = matmul(amat, wmat)
         gram = matmul(atil, transpose(atil))
         allocate(delta(m,n))
         do i = 1, m
            allocate(rhs(m))
            rhs = 0.0_dp
            rhs(i) = 1.0_dp
            call solve_spd(gram, rhs, sol, info)
            if (info /= 0) then
               status = coneproj_invalid_input
               return
            end if
            tmpv = matmul(transpose(atil), sol)
            delta(i,:) = matmul(wmat, tmpv)
            deallocate(rhs, sol, tmpv)
         end do
      else
         gram = matmul(amat, transpose(amat))
         allocate(delta(m,n))
         do i = 1, m
            rhs = amat(i,:)
            ! Solve (A A') z = e_i, then delta_i = A' z.
            deallocate(rhs)
            allocate(rhs(m))
            rhs = 0.0_dp
            rhs(i) = 1.0_dp
            call solve_spd(gram, rhs, sol, info)
            if (info /= 0) then
               status = coneproj_invalid_input
               return
            end if
            delta(i,:) = matmul(transpose(amat), sol)
            deallocate(rhs, sol)
         end do
      end if

      if (shape == shape_convex .or. shape == shape_concave) then
         allocate(pr(n,2), pg(2,2))
         pr(:,1) = 1.0_dp
         pr(:,2) = x
         pg = matmul(transpose(pr), pr)
         do i = 1, m
            rhs = matmul(transpose(pr), delta(i,:))
            call solve_spd(pg, rhs, sol, info)
            if (info /= 0) then
               status = coneproj_invalid_input
               return
            end if
            delta(i,:) = delta(i,:) - matmul(pr, sol)
            deallocate(rhs, sol)
         end do
      else
         do i = 1, m
            delta(i,:) = delta(i,:) - sum(delta(i,:)) / real(n,dp)
         end do
      end if
      do i = 1, m
         nrm = sqrt(dot_product(delta(i,:), delta(i,:)))
         if (nrm <= tiny(1.0_dp)) then
            status = coneproj_invalid_input
            return
         end if
         delta(i,:) = delta(i,:) / nrm
      end do
   end subroutine make_delta

   subroutine check_irreducible(edges, keep, reducible, equal_edges, status)
      real(dp), intent(in) :: edges(:,:)
      integer, allocatable, intent(out) :: keep(:), reducible(:), equal_edges(:)
      integer, intent(out) :: status
      real(dp), allocatable :: others(:,:), hd(:)
      logical, allocatable :: removed(:), equals(:)
      type(cone_result) :: ans
      integer :: m, n, i, j, k, no
      real(dp), parameter :: tol = 1.0e-8_dp

      n = size(edges,1)
      m = size(edges,2)
      status = coneproj_success
      if (n < 1 .or. m < 1) then
         status = coneproj_invalid_input
         allocate(keep(0), reducible(0), equal_edges(0))
         return
      end if
      allocate(removed(m), equals(m))
      removed = .false.; equals = .false.
      do i = 1, m
         no = m - 1
         if (no < 1) exit
         allocate(others(n,no), hd(n))
         hd = edges(:,i)
         k = 0
         do j = 1, m
            if (j /= i .and. .not. removed(j)) then
               k = k + 1
               if (k <= no) others(:,k) = edges(:,j)
            end if
         end do
         if (k == 0) then
            deallocate(others, hd)
            cycle
         end if
         if (k < no) others = others(:,1:k)
         call cone_b(hd, others, ans)
         if (ans%status == coneproj_success .and. maxval(abs(ans%fit-hd)) <= tol) then
            removed(i) = .true.
         else
            call cone_b(-hd, others, ans)
            if (ans%status == coneproj_success .and. maxval(abs(ans%fit+hd)) <= tol) equals(i) = .true.
         end if
         deallocate(others, hd)
      end do
      call mask_indices(.not. removed, keep)
      call mask_indices(removed, reducible)
      call mask_indices(equals, equal_edges)
   end subroutine check_irreducible

   subroutine sorted_unique(x, xu, first_idx)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: xu(:)
      integer, allocatable, intent(out) :: first_idx(:)
      real(dp), allocatable :: xs(:)
      integer :: i, j, n, k, imin
      real(dp) :: tmp
      n = size(x)
      allocate(xs(n))
      xs = x
      do i = 1, n-1
         imin = i
         do j = i+1, n
            if (xs(j) < xs(imin)) imin = j
         end do
         if (imin /= i) then
            tmp = xs(i); xs(i) = xs(imin); xs(imin) = tmp
         end if
      end do
      k = 1
      do i = 2, n
         if (abs(xs(i)-xs(k)) > 1.0e-12_dp * max(1.0_dp,abs(xs(k)))) then
            k = k + 1
            xs(k) = xs(i)
         end if
      end do
      allocate(xu(k), first_idx(k))
      xu = xs(1:k)
      do i = 1, k
         first_idx(i) = 1
         do j = 1, n
            if (abs(x(j)-xu(i)) <= 1.0e-8_dp) then
               first_idx(i) = j
               exit
            end if
         end do
      end do
   end subroutine sorted_unique

   subroutine mask_indices(mask, idx)
      logical, intent(in) :: mask(:)
      integer, allocatable, intent(out) :: idx(:)
      integer :: i, k
      allocate(idx(count(mask)))
      k = 0
      do i = 1, size(mask)
         if (mask(i)) then
            k = k + 1
            idx(k) = i
         end if
      end do
   end subroutine mask_indices

end module coneproj_shape
