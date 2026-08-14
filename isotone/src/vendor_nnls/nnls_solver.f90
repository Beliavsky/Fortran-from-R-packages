! SPDX-License-Identifier: GPL-2.0-or-later
module nnls_solver
   use nnls_kinds, only : dp
   use nnls_linalg, only : least_squares_qr, norm2_stable
   implicit none
   private

   integer, parameter, public :: NNLS_SUCCESS = 1
   integer, parameter, public :: NNLS_BAD_DIMENSIONS = 2
   integer, parameter, public :: NNLS_ITERATION_LIMIT = 3

   type, public :: nnls_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: dual(:)
      integer, allocatable :: passive(:)
      integer, allocatable :: bound(:)
      real(dp) :: rnorm = 0.0_dp
      real(dp) :: deviance = 0.0_dp
      integer :: mode = NNLS_BAD_DIMENSIONS
      integer :: nsetp = 0
      integer :: iterations = 0
   end type nnls_result

   public :: nnls_fit, nnnpls_fit, nnls_solve

   interface nnls_solve
      module procedure nnls_fit
      module procedure nnnpls_fit
   end interface nnls_solve
contains
   subroutine nnls_fit(a, b, result, max_iter)
      real(dp), intent(in) :: a(:,:), b(:)
      type(nnls_result), intent(out) :: result
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: signs(:)
      integer :: n

      n = size(a,2)
      allocate(signs(n))
      signs = 1.0_dp
      call signed_nnls(a, b, signs, result, max_iter)
   end subroutine nnls_fit

   subroutine nnnpls_fit(a, b, con, result, max_iter)
      real(dp), intent(in) :: a(:,:), b(:), con(:)
      type(nnls_result), intent(out) :: result
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: signs(:)
      integer :: j, n

      n = size(a,2)
      call init_empty_result(result, size(a,1), n)
      if (size(con) /= n) then
         result%mode = NNLS_BAD_DIMENSIONS
         return
      end if
      allocate(signs(n))
      do j = 1, n
         if (con(j) < 0.0_dp) then
            signs(j) = -1.0_dp
         else
            signs(j) = 1.0_dp
         end if
      end do
      call signed_nnls(a, b, signs, result, max_iter)
   end subroutine nnnpls_fit

   subroutine signed_nnls(a, b, signs, result, max_iter)
      real(dp), intent(in) :: a(:,:), b(:), signs(:)
      type(nnls_result), intent(out) :: result
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: aw(:,:), z(:), w(:), r(:), zp(:)
      logical, allocatable :: passive(:), rejected(:)
      integer, allocatable :: order(:), pidx(:)
      integer :: m, n, j, jmax, np, iter, itmax, k, move_pos
      real(dp) :: wmax, alpha, t, wtol
      logical :: rank_ok, all_positive

      m = size(a,1)
      n = size(a,2)
      call init_empty_result(result, m, n)
      if (m <= 0 .or. n <= 0 .or. size(b) /= m .or. size(signs) /= n) then
         result%mode = NNLS_BAD_DIMENSIONS
         return
      end if
      if (present(max_iter)) then
         itmax = max_iter
      else
         itmax = 3 * n
      end if
      if (itmax < 0) itmax = 0

      allocate(aw(m,n), z(n), w(n), r(m), passive(n), rejected(n), order(n))
      do j = 1, n
         aw(:,j) = signs(j) * a(:,j)
      end do
      z = 0.0_dp
      passive = .false.
      order = [(j, j=1,n)]
      iter = 0
      wtol = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aw))) * &
         max(1.0_dp, norm2_stable(b))

      main_loop: do
         r = b - matmul(aw, z)
         w = matmul(transpose(aw), r)
         rejected = .false.

         candidate_loop: do
            wmax = 0.0_dp
            jmax = 0
            do j = 1, n
               if (.not. passive(j) .and. .not. rejected(j)) then
                  if (w(j) > wmax) then
                     wmax = w(j)
                     jmax = j
                  end if
               end if
            end do
            if (jmax == 0 .or. wmax <= wtol) exit main_loop
            if (count(passive) >= m) exit main_loop

            passive(jmax) = .true.
            call move_into_passive_order(order, passive, jmax)
            call solve_passive(aw, b, passive, order, zp, pidx, rank_ok)
            if (.not. rank_ok) then
               passive(jmax) = .false.
               rejected(jmax) = .true.
               call normalize_order(order, passive)
               cycle candidate_loop
            end if
            exit candidate_loop
         end do candidate_loop

         secondary_loop: do
            iter = iter + 1
            if (iter > itmax) then
               result%mode = NNLS_ITERATION_LIMIT
               exit main_loop
            end if
            np = size(pidx)
            all_positive = all(zp > 0.0_dp)
            if (all_positive) then
               do k = 1, np
                  z(pidx(k)) = zp(k)
               end do
               exit secondary_loop
            end if

            alpha = 2.0_dp
            move_pos = 0
            do k = 1, np
               j = pidx(k)
               if (zp(k) <= 0.0_dp) then
                  t = -z(j) / (zp(k) - z(j))
                  if (alpha > t) then
                     alpha = t
                     move_pos = k
                  end if
               end if
            end do
            if (move_pos == 0 .or. alpha > 1.0_dp) then
               do k = 1, np
                  z(pidx(k)) = max(0.0_dp, zp(k))
               end do
               exit secondary_loop
            end if
            do k = 1, np
               j = pidx(k)
               z(j) = z(j) + alpha * (zp(k) - z(j))
            end do
            do k = np, 1, -1
               j = pidx(k)
               if (z(j) <= wtol) then
                  z(j) = 0.0_dp
                  passive(j) = .false.
               end if
            end do
            call normalize_order(order, passive)
            call solve_passive(aw, b, passive, order, zp, pidx, rank_ok)
            if (.not. rank_ok) then
               result%mode = NNLS_ITERATION_LIMIT
               exit main_loop
            end if
         end do secondary_loop
      end do main_loop

      result%iterations = iter
      if (result%mode /= NNLS_ITERATION_LIMIT) result%mode = NNLS_SUCCESS
      result%x = signs * z
      result%fitted = matmul(a, result%x)
      result%residuals = b - result%fitted
      result%rnorm = norm2_stable(result%residuals)
      result%deviance = result%rnorm**2
      result%dual = signs * matmul(transpose(aw), b - matmul(aw,z))
      call fill_sets(order, passive, result)
   end subroutine signed_nnls

   subroutine solve_passive(a, b, passive, order, zp, pidx, rank_ok)
      real(dp), intent(in) :: a(:,:), b(:)
      logical, intent(in) :: passive(:)
      integer, intent(in) :: order(:)
      real(dp), allocatable, intent(out) :: zp(:)
      integer, allocatable, intent(out) :: pidx(:)
      logical, intent(out) :: rank_ok
      real(dp), allocatable :: ap(:,:)
      integer :: np, k, j

      np = count(passive)
      allocate(zp(np), pidx(np))
      if (np == 0) then
         rank_ok = .true.
         return
      end if
      allocate(ap(size(a,1),np))
      k = 0
      do j = 1, size(order)
         if (passive(order(j))) then
            k = k + 1
            pidx(k) = order(j)
            ap(:,k) = a(:,order(j))
         end if
      end do
      call least_squares_qr(ap, b, zp, rank_ok)
   end subroutine solve_passive

   subroutine move_into_passive_order(order, passive, idx)
      integer, intent(inout) :: order(:)
      logical, intent(in) :: passive(:)
      integer, intent(in) :: idx
      integer :: pos, target, tmp
      target = count(passive)
      pos = find_index(order, idx)
      if (pos > target .and. target > 0) then
         tmp = order(target)
         order(target) = order(pos)
         order(pos) = tmp
      end if
      call normalize_order(order, passive)
   end subroutine move_into_passive_order

   subroutine normalize_order(order, passive)
      integer, intent(inout) :: order(:)
      logical, intent(in) :: passive(:)
      integer, allocatable :: tmp(:)
      integer :: i, k
      allocate(tmp(size(order)))
      k = 0
      do i = 1, size(order)
         if (passive(order(i))) then
            k = k + 1
            tmp(k) = order(i)
         end if
      end do
      do i = 1, size(order)
         if (.not. passive(order(i))) then
            k = k + 1
            tmp(k) = order(i)
         end if
      end do
      order = tmp
   end subroutine normalize_order

   integer function find_index(a, value) result(pos)
      integer, intent(in) :: a(:), value
      integer :: i
      pos = 0
      do i = 1, size(a)
         if (a(i) == value) then
            pos = i
            return
         end if
      end do
   end function find_index

   subroutine fill_sets(order, passive, result)
      integer, intent(in) :: order(:)
      logical, intent(in) :: passive(:)
      type(nnls_result), intent(inout) :: result
      integer :: np, nb, i, kp, kb
      np = count(passive)
      nb = size(passive) - np
      if (allocated(result%passive)) deallocate(result%passive)
      if (allocated(result%bound)) deallocate(result%bound)
      allocate(result%passive(np), result%bound(nb))
      kp = 0
      kb = 0
      do i = 1, size(order)
         if (passive(order(i))) then
            kp = kp + 1
            result%passive(kp) = order(i)
         else
            kb = kb + 1
            result%bound(kb) = order(i)
         end if
      end do
      result%nsetp = np
   end subroutine fill_sets

   subroutine init_empty_result(result, m, n)
      type(nnls_result), intent(out) :: result
      integer, intent(in) :: m, n
      integer :: i
      allocate(result%x(max(0,n)), result%fitted(max(0,m)), &
         result%residuals(max(0,m)), result%dual(max(0,n)))
      allocate(result%passive(0), result%bound(max(0,n)))
      result%x = 0.0_dp
      result%fitted = 0.0_dp
      result%residuals = 0.0_dp
      result%dual = 0.0_dp
      if (n > 0) result%bound = [(i, i=1,n)]
      result%rnorm = 0.0_dp
      result%deviance = 0.0_dp
      result%mode = NNLS_BAD_DIMENSIONS
      result%nsetp = 0
      result%iterations = 0
   end subroutine init_empty_result
end module nnls_solver
