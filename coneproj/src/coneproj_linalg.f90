! SPDX-License-Identifier: GPL-2.0-or-later
module coneproj_linalg
   use coneproj_kinds, only : dp
   implicit none
   private
   public :: solve_linear, solve_spd, cholesky_upper, solve_upper, solve_upper_transpose
   public :: inverse_upper, least_squares, matrix_rank, column_basis, inverse_matrix

contains

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:,:), bb(:), rowtmp(:)
      real(dp) :: piv, factor, scale
      integer :: n, i, k, p

      n = size(a, 1)
      info = 0
      allocate(x(n))
      x = 0.0_dp
      if (size(a, 2) /= n .or. size(b) /= n) then
         info = 1
         return
      end if
      allocate(aa(n,n), bb(n), rowtmp(n))
      aa = a
      bb = b
      scale = max(1.0_dp, maxval(abs(aa)))
      do k = 1, n
         p = k - 1 + maxloc(abs(aa(k:n,k)), dim=1)
         if (abs(aa(p,k)) <= 100.0_dp * epsilon(1.0_dp) * scale) then
            info = 2
            return
         end if
         if (p /= k) then
            rowtmp = aa(k,:)
            aa(k,:) = aa(p,:)
            aa(p,:) = rowtmp
            piv = bb(k)
            bb(k) = bb(p)
            bb(p) = piv
         end if
         piv = aa(k,k)
         do i = k + 1, n
            factor = aa(i,k) / piv
            aa(i,k:n) = aa(i,k:n) - factor * aa(k,k:n)
            bb(i) = bb(i) - factor * bb(k)
         end do
      end do
      do i = n, 1, -1
         if (i < n) then
            x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n))) / aa(i,i)
         else
            x(i) = bb(i) / aa(i,i)
         end if
      end do
   end subroutine solve_linear

   subroutine cholesky_upper(a, u, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: u(:,:)
      integer, intent(out) :: info
      real(dp) :: s
      integer :: n, i, j, k

      n = size(a,1)
      allocate(u(n,n))
      u = 0.0_dp
      info = 0
      if (size(a,2) /= n) then
         info = 1
         return
      end if
      do j = 1, n
         do i = 1, j
            s = a(i,j)
            do k = 1, i - 1
               s = s - u(k,i) * u(k,j)
            end do
            if (i == j) then
               if (s <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(a(i,i)))) then
                  info = i
                  return
               end if
               u(i,j) = sqrt(s)
            else
               u(i,j) = s / u(i,i)
            end if
         end do
      end do
   end subroutine cholesky_upper

   subroutine solve_upper(u, b, x, info)
      real(dp), intent(in) :: u(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      integer :: n, i
      n = size(u,1)
      allocate(x(n))
      x = 0.0_dp
      info = 0
      if (size(u,2) /= n .or. size(b) /= n) then
         info = 1
         return
      end if
      do i = n, 1, -1
         if (abs(u(i,i)) <= tiny(1.0_dp)) then
            info = i
            return
         end if
         if (i < n) then
            x(i) = (b(i) - dot_product(u(i,i+1:n), x(i+1:n))) / u(i,i)
         else
            x(i) = b(i) / u(i,i)
         end if
      end do
   end subroutine solve_upper

   subroutine solve_upper_transpose(u, b, x, info)
      real(dp), intent(in) :: u(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      integer :: n, i
      n = size(u,1)
      allocate(x(n))
      x = 0.0_dp
      info = 0
      if (size(u,2) /= n .or. size(b) /= n) then
         info = 1
         return
      end if
      do i = 1, n
         if (abs(u(i,i)) <= tiny(1.0_dp)) then
            info = i
            return
         end if
         if (i > 1) then
            x(i) = (b(i) - dot_product(u(1:i-1,i), x(1:i-1))) / u(i,i)
         else
            x(i) = b(i) / u(i,i)
         end if
      end do
   end subroutine solve_upper_transpose

   subroutine solve_spd(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: u(:,:), z(:)
      call cholesky_upper(a, u, info)
      if (info /= 0) then
         allocate(x(size(b)))
         x = 0.0_dp
         return
      end if
      call solve_upper_transpose(u, b, z, info)
      if (info /= 0) then
         allocate(x(size(b)))
         x = 0.0_dp
         return
      end if
      call solve_upper(u, z, x, info)
   end subroutine solve_spd

   subroutine inverse_upper(u, ui, info)
      real(dp), intent(in) :: u(:,:)
      real(dp), allocatable, intent(out) :: ui(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: e(:), col(:)
      integer :: n, j, istat
      n = size(u,1)
      allocate(ui(n,n), e(n))
      ui = 0.0_dp
      info = 0
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_upper(u, e, col, istat)
         if (istat /= 0) then
            info = istat
            return
         end if
         ui(:,j) = col
         deallocate(col)
      end do
   end subroutine inverse_upper

   subroutine inverse_matrix(a, ai, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ai(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: e(:), col(:)
      integer :: n, j, istat
      n = size(a,1)
      allocate(ai(n,n), e(n))
      ai = 0.0_dp
      info = 0
      if (size(a,2) /= n) then
         info = 1
         return
      end if
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear(a, e, col, istat)
         if (istat /= 0) then
            info = istat
            return
         end if
         ai(:,j) = col
         deallocate(col)
      end do
   end subroutine inverse_matrix

   subroutine least_squares(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: ata(:,:), atb(:)
      real(dp) :: ridge, scale
      integer :: p, k

      p = size(a,2)
      allocate(ata(p,p), atb(p))
      ata = matmul(transpose(a), a)
      atb = matmul(transpose(a), b)
      call solve_spd(ata, atb, x, info)
      if (info == 0) return
      scale = max(1.0_dp, maxval(abs(ata)))
      ridge = 100.0_dp * epsilon(1.0_dp) * scale
      do k = 1, 8
         ata = matmul(transpose(a), a)
         ata = ata + ridge * identity_matrix(p)
         call solve_spd(ata, atb, x, info)
         if (info == 0) return
         ridge = ridge * 100.0_dp
      end do
   end subroutine least_squares

   function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine column_basis(a, q, rank, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: q(:,:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: work(:,:), v(:)
      real(dp) :: nrm, eps
      integer :: n, p, j, k

      n = size(a,1)
      p = size(a,2)
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      allocate(work(n,max(1,p)), v(n))
      work = 0.0_dp
      rank = 0
      do j = 1, p
         v = a(:,j)
         do k = 1, rank
            v = v - dot_product(work(:,k), a(:,j)) * work(:,k)
         end do
         ! Reorthogonalize for stability.
         do k = 1, rank
            v = v - dot_product(work(:,k), v) * work(:,k)
         end do
         nrm = sqrt(max(0.0_dp, dot_product(v,v)))
         if (nrm > eps) then
            rank = rank + 1
            work(:,rank) = v / nrm
         end if
      end do
      allocate(q(n,rank))
      if (rank > 0) q = work(:,1:rank)
   end subroutine column_basis

   integer function matrix_rank(a, tol) result(rank)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: q(:,:)
      call column_basis(a, q, rank, tol)
   end function matrix_rank

end module coneproj_linalg
