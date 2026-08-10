! SPDX-License-Identifier: GPL-3.0-only
module matrix_advanced
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_shape, matrix_err_singular
   use matrix_dense, only : eye
   use matrix_decompositions, only : singular_value_decomposition, inverse_matrix, solve_linear
   implicit none
   private
   public :: pseudoinverse, null_space, schur_complement
   public :: ldlt_factor, ldlt_solve

contains

   subroutine pseudoinverse(a, ap, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ap(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: u(:,:), s(:), vt(:,:), sinv(:,:)
      real(dp) :: eps
      integer :: i, k
      call singular_value_decomposition(a, u, s, vt, info)
      if (info /= matrix_success) then
         allocate(ap(0, 0))
         return
      end if
      k = size(s)
      allocate(sinv(k, k), source=0.0_dp)
      if (k == 0) then
         eps = 0.0_dp
      else
         eps = real(max(size(a, 1), size(a, 2)), dp) * epsilon(1.0_dp) * maxval(s)
      end if
      if (present(tol)) eps = max(0.0_dp, tol)
      do i = 1, k
         if (s(i) > eps) sinv(i, i) = 1.0_dp / s(i)
      end do
      ap = matmul(transpose(vt), matmul(sinv, transpose(u)))
      info = matrix_success
   end subroutine pseudoinverse

   subroutine null_space(a, basis, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: basis(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: u(:,:), s(:), vt(:,:)
      real(dp) :: eps
      integer :: r, n
      call singular_value_decomposition(a, u, s, vt, info)
      if (info /= matrix_success) then
         allocate(basis(0, 0))
         return
      end if
      n = size(a, 2)
      if (size(s) == 0) then
         eps = 0.0_dp
      else
         eps = real(max(size(a, 1), size(a, 2)), dp) * epsilon(1.0_dp) * maxval(s)
      end if
      if (present(tol)) eps = max(0.0_dp, tol)
      r = count(s > eps)
      if (r >= n) then
         allocate(basis(n, 0))
      else if (size(vt, 1) == n) then
         basis = transpose(vt(r + 1:n, :))
      else
         ! The compact SVD cannot expose all null vectors for a wide matrix.
         call null_space_via_projector(a, basis, info, eps)
      end if
   end subroutine null_space

   subroutine null_space_via_projector(a, basis, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: basis(:,:)
      integer, intent(out) :: info
      real(dp), intent(in) :: tol
      real(dp), allocatable :: ap(:,:), projector(:,:), q(:,:), candidate(:)
      integer :: n, i, j, count_basis
      n = size(a, 2)
      call pseudoinverse(a, ap, info, tol)
      if (info /= matrix_success) then
         allocate(basis(0, 0))
         return
      end if
      projector = eye(n) - matmul(ap, a)
      allocate(q(n, n), source=0.0_dp)
      allocate(candidate(n))
      count_basis = 0
      do j = 1, n
         candidate = projector(:, j)
         do i = 1, count_basis
            candidate = candidate - dot_product(q(:, i), candidate) * q(:, i)
         end do
         if (sqrt(dot_product(candidate, candidate)) > tol) then
            count_basis = count_basis + 1
            q(:, count_basis) = candidate / sqrt(dot_product(candidate, candidate))
         end if
      end do
      basis = q(:, :count_basis)
   end subroutine null_space_via_projector

   subroutine schur_complement(a, nleading, s, info, leading)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: nleading
      real(dp), allocatable, intent(out) :: s(:,:)
      integer, intent(out) :: info
      logical, intent(in), optional :: leading
      logical :: eliminate_leading
      integer :: n, k
      real(dp), allocatable :: x(:,:)
      if (size(a, 1) /= size(a, 2)) then
         allocate(s(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      k = nleading
      if (k < 0 .or. k > n) then
         allocate(s(0, 0))
         info = matrix_err_shape
         return
      end if
      eliminate_leading = .true.
      if (present(leading)) eliminate_leading = leading
      if (eliminate_leading) then
         call solve_linear(a(:k, :k), a(:k, k + 1:), x, info)
         if (info /= matrix_success) return
         s = a(k + 1:, k + 1:) - matmul(a(k + 1:, :k), x)
      else
         call solve_linear(a(k + 1:, k + 1:), a(k + 1:, :k), x, info)
         if (info /= matrix_success) return
         s = a(:k, :k) - matmul(a(:k, k + 1:), x)
      end if
   end subroutine schur_complement

   subroutine ldlt_factor(a, l, d, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:), d(:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      real(dp) :: eps, value
      integer :: n, i, j, k
      if (size(a, 1) /= size(a, 2)) then
         allocate(l(0, 0), d(0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      allocate(l(n, n), source=0.0_dp)
      allocate(d(n), source=0.0_dp)
      eps = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))
      if (present(tol)) eps = max(0.0_dp, tol)
      do i = 1, n
         l(i, i) = 1.0_dp
      end do
      do j = 1, n
         value = a(j, j)
         do k = 1, j - 1
            value = value - l(j, k) * l(j, k) * d(k)
         end do
         if (abs(value) <= eps) then
            info = matrix_err_singular
            return
         end if
         d(j) = value
         do i = j + 1, n
            value = a(i, j)
            do k = 1, j - 1
               value = value - l(i, k) * l(j, k) * d(k)
            end do
            l(i, j) = value / d(j)
         end do
      end do
      info = matrix_success
   end subroutine ldlt_factor

   subroutine ldlt_solve(l, d, b, x, info)
      real(dp), intent(in) :: l(:,:), d(:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      integer :: n, i, k
      if (size(l, 1) /= size(l, 2) .or. size(d) /= size(l, 1) .or. &
         size(b, 1) /= size(l, 1)) then
         allocate(x(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(l, 1)
      x = b
      do i = 2, n
         do k = 1, i - 1
            x(i, :) = x(i, :) - l(i, k) * x(k, :)
         end do
      end do
      do i = 1, n
         if (abs(d(i)) <= tiny(1.0_dp)) then
            info = matrix_err_singular
            return
         end if
         x(i, :) = x(i, :) / d(i)
      end do
      do i = n - 1, 1, -1
         do k = i + 1, n
            x(i, :) = x(i, :) - l(k, i) * x(k, :)
         end do
      end do
      info = matrix_success
   end subroutine ldlt_solve

end module matrix_advanced
