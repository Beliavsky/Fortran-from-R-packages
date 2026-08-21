! SPDX-License-Identifier: GPL-3.0-only
module matrix_functions
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_shape, matrix_err_invalid, &
      matrix_err_singular, matrix_err_convergence
   use matrix_dense, only : eye, one_norm, frobenius_norm, symmpart
   use matrix_decompositions, only : inverse_matrix, symmetric_eigen, solve_linear
   implicit none
   private
   public :: matrix_exponential, matrix_power, matrix_sqrt_sym
   public :: near_positive_definite, condition_number_1, reciprocal_condition_1
   public :: projection_psd

contains

   subroutine matrix_exponential(a, expa, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: expa(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: scaled(:,:), term(:,:), id(:,:)
      real(dp) :: nrm, scale
      integer :: n, s, k
      if (size(a, 1) /= size(a, 2)) then
         allocate(expa(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      if (n == 0) then
         allocate(expa(0, 0))
         info = matrix_success
         return
      end if
      nrm = one_norm(a)
      if (nrm <= 0.5_dp) then
         s = 0
      else
         s = max(0, ceiling(log(nrm / 0.5_dp) / log(2.0_dp)))
      end if
      scale = 2.0_dp ** real(s, dp)
      scaled = a / scale
      id = eye(n)
      expa = id
      term = id
      do k = 1, 80
         term = matmul(term, scaled) / real(k, dp)
         expa = expa + term
         if (frobenius_norm(term) <= 10.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, frobenius_norm(expa))) exit
      end do
      if (k > 80) then
         info = matrix_err_convergence
         return
      end if
      do k = 1, s
         expa = matmul(expa, expa)
      end do
      info = matrix_success
   end subroutine matrix_exponential

   subroutine matrix_power(a, p, ap, info)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: p
      real(dp), allocatable, intent(out) :: ap(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: base(:,:)
      integer :: n, exponent
      if (size(a, 1) /= size(a, 2)) then
         allocate(ap(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      ap = eye(n)
      if (p == 0) then
         info = matrix_success
         return
      end if
      if (p < 0) then
         call inverse_matrix(a, base, info)
         if (info /= matrix_success) return
         exponent = -p
      else
         base = a
         exponent = p
      end if
      do while (exponent > 0)
         if (mod(exponent, 2) == 1) ap = matmul(ap, base)
         exponent = exponent / 2
         if (exponent > 0) base = matmul(base, base)
      end do
      info = matrix_success
   end subroutine matrix_power

   subroutine matrix_sqrt_sym(a, root, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: root(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: values(:), vectors(:,:), d(:,:)
      real(dp) :: eps
      integer :: i, n
      if (size(a, 1) /= size(a, 2)) then
         allocate(root(0, 0))
         info = matrix_err_shape
         return
      end if
      call symmetric_eigen(a, values, vectors, info)
      if (info /= matrix_success) then
         allocate(root(0, 0))
         return
      end if
      eps = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))
      if (present(tol)) eps = max(0.0_dp, tol)
      if (any(values < -eps)) then
         allocate(root(0, 0))
         info = matrix_err_invalid
         return
      end if
      n = size(values)
      allocate(d(n, n), source=0.0_dp)
      do i = 1, n
         d(i, i) = sqrt(max(0.0_dp, values(i)))
      end do
      root = matmul(vectors, matmul(d, transpose(vectors)))
      info = matrix_success
   end subroutine matrix_sqrt_sym

   subroutine projection_psd(a, out, info, min_eigenvalue)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: out(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: min_eigenvalue
      real(dp), allocatable :: values(:), vectors(:,:), d(:,:)
      real(dp) :: floor_value
      integer :: i, n
      if (size(a, 1) /= size(a, 2)) then
         allocate(out(0, 0))
         info = matrix_err_shape
         return
      end if
      call symmetric_eigen(0.5_dp * (a + transpose(a)), values, vectors, info)
      if (info /= matrix_success) then
         allocate(out(0, 0))
         return
      end if
      floor_value = 0.0_dp
      if (present(min_eigenvalue)) floor_value = min_eigenvalue
      n = size(values)
      allocate(d(n, n), source=0.0_dp)
      do i = 1, n
         d(i, i) = max(values(i), floor_value)
      end do
      out = matmul(vectors, matmul(d, transpose(vectors)))
      out = 0.5_dp * (out + transpose(out))
      info = matrix_success
   end subroutine projection_psd

   subroutine near_positive_definite(a, out, info, keep_diagonal, corr, tol, max_iter)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: out(:,:)
      integer, intent(out) :: info
      logical, intent(in), optional :: keep_diagonal, corr
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: y(:,:), x(:,:), r(:,:), delta_s(:,:), previous(:,:), original_diag(:)
      real(dp) :: eps, rel_change, eigen_floor
      integer :: n, i, iter, maxit, pinfo
      logical :: preserve_diag, correlation
      if (size(a, 1) /= size(a, 2)) then
         allocate(out(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      preserve_diag = .false.
      correlation = .false.
      if (present(keep_diagonal)) preserve_diag = keep_diagonal
      if (present(corr)) correlation = corr
      eps = 1.0e-8_dp
      if (present(tol)) eps = max(tol, epsilon(1.0_dp))
      maxit = 100
      if (present(max_iter)) maxit = max_iter
      allocate(original_diag(n))
      do i = 1, n
         original_diag(i) = a(i, i)
      end do
      y = 0.5_dp * (a + transpose(a))
      allocate(delta_s(n, n), source=0.0_dp)
      info = matrix_err_convergence
      do iter = 1, maxit
         previous = y
         r = y - delta_s
         call projection_psd(r, x, pinfo)
         if (pinfo /= matrix_success) return
         delta_s = x - r
         y = x
         if (correlation) then
            do i = 1, n
               y(i, i) = 1.0_dp
            end do
         else if (preserve_diag) then
            do i = 1, n
               y(i, i) = original_diag(i)
            end do
         end if
         rel_change = frobenius_norm(y - previous) / max(1.0_dp, frobenius_norm(y))
         if (rel_change <= eps) then
            info = matrix_success
            exit
         end if
      end do
      eigen_floor = eps * max(1.0_dp, maxval(abs(y)))
      call projection_psd(y, out, pinfo, min_eigenvalue=eigen_floor)
      if (pinfo /= matrix_success) then
         info = pinfo
         return
      end if
      if (correlation) then
         call normalize_correlation(out)
      else if (preserve_diag) then
         do i = 1, n
            out(i, i) = original_diag(i)
         end do
      end if
   end subroutine near_positive_definite

   subroutine normalize_correlation(a)
      real(dp), intent(inout) :: a(:,:)
      real(dp), allocatable :: d(:)
      integer :: i, j, n
      n = size(a, 1)
      allocate(d(n))
      do i = 1, n
         d(i) = sqrt(max(a(i, i), tiny(1.0_dp)))
      end do
      do j = 1, n
         do i = 1, n
            a(i, j) = a(i, j) / (d(i) * d(j))
         end do
      end do
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end subroutine normalize_correlation

   subroutine condition_number_1(a, cond, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: cond
      integer, intent(out) :: info
      real(dp), allocatable :: ainv(:,:)
      if (size(a, 1) /= size(a, 2)) then
         cond = huge(1.0_dp)
         info = matrix_err_shape
         return
      end if
      call inverse_matrix(a, ainv, info)
      if (info /= matrix_success) then
         cond = huge(1.0_dp)
         return
      end if
      cond = one_norm(a) * one_norm(ainv)
   end subroutine condition_number_1

   subroutine reciprocal_condition_1(a, rcond, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: rcond
      integer, intent(out) :: info
      real(dp) :: cond
      call condition_number_1(a, cond, info)
      if (info /= matrix_success .or. cond <= 0.0_dp .or. cond >= huge(1.0_dp)) then
         rcond = 0.0_dp
      else
         rcond = 1.0_dp / cond
      end if
   end subroutine reciprocal_condition_1

end module matrix_functions
