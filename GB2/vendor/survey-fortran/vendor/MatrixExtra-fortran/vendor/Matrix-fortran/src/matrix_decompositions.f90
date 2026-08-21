! SPDX-License-Identifier: GPL-3.0-only
module matrix_decompositions
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_shape, matrix_err_singular, &
      matrix_err_not_posdef, matrix_err_convergence
   use matrix_dense, only : eye, frobenius_norm
   implicit none
   private
   public :: lu_factor, lu_solve, solve_linear, inverse_matrix
   public :: determinant, log_determinant
   public :: cholesky_factor, cholesky_solve
   public :: qr_factor, least_squares
   public :: symmetric_eigen, singular_value_decomposition, rank_matrix

contains

   subroutine lu_factor(a, lu, piv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: lu(:,:)
      integer, allocatable, intent(out) :: piv(:)
      integer, intent(out) :: info
      integer :: n, i, j, k, p
      real(dp) :: maxv
      real(dp), allocatable :: row(:)
      if (size(a, 1) /= size(a, 2)) then
         allocate(lu(0, 0), piv(0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      lu = a
      allocate(piv(n), row(n))
      piv = [(i, i = 1, n)]
      info = matrix_success
      do k = 1, n - 1
         p = k
         maxv = abs(lu(k, k))
         do i = k + 1, n
            if (abs(lu(i, k)) > maxv) then
               maxv = abs(lu(i, k))
               p = i
            end if
         end do
         if (maxv <= tiny(1.0_dp)) then
            info = matrix_err_singular
            return
         end if
         if (p /= k) then
            row = lu(k, :)
            lu(k, :) = lu(p, :)
            lu(p, :) = row
            i = piv(k)
            piv(k) = piv(p)
            piv(p) = i
         end if
         do i = k + 1, n
            lu(i, k) = lu(i, k) / lu(k, k)
            do j = k + 1, n
               lu(i, j) = lu(i, j) - lu(i, k) * lu(k, j)
            end do
         end do
      end do
      if (n > 0 .and. abs(lu(n, n)) <= tiny(1.0_dp)) info = matrix_err_singular
   end subroutine lu_factor

   subroutine lu_solve(lu, piv, b, x, info)
      real(dp), intent(in) :: lu(:,:), b(:,:)
      integer, intent(in) :: piv(:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      integer :: n, nrhs, i, k
      if (size(lu, 1) /= size(lu, 2) .or. size(piv) /= size(lu, 1) .or. &
         size(b, 1) /= size(lu, 1)) then
         allocate(x(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(lu, 1)
      nrhs = size(b, 2)
      allocate(x(n, nrhs))
      do i = 1, n
         x(i, :) = b(piv(i), :)
      end do
      do i = 2, n
         do k = 1, i - 1
            x(i, :) = x(i, :) - lu(i, k) * x(k, :)
         end do
      end do
      do i = n, 1, -1
         if (abs(lu(i, i)) <= tiny(1.0_dp)) then
            info = matrix_err_singular
            return
         end if
         do k = i + 1, n
            x(i, :) = x(i, :) - lu(i, k) * x(k, :)
         end do
         x(i, :) = x(i, :) / lu(i, i)
      end do
      info = matrix_success
   end subroutine lu_solve

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: lu(:,:)
      integer, allocatable :: piv(:)
      call lu_factor(a, lu, piv, info)
      if (info /= matrix_success) then
         allocate(x(0, 0))
         return
      end if
      call lu_solve(lu, piv, b, x, info)
   end subroutine solve_linear

   subroutine inverse_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: id(:,:)
      if (size(a, 1) /= size(a, 2)) then
         allocate(ainv(0, 0))
         info = matrix_err_shape
         return
      end if
      id = eye(size(a, 1))
      call solve_linear(a, id, ainv, info)
   end subroutine inverse_matrix

   subroutine determinant(a, det, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: det
      integer, intent(out) :: info
      real(dp), allocatable :: lu(:,:)
      integer, allocatable :: piv(:)
      integer :: i, swaps
      call lu_factor(a, lu, piv, info)
      if (info /= matrix_success) then
         det = 0.0_dp
         return
      end if
      det = 1.0_dp
      swaps = 0
      do i = 1, size(lu, 1)
         det = det * lu(i, i)
         if (piv(i) /= i) swaps = swaps + 1
      end do
      ! The cycle parity, not simply the number of displaced entries, sets the sign.
      swaps = permutation_swaps(piv)
      if (mod(swaps, 2) /= 0) det = -det
   end subroutine determinant

   subroutine log_determinant(a, logabsdet, sign_det, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: logabsdet, sign_det
      integer, intent(out) :: info
      real(dp), allocatable :: lu(:,:)
      integer, allocatable :: piv(:)
      integer :: i, swaps
      call lu_factor(a, lu, piv, info)
      if (info /= matrix_success) then
         logabsdet = -huge(1.0_dp)
         sign_det = 0.0_dp
         return
      end if
      logabsdet = 0.0_dp
      sign_det = 1.0_dp
      do i = 1, size(lu, 1)
         if (abs(lu(i, i)) <= tiny(1.0_dp)) then
            info = matrix_err_singular
            logabsdet = -huge(1.0_dp)
            sign_det = 0.0_dp
            return
         end if
         logabsdet = logabsdet + log(abs(lu(i, i)))
         if (lu(i, i) < 0.0_dp) sign_det = -sign_det
      end do
      swaps = permutation_swaps(piv)
      if (mod(swaps, 2) /= 0) sign_det = -sign_det
   end subroutine log_determinant

   integer function permutation_swaps(p) result(nswap)
      integer, intent(in) :: p(:)
      logical, allocatable :: seen(:)
      integer :: i, j, cycle_len
      allocate(seen(size(p)), source=.false.)
      nswap = 0
      do i = 1, size(p)
         if (.not. seen(i)) then
            j = i
            cycle_len = 0
            do while (.not. seen(j))
               seen(j) = .true.
               cycle_len = cycle_len + 1
               j = p(j)
            end do
            if (cycle_len > 0) nswap = nswap + cycle_len - 1
         end if
      end do
   end function permutation_swaps

   subroutine cholesky_factor(a, l, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer :: n, i, j, k
      real(dp) :: s, eps
      if (size(a, 1) /= size(a, 2)) then
         allocate(l(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      eps = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))
      if (present(tol)) eps = max(0.0_dp, tol)
      allocate(l(n, n), source=0.0_dp)
      do j = 1, n
         s = a(j, j)
         do k = 1, j - 1
            s = s - l(j, k) * l(j, k)
         end do
         if (s <= eps) then
            info = matrix_err_not_posdef
            return
         end if
         l(j, j) = sqrt(s)
         do i = j + 1, n
            s = a(i, j)
            do k = 1, j - 1
               s = s - l(i, k) * l(j, k)
            end do
            l(i, j) = s / l(j, j)
         end do
      end do
      info = matrix_success
   end subroutine cholesky_factor

   subroutine cholesky_solve(l, b, x, info)
      real(dp), intent(in) :: l(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      integer :: n, i, k
      if (size(l, 1) /= size(l, 2) .or. size(b, 1) /= size(l, 1)) then
         allocate(x(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(l, 1)
      x = b
      do i = 1, n
         if (abs(l(i, i)) <= tiny(1.0_dp)) then
            info = matrix_err_singular
            return
         end if
         do k = 1, i - 1
            x(i, :) = x(i, :) - l(i, k) * x(k, :)
         end do
         x(i, :) = x(i, :) / l(i, i)
      end do
      do i = n, 1, -1
         do k = i + 1, n
            x(i, :) = x(i, :) - l(k, i) * x(k, :)
         end do
         x(i, :) = x(i, :) / l(i, i)
      end do
      info = matrix_success
   end subroutine cholesky_solve

   subroutine qr_factor(a, q, r, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: q(:,:), r(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer :: m, n, i, j
      real(dp) :: nrm, eps
      m = size(a, 1)
      n = size(a, 2)
      allocate(q(m, n), r(n, n), source=0.0_dp)
      eps = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, frobenius_norm(a))
      if (present(tol)) eps = max(0.0_dp, tol)
      do j = 1, n
         q(:, j) = a(:, j)
         do i = 1, j - 1
            r(i, j) = dot_product(q(:, i), q(:, j))
            q(:, j) = q(:, j) - r(i, j) * q(:, i)
         end do
         ! Reorthogonalization improves modified Gram-Schmidt robustness.
         do i = 1, j - 1
            nrm = dot_product(q(:, i), q(:, j))
            r(i, j) = r(i, j) + nrm
            q(:, j) = q(:, j) - nrm * q(:, i)
         end do
         r(j, j) = sqrt(dot_product(q(:, j), q(:, j)))
         if (r(j, j) <= eps) then
            info = matrix_err_singular
            return
         end if
         q(:, j) = q(:, j) / r(j, j)
      end do
      info = matrix_success
   end subroutine qr_factor

   subroutine least_squares(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: q(:,:), r(:,:), y(:,:)
      integer :: i, k, n
      if (size(a, 1) /= size(b, 1) .or. size(a, 1) < size(a, 2)) then
         allocate(x(0, 0))
         info = matrix_err_shape
         return
      end if
      call qr_factor(a, q, r, info)
      if (info /= matrix_success) then
         allocate(x(0, 0))
         return
      end if
      y = matmul(transpose(q), b)
      n = size(a, 2)
      x = y
      do i = n, 1, -1
         do k = i + 1, n
            x(i, :) = x(i, :) - r(i, k) * x(k, :)
         end do
         x(i, :) = x(i, :) / r(i, i)
      end do
   end subroutine least_squares

   subroutine symmetric_eigen(a, values, vectors, info, tol, max_iter)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: b(:,:)
      real(dp) :: app, aqq, apq, tau, t, c, s, eps, off
      integer :: n, p, q, i, iter, maxit
      if (size(a, 1) /= size(a, 2)) then
         allocate(values(0), vectors(0, 0))
         info = matrix_err_shape
         return
      end if
      n = size(a, 1)
      b = 0.5_dp * (a + transpose(a))
      vectors = eye(n)
      eps = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(b)))
      if (present(tol)) eps = max(0.0_dp, tol)
      maxit = max(30, 50 * n * n)
      if (present(max_iter)) maxit = max_iter
      info = matrix_err_convergence
      do iter = 1, maxit
         call largest_offdiag(b, p, q, off)
         if (off <= eps) then
            info = matrix_success
            exit
         end if
         app = b(p, p)
         aqq = b(q, q)
         apq = b(p, q)
         tau = (aqq - app) / (2.0_dp * apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
         else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
         end if
         c = 1.0_dp / sqrt(1.0_dp + t * t)
         s = t * c
         do i = 1, n
            if (i /= p .and. i /= q) then
               app = b(i, p)
               aqq = b(i, q)
               b(i, p) = c * app - s * aqq
               b(p, i) = b(i, p)
               b(i, q) = s * app + c * aqq
               b(q, i) = b(i, q)
            end if
         end do
         app = b(p, p)
         aqq = b(q, q)
         apq = b(p, q)
         b(p, p) = c * c * app - 2.0_dp * c * s * apq + s * s * aqq
         b(q, q) = s * s * app + 2.0_dp * c * s * apq + c * c * aqq
         b(p, q) = 0.0_dp
         b(q, p) = 0.0_dp
         do i = 1, n
            app = vectors(i, p)
            aqq = vectors(i, q)
            vectors(i, p) = c * app - s * aqq
            vectors(i, q) = s * app + c * aqq
         end do
      end do
      allocate(values(n))
      do i = 1, n
         values(i) = b(i, i)
      end do
      call sort_eigenpairs(values, vectors, descending=.false.)
   end subroutine symmetric_eigen

   subroutine largest_offdiag(a, p, q, value)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out) :: p, q
      real(dp), intent(out) :: value
      integer :: i, j
      p = 1
      q = min(2, size(a, 1))
      value = 0.0_dp
      do j = 2, size(a, 2)
         do i = 1, j - 1
            if (abs(a(i, j)) > value) then
               value = abs(a(i, j))
               p = i
               q = j
            end if
         end do
      end do
   end subroutine largest_offdiag

   subroutine sort_eigenpairs(values, vectors, descending)
      real(dp), intent(inout) :: values(:), vectors(:,:)
      logical, intent(in) :: descending
      integer :: i, j, k
      real(dp) :: tmp
      real(dp), allocatable :: col(:)
      allocate(col(size(vectors, 1)))
      do i = 1, size(values) - 1
         k = i
         do j = i + 1, size(values)
            if (descending .neqv. (values(j) < values(k))) k = j
         end do
         if (k /= i) then
            tmp = values(i)
            values(i) = values(k)
            values(k) = tmp
            col = vectors(:, i)
            vectors(:, i) = vectors(:, k)
            vectors(:, k) = col
         end if
      end do
   end subroutine sort_eigenpairs

   subroutine singular_value_decomposition(a, u, s, vt, info, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: u(:,:), s(:), vt(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: ata(:,:), values(:), v(:,:), q(:,:), r(:,:)
      real(dp) :: eps
      integer :: m, n, k, i, qinfo
      m = size(a, 1)
      n = size(a, 2)
      k = min(m, n)
      ata = matmul(transpose(a), a)
      call symmetric_eigen(ata, values, v, info)
      if (info /= matrix_success) then
         allocate(u(0, 0), s(0), vt(0, 0))
         return
      end if
      call sort_eigenpairs(values, v, descending=.true.)
      allocate(s(k), u(m, k), vt(k, n))
      eps = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, sqrt(maxval(max(values, 0.0_dp))))
      if (present(tol)) eps = max(0.0_dp, tol)
      do i = 1, k
         s(i) = sqrt(max(0.0_dp, values(i)))
         vt(i, :) = v(:, i)
         if (s(i) > eps) then
            u(:, i) = matmul(a, v(:, i)) / s(i)
         else
            u(:, i) = 0.0_dp
            if (i <= m) u(i, i) = 1.0_dp
         end if
      end do
      ! Reorthogonalize U to limit loss of orthogonality for repeated singular values.
      call qr_factor(u, q, r, qinfo, tol=eps * 0.01_dp)
      if (qinfo == matrix_success) u = q
      info = matrix_success
   end subroutine singular_value_decomposition

   integer function rank_matrix(a, tol, info) result(rnk)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: info
      real(dp), allocatable :: u(:,:), s(:), vt(:,:)
      real(dp) :: eps
      integer :: istat
      call singular_value_decomposition(a, u, s, vt, istat)
      if (istat /= matrix_success) then
         rnk = 0
      else
         if (size(s) == 0) then
            eps = 0.0_dp
         else
            eps = real(max(size(a, 1), size(a, 2)), dp) * epsilon(1.0_dp) * maxval(s)
         end if
         if (present(tol)) eps = max(0.0_dp, tol)
         rnk = count(s > eps)
      end if
      if (present(info)) info = istat
   end function rank_matrix

end module matrix_decompositions
