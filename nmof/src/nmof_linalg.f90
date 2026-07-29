! SPDX-License-Identifier: GPL-3.0-only
module nmof_linalg
   use nmof_kinds, only: dp
   implicit none
   private
   public :: solve_linear, solve_linear_matrix, invert_matrix, solve_spd
   public :: eigen_symmetric, cholesky_lower, identity_matrix, vector_norm2
   public :: covariance_matrix, column_means, column_sds, matrix_rank_subset

   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         integer, intent(out) :: ipiv(*)
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
         integer, intent(out) :: info
      end subroutine dgesv
      subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
         integer, intent(out) :: info
      end subroutine dposv
      subroutine dpotrf(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(out) :: info
      end subroutine dpotrf
      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*), work(*)
         integer, intent(out) :: info
      end subroutine dsyev
      subroutine dgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info)
         import dp
         integer, intent(in) :: m, n, lda, lwork
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(inout) :: jpvt(*)
         real(dp), intent(out) :: tau(*), work(*)
         integer, intent(out) :: info
      end subroutine dgeqp3
   end interface
contains
   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure real(dp) function vector_norm2(x) result(v)
      real(dp), intent(in) :: x(:)
      v = sqrt(max(0.0_dp, dot_product(x, x)))
   end function vector_norm2

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: ac(:, :), bc(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(ac(n, n), bc(n, 1), ipiv(n))
      ac = a
      bc(:, 1) = b
      call dgesv(n, 1, ac, n, ipiv, bc, n, info)
      if (info == 0) x = bc(:, 1)
   end subroutine solve_linear

   subroutine solve_linear_matrix(a, b, x, info)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(out) :: x(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: ac(:, :), bc(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n, nrhs
      n = size(a, 1)
      nrhs = size(b, 2)
      if (size(a, 2) /= n .or. size(b, 1) /= n .or. any(shape(x) /= shape(b))) then
         info = -1
         return
      end if
      allocate(ac(n, n), bc(n, nrhs), ipiv(n))
      ac = a
      bc = b
      call dgesv(n, nrhs, ac, n, ipiv, bc, n, info)
      if (info == 0) x = bc
   end subroutine solve_linear_matrix

   subroutine invert_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: ainv(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: b(:, :)
      integer :: n
      n = size(a, 1)
      if (size(a, 2) /= n .or. any(shape(ainv) /= [n, n])) then
         info = -1
         return
      end if
      b = identity_matrix(n)
      call solve_linear_matrix(a, b, ainv, info)
   end subroutine invert_matrix

   subroutine solve_spd(a, b, x, info)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: ac(:, :), bc(:, :)
      integer :: n
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(ac(n, n), bc(n, 1))
      ac = a
      bc(:, 1) = b
      call dposv('U', n, 1, ac, n, bc, n, info)
      if (info == 0) x = bc(:, 1)
   end subroutine solve_spd

   subroutine cholesky_lower(a, l, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: l(:, :)
      integer, intent(out) :: info
      integer :: n, i, j
      n = size(a, 1)
      if (size(a, 2) /= n .or. any(shape(l) /= [n, n])) then
         info = -1
         return
      end if
      l = a
      call dpotrf('L', n, l, n, info)
      if (info == 0) then
         do j = 2, n
            do i = 1, j - 1
               l(i, j) = 0.0_dp
            end do
         end do
      end if
   end subroutine cholesky_lower

   subroutine eigen_symmetric(a, values, vectors, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: values(:), vectors(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: work(:), ac(:, :)
      real(dp) :: query(1)
      integer :: n, lwork
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(values) /= n .or. any(shape(vectors) /= [n, n])) then
         info = -1
         return
      end if
      allocate(ac(n, n))
      ac = a
      call dsyev('V', 'U', n, ac, n, values, query, -1, info)
      if (info /= 0) return
      lwork = max(1, int(query(1)))
      allocate(work(lwork))
      call dsyev('V', 'U', n, ac, n, values, work, lwork, info)
      if (info == 0) vectors = ac
   end subroutine eigen_symmetric

   pure function column_means(x) result(mu)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: mu(size(x, 2))
      if (size(x, 1) > 0) then
         mu = sum(x, dim=1) / real(size(x, 1), dp)
      else
         mu = 0.0_dp
      end if
   end function column_means

   pure function column_sds(x) result(sd)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: sd(size(x, 2)), mu(size(x, 2))
      integer :: n, j
      n = size(x, 1)
      mu = column_means(x)
      sd = 0.0_dp
      if (n > 1) then
         do j = 1, size(x, 2)
            sd(j) = sqrt(sum((x(:, j) - mu(j))**2) / real(n - 1, dp))
         end do
      end if
   end function column_sds

   pure function covariance_matrix(x) result(cov)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: cov(size(x, 2), size(x, 2))
      real(dp) :: mu(size(x, 2))
      integer :: n, i, j
      n = size(x, 1)
      mu = column_means(x)
      cov = 0.0_dp
      if (n > 1) then
         do j = 1, size(x, 2)
            do i = 1, j
               cov(i, j) = dot_product(x(:, i) - mu(i), x(:, j) - mu(j)) / real(n - 1, dp)
               cov(j, i) = cov(i, j)
            end do
         end do
      end if
   end function covariance_matrix

   subroutine matrix_rank_subset(x, columns, multiplier, rank, info, tol)
      real(dp), intent(in) :: x(:, :)
      integer, allocatable, intent(out) :: columns(:)
      real(dp), allocatable, intent(out) :: multiplier(:, :)
      integer, intent(out) :: rank, info
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: a(:, :), tau(:), work(:), b(:, :), sol(:, :)
      real(dp) :: query(1), threshold
      integer, allocatable :: jpvt(:)
      integer :: m, n, lwork, i, j, k, linfo
      m = size(x, 1)
      n = size(x, 2)
      allocate(a(m, n), jpvt(n), tau(min(m, n)))
      a = x
      jpvt = 0
      call dgeqp3(m, n, a, m, jpvt, tau, query, -1, info)
      if (info /= 0) return
      lwork = max(1, int(query(1)))
      allocate(work(lwork))
      a = x
      jpvt = 0
      call dgeqp3(m, n, a, m, jpvt, tau, work, lwork, info)
      if (info /= 0) return
      threshold = max(m, n) * epsilon(1.0_dp) * abs(a(1, 1))
      if (present(tol)) threshold = tol
      rank = 0
      do i = 1, min(m, n)
         if (abs(a(i, i)) > threshold) rank = rank + 1
      end do
      allocate(columns(rank), multiplier(rank, n))
      columns = jpvt(1:rank)
      multiplier = 0.0_dp
      do i = 1, rank
         multiplier(i, columns(i)) = 1.0_dp
      end do
      if (rank < n .and. rank > 0) then
         allocate(b(m, n - rank), sol(rank, n - rank))
         b = x(:, jpvt(rank + 1:n))
         call least_squares_normal(x(:, columns), b, sol, linfo)
         if (linfo /= 0) then
            info = linfo
            return
         end if
         do j = 1, n - rank
            k = jpvt(rank + j)
            multiplier(:, k) = sol(:, j)
         end do
      end if
      info = 0
   contains
      subroutine least_squares_normal(a0, b0, x0, ierr)
         real(dp), intent(in) :: a0(:, :), b0(:, :)
         real(dp), intent(out) :: x0(:, :)
         integer, intent(out) :: ierr
         real(dp), allocatable :: ata(:, :), atb(:, :)
         ata = matmul(transpose(a0), a0)
         atb = matmul(transpose(a0), b0)
         call solve_linear_matrix(ata, atb, x0, ierr)
      end subroutine least_squares_normal
   end subroutine matrix_rank_subset
end module nmof_linalg
