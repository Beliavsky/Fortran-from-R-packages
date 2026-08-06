! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_linalg
   use mgcv_kinds, only : dp
   implicit none
   private
   public :: cholesky_upper, spd_solve, spd_inverse, logdet_spd
   public :: solve_linear, jacobi_eigen, symmetric_root, matrix_rank
   public :: kronecker_product, trace_product, symmetrize

contains

   subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:, :)
      integer :: i, j, n
      n = min(size(a, 1), size(a, 2))
      do j = 1, n
         do i = j + 1, n
            a(i, j) = 0.5_dp * (a(i, j) + a(j, i))
            a(j, i) = a(i, j)
         end do
      end do
   end subroutine symmetrize

   subroutine cholesky_upper(a, r, status, jitter)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: r(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: jitter
      real(dp), allocatable :: work(:, :)
      real(dp) :: s, add
      integer :: i, j, k, n, attempt

      status = 0
      if (size(a, 1) /= size(a, 2)) then
         allocate(r(0, 0)); status = 1; return
      end if
      n = size(a, 1)
      allocate(work(n, n), r(n, n))
      work = a
      call symmetrize(work)
      add = 0.0_dp
      if (present(jitter)) add = max(0.0_dp, jitter)

      do attempt = 0, 7
         r = 0.0_dp
         status = 0
         do j = 1, n
            s = work(j, j) + add
            do k = 1, j - 1
               s = s - r(k, j) * r(k, j)
            end do
            if (s <= epsilon(1.0_dp) * max(1.0_dp, abs(work(j, j)))) then
               status = 2
               exit
            end if
            r(j, j) = sqrt(s)
            do i = j + 1, n
               s = work(j, i)
               do k = 1, j - 1
                  s = s - r(k, j) * r(k, i)
               end do
               r(j, i) = s / r(j, j)
            end do
         end do
         if (status == 0) return
         if (add <= 0.0_dp) then
            add = max(1.0e-12_dp, 1.0e-12_dp * maxval(abs([(work(i, i), i=1,n)])))
         else
            add = 10.0_dp * add
         end if
      end do
   end subroutine cholesky_upper

   subroutine spd_solve(a, b, x, status, jitter)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: jitter
      real(dp), allocatable :: r(:, :), y(:, :)
      integer :: i, j, k, n, nrhs

      if (size(a, 1) /= size(a, 2) .or. size(b, 1) /= size(a, 1)) then
         allocate(x(0, 0)); status = 1; return
      end if
      call cholesky_upper(a, r, status, jitter)
      if (status /= 0) then
         allocate(x(0, 0)); return
      end if
      n = size(a, 1); nrhs = size(b, 2)
      allocate(y(n, nrhs), x(n, nrhs))
      y = b
      do j = 1, nrhs
         do i = 1, n
            do k = 1, i - 1
               y(i, j) = y(i, j) - r(k, i) * y(k, j)
            end do
            y(i, j) = y(i, j) / r(i, i)
         end do
         x(:, j) = y(:, j)
         do i = n, 1, -1
            do k = i + 1, n
               x(i, j) = x(i, j) - r(i, k) * x(k, j)
            end do
            x(i, j) = x(i, j) / r(i, i)
         end do
      end do
   end subroutine spd_solve

   subroutine spd_inverse(a, ainv, status, jitter)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: jitter
      real(dp), allocatable :: eye(:, :)
      integer :: i, n
      if (size(a, 1) /= size(a, 2)) then
         allocate(ainv(0, 0)); status = 1; return
      end if
      n = size(a, 1)
      allocate(eye(n, n)); eye = 0.0_dp
      do i = 1, n
         eye(i, i) = 1.0_dp
      end do
      call spd_solve(a, eye, ainv, status, jitter)
      if (status == 0) call symmetrize(ainv)
   end subroutine spd_inverse

   function logdet_spd(a, status, jitter) result(value)
      real(dp), intent(in) :: a(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: jitter
      real(dp) :: value
      real(dp), allocatable :: r(:, :)
      integer :: i
      call cholesky_upper(a, r, status, jitter)
      if (status /= 0) then
         value = huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, size(r, 1)
         value = value + 2.0_dp * log(r(i, i))
      end do
   end function logdet_spd

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: aa(:, :), bb(:, :), row(:)
      real(dp) :: pivot, factor
      integer :: i, k, p, n, nrhs

      status = 0
      if (size(a, 1) /= size(a, 2) .or. size(b, 1) /= size(a, 1)) then
         allocate(x(0, 0)); status = 1; return
      end if
      n = size(a, 1); nrhs = size(b, 2)
      allocate(aa(n, n), bb(n, nrhs), row(max(n, nrhs)), x(n, nrhs))
      aa = a; bb = b
      do k = 1, n
         p = k
         do i = k + 1, n
            if (abs(aa(i, k)) > abs(aa(p, k))) p = i
         end do
         if (abs(aa(p, k)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aa)))) then
            status = 2; x = 0.0_dp; return
         end if
         if (p /= k) then
            row(1:n) = aa(k, :); aa(k, :) = aa(p, :); aa(p, :) = row(1:n)
            row(1:nrhs) = bb(k, :); bb(k, :) = bb(p, :); bb(p, :) = row(1:nrhs)
         end if
         pivot = aa(k, k)
         aa(k, k:n) = aa(k, k:n) / pivot
         bb(k, :) = bb(k, :) / pivot
         do i = 1, n
            if (i == k) cycle
            factor = aa(i, k)
            if (abs(factor) <= tiny(1.0_dp)) cycle
            aa(i, k:n) = aa(i, k:n) - factor * aa(k, k:n)
            bb(i, :) = bb(i, :) - factor * bb(k, :)
         end do
      end do
      x = bb
   end subroutine solve_linear

   subroutine jacobi_eigen(a, values, vectors, status, tol, max_sweeps)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_sweeps
      real(dp), allocatable :: d(:, :)
      real(dp) :: threshold, app, aqq, apq, tau, t, c, s, dip, diq, vip, viq
      integer :: n, sweep, p, q, i, maxit, idx
      integer, allocatable :: order(:)
      real(dp), allocatable :: tmpv(:), vals(:)

      status = 0
      if (size(a, 1) /= size(a, 2)) then
         allocate(values(0), vectors(0, 0)); status = 1; return
      end if
      n = size(a, 1)
      allocate(d(n, n), vectors(n, n), values(n))
      d = a; call symmetrize(d)
      vectors = 0.0_dp
      do i = 1, n
         vectors(i, i) = 1.0_dp
      end do
      threshold = 1.0e-13_dp; if (present(tol)) threshold = tol
      maxit = max(20, 10 * n * n); if (present(max_sweeps)) maxit = max_sweeps
      do sweep = 1, maxit
         apq = 0.0_dp; p = 1; q = min(2, n)
         do i = 1, n - 1
            do idx = i + 1, n
               if (abs(d(i, idx)) > abs(apq)) then
                  apq = d(i, idx); p = i; q = idx
               end if
            end do
         end do
         if (abs(apq) <= threshold * max(1.0_dp, maxval(abs(d)))) exit
         app = d(p, p); aqq = d(q, q)
         tau = (aqq - app) / (2.0_dp * apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
         else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
         end if
         c = 1.0_dp / sqrt(1.0_dp + t * t); s = t * c
         do i = 1, n
            if (i /= p .and. i /= q) then
               dip = d(i, p); diq = d(i, q)
               d(i, p) = c * dip - s * diq; d(p, i) = d(i, p)
               d(i, q) = s * dip + c * diq; d(q, i) = d(i, q)
            end if
            vip = vectors(i, p); viq = vectors(i, q)
            vectors(i, p) = c * vip - s * viq
            vectors(i, q) = s * vip + c * viq
         end do
         d(p, p) = app - t * apq
         d(q, q) = aqq + t * apq
         d(p, q) = 0.0_dp; d(q, p) = 0.0_dp
      end do
      if (sweep > maxit) status = 2
      do i = 1, n
         values(i) = d(i, i)
      end do
      allocate(order(n), vals(n), tmpv(n))
      order = [(i, i=1,n)]
      do i = 1, n - 1
         idx = i
         do q = i + 1, n
            if (values(order(q)) > values(order(idx))) idx = q
         end do
         if (idx /= i) then
            p = order(i); order(i) = order(idx); order(idx) = p
         end if
      end do
      vals = values(order)
      d = vectors
      do i = 1, n
         vectors(:, i) = d(:, order(i))
      end do
      values = vals
   end subroutine jacobi_eigen

   subroutine symmetric_root(a, root, rank, status, tol)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: root(:, :)
      integer, intent(out) :: rank, status
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: vals(:), vecs(:, :)
      real(dp) :: cutoff
      integer :: i, n
      call jacobi_eigen(a, vals, vecs, status)
      if (status /= 0) then
         allocate(root(0, 0)); rank = 0; return
      end if
      n = size(vals)
      cutoff = max(1.0e-12_dp, sqrt(epsilon(1.0_dp))) * max(1.0_dp, abs(vals(1)))
      if (present(tol)) cutoff = tol * max(1.0_dp, abs(vals(1)))
      rank = count(vals > cutoff)
      allocate(root(n, rank))
      do i = 1, rank
         root(:, i) = vecs(:, i) * sqrt(max(0.0_dp, vals(i)))
      end do
   end subroutine symmetric_root

   integer function matrix_rank(a, tol) result(rank)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: ata(:, :), vals(:), vecs(:, :)
      real(dp) :: cutoff
      integer :: status
      ata = matmul(transpose(a), a)
      call jacobi_eigen(ata, vals, vecs, status)
      if (status /= 0 .or. size(vals) == 0) then
         rank = 0; return
      end if
      cutoff = max(size(a, 1), size(a, 2)) * epsilon(1.0_dp) * max(1.0_dp, sqrt(max(0.0_dp, vals(1))))
      if (present(tol)) cutoff = tol
      rank = count(sqrt(max(0.0_dp, vals)) > cutoff)
   end function matrix_rank

   function kronecker_product(a, b) result(c)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable :: c(:, :)
      integer :: i, j, nr, nc
      nr = size(b, 1); nc = size(b, 2)
      allocate(c(size(a, 1) * nr, size(a, 2) * nc))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            c((i - 1) * nr + 1:i * nr, (j - 1) * nc + 1:j * nc) = a(i, j) * b
         end do
      end do
   end function kronecker_product

   real(dp) function trace_product(a, b) result(value)
      real(dp), intent(in) :: a(:, :), b(:, :)
      integer :: i, j
      value = 0.0_dp
      do j = 1, min(size(a, 2), size(b, 1))
         do i = 1, min(size(a, 1), size(b, 2))
            value = value + a(i, j) * b(j, i)
         end do
      end do
   end function trace_product

end module mgcv_linalg
