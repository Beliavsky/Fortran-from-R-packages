! SPDX-License-Identifier: GPL-3.0-or-later
! Eigenvalue and singular-value routines translated from sanic.
module sanic_spectral
   use sanic_kinds, only: dp
   use sanic_direct, only: dense_is_symmetric
   implicit none
   private

   integer, parameter, public :: sanic_svd_bdc = 1
   integer, parameter, public :: sanic_svd_jacobi = 2

   type, public :: eigen_result
      complex(dp), allocatable :: values(:)
      complex(dp), allocatable :: vectors(:, :)
      integer :: info = 0
      logical :: symmetric = .false.
   end type eigen_result

   type, public :: svd_result
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: v(:, :)
      integer :: info = 0
      logical :: thin = .true.
      logical :: vectors = .true.
      integer :: method = sanic_svd_bdc
   end type svd_result

   public :: eigen2
   public :: svd2
   public :: symmetric_eigen_real
   public :: general_eigen_complex

contains

   pure subroutine eigen2(a, out, symmetric, vectors, tol, max_iter)
      !! Compute symmetric or general eigenvalues and optional eigenvectors.
      real(dp), intent(in) :: a(:, :) !! Square dense matrix whose eigenpairs are requested.
      type(eigen_result), intent(out) :: out !! Eigenvalues, optional vectors, symmetry flag, and status.
      logical, intent(in), optional :: symmetric !! If absent, symmetry is detected at machine precision.
      logical, intent(in), optional :: vectors !! Whether eigenvectors are returned; default true.
      real(dp), intent(in), optional :: tol !! Positive convergence tolerance; default scales machine epsilon.
      integer, intent(in), optional :: max_iter !! Maximum QR/Jacobi iteration budget; default depends on matrix order.
      logical :: compute_vectors, use_symmetric
      real(dp), allocatable :: real_values(:), real_vectors(:, :)
      complex(dp), allocatable :: complex_values(:), complex_vectors(:, :)
      integer :: status

      if (size(a, 1) /= size(a, 2) .or. size(a, 1) < 1) then
         out%info = 1
         allocate(out%values(0), out%vectors(0, 0))
         return
      end if
      use_symmetric = dense_is_symmetric(a)
      if (present(symmetric)) use_symmetric = symmetric
      compute_vectors = .true.
      if (present(vectors)) compute_vectors = vectors
      out%symmetric = use_symmetric

      if (use_symmetric) then
         call symmetric_eigen_real(a, real_values, real_vectors, status, tol, max_iter)
         out%info = status
         out%values = cmplx(real_values, 0.0_dp, kind=dp)
         if (compute_vectors) then
            out%vectors = cmplx(real_vectors, 0.0_dp, kind=dp)
         else
            allocate(out%vectors(0, 0))
         end if
      else
         call general_eigen_complex(a, complex_values, complex_vectors, status, tol, max_iter, compute_vectors)
         out%info = status
         out%values = complex_values
         if (compute_vectors) then
            out%vectors = complex_vectors
         else
            allocate(out%vectors(0, 0))
         end if
      end if
   end subroutine eigen2

   pure subroutine symmetric_eigen_real(a, values, vectors, info, tol, max_iter)
      !! Self-adjoint eigensolver using Householder tridiagonalization and implicit QL, with Jacobi fallback.
      real(dp), intent(in) :: a(:, :) !! Real square symmetric matrix.
      real(dp), allocatable, intent(out) :: values(:) !! Eigenvalues sorted in ascending order.
      real(dp), allocatable, intent(out) :: vectors(:, :) !! Orthonormal eigenvectors in matching columns.
      integer, intent(out) :: info !! Status: zero converged, one invalid shape, two if both iterative methods fail.
      real(dp), intent(in), optional :: tol !! Relative tridiagonal deflation tolerance.
      integer, intent(in), optional :: max_iter !! Maximum implicit-QL steps per eigenvalue before Jacobi fallback.
      real(dp), allocatable :: diagonal(:), offdiagonal(:), q(:, :), work(:, :)
      integer :: status

      if (size(a, 1) /= size(a, 2) .or. size(a, 1) < 1) then
         allocate(values(0), vectors(0, 0))
         info = 1
         return
      end if
      call selfadjoint_from_lower(a, work)
      call symmetric_tridiagonalize(work, diagonal, offdiagonal, q)
      call tridiagonal_ql(diagonal, offdiagonal, q, status, tol, max_iter)
      if (status == 0) then
         values = diagonal
         vectors = q
         call sort_real_eigenpairs(values, vectors, ascending=.true.)
         info = 0
         return
      end if
      call symmetric_eigen_jacobi(work, values, vectors, info, tol, max_iter)
   end subroutine symmetric_eigen_real

   pure subroutine selfadjoint_from_lower(a, symmetric_a)
      !! Form the self-adjoint matrix represented by the lower triangle of A, matching Eigen semantics.
      real(dp), intent(in) :: a(:, :) !! Square input whose diagonal and lower triangle define the self-adjoint matrix.
      real(dp), allocatable, intent(out) :: symmetric_a(:, :) !! Symmetric matrix formed by reflecting A's lower triangle.
      integer :: i, j, n

      n = size(a, 1)
      allocate(symmetric_a(n, n))
      symmetric_a = 0.0_dp
      do j = 1, n
         symmetric_a(j, j) = a(j, j)
         do i = j + 1, n
            symmetric_a(i, j) = a(i, j)
            symmetric_a(j, i) = a(i, j)
         end do
      end do
   end subroutine selfadjoint_from_lower

   pure subroutine symmetric_tridiagonalize(a, diagonal, offdiagonal, q)
      !! Reduce a real symmetric matrix to tridiagonal form while accumulating the orthogonal similarity transform.
      real(dp), intent(in) :: a(:, :) !! Real symmetric square matrix to reduce.
      real(dp), allocatable, intent(out) :: diagonal(:) !! Tridiagonal diagonal entries.
      real(dp), allocatable, intent(out) :: offdiagonal(:) !! Subdiagonal entries stored with offdiagonal(1)=0.
      real(dp), allocatable, intent(out) :: q(:, :) !! Orthogonal transformation whose columns map tridiagonal eigenvectors.
      real(dp) :: f, g, h, hh, scale
      integer :: i, j, k, l, n

      n = size(a, 1)
      q = a
      allocate(diagonal(n), offdiagonal(n))
      diagonal = 0.0_dp
      offdiagonal = 0.0_dp
      do i = n, 2, -1
         l = i - 1
         h = 0.0_dp
         scale = 0.0_dp
         if (l > 1) then
            scale = sum(abs(q(i, 1:l)))
            if (scale <= tiny(1.0_dp)) then
               offdiagonal(i) = q(i, l)
            else
               q(i, 1:l) = q(i, 1:l)/scale
               h = dot_product(q(i, 1:l), q(i, 1:l))
               f = q(i, l)
               g = -sign(sqrt(h), f)
               offdiagonal(i) = scale*g
               h = h - f*g
               q(i, l) = f - g
               offdiagonal(1:l) = 0.0_dp
               do j = 1, l
                  q(j, i) = q(i, j)/h
                  g = dot_product(q(j, 1:j), q(i, 1:j))
                  if (j < l) g = g + dot_product(q(j + 1:l, j), q(i, j + 1:l))
                  offdiagonal(j) = g/h
               end do
               f = dot_product(offdiagonal(1:l), q(i, 1:l))
               hh = f/(2.0_dp*h)
               do j = 1, l
                  f = q(i, j)
                  g = offdiagonal(j) - hh*f
                  offdiagonal(j) = g
                  do k = 1, j
                     q(j, k) = q(j, k) - f*offdiagonal(k) - g*q(i, k)
                  end do
               end do
            end if
         else
            offdiagonal(i) = q(i, l)
         end if
         diagonal(i) = h
      end do
      diagonal(1) = 0.0_dp
      offdiagonal(1) = 0.0_dp

      do i = 1, n
         l = i - 1
         if (abs(diagonal(i)) > tiny(1.0_dp) .and. l > 0) then
            do j = 1, l
               g = dot_product(q(i, 1:l), q(1:l, j))
               q(1:l, j) = q(1:l, j) - g*q(1:l, i)
            end do
         end if
         diagonal(i) = q(i, i)
         q(i, i) = 1.0_dp
         if (l > 0) then
            q(i, 1:l) = 0.0_dp
            q(1:l, i) = 0.0_dp
         end if
      end do
   end subroutine symmetric_tridiagonalize

   pure subroutine tridiagonal_ql(diagonal, offdiagonal, vectors, info, tol, max_iter)
      !! Diagonalize a symmetric tridiagonal matrix by implicit-shift QL while accumulating eigenvectors.
      real(dp), intent(inout) :: diagonal(:) !! Tridiagonal diagonal entries, replaced by eigenvalues.
      real(dp), intent(inout) :: offdiagonal(:) !! Subdiagonal workspace from symmetric_tridiagonalize.
      real(dp), intent(inout) :: vectors(:, :) !! Orthogonal transformation updated to eigenvectors.
      integer, intent(out) :: info !! Status: zero success, two if an eigenvalue exceeds its QL iteration budget.
      real(dp), intent(in), optional :: tol !! Relative deflation tolerance; default scales machine epsilon.
      integer, intent(in), optional :: max_iter !! Maximum QL steps for each active eigenvalue; default 100.
      real(dp) :: b, c, f, g, p, r, relative_tol, ss
      integer :: i, iteration, k, l, limit, m, n
      logical :: restart_iteration

      n = size(diagonal)
      if (n <= 1) then
         info = 0
         return
      end if
      offdiagonal(1:n - 1) = offdiagonal(2:n)
      offdiagonal(n) = 0.0_dp
      relative_tol = epsilon(1.0_dp)*real(max(1, n), dp)
      if (present(tol)) then
         if (tol > 0.0_dp) relative_tol = tol
      end if
      limit = 100
      if (present(max_iter)) then
         if (max_iter > 0) limit = max_iter
      end if

      info = 0
      do l = 1, n
         iteration = 0
         do
            m = l
            do while (m < n)
               if (abs(offdiagonal(m)) <= relative_tol*(1.0_dp + abs(diagonal(m)) + abs(diagonal(m + 1)))) exit
               m = m + 1
            end do
            if (m == l) exit
            iteration = iteration + 1
            if (iteration > limit) then
               info = 2
               return
            end if

            g = (diagonal(l + 1) - diagonal(l))/(2.0_dp*offdiagonal(l))
            r = hypot(g, 1.0_dp)
            g = diagonal(m) - diagonal(l) + offdiagonal(l)/(g + sign(r, g))
            ss = 1.0_dp
            c = 1.0_dp
            p = 0.0_dp
            restart_iteration = .false.
            do i = m - 1, l, -1
               f = ss*offdiagonal(i)
               b = c*offdiagonal(i)
               r = hypot(f, g)
               offdiagonal(i + 1) = r
               if (r <= tiny(1.0_dp)) then
                  diagonal(i + 1) = diagonal(i + 1) - p
                  offdiagonal(m) = 0.0_dp
                  restart_iteration = .true.
                  exit
               end if
               ss = f/r
               c = g/r
               g = diagonal(i + 1) - p
               r = (diagonal(i) - g)*ss + 2.0_dp*c*b
               p = ss*r
               diagonal(i + 1) = g + p
               g = c*r - b
               do k = 1, n
                  f = vectors(k, i + 1)
                  vectors(k, i + 1) = ss*vectors(k, i) + c*f
                  vectors(k, i) = c*vectors(k, i) - ss*f
               end do
            end do
            if (restart_iteration) cycle
            diagonal(l) = diagonal(l) - p
            offdiagonal(l) = g
            offdiagonal(m) = 0.0_dp
         end do
      end do
   end subroutine tridiagonal_ql

   pure subroutine symmetric_eigen_jacobi(a, values, vectors, info, tol, max_iter)
      !! Fallback Jacobi eigensolver for a real symmetric matrix, returning ascending eigenvalues.
      real(dp), intent(in) :: a(:, :) !! Real square symmetric matrix.
      real(dp), allocatable, intent(out) :: values(:) !! Eigenvalues sorted in ascending order.
      real(dp), allocatable, intent(out) :: vectors(:, :) !! Orthonormal eigenvectors in matching columns.
      integer, intent(out) :: info !! Status: zero converged, one invalid shape, two iteration limit reached.
      real(dp), intent(in), optional :: tol !! Relative off-diagonal convergence tolerance.
      integer, intent(in), optional :: max_iter !! Maximum Jacobi rotations; default 100*n*n.
      real(dp), allocatable :: work(:, :)
      real(dp) :: app, apq, aqq, c, scale, s, t, tau, threshold, vip, viq, wip, wiq
      integer :: i, iteration, limit, n, p, q

      n = size(a, 1)
      if (size(a, 2) /= n .or. n < 1) then
         allocate(values(0), vectors(0, 0))
         info = 1
         return
      end if
      work = 0.5_dp*(a + transpose(a))
      allocate(vectors(n, n))
      vectors = 0.0_dp
      do i = 1, n
         vectors(i, i) = 1.0_dp
      end do
      scale = max(1.0_dp, maxval(abs(work)))
      threshold = epsilon(1.0_dp)*scale*real(max(1, n), dp)
      if (present(tol)) then
         if (tol > 0.0_dp) threshold = tol*scale
      end if
      limit = 100*n*n
      if (present(max_iter)) then
         if (max_iter > 0) limit = max_iter
      end if

      info = 2
      do iteration = 1, limit
         if (n == 1) then
            info = 0
            exit
         end if
         call largest_offdiagonal(work, p, q, apq)
         if (abs(apq) <= threshold) then
            info = 0
            exit
         end if
         app = work(p, p)
         aqq = work(q, q)
         tau = (aqq - app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau + sqrt(1.0_dp + tau*tau))
         else
            t = -1.0_dp/(-tau + sqrt(1.0_dp + tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp + t*t)
         s = t*c

         do i = 1, n
            if (i /= p .and. i /= q) then
               wip = work(i, p)
               wiq = work(i, q)
               work(i, p) = c*wip - s*wiq
               work(p, i) = work(i, p)
               work(i, q) = s*wip + c*wiq
               work(q, i) = work(i, q)
            end if
            vip = vectors(i, p)
            viq = vectors(i, q)
            vectors(i, p) = c*vip - s*viq
            vectors(i, q) = s*vip + c*viq
         end do
         work(p, p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         work(q, q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         work(p, q) = 0.0_dp
         work(q, p) = 0.0_dp
      end do

      allocate(values(n))
      do i = 1, n
         values(i) = work(i, i)
      end do
      call sort_real_eigenpairs(values, vectors, ascending=.true.)
   end subroutine symmetric_eigen_jacobi

   pure subroutine general_eigen_complex(a, values, vectors, info, tol, max_iter, compute_vectors)
      !! Shifted complex QR eigensolver with inverse iteration for right eigenvectors.
      real(dp), intent(in) :: a(:, :) !! Real square matrix; complex eigenpairs are supported.
      complex(dp), allocatable, intent(out) :: values(:) !! Complex eigenvalue estimates.
      complex(dp), allocatable, intent(out) :: vectors(:, :) !! Right eigenvectors, or zero-size when not requested.
      integer, intent(out) :: info !! Status: zero converged, one invalid shape, two QR convergence limit reached.
      real(dp), intent(in), optional :: tol !! Relative QR deflation tolerance.
      integer, intent(in), optional :: max_iter !! Iteration budget per active trailing block.
      logical, intent(in), optional :: compute_vectors !! Whether to compute right eigenvectors; default true.
      complex(dp), allocatable :: h(:, :), q(:, :), r(:, :), block(:, :)
      real(dp), allocatable :: h_real(:, :)
      complex(dp) :: mu
      real(dp) :: scale, threshold
      integer :: active, i, iteration, limit, n, status
      logical :: want_vectors

      n = size(a, 1)
      if (size(a, 2) /= n .or. n < 1) then
         allocate(values(0), vectors(0, 0))
         info = 1
         return
      end if
      call hessenberg_reduce_real(a, h_real)
      h = cmplx(h_real, 0.0_dp, kind=dp)
      scale = max(1.0_dp, maxval(abs(a)))
      threshold = epsilon(1.0_dp)*scale*real(max(1, n), dp)
      if (present(tol)) then
         if (tol > 0.0_dp) threshold = tol*scale
      end if
      limit = max(200, 100*n*n)
      if (present(max_iter)) then
         if (max_iter > 0) limit = max_iter
      end if
      info = 0

      active = n
      do while (active > 1)
         iteration = 0
         do while (abs(h(active, active - 1)) > &
                   threshold*(1.0_dp + abs(h(active - 1, active - 1)) + abs(h(active, active))))
            iteration = iteration + 1
            if (iteration > limit) then
               info = 2
               exit
            end if
            mu = wilkinson_shift_complex(h(active - 1, active - 1), h(active - 1, active), &
                                         h(active, active - 1), h(active, active))
            block = h(1:active, 1:active)
            do i = 1, active
               block(i, i) = block(i, i) - mu
            end do
            call complex_qr(block, q, r, status)
            if (status /= 0) then
               mu = h(active, active) + cmplx(0.0_dp, threshold, kind=dp)
               block = h(1:active, 1:active)
               do i = 1, active
                  block(i, i) = block(i, i) - mu
               end do
               call complex_qr(block, q, r, status)
               if (status /= 0) then
                  info = 2
                  exit
               end if
            end if
            h(1:active, 1:active) = matmul(r, q)
            do i = 1, active
               h(i, i) = h(i, i) + mu
            end do
            call zero_below_hessenberg(h, active, threshold)
         end do
         if (info /= 0) exit
         h(active, active - 1) = cmplx(0.0_dp, 0.0_dp, kind=dp)
         active = active - 1
      end do

      allocate(values(n))
      do i = 1, n
         values(i) = h(i, i)
      end do
      want_vectors = .true.
      if (present(compute_vectors)) want_vectors = compute_vectors
      if (want_vectors) then
         call inverse_iteration_vectors(a, values, vectors)
      else
         allocate(vectors(0, 0))
      end if
   end subroutine general_eigen_complex

   pure subroutine hessenberg_reduce_real(a, h)
      !! Reduce a real square matrix to upper Hessenberg form by Householder similarity transformations.
      real(dp), intent(in) :: a(:, :) !! Real square matrix to reduce.
      real(dp), allocatable, intent(out) :: h(:, :) !! Orthogonally similar upper-Hessenberg matrix.
      real(dp), allocatable :: reflector(:)
      real(dp) :: alpha, beta, norm_value, projection, reflector_norm, threshold
      integer :: i, j, k, n

      n = size(a, 1)
      allocate(h(n, n))
      h = a
      if (n <= 2) return
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(1, n), dp)
      do k = 1, n - 2
         reflector = h(k + 1:n, k)
         norm_value = sqrt(max(0.0_dp, dot_product(reflector, reflector)))
         if (norm_value <= threshold) cycle
         alpha = -sign(norm_value, reflector(1))
         reflector(1) = reflector(1) - alpha
         reflector_norm = dot_product(reflector, reflector)
         if (reflector_norm <= tiny(1.0_dp)) cycle
         beta = 2.0_dp/reflector_norm

         do j = k, n
            projection = beta*dot_product(reflector, h(k + 1:n, j))
            h(k + 1:n, j) = h(k + 1:n, j) - reflector*projection
         end do
         do i = 1, n
            projection = beta*dot_product(h(i, k + 1:n), reflector)
            h(i, k + 1:n) = h(i, k + 1:n) - projection*reflector
         end do
         h(k + 1, k) = alpha
         if (k + 2 <= n) h(k + 2:n, k) = 0.0_dp
      end do
   end subroutine hessenberg_reduce_real

   pure subroutine zero_below_hessenberg(h, active, threshold)
      !! Remove roundoff fill below the first subdiagonal of an active Hessenberg block.
      complex(dp), intent(inout) :: h(:, :) !! Complex Hessenberg iterate to clean in place.
      integer, intent(in) :: active !! Leading active block order.
      real(dp), intent(in) :: threshold !! Scale-aware threshold used to discard roundoff fill.
      integer :: i, j

      do j = 1, active - 2
         do i = j + 2, active
            if (abs(h(i, j)) <= 100.0_dp*threshold) h(i, j) = cmplx(0.0_dp, 0.0_dp, kind=dp)
         end do
      end do
   end subroutine zero_below_hessenberg

   pure subroutine svd2(a, out, method, vectors, thin, tol, max_iter)
      !! Compute singular values and optional thin/full singular vectors with distinct BDC-like and Jacobi backends.
      real(dp), intent(in) :: a(:, :) !! Dense matrix to decompose.
      type(svd_result), intent(out) :: out !! Singular values/vectors and status information.
      integer, intent(in), optional :: method !! sanic_svd_bdc selects bidiagonal QR; sanic_svd_jacobi selects Jacobi SVD.
      logical, intent(in), optional :: vectors !! Whether singular vectors are returned; default true.
      logical, intent(in), optional :: thin !! Whether vectors use economy shapes; default true.
      real(dp), intent(in), optional :: tol !! Relative pair-orthogonality tolerance for Jacobi sweeps.
      integer, intent(in), optional :: max_iter !! Maximum Jacobi sweeps; default is 100.
      real(dp), allocatable :: tall_u(:, :), tall_v(:, :)
      integer :: m, n, r, status
      logical :: thin_vectors, want_vectors

      m = size(a, 1)
      n = size(a, 2)
      r = min(m, n)
      out%method = sanic_svd_bdc
      if (present(method)) out%method = method
      want_vectors = .true.
      if (present(vectors)) want_vectors = vectors
      thin_vectors = .true.
      if (present(thin)) thin_vectors = thin
      out%vectors = want_vectors
      out%thin = thin_vectors
      if (m < 1 .or. n < 1 .or. out%method < sanic_svd_bdc .or. out%method > sanic_svd_jacobi) then
         out%info = 1
         allocate(out%d(0), out%u(0, 0), out%v(0, 0))
         return
      end if

      if (m >= n) then
         if (out%method == sanic_svd_bdc) then
            call bidiagonal_qr_tall(a, out%d, tall_u, tall_v, status, tol, max_iter)
         else
            call jacobi_svd_tall(a, out%d, tall_u, tall_v, status, tol, max_iter)
         end if
         if (status /= 0) then
            out%info = status
            allocate(out%u(0, 0), out%v(0, 0))
            return
         end if
         if (.not. want_vectors) then
            allocate(out%u(0, 0), out%v(0, 0))
            out%info = 0
            return
         end if
         if (thin_vectors) then
            out%u = tall_u
         else
            allocate(out%u(m, m))
            out%u = 0.0_dp
            out%u(:, 1:r) = tall_u
            if (m > r) call complete_basis(out%u, r)
         end if
         out%v = tall_v
      else
         if (out%method == sanic_svd_bdc) then
            call bidiagonal_qr_tall(transpose(a), out%d, tall_u, tall_v, status, tol, max_iter)
         else
            call jacobi_svd_tall(transpose(a), out%d, tall_u, tall_v, status, tol, max_iter)
         end if
         if (status /= 0) then
            out%info = status
            allocate(out%u(0, 0), out%v(0, 0))
            return
         end if
         if (.not. want_vectors) then
            allocate(out%u(0, 0), out%v(0, 0))
            out%info = 0
            return
         end if
         out%u = tall_v
         if (thin_vectors) then
            out%v = tall_u
         else
            allocate(out%v(n, n))
            out%v = 0.0_dp
            out%v(:, 1:r) = tall_u
            if (n > r) call complete_basis(out%v, r)
         end if
      end if
      out%info = 0
   end subroutine svd2

   pure subroutine bidiagonal_qr_tall(a, singular_values, u, v, info, tol, max_iter)
      !! Compute a tall-matrix SVD by Householder bidiagonalization followed by implicit QR iteration.
      real(dp), intent(in) :: a(:, :) !! Tall or square m-by-n matrix with m >= n.
      real(dp), allocatable, intent(out) :: singular_values(:) !! Singular values in nonincreasing order.
      real(dp), allocatable, intent(out) :: u(:, :) !! Thin left singular vectors with shape m-by-n.
      real(dp), allocatable, intent(out) :: v(:, :) !! Full n-by-n right singular-vector matrix.
      integer, intent(out) :: info !! Status: zero converged, one invalid shape, two QR iteration limit reached.
      real(dp), intent(in), optional :: tol !! Relative bidiagonal deflation tolerance; default is machine epsilon.
      integer, intent(in), optional :: max_iter !! Maximum implicit-QR steps per singular value; default is 100.
      real(dp), allocatable :: rv1(:), work(:, :)
      real(dp) :: anorm, c, f, g, h, relative_tol, s_rotation, scale, x, y, z
      integer :: i, iteration, j, k, l, limit, m, n, nm
      logical :: cancellation

      m = size(a, 1)
      n = size(a, 2)
      if (m < n .or. n < 1) then
         allocate(singular_values(0), u(0, 0), v(0, 0))
         info = 1
         return
      end if
      work = a
      allocate(singular_values(n), rv1(n), v(n, n))
      singular_values = 0.0_dp
      rv1 = 0.0_dp
      v = 0.0_dp
      relative_tol = epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) relative_tol = tol
      end if
      limit = 100
      if (present(max_iter)) then
         if (max_iter > 0) limit = max_iter
      end if

      g = 0.0_dp
      scale = 0.0_dp
      anorm = 0.0_dp
      do i = 1, n
         l = i + 1
         rv1(i) = scale*g
         g = 0.0_dp
         s_rotation = 0.0_dp
         scale = sum(abs(work(i:m, i)))
         if (scale > 0.0_dp) then
            work(i:m, i) = work(i:m, i)/scale
            s_rotation = dot_product(work(i:m, i), work(i:m, i))
            f = work(i, i)
            g = -sign(sqrt(s_rotation), f)
            h = f*g - s_rotation
            work(i, i) = f - g
            if (i < n) then
               do j = l, n
                  s_rotation = dot_product(work(i:m, i), work(i:m, j))
                  f = s_rotation/h
                  work(i:m, j) = work(i:m, j) + f*work(i:m, i)
               end do
            end if
            work(i:m, i) = scale*work(i:m, i)
         end if
         singular_values(i) = scale*g

         g = 0.0_dp
         s_rotation = 0.0_dp
         scale = 0.0_dp
         if (i < n) then
            scale = sum(abs(work(i, l:n)))
            if (scale > 0.0_dp) then
               work(i, l:n) = work(i, l:n)/scale
               s_rotation = dot_product(work(i, l:n), work(i, l:n))
               f = work(i, l)
               g = -sign(sqrt(s_rotation), f)
               h = f*g - s_rotation
               work(i, l) = f - g
               rv1(l:n) = work(i, l:n)/h
               if (i < m) then
                  do j = l, m
                     s_rotation = dot_product(work(j, l:n), work(i, l:n))
                     work(j, l:n) = work(j, l:n) + s_rotation*rv1(l:n)
                  end do
               end if
               work(i, l:n) = scale*work(i, l:n)
            end if
         end if
         anorm = max(anorm, abs(singular_values(i)) + abs(rv1(i)))
      end do

      g = 0.0_dp
      l = n + 1
      do i = n, 1, -1
         if (i < n) then
            if (abs(g) > tiny(1.0_dp)) then
               v(l:n, i) = (work(i, l:n)/work(i, l))/g
               do j = l, n
                  s_rotation = dot_product(work(i, l:n), v(l:n, j))
                  v(l:n, j) = v(l:n, j) + s_rotation*v(l:n, i)
               end do
            end if
            v(i, l:n) = 0.0_dp
            v(l:n, i) = 0.0_dp
         end if
         v(i, i) = 1.0_dp
         g = rv1(i)
         l = i
      end do

      do i = n, 1, -1
         l = i + 1
         g = singular_values(i)
         if (i < n) work(i, l:n) = 0.0_dp
         if (abs(g) > tiny(1.0_dp)) then
            g = 1.0_dp/g
            if (i < n) then
               do j = l, n
                  s_rotation = dot_product(work(l:m, i), work(l:m, j))
                  f = (s_rotation/work(i, i))*g
                  work(i:m, j) = work(i:m, j) + f*work(i:m, i)
               end do
            end if
            work(i:m, i) = g*work(i:m, i)
         else
            work(i:m, i) = 0.0_dp
         end if
         work(i, i) = work(i, i) + 1.0_dp
      end do

      info = 0
      do k = n, 1, -1
         do iteration = 1, limit
            cancellation = .true.
            l = k
            do while (l >= 1)
               if (abs(rv1(l)) <= relative_tol*max(1.0_dp, anorm)) then
                  cancellation = .false.
                  exit
               end if
               if (l == 1) exit
               nm = l - 1
               if (abs(singular_values(nm)) <= relative_tol*max(1.0_dp, anorm)) exit
               l = l - 1
            end do

            if (cancellation) then
               c = 0.0_dp
               s_rotation = 1.0_dp
               nm = l - 1
               do i = l, k
                  f = s_rotation*rv1(i)
                  rv1(i) = c*rv1(i)
                  if (abs(f) <= relative_tol*max(1.0_dp, anorm)) exit
                  g = singular_values(i)
                  h = hypot(f, g)
                  singular_values(i) = h
                  if (h <= tiny(1.0_dp)) cycle
                  h = 1.0_dp/h
                  c = g*h
                  s_rotation = -f*h
                  do j = 1, m
                     y = work(j, nm)
                     z = work(j, i)
                     work(j, nm) = y*c + z*s_rotation
                     work(j, i) = z*c - y*s_rotation
                  end do
               end do
            end if

            z = singular_values(k)
            if (l == k) then
               if (z < 0.0_dp) then
                  singular_values(k) = -z
                  v(:, k) = -v(:, k)
               end if
               exit
            end if
            if (iteration == limit) then
               info = 2
               allocate(u(0, 0))
               return
            end if

            nm = k - 1
            x = singular_values(l)
            y = singular_values(nm)
            g = rv1(nm)
            h = rv1(k)
            if (abs(h*y) <= tiny(1.0_dp)) then
               f = 0.0_dp
            else
               f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0_dp*h*y)
            end if
            g = hypot(f, 1.0_dp)
            if (abs(x) <= tiny(1.0_dp)) then
               f = 0.0_dp
            else
               f = ((x - z)*(x + z) + h*(y/(f + sign(g, f)) - h))/x
            end if
            c = 1.0_dp
            s_rotation = 1.0_dp
            do j = l, nm
               i = j + 1
               g = rv1(i)
               y = singular_values(i)
               h = s_rotation*g
               g = c*g
               z = hypot(f, h)
               rv1(j) = z
               if (z > tiny(1.0_dp)) then
                  c = f/z
                  s_rotation = h/z
               else
                  c = 1.0_dp
                  s_rotation = 0.0_dp
               end if
               f = x*c + g*s_rotation
               g = g*c - x*s_rotation
               h = y*s_rotation
               y = y*c
               do i = 1, n
                  x = v(i, j)
                  z = v(i, j + 1)
                  v(i, j) = x*c + z*s_rotation
                  v(i, j + 1) = z*c - x*s_rotation
               end do
               z = hypot(f, h)
               singular_values(j) = z
               if (z > tiny(1.0_dp)) then
                  z = 1.0_dp/z
                  c = f*z
                  s_rotation = h*z
               else
                  c = 1.0_dp
                  s_rotation = 0.0_dp
               end if
               f = c*g + s_rotation*y
               x = c*y - s_rotation*g
               do i = 1, m
                  y = work(i, j)
                  z = work(i, j + 1)
                  work(i, j) = y*c + z*s_rotation
                  work(i, j + 1) = z*c - y*s_rotation
               end do
            end do
            rv1(l) = 0.0_dp
            rv1(k) = f
            singular_values(k) = x
         end do
      end do

      allocate(u(m, n))
      u = work(:, 1:n)
      call sort_svd_columns(singular_values, u, v)
   end subroutine bidiagonal_qr_tall

   pure subroutine jacobi_svd_tall(a, singular_values, u, v, info, tol, max_iter)
      !! Compute a Jacobi SVD after column-pivoted Householder QR preconditioning, matching Eigen's default strategy.
      real(dp), intent(in) :: a(:, :) !! Tall or square m-by-n matrix with m >= n.
      real(dp), allocatable, intent(out) :: singular_values(:) !! Singular values in nonincreasing order.
      real(dp), allocatable, intent(out) :: u(:, :) !! Thin left singular vectors with shape m-by-n.
      real(dp), allocatable, intent(out) :: v(:, :) !! Full n-by-n right singular-vector matrix.
      integer, intent(out) :: info !! Status: zero converged, one invalid shape, two Jacobi sweep limit reached.
      real(dp), intent(in), optional :: tol !! Relative two-sided Jacobi convergence tolerance.
      integer, intent(in), optional :: max_iter !! Maximum number of two-sided Jacobi sweeps; default 100.
      integer, allocatable :: perm(:)
      real(dp), allocatable :: qthin(:, :), r(:, :), square_u(:, :), square_v(:, :)
      integer :: i, n, status

      n = size(a, 2)
      if (size(a, 1) < n .or. n < 1) then
         allocate(singular_values(0), u(0, 0), v(0, 0))
         info = 1
         return
      end if

      call column_pivot_qr_precondition(a, qthin, r, perm, status)
      if (status /= 0) then
         call one_sided_jacobi_tall(a, singular_values, u, v, info, tol, max_iter)
         return
      end if
      call two_sided_jacobi_square(r, singular_values, square_u, square_v, status, tol, max_iter)
      if (status /= 0) then
         call one_sided_jacobi_tall(a, singular_values, u, v, info, tol, max_iter)
         return
      end if

      u = matmul(qthin, square_u)
      allocate(v(n, n))
      v = 0.0_dp
      do i = 1, n
         v(perm(i), :) = square_v(i, :)
      end do
      info = 0
   end subroutine jacobi_svd_tall

   pure subroutine column_pivot_qr_precondition(a, qthin, r, perm, info)
      !! Form a column-pivoted Householder QR preconditioner before the square two-sided Jacobi sweeps.
      real(dp), intent(in) :: a(:, :) !! Tall or square m-by-n matrix with m >= n.
      real(dp), allocatable, intent(out) :: qthin(:, :) !! First n columns of the orthogonal Q factor.
      real(dp), allocatable, intent(out) :: r(:, :) !! Leading n-by-n upper-triangular R factor.
      integer, allocatable, intent(out) :: perm(:) !! New-index to original-index column permutation.
      integer, intent(out) :: info !! Status: zero success, one invalid shape.
      real(dp), allocatable :: beta(:), reflectors(:, :), temporary(:), vector(:), work(:, :)
      real(dp) :: alpha, denominator, norm_value, projection, threshold
      integer :: i, j, k, m, n, pivot_col, temporary_index

      m = size(a, 1)
      n = size(a, 2)
      if (m < n .or. n < 1) then
         allocate(qthin(0, 0), r(0, 0), perm(0))
         info = 1
         return
      end if
      work = a
      allocate(beta(n), reflectors(m, n), perm(n))
      beta = 0.0_dp
      reflectors = 0.0_dp
      perm = [(i, i=1, n)]
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(m, n), dp)

      do k = 1, n
         pivot_col = k
         norm_value = -1.0_dp
         do j = k, n
            denominator = dot_product(work(k:m, j), work(k:m, j))
            if (denominator > norm_value) then
               norm_value = denominator
               pivot_col = j
            end if
         end do
         if (pivot_col /= k) then
            temporary = work(:, k)
            work(:, k) = work(:, pivot_col)
            work(:, pivot_col) = temporary
            temporary_index = perm(k)
            perm(k) = perm(pivot_col)
            perm(pivot_col) = temporary_index
         end if

         vector = work(k:m, k)
         norm_value = sqrt(max(0.0_dp, dot_product(vector, vector)))
         if (norm_value <= threshold) cycle
         alpha = -sign(norm_value, vector(1))
         vector(1) = vector(1) - alpha
         denominator = dot_product(vector, vector)
         if (denominator <= tiny(1.0_dp)) cycle
         beta(k) = 2.0_dp/denominator
         reflectors(k:m, k) = vector
         do j = k, n
            projection = beta(k)*dot_product(vector, work(k:m, j))
            work(k:m, j) = work(k:m, j) - vector*projection
         end do
         work(k, k) = alpha
         if (k < m) work(k + 1:m, k) = 0.0_dp
      end do

      r = work(1:n, 1:n)
      allocate(qthin(m, n))
      qthin = 0.0_dp
      do i = 1, n
         qthin(i, i) = 1.0_dp
      end do
      do k = n, 1, -1
         if (abs(beta(k)) <= tiny(1.0_dp)) cycle
         do j = 1, n
            projection = beta(k)*dot_product(reflectors(k:m, k), qthin(k:m, j))
            qthin(k:m, j) = qthin(k:m, j) - reflectors(k:m, k)*projection
         end do
      end do
      info = 0
   end subroutine column_pivot_qr_precondition

   pure subroutine two_sided_jacobi_square(a, singular_values, u, v, info, tol, max_iter)
      !! Diagonalize a real square matrix by two-sided Jacobi SVD sweeps.
      real(dp), intent(in) :: a(:, :) !! Real square matrix after rectangular QR preconditioning.
      real(dp), allocatable, intent(out) :: singular_values(:) !! Singular values in nonincreasing order.
      real(dp), allocatable, intent(out) :: u(:, :) !! Orthogonal left singular-vector matrix.
      real(dp), allocatable, intent(out) :: v(:, :) !! Orthogonal right singular-vector matrix.
      integer, intent(out) :: info !! Status: zero converged, one invalid shape, two sweep limit reached.
      real(dp), intent(in), optional :: tol !! Relative off-diagonal convergence tolerance.
      integer, intent(in), optional :: max_iter !! Maximum complete pair sweeps; default 100.
      real(dp), allocatable :: column_p(:), column_q(:), row_p(:), row_q(:), work(:, :)
      real(dp) :: local_u(2, 2), local_v(2, 2), pair_tolerance, scale
      integer :: i, limit, n, p, q, sweep
      logical :: changed

      n = size(a, 1)
      if (size(a, 2) /= n .or. n < 1) then
         allocate(singular_values(0), u(0, 0), v(0, 0))
         info = 1
         return
      end if
      work = a
      allocate(u(n, n), v(n, n))
      u = 0.0_dp
      v = 0.0_dp
      do i = 1, n
         u(i, i) = 1.0_dp
         v(i, i) = 1.0_dp
      end do
      pair_tolerance = 10.0_dp*epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) pair_tolerance = tol
      end if
      limit = 100
      if (present(max_iter)) then
         if (max_iter > 0) limit = max_iter
      end if

      info = 2
      do sweep = 1, limit
         changed = .false.
         do p = 1, n - 1
            do q = p + 1, n
               scale = max(1.0_dp, hypot(work(p, p), work(q, q)))
               if (hypot(work(p, q), work(q, p)) <= pair_tolerance*scale) cycle
               call svd_rotation_2x2(work(p, p), work(p, q), work(q, p), work(q, q), local_u, local_v)

               row_p = work(p, :)
               row_q = work(q, :)
               work(p, :) = local_u(1, 1)*row_p + local_u(2, 1)*row_q
               work(q, :) = local_u(1, 2)*row_p + local_u(2, 2)*row_q
               column_p = work(:, p)
               column_q = work(:, q)
               work(:, p) = local_v(1, 1)*column_p + local_v(2, 1)*column_q
               work(:, q) = local_v(1, 2)*column_p + local_v(2, 2)*column_q

               column_p = u(:, p)
               column_q = u(:, q)
               u(:, p) = local_u(1, 1)*column_p + local_u(2, 1)*column_q
               u(:, q) = local_u(1, 2)*column_p + local_u(2, 2)*column_q
               column_p = v(:, p)
               column_q = v(:, q)
               v(:, p) = local_v(1, 1)*column_p + local_v(2, 1)*column_q
               v(:, q) = local_v(1, 2)*column_p + local_v(2, 2)*column_q
               changed = .true.
            end do
         end do
         if (.not. changed) then
            info = 0
            exit
         end if
      end do
      if (info /= 0) then
         allocate(singular_values(0))
         return
      end if

      allocate(singular_values(n))
      do i = 1, n
         if (work(i, i) < 0.0_dp) u(:, i) = -u(:, i)
         singular_values(i) = abs(work(i, i))
      end do
      call sort_svd_columns(singular_values, u, v)
   end subroutine two_sided_jacobi_square

   pure subroutine svd_rotation_2x2(a11, a12, a21, a22, u, v)
      !! Compute orthogonal 2-by-2 left/right rotations whose product diagonalizes a real 2-by-2 block.
      real(dp), intent(in) :: a11 !! Upper-left block entry.
      real(dp), intent(in) :: a12 !! Upper-right block entry.
      real(dp), intent(in) :: a21 !! Lower-left block entry.
      real(dp), intent(in) :: a22 !! Lower-right block entry.
      real(dp), intent(out) :: u(2, 2) !! Left orthogonal rotation.
      real(dp), intent(out) :: v(2, 2) !! Right orthogonal rotation.
      real(dp) :: b11, b12, b21, b22, c, cosine, dot_value, norm_value, s_value, sine, theta

      b11 = a11*a11 + a21*a21
      b12 = a11*a12 + a21*a22
      b22 = a12*a12 + a22*a22
      theta = 0.5_dp*atan2(2.0_dp*b12, b11 - b22)
      cosine = cos(theta)
      sine = sin(theta)
      v(1, 1) = cosine
      v(1, 2) = -sine
      v(2, 1) = sine
      v(2, 2) = cosine

      b11 = a11*cosine + a12*sine
      b21 = a21*cosine + a22*sine
      b12 = -a11*sine + a12*cosine
      b22 = -a21*sine + a22*cosine
      norm_value = hypot(b11, b21)
      if (norm_value > tiny(1.0_dp)) then
         u(1, 1) = b11/norm_value
         u(2, 1) = b21/norm_value
      else
         u(1, 1) = 1.0_dp
         u(2, 1) = 0.0_dp
      end if
      dot_value = u(1, 1)*b12 + u(2, 1)*b22
      b12 = b12 - dot_value*u(1, 1)
      b22 = b22 - dot_value*u(2, 1)
      norm_value = hypot(b12, b22)
      if (norm_value > tiny(1.0_dp)) then
         u(1, 2) = b12/norm_value
         u(2, 2) = b22/norm_value
      else
         u(1, 2) = -u(2, 1)
         u(2, 2) = u(1, 1)
      end if
      s_value = u(1, 2)*(-a11*sine + a12*cosine) + u(2, 2)*(-a21*sine + a22*cosine)
      if (s_value < 0.0_dp) u(:, 2) = -u(:, 2)
      c = u(1, 1)*(a11*cosine + a12*sine) + u(2, 1)*(a21*cosine + a22*sine)
      if (c < 0.0_dp) u(:, 1) = -u(:, 1)
   end subroutine svd_rotation_2x2

   pure subroutine one_sided_jacobi_tall(a, singular_values, u, v, info, tol, max_iter)
      !! Direct one-sided Jacobi SVD for a matrix with at least as many rows as columns.
      real(dp), intent(in) :: a(:, :) !! Tall or square m-by-n matrix with m >= n.
      real(dp), allocatable, intent(out) :: singular_values(:) !! Singular values in nonincreasing order.
      real(dp), allocatable, intent(out) :: u(:, :) !! Thin left singular vectors with shape m-by-n.
      real(dp), allocatable, intent(out) :: v(:, :) !! Full n-by-n right singular-vector matrix.
      integer, intent(out) :: info !! Status: zero converged, one invalid shape, two sweep limit reached.
      real(dp), intent(in), optional :: tol !! Relative column-orthogonality tolerance; default is 10 machine epsilons.
      integer, intent(in), optional :: max_iter !! Maximum number of complete Jacobi sweeps; default 100.
      real(dp), allocatable :: work(:, :), column_p(:), column_q(:), v_p(:), v_q(:)
      real(dp) :: alpha, beta, c, gamma, pair_tolerance, s_rotation, scale, t, tau
      integer :: j, limit, n, p, q, sweep
      logical :: changed

      n = size(a, 2)
      if (size(a, 1) < n .or. n < 1) then
         allocate(singular_values(0), u(0, 0), v(0, 0))
         info = 1
         return
      end if
      work = a
      allocate(v(n, n))
      v = 0.0_dp
      do j = 1, n
         v(j, j) = 1.0_dp
      end do
      pair_tolerance = 10.0_dp*epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) pair_tolerance = tol
      end if
      limit = 100
      if (present(max_iter)) then
         if (max_iter > 0) limit = max_iter
      end if

      info = 2
      do sweep = 1, limit
         changed = .false.
         do p = 1, n - 1
            do q = p + 1, n
               alpha = dot_product(work(:, p), work(:, p))
               beta = dot_product(work(:, q), work(:, q))
               gamma = dot_product(work(:, p), work(:, q))
               scale = sqrt(max(0.0_dp, alpha*beta))
               if (scale <= tiny(1.0_dp)) cycle
               if (abs(gamma) <= pair_tolerance*scale) cycle
               tau = (beta - alpha)/(2.0_dp*gamma)
               t = sign(1.0_dp, tau)/(abs(tau) + sqrt(1.0_dp + tau*tau))
               c = 1.0_dp/sqrt(1.0_dp + t*t)
               s_rotation = c*t
               column_p = work(:, p)
               column_q = work(:, q)
               work(:, p) = c*column_p - s_rotation*column_q
               work(:, q) = s_rotation*column_p + c*column_q
               v_p = v(:, p)
               v_q = v(:, q)
               v(:, p) = c*v_p - s_rotation*v_q
               v(:, q) = s_rotation*v_p + c*v_q
               changed = .true.
            end do
         end do
         if (.not. changed) then
            info = 0
            exit
         end if
      end do
      if (info /= 0) then
         allocate(singular_values(0), u(0, 0))
         return
      end if

      allocate(singular_values(n), u(size(a, 1), n))
      do j = 1, n
         singular_values(j) = sqrt(max(0.0_dp, dot_product(work(:, j), work(:, j))))
      end do
      call sort_svd_columns(singular_values, work, v)
      u = 0.0_dp
      scale = max(1.0_dp, maxval(singular_values))
      do j = 1, n
         if (singular_values(j) > epsilon(1.0_dp)*scale*real(max(size(a, 1), n), dp)) then
            u(:, j) = work(:, j)/singular_values(j)
         end if
         call normalize_with_completion(u, j)
      end do
   end subroutine one_sided_jacobi_tall

   pure subroutine sort_svd_columns(singular_values, work, v)
      !! Sort singular values descending and apply the same swaps to transformed and right-vector columns.
      real(dp), intent(inout) :: singular_values(:) !! Singular values to sort in nonincreasing order.
      real(dp), intent(inout) :: work(:, :) !! Rotated A columns paired with the singular values.
      real(dp), intent(inout) :: v(:, :) !! Right singular-vector columns paired with the singular values.
      real(dp), allocatable :: temporary(:)
      real(dp) :: value
      integer :: i, j, selected

      do i = 1, size(singular_values) - 1
         selected = i
         do j = i + 1, size(singular_values)
            if (singular_values(j) > singular_values(selected)) selected = j
         end do
         if (selected /= i) then
            value = singular_values(i)
            singular_values(i) = singular_values(selected)
            singular_values(selected) = value
            temporary = work(:, i)
            work(:, i) = work(:, selected)
            work(:, selected) = temporary
            temporary = v(:, i)
            v(:, i) = v(:, selected)
            v(:, selected) = temporary
         end if
      end do
   end subroutine sort_svd_columns

   pure subroutine largest_offdiagonal(a, p, q, value)
      !! Locate the largest-magnitude strict upper-triangular entry.
      real(dp), intent(in) :: a(:, :) !! Square symmetric work matrix.
      integer, intent(out) :: p !! Row index of the selected off-diagonal entry.
      integer, intent(out) :: q !! Column index of the selected off-diagonal entry.
      real(dp), intent(out) :: value !! Signed selected matrix entry.
      real(dp) :: largest
      integer :: i, j, n

      n = size(a, 1)
      p = 1
      q = min(2, n)
      largest = -1.0_dp
      value = 0.0_dp
      do j = 2, n
         do i = 1, j - 1
            if (abs(a(i, j)) > largest) then
               largest = abs(a(i, j))
               value = a(i, j)
               p = i
               q = j
            end if
         end do
      end do
   end subroutine largest_offdiagonal

   pure subroutine sort_real_eigenpairs(values, vectors, ascending)
      !! Sort real eigenvalues and apply the same permutation to eigenvector columns.
      real(dp), intent(inout) :: values(:) !! Eigenvalues to sort in place.
      real(dp), intent(inout) :: vectors(:, :) !! Eigenvector columns paired with values.
      logical, intent(in) :: ascending !! True for ascending, false for descending order.
      real(dp), allocatable :: tmp(:)
      real(dp) :: value_tmp
      integer :: i, j, selected

      do i = 1, size(values) - 1
         selected = i
         do j = i + 1, size(values)
            if ((ascending .and. values(j) < values(selected)) .or. &
                (.not. ascending .and. values(j) > values(selected))) selected = j
         end do
         if (selected /= i) then
            value_tmp = values(i)
            values(i) = values(selected)
            values(selected) = value_tmp
            tmp = vectors(:, i)
            vectors(:, i) = vectors(:, selected)
            vectors(:, selected) = tmp
         end if
      end do
   end subroutine sort_real_eigenpairs

   pure function wilkinson_shift_complex(a, b, c, d) result(mu)
      !! Return the eigenvalue of a trailing 2-by-2 block closest to d.
      complex(dp), intent(in) :: a !! Upper-left entry of the trailing 2-by-2 block.
      complex(dp), intent(in) :: b !! Upper-right entry of the trailing 2-by-2 block.
      complex(dp), intent(in) :: c !! Lower-left entry of the trailing 2-by-2 block.
      complex(dp), intent(in) :: d !! Lower-right entry of the trailing 2-by-2 block.
      complex(dp) :: mu, disc, root1, root2, trace_half

      trace_half = 0.5_dp*(a + d)
      disc = sqrt(0.25_dp*(a - d)*(a - d) + b*c)
      root1 = trace_half + disc
      root2 = trace_half - disc
      if (abs(root1 - d) <= abs(root2 - d)) then
         mu = root1
      else
         mu = root2
      end if
   end function wilkinson_shift_complex

   pure subroutine complex_qr(a, q, r, info)
      !! Rank-tolerant complex Householder QR factorization used by the nonsymmetric eigen iteration.
      complex(dp), intent(in) :: a(:, :) !! Complex square matrix to factor.
      complex(dp), allocatable, intent(out) :: q(:, :) !! Unitary factor.
      complex(dp), allocatable, intent(out) :: r(:, :) !! Upper-triangular factor.
      integer, intent(out) :: info !! Status: zero success, one invalid nonsquare shape.
      complex(dp), allocatable :: reflector(:)
      complex(dp) :: alpha, phase, projection
      real(dp) :: beta, norm_value, reflector_norm, threshold
      integer :: i, j, k, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(q(0, 0), r(0, 0))
         info = 1
         return
      end if
      r = a
      allocate(q(n, n))
      q = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do i = 1, n
         q(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
      end do
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(1, n), dp)

      do k = 1, n
         reflector = r(k:n, k)
         norm_value = sqrt(sum(abs(reflector)**2))
         if (norm_value <= threshold) cycle
         if (abs(reflector(1)) > threshold) then
            phase = reflector(1)/abs(reflector(1))
         else
            phase = cmplx(1.0_dp, 0.0_dp, kind=dp)
         end if
         alpha = -phase*norm_value
         reflector(1) = reflector(1) - alpha
         reflector_norm = real(dot_product(reflector, reflector), dp)
         if (reflector_norm <= tiny(1.0_dp)) cycle
         beta = 2.0_dp/reflector_norm

         do j = k, n
            projection = beta*dot_product(reflector, r(k:n, j))
            r(k:n, j) = r(k:n, j) - reflector*projection
         end do
         if (k < n) r(k + 1:n, k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
         r(k, k) = alpha
         do i = 1, n
            projection = beta*sum(q(i, k:n)*reflector)
            q(i, k:n) = q(i, k:n) - projection*conjg(reflector)
         end do
      end do
      info = 0
   end subroutine complex_qr

   pure subroutine inverse_iteration_vectors(a, values, vectors)
      !! Estimate right eigenvectors for supplied eigenvalues by shifted inverse iteration.
      real(dp), intent(in) :: a(:, :) !! Original real square matrix.
      complex(dp), intent(in) :: values(:) !! Approximate eigenvalues.
      complex(dp), allocatable, intent(out) :: vectors(:, :) !! Normalized right eigenvectors.
      complex(dp), allocatable :: matrix(:, :), rhs(:), next(:)
      complex(dp) :: shift
      real(dp) :: delta, norm_value, scale
      integer :: i, iteration, n, status

      n = size(a, 1)
      allocate(vectors(n, n))
      vectors = cmplx(0.0_dp, 0.0_dp, kind=dp)
      scale = max(1.0_dp, maxval(abs(a)))
      delta = sqrt(epsilon(1.0_dp))*scale
      do i = 1, n
         allocate(rhs(n))
         rhs = cmplx(1.0_dp, 0.0_dp, kind=dp)
         rhs = rhs/sqrt(real(n, dp))
         shift = values(i) + cmplx(delta, delta*real(i, dp)/real(max(1, n), dp), kind=dp)
         do iteration = 1, 8
            matrix = cmplx(a, 0.0_dp, kind=dp)
            matrix = matrix - shift*complex_identity(n)
            call complex_linear_solve(matrix, rhs, next, status)
            if (status /= 0) exit
            norm_value = sqrt(max(0.0_dp, real(dot_product(next, next), dp)))
            if (norm_value <= tiny(1.0_dp)) exit
            rhs = next/norm_value
         end do
         norm_value = sqrt(max(0.0_dp, real(dot_product(rhs, rhs), dp)))
         if (norm_value > tiny(1.0_dp)) vectors(:, i) = rhs/norm_value
         deallocate(rhs)
      end do
   end subroutine inverse_iteration_vectors

   pure subroutine complex_linear_solve(a, b, x, info)
      !! Solve one complex square linear system with partial row pivoting.
      complex(dp), intent(in) :: a(:, :) !! Complex square coefficient matrix.
      complex(dp), intent(in) :: b(:) !! Complex RHS vector.
      complex(dp), allocatable, intent(out) :: x(:) !! Solution vector.
      integer, intent(out) :: info !! Status: zero success, two numerical singularity.
      complex(dp), allocatable :: u(:, :), rhs(:), tmp_row(:)
      complex(dp) :: factor, tmp_value
      real(dp) :: threshold
      integer :: i, k, n, p

      n = size(a, 1)
      u = a
      rhs = b
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(1, n), dp)
      do k = 1, n - 1
         p = k - 1 + maxloc(abs(u(k:n, k)), dim=1)
         if (abs(u(p, k)) <= threshold) then
            allocate(x(0))
            info = 2
            return
         end if
         if (p /= k) then
            tmp_row = u(k, :)
            u(k, :) = u(p, :)
            u(p, :) = tmp_row
            tmp_value = rhs(k)
            rhs(k) = rhs(p)
            rhs(p) = tmp_value
         end if
         do i = k + 1, n
            factor = u(i, k)/u(k, k)
            u(i, k:n) = u(i, k:n) - factor*u(k, k:n)
            rhs(i) = rhs(i) - factor*rhs(k)
         end do
      end do
      if (n > 0 .and. abs(u(n, n)) <= threshold) then
         allocate(x(0))
         info = 2
         return
      end if
      allocate(x(n))
      x = rhs
      do i = n, 1, -1
         if (i < n) x(i) = x(i) - dot_product(conjg(u(i, i + 1:n)), x(i + 1:n))
         x(i) = x(i)/u(i, i)
      end do
      info = 0
   end subroutine complex_linear_solve

   pure function complex_identity(n) result(identity)
      !! Construct an n-by-n complex identity matrix.
      integer, intent(in) :: n !! Requested nonnegative identity order.
      complex(dp), allocatable :: identity(:, :)
      integer :: i

      allocate(identity(n, n))
      identity = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do i = 1, n
         identity(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
      end do
   end function complex_identity

   pure subroutine normalize_with_completion(matrix, column)
      !! Normalize one column, replacing a near-zero candidate by an orthogonal canonical vector.
      real(dp), intent(inout) :: matrix(:, :) !! Basis matrix whose earlier columns are treated as fixed.
      integer, intent(in) :: column !! One-based column to normalize/complete.
      real(dp) :: norm_value
      integer :: candidate, j

      do j = 1, column - 1
         matrix(:, column) = matrix(:, column) - dot_product(matrix(:, j), matrix(:, column))*matrix(:, j)
      end do
      norm_value = sqrt(max(0.0_dp, dot_product(matrix(:, column), matrix(:, column))))
      if (norm_value <= sqrt(epsilon(1.0_dp))) then
         do candidate = 1, size(matrix, 1)
            matrix(:, column) = 0.0_dp
            matrix(candidate, column) = 1.0_dp
            do j = 1, column - 1
               matrix(:, column) = matrix(:, column) - dot_product(matrix(:, j), matrix(:, column))*matrix(:, j)
            end do
            norm_value = sqrt(max(0.0_dp, dot_product(matrix(:, column), matrix(:, column))))
            if (norm_value > sqrt(epsilon(1.0_dp))) exit
         end do
      end if
      if (norm_value > tiny(1.0_dp)) matrix(:, column) = matrix(:, column)/norm_value
   end subroutine normalize_with_completion

   pure subroutine complete_basis(matrix, filled)
      !! Fill columns after filled with deterministic orthonormal completion vectors.
      real(dp), intent(inout) :: matrix(:, :) !! Matrix whose leading filled columns already form a basis.
      integer, intent(in) :: filled !! Number of leading columns considered initialized.
      integer :: j

      do j = filled + 1, size(matrix, 2)
         matrix(:, j) = 0.0_dp
         call normalize_with_completion(matrix, j)
      end do
   end subroutine complete_basis


end module sanic_spectral
