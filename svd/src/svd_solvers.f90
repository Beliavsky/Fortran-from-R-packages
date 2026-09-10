! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of PROPACK/nuTRLan-facing APIs from R package svd 0.5.8.
module svd_solvers
   use svd_kinds, only : dp
   use rspectra, only : dense_operator, eigs_opts, eigs_sym, eigs_sym_result
   use rspectra, only : linear_operator, make_dense_operator, svds, svds_opts, svds_result
   use r_linalg, only : complex_thin_svd
   implicit none
   private

   integer, parameter, public :: solver_success = 0
   integer, parameter, public :: solver_invalid_input = -2001
   integer, parameter, public :: solver_invalid_warm_start = -2002

   type, public :: svd_options
      !! Common translation of the named R option-list entries used by PROPACK and nuTRLan wrappers.
      integer :: kmax = 0
      integer :: dim = 0
      integer :: p = 0
      integer :: maxiter = 0
      integer :: verbose = 0
      real(dp) :: tol = -1.0_dp
   end type svd_options

   type, public :: propack_svd_result
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: v(:, :)
      integer :: nconv = 0
      integer :: niter = 0
      integer :: nops = 0
      integer :: info = 0
   end type propack_svd_result

   type, public :: trlan_svd_result
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: u(:, :)
      integer :: nconv = 0
      integer :: niter = 0
      integer :: nops = 0
      integer :: info = 0
   end type trlan_svd_result

   type, public :: ztrlan_svd_result
      real(dp), allocatable :: d(:)
      complex(dp), allocatable :: u(:, :)
      integer :: nconv = 0
      integer :: niter = 0
      integer :: nops = 0
      integer :: info = 0
   end type ztrlan_svd_result

   type, public :: trlan_eigen_result
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: u(:, :)
      integer :: nconv = 0
      integer :: niter = 0
      integer :: nops = 0
      integer :: info = 0
   end type trlan_eigen_result

   type, extends(linear_operator) :: left_normal_operator
      class(linear_operator), pointer :: base => null()
      real(dp), allocatable :: work(:)
   contains
      procedure :: prod => left_normal_prod
      procedure :: tprod => left_normal_tprod
   end type left_normal_operator

   public :: propack_svd_dense
   public :: propack_svd_operator
   public :: trlan_eigen_dense
   public :: trlan_eigen_operator
   public :: trlan_svd_dense
   public :: trlan_svd_operator
   public :: ztrlan_svd_complex

   interface propack_svd
      module procedure propack_svd_dense
      module procedure propack_svd_operator
   end interface propack_svd
   public :: propack_svd

   interface trlan_svd
      module procedure trlan_svd_dense
      module procedure trlan_svd_operator
   end interface trlan_svd
   public :: trlan_svd

   interface trlan_eigen
      module procedure trlan_eigen_dense
      module procedure trlan_eigen_operator
   end interface trlan_eigen
   public :: trlan_eigen

contains

   function propack_svd_dense(a, neig, options) result(out)
      !! Computes leading real singular triplets for a dense matrix using the shared RSpectra backend.
      real(dp), intent(in) :: a(:, :) !! Dense input matrix with shape `(m, n)`.
      integer, intent(in), optional :: neig !! Requested leading singular triplets; defaults to `min(m,n)` and is clipped there.
      type(svd_options), intent(in), optional :: options !! PROPACK controls for subspace size, tolerance, and restarts.
      type(propack_svd_result) :: out
      type(svds_opts) :: controls
      type(svds_result) :: result
      type(svd_options) :: local
      integer :: k, wd

      wd = min(size(a, 1), size(a, 2))
      k = wd
      if (present(neig)) k = min(neig, wd)
      if (wd < 1 .or. k < 1) then
         call empty_propack(out, size(a, 1), size(a, 2), solver_invalid_input)
         return
      end if

      local = svd_options()
      if (present(options)) local = options
      controls = propack_controls(local, k, wd)
      result = svds(a, k, nu=k, nv=k, opts=controls)
      call copy_propack_result(result, out)
   end function propack_svd_dense

   function propack_svd_operator(op, neig, options) result(out)
      !! Computes leading real singular triplets for a matrix-free operator using the shared iterative SVD backend.
      class(linear_operator), intent(inout) :: op !! Matrix-free operator providing both forward and transpose products.
      integer, intent(in), optional :: neig !! Leading triplet count; defaults to and is clipped at `min(nrow,ncol)`.
      type(svd_options), intent(in), optional :: options !! PROPACK-style controls mapped onto RSpectra/ARPACK controls.
      type(propack_svd_result) :: out
      real(dp), allocatable :: a(:, :)
      type(svds_opts) :: controls
      type(svds_result) :: result
      type(svd_options) :: local
      integer :: k, wd

      wd = min(op%nrow, op%ncol)
      k = wd
      if (present(neig)) k = min(neig, wd)
      if (wd < 1 .or. k < 1) then
         call empty_propack(out, op%nrow, op%ncol, solver_invalid_input)
         return
      end if

      local = svd_options()
      if (present(options)) local = options
      if (k == wd .or. wd < 3) then
         call materialize_operator(op, a)
         out = propack_svd_dense(a, k, local)
         out%nops = out%nops + op%ncol
         return
      end if

      controls = propack_controls(local, k, wd)
      result = svds(op, k, nu=k, nv=k, opts=controls)
      call copy_propack_result(result, out)
   end function propack_svd_operator

   function trlan_svd_dense(a, neig, options, lambda_start, u_start) result(out)
      !! Computes leading singular values and left vectors through the left normal operator, matching nuTRLan's R result shape.
      real(dp), intent(in) :: a(:, :) !! Dense real input matrix with shape `(m, n)`.
      integer, intent(in), optional :: neig !! Requested singular pairs; defaults to `min(m,n)` and is clipped there.
      type(svd_options), intent(in), optional :: options !! nuTRLan controls for Krylov size, tolerance, and product budget.
      real(dp), intent(in), optional :: lambda_start(:) !! Prior singular values paired with columns of `u_start`.
      real(dp), intent(in), optional :: u_start(:, :) !! Prior left vectors; one matched vector seeds the shared solver.
      type(trlan_svd_result) :: out
      type(dense_operator) :: op

      op = make_dense_operator(a)
      out = trlan_svd_operator(op, neig, options, lambda_start, u_start)
   end function trlan_svd_dense

   function trlan_svd_operator(op, neig, options, lambda_start, u_start) result(out)
      !! Computes leading singular values/left vectors for a matrix-free operator via `A*A^T` Lanczos iterations.
      class(linear_operator), target, intent(inout) :: op !! Matrix-free real operator representing `A`.
      integer, intent(in), optional :: neig !! Requested singular pairs; defaults to `min(nrow,ncol)` and is clipped there.
      type(svd_options), intent(in), optional :: options !! nuTRLan-style controls for subspace size, tolerance, and iterations.
      real(dp), intent(in), optional :: lambda_start(:) !! Prior singular values paired with `u_start`; values are not locked.
      real(dp), intent(in), optional :: u_start(:, :) !! Prior left vectors; the last match seeds the ARPACK residual.
      type(trlan_svd_result) :: out
      type(left_normal_operator) :: normal
      type(eigs_opts) :: controls
      type(eigs_sym_result) :: result
      type(svd_options) :: local
      real(dp), allocatable :: a(:, :), normal_matrix(:, :)
      integer :: k, wd, warm_status

      wd = min(op%nrow, op%ncol)
      k = wd
      if (present(neig)) k = min(neig, wd)
      if (wd < 1 .or. k < 1) then
         call empty_trlan_svd(out, op%nrow, solver_invalid_input)
         return
      end if

      local = svd_options()
      if (present(options)) local = options
      controls = trlan_controls(local, k, op%nrow, min(op%nrow, wd + 1))
      call apply_real_warm_start(controls, op%nrow, lambda_start, u_start, warm_status)
      if (warm_status /= solver_success) then
         call empty_trlan_svd(out, op%nrow, warm_status)
         return
      end if

      if (k == op%nrow) then
         call materialize_operator(op, a)
         normal_matrix = matmul(a, transpose(a))
         result = eigs_sym(normal_matrix, k, which='LA', opts=controls)
         result%nops = result%nops + op%ncol
      else
         normal%nrow = op%nrow
         normal%ncol = op%nrow
         normal%base => op
         allocate(normal%work(op%ncol))
         result = eigs_sym(normal, k, which='LA', opts=controls)
      end if
      call copy_trlan_svd_result(result, out)
   end function trlan_svd_operator

   function trlan_eigen_dense(a, neig, options, lambda_start, u_start) result(out)
      !! Computes leading eigenpairs of a dense real symmetric matrix using the shared RSpectra backend.
      real(dp), intent(in) :: a(:, :) !! Dense square matrix; the upstream algorithm assumes a symmetric operator.
      integer, intent(in), optional :: neig !! Requested leading eigenpairs; defaults to the matrix order and is clipped there.
      type(svd_options), intent(in), optional :: options !! nuTRLan-style subspace, tolerance, and iteration controls.
      real(dp), intent(in), optional :: lambda_start(:) !! Prior eigenvalues paired with `u_start`; values are not locked.
      real(dp), intent(in), optional :: u_start(:, :) !! Prior eigenvectors; the last matched vector seeds the iterative path.
      type(trlan_eigen_result) :: out
      type(eigs_opts) :: controls
      type(eigs_sym_result) :: result
      type(svd_options) :: local
      integer :: k, n, warm_status

      n = size(a, 1)
      if (size(a, 2) /= n .or. n < 1) then
         call empty_trlan_eigen(out, n, solver_invalid_input)
         return
      end if
      k = n
      if (present(neig)) k = min(neig, n)
      if (k < 1) then
         call empty_trlan_eigen(out, n, solver_invalid_input)
         return
      end if

      local = svd_options()
      if (present(options)) local = options
      controls = trlan_controls(local, k, n)
      call apply_real_warm_start(controls, n, lambda_start, u_start, warm_status)
      if (warm_status /= solver_success) then
         call empty_trlan_eigen(out, n, warm_status)
         return
      end if

      result = eigs_sym(a, k, which='LA', opts=controls)
      call copy_trlan_eigen_result(result, out)
   end function trlan_eigen_dense

   function trlan_eigen_operator(op, neig, options, lambda_start, u_start) result(out)
      !! Computes leading eigenpairs of a square real matrix-free symmetric operator.
      class(linear_operator), intent(inout) :: op !! Square operator whose forward product represents the symmetric matrix.
      integer, intent(in), optional :: neig !! Requested leading eigenpairs; defaults to the operator order and is clipped there.
      type(svd_options), intent(in), optional :: options !! nuTRLan-style subspace, tolerance, and iteration controls.
      real(dp), intent(in), optional :: lambda_start(:) !! Prior eigenvalues paired with `u_start`; values are not locked.
      real(dp), intent(in), optional :: u_start(:, :) !! Prior eigenvectors; the last matched vector seeds the iterative path.
      type(trlan_eigen_result) :: out
      real(dp), allocatable :: a(:, :)
      type(eigs_opts) :: controls
      type(eigs_sym_result) :: result
      type(svd_options) :: local
      integer :: k, n, warm_status

      n = op%nrow
      if (op%ncol /= n .or. n < 1) then
         call empty_trlan_eigen(out, n, solver_invalid_input)
         return
      end if
      k = n
      if (present(neig)) k = min(neig, n)
      if (k < 1) then
         call empty_trlan_eigen(out, n, solver_invalid_input)
         return
      end if

      local = svd_options()
      if (present(options)) local = options
      controls = trlan_controls(local, k, n)
      call apply_real_warm_start(controls, n, lambda_start, u_start, warm_status)
      if (warm_status /= solver_success) then
         call empty_trlan_eigen(out, n, warm_status)
         return
      end if

      if (k == n) then
         call materialize_operator(op, a)
         result = eigs_sym(a, k, which='LA', opts=controls)
         result%nops = result%nops + n
      else
         result = eigs_sym(op, k, which='LA', opts=controls)
      end if
      call copy_trlan_eigen_result(result, out)
   end function trlan_eigen_operator

   function ztrlan_svd_complex(a, neig, options, lambda_start, u_start) result(out)
      !! Computes the conventional complex SVD and returns leading singular values and left vectors.
      complex(dp), intent(in) :: a(:, :) !! Dense complex matrix with shape `(m, n)`.
      integer, intent(in), optional :: neig !! Requested leading singular pairs; defaults to `min(m,n)` and is clipped there.
      type(svd_options), intent(in), optional :: options !! Compatibility controls; dense LAPACK SVD ignores Lanczos settings.
      real(dp), intent(in), optional :: lambda_start(:) !! Prior singular values; dense SVD validates but does not reuse them.
      complex(dp), intent(in), optional :: u_start(:, :) !! Prior left vectors; dense SVD validates but does not reuse them.
      type(ztrlan_svd_result) :: out
      complex(dp), allocatable :: u(:, :), vt(:, :)
      real(dp), allocatable :: singular_values(:)
      integer :: info, k, wd, warm_status

      wd = min(size(a, 1), size(a, 2))
      k = wd
      if (present(neig)) k = min(neig, wd)
      if (wd < 1 .or. k < 1) then
         call empty_ztrlan_svd(out, size(a, 1), solver_invalid_input)
         return
      end if

      call validate_complex_warm_start(size(a, 1), lambda_start, u_start, warm_status)
      if (warm_status /= solver_success) then
         call empty_ztrlan_svd(out, size(a, 1), warm_status)
         return
      end if
      if (present(options)) then
         if (options%tol < -1.0_dp) then
            call empty_ztrlan_svd(out, size(a, 1), solver_invalid_input)
            return
         end if
      end if

      call complex_thin_svd(a, u, singular_values, vt, info)
      if (info /= 0) then
         call empty_ztrlan_svd(out, size(a, 1), info)
         return
      end if
      allocate(out%d(k), out%u(size(a, 1), k))
      out%d = singular_values(1:k)
      out%u = u(:, 1:k)
      out%nconv = k
      out%niter = 0
      out%nops = 0
      out%info = solver_success
   end function ztrlan_svd_complex

   pure function propack_controls(options, k, wd) result(controls)
      !! Maps PROPACK's Krylov/restart controls onto RSpectra's iterative-SVD options.
      type(svd_options), intent(in) :: options !! User controls following the upstream named-list fields.
      integer, intent(in) :: k !! Number of requested singular triplets.
      integer, intent(in) :: wd !! Minimum input dimension, which bounds the iterative subspace.
      type(svds_opts) :: controls
      integer :: ncv, restart_count

      ncv = 5 * k
      if (options%kmax > 0) ncv = options%kmax
      if (options%dim > 0) ncv = options%dim
      if (k < wd) ncv = min(wd, max(k + 1, ncv))
      if (k == wd) ncv = wd
      controls%ncv = ncv

      restart_count = 10
      if (options%maxiter > 0) restart_count = options%maxiter
      controls%maxitr = max(1, restart_count * max(1, ncv))
      controls%tol = 1.0e-12_dp
      if (options%tol >= 0.0_dp) controls%tol = options%tol
   end function propack_controls

   pure function trlan_controls(options, k, n, subspace_limit) result(controls)
      !! Maps nuTRLan's subspace/tolerance/matrix-vector budget controls onto symmetric RSpectra options.
      type(svd_options), intent(in) :: options !! User controls following the upstream named-list fields.
      integer, intent(in) :: k !! Number of requested eigenpairs or singular pairs.
      integer, intent(in) :: n !! Dimension of the symmetric operator seen by Lanczos.
      integer, intent(in), optional :: subspace_limit !! Optional upstream-style upper bound on the Krylov basis size.
      type(eigs_opts) :: controls
      integer :: limit, ncv

      limit = n
      if (present(subspace_limit)) limit = min(n, max(k + 1, subspace_limit))
      ncv = 5 * k
      if (options%kmax > 0) ncv = options%kmax
      if (k < n) ncv = min(limit, max(k + 1, ncv))
      if (k == n) ncv = n
      controls%ncv = ncv
      controls%maxitr = max(1, k * n)
      if (options%maxiter > 0) controls%maxitr = options%maxiter
      controls%tol = sqrt(epsilon(1.0_dp))
      if (options%tol >= 0.0_dp) controls%tol = options%tol
      controls%retvec = .true.
   end function trlan_controls

   pure subroutine apply_real_warm_start(controls, n, lambda_start, u_start, info)
      !! Validates the paired upstream warm start and maps one supplied vector to ARPACK's initial residual.
      type(eigs_opts), intent(inout) :: controls !! Solver controls whose optional initial vector is populated.
      integer, intent(in) :: n !! Required row dimension of supplied eigen/singular vectors.
      real(dp), intent(in), optional :: lambda_start(:) !! Prior values whose length selects the matched warm-start count.
      real(dp), intent(in), optional :: u_start(:, :) !! Previously converged vectors with `n` rows.
      integer, intent(out) :: info !! Zero when absent/valid; negative for inconsistent supplied warm-start shapes.
      real(dp) :: norm2
      integer :: j, matched

      info = solver_success
      if (.not. present(lambda_start) .or. .not. present(u_start)) return
      if (size(u_start, 1) /= n) then
         info = solver_invalid_warm_start
         return
      end if
      matched = min(size(lambda_start), size(u_start, 2))
      if (matched < 1) return
      j = matched
      norm2 = sqrt(sum(u_start(:, j) * u_start(:, j)))
      if (norm2 <= sqrt(tiny(1.0_dp))) return
      controls%initvec = u_start(:, j) / norm2
   end subroutine apply_real_warm_start

   pure subroutine validate_complex_warm_start(n, lambda_start, u_start, info)
      !! Validates optional complex nuTRLan warm-start arrays even though dense LAPACK cannot reuse them.
      integer, intent(in) :: n !! Required row dimension of the complex left-vector matrix.
      real(dp), intent(in), optional :: lambda_start(:) !! Previously converged singular values.
      complex(dp), intent(in), optional :: u_start(:, :) !! Previously converged complex left vectors.
      integer, intent(out) :: info !! Zero when absent/shape-compatible; negative for inconsistent supplied shapes.

      info = solver_success
      if (.not. present(lambda_start) .or. .not. present(u_start)) return
      if (size(u_start, 1) /= n) then
         info = solver_invalid_warm_start
         return
      end if
      if (min(size(lambda_start), size(u_start, 2)) < 0) info = solver_invalid_warm_start
   end subroutine validate_complex_warm_start

   subroutine left_normal_prod(self, x, y)
      !! Applies `A*A^T` without materializing either `A` or the normal matrix.
      class(left_normal_operator), intent(inout) :: self !! Left-normal wrapper around the base matrix-free operator.
      real(dp), intent(in) :: x(:) !! Input vector with length equal to the base operator row count.
      real(dp), intent(out) :: y(:) !! Output vector containing `A*A^T*x`.

      if (.not. associated(self%base)) error stop "svd_solvers: missing normal-operator base"
      call self%base%tprod(x, self%work)
      call self%base%prod(self%work, y)
   end subroutine left_normal_prod

   subroutine left_normal_tprod(self, x, y)
      !! Applies the transpose of the symmetric left-normal operator, identical to its forward product.
      class(left_normal_operator), intent(inout) :: self !! Left-normal wrapper around the base matrix-free operator.
      real(dp), intent(in) :: x(:) !! Input vector with length equal to the base operator row count.
      real(dp), intent(out) :: y(:) !! Output vector containing `A*A^T*x`.

      call self%prod(x, y)
   end subroutine left_normal_tprod

   subroutine materialize_operator(op, a)
      !! Materializes a generic real linear operator by applying it to every canonical basis vector.
      class(linear_operator), intent(inout) :: op !! Operator to materialize; must implement a valid forward product.
      real(dp), allocatable, intent(out) :: a(:, :) !! Allocated dense representation with shape `(nrow,ncol)`.
      real(dp), allocatable :: basis(:)
      integer :: j

      allocate(a(op%nrow, op%ncol), basis(op%ncol))
      do j = 1, op%ncol
         basis = 0.0_dp
         basis(j) = 1.0_dp
         call op%prod(basis, a(:, j))
      end do
   end subroutine materialize_operator

   pure subroutine copy_propack_result(source, target)
      !! Copies the shared SVD backend result into the package-compatible PROPACK result type.
      type(svds_result), intent(in) :: source !! RSpectra SVD result to translate.
      type(propack_svd_result), intent(out) :: target !! Package result receiving singular values, vectors, and diagnostics.

      if (allocated(source%d)) target%d = source%d
      if (allocated(source%u)) target%u = source%u
      if (allocated(source%v)) target%v = source%v
      target%nconv = source%nconv
      target%niter = source%niter
      target%nops = source%nops
      target%info = source%info
   end subroutine copy_propack_result

   pure subroutine copy_trlan_svd_result(source, target)
      !! Converts eigenpairs of `A*A^T` to nuTRLan-compatible singular values and left vectors.
      type(eigs_sym_result), intent(in) :: source !! Shared symmetric eigensolver result for the left normal operator.
      type(trlan_svd_result), intent(out) :: target !! Result receiving square-root eigenvalues and left vectors.

      if (allocated(source%values)) target%d = sqrt(max(source%values, 0.0_dp))
      if (allocated(source%vectors)) target%u = source%vectors
      target%nconv = source%nconv
      target%niter = source%niter
      target%nops = 2 * source%nops
      target%info = source%info
   end subroutine copy_trlan_svd_result

   pure subroutine copy_trlan_eigen_result(source, target)
      !! Copies a shared symmetric-eigensolver result into the nuTRLan-compatible eigen result type.
      type(eigs_sym_result), intent(in) :: source !! Shared symmetric eigensolver result.
      type(trlan_eigen_result), intent(out) :: target !! Result receiving leading eigenvalues, eigenvectors, and diagnostics.

      if (allocated(source%values)) target%d = source%values
      if (allocated(source%vectors)) target%u = source%vectors
      target%nconv = source%nconv
      target%niter = source%niter
      target%nops = source%nops
      target%info = source%info
   end subroutine copy_trlan_eigen_result

   pure subroutine empty_propack(out, m, n, info)
      !! Initializes an empty PROPACK-style result after input or backend failure.
      type(propack_svd_result), intent(out) :: out !! Result object to initialize.
      integer, intent(in) :: m !! Row dimension used for the zero-column left-vector array.
      integer, intent(in) :: n !! Column dimension used for the zero-column right-vector array.
      integer, intent(in) :: info !! Status code to store in the result.

      allocate(out%d(0), out%u(max(0, m), 0), out%v(max(0, n), 0))
      out%info = info
   end subroutine empty_propack

   pure subroutine empty_trlan_svd(out, m, info)
      !! Initializes an empty real nuTRLan-SVD result after input or backend failure.
      type(trlan_svd_result), intent(out) :: out !! Result object to initialize.
      integer, intent(in) :: m !! Row dimension used for the zero-column left-vector array.
      integer, intent(in) :: info !! Status code to store in the result.

      allocate(out%d(0), out%u(max(0, m), 0))
      out%info = info
   end subroutine empty_trlan_svd

   pure subroutine empty_ztrlan_svd(out, m, info)
      !! Initializes an empty complex nuTRLan-SVD result after input or backend failure.
      type(ztrlan_svd_result), intent(out) :: out !! Result object to initialize.
      integer, intent(in) :: m !! Row dimension used for the zero-column complex left-vector array.
      integer, intent(in) :: info !! Status code to store in the result.

      allocate(out%d(0), out%u(max(0, m), 0))
      out%info = info
   end subroutine empty_ztrlan_svd

   pure subroutine empty_trlan_eigen(out, n, info)
      !! Initializes an empty nuTRLan-eigen result after input or backend failure.
      type(trlan_eigen_result), intent(out) :: out !! Result object to initialize.
      integer, intent(in) :: n !! Matrix order used for the zero-column eigenvector array.
      integer, intent(in) :: info !! Status code to store in the result.

      allocate(out%d(0), out%u(max(0, n), 0))
      out%info = info
   end subroutine empty_trlan_eigen

end module svd_solvers
