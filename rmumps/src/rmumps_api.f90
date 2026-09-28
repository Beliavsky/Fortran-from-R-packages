! Modern Fortran translation of computational behavior from the R package rmumps.
! Upstream rmumps: Copyright 2015-2026 INRAE/INSA/CNRS, GPL-2.0-or-later.
! This translated file was substantially rewritten on 2026-09-22.
! See LICENSE, NOTICE.md, and UPSTREAM_COPYRIGHTS for license and provenance details.
module rmumps_api
   use rmumps_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: RMUMPS_PERM_AMD = 0
   integer, parameter, public :: RMUMPS_PERM_AMF = 2
   integer, parameter, public :: RMUMPS_PERM_SCOTCH = 3
   integer, parameter, public :: RMUMPS_PERM_PORD = 4
   integer, parameter, public :: RMUMPS_PERM_METIS = 5
   integer, parameter, public :: RMUMPS_PERM_QAMD = 6
   integer, parameter, public :: RMUMPS_PERM_AUTO = 7

   type, public :: rmumps_matrix
      integer :: n = 0
      integer :: sym = 0
      integer :: permutation = RMUMPS_PERM_AUTO
      logical :: symbolic_ready = .false.
      logical :: numeric_ready = .false.
      integer :: swap_count = 0
      integer :: last_info = 0
      integer :: icntl(33) = 0
      real(dp) :: cntl(5) = 0.0_dp
      integer :: keep(500) = 0
      logical :: copy = .true.
      integer :: ncore = 1
      integer, allocatable :: irn(:)
      integer, allocatable :: jcn(:)
      integer, allocatable :: piv(:)
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: dense(:, :)
      real(dp), allocatable :: lu(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: mrhs(:, :)
      integer, allocatable :: perm_in(:)
   end type rmumps_matrix

   type, public :: rmumps_infos
      integer :: info(27) = 0
      real(dp) :: rinfo(3) = 0.0_dp
      integer :: infog(34) = 0
      real(dp) :: rinfog(13) = 0.0_dp
   end type rmumps_infos

   interface rmumps_solve
      module procedure rmumps_solve_vector
      module procedure rmumps_solve_matrix
      module procedure rmumps_solve_sparse_rhs
   end interface rmumps_solve

   interface rmumps_solve_inplace
      module procedure rmumps_solve_vector_inplace
      module procedure rmumps_solve_matrix_inplace
   end interface rmumps_solve_inplace

   interface rmumps_solvet
      module procedure rmumps_solvet_vector
      module procedure rmumps_solvet_matrix
      module procedure rmumps_solvet_sparse_rhs
   end interface rmumps_solvet

   public :: dp
   public :: rmumps_init_triplet
   public :: rmumps_init_triplet_zero_based
   public :: rmumps_symbolic
   public :: rmumps_numeric
   public :: rmumps_solve
   public :: rmumps_solve_vector
   public :: rmumps_solve_matrix
   public :: rmumps_solve_sparse_rhs
   public :: rmumps_solve_inplace
   public :: rmumps_solve_vector_inplace
   public :: rmumps_solve_matrix_inplace
   public :: rmumps_solvet
   public :: rmumps_solvet_vector
   public :: rmumps_solvet_matrix
   public :: rmumps_solvet_sparse_rhs
   public :: rmumps_inverse
   public :: rmumps_inverse_transpose
   public :: rmumps_det
   public :: rmumps_log2det
   public :: rmumps_determinant
   public :: rmumps_set_values
   public :: rmumps_get_triplet
   public :: rmumps_set_permutation
   public :: rmumps_get_permutation
   public :: rmumps_get_sym
   public :: rmumps_set_copy
   public :: rmumps_get_copy
   public :: rmumps_set_ncore
   public :: rmumps_get_ncore
   public :: rmumps_set_rhs
   public :: rmumps_get_rhs
   public :: rmumps_set_mrhs
   public :: rmumps_get_mrhs
   public :: rmumps_set_perm_in
   public :: rmumps_set_icntl
   public :: rmumps_get_icntl
   public :: rmumps_set_cntl
   public :: rmumps_get_cntl
   public :: rmumps_set_keep
   public :: rmumps_get_keep
   public :: rmumps_get_infos
   public :: rmumps_dim
   public :: rmumps_nrow
   public :: rmumps_ncol

contains

   subroutine rmumps_init_triplet(a, i, j, x, n, sym, copy, info)
      type(rmumps_matrix), intent(out) :: a !! Matrix object initialized from a sparse triplet.
      integer, intent(in) :: i(:) !! One-based row indices; each entry must lie in 1:n.
      integer, intent(in) :: j(:) !! One-based column indices; each entry must lie in 1:n.
      real(dp), intent(in) :: x(:) !! Triplet values corresponding elementwise to i and j.
      integer, intent(in) :: n !! Matrix order; must be nonnegative and defines an n by n system.
      integer, intent(in), optional :: sym !! Symmetry mode: 0 general, 1 or 2 symmetric using one stored triangle; default 0.
      logical, intent(in), optional :: copy !! Compatibility copy flag stored in the object; default true.
      integer, intent(out), optional :: info !! Status: 0 success, negative for invalid dimensions, indices, or symmetry mode.
      integer :: local_info
      integer :: k
      integer :: s
      logical :: copy_values

      local_info = 0
      s = 0
      if (present(sym)) s = sym
      copy_values = .true.
      if (present(copy)) copy_values = copy

      if (n < 0) local_info = -1
      if (size(i) /= size(j) .or. size(i) /= size(x)) local_info = -2
      if (s < 0 .or. s > 2) local_info = -3

      if (local_info == 0) then
         do k = 1, size(i)
            if (i(k) < 1 .or. i(k) > n .or. j(k) < 1 .or. j(k) > n) then
               local_info = -4
               exit
            end if
         end do
      end if

      if (local_info /= 0) then
         a%last_info = local_info
         if (present(info)) info = local_info
         return
      end if

      a%n = n
      a%sym = s
      a%permutation = RMUMPS_PERM_AUTO
      a%symbolic_ready = .false.
      a%numeric_ready = .false.
      a%swap_count = 0
      a%last_info = 0
      a%copy = copy_values
      a%ncore = 1
      call initialize_controls(a)
      allocate(a%irn(size(i)), a%jcn(size(j)), a%values(size(x)))
      a%irn = i
      a%jcn = j
      a%values = x
      call assemble_dense(a)

      if (present(info)) info = 0
   end subroutine rmumps_init_triplet

   subroutine rmumps_init_triplet_zero_based(a, i0, j0, x, n, sym, copy, info)
      type(rmumps_matrix), intent(out) :: a !! Matrix object initialized from an upstream-style zero-based sparse triplet.
      integer, intent(in) :: i0(:) !! Zero-based row indices, matching Rcpp dgTMatrix slot conventions.
      integer, intent(in) :: j0(:) !! Zero-based column indices, matching Rcpp dgTMatrix slot conventions.
      real(dp), intent(in) :: x(:) !! Triplet values corresponding elementwise to i0 and j0.
      integer, intent(in) :: n !! Matrix order; must be nonnegative and defines an n by n system.
      integer, intent(in), optional :: sym !! Symmetry mode: 0 general, 1 or 2 symmetric using one stored triangle; default 0.
      logical, intent(in), optional :: copy !! Compatibility copy flag stored in the object; default true.
      integer, intent(out), optional :: info !! Status returned by one-based initialization after index conversion.
      integer, allocatable :: i(:)
      integer, allocatable :: j(:)
      integer :: local_info
      integer :: s
      logical :: copy_values

      s = 0
      if (present(sym)) s = sym
      copy_values = .true.
      if (present(copy)) copy_values = copy
      allocate(i(size(i0)), j(size(j0)))
      i = i0 + 1
      j = j0 + 1
      call rmumps_init_triplet(a, i, j, x, n, sym=s, copy=copy_values, info=local_info)
      if (present(info)) info = local_info
   end subroutine rmumps_init_triplet_zero_based

   subroutine rmumps_symbolic(a, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose symbolic-ready state is to be established.
      integer, intent(out), optional :: info !! Status: 0 when matrix storage is initialized, otherwise -1.

      if (.not. allocated(a%dense)) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      a%symbolic_ready = .true.
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_symbolic

   subroutine rmumps_numeric(a, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object to factor numerically; LU factors are cached in the object.
      integer, intent(out), optional :: info !! Status: 0 success, -10 for a singular factorization, or -1 if uninitialized.
      integer :: local_info

      if (.not. allocated(a%dense)) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if

      if (.not. a%symbolic_ready) a%symbolic_ready = .true.
      if (allocated(a%lu)) deallocate(a%lu)
      if (allocated(a%piv)) deallocate(a%piv)
      allocate(a%lu(a%n, a%n), a%piv(a%n))
      a%lu = a%dense
      call lu_factor(a%lu, a%piv, a%swap_count, local_info)
      if (local_info > 0) local_info = -10
      a%numeric_ready = local_info == 0
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_numeric

   subroutine rmumps_solve_vector(a, b, x, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object; its cached LU factorization is created or reused.
      real(dp), intent(in) :: b(:) !! Dense right-hand-side vector of length n.
      real(dp), intent(out) :: x(:) !! Solution vector of length n; undefined when info is nonzero.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      integer :: local_info

      local_info = 0
      if (size(b) /= a%n .or. size(x) /= a%n) local_info = -2
      if (local_info == 0 .and. .not. a%numeric_ready) call rmumps_numeric(a, local_info)
      if (local_info == 0) then
         x = b
         if (a%icntl(9) == 2) then
            call lu_solve_transpose_vector(a%lu, a%piv, x)
         else
            call lu_solve_vector(a%lu, a%piv, x)
         end if
         a%rhs = x
         a%icntl(20) = 0
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solve_vector

   subroutine rmumps_solve_matrix(a, b, x, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object; its cached LU factorization is created or reused.
      real(dp), intent(in) :: b(:, :) !! Dense right-hand-side matrix with n rows and zero or more columns.
      real(dp), intent(out) :: x(:, :) !! Solution matrix with the same shape as b; undefined when info is nonzero.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      integer :: local_info
      integer :: k

      local_info = 0
      if (size(b, 1) /= a%n) local_info = -2
      if (any(shape(x) /= shape(b))) local_info = -3
      if (local_info == 0 .and. .not. a%numeric_ready) call rmumps_numeric(a, local_info)
      if (local_info == 0) then
         x = b
         do k = 1, size(x, 2)
            if (a%icntl(9) == 2) then
               call lu_solve_transpose_vector(a%lu, a%piv, x(:, k))
            else
               call lu_solve_vector(a%lu, a%piv, x(:, k))
            end if
         end do
         a%mrhs = x
         a%icntl(20) = 0
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solve_matrix

   subroutine rmumps_solve_sparse_rhs(a, i, j, v, nrow, ncol, x, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object; its cached LU factorization is created or reused.
      integer, intent(in) :: i(:) !! One-based row indices of sparse right-hand-side entries.
      integer, intent(in) :: j(:) !! One-based column indices of sparse right-hand-side entries.
      real(dp), intent(in) :: v(:) !! Sparse right-hand-side values corresponding elementwise to i and j.
      integer, intent(in) :: nrow !! Sparse right-hand-side row count; normally equal to matrix order n.
      integer, intent(in) :: ncol !! Sparse right-hand-side column count; zero with nrow=0 requests the inverse.
      real(dp), intent(out) :: x(:, :) !! Dense solutions, shape n by ncol; for a 0 by 0 rhs, shape n by n.
      integer, intent(out), optional :: info !! Status: 0 success, negative sparse-input error, or -10 for a singular factorization.
      real(dp), allocatable :: b(:, :)
      integer :: local_info

      call validate_sparse_rhs(a, i, j, v, nrow, ncol, x, local_info)
      if (local_info /= 0) then
         a%last_info = local_info
         if (present(info)) info = local_info
         return
      end if

      if (nrow == 0 .and. ncol == 0) then
         call rmumps_inverse(a, x, local_info)
      else
         allocate(b(nrow, ncol))
         call assemble_sparse_rhs(i, j, v, b)
         call rmumps_solve_matrix(a, b, x, local_info)
         if (local_info == 0) a%icntl(20) = 1
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solve_sparse_rhs

   subroutine rmumps_solve_vector_inplace(a, b, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object; its cached LU factorization is created or reused.
      real(dp), intent(inout) :: b(:) !! Dense right-hand side of length n, overwritten by the solution.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      real(dp), allocatable :: x(:)
      integer :: local_info

      allocate(x(size(b)))
      call rmumps_solve_vector(a, b, x, local_info)
      if (local_info == 0) b = x
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solve_vector_inplace

   subroutine rmumps_solve_matrix_inplace(a, b, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object; its cached LU factorization is created or reused.
      real(dp), intent(inout) :: b(:, :) !! Dense right-hand-side matrix with n rows, overwritten by the solutions.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      real(dp), allocatable :: x(:, :)
      integer :: local_info

      allocate(x(size(b, 1), size(b, 2)))
      call rmumps_solve_matrix(a, b, x, local_info)
      if (local_info == 0) b = x
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solve_matrix_inplace

   subroutine rmumps_solvet_vector(a, b, x, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose transpose defines the system; cached LU factors are reused.
      real(dp), intent(in) :: b(:) !! Dense right-hand-side vector of length n.
      real(dp), intent(out) :: x(:) !! Solution of transpose(A) x = b, length n.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      integer :: local_info

      a%icntl(9) = 2
      call rmumps_solve_vector(a, b, x, local_info)
      a%icntl(9) = 1
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solvet_vector

   subroutine rmumps_solvet_matrix(a, b, x, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose transpose defines the system; cached LU factors are reused.
      real(dp), intent(in) :: b(:, :) !! Dense right-hand-side matrix with n rows and zero or more columns.
      real(dp), intent(out) :: x(:, :) !! Solutions of transpose(A) X = B with the same shape as b.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      integer :: local_info

      a%icntl(9) = 2
      call rmumps_solve_matrix(a, b, x, local_info)
      a%icntl(9) = 1
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solvet_matrix

   subroutine rmumps_solvet_sparse_rhs(a, i, j, v, nrow, ncol, x, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose transpose defines the system; cached LU factors are reused.
      integer, intent(in) :: i(:) !! One-based row indices of sparse right-hand-side entries.
      integer, intent(in) :: j(:) !! One-based column indices of sparse right-hand-side entries.
      real(dp), intent(in) :: v(:) !! Sparse right-hand-side values corresponding elementwise to i and j.
      integer, intent(in) :: nrow !! Sparse right-hand-side row count; normally equal to matrix order n.
      integer, intent(in) :: ncol !! Sparse right-hand-side column count; zero with nrow=0 requests inverse transpose.
      real(dp), intent(out) :: x(:, :) !! Dense transpose-system solutions; for a 0 by 0 rhs, shape n by n.
      integer, intent(out), optional :: info !! Status: 0 success, negative sparse-input error, or -10 for a singular factorization.
      real(dp), allocatable :: b(:, :)
      integer :: local_info

      call validate_sparse_rhs(a, i, j, v, nrow, ncol, x, local_info)
      if (local_info /= 0) then
         a%last_info = local_info
         if (present(info)) info = local_info
         return
      end if

      if (nrow == 0 .and. ncol == 0) then
         call rmumps_inverse_transpose(a, x, local_info)
      else
         allocate(b(nrow, ncol))
         call assemble_sparse_rhs(i, j, v, b)
         call rmumps_solvet_matrix(a, b, x, local_info)
         if (local_info == 0) a%icntl(20) = 1
      end if
      a%icntl(9) = 1
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_solvet_sparse_rhs

   subroutine rmumps_inverse(a, ainv, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object to invert through its cached LU factorization.
      real(dp), intent(out) :: ainv(:, :) !! Dense n by n inverse matrix; undefined when info is nonzero.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      integer :: j
      integer :: local_info

      local_info = 0
      if (size(ainv, 1) /= a%n .or. size(ainv, 2) /= a%n) local_info = -2
      if (local_info == 0 .and. .not. a%numeric_ready) call rmumps_numeric(a, local_info)
      if (local_info == 0) then
         ainv = 0.0_dp
         do j = 1, a%n
            ainv(j, j) = 1.0_dp
            call lu_solve_vector(a%lu, a%piv, ainv(:, j))
         end do
         a%mrhs = ainv
         a%icntl(20) = 1
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_inverse

   subroutine rmumps_inverse_transpose(a, ainvt, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose inverse transpose is requested using cached LU factors.
      real(dp), intent(out) :: ainvt(:, :) !! Dense n by n inverse of transpose(A); undefined when info is nonzero.
      integer, intent(out), optional :: info !! Status: 0 success, negative shape error, or -10 for a singular factorization.
      integer :: j
      integer :: local_info

      local_info = 0
      if (size(ainvt, 1) /= a%n .or. size(ainvt, 2) /= a%n) local_info = -2
      if (local_info == 0 .and. .not. a%numeric_ready) call rmumps_numeric(a, local_info)
      if (local_info == 0) then
         ainvt = 0.0_dp
         do j = 1, a%n
            ainvt(j, j) = 1.0_dp
            call lu_solve_transpose_vector(a%lu, a%piv, ainvt(:, j))
         end do
         a%mrhs = ainvt
         a%icntl(20) = 1
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_inverse_transpose

   function rmumps_det(a, info) result(d)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose determinant is requested; factorization is cached.
      integer, intent(out), optional :: info !! Status: 0 success, -10 if singular, or negative initialization error.
      real(dp) :: d
      integer :: k
      integer :: local_info

      local_info = 0
      if (.not. a%numeric_ready) call rmumps_numeric(a, local_info)
      if (local_info /= 0) then
         d = 0.0_dp
      else
         d = 1.0_dp
         do k = 1, a%n
            d = d * a%lu(k, k)
         end do
         if (mod(a%swap_count, 2) /= 0) d = -d
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end function rmumps_det

   function rmumps_log2det(a, sign_det, info) result(logabsdet)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose base-2 log absolute determinant is requested.
      real(dp), intent(out) :: sign_det !! Determinant sign: -1, 0, or +1; zero for a singular matrix.
      integer, intent(out), optional :: info !! Status: 0 success, -10 if singular, or negative initialization error.
      real(dp) :: logabsdet
      integer :: local_info

      local_info = 0
      if (.not. a%numeric_ready) call rmumps_numeric(a, local_info)
      if (local_info /= 0) then
         sign_det = 0.0_dp
         logabsdet = -huge(1.0_dp)
      else
         call factor_log2det(a%lu, a%swap_count, sign_det, logabsdet)
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end function rmumps_log2det

   subroutine rmumps_determinant(a, modulus, sign_det, logarithm, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose R-style determinant decomposition is requested.
      real(dp), intent(out) :: modulus !! Absolute determinant, or natural log absolute determinant when logarithm is true.
      real(dp), intent(out) :: sign_det !! Determinant sign: -1, 0, or +1; zero for a singular matrix.
      logical, intent(in), optional :: logarithm !! If true, return log(abs(det(A))); default true as in R determinant().
      integer, intent(out), optional :: info !! Status: 0 success, -10 if singular, or negative initialization error.
      logical :: use_log
      real(dp) :: log2abs
      integer :: local_info

      use_log = .true.
      if (present(logarithm)) use_log = logarithm
      log2abs = rmumps_log2det(a, sign_det, local_info)
      if (local_info /= 0) then
         if (use_log) then
            modulus = -huge(1.0_dp)
         else
            modulus = 0.0_dp
         end if
      else if (use_log) then
         modulus = log2abs * log(2.0_dp)
      else
         modulus = exp(log2abs * log(2.0_dp))
      end if
      a%last_info = local_info
      if (present(info)) info = local_info
   end subroutine rmumps_determinant

   subroutine rmumps_set_values(a, x, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose existing triplet values are replaced in pattern order.
      real(dp), intent(in) :: x(:) !! Replacement values; length must equal the stored triplet nonzero count.
      integer, intent(out), optional :: info !! Status: 0 success or -1 for a size mismatch or uninitialized pattern.

      if (.not. allocated(a%values) .or. size(x) /= size(a%values)) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      a%values = x
      call assemble_dense(a)
      call invalidate_factorization(a, invalidate_symbolic=.false.)
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_values

   subroutine rmumps_get_triplet(a, i, j, x)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose stored sparse triplet is copied to outputs.
      integer, allocatable, intent(out) :: i(:) !! Allocated one-based row-index vector copied from the stored triplet.
      integer, allocatable, intent(out) :: j(:) !! Allocated one-based column-index vector copied from the stored triplet.
      real(dp), allocatable, intent(out) :: x(:) !! Allocated value vector copied from the stored triplet.

      if (.not. allocated(a%values)) then
         allocate(i(0), j(0), x(0))
         return
      end if
      allocate(i(size(a%irn)), j(size(a%jcn)), x(size(a%values)))
      i = a%irn
      j = a%jcn
      x = a%values
   end subroutine rmumps_get_triplet

   subroutine rmumps_set_permutation(a, permutation, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose requested ordering code is stored.
      integer, intent(in) :: permutation !! Ordering code: 0, 2, 3, 4, 5, 6, or 7 matching rmumps constants.
      integer, intent(out), optional :: info !! Status: 0 success or -1 for a code rejected by upstream rmumps rules.

      if (permutation < RMUMPS_PERM_AMD .or. permutation > RMUMPS_PERM_AUTO .or. permutation == 1) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      if (permutation /= a%permutation) call invalidate_factorization(a, invalidate_symbolic=.true.)
      a%permutation = permutation
      a%icntl(7) = permutation
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_permutation

   pure integer function rmumps_get_permutation(a) result(permutation)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose stored ordering code is queried.

      permutation = a%permutation
   end function rmumps_get_permutation

   pure integer function rmumps_get_sym(a) result(sym)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose symmetry mode is queried: 0 general, 1 or 2 symmetric.

      sym = a%sym
   end function rmumps_get_sym

   pure subroutine rmumps_set_copy(a, copy)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose R-wrapper copy flag is stored for compatibility.
      logical, intent(in) :: copy !! Requested copy behavior; Fortran allocatable assignments still have value semantics.

      a%copy = copy
   end subroutine rmumps_set_copy

   pure logical function rmumps_get_copy(a) result(copy)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose stored R-wrapper copy flag is queried.

      copy = a%copy
   end function rmumps_get_copy

   pure subroutine rmumps_set_ncore(a, ncore, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose requested core count is stored for compatibility.
      integer, intent(in) :: ncore !! Requested positive core count; the current serial dense backend does not use it.
      integer, intent(out), optional :: info !! Status: 0 success or -1 when ncore is less than one.

      if (ncore < 1) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      a%ncore = ncore
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_ncore

   pure integer function rmumps_get_ncore(a) result(ncore)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose stored requested core count is queried.

      ncore = a%ncore
   end function rmumps_get_ncore

   pure subroutine rmumps_set_rhs(a, b, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose single dense right-hand-side state is replaced.
      real(dp), intent(in) :: b(:) !! Dense right-hand-side vector; its length must equal the matrix order.
      integer, intent(out), optional :: info !! Status: 0 success or -1 when the vector length differs from the matrix order.

      if (size(b) /= a%n) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      a%rhs = b
      a%icntl(20) = 0
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_rhs

   pure subroutine rmumps_get_rhs(a, b)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose stored single dense right-hand side is copied out.
      real(dp), allocatable, intent(out) :: b(:) !! Allocated stored vector, or a zero-length vector when none is stored.

      if (allocated(a%rhs)) then
         allocate(b(size(a%rhs)))
         b = a%rhs
      else
         allocate(b(0))
      end if
   end subroutine rmumps_get_rhs

   pure subroutine rmumps_set_mrhs(a, b, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose multiple dense right-hand-side state is replaced.
      real(dp), intent(in) :: b(:, :) !! Dense right-hand-side matrix; its row count must equal the matrix order.
      integer, intent(out), optional :: info !! Status: 0 success or -1 when the row count differs from the matrix order.

      if (size(b, 1) /= a%n) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      a%mrhs = b
      a%icntl(20) = 0
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_mrhs

   pure subroutine rmumps_get_mrhs(a, b)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose stored multiple dense right-hand sides are copied out.
      real(dp), allocatable, intent(out) :: b(:, :) !! Allocated stored matrix, or a zero-by-zero matrix when none is stored.

      if (allocated(a%mrhs)) then
         allocate(b(size(a%mrhs, 1), size(a%mrhs, 2)))
         b = a%mrhs
      else
         allocate(b(0, 0))
      end if
   end subroutine rmumps_get_mrhs

   pure subroutine rmumps_set_perm_in(a, permutation, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object receiving a user-supplied symmetric permutation request.
      integer, intent(in) :: permutation(:) !! One-based permutation of 1:n; stored but not applied by the dense fallback.
      integer, intent(out), optional :: info !! Status: 0 success, -1 size mismatch, or -2 invalid/duplicate permutation entry.
      logical, allocatable :: seen(:)
      integer :: k

      if (size(permutation) /= a%n) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      allocate(seen(a%n))
      seen = .false.
      do k = 1, a%n
         if (permutation(k) < 1 .or. permutation(k) > a%n) then
            a%last_info = -2
            if (present(info)) info = -2
            return
         end if
         if (seen(permutation(k))) then
            a%last_info = -2
            if (present(info)) info = -2
            return
         end if
         seen(permutation(k)) = .true.
      end do

      a%perm_in = permutation
      a%icntl(7) = 1
      a%icntl(28) = 1
      a%permutation = 1
      call invalidate_factorization(a, invalidate_symbolic=.true.)
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_perm_in

   subroutine rmumps_set_icntl(a, values, indices, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose MUMPS-style integer control state is updated.
      integer, intent(in) :: values(:) !! Integer control values corresponding elementwise to indices.
      integer, intent(in) :: indices(:) !! One-based ICNTL positions; values outside 1:33 are ignored as upstream does.
      integer, intent(out), optional :: info !! Status: 0 success or -1 when values and indices have different lengths.
      integer :: k

      if (size(values) /= size(indices)) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      do k = 1, size(indices)
         if (indices(k) < 1 .or. indices(k) > 33) cycle
         a%icntl(indices(k)) = values(k)
         if (indices(k) == 7) a%permutation = values(k)
      end do
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_icntl

   pure function rmumps_get_icntl(a) result(values)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose 33-entry MUMPS-style integer control state is queried.
      integer :: values(33)

      values = a%icntl
   end function rmumps_get_icntl

   subroutine rmumps_set_cntl(a, values, indices, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose MUMPS-style real control state is updated.
      real(dp), intent(in) :: values(:) !! Real control values corresponding elementwise to indices.
      integer, intent(in) :: indices(:) !! One-based CNTL positions; values outside 1:5 are ignored as upstream does.
      integer, intent(out), optional :: info !! Status: 0 success or -1 when values and indices have different lengths.
      integer :: k

      if (size(values) /= size(indices)) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      do k = 1, size(indices)
         if (indices(k) < 1 .or. indices(k) > 5) cycle
         a%cntl(indices(k)) = values(k)
      end do
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_cntl

   pure function rmumps_get_cntl(a) result(values)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose five-entry MUMPS-style real control state is queried.
      real(dp) :: values(5)

      values = a%cntl
   end function rmumps_get_cntl

   subroutine rmumps_set_keep(a, values, indices, info)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose MUMPS-style KEEP state is updated.
      integer, intent(in) :: values(:) !! Integer KEEP values corresponding elementwise to indices.
      integer, intent(in) :: indices(:) !! One-based KEEP positions; values outside 1:500 are ignored as upstream does.
      integer, intent(out), optional :: info !! Status: 0 success or -1 when values and indices have different lengths.
      integer :: k

      if (size(values) /= size(indices)) then
         a%last_info = -1
         if (present(info)) info = -1
         return
      end if
      do k = 1, size(indices)
         if (indices(k) < 1 .or. indices(k) > 500) cycle
         a%keep(indices(k)) = values(k)
      end do
      a%last_info = 0
      if (present(info)) info = 0
   end subroutine rmumps_set_keep

   pure function rmumps_get_keep(a) result(values)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose 500-entry MUMPS-style KEEP state is queried.
      integer :: values(500)

      values = a%keep
   end function rmumps_get_keep

   subroutine rmumps_get_infos(a, infos)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose available backend status is mapped into MUMPS-shaped arrays.
      type(rmumps_infos), intent(out) :: infos !! MUMPS-shaped status arrays; unsupported entries are returned as zero.
      real(dp) :: log2abs
      real(dp) :: sign_det
      real(dp) :: exponent_real
      integer :: exponent_value

      infos%info = 0
      infos%rinfo = 0.0_dp
      infos%infog = 0
      infos%rinfog = 0.0_dp
      infos%info(1) = a%last_info
      infos%infog(1) = a%last_info

      if (.not. a%numeric_ready) return
      call factor_log2det(a%lu, a%swap_count, sign_det, log2abs)
      exponent_real = floor(log2abs)
      exponent_value = int(exponent_real)
      infos%infog(34) = exponent_value
      infos%rinfog(12) = sign_det * exp((log2abs - exponent_real) * log(2.0_dp))
   end subroutine rmumps_get_infos

   pure function rmumps_dim(a) result(dims)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose two-dimensional square shape is queried.
      integer :: dims(2)

      dims = [a%n, a%n]
   end function rmumps_dim

   pure integer function rmumps_nrow(a) result(n)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose row count is queried.

      n = a%n
   end function rmumps_nrow

   pure integer function rmumps_ncol(a) result(n)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose column count is queried.

      n = a%n
   end function rmumps_ncol

   subroutine initialize_controls(a)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose control arrays are reset to upstream rmumps defaults.
      integer :: integer_bytes

      a%icntl = 0
      a%cntl = 0.0_dp
      a%keep = 0
      a%icntl(1) = -1
      a%icntl(2) = -1
      a%icntl(3) = -1
      a%icntl(4) = 0
      a%icntl(5) = 0
      a%icntl(6) = 0
      a%icntl(7) = RMUMPS_PERM_AUTO
      a%icntl(8) = 77
      a%icntl(9) = 1
      a%icntl(10) = 0
      a%icntl(11) = 0
      a%icntl(12) = 0
      a%icntl(13) = 0
      a%icntl(14) = 50
      a%icntl(18) = 0
      a%icntl(19) = 0
      a%icntl(20) = 0
      a%icntl(21) = 0
      a%icntl(22) = 0
      a%icntl(23) = 0
      a%icntl(24) = 0
      a%icntl(25) = 0
      a%icntl(26) = 0
      a%icntl(27) = -8
      a%icntl(28) = 0
      a%icntl(29) = 0
      a%icntl(30) = 0
      a%icntl(31) = 0
      a%icntl(32) = 0
      a%icntl(33) = 1
      integer_bytes = storage_size(0) / 8
      a%keep(10) = 8 / integer_bytes
      a%keep(34) = integer_bytes
      a%keep(40) = 0
   end subroutine initialize_controls

   pure subroutine invalidate_factorization(a, invalidate_symbolic)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose cached analysis/factorization state is invalidated.
      logical, intent(in) :: invalidate_symbolic !! If true, symbolic-ready state is also cleared in addition to numeric state.

      if (invalidate_symbolic) a%symbolic_ready = .false.
      a%numeric_ready = .false.
      a%swap_count = 0
      if (allocated(a%lu)) deallocate(a%lu)
      if (allocated(a%piv)) deallocate(a%piv)
   end subroutine invalidate_factorization

   pure subroutine factor_log2det(lu, swap_count, sign_det, log2abs)
      real(dp), intent(in) :: lu(:, :) !! Compact LU factors whose diagonal determines determinant magnitude and sign.
      integer, intent(in) :: swap_count !! Number of row swaps used by the LU factorization.
      real(dp), intent(out) :: sign_det !! Determinant sign, either -1 or +1 for a nonsingular factorization.
      real(dp), intent(out) :: log2abs !! Base-2 logarithm of the absolute determinant.
      integer :: k

      sign_det = 1.0_dp
      if (mod(swap_count, 2) /= 0) sign_det = -1.0_dp
      log2abs = 0.0_dp
      do k = 1, size(lu, 1)
         if (lu(k, k) < 0.0_dp) sign_det = -sign_det
         log2abs = log2abs + log(abs(lu(k, k))) / log(2.0_dp)
      end do
   end subroutine factor_log2det

   subroutine assemble_dense(a)
      type(rmumps_matrix), intent(inout) :: a !! Matrix object whose triplet data are assembled into dense storage.
      integer :: k
      integer :: row
      integer :: col

      if (allocated(a%dense)) deallocate(a%dense)
      allocate(a%dense(a%n, a%n))
      a%dense = 0.0_dp
      do k = 1, size(a%values)
         row = a%irn(k)
         col = a%jcn(k)
         a%dense(row, col) = a%dense(row, col) + a%values(k)
         if (a%sym /= 0 .and. row /= col) then
            a%dense(col, row) = a%dense(col, row) + a%values(k)
         end if
      end do
   end subroutine assemble_dense

   pure subroutine assemble_sparse_rhs(i, j, v, b)
      integer, intent(in) :: i(:) !! One-based sparse right-hand-side row indices already validated against b.
      integer, intent(in) :: j(:) !! One-based sparse right-hand-side column indices already validated against b.
      real(dp), intent(in) :: v(:) !! Sparse right-hand-side values corresponding elementwise to i and j.
      real(dp), intent(out) :: b(:, :) !! Dense right-hand-side matrix assembled by summing duplicate triplet entries.
      integer :: k

      b = 0.0_dp
      do k = 1, size(v)
         b(i(k), j(k)) = b(i(k), j(k)) + v(k)
      end do
   end subroutine assemble_sparse_rhs

   pure subroutine validate_sparse_rhs(a, i, j, v, nrow, ncol, x, info)
      type(rmumps_matrix), intent(in) :: a !! Matrix object whose order constrains the sparse right-hand side.
      integer, intent(in) :: i(:) !! One-based row indices to validate.
      integer, intent(in) :: j(:) !! One-based column indices to validate.
      real(dp), intent(in) :: v(:) !! Sparse values whose length must match i and j.
      integer, intent(in) :: nrow !! Declared sparse right-hand-side row count.
      integer, intent(in) :: ncol !! Declared sparse right-hand-side column count.
      real(dp), intent(out) :: x(:, :) !! Output array whose required shape is validated but whose values are not set here.
      integer, intent(out) :: info !! Status: 0 valid, -2 dimensions, -3 output shape, -4 lengths, or -5 indices.
      integer :: k

      info = 0
      if (nrow < 0 .or. ncol < 0) then
         info = -2
         return
      end if
      if (size(i) /= size(j) .or. size(i) /= size(v)) then
         info = -4
         return
      end if

      if (nrow == 0 .and. ncol == 0) then
         if (size(i) /= 0) then
            info = -4
            return
         end if
         if (size(x, 1) /= a%n .or. size(x, 2) /= a%n) info = -3
         return
      end if

      if (nrow /= a%n) then
         info = -2
         return
      end if
      if (size(x, 1) /= a%n .or. size(x, 2) /= ncol) then
         info = -3
         return
      end if
      do k = 1, size(i)
         if (i(k) < 1 .or. i(k) > nrow .or. j(k) < 1 .or. j(k) > ncol) then
            info = -5
            return
         end if
      end do
   end subroutine validate_sparse_rhs

   pure subroutine lu_factor(a, piv, swap_count, info)
      real(dp), intent(inout) :: a(:, :) !! Square matrix overwritten by compact LU factors with unit lower diagonal implied.
      integer, intent(out) :: piv(:) !! Pivot row chosen at each elimination step; length must equal matrix order.
      integer, intent(out) :: swap_count !! Number of row interchanges, used to recover determinant sign.
      integer, intent(out) :: info !! Status: 0 success, positive first singular pivot index, or negative shape mismatch.
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: p
      real(dp) :: pivot_abs
      real(dp) :: temp
      real(dp) :: tol

      n = size(a, 1)
      info = 0
      swap_count = 0
      if (size(a, 2) /= n .or. size(piv) /= n) then
         info = -1
         return
      end if
      if (n == 0) return
      tol = epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a))) * real(max(1, n), dp)

      do k = 1, n
         p = k
         pivot_abs = abs(a(k, k))
         do i = k + 1, n
            if (abs(a(i, k)) > pivot_abs) then
               p = i
               pivot_abs = abs(a(i, k))
            end if
         end do
         piv(k) = p
         if (pivot_abs <= tol) then
            info = k
            return
         end if
         if (p /= k) then
            do j = 1, n
               temp = a(k, j)
               a(k, j) = a(p, j)
               a(p, j) = temp
            end do
            swap_count = swap_count + 1
         end if
         do i = k + 1, n
            a(i, k) = a(i, k) / a(k, k)
            do j = k + 1, n
               a(i, j) = a(i, j) - a(i, k) * a(k, j)
            end do
         end do
      end do
   end subroutine lu_factor

   pure subroutine lu_solve_vector(lu, piv, b)
      real(dp), intent(in) :: lu(:, :) !! Compact LU factors produced by lu_factor.
      integer, intent(in) :: piv(:) !! Pivot rows produced by lu_factor.
      real(dp), intent(inout) :: b(:) !! Right-hand side overwritten by the solution of A x = b.
      integer :: i
      integer :: k
      integer :: n
      integer :: p
      real(dp) :: temp

      n = size(b)
      do k = 1, n
         p = piv(k)
         if (p /= k) then
            temp = b(k)
            b(k) = b(p)
            b(p) = temp
         end if
      end do
      do i = 2, n
         do k = 1, i - 1
            b(i) = b(i) - lu(i, k) * b(k)
         end do
      end do
      do i = n, 1, -1
         do k = i + 1, n
            b(i) = b(i) - lu(i, k) * b(k)
         end do
         b(i) = b(i) / lu(i, i)
      end do
   end subroutine lu_solve_vector

   pure subroutine lu_solve_transpose_vector(lu, piv, b)
      real(dp), intent(in) :: lu(:, :) !! Compact LU factors satisfying P A = L U.
      integer, intent(in) :: piv(:) !! Pivot rows produced by lu_factor and encoding the row permutation P.
      real(dp), intent(inout) :: b(:) !! Right-hand side overwritten by the solution of transpose(A) x = b.
      integer :: i
      integer :: k
      integer :: n
      integer :: p
      real(dp) :: temp

      n = size(b)

      ! U^T y = b.
      do i = 1, n
         do k = 1, i - 1
            b(i) = b(i) - lu(k, i) * b(k)
         end do
         b(i) = b(i) / lu(i, i)
      end do

      ! L^T z = y, where L has an implied unit diagonal.
      do i = n - 1, 1, -1
         do k = i + 1, n
            b(i) = b(i) - lu(k, i) * b(k)
         end do
      end do

      ! P x = z, so x = P^T z. Undo the pivot swaps in reverse order.
      do k = n, 1, -1
         p = piv(k)
         if (p /= k) then
            temp = b(k)
            b(k) = b(p)
            b(p) = temp
         end if
      end do
   end subroutine lu_solve_transpose_vector

end module rmumps_api
