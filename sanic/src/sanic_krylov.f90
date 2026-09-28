! SPDX-License-Identifier: GPL-3.0-or-later
! Arnoldi and Lanczos iterations translated from sanic.
module sanic_krylov
   use sanic_kinds, only: dp
   use sanic_sparse, only: sparse_matrix, sparse_matvec, sparse_is_symmetric
   use sanic_spectral, only: symmetric_eigen_real, general_eigen_complex
   implicit none
   private

   type, public :: krylov_result
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: diagonal(:)
      real(dp), allocatable :: subdiagonal(:)
      real(dp), allocatable :: q(:, :)
      complex(dp), allocatable :: values(:)
      integer :: iterations = 0
      integer :: info = 0
      logical :: symmetric = .false.
   end type krylov_result

   public :: arnoldi_sparse
   public :: lanczos_sparse

contains

   pure subroutine arnoldi_sparse(a, out, b, symmetric, iter, tol, eigen, orthogonalise)
      !! Compute an Arnoldi basis, dispatching to Lanczos when symmetric is true or detected.
      type(sparse_matrix), intent(in) :: a !! Sparse square matrix whose Krylov projection is requested.
      type(krylov_result), intent(out) :: out !! Hessenberg/tridiagonal projection, basis, eigenvalues, and status.
      real(dp), intent(in), optional :: b(:) !! Nonzero starting vector; absent uses deterministic entries 1,2,...,n.
      logical, intent(in), optional :: symmetric !! If absent, sparse symmetry is checked at machine precision.
      integer, intent(in), optional :: iter !! Requested Krylov dimension in 1..a%nrow; default a%nrow.
      real(dp), intent(in), optional :: tol !! Breakdown tolerance; zero/absent uses machine epsilon.
      logical, intent(in), optional :: eigen !! Whether projected eigenvalues are computed; default true.
      logical, intent(in), optional :: orthogonalise !! Lanczos full-reorthogonalization flag; ignored for nonsymmetric Arnoldi.
      logical :: use_symmetric

      use_symmetric = sparse_is_symmetric(a)
      if (present(symmetric)) use_symmetric = symmetric
      if (use_symmetric) then
         call lanczos_sparse(a, out, b, iter, tol, eigen, orthogonalise)
      else
         call arnoldi_nonsymmetric(a, out, b, iter, tol, eigen)
      end if
   end subroutine arnoldi_sparse

   pure subroutine lanczos_sparse(a, out, b, iter, tol, eigen, orthogonalise)
      !! Compute the Lanczos tridiagonal projection of a real symmetric sparse matrix.
      type(sparse_matrix), intent(in) :: a !! Sparse square symmetric matrix.
      type(krylov_result), intent(out) :: out !! Tridiagonal coefficients, basis, eigenvalues, and status.
      real(dp), intent(in), optional :: b(:) !! Nonzero starting vector; absent uses deterministic entries 1,2,...,n.
      integer, intent(in), optional :: iter !! Requested Krylov dimension in 1..a%nrow; default a%nrow.
      real(dp), intent(in), optional :: tol !! Breakdown tolerance; zero/absent uses machine epsilon.
      logical, intent(in), optional :: eigen !! Whether projected eigenvalues are computed; default true.
      logical, intent(in), optional :: orthogonalise !! If true, apply full modified Gram-Schmidt reorthogonalization.
      real(dp), allocatable :: d(:), q(:, :), s(:), start(:), v(:), tri(:, :), eigval(:), eigvec(:, :)
      real(dp) :: norm_start, norm_v, tolerance
      integer :: actual, i, j, requested, status
      logical :: compute_eigen, reorthogonalize

      out%symmetric = .true.
      if (a%nrow /= a%ncol .or. a%nrow < 1) then
         out%info = 1
         call allocate_empty_krylov(out)
         return
      end if
      requested = a%nrow
      if (present(iter)) requested = iter
      if (requested < 1 .or. requested > a%nrow) then
         out%info = 1
         call allocate_empty_krylov(out)
         return
      end if
      if (present(tol)) then
         if (tol < 0.0_dp) then
            out%info = 1
            call allocate_empty_krylov(out)
            return
         end if
      end if
      tolerance = epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) tolerance = tol
      end if
      compute_eigen = .true.
      if (present(eigen)) compute_eigen = eigen
      reorthogonalize = .true.
      if (present(orthogonalise)) reorthogonalize = orthogonalise

      call starting_vector(a%nrow, b, start, status)
      if (status /= 0) then
         out%info = 1
         call allocate_empty_krylov(out)
         return
      end if
      norm_start = sqrt(dot_product(start, start))
      allocate(q(a%nrow, requested), d(requested), s(max(0, requested - 1)))
      q = 0.0_dp
      d = 0.0_dp
      s = 0.0_dp
      q(:, 1) = start/norm_start

      call sparse_matvec(a, q(:, 1), v, status)
      if (status /= 0) then
         out%info = 1
         call allocate_empty_krylov(out)
         return
      end if
      d(1) = dot_product(q(:, 1), v)
      v = v - d(1)*q(:, 1)
      actual = 1

      do i = 2, requested
         if (reorthogonalize) then
            do j = 1, i - 1
               v = v - dot_product(q(:, j), v)*q(:, j)
            end do
         end if
         norm_v = sqrt(max(0.0_dp, dot_product(v, v)))
         if (norm_v < tolerance) exit
         s(i - 1) = norm_v
         q(:, i) = v/norm_v
         actual = i
         call sparse_matvec(a, q(:, i), v, status)
         if (status /= 0) then
            out%info = 1
            call allocate_empty_krylov(out)
            return
         end if
         d(i) = dot_product(q(:, i), v)
         v = v - d(i)*q(:, i) - s(i - 1)*q(:, i - 1)
      end do

      out%iterations = actual
      out%q = q
      out%diagonal = d
      out%subdiagonal = s
      allocate(out%h(0, 0))

      if (compute_eigen) then
         allocate(tri(requested, requested))
         tri = 0.0_dp
         do i = 1, requested
            tri(i, i) = out%diagonal(i)
         end do
         do i = 1, requested - 1
            tri(i, i + 1) = out%subdiagonal(i)
            tri(i + 1, i) = out%subdiagonal(i)
         end do
         call symmetric_eigen_real(tri, eigval, eigvec, status)
         out%values = cmplx(eigval, 0.0_dp, kind=dp)
         if (status /= 0) then
            out%info = 2
         else
            out%info = 0
         end if
      else
         allocate(out%values(0))
         out%info = 0
      end if
   end subroutine lanczos_sparse

   pure subroutine arnoldi_nonsymmetric(a, out, b, iter, tol, eigen)
      !! Compute a modified-Gram-Schmidt Arnoldi Hessenberg projection.
      type(sparse_matrix), intent(in) :: a !! Sparse square nonsymmetric matrix.
      type(krylov_result), intent(out) :: out !! Hessenberg projection, basis, projected eigenvalues, and status.
      real(dp), intent(in), optional :: b(:) !! Nonzero starting vector; absent uses deterministic entries 1,2,...,n.
      integer, intent(in), optional :: iter !! Requested Krylov dimension in 1..a%nrow; default a%nrow.
      real(dp), intent(in), optional :: tol !! Breakdown tolerance; zero/absent uses machine epsilon.
      logical, intent(in), optional :: eigen !! Whether projected eigenvalues are computed; default true.
      real(dp), allocatable :: h(:, :), q(:, :), start(:), v(:)
      complex(dp), allocatable :: eigval(:), eigvec(:, :)
      real(dp) :: norm_start, norm_v, tolerance
      integer :: actual, i, j, requested, status
      logical :: compute_eigen

      out%symmetric = .false.
      if (a%nrow /= a%ncol .or. a%nrow < 1) then
         out%info = 1
         call allocate_empty_krylov(out)
         return
      end if
      requested = a%nrow
      if (present(iter)) requested = iter
      if (requested < 1 .or. requested > a%nrow) then
         out%info = 1
         call allocate_empty_krylov(out)
         return
      end if
      if (present(tol)) then
         if (tol < 0.0_dp) then
            out%info = 1
            call allocate_empty_krylov(out)
            return
         end if
      end if
      tolerance = epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) tolerance = tol
      end if
      compute_eigen = .true.
      if (present(eigen)) compute_eigen = eigen

      call starting_vector(a%nrow, b, start, status)
      if (status /= 0) then
         out%info = 1
         call allocate_empty_krylov(out)
         return
      end if
      norm_start = sqrt(dot_product(start, start))
      allocate(q(a%nrow, requested), h(requested, requested))
      q = 0.0_dp
      h = 0.0_dp
      q(:, 1) = start/norm_start
      actual = 0

      do i = 1, requested
         call sparse_matvec(a, q(:, i), v, status)
         if (status /= 0) then
            out%info = 1
            call allocate_empty_krylov(out)
            return
         end if
         do j = 1, i
            h(j, i) = dot_product(q(:, j), v)
            v = v - h(j, i)*q(:, j)
         end do
         actual = i
         norm_v = sqrt(max(0.0_dp, dot_product(v, v)))
         if (norm_v < tolerance) exit
         if (i < requested) then
            h(i + 1, i) = norm_v
            q(:, i + 1) = v/norm_v
         end if
      end do

      out%iterations = actual
      out%h = h
      out%q = q
      allocate(out%diagonal(0), out%subdiagonal(0))
      if (compute_eigen) then
         call general_eigen_complex(out%h, eigval, eigvec, status, compute_vectors=.false.)
         out%values = eigval
         if (status /= 0) then
            out%info = 2
         else
            out%info = 0
         end if
      else
         allocate(out%values(0))
         out%info = 0
      end if
   end subroutine arnoldi_nonsymmetric

   pure subroutine starting_vector(n, b, start, info)
      !! Validate a supplied start vector or create a deterministic nonzero default.
      integer, intent(in) :: n !! Required starting-vector length.
      real(dp), intent(in), optional :: b(:) !! Optional user starting vector of length n.
      real(dp), allocatable, intent(out) :: start(:) !! Validated or deterministic starting vector.
      integer, intent(out) :: info !! Status: zero success, one for wrong length or zero norm.
      integer :: i

      info = 0
      allocate(start(n))
      if (present(b)) then
         if (size(b) /= n) then
            info = 1
            start = 0.0_dp
            return
         end if
         start = b
      else
         do i = 1, n
            start(i) = real(i, dp)
         end do
      end if
      if (dot_product(start, start) <= tiny(1.0_dp)) info = 1
   end subroutine starting_vector

   pure subroutine allocate_empty_krylov(out)
      !! Allocate zero-size arrays for an invalid Krylov request.
      type(krylov_result), intent(inout) :: out !! Result object receiving zero-size array components.

      if (.not. allocated(out%h)) allocate(out%h(0, 0))
      if (.not. allocated(out%diagonal)) allocate(out%diagonal(0))
      if (.not. allocated(out%subdiagonal)) allocate(out%subdiagonal(0))
      if (.not. allocated(out%q)) allocate(out%q(0, 0))
      if (.not. allocated(out%values)) allocate(out%values(0))
      out%iterations = 0
   end subroutine allocate_empty_krylov

end module sanic_krylov
