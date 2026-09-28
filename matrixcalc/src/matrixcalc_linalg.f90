module matrixcalc_linalg
   use r_kinds, only : dp
   use r_linalg, only : solve_system, symmetric_eigen
   implicit none
   private

   type, public :: lu_result
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: u(:, :)
   end type lu_result

   public :: is_square_matrix, is_symmetric_matrix, is_skew_symmetric_matrix
   public :: is_diagonal_matrix, is_idempotent_matrix
   public :: is_singular_matrix, is_non_singular_matrix
   public :: is_positive_definite, is_positive_semi_definite
   public :: is_negative_definite, is_negative_semi_definite, is_indefinite
   public :: lu_decomposition, matrix_inverse, matrix_power, matrix_rank
   public :: spectral_norm, svd_inverse

   interface spectral_norm
      module procedure spectral_norm_matrix
      module procedure spectral_norm_vector
   end interface spectral_norm

contains

   pure logical function is_square_matrix(x) result(answer)
      real(dp), intent(in) :: x(:, :) !! Matrix whose row and column dimensions are compared.

      answer = size(x, 1) == size(x, 2)
   end function is_square_matrix

   pure logical function is_symmetric_matrix(x) result(answer)
      real(dp), intent(in) :: x(:, :) !! Matrix tested for exact symmetry, matching the upstream function.

      if (.not. is_square_matrix(x)) then
         answer = .false.
      else
         answer = all(abs(x - transpose(x)) <= 0.0_dp)
      end if
   end function is_symmetric_matrix

   pure logical function is_skew_symmetric_matrix(x, tol) result(answer)
      real(dp), intent(in) :: x(:, :) !! Square matrix tested for skew symmetry.
      real(dp), intent(in), optional :: tol !! Absolute tolerance; default is 1e-8.
      real(dp) :: tolerance

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      if (.not. is_square_matrix(x)) then
         answer = .false.
      else
         answer = all(abs(x + transpose(x)) < tolerance)
      end if
   end function is_skew_symmetric_matrix

   pure logical function is_diagonal_matrix(x, tol) result(answer)
      real(dp), intent(in) :: x(:, :) !! Square matrix tested for zero off-diagonal entries.
      real(dp), intent(in), optional :: tol !! Absolute off-diagonal tolerance; default is 1e-8.
      real(dp) :: tolerance
      integer :: i
      integer :: j

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      if (.not. is_square_matrix(x)) then
         answer = .false.
         return
      end if
      answer = .true.
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (i /= j .and. abs(x(i, j)) >= tolerance) then
               answer = .false.
               return
            end if
         end do
      end do
   end function is_diagonal_matrix

   pure logical function is_idempotent_matrix(x, tol) result(answer)
      real(dp), intent(in) :: x(:, :) !! Square matrix tested for X*X = X.
      real(dp), intent(in), optional :: tol !! Absolute elementwise tolerance; default is 1e-8.
      real(dp) :: tolerance

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      if (.not. is_square_matrix(x)) then
         answer = .false.
      else
         answer = all(abs(matmul(x, x) - x) < tolerance)
      end if
   end function is_idempotent_matrix

   logical function is_singular_matrix(x, tol, status) result(answer)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose determinant magnitude is tested.
      real(dp), intent(in), optional :: tol !! Determinant tolerance; default is 1e-8.
      integer, intent(out), optional :: status !! Zero on success; nonzero for invalid dimensions.
      real(dp) :: tolerance
      real(dp) :: det_value
      integer :: info

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call determinant_value(x, det_value, info)
      if (present(status)) status = info
      if (info /= 0) then
         answer = .true.
      else
         answer = abs(det_value) < tolerance
      end if
   end function is_singular_matrix

   logical function is_non_singular_matrix(x, tol, status) result(answer)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose determinant magnitude is tested.
      real(dp), intent(in), optional :: tol !! Determinant tolerance; default is 1e-8.
      integer, intent(out), optional :: status !! Zero on success; nonzero for invalid dimensions.
      real(dp) :: tolerance
      real(dp) :: det_value
      integer :: info

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call determinant_value(x, det_value, info)
      if (present(status)) status = info
      if (info /= 0) then
         answer = .false.
      else
         answer = abs(det_value) >= tolerance
      end if
   end function is_non_singular_matrix

   logical function is_positive_definite(x, tol, status) result(answer)
      real(dp), intent(in) :: x(:, :) !! Symmetric matrix whose eigenvalues are tested.
      real(dp), intent(in), optional :: tol !! Eigenvalue zeroing tolerance; default is 1e-8.
      integer, intent(out), optional :: status !! Zero on success; nonzero for nonsymmetry/eigensolver failure.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: tolerance
      integer :: info

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call symmetric_values(x, tolerance, values, vectors, info)
      if (present(status)) status = info
      if (info /= 0) then
         answer = .false.
      else
         answer = all(values > 0.0_dp)
      end if
   end function is_positive_definite

   logical function is_positive_semi_definite(x, tol, status) result(answer)
      real(dp), intent(in) :: x(:, :) !! Symmetric matrix whose eigenvalues are tested.
      real(dp), intent(in), optional :: tol !! Eigenvalue zeroing tolerance; default is 1e-8.
      integer, intent(out), optional :: status !! Zero on success; nonzero for nonsymmetry/eigensolver failure.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: tolerance
      integer :: info

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call symmetric_values(x, tolerance, values, vectors, info)
      if (present(status)) status = info
      if (info /= 0) then
         answer = .false.
      else
         answer = all(values >= 0.0_dp)
      end if
   end function is_positive_semi_definite

   logical function is_negative_definite(x, tol, status) result(answer)
      real(dp), intent(in) :: x(:, :) !! Symmetric matrix whose eigenvalues are tested.
      real(dp), intent(in), optional :: tol !! Eigenvalue zeroing tolerance; default is 1e-8.
      integer, intent(out), optional :: status !! Zero on success; nonzero for nonsymmetry/eigensolver failure.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: tolerance
      integer :: info

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call symmetric_values(x, tolerance, values, vectors, info)
      if (present(status)) status = info
      if (info /= 0) then
         answer = .false.
      else
         answer = all(values < 0.0_dp)
      end if
   end function is_negative_definite

   logical function is_negative_semi_definite(x, tol, status) result(answer)
      real(dp), intent(in) :: x(:, :) !! Symmetric matrix whose eigenvalues are tested.
      real(dp), intent(in), optional :: tol !! Eigenvalue zeroing tolerance; default is 1e-8.
      integer, intent(out), optional :: status !! Zero on success; nonzero for nonsymmetry/eigensolver failure.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: tolerance
      integer :: info

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call symmetric_values(x, tolerance, values, vectors, info)
      if (present(status)) status = info
      if (info /= 0) then
         answer = .false.
      else
         answer = all(values <= 0.0_dp)
      end if
   end function is_negative_semi_definite

   logical function is_indefinite(x, tol, status) result(answer)
      real(dp), intent(in) :: x(:, :) !! Symmetric matrix tested for both positive and negative eigenvalues.
      real(dp), intent(in), optional :: tol !! Eigenvalue zeroing tolerance; default is 1e-8.
      integer, intent(out), optional :: status !! Zero on success; nonzero for nonsymmetry/eigensolver failure.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: tolerance
      integer :: info

      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call symmetric_values(x, tolerance, values, vectors, info)
      if (present(status)) status = info
      if (info /= 0) then
         answer = .false.
      else
         answer = any(values > 0.0_dp) .and. any(values < 0.0_dp)
      end if
   end function is_indefinite

   function lu_decomposition(x, status) result(factors)
      real(dp), intent(in) :: x(:, :) !! Square matrix factored without pivoting, matching upstream Doolittle LU.
      integer, intent(out), optional :: status !! Zero on success; 1 for nonsquare input, 2 for a zero pivot.
      type(lu_result) :: factors
      integer :: n
      integer :: i
      integer :: j
      integer :: k
      real(dp) :: scale
      real(dp) :: pivot_tol
      integer :: info

      n = size(x, 1)
      allocate(factors%l(n, n), factors%u(n, n))
      factors%l = 0.0_dp
      factors%u = 0.0_dp
      info = 0
      if (size(x, 2) /= n) then
         info = 1
         if (present(status)) status = info
         return
      end if
      do i = 1, n
         factors%l(i, i) = 1.0_dp
      end do
      scale = max(1.0_dp, maxval(abs(x)))
      pivot_tol = epsilon(1.0_dp) * scale
      do i = 1, n
         do j = i, n
            factors%u(i, j) = x(i, j)
            do k = 1, i - 1
               factors%u(i, j) = factors%u(i, j) - factors%l(i, k) * factors%u(k, j)
            end do
         end do
         if (i < n) then
            if (abs(factors%u(i, i)) <= pivot_tol) then
               info = 2
               if (present(status)) status = info
               return
            end if
            do j = i + 1, n
               factors%l(j, i) = x(j, i)
               do k = 1, i - 1
                  factors%l(j, i) = factors%l(j, i) - factors%l(j, k) * factors%u(k, i)
               end do
               factors%l(j, i) = factors%l(j, i) / factors%u(i, i)
            end do
         end if
      end do
      if (present(status)) status = info
   end function lu_decomposition

   function matrix_inverse(x, status) result(inv)
      real(dp), intent(in) :: x(:, :) !! Nonsingular square matrix to invert.
      integer, intent(out), optional :: status !! Zero on success; nonzero for invalid dimensions or solve failure.
      real(dp), allocatable :: inv(:, :)
      real(dp), allocatable :: identity(:, :)
      integer :: n
      integer :: i
      integer :: info

      n = size(x, 1)
      allocate(inv(n, n), identity(n, n))
      inv = 0.0_dp
      identity = 0.0_dp
      if (size(x, 2) /= n) then
         if (present(status)) status = 1
         return
      end if
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      call solve_system(x, identity, inv, info)
      if (present(status)) status = info
   end function matrix_inverse

   function matrix_power(x, k, status) result(y)
      real(dp), intent(in) :: x(:, :) !! Square matrix raised to an integer power.
      integer, intent(in) :: k !! Integer exponent; negative values use the matrix inverse.
      integer, intent(out), optional :: status !! Zero on success; nonzero if inversion or dimensions fail.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: base(:, :)
      real(dp), allocatable :: work(:, :)
      integer :: n
      integer :: exponent
      integer :: i
      integer :: info

      n = size(x, 1)
      allocate(y(n, n))
      y = 0.0_dp
      if (size(x, 2) /= n) then
         if (present(status)) status = 1
         return
      end if
      do i = 1, n
         y(i, i) = 1.0_dp
      end do
      info = 0
      if (k < 0) then
         base = matrix_inverse(x, info)
         if (info /= 0) then
            if (present(status)) status = info
            return
         end if
         exponent = -k
      else
         allocate(base(n, n))
         base = x
         exponent = k
      end if
      do while (exponent > 0)
         if (mod(exponent, 2) == 1) then
            work = matmul(y, base)
            y = work
         end if
         exponent = exponent / 2
         if (exponent > 0) then
            work = matmul(base, base)
            base = work
         end if
      end do
      if (present(status)) status = info
   end function matrix_power

   integer function matrix_rank(x, status, tol) result(rank_value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose numerical rank is estimated from Gram-matrix eigenvalues.
      integer, intent(out), optional :: status !! Zero on success; otherwise the symmetric eigensolver status.
      real(dp), intent(in), optional :: tol !! Absolute singular-value cutoff; default is a scale-aware machine tolerance.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp), allocatable :: gram(:, :)
      real(dp) :: cutoff
      real(dp) :: smax
      integer :: info

      if (size(x, 1) >= size(x, 2)) then
         gram = matmul(transpose(x), x)
      else
         gram = matmul(x, transpose(x))
      end if
      call symmetric_eigen(gram, values, vectors, info)
      if (present(status)) status = info
      if (info /= 0 .or. size(values) == 0) then
         rank_value = 0
         return
      end if
      values = max(values, 0.0_dp)
      smax = sqrt(maxval(values))
      cutoff = real(max(size(x, 1), size(x, 2)), dp) * epsilon(1.0_dp) * smax
      if (present(tol)) cutoff = tol
      rank_value = count(sqrt(values) > cutoff)
   end function matrix_rank

   real(dp) function spectral_norm_matrix(x, status) result(value)
      real(dp), intent(in) :: x(:, :) !! Matrix whose largest singular value is returned.
      integer, intent(out), optional :: status !! Zero on success; otherwise the symmetric eigensolver status.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp), allocatable :: gram(:, :)
      integer :: info

      if (size(x, 1) >= size(x, 2)) then
         gram = matmul(transpose(x), x)
      else
         gram = matmul(x, transpose(x))
      end if
      call symmetric_eigen(gram, values, vectors, info)
      if (present(status)) status = info
      if (info /= 0 .or. size(values) == 0) then
         value = 0.0_dp
      else
         value = sqrt(max(0.0_dp, maxval(values)))
      end if
   end function spectral_norm_matrix

   pure real(dp) function spectral_norm_vector(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose Euclidean norm is returned.

      value = sqrt(sum(x * x))
   end function spectral_norm_vector

   function svd_inverse(x, status, tol) result(xinv)
      real(dp), intent(in) :: x(:, :) !! Full-rank matrix whose Moore-Penrose inverse is constructed from a Gram eigendecomposition.
      integer, intent(out), optional :: status !! Zero on success; 2 when a singular value is effectively zero.
      real(dp), intent(in), optional :: tol !! Singular-value cutoff; default is a scale-aware machine tolerance.
      real(dp), allocatable :: xinv(:, :)
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp), allocatable :: gram(:, :)
      real(dp), allocatable :: gram_inv(:, :)
      real(dp) :: cutoff
      real(dp) :: smax
      real(dp) :: lambda
      integer :: info
      integer :: i
      integer :: j
      integer :: q
      integer :: ng

      allocate(xinv(size(x, 2), size(x, 1)))
      xinv = 0.0_dp
      if (size(x, 1) >= size(x, 2)) then
         gram = matmul(transpose(x), x)
      else
         gram = matmul(x, transpose(x))
      end if
      call symmetric_eigen(gram, values, vectors, info)
      if (info /= 0) then
         if (present(status)) status = info
         return
      end if
      values = max(values, 0.0_dp)
      smax = sqrt(maxval(values))
      cutoff = real(max(size(x, 1), size(x, 2)), dp) * epsilon(1.0_dp) * smax
      if (present(tol)) cutoff = tol
      if (any(sqrt(values) <= cutoff)) then
         if (present(status)) status = 2
         return
      end if

      ng = size(values)
      allocate(gram_inv(ng, ng))
      gram_inv = 0.0_dp
      do q = 1, ng
         lambda = values(q)
         do j = 1, ng
            do i = 1, ng
               gram_inv(i, j) = gram_inv(i, j) + vectors(i, q) * vectors(j, q) / lambda
            end do
         end do
      end do
      if (size(x, 1) >= size(x, 2)) then
         xinv = matmul(gram_inv, transpose(x))
      else
         xinv = matmul(transpose(x), gram_inv)
      end if
      if (present(status)) status = 0
   end function svd_inverse

   subroutine symmetric_values(x, tol, values, vectors, info)
      real(dp), intent(in) :: x(:, :) !! Candidate symmetric matrix.
      real(dp), intent(in) :: tol !! Absolute value below which eigenvalues are set to zero.
      real(dp), allocatable, intent(out) :: values(:) !! Eigenvalues after tolerance zeroing.
      real(dp), allocatable, intent(out) :: vectors(:, :) !! Corresponding orthonormal eigenvectors.
      integer, intent(out) :: info !! Zero on success; 1 for nonsymmetry, otherwise eigensolver status.
      integer :: i

      if (.not. is_symmetric_matrix(x)) then
         allocate(values(0), vectors(0, 0))
         info = 1
         return
      end if
      call symmetric_eigen(x, values, vectors, info)
      if (info /= 0) return
      do i = 1, size(values)
         if (abs(values(i)) < tol) values(i) = 0.0_dp
      end do
   end subroutine symmetric_values

   subroutine determinant_value(x, value, info)
      real(dp), intent(in) :: x(:, :) !! Square matrix whose determinant is evaluated by pivoted elimination.
      real(dp), intent(out) :: value !! Determinant value; zero for a numerically singular matrix.
      integer, intent(out) :: info !! Zero on success; 1 for nonsquare input.
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: tmp(:)
      real(dp) :: factor
      real(dp) :: scale
      real(dp) :: pivot_tol
      integer :: n
      integer :: i
      integer :: j
      integer :: pivot
      integer :: sign_value

      n = size(x, 1)
      value = 0.0_dp
      if (size(x, 2) /= n) then
         info = 1
         return
      end if
      allocate(a(n, n), tmp(n))
      a = x
      scale = max(1.0_dp, maxval(abs(a)))
      pivot_tol = tiny(1.0_dp) * scale
      sign_value = 1
      value = 1.0_dp
      do i = 1, n
         pivot = i - 1 + maxloc(abs(a(i:n, i)), dim=1)
         if (abs(a(pivot, i)) <= pivot_tol) then
            value = 0.0_dp
            info = 0
            return
         end if
         if (pivot /= i) then
            tmp = a(i, :)
            a(i, :) = a(pivot, :)
            a(pivot, :) = tmp
            sign_value = -sign_value
         end if
         value = value * a(i, i)
         do j = i + 1, n
            factor = a(j, i) / a(i, i)
            a(j, i:n) = a(j, i:n) - factor * a(i, i:n)
         end do
      end do
      value = real(sign_value, dp) * value
      info = 0
   end subroutine determinant_value

end module matrixcalc_linalg
