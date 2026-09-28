module optr_model
   use optr_kinds, only : dp
   use optr_status, only : optr_success, optr_invalid_shape, optr_invalid_argument
   use optr_linear, only : gaussian_solve, lu_solve, cholesky_solve, invert_matrix
   use optr_iterative, only : gauss_seidel, conjugate_gradient
   implicit none
   private

   public :: hat_matrix
   public :: optr_fit_design
   public :: optr_solve
   public :: predict_linear

   interface optr_solve
      module procedure optr_solve_vector
      module procedure optr_solve_matrix
   end interface optr_solve

   interface optr_fit_design
      module procedure optr_fit_design_vector
      module procedure optr_fit_design_matrix
   end interface optr_fit_design

   interface predict_linear
      module procedure predict_linear_vector
      module procedure predict_linear_matrix
   end interface predict_linear

contains

   pure subroutine optr_solve_vector(a, b, method, x, status, tol, max_iter)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A x = b`.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      character(len=*), intent(in) :: method !! Solver name: `gauss`, `LU`, `gaussseidel`, `cgm`, `choleski`, or `cholesky`.
      real(dp), allocatable, intent(out) :: x(:) !! Solution or final iterate.
      integer, intent(out) :: status !! Zero on success; nonzero on invalid input, factorization failure, or nonconvergence.
      real(dp), optional, intent(in) :: tol !! Solver tolerance; defaults to `1.0e-7_dp`.
      integer, optional, intent(in) :: max_iter !! Iteration cap for iterative methods; defaults to 500.

      real(dp), allocatable :: convergence(:)
      real(dp) :: local_tol
      integer :: iterations
      integer :: local_max_iter

      local_tol = 1.0e-7_dp
      if (present(tol)) local_tol = tol
      local_max_iter = 500
      if (present(max_iter)) local_max_iter = max_iter

      select case (trim(method))
      case ("gauss", "GAUSS")
         call gaussian_solve(a, b, x, status, local_tol)
      case ("LU", "lu", "Lu", "lU")
         call lu_solve(a, b, x, status, local_tol)
      case ("gaussseidel", "GAUSSSEIDEL")
         call gauss_seidel(a, b, x, status, iterations, tol=local_tol, max_iter=local_max_iter)
      case ("cgm", "CGM")
         call conjugate_gradient(a, b, x, status, iterations, convergence, tol=local_tol, max_iter=local_max_iter)
      case ("choleski", "CHOLESKI", "cholesky", "CHOLESKY")
         call cholesky_solve(a, b, x, status, local_tol)
      case default
         status = optr_invalid_argument
         allocate(x(0))
      end select
   end subroutine optr_solve_vector

   pure subroutine optr_solve_matrix(a, b, method, x, status, tol, max_iter)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A X = B`.
      real(dp), intent(in) :: b(:, :) !! Right-hand-side matrix with one or more columns.
      character(len=*), intent(in) :: method !! Solver name: `gauss`, `LU`, `gaussseidel`, `cgm`, `choleski`, or `cholesky`.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix or final iterates.
      integer, intent(out) :: status !! Zero on success; nonzero on invalid input, factorization failure, or nonconvergence.
      real(dp), optional, intent(in) :: tol !! Solver tolerance; defaults to `1.0e-7_dp`.
      integer, optional, intent(in) :: max_iter !! Iteration cap for iterative methods; defaults to 500.

      real(dp), allocatable :: column(:)
      real(dp) :: local_tol
      integer :: j
      integer :: local_max_iter
      integer :: n

      if (size(b, 1) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0, 0))
         return
      end if

      local_tol = 1.0e-7_dp
      if (present(tol)) local_tol = tol
      local_max_iter = 500
      if (present(max_iter)) local_max_iter = max_iter

      select case (trim(method))
      case ("gauss", "GAUSS")
         call gaussian_solve(a, b, x, status, local_tol)
      case ("LU", "lu", "Lu", "lU")
         call lu_solve(a, b, x, status, local_tol)
      case ("choleski", "CHOLESKI", "cholesky", "CHOLESKY")
         call cholesky_solve(a, b, x, status, local_tol)
      case ("gaussseidel", "GAUSSSEIDEL", "cgm", "CGM")
         n = size(a, 1)
         allocate(x(n, size(b, 2)))
         do j = 1, size(b, 2)
            call optr_solve_vector(a, b(:, j), method, column, status, local_tol, local_max_iter)
            if (status /= optr_success) return
            x(:, j) = column
         end do
      case default
         status = optr_invalid_argument
         allocate(x(0, 0))
      end select
   end subroutine optr_solve_matrix

   pure subroutine optr_fit_design_matrix(design, response, beta, status, method, weights, tol, max_iter)
      real(dp), intent(in) :: design(:, :) !! Numeric design matrix with observations in rows and coefficients in columns.
      real(dp), intent(in) :: response(:, :) !! Response matrix with the same number of observations as `design`.
      real(dp), allocatable, intent(out) :: beta(:, :) !! Fitted coefficient matrix, one column per response.
      integer, intent(out) :: status !! Zero on success; nonzero on invalid weights/shapes or solver failure.
      character(len=*), optional, intent(in) :: method !! Solver used for square or normal-equation systems; defaults to `LU`.
      real(dp), optional, intent(in) :: weights(:) !! Nonnegative observation weights; omitted for unweighted fitting.
      real(dp), optional, intent(in) :: tol !! Solver tolerance; defaults to `1.0e-7_dp`.
      integer, optional, intent(in) :: max_iter !! Iteration cap for iterative methods; defaults to 500.

      real(dp), allocatable :: lhs(:, :)
      real(dp), allocatable :: response_weighted(:, :)
      real(dp), allocatable :: rhs(:, :)
      real(dp), allocatable :: design_weighted(:, :)
      real(dp) :: local_tol
      real(dp) :: sw
      character(len=16) :: local_method
      integer :: i
      integer :: local_max_iter
      integer :: n
      integer :: p

      n = size(design, 1)
      p = size(design, 2)
      if (size(response, 1) /= n .or. n < 1 .or. p < 1) then
         status = optr_invalid_shape
         allocate(beta(0, 0))
         return
      end if

      local_method = "LU"
      if (present(method)) local_method = method
      local_tol = 1.0e-7_dp
      if (present(tol)) local_tol = tol
      local_max_iter = 500
      if (present(max_iter)) local_max_iter = max_iter

      allocate(design_weighted(n, p), response_weighted(n, size(response, 2)))
      design_weighted = design
      response_weighted = response
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = optr_invalid_argument
            allocate(beta(0, 0))
            return
         end if
         do i = 1, n
            sw = sqrt(weights(i))
            design_weighted(i, :) = sw * design_weighted(i, :)
            response_weighted(i, :) = sw * response_weighted(i, :)
         end do
      end if

      if (n == p) then
         call optr_solve_matrix(design_weighted, response_weighted, trim(local_method), beta, status, &
            local_tol, local_max_iter)
      else
         allocate(lhs(p, p), rhs(p, size(response, 2)))
         lhs = matmul(transpose(design_weighted), design_weighted)
         rhs = matmul(transpose(design_weighted), response_weighted)
         call optr_solve_matrix(lhs, rhs, trim(local_method), beta, status, local_tol, local_max_iter)
      end if
   end subroutine optr_fit_design_matrix

   pure subroutine optr_fit_design_vector(design, response, beta, status, method, weights, tol, max_iter)
      real(dp), intent(in) :: design(:, :) !! Numeric design matrix with observations in rows and coefficients in columns.
      real(dp), intent(in) :: response(:) !! Response vector with one element per observation.
      real(dp), allocatable, intent(out) :: beta(:) !! Fitted coefficient vector.
      integer, intent(out) :: status !! Zero on success; nonzero on invalid weights/shapes or solver failure.
      character(len=*), optional, intent(in) :: method !! Solver used for square or normal-equation systems; defaults to `LU`.
      real(dp), optional, intent(in) :: weights(:) !! Nonnegative observation weights; omitted for unweighted fitting.
      real(dp), optional, intent(in) :: tol !! Solver tolerance; defaults to `1.0e-7_dp`.
      integer, optional, intent(in) :: max_iter !! Iteration cap for iterative methods; defaults to 500.

      real(dp), allocatable :: beta_matrix(:, :)
      real(dp), allocatable :: response_matrix(:, :)
      real(dp) :: local_tol
      character(len=16) :: local_method
      integer :: local_max_iter

      if (size(response) /= size(design, 1)) then
         status = optr_invalid_shape
         allocate(beta(0))
         return
      end if
      allocate(response_matrix(size(response), 1))
      response_matrix(:, 1) = response
      local_method = "LU"
      if (present(method)) local_method = method
      local_tol = 1.0e-7_dp
      if (present(tol)) local_tol = tol
      local_max_iter = 500
      if (present(max_iter)) local_max_iter = max_iter

      if (present(weights)) then
         call optr_fit_design_matrix(design, response_matrix, beta_matrix, status, trim(local_method), &
            weights, local_tol, local_max_iter)
      else
         call optr_fit_design_matrix(design, response_matrix, beta_matrix, status, trim(local_method), &
            tol=local_tol, max_iter=local_max_iter)
      end if
      if (status /= optr_success) then
         allocate(beta(0))
         return
      end if
      allocate(beta(size(beta_matrix, 1)))
      beta = beta_matrix(:, 1)
   end subroutine optr_fit_design_vector

   pure subroutine predict_linear_vector(design, beta, fitted, status)
      real(dp), intent(in) :: design(:, :) !! New design matrix with the same number of columns as coefficients.
      real(dp), intent(in) :: beta(:) !! Coefficient vector.
      real(dp), allocatable, intent(out) :: fitted(:) !! Predicted response vector `design * beta`.
      integer, intent(out) :: status !! Zero on success; nonzero when matrix and coefficient dimensions differ.

      if (size(design, 2) /= size(beta)) then
         status = optr_invalid_shape
         allocate(fitted(0))
         return
      end if
      allocate(fitted(size(design, 1)))
      fitted = matmul(design, beta)
      status = optr_success
   end subroutine predict_linear_vector

   pure subroutine predict_linear_matrix(design, beta, fitted, status)
      real(dp), intent(in) :: design(:, :) !! New design matrix with the same number of columns as coefficient rows.
      real(dp), intent(in) :: beta(:, :) !! Coefficient matrix with one column per response.
      real(dp), allocatable, intent(out) :: fitted(:, :) !! Predicted response matrix `design * beta`.
      integer, intent(out) :: status !! Zero on success; nonzero when matrix and coefficient dimensions differ.

      if (size(design, 2) /= size(beta, 1)) then
         status = optr_invalid_shape
         allocate(fitted(0, 0))
         return
      end if
      allocate(fitted(size(design, 1), size(beta, 2)))
      fitted = matmul(design, beta)
      status = optr_success
   end subroutine predict_linear_matrix

   pure subroutine hat_matrix(design, h, status, covmat, tol)
      real(dp), intent(in) :: design(:, :) !! Design matrix `X` with observations in rows.
      real(dp), allocatable, intent(out) :: h(:, :) !! Projection/hat matrix, generalized when `covmat` is supplied.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or a singular covariance/Gram matrix.
      real(dp), optional, intent(in) :: covmat(:, :) !! Observation covariance matrix; omitted for ordinary least squares.
      real(dp), optional, intent(in) :: tol !! Inversion pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: cov_inverse(:, :)
      real(dp), allocatable :: gram(:, :)
      real(dp), allocatable :: gram_inverse(:, :)
      real(dp) :: local_tol
      integer :: n

      n = size(design, 1)
      if (n < 1 .or. size(design, 2) < 1) then
         status = optr_invalid_shape
         allocate(h(0, 0))
         return
      end if
      local_tol = 1.0e-7_dp
      if (present(tol)) local_tol = tol

      if (present(covmat)) then
         if (size(covmat, 1) /= n .or. size(covmat, 2) /= n) then
            status = optr_invalid_shape
            allocate(h(0, 0))
            return
         end if
         call invert_matrix(covmat, cov_inverse, status, local_tol)
         if (status /= optr_success) then
            allocate(h(0, 0))
            return
         end if
         allocate(gram(size(design, 2), size(design, 2)))
         gram = matmul(transpose(design), matmul(cov_inverse, design))
         call invert_matrix(gram, gram_inverse, status, local_tol)
         if (status /= optr_success) then
            allocate(h(0, 0))
            return
         end if
         allocate(h(n, n))
         h = matmul(design, matmul(gram_inverse, matmul(transpose(design), cov_inverse)))
      else
         allocate(gram(size(design, 2), size(design, 2)))
         gram = matmul(transpose(design), design)
         call invert_matrix(gram, gram_inverse, status, local_tol)
         if (status /= optr_success) then
            allocate(h(0, 0))
            return
         end if
         allocate(h(n, n))
         h = matmul(design, matmul(gram_inverse, transpose(design)))
      end if
      status = optr_success
   end subroutine hat_matrix

end module optr_model
