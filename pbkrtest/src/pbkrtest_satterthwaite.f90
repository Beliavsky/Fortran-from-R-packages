! SPDX-License-Identifier: GPL-2.0-or-later
module pbkrtest_satterthwaite
   use r_kinds, only : dp
   use r_distributions, only : r_pf
   use r_linalg, only : symmetric_eigen
   use numderiv, only : hessian, jacobian, nd_success
   use pbkrtest_types, only : auxiliary_callbacks_t, auxiliary_result_t, pbkr_invalid_argument, &
      pbkr_invalid_shape, pbkr_linalg_failure, pbkr_numderiv_failure, pbkr_success, &
      satterthwaite_result_t
   use pbkrtest_utils, only : qform
   implicit none
   private
   public :: compute_auxiliary_numeric
   public :: get_fstat_ddf
   public :: satterthwaite_test

contains

   subroutine compute_auxiliary_numeric(varpar_opt, p, callbacks, result, status, tolerance)
      real(dp), intent(in) :: varpar_opt(:) !! Variance parameters where numerical derivatives are evaluated.
      integer, intent(in) :: p !! Number of fixed-effect coefficients, determining the `p` by `p` covariance matrix shape.
      type(auxiliary_callbacks_t), intent(in) :: callbacks !! Deviance and covariance callbacks for `numDeriv`.
      type(auxiliary_result_t), intent(out) :: result !! Numerical derivative matrices and Hessian diagnostics.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), intent(in), optional :: tolerance !! Relative Hessian eigenvalue threshold for mode classification.
      real(dp), allocatable :: eig_values(:)
      real(dp), allocatable :: eig_vectors(:, :)
      real(dp), allocatable :: h_inverse(:, :)
      real(dp), allocatable :: jac(:, :)
      real(dp) :: tol
      integer :: h_status
      integer :: j_status
      integer :: k
      integer :: nvar
      integer :: q

      status = pbkr_success
      if (size(varpar_opt) == 0 .or. p <= 0) then
         status = pbkr_invalid_argument
         return
      end if
      if (.not. associated(callbacks%deviance) .or. .not. associated(callbacks%covbeta_vector)) then
         status = pbkr_invalid_argument
         return
      end if
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      if (tol < 0.0_dp) then
         status = pbkr_invalid_argument
         return
      end if

      call hessian(callbacks%deviance, varpar_opt, result%hessian, status=h_status)
      if (h_status /= nd_success) then
         status = pbkr_numderiv_failure
         return
      end if
      call symmetric_eigen(result%hessian, eig_values, eig_vectors, h_status, descending=.true.)
      if (h_status /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      result%hessian_eigenvalues = eig_values
      result%negative_eigenvalues = count(eig_values < -tol)
      result%near_zero_eigenvalues = count(eig_values > -tol .and. eig_values < tol)
      q = count(eig_values > tol)
      nvar = size(varpar_opt)
      allocate(h_inverse(nvar, nvar))
      h_inverse = 0.0_dp
      do k = 1, q
         h_inverse = h_inverse + outer_product(eig_vectors(:, k), eig_vectors(:, k)) / eig_values(k)
      end do
      allocate(result%vcov_varpar(nvar, nvar))
      result%vcov_varpar = 2.0_dp * h_inverse

      call jacobian(callbacks%covbeta_vector, varpar_opt, jac, status=j_status)
      if (j_status /= nd_success) then
         status = pbkr_numderiv_failure
         return
      end if
      if (size(jac, 1) /= p * p .or. size(jac, 2) /= nvar) then
         status = pbkr_invalid_shape
         return
      end if
      allocate(result%jacobian(p, p, nvar))
      do k = 1, nvar
         result%jacobian(:, :, k) = reshape(jac(:, k), [p, p])
      end do
   end subroutine compute_auxiliary_numeric

   subroutine satterthwaite_test(beta, beta_h, l, vcov_beta, vcov_varpar, jacobian_matrices, &
      result, status, eigen_tolerance)
      real(dp), intent(in) :: beta(:) !! Estimated fixed-effect coefficient vector of length `p`.
      real(dp), intent(in) :: beta_h(:) !! Null-hypothesis coefficient vector of length `p`.
      real(dp), intent(in) :: l(:, :) !! Restriction matrix with shape `(r,p)`.
      real(dp), intent(in) :: vcov_beta(:, :) !! Fixed-effect covariance matrix with shape `(p,p)`.
      real(dp), intent(in) :: vcov_varpar(:, :) !! Variance-parameter covariance matrix with shape `(m,m)`.
      real(dp), intent(in) :: jacobian_matrices(:, :, :) !! Derivatives of `vcov_beta`, with shape `(p,p,m)`.
      type(satterthwaite_result_t), intent(out) :: result !! Satterthwaite statistic, df, p-value, and component df.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), intent(in), optional :: eigen_tolerance !! Relative threshold for contrast-covariance eigenvalues.
      real(dp), allocatable :: beta_diff(:)
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: grad(:, :)
      real(dp), allocatable :: p_vectors(:, :)
      real(dp), allocatable :: ptl(:, :)
      real(dp), allocatable :: t2(:)
      real(dp), allocatable :: v(:, :)
      real(dp) :: eps_value
      real(dp) :: tol
      integer :: i
      integer :: info
      integer :: k
      integer :: m
      integer :: p
      integer :: qq
      integer :: r

      status = pbkr_success
      p = size(beta)
      r = size(l, 1)
      m = size(vcov_varpar, 1)
      if (p <= 0 .or. r <= 0 .or. m <= 0) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(beta_h) /= p .or. size(l, 2) /= p) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(vcov_beta, 1) /= p .or. size(vcov_beta, 2) /= p) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(vcov_varpar, 2) /= m .or. size(jacobian_matrices, 1) /= p .or. &
          size(jacobian_matrices, 2) /= p .or. size(jacobian_matrices, 3) /= m) then
         status = pbkr_invalid_shape
         return
      end if

      allocate(v(r, r))
      v = matmul(matmul(l, vcov_beta), transpose(l))
      call symmetric_eigen(v, d, p_vectors, info, descending=.true.)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      eps_value = sqrt(epsilon(1.0_dp))
      if (present(eigen_tolerance)) eps_value = eigen_tolerance
      tol = max(eps_value * d(1), 0.0_dp)
      qq = count(d > tol)
      if (qq <= 0) then
         status = pbkr_invalid_argument
         return
      end if

      allocate(ptl(qq, p), beta_diff(p), t2(qq), grad(qq, m), result%nu(qq))
      ptl = matmul(transpose(p_vectors(:, 1:qq)), l)
      beta_diff = beta - beta_h
      do i = 1, qq
         t2(i) = dot_product(ptl(i, :), beta_diff) ** 2 / d(i)
      end do
      result%f_stat = sum(t2) / real(qq, dp)

      do i = 1, qq
         do k = 1, m
            grad(i, k) = qform(ptl(i, :), jacobian_matrices(:, :, k))
         end do
         result%nu(i) = 2.0_dp * d(i) ** 2 / qform(grad(i, :), vcov_varpar)
      end do
      result%ndf = qq
      result%ddf = get_fstat_ddf(result%nu)
      result%p_value = r_pf(result%f_stat, real(qq, dp), result%ddf, lower_tail=.false.)
   end subroutine satterthwaite_test

   pure real(dp) function get_fstat_ddf(nu, tolerance) result(ddf)
      real(dp), intent(in) :: nu(:) !! Denominator degrees of freedom from the independent component t-statistics.
      real(dp), intent(in), optional :: tolerance !! Difference threshold for treating all component degrees of freedom as equal.
      real(dp) :: e_value
      real(dp) :: tol
      integer :: n

      n = size(nu)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      if (n == 1) then
         ddf = nu(1)
         return
      end if
      if (n <= 0) then
         ddf = 0.0_dp
         return
      end if
      if (all(abs(nu(2:n) - nu(1:n - 1)) < tol)) then
         ddf = sum(nu) / real(n, dp)
         return
      end if
      if (any(nu <= 2.0_dp)) then
         ddf = 2.0_dp
         return
      end if
      e_value = sum(nu / (nu - 2.0_dp))
      ddf = 2.0_dp * e_value / (e_value - real(n, dp))
   end function get_fstat_ddf

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector of the outer product.
      real(dp), intent(in) :: y(:) !! Right vector of the outer product.
      real(dp) :: a(size(x), size(y))
      integer :: i
      integer :: j

      do j = 1, size(y)
         do i = 1, size(x)
            a(i, j) = x(i) * y(j)
         end do
      end do
   end function outer_product

end module pbkrtest_satterthwaite
