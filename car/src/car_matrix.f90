module car_matrix
   use car_kinds, only: dp
   use car_utils, only: matrix_determinant, matrix_inverse, chi_square_survival, f_survival
   implicit none
   private

   public :: delta_method_gradient
   public :: gvif_from_covariance
   public :: hccm_matrix
   public :: linear_hypothesis_test
   public :: sigma_hat
   public :: wcrossprod
   public :: wcrossprod_matrix_weight
   public :: wcrossprod_vector_weight

contains

   pure function wcrossprod(x, y) result(result_matrix)
      real(dp), intent(in) :: x(:, :) !! Left data matrix with observations in rows.
      real(dp), intent(in), optional :: y(:, :) !! Optional right data matrix with the same number of rows as x.
      real(dp), allocatable :: result_matrix(:, :)

      if (present(y)) then
         if (size(y, 1) /= size(x, 1)) then
            allocate(result_matrix(0, 0))
            return
         end if
         result_matrix = matmul(transpose(x), y)
      else
         result_matrix = matmul(transpose(x), x)
      end if
   end function wcrossprod

   pure function wcrossprod_vector_weight(x, w, y) result(result_matrix)
      real(dp), intent(in) :: x(:, :) !! Left data matrix with observations in rows.
      real(dp), intent(in) :: w(:) !! Observation weights; length must equal the number of rows of x.
      real(dp), intent(in), optional :: y(:, :) !! Optional right data matrix with the same number of rows as x.
      real(dp), allocatable :: result_matrix(:, :)
      real(dp), allocatable :: wy(:, :)
      integer :: j

      if (size(w) /= size(x, 1)) then
         allocate(result_matrix(0, 0))
         return
      end if
      if (present(y)) then
         if (size(y, 1) /= size(x, 1)) then
            allocate(result_matrix(0, 0))
            return
         end if
         allocate(wy(size(y, 1), size(y, 2)))
         do j = 1, size(y, 2)
            wy(:, j) = w * y(:, j)
         end do
         result_matrix = matmul(transpose(x), wy)
      else
         allocate(wy(size(x, 1), size(x, 2)))
         do j = 1, size(x, 2)
            wy(:, j) = w * x(:, j)
         end do
         result_matrix = matmul(transpose(x), wy)
      end if
   end function wcrossprod_vector_weight

   pure function wcrossprod_matrix_weight(x, w, y) result(result_matrix)
      real(dp), intent(in) :: x(:, :) !! Left data matrix with observations in rows.
      real(dp), intent(in) :: w(:, :) !! Square observation-weight matrix with order equal to the number of rows of x.
      real(dp), intent(in), optional :: y(:, :) !! Optional right data matrix with the same number of rows as x.
      real(dp), allocatable :: result_matrix(:, :)

      if (size(w, 1) /= size(w, 2) .or. size(w, 1) /= size(x, 1)) then
         allocate(result_matrix(0, 0))
         return
      end if
      if (present(y)) then
         if (size(y, 1) /= size(x, 1)) then
            allocate(result_matrix(0, 0))
            return
         end if
         result_matrix = matmul(transpose(x), matmul(w, y))
      else
         result_matrix = matmul(transpose(x), matmul(w, x))
      end if
   end function wcrossprod_matrix_weight

   pure function sigma_hat(residuals, residual_df) result(value)
      real(dp), intent(in) :: residuals(:) !! Model residuals used to estimate the residual standard deviation.
      integer, intent(in) :: residual_df !! Positive residual degrees of freedom.
      real(dp) :: value

      if (residual_df <= 0) then
         value = -1.0_dp
      else
         value = sqrt(sum(residuals**2) / real(residual_df, dp))
      end if
   end function sigma_hat

   subroutine hccm_matrix(x, residuals, cov_unscaled, leverage, covariance, status, hc_type, weights, residual_df)
      real(dp), intent(in) :: x(:, :) !! Model matrix after removal of aliased columns, with observations in rows.
      real(dp), intent(in) :: residuals(:) !! Regression residuals corresponding to rows of x.
      real(dp), intent(in) :: cov_unscaled(:, :) !! Unscaled coefficient covariance matrix, ordinarily inverse(X'WX).
      real(dp), intent(in) :: leverage(:) !! Hat values for observations, each normally in [0, 1).
      real(dp), intent(out) :: covariance(:, :) !! Heteroscedasticity-consistent coefficient covariance matrix.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, HC type, or leverage.
      character(len=*), intent(in), optional :: hc_type !! Estimator name: hc0, hc1, hc2, hc3, or hc4; default hc3.
      real(dp), intent(in), optional :: weights(:) !! Optional prior regression weights; default is one for every row.
      integer, intent(in), optional :: residual_df !! Residual degrees of freedom required by HC1; default n-p.
      character(len=8) :: kind
      real(dp), allocatable :: e(:), factor(:), meat(:, :), scaled_x(:, :)
      integer :: i, n, p, j, df

      status = 0
      covariance = 0.0_dp
      n = size(x, 1)
      p = size(x, 2)
      if (size(residuals) /= n .or. size(leverage) /= n) then
         status = 1
         return
      end if
      if (size(cov_unscaled, 1) /= p .or. size(cov_unscaled, 2) /= p) then
         status = 1
         return
      end if
      if (size(covariance, 1) /= p .or. size(covariance, 2) /= p) then
         status = 1
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n) then
            status = 1
            return
         end if
      end if

      kind = "hc3"
      if (present(hc_type)) kind = lowercase(adjustl(hc_type))
      allocate(e(n), factor(n), scaled_x(n, p))
      e = residuals
      if (present(weights)) e = e * weights

      select case (trim(kind))
      case ("hc0")
         factor = 1.0_dp
      case ("hc1")
         df = n - p
         if (present(residual_df)) df = residual_df
         if (df <= 0) then
            status = 2
            return
         end if
         factor = real(df, dp) / real(n, dp)
      case ("hc2")
         if (any(leverage >= 1.0_dp)) then
            status = 3
            return
         end if
         factor = 1.0_dp - leverage
      case ("hc3")
         if (any(leverage >= 1.0_dp)) then
            status = 3
            return
         end if
         factor = (1.0_dp - leverage)**2
      case ("hc4")
         if (p <= 0) then
            status = 2
            return
         end if
         do i = 1, n
            factor(i) = (1.0_dp - leverage(i))**min(4.0_dp, real(n, dp) * leverage(i) / real(p, dp))
            if (factor(i) <= tiny(1.0_dp)) then
               status = 3
               return
            end if
         end do
      case default
         status = 4
         return
      end select

      do j = 1, p
         scaled_x(:, j) = x(:, j) * (e**2 / factor)
      end do
      meat = matmul(transpose(x), scaled_x)
      covariance = matmul(cov_unscaled, matmul(meat, cov_unscaled))
   end subroutine hccm_matrix

   subroutine gvif_from_covariance(vcov, assign, gvif, degrees_freedom, adjusted_gvif, status)
      real(dp), intent(in) :: vcov(:, :) !! Coefficient covariance matrix with intercept already removed when applicable.
      integer, intent(in) :: assign(:) !! Positive term number assigned to each covariance-matrix coefficient.
      real(dp), allocatable, intent(out) :: gvif(:) !! Generalized variance-inflation factor for each term number.
      integer, allocatable, intent(out) :: degrees_freedom(:) !! Number of coefficients assigned to each term.
      real(dp), allocatable, intent(out) :: adjusted_gvif(:) !! GVIF raised to 1/(2*df) for each term.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid input or singular correlation matrices.
      real(dp), allocatable :: corr(:, :), sub(:, :), comp(:, :)
      real(dp), allocatable :: scale(:)
      logical, allocatable :: keep(:)
      real(dp) :: det_all, det_sub, det_comp
      integer :: i, j, k, n, n_terms, ns, nc, det_status

      status = 0
      n = size(vcov, 1)
      if (size(vcov, 2) /= n .or. size(assign) /= n .or. n == 0 .or. any(assign <= 0)) then
         status = 1
         allocate(gvif(0), degrees_freedom(0), adjusted_gvif(0))
         return
      end if
      allocate(scale(n), corr(n, n))
      do i = 1, n
         if (vcov(i, i) <= 0.0_dp) then
            status = 2
            allocate(gvif(0), degrees_freedom(0), adjusted_gvif(0))
            return
         end if
         scale(i) = sqrt(vcov(i, i))
      end do
      do j = 1, n
         do i = 1, n
            corr(i, j) = vcov(i, j) / (scale(i) * scale(j))
         end do
      end do
      call matrix_determinant(corr, det_all, det_status)
      if (det_status /= 0 .or. abs(det_all) <= 100.0_dp * epsilon(1.0_dp)) then
         status = 3
         allocate(gvif(0), degrees_freedom(0), adjusted_gvif(0))
         return
      end if

      n_terms = maxval(assign)
      allocate(gvif(n_terms), degrees_freedom(n_terms), adjusted_gvif(n_terms), keep(n))
      gvif = 0.0_dp
      degrees_freedom = 0
      adjusted_gvif = 0.0_dp
      do k = 1, n_terms
         keep = assign == k
         ns = count(keep)
         nc = n - ns
         if (ns == 0) then
            status = 4
            return
         end if
         degrees_freedom(k) = ns
         allocate(sub(ns, ns))
         call extract_principal(corr, keep, .true., sub)
         call matrix_determinant(sub, det_sub, det_status)
         if (det_status /= 0) then
            status = 3
            return
         end if
         deallocate(sub)
         if (nc > 0) then
            allocate(comp(nc, nc))
            call extract_principal(corr, keep, .false., comp)
            call matrix_determinant(comp, det_comp, det_status)
            if (det_status /= 0) then
               status = 3
               return
            end if
            deallocate(comp)
         else
            det_comp = 1.0_dp
         end if
         gvif(k) = det_sub * det_comp / det_all
         adjusted_gvif(k) = gvif(k)**(1.0_dp / (2.0_dp * real(ns, dp)))
      end do
   end subroutine gvif_from_covariance

   pure subroutine linear_hypothesis_test(beta, vcov, hypothesis, rhs, statistic, p_value, status, residual_df, &
                                     estimated_value, hypothesis_vcov)
      real(dp), intent(in) :: beta(:) !! Estimated coefficient vector.
      real(dp), intent(in) :: vcov(:, :) !! Covariance matrix for beta.
      real(dp), intent(in) :: hypothesis(:, :) !! Rows of the numeric linear-hypothesis matrix L.
      real(dp), intent(in) :: rhs(:) !! Right-hand-side vector in hypotheses L beta = rhs.
      real(dp), intent(out) :: statistic !! Chi-square statistic, or F statistic when residual_df is supplied.
      real(dp), intent(out) :: p_value !! Upper-tail probability for the returned statistic.
      integer, intent(out) :: status !! Zero on success; nonzero for bad dimensions or singular hypothesis covariance.
      integer, intent(in), optional :: residual_df !! Positive residual degrees of freedom selects an F test.
      real(dp), allocatable, intent(out), optional :: estimated_value(:) !! Optional estimated L beta - rhs vector.
      real(dp), allocatable, intent(out), optional :: hypothesis_vcov(:, :) !! Optional covariance L V L' of the hypotheses.
      real(dp), allocatable :: value(:), hv(:, :), hv_inv(:, :)
      integer :: q, p, inv_status
      real(dp) :: ssh

      statistic = 0.0_dp
      p_value = 1.0_dp
      status = 0
      p = size(beta)
      q = size(hypothesis, 1)
      if (size(vcov, 1) /= p .or. size(vcov, 2) /= p .or. size(hypothesis, 2) /= p .or. &
          size(rhs) /= q .or. q <= 0) then
         status = 1
         return
      end if
      value = matmul(hypothesis, beta) - rhs
      hv = matmul(hypothesis, matmul(vcov, transpose(hypothesis)))
      allocate(hv_inv(q, q))
      call matrix_inverse(hv, hv_inv, inv_status)
      if (inv_status /= 0) then
         status = 2
         return
      end if
      ssh = dot_product(value, matmul(hv_inv, value))
      if (present(residual_df)) then
         if (residual_df > 0) then
            statistic = ssh / real(q, dp)
            p_value = f_survival(statistic, q, residual_df)
         else
            status = 3
            return
         end if
      else
         statistic = ssh
         p_value = chi_square_survival(statistic, q)
      end if
      if (present(estimated_value)) estimated_value = value
      if (present(hypothesis_vcov)) hypothesis_vcov = hv
   end subroutine linear_hypothesis_test

   subroutine delta_method_gradient(estimate, gradient, vcov, transformed_estimate, standard_error, status)
      real(dp), intent(in) :: estimate !! Value of a scalar transformed parameter function evaluated at the estimates.
      real(dp), intent(in) :: gradient(:) !! Gradient of the transformed parameter function with respect to coefficients.
      real(dp), intent(in) :: vcov(:, :) !! Covariance matrix of the original coefficient estimates.
      real(dp), intent(out) :: transformed_estimate !! Returned transformed estimate, equal to the supplied estimate.
      real(dp), intent(out) :: standard_error !! Delta-method standard error sqrt(g' V g).
      integer, intent(out) :: status !! Zero on success; nonzero for dimension mismatch or negative numerical variance.
      real(dp) :: variance

      transformed_estimate = estimate
      standard_error = 0.0_dp
      status = 0
      if (size(vcov, 1) /= size(gradient) .or. size(vcov, 2) /= size(gradient)) then
         status = 1
         return
      end if
      variance = dot_product(gradient, matmul(vcov, gradient))
      if (variance < -100.0_dp * epsilon(1.0_dp)) then
         status = 2
         return
      end if
      standard_error = sqrt(max(0.0_dp, variance))
   end subroutine delta_method_gradient

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text !! Text converted to lowercase for option matching.
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

   pure subroutine extract_principal(matrix, mask, selected, result_matrix)
      real(dp), intent(in) :: matrix(:, :) !! Square source matrix from which a principal submatrix is extracted.
      logical, intent(in) :: mask(:) !! Mask identifying the selected set of rows and columns.
      logical, intent(in) :: selected !! True extracts masked entries; false extracts their complement.
      real(dp), intent(out) :: result_matrix(:, :) !! Principal submatrix in original coefficient order.
      integer :: i, j, ii, jj

      ii = 0
      do i = 1, size(mask)
         if (mask(i) .eqv. selected) then
            ii = ii + 1
            jj = 0
            do j = 1, size(mask)
               if (mask(j) .eqv. selected) then
                  jj = jj + 1
                  result_matrix(ii, jj) = matrix(i, j)
               end if
            end do
         end if
      end do
   end subroutine extract_principal

end module car_matrix
