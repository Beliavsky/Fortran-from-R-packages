! SPDX-License-Identifier: GPL-2.0-or-later
module pbkrtest_kr
   use r_kinds, only : dp
   use r_distributions, only : r_pf
   use r_linalg, only : inverse_matrix, numerical_rank, solve_system, spd_inverse_logdet, &
      symmetric_eigen, symmetric_eigenvalues, symmetrize
   use pbkrtest_types, only : kr_result_t, pbkr_invalid_argument, pbkr_invalid_shape, &
      pbkr_linalg_failure, pbkr_success, vcov_adjustment_t
   use pbkrtest_utils, only : div_zero, pair_index_upper, trace_matrix
   implicit none
   private
   public :: ddf_lb_scalar
   public :: kr_adjust
   public :: lb_ddf
   public :: vcov_adjust_kr

contains

   subroutine vcov_adjust_kr(phi, sigma, g, x, result, status)
      real(dp), intent(in) :: phi(:, :) !! Unadjusted covariance matrix of fixed effects, with shape `(p,p)`.
      real(dp), intent(in) :: sigma(:, :) !! Marginal response covariance matrix, with shape `(n,n)` and positive definite.
      real(dp), intent(in) :: g(:, :, :) !! Covariance-component derivative matrices, with shape `(n,n,m)`.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect model matrix, with shape `(n,p)`.
      type(vcov_adjustment_t), intent(out) :: result !! Adjusted covariance and KR auxiliary matrices.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: hh(:, :, :)
      real(dp), allocatable :: ie_inverse(:, :)
      real(dp), allocatable :: ktrace(:, :)
      real(dp), allocatable :: oo(:, :, :)
      real(dp), allocatable :: phi_p_i(:, :)
      real(dp), allocatable :: qq(:, :, :)
      real(dp), allocatable :: sigma_inverse(:, :)
      real(dp), allocatable :: tt(:, :)
      real(dp), allocatable :: uu(:, :)
      real(dp) :: ignored_logdet
      integer :: i
      integer :: info
      integer :: j
      integer :: m
      integer :: n
      integer :: nq
      integer :: p
      integer :: qindex

      status = pbkr_success
      n = size(sigma, 1)
      p = size(phi, 1)
      m = size(g, 3)
      if (n <= 0 .or. p <= 0 .or. m <= 0) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(sigma, 2) /= n .or. size(x, 1) /= n .or. size(x, 2) /= p) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(phi, 2) /= p .or. size(g, 1) /= n .or. size(g, 2) /= n) then
         status = pbkr_invalid_shape
         return
      end if

      call spd_inverse_logdet(sigma, sigma_inverse, ignored_logdet, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if

      nq = m * (m + 1) / 2
      allocate(tt(n, p), hh(n, n, m), oo(n, p, m))
      allocate(result%p_matrices(p, p, m), qq(p, p, nq))
      allocate(ktrace(m, m), result%information(m, m), uu(p, p))
      tt = matmul(sigma_inverse, x)
      do i = 1, m
         hh(:, :, i) = matmul(g(:, :, i), sigma_inverse)
         oo(:, :, i) = matmul(hh(:, :, i), x)
      end do

      do i = 1, m
         result%p_matrices(:, :, i) = symmetrize(-matmul(transpose(oo(:, :, i)), tt))
         do j = i, m
            qindex = pair_index_upper(i, j, m)
            qq(:, :, qindex) = matmul(matmul(transpose(oo(:, :, i)), sigma_inverse), oo(:, :, j))
         end do
      end do

      ktrace = 0.0_dp
      do i = 1, m
         do j = i, m
            ktrace(i, j) = sum(transpose(hh(:, :, i)) * hh(:, :, j))
            ktrace(j, i) = ktrace(i, j)
         end do
      end do

      allocate(phi_p_i(p, p))
      result%information = 0.0_dp
      do i = 1, m
         phi_p_i = matmul(phi, result%p_matrices(:, :, i))
         do j = i, m
            qindex = pair_index_upper(i, j, m)
            result%information(i, j) = ktrace(i, j) - &
               2.0_dp * sum(phi * qq(:, :, qindex)) + &
               sum(phi_p_i * matmul(result%p_matrices(:, :, j), phi))
            result%information(j, i) = result%information(i, j)
         end do
      end do

      call symmetric_eigenvalues(result%information, eigenvalues, info)
      if (info /= 0 .or. size(eigenvalues) /= m) then
         status = pbkr_linalg_failure
         return
      end if
      result%condition = minval(abs(eigenvalues))
      if (result%condition > 1.0e-10_dp) then
         call inverse_matrix(result%information, ie_inverse, info)
      else
         call symmetric_pseudoinverse(result%information, ie_inverse, info)
      end if
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      allocate(result%w(m, m))
      result%w = symmetrize(2.0_dp * ie_inverse)

      uu = 0.0_dp
      do i = 1, m - 1
         do j = i + 1, m
            qindex = pair_index_upper(i, j, m)
            uu = uu + result%w(i, j) * &
               (qq(:, :, qindex) - matmul(matmul(result%p_matrices(:, :, i), phi), &
               result%p_matrices(:, :, j)))
         end do
      end do
      uu = uu + transpose(uu)
      do i = 1, m
         qindex = pair_index_upper(i, i, m)
         uu = uu + result%w(i, i) * &
            (qq(:, :, qindex) - matmul(matmul(result%p_matrices(:, :, i), phi), &
            result%p_matrices(:, :, i)))
      end do

      allocate(result%phi_adjusted(p, p))
      result%phi_adjusted = phi + 2.0_dp * matmul(matmul(phi, uu), phi)
      result%phi_adjusted = symmetrize(result%phi_adjusted)
   end subroutine vcov_adjust_kr

   subroutine kr_adjust(phi_adjusted, phi, l, beta, beta_h, p_matrices, w, result, status)
      real(dp), intent(in) :: phi_adjusted(:, :) !! Kenward-Roger adjusted fixed-effect covariance matrix, shape `(p,p)`.
      real(dp), intent(in) :: phi(:, :) !! Unadjusted fixed-effect covariance matrix, shape `(p,p)`.
      real(dp), intent(in) :: l(:, :) !! Full-row-rank restriction matrix, shape `(r,p)`.
      real(dp), intent(in) :: beta(:) !! Estimated fixed-effect coefficients, length `p`.
      real(dp), intent(in) :: beta_h(:) !! Null-hypothesis coefficient vector subtracted from `beta`, length `p`.
      real(dp), intent(in) :: p_matrices(:, :, :) !! `P` matrices returned by `vcov_adjust_kr`, shape `(p,p,m)`.
      real(dp), intent(in) :: w(:, :) !! Covariance-parameter variance matrix returned by `vcov_adjust_kr`, shape `(m,m)`.
      type(kr_result_t), intent(out) :: result !! Kenward-Roger numerator/denominator df, F statistics, p-values, and auxiliaries.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: beta_diff(:)
      real(dp), allocatable :: lb2(:)
      real(dp), allocatable :: solution(:, :)
      real(dp), allocatable :: solution_vector(:)
      real(dp), allocatable :: theta(:, :)
      real(dp) :: b_value
      real(dp) :: c1
      real(dp) :: c2
      real(dp) :: c3
      real(dp) :: g_value
      real(dp) :: qreal
      real(dp) :: wald
      integer :: info
      integer :: m
      integer :: p
      integer :: q
      integer :: rows_l

      status = pbkr_success
      p = size(phi, 1)
      rows_l = size(l, 1)
      m = size(p_matrices, 3)
      if (p <= 0 .or. rows_l <= 0 .or. m <= 0) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(phi, 2) /= p .or. size(phi_adjusted, 1) /= p .or. &
          size(phi_adjusted, 2) /= p) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(l, 2) /= p .or. size(beta) /= p .or. size(beta_h) /= p) then
         status = pbkr_invalid_shape
         return
      end if
      if (size(p_matrices, 1) /= p .or. size(p_matrices, 2) /= p .or. &
          size(w, 1) /= m .or. size(w, 2) /= m) then
         status = pbkr_invalid_shape
         return
      end if

      call numerical_rank(l, q, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      if (q <= 0 .or. q /= rows_l) then
         status = pbkr_invalid_argument
         return
      end if
      qreal = real(q, dp)

      allocate(a(q, q), solution(q, p), theta(p, p))
      a = matmul(matmul(l, phi), transpose(l))
      call solve_system(a, l, solution, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      theta = matmul(transpose(l), solution)
      call accumulate_a1_a2(theta, phi, p_matrices, w, result%a1, result%a2)

      b_value = (result%a1 + 6.0_dp * result%a2) / (2.0_dp * qreal)
      g_value = ((qreal + 1.0_dp) * result%a1 - (qreal + 4.0_dp) * result%a2) / &
         ((qreal + 2.0_dp) * result%a2)
      c1 = g_value / (3.0_dp * qreal + 2.0_dp * (1.0_dp - g_value))
      c2 = (qreal - g_value) / (3.0_dp * qreal + 2.0_dp * (1.0_dp - g_value))
      c3 = (qreal + 2.0_dp - g_value) / (3.0_dp * qreal + 2.0_dp * (1.0_dp - g_value))
      result%v0 = 1.0_dp + c1 * b_value
      result%v1 = 1.0_dp - c2 * b_value
      result%v2 = 1.0_dp - c3 * b_value
      if (abs(result%v0) < 1.0e-10_dp) result%v0 = 0.0_dp
      result%rho = (div_zero(1.0_dp - result%a2 / qreal, result%v1) ** 2) * &
         result%v0 / (qreal * result%v2)
      result%ddf = 4.0_dp + (qreal + 2.0_dp) / (qreal * result%rho - 1.0_dp)
      if (abs(result%ddf - 2.0_dp) < 1.0e-2_dp) then
         result%f_scaling = 1.0_dp
      else
         result%f_scaling = result%ddf * (1.0_dp - result%a2 / qreal) / (result%ddf - 2.0_dp)
      end if

      allocate(beta_diff(p), lb2(q), solution_vector(q))
      beta_diff = beta - beta_h
      lb2 = matmul(l, beta_diff)
      a = matmul(matmul(l, phi_adjusted), transpose(l))
      call solve_system(a, lb2, solution_vector, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      wald = dot_product(lb2, solution_vector)
      a = matmul(matmul(l, phi), transpose(l))
      call solve_system(a, lb2, solution_vector, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      result%wald_unadjusted = dot_product(lb2, solution_vector)

      result%ndf = q
      result%f_stat_unscaled = wald / qreal
      result%p_value_unscaled = r_pf(result%f_stat_unscaled, qreal, result%ddf, lower_tail=.false.)
      result%f_stat = result%f_scaling * result%f_stat_unscaled
      result%p_value = r_pf(result%f_stat, qreal, result%ddf, lower_tail=.false.)
   end subroutine kr_adjust

   pure subroutine lb_ddf(l, v0, p_matrices, w, ddf, status)
      real(dp), intent(in) :: l(:, :) !! Restriction matrix with `p` columns; upstream uses its row count as numerator df.
      real(dp), intent(in) :: v0(:, :) !! Unadjusted fixed-effect covariance matrix, shape `(p,p)`.
      real(dp), intent(in) :: p_matrices(:, :, :) !! Kenward-Roger `P` matrices, shape `(p,p,m)`.
      real(dp), intent(in) :: w(:, :) !! Kenward-Roger covariance-parameter variance matrix, shape `(m,m)`.
      real(dp), intent(out) :: ddf !! Kenward-Roger denominator degrees of freedom for the restriction.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: solution(:, :)
      real(dp), allocatable :: theta(:, :)
      real(dp) :: a1
      real(dp) :: a2
      real(dp) :: b_value
      real(dp) :: c1
      real(dp) :: c2
      real(dp) :: c3
      real(dp) :: g_value
      real(dp) :: qreal
      real(dp) :: rho
      real(dp) :: v0s
      real(dp) :: v1
      real(dp) :: v2
      integer :: info
      integer :: m
      integer :: p
      integer :: q

      status = pbkr_success
      p = size(v0, 1)
      q = size(l, 1)
      m = size(p_matrices, 3)
      if (p <= 0 .or. q <= 0 .or. m <= 0) then
         status = pbkr_invalid_shape
         ddf = 0.0_dp
         return
      end if
      if (size(v0, 2) /= p .or. size(l, 2) /= p) then
         status = pbkr_invalid_shape
         ddf = 0.0_dp
         return
      end if
      if (size(p_matrices, 1) /= p .or. size(p_matrices, 2) /= p .or. &
          size(w, 1) /= m .or. size(w, 2) /= m) then
         status = pbkr_invalid_shape
         ddf = 0.0_dp
         return
      end if

      allocate(a(q, q), solution(q, p), theta(p, p))
      a = matmul(matmul(l, v0), transpose(l))
      call solve_system(a, l, solution, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         ddf = 0.0_dp
         return
      end if
      theta = matmul(transpose(l), solution)
      call accumulate_a1_a2(theta, v0, p_matrices, w, a1, a2)

      qreal = real(q, dp)
      b_value = (a1 + 6.0_dp * a2) / (2.0_dp * qreal)
      g_value = ((qreal + 1.0_dp) * a1 - (qreal + 4.0_dp) * a2) / ((qreal + 2.0_dp) * a2)
      c1 = g_value / (3.0_dp * qreal + 2.0_dp * (1.0_dp - g_value))
      c2 = (qreal - g_value) / (3.0_dp * qreal + 2.0_dp * (1.0_dp - g_value))
      c3 = (qreal + 2.0_dp - g_value) / (3.0_dp * qreal + 2.0_dp * (1.0_dp - g_value))
      v0s = 1.0_dp + c1 * b_value
      v1 = 1.0_dp - c2 * b_value
      v2 = 1.0_dp - c3 * b_value
      if (abs(v0s) < 1.0e-10_dp) v0s = 0.0_dp
      rho = (div_zero(1.0_dp - a2 / qreal, v1) ** 2) * v0s / (qreal * v2)
      ddf = 4.0_dp + (qreal + 2.0_dp) / (qreal * rho - 1.0_dp)
   end subroutine lb_ddf

   pure subroutine ddf_lb_scalar(vv0, lcoef, p_matrices, w, ddf, status)
      real(dp), intent(in) :: vv0(:, :) !! Unadjusted fixed-effect covariance matrix, shape `(p,p)`.
      real(dp), intent(in) :: lcoef(:) !! Single linear contrast coefficient vector, length `p`.
      real(dp), intent(in) :: p_matrices(:, :, :) !! Kenward-Roger `P` matrices, shape `(p,p,m)`.
      real(dp), intent(in) :: w(:, :) !! Kenward-Roger covariance-parameter variance matrix, shape `(m,m)`.
      real(dp), intent(out) :: ddf !! Scalar-contrast Kenward-Roger denominator degrees of freedom.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: theta(:, :)
      real(dp) :: a1
      real(dp) :: a2
      real(dp) :: b_value
      real(dp) :: c1
      real(dp) :: c2
      real(dp) :: c3
      real(dp) :: g_value
      real(dp) :: rho
      real(dp) :: v0s
      real(dp) :: v1
      real(dp) :: v2
      real(dp) :: vlb
      integer :: i
      integer :: j
      integer :: m
      integer :: p

      status = pbkr_success
      p = size(vv0, 1)
      m = size(p_matrices, 3)
      if (p <= 0 .or. m <= 0 .or. size(vv0, 2) /= p .or. size(lcoef) /= p) then
         status = pbkr_invalid_shape
         ddf = 0.0_dp
         return
      end if
      if (size(p_matrices, 1) /= p .or. size(p_matrices, 2) /= p .or. &
          size(w, 1) /= m .or. size(w, 2) /= m) then
         status = pbkr_invalid_shape
         ddf = 0.0_dp
         return
      end if
      vlb = dot_product(lcoef, matmul(vv0, lcoef))
      if (abs(vlb) <= tiny(1.0_dp)) then
         status = pbkr_invalid_argument
         ddf = 0.0_dp
         return
      end if
      allocate(theta(p, p))
      do i = 1, p
         do j = 1, p
            theta(i, j) = lcoef(i) * lcoef(j) / vlb
         end do
      end do
      call accumulate_a1_a2(theta, vv0, p_matrices, w, a1, a2)

      b_value = (a1 + 6.0_dp * a2) / 2.0_dp
      g_value = (2.0_dp * a1 - 5.0_dp * a2) / (3.0_dp * a2)
      c1 = g_value / (3.0_dp + 2.0_dp * (1.0_dp - g_value))
      c2 = (1.0_dp - g_value) / (3.0_dp + 2.0_dp * (1.0_dp - g_value))
      c3 = (3.0_dp - g_value) / (3.0_dp + 2.0_dp * (1.0_dp - g_value))
      v0s = 1.0_dp + c1 * b_value
      v1 = 1.0_dp - c2 * b_value
      v2 = 1.0_dp - c3 * b_value
      if (abs(v0s) < 1.0e-10_dp) v0s = 0.0_dp
      rho = div_zero(1.0_dp - a2, v1) ** 2 * v0s / v2
      ddf = 4.0_dp + 3.0_dp / (rho - 1.0_dp)
   end subroutine ddf_lb_scalar

   pure subroutine accumulate_a1_a2(theta, covariance, p_matrices, w, a1, a2)
      real(dp), intent(in) :: theta(:, :) !! Projection-like matrix used in the Kenward-Roger moment calculation.
      real(dp), intent(in) :: covariance(:, :) !! Fixed-effect covariance matrix multiplied around each `P` matrix.
      real(dp), intent(in) :: p_matrices(:, :, :) !! Kenward-Roger `P` matrices, shape `(p,p,m)`.
      real(dp), intent(in) :: w(:, :) !! Covariance-parameter variance matrix, shape `(m,m)`.
      real(dp), intent(out) :: a1 !! First Kenward-Roger moment correction accumulator.
      real(dp), intent(out) :: a2 !! Second Kenward-Roger moment correction accumulator.
      real(dp), allocatable :: theta_covariance(:, :)
      real(dp), allocatable :: ui(:, :)
      real(dp), allocatable :: uj(:, :)
      real(dp) :: multiplier
      integer :: i
      integer :: j
      integer :: m
      integer :: p

      p = size(covariance, 1)
      m = size(p_matrices, 3)
      allocate(theta_covariance(p, p), ui(p, p), uj(p, p))
      theta_covariance = matmul(theta, covariance)
      a1 = 0.0_dp
      a2 = 0.0_dp
      do i = 1, m
         do j = i, m
            if (i == j) then
               multiplier = 1.0_dp
            else
               multiplier = 2.0_dp
            end if
            ui = matmul(matmul(theta_covariance, p_matrices(:, :, i)), covariance)
            uj = matmul(matmul(theta_covariance, p_matrices(:, :, j)), covariance)
            a1 = a1 + multiplier * w(i, j) * trace_matrix(ui) * trace_matrix(uj)
            a2 = a2 + multiplier * w(i, j) * sum(ui * transpose(uj))
         end do
      end do
   end subroutine accumulate_a1_a2

   subroutine symmetric_pseudoinverse(a, inverse, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix for which a Moore-Penrose inverse is required.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Allocated symmetric pseudo-inverse with the same shape as `a`.
      integer, intent(out) :: info !! Zero on success or the underlying eigensolver status.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp) :: scale
      real(dp) :: threshold
      integer :: i
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         allocate(inverse(0, 0))
         info = -1
         return
      end if
      call symmetric_eigen(a, values, vectors, info)
      if (info /= 0) then
         allocate(inverse(0, 0))
         return
      end if
      allocate(inverse(n, n))
      inverse = 0.0_dp
      if (n == 0) return
      scale = maxval(abs(values))
      threshold = sqrt(epsilon(1.0_dp)) * scale
      do i = 1, n
         if (abs(values(i)) > threshold) then
            inverse = inverse + outer_product(vectors(:, i), vectors(:, i)) / values(i)
         end if
      end do
      inverse = symmetrize(inverse)
   end subroutine symmetric_pseudoinverse

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

end module pbkrtest_kr
