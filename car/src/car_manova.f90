module car_manova
   use car_kinds, only: dp
   use car_utils, only: f_survival_real, matrix_determinant, matrix_inverse
   implicit none
   private

   public :: manova_test_statistics
   public :: manova_type3_ssp

contains

   pure subroutine manova_type3_ssp(x, y, assign, hypothesis_ssp, error_ssp, term_df, error_df, status, weights)
      real(dp), intent(in) :: x(:, :) !! Full-rank between-subjects design matrix with observations in rows.
      real(dp), intent(in) :: y(:, :) !! Multivariate response matrix with observations in rows and responses in columns.
      integer, intent(in) :: assign(:) !! Term assignment for each coefficient; zero denotes an intercept and positive integers denote terms.
      real(dp), allocatable, intent(out) :: hypothesis_ssp(:, :, :) !! Type III hypothesis SSP matrix for each positive term number.
      real(dp), allocatable, intent(out) :: error_ssp(:, :) !! Residual sums-of-squares-and-products matrix of the full model.
      integer, allocatable, intent(out) :: term_df(:) !! Hypothesis degrees of freedom for each term.
      integer, intent(out) :: error_df !! Residual degrees of freedom n - p for the full-rank multivariate fit.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, weights, or singular matrices.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative case weights; omitted means unit weights.
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: c(:, :)
      real(dp), allocatable :: hypothesis(:, :)
      real(dp), allocatable :: m(:, :)
      real(dp), allocatable :: m_inverse(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: xtwx(:, :)
      real(dp), allocatable :: xtwx_inverse(:, :)
      real(dp), allocatable :: xtwy(:, :)
      real(dp) :: weight
      integer :: i
      integer :: inv_status
      integer :: j
      integer :: k
      integer :: n
      integer :: n_response
      integer :: n_terms
      integer :: p
      integer :: q

      n = size(x, 1)
      p = size(x, 2)
      n_response = size(y, 2)
      error_df = n - p
      status = 0
      if (n <= 0 .or. p <= 0 .or. n_response <= 0 .or. size(y, 1) /= n .or. size(assign) /= p .or. &
          any(assign < 0) .or. error_df <= 0) then
         status = 1
         allocate(hypothesis_ssp(0, 0, 0), error_ssp(0, 0), term_df(0))
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) then
            status = 2
            allocate(hypothesis_ssp(0, 0, 0), error_ssp(0, 0), term_df(0))
            return
         end if
      end if
      n_terms = maxval(assign)
      if (n_terms <= 0) then
         status = 3
         allocate(hypothesis_ssp(0, 0, 0), error_ssp(0, 0), term_df(0))
         return
      end if
      do i = 1, n_terms
         if (count(assign == i) == 0) then
            status = 4
            allocate(hypothesis_ssp(0, 0, 0), error_ssp(0, 0), term_df(0))
            return
         end if
      end do

      allocate(xtwx(p, p), xtwx_inverse(p, p), xtwy(p, n_response), beta(p, n_response), residuals(n, n_response))
      xtwx = 0.0_dp
      xtwy = 0.0_dp
      do i = 1, n
         if (present(weights)) then
            weight = weights(i)
         else
            weight = 1.0_dp
         end if
         do j = 1, p
            xtwy(j, :) = xtwy(j, :) + weight * x(i, j) * y(i, :)
            do k = 1, p
               xtwx(j, k) = xtwx(j, k) + weight * x(i, j) * x(i, k)
            end do
         end do
      end do
      call matrix_inverse(xtwx, xtwx_inverse, inv_status)
      if (inv_status /= 0) then
         status = 5
         allocate(hypothesis_ssp(0, 0, 0), error_ssp(0, 0), term_df(0))
         return
      end if
      beta = matmul(xtwx_inverse, xtwy)
      residuals = y - matmul(x, beta)

      allocate(error_ssp(n_response, n_response), hypothesis_ssp(n_response, n_response, n_terms), term_df(n_terms))
      error_ssp = 0.0_dp
      hypothesis_ssp = 0.0_dp
      term_df = 0
      do i = 1, n
         if (present(weights)) then
            weight = weights(i)
         else
            weight = 1.0_dp
         end if
         do j = 1, n_response
            do k = 1, n_response
               error_ssp(j, k) = error_ssp(j, k) + weight * residuals(i, j) * residuals(i, k)
            end do
         end do
      end do

      do i = 1, n_terms
         q = count(assign == i)
         term_df(i) = q
         allocate(hypothesis(q, p), c(q, n_response), m(q, q), m_inverse(q, q))
         hypothesis = 0.0_dp
         j = 0
         do k = 1, p
            if (assign(k) == i) then
               j = j + 1
               hypothesis(j, k) = 1.0_dp
            end if
         end do
         c = matmul(hypothesis, beta)
         m = matmul(hypothesis, matmul(xtwx_inverse, transpose(hypothesis)))
         call matrix_inverse(m, m_inverse, inv_status)
         if (inv_status /= 0) then
            status = 10 + i
            return
         end if
         hypothesis_ssp(:, :, i) = matmul(transpose(c), matmul(m_inverse, c))
         deallocate(hypothesis, c, m, m_inverse)
      end do
   end subroutine manova_type3_ssp

   pure subroutine manova_test_statistics(hypothesis_ssp, error_ssp, hypothesis_df, error_df, &
                                          pillai, pillai_f, pillai_num_df, pillai_den_df, pillai_p, &
                                          wilks, wilks_f, wilks_num_df, wilks_den_df, wilks_p, &
                                          hotelling_lawley, hotelling_f, hotelling_num_df, &
                                          hotelling_den_df, hotelling_p, status)
      real(dp), intent(in) :: hypothesis_ssp(:, :) !! Hypothesis sums-of-squares-and-products matrix for one MANOVA term.
      real(dp), intent(in) :: error_ssp(:, :) !! Residual sums-of-squares-and-products matrix from the full model.
      integer, intent(in) :: hypothesis_df !! Positive hypothesis degrees of freedom for the tested term.
      integer, intent(in) :: error_df !! Positive residual degrees of freedom of the multivariate linear model.
      real(dp), intent(out) :: pillai !! Pillai trace statistic.
      real(dp), intent(out) :: pillai_f !! Approximate F statistic corresponding to the Pillai trace.
      real(dp), intent(out) :: pillai_num_df !! Approximate numerator degrees of freedom for the Pillai F test.
      real(dp), intent(out) :: pillai_den_df !! Approximate denominator degrees of freedom for the Pillai F test.
      real(dp), intent(out) :: pillai_p !! Upper-tail probability for the Pillai approximate F statistic.
      real(dp), intent(out) :: wilks !! Wilks lambda statistic.
      real(dp), intent(out) :: wilks_f !! Approximate F statistic corresponding to Wilks lambda.
      real(dp), intent(out) :: wilks_num_df !! Approximate numerator degrees of freedom for the Wilks F test.
      real(dp), intent(out) :: wilks_den_df !! Approximate denominator degrees of freedom for the Wilks F test.
      real(dp), intent(out) :: wilks_p !! Upper-tail probability for the Wilks approximate F statistic.
      real(dp), intent(out) :: hotelling_lawley !! Hotelling-Lawley trace statistic.
      real(dp), intent(out) :: hotelling_f !! Approximate F statistic corresponding to the Hotelling-Lawley trace.
      real(dp), intent(out) :: hotelling_num_df !! Approximate numerator degrees of freedom for the Hotelling-Lawley F test.
      real(dp), intent(out) :: hotelling_den_df !! Approximate denominator degrees of freedom for the Hotelling-Lawley F test.
      real(dp), intent(out) :: hotelling_p !! Upper-tail probability for the Hotelling-Lawley approximate F statistic.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or singular SSP matrices.
      real(dp), allocatable :: e_inverse(:, :)
      real(dp), allocatable :: total_inverse(:, :)
      real(dp), allocatable :: total_ssp(:, :)
      real(dp) :: det_e
      real(dp) :: det_total
      real(dp) :: m
      real(dp) :: n_value
      real(dp) :: s_value
      real(dp) :: tmp1
      real(dp) :: tmp2
      real(dp) :: tmp3
      real(dp) :: tmp3_arg
      integer :: det_status
      integer :: inv_status
      integer :: p

      pillai = 0.0_dp
      pillai_f = 0.0_dp
      pillai_num_df = 0.0_dp
      pillai_den_df = 0.0_dp
      pillai_p = 1.0_dp
      wilks = 0.0_dp
      wilks_f = 0.0_dp
      wilks_num_df = 0.0_dp
      wilks_den_df = 0.0_dp
      wilks_p = 1.0_dp
      hotelling_lawley = 0.0_dp
      hotelling_f = 0.0_dp
      hotelling_num_df = 0.0_dp
      hotelling_den_df = 0.0_dp
      hotelling_p = 1.0_dp
      status = 0
      p = size(error_ssp, 1)
      if (p <= 0 .or. size(error_ssp, 2) /= p .or. size(hypothesis_ssp, 1) /= p .or. &
          size(hypothesis_ssp, 2) /= p .or. hypothesis_df <= 0 .or. error_df <= 0) then
         status = 1
         return
      end if

      allocate(e_inverse(p, p), total_ssp(p, p), total_inverse(p, p))
      total_ssp = error_ssp + hypothesis_ssp
      call matrix_inverse(error_ssp, e_inverse, inv_status)
      if (inv_status /= 0) then
         status = 2
         return
      end if
      call matrix_inverse(total_ssp, total_inverse, inv_status)
      if (inv_status /= 0) then
         status = 3
         return
      end if
      call matrix_determinant(error_ssp, det_e, det_status)
      if (det_status /= 0) then
         status = 4
         return
      end if
      call matrix_determinant(total_ssp, det_total, det_status)
      if (det_status /= 0 .or. det_total <= 0.0_dp .or. det_e < 0.0_dp) then
         status = 4
         return
      end if

      pillai = matrix_trace(matmul(hypothesis_ssp, total_inverse))
      wilks = det_e / det_total
      hotelling_lawley = matrix_trace(matmul(e_inverse, hypothesis_ssp))

      s_value = real(min(p, hypothesis_df), dp)
      n_value = 0.5_dp * real(error_df - p - 1, dp)
      m = 0.5_dp * (real(abs(p - hypothesis_df), dp) - 1.0_dp)
      tmp1 = 2.0_dp * m + s_value + 1.0_dp
      tmp2 = 2.0_dp * n_value + s_value + 1.0_dp
      pillai_num_df = s_value * tmp1
      pillai_den_df = s_value * tmp2
      if (s_value - pillai > 0.0_dp .and. pillai_num_df > 0.0_dp .and. pillai_den_df > 0.0_dp) then
         pillai_f = (tmp2 / tmp1 * pillai) / (s_value - pillai)
         pillai_p = f_survival_real(pillai_f, pillai_num_df, pillai_den_df)
      else
         status = 5
         return
      end if

      tmp1 = real(error_df, dp) - 0.5_dp * real(p - hypothesis_df + 1, dp)
      tmp2 = (real(p * hypothesis_df, dp) - 2.0_dp) / 4.0_dp
      tmp3_arg = real(p * p + hypothesis_df * hypothesis_df - 5, dp)
      if (tmp3_arg > 0.0_dp) then
         tmp3 = sqrt((real((p * hypothesis_df)**2 - 4, dp)) / tmp3_arg)
      else
         tmp3 = 1.0_dp
      end if
      wilks_num_df = real(p * hypothesis_df, dp)
      wilks_den_df = tmp1 * tmp3 - 2.0_dp * tmp2
      if (wilks > 0.0_dp .and. wilks_num_df > 0.0_dp .and. wilks_den_df > 0.0_dp) then
         wilks_f = ((wilks**(-1.0_dp / tmp3) - 1.0_dp) * wilks_den_df) / real(p * hypothesis_df, dp)
         wilks_p = f_survival_real(wilks_f, wilks_num_df, wilks_den_df)
      else
         status = 6
         return
      end if

      tmp1 = 2.0_dp * m + s_value + 1.0_dp
      tmp2 = 2.0_dp * (s_value * n_value + 1.0_dp)
      hotelling_num_df = s_value * tmp1
      hotelling_den_df = tmp2
      if (hotelling_num_df > 0.0_dp .and. hotelling_den_df > 0.0_dp) then
         hotelling_f = (tmp2 * hotelling_lawley) / (s_value * s_value * tmp1)
         hotelling_p = f_survival_real(hotelling_f, hotelling_num_df, hotelling_den_df)
      else
         status = 7
      end if
   end subroutine manova_test_statistics

   pure function matrix_trace(a) result(trace_value)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose diagonal sum is returned.
      real(dp) :: trace_value
      integer :: i

      trace_value = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         trace_value = trace_value + a(i, i)
      end do
   end function matrix_trace

end module car_manova
