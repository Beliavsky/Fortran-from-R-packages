program test_car
   use car
   implicit none

   integer :: failures

   failures = 0
   call test_transforms(failures)
   call test_contrasts_and_geometry(failures)
   call test_kernels(failures)
   call test_matrix_statistics(failures)
   call test_diagnostics(failures)
   call test_data_helpers(failures)
   call test_power_fits(failures)
   call test_box_tidwell(failures)
   call test_po_test(failures)
   call test_lm_helpers(failures)
   call test_anova_lm(failures)
   call test_manova(failures)
   call test_remaining_helpers(failures)

   if (failures /= 0) then
      write(*, '(a,i0)') 'test_car: FAIL, failures = ', failures
      error stop 1
   end if
   write(*, '(a)') 'test_car: PASS'

contains

   subroutine test_transforms(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: u(3), z(3), back(3), p(3), logits(3), bcv(3)
      integer :: status, i

      u = [1.0_dp, 2.0_dp, 4.0_dp]
      call bc_power(u, 0.0_dp, z, status)
      call require(status == 0, 'bcPower status', failures)
      call require_all_close(z, log(u), 2.0e-12_dp, 'bcPower log limit', failures)

      call bc_power(u, 1.0_dp, z, status)
      call require_all_close(z, u - 1.0_dp, 2.0e-12_dp, 'bcPower lambda one', failures)

      u = [0.0_dp, 1.0_dp, -1.0_dp]
      call yj_power(u, 1.0_dp, z, status)
      call require(status == 0, 'yjPower status', failures)
      call require_all_close(z, u, 2.0e-12_dp, 'yjPower identity', failures)

      u = [0.25_dp, 1.0_dp, 3.0_dp]
      call bcn_power(u, 0.4_dp, 0.7_dp, z, status)
      call require(status == 0, 'bcnPower status', failures)
      do i = 1, 3
         back(i) = bcn_power_inverse(z(i), 0.4_dp, 0.7_dp)
      end do
      call require_all_close(back, u, 2.0e-11_dp, 'bcnPower inverse', failures)

      p = [0.25_dp, 0.5_dp, 0.75_dp]
      call logit_transform(p, logits, status, percents=.false., adjust=0.0_dp)
      call require(status == 0, 'logit status', failures)
      call require_all_close(logits, log(p / (1.0_dp - p)), 2.0e-12_dp, 'logit values', failures)

      u = [1.0_dp, 2.0_dp, 4.0_dp]
      call box_cox_variable(u, bcv, status)
      call require(status == 0, 'boxCoxVariable status', failures)
      call require_close(sum(bcv), sum(u * (log(u / exp(sum(log(u)) / 3.0_dp)) - 1.0_dp)), &
                         2.0e-12_dp, 'boxCoxVariable values', failures)
   end subroutine test_transforms

   subroutine test_contrasts_and_geometry(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp), allocatable :: contrast(:, :), xy(:, :)
      real(dp) :: center(2), shape(2, 2)
      integer :: status

      call treatment_contrast(3, contrast, status)
      call require(status == 0, 'treatment contrast status', failures)
      call require_all_close(reshape(contrast, [6]), [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
                             1.0e-14_dp, 'treatment contrast values', failures)

      call sum_contrast(3, contrast, status)
      call require_all_close(reshape(contrast, [6]), [1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, -1.0_dp], &
                             1.0e-14_dp, 'sum contrast values', failures)

      call helmert_contrast(3, contrast, status)
      call require_all_close(reshape(contrast, [6]), [-1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, 2.0_dp], &
                             1.0e-14_dp, 'Helmert contrast values', failures)

      center = [2.0_dp, -1.0_dp]
      shape = reshape([4.0_dp, 0.0_dp, 0.0_dp, 9.0_dp], [2, 2])
      call ellipse_coordinates(center, shape, 1.0_dp, 4, xy, status)
      call require(status == 0, 'ellipse status', failures)
      call require_close(xy(1, 1), 4.0_dp, 1.0e-12_dp, 'ellipse first x', failures)
      call require_close(xy(1, 2), -1.0_dp, 1.0e-12_dp, 'ellipse first y', failures)
      call require_all_close(xy(1, :), xy(5, :), 1.0e-12_dp, 'ellipse closure', failures)
   end subroutine test_contrasts_and_geometry

   subroutine test_kernels(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: x(5), bw, used
      real(dp), allocatable :: grid(:), density(:), initial(:)
      integer :: status

      call require_close(epanechnikov_kernel(0.0_dp), 0.75_dp, 1.0e-14_dp, 'Epanechnikov center', failures)
      call require_close(epanechnikov_kernel(2.0_dp), 0.0_dp, 1.0e-14_dp, 'Epanechnikov support', failures)
      call require_close(biweight_kernel(0.0_dp), 15.0_dp / 16.0_dp, 1.0e-14_dp, 'biweight center', failures)

      x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      call bw_nrd0(x, bw, status)
      call require(status == 0 .and. bw > 0.0_dp, 'bw.nrd0', failures)
      call adaptive_kernel(x, grid, density, initial, used, status, n_points=31, from=0.0_dp, to=6.0_dp)
      call require(status == 0, 'adaptiveKernel status', failures)
      call require(size(grid) == 31 .and. size(density) == 31, 'adaptiveKernel size', failures)
      call require_close(grid(1), 0.0_dp, 1.0e-14_dp, 'adaptiveKernel lower grid', failures)
      call require_close(grid(31), 6.0_dp, 1.0e-14_dp, 'adaptiveKernel upper grid', failures)
      call require(all(density >= 0.0_dp), 'adaptiveKernel nonnegative', failures)
   end subroutine test_kernels

   subroutine test_matrix_statistics(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: x(3, 2), cp_expected(2, 2), vcov(2, 2), leverage(3), residuals(3)
      real(dp) :: hc(2, 2), beta(2), l(1, 2), rhs(1), statistic, p_value, transformed, se
      real(dp), allocatable :: cp(:, :), gvif(:), adjusted(:)
      integer, allocatable :: gvif_df(:)
      integer :: assign(2), status

      x = reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp], [3, 2])
      cp = wcrossprod(x)
      cp_expected = reshape([3.0_dp, 3.0_dp, 3.0_dp, 5.0_dp], [2, 2])
      call require_all_close(reshape(cp, [4]), reshape(cp_expected, [4]), 1.0e-14_dp, 'wcrossprod', failures)

      vcov = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
      assign = [1, 2]
      call gvif_from_covariance(vcov, assign, gvif, gvif_df, adjusted, status)
      call require(status == 0, 'GVIF status', failures)
      call require_all_close(gvif, [4.0_dp / 3.0_dp, 4.0_dp / 3.0_dp], 2.0e-12_dp, 'GVIF values', failures)

      vcov = reshape([5.0_dp / 6.0_dp, -0.5_dp, -0.5_dp, 0.5_dp], [2, 2])
      leverage = [5.0_dp / 6.0_dp, 1.0_dp / 3.0_dp, 5.0_dp / 6.0_dp]
      residuals = [1.0_dp, -1.0_dp, 0.5_dp]
      call hccm_matrix(x, residuals, vcov, leverage, hc, status, hc_type='hc0')
      call require(status == 0, 'hccm status', failures)
      call require_all_close(reshape(hc, [4]), [0.8125_dp, -0.4375_dp, -0.4375_dp, 0.3125_dp], &
                             2.0e-12_dp, 'hccm HC0', failures)

      beta = [1.0_dp, 2.0_dp]
      vcov = reshape([0.25_dp, 0.0_dp, 0.0_dp, 0.25_dp], [2, 2])
      l = reshape([0.0_dp, 1.0_dp], [1, 2])
      rhs = 0.0_dp
      call linear_hypothesis_test(beta, vcov, l, rhs, statistic, p_value, status)
      call require(status == 0, 'linear hypothesis status', failures)
      call require_close(statistic, 16.0_dp, 1.0e-12_dp, 'linear hypothesis chi-square', failures)
      call require(p_value < 1.0e-4_dp, 'linear hypothesis p-value', failures)

      call delta_method_gradient(3.0_dp, [1.0_dp, 2.0_dp], vcov, transformed, se, status)
      call require(status == 0, 'deltaMethod status', failures)
      call require_close(transformed, 3.0_dp, 1.0e-14_dp, 'deltaMethod estimate', failures)
      call require_close(se, sqrt(1.25_dp), 1.0e-12_dp, 'deltaMethod standard error', failures)

      call require_close(sigma_hat([1.0_dp, -1.0_dp, 2.0_dp], 2), sqrt(3.0_dp), 1.0e-12_dp, &
                         'sigmaHat', failures)
   end subroutine test_matrix_statistics

   subroutine test_diagnostics(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp), allocatable :: dw(:), ac(:), p(:), bp(:)
      real(dp) :: stat, pval
      integer :: status, df1, df2

      call durbin_watson([1.0_dp, 2.0_dp, 1.0_dp], 1, dw, ac, status)
      call require(status == 0, 'Durbin-Watson status', failures)
      call require_close(dw(1), 1.0_dp / 3.0_dp, 1.0e-12_dp, 'Durbin-Watson value', failures)
      call require_close(ac(1), 2.0_dp / 3.0_dp, 1.0e-12_dp, 'Durbin-Watson autocorrelation', failures)

      call levene_test([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], [1, 1, 1, 2, 2, 2], &
                       stat, pval, df1, df2, status)
      call require(status == 0 .and. df1 == 1 .and. df2 == 4, 'Levene status and df', failures)
      call require_close(stat, 0.8_dp, 1.0e-12_dp, 'Levene statistic', failures)

      call ncv_test([1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp], [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], 2, &
                    stat, pval, status)
      call require(status == 0, 'ncvTest status', failures)
      call require_close(stat, 0.576_dp, 2.0e-12_dp, 'ncvTest statistic', failures)

      call outlier_bonferroni([0.0_dp, 2.0_dp], 10, p, bp, status)
      call require(status == 0, 'outlierTest status', failures)
      call require_close(p(1), 1.0_dp, 2.0e-12_dp, 'outlierTest zero residual p', failures)
      call require(bp(2) >= p(2) .and. bp(2) <= 1.0_dp, 'outlierTest Bonferroni', failures)
   end subroutine test_diagnostics

   subroutine test_data_helpers(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      character(len=8) :: object_names(4)
      character(len=8) :: query_names(3)
      character(len=5) :: string_data(4, 2)
      character(len=:), allocatable :: levels(:, :)
      logical, allocatable :: converted(:)
      logical, allocatable :: matched(:)
      real(dp) :: coefficient_matrix(2, 2)
      real(dp) :: covariance(2, 2, 2)
      real(dp), allocatable :: grouped(:)
      real(dp), allocatable :: p_values(:, :)
      real(dp), allocatable :: result(:)
      real(dp), allocatable :: std_errors(:, :)
      real(dp), allocatable :: z_values(:, :)
      integer, allocatable :: codes(:, :)
      integer, allocatable :: counts(:)
      integer, allocatable :: indices(:)
      integer, allocatable :: n_levels(:)
      integer, allocatable :: sampled(:)
      integer :: status

      call recode_numeric([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1.0_dp, 3.0_dp], [2.0_dp, 4.0_dp], &
                          [10.0_dp, 20.0_dp], result, status)
      call require(status == 0, 'recode status', failures)
      call require_all_close(result, [10.0_dp, 10.0_dp, 20.0_dp, 20.0_dp], 1.0e-14_dp, 'recode values', failures)

      call tapply_numeric([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], [1, 1, 2, 2], 'mean', grouped, counts, status)
      call require(status == 0, 'Tapply status', failures)
      call require_all_close(grouped, [2.0_dp, 3.0_dp], 1.0e-14_dp, 'Tapply means', failures)
      call require(all(counts == [2, 2]), 'Tapply counts', failures)

      object_names = ['alpha   ', 'beta    ', 'gamma   ', 'delta   ']
      query_names = ['delta   ', 'missing ', 'beta    ']
      call which_names(object_names, query_names, indices, matched, status)
      call require(status == 0, 'whichNames status', failures)
      call require(all(matched .eqv. [.true., .false., .true.]), 'whichNames matched flags', failures)
      call require(size(indices) == 2 .and. all(indices == [4, 2]), 'whichNames positions', failures)

      call some_indices(6, 3, sampled, status, random_values=[0.9_dp, 0.0_dp, 0.4_dp])
      call require(status == 0, 'some status', failures)
      call require(size(sampled) == 3 .and. all(sampled == [2, 4, 6]), 'some sampled positions', failures)

      string_data(:, 1) = ['b    ', 'a    ', 'b    ', 'c    ']
      string_data(:, 2) = ['u    ', 'v    ', 'w    ', 'x    ']
      call strings2factors_codes(string_data, codes, levels, n_levels, converted, status)
      call require(status == 0, 'strings2factors status', failures)
      call require(all(converted .eqv. [.true., .false.]), 'strings2factors converted columns', failures)
      call require(all(n_levels == [3, 0]), 'strings2factors level counts', failures)
      call require(all(codes(:, 1) == [2, 1, 2, 3]), 'strings2factors codes', failures)
      call require(trim(levels(1, 1)) == 'a' .and. trim(levels(2, 1)) == 'b' .and. &
                   trim(levels(3, 1)) == 'c', 'strings2factors levels', failures)

      coefficient_matrix(:, 1) = [1.0_dp, 2.0_dp]
      coefficient_matrix(:, 2) = [1.5_dp, 3.0_dp]
      covariance = 0.0_dp
      covariance(1, 1, 1) = 0.25_dp
      covariance(2, 2, 1) = 1.0_dp
      covariance(1, 1, 2) = 0.09_dp
      covariance(2, 2, 2) = 0.25_dp
      call compare_coefs_numeric(coefficient_matrix, covariance, std_errors, z_values, p_values, status)
      call require(status == 0, 'compareCoefs status', failures)
      call require_all_close(std_errors(:, 1), [0.5_dp, 1.0_dp], 1.0e-14_dp, 'compareCoefs model one SE', failures)
      call require_all_close(std_errors(:, 2), [0.3_dp, 0.5_dp], 1.0e-14_dp, 'compareCoefs model two SE', failures)
      call require_all_close(z_values(:, 1), [2.0_dp, 2.0_dp], 1.0e-14_dp, 'compareCoefs model one z', failures)
      call require_all_close(z_values(:, 2), [5.0_dp, 6.0_dp], 1.0e-14_dp, 'compareCoefs model two z', failures)
      call require(all(p_values >= 0.0_dp .and. p_values <= 1.0_dp), 'compareCoefs p values', failures)
   end subroutine test_data_helpers

   subroutine test_power_fits(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: design(8, 2), y(8), predictor(8), lambda, objective, lr, p_value
      real(dp), allocatable :: profile(:)
      integer :: i, status

      do i = 1, 8
         design(i, 1) = 1.0_dp
         design(i, 2) = real(i - 1, dp)
      end do
      y = [1.10_dp, 1.35_dp, 1.82_dp, 2.30_dp, 3.10_dp, 4.05_dp, 5.60_dp, 7.10_dp]
      call power_transform_fit(design, y, 'bcPower', lambda, objective, status)
      call require(status == 0 .and. lambda >= -3.0_dp .and. lambda <= 3.0_dp, 'powerTransform fit', failures)
      call power_transform_test(design, y, 'bcPower', lambda, 1.0_dp, lr, p_value, status)
      call require(status == 0 .and. lr >= 0.0_dp .and. p_value >= 0.0_dp .and. p_value <= 1.0_dp, &
                   'testTransform', failures)
      call box_cox_profile(design, y, [-1.0_dp, 0.0_dp, 1.0_dp], profile, status)
      call require(status == 0 .and. size(profile) == 3, 'boxCox profile', failures)

      predictor = [1.0_dp, 1.5_dp, 2.2_dp, 3.0_dp, 4.5_dp, 6.5_dp, 9.0_dp, 12.0_dp]
      y = 1.0_dp + 2.0_dp * log(predictor) + [0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp, -0.02_dp, 0.03_dp, -0.01_dp, -0.04_dp]
      call inv_transform_estimate(predictor, y, 'bcPower', lambda, status)
      call require(status == 0 .and. abs(lambda) < 1.0_dp, 'invTranEstimate', failures)
   end subroutine test_power_fits


   subroutine test_box_tidwell(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: x1(12, 1)
      real(dp) :: y(12)
      real(dp) :: omnibus_f
      real(dp) :: omnibus_p
      real(dp), allocatable :: powers(:)
      real(dp), allocatable :: p_values(:)
      real(dp), allocatable :: score_t(:)
      integer :: i
      integer :: iterations
      integer :: status
      logical :: converged

      do i = 1, 12
         x1(i, 1) = real(i, dp)
         y(i) = 1.0_dp + 2.5_dp * sqrt(x1(i, 1)) + 0.01_dp * real(mod(i, 3) - 1, dp)
      end do
      call box_tidwell(y, x1, powers, score_t, p_values, iterations, status, omnibus_f=omnibus_f, &
                       omnibus_p=omnibus_p, converged=converged)
      call require(status == 0 .and. converged, 'boxTidwell status', failures)
      call require(size(powers) == 1 .and. iterations > 0, 'boxTidwell shape and iterations', failures)
      call require_close(powers(1), 0.4997835203448151_dp, 5.0e-10_dp, 'boxTidwell power', failures)
      call require(p_values(1) >= 0.0_dp .and. p_values(1) <= 1.0_dp, 'boxTidwell score p-value', failures)
      call require(omnibus_f >= 0.0_dp .and. omnibus_p >= 0.0_dp .and. omnibus_p <= 1.0_dp, &
                   'boxTidwell omnibus test', failures)
   end subroutine test_box_tidwell

   subroutine test_po_test(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: x(15, 2)
      real(dp) :: overall_chisq
      real(dp) :: overall_p
      real(dp), allocatable :: beta_threshold(:, :)
      real(dp), allocatable :: predictor_chisq(:)
      real(dp), allocatable :: predictor_p(:)
      integer :: overall_df
      integer :: predictor_df
      integer :: status
      integer :: y_level(15)

      x(:, 1) = 1.0_dp
      x(:, 2) = [-2.0_dp, -2.0_dp, -2.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                  1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
      y_level = [1, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 3]
      call po_test(x, y_level, overall_chisq, overall_df, overall_p, predictor_chisq, predictor_p, &
                   predictor_df, status, beta_threshold=beta_threshold)
      call require(status == 0, 'poTest status', failures)
      call require(overall_df == 1 .and. predictor_df == 1, 'poTest degrees of freedom', failures)
      call require(size(beta_threshold, 1) == 1 .and. size(beta_threshold, 2) == 2, 'poTest beta shape', failures)
      call require_close(overall_chisq, 0.0_dp, 1.0e-12_dp, 'poTest symmetric slopes', failures)
      call require(overall_p >= 0.0_dp .and. overall_p <= 1.0_dp, 'poTest overall p-value', failures)
      call require(size(predictor_chisq) == 1 .and. size(predictor_p) == 1, 'poTest predictor shape', failures)
   end subroutine test_po_test

   subroutine test_lm_helpers(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: x(5, 2)
      real(dp) :: x_new(2, 2)
      real(dp) :: y(5)
      real(dp) :: sigma
      real(dp) :: r_squared
      real(dp) :: adjusted_r_squared
      real(dp) :: f_statistic
      real(dp) :: f_p_value
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: t_values(:)
      real(dp), allocatable :: p_values(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: fit(:)
      real(dp), allocatable :: pred_se(:)
      real(dp), allocatable :: pred_lower(:)
      real(dp), allocatable :: pred_upper(:)
      real(dp), allocatable :: boot(:, :)
      integer :: indices(5, 2)
      integer :: residual_df
      integer :: status

      x(:, 1) = 1.0_dp
      x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      y = [1.1_dp, 2.9_dp, 5.2_dp, 6.8_dp, 9.1_dp]
      call lm_summary(x, y, beta, se, t_values, p_values, sigma, r_squared, adjusted_r_squared, &
                      f_statistic, f_p_value, residual_df, status, covariance=covariance)
      call require(status == 0 .and. residual_df == 3, 'S lm numeric status', failures)
      call require_all_close(beta, [1.04_dp, 1.99_dp], 2.0e-12_dp, 'S lm coefficients', failures)
      call require_close(sigma, 0.18885620632287084_dp, 2.0e-12_dp, 'S lm sigma', failures)
      call require_close(r_squared, 0.9973053289009771_dp, 2.0e-12_dp, 'S lm r squared', failures)
      call require_close(adjusted_r_squared, 0.9964071052013028_dp, 2.0e-12_dp, &
                         'S lm adjusted r squared', failures)
      call require_close(f_statistic, 1110.3084112149502_dp, 2.0e-9_dp, 'S lm F statistic', failures)
      call require(f_p_value < 1.0e-4_dp .and. all(p_values < 0.01_dp), 'S lm p values', failures)

      call confint_lm(beta, covariance, residual_df, 0.95_dp, lower, upper, status)
      call require(status == 0, 'Confint lm status', failures)
      call require_all_close(0.5_dp * (lower + upper), beta, 2.0e-12_dp, 'Confint lm centers', failures)
      call require(all(lower < beta) .and. all(upper > beta), 'Confint lm bounds', failures)

      x_new = reshape([1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp], [2, 2])
      call predict_lm(beta, covariance, x_new, fit, pred_se, pred_lower, pred_upper, status, &
                      level=0.95_dp, residual_df=residual_df)
      call require(status == 0, 'Predict lm status', failures)
      call require_all_close(fit, [5.02_dp, 9.0_dp], 2.0e-12_dp, 'Predict lm fits', failures)
      call require(all(pred_se > 0.0_dp), 'Predict lm standard errors', failures)
      call require(all(pred_lower < fit) .and. all(pred_upper > fit), 'Predict lm confidence bounds', failures)

      indices(:, 1) = [1, 2, 3, 4, 5]
      indices(:, 2) = [1, 1, 3, 4, 5]
      call lm_bootstrap_case(x, y, indices, boot, status)
      call require(status == 0, 'Boot lm case status', failures)
      call require(size(boot, 1) == 2 .and. size(boot, 2) == 2, 'Boot lm case shape', failures)
      call require_all_close(boot(1, :), beta, 2.0e-12_dp, 'Boot lm identity resample', failures)
   end subroutine test_lm_helpers

   subroutine test_anova_lm(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: x(5, 2)
      real(dp) :: y(5)
      real(dp) :: residual_ss
      real(dp), allocatable :: term_ss(:)
      real(dp), allocatable :: f_values(:)
      real(dp), allocatable :: p_values(:)
      integer, allocatable :: term_df(:)
      integer :: assign(2)
      integer :: residual_df
      integer :: status

      x(:, 1) = 1.0_dp
      x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      y = [1.1_dp, 2.9_dp, 5.2_dp, 6.8_dp, 9.1_dp]
      assign = [0, 1]
      call anova_lm_type3(x, y, assign, term_ss, term_df, f_values, p_values, residual_ss, residual_df, status)
      call require(status == 0 .and. residual_df == 3, 'Anova lm Type III status', failures)
      call require(size(term_ss) == 1 .and. term_df(1) == 1, 'Anova lm Type III shape', failures)
      call require_close(term_ss(1), 39.601_dp, 2.0e-10_dp, 'Anova lm Type III sum squares', failures)
      call require_close(residual_ss, 0.107_dp, 2.0e-12_dp, 'Anova lm residual sum squares', failures)
      call require_close(f_values(1), 1110.3084112149502_dp, 2.0e-9_dp, 'Anova lm Type III F', failures)
      call require(p_values(1) < 1.0e-4_dp, 'Anova lm Type III p value', failures)
   end subroutine test_anova_lm

   subroutine test_manova(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      real(dp) :: x(5, 2)
      real(dp) :: y(5, 2)
      real(dp), allocatable :: h(:, :, :)
      real(dp), allocatable :: e(:, :)
      integer, allocatable :: term_df(:)
      integer :: assign(2)
      integer :: error_df
      integer :: status
      real(dp) :: pillai
      real(dp) :: pillai_f
      real(dp) :: pillai_num_df
      real(dp) :: pillai_den_df
      real(dp) :: pillai_p
      real(dp) :: wilks
      real(dp) :: wilks_f
      real(dp) :: wilks_num_df
      real(dp) :: wilks_den_df
      real(dp) :: wilks_p
      real(dp) :: hotelling
      real(dp) :: hotelling_f
      real(dp) :: hotelling_num_df
      real(dp) :: hotelling_den_df
      real(dp) :: hotelling_p

      x(:, 1) = 1.0_dp
      x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      y(:, 1) = [1.1_dp, 2.9_dp, 5.2_dp, 6.8_dp, 9.1_dp]
      y(:, 2) = [2.0_dp, 4.1_dp, 6.0_dp, 8.2_dp, 10.1_dp]
      assign = [0, 1]
      call manova_type3_ssp(x, y, assign, h, e, term_df, error_df, status)
      call require(status == 0 .and. error_df == 3, 'Manova Type III status', failures)
      call require(size(h, 3) == 1 .and. term_df(1) == 1, 'Manova Type III shape', failures)
      call require_all_close(reshape(e, [4]), [0.107_dp, -0.045_dp, -0.045_dp, 0.019_dp], &
                             2.0e-12_dp, 'Manova error SSP', failures)
      call require_all_close(reshape(h(:, :, 1), [4]), [39.601_dp, 40.397_dp, 40.397_dp, 41.209_dp], &
                             2.0e-10_dp, 'Manova hypothesis SSP', failures)
      call manova_test_statistics(h(:, :, 1), e, term_df(1), error_df, pillai, pillai_f, &
                                  pillai_num_df, pillai_den_df, pillai_p, wilks, wilks_f, &
                                  wilks_num_df, wilks_den_df, wilks_p, hotelling, hotelling_f, &
                                  hotelling_num_df, hotelling_den_df, hotelling_p, status)
      call require(status == 0, 'Manova test statistics status', failures)
      call require_close(pillai, 0.9999990906528489_dp, 2.0e-10_dp, 'Manova Pillai', failures)
      call require_close(wilks, 9.093471796596996e-7_dp, 2.0e-13_dp, 'Manova Wilks', failures)
      call require_close(hotelling, 1099688.9999999767_dp, 2.0e-6_dp, 'Manova Hotelling-Lawley', failures)
      call require(pillai_p < 1.0e-5_dp .and. wilks_p < 1.0e-5_dp .and. hotelling_p < 1.0e-5_dp, &
                   'Manova multivariate p values', failures)
   end subroutine test_manova

   subroutine test_remaining_helpers(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count updated by this test group.
      character(len=12) :: coefficient_names(4)
      character(len=40) :: hypotheses(2)
      real(dp) :: x(5, 2)
      real(dp) :: y(5)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: hypothesis_matrix(:, :)
      real(dp), allocatable :: rhs(:)
      integer, allocatable :: matches(:)
      integer :: status

      x(:, 1) = 1.0_dp
      x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      y = [1.1_dp, 2.9_dp, 5.2_dp, 6.8_dp, 9.1_dp]
      call lm_refit_rows(x, y, [1, 2, 3, 4, 5], beta, status)
      call require(status == 0, 'nextBoot lm refit status', failures)
      call require_all_close(beta, [1.04_dp, 1.99_dp], 2.0e-12_dp, 'nextBoot lm refit coefficients', failures)

      coefficient_names = ['(Intercept) ', 'age         ', 'age:sex     ', 'income      ']
      call match_coefficient_names(coefficient_names, 'age', matches, status)
      call require(status == 0 .and. size(matches) == 2, 'matchCoefs status and size', failures)
      call require(all(matches == [2, 3]), 'matchCoefs literal matches', failures)

      coefficient_names(1:3) = ['(Intercept) ', 'x           ', 'z           ']
      hypotheses(1) = '2*x - z = 3'
      hypotheses(2) = 'x + z'
      call make_hypothesis(coefficient_names(1:3), hypotheses, hypothesis_matrix, rhs, status)
      call require(status == 0, 'makeHypothesis parser status', failures)
      call require_all_close(hypothesis_matrix(1, :), [0.0_dp, 2.0_dp, -1.0_dp], 1.0e-14_dp, &
                             'makeHypothesis first row', failures)
      call require_close(rhs(1), 3.0_dp, 1.0e-14_dp, 'makeHypothesis first rhs', failures)
      call require_all_close(hypothesis_matrix(2, :), [0.0_dp, 1.0_dp, 1.0_dp], 1.0e-14_dp, &
                             'makeHypothesis second row', failures)
      call require_close(rhs(2), 0.0_dp, 1.0e-14_dp, 'makeHypothesis second rhs', failures)
   end subroutine test_remaining_helpers

   subroutine require(condition, label, failures)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: label !! Short description printed when the assertion fails.
      integer, intent(inout) :: failures !! Running failure count incremented for a failed assertion.

      if (.not. condition) then
         failures = failures + 1
         write(*, '(a,a)') 'FAIL: ', trim(label)
      end if
   end subroutine require

   subroutine require_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tolerance !! Absolute comparison tolerance.
      character(len=*), intent(in) :: label !! Short description printed when comparison fails.
      integer, intent(inout) :: failures !! Running failure count incremented for a failed comparison.

      call require(abs(actual - expected) <= tolerance, label, failures)
   end subroutine require_close

   subroutine require_all_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual(:) !! Computed vector under test.
      real(dp), intent(in) :: expected(:) !! Reference vector with the same size as actual.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute elementwise difference.
      character(len=*), intent(in) :: label !! Short description printed when comparison fails.
      integer, intent(inout) :: failures !! Running failure count incremented for a failed comparison.

      call require(size(actual) == size(expected), trim(label)//' size', failures)
      if (size(actual) == size(expected)) then
         call require(maxval(abs(actual - expected)) <= tolerance, label, failures)
      end if
   end subroutine require_all_close

end program test_car
