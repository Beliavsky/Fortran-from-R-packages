program test_pbkrtest
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use pbkrtest
   implicit none

   call test_sigma_g()
   call test_column_spaces()
   call test_hypothesis_helpers()
   call test_kr_kernels()
   call test_satterthwaite_kernels()
   call test_auxiliary_derivatives()
   print '(a)', 'pbkrtest tests passed'

contains

   subroutine test_sigma_g()
      type(random_sigma_term_t), allocatable :: terms(:)
      type(sigma_g_result_t) :: result
      real(dp) :: expected_sigma(4, 4)
      real(dp) :: expected_gamma(4)
      integer :: status

      allocate(terms(1))
      terms(1)%n_levels = 2
      allocate(terms(1)%z(4, 4), terms(1)%covariance(2, 2))
      terms(1)%z(1, :) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
      terms(1)%z(2, :) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
      terms(1)%z(3, :) = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
      terms(1)%z(4, :) = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
      terms(1)%covariance = reshape([2.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2])
      call build_sigma_g(terms, 0.5_dp, result, status)
      call assert_int(status, pbkr_success, 'build_sigma_g status')
      expected_gamma = [2.0_dp, 0.3_dp, 1.0_dp, 0.5_dp]
      expected_sigma(1, :) = [4.5_dp, 0.6_dp, 2.3_dp, 2.3_dp]
      expected_sigma(2, :) = [0.6_dp, 2.5_dp, 1.3_dp, 1.3_dp]
      expected_sigma(3, :) = [2.3_dp, 1.3_dp, 4.1_dp, 0.0_dp]
      expected_sigma(4, :) = [2.3_dp, 1.3_dp, 0.0_dp, 4.1_dp]
      call assert_vector_close(result%gamma, expected_gamma, 1.0e-13_dp, 'Sigma gamma')
      call assert_matrix_close(result%sigma, expected_sigma, 1.0e-13_dp, 'Sigma matrix')
      call assert_matrix_close(result%g(:, :, 4), identity_matrix(4), 1.0e-13_dp, 'residual G')
   end subroutine test_sigma_g

   subroutine test_column_spaces()
      real(dp), allocatable :: basis(:, :)
      real(dp), allocatable :: l_full(:, :)
      real(dp), allocatable :: l_recovered(:, :)
      real(dp), allocatable :: x_small(:, :)
      real(dp) :: empty_w(3, 0)
      real(dp) :: l(1, 2)
      real(dp) :: l_redundant(2, 2)
      real(dp) :: x(4, 2)
      real(dp) :: x_intercept(4, 1)
      integer :: relationship
      integer :: status

      x(:, 1) = 1.0_dp
      x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      x_intercept(:, 1) = 1.0_dp
      l = reshape([0.0_dp, 1.0_dp], [1, 2])
      call compare_column_space(x, x_intercept, relationship, status)
      call assert_int(status, pbkr_success, 'compare column status')
      call assert_int(relationship, 1, 'nested column relationship')
      call compare_column_space(x_intercept, x, relationship, status)
      call assert_int(relationship, 0, 'reverse nested column relationship')
      call compare_column_space(x, x, relationship, status)
      call assert_int(relationship, -1, 'equal column relationship')

      call orthogonal_complement(empty_w, basis, status)
      call assert_int(status, pbkr_success, 'empty orthogonal complement status')
      call assert_matrix_close(basis, identity_matrix(3), 1.0e-13_dp, 'empty complement identity')

      call orthogonal_complement(transpose(l), basis, status)
      call assert_int(status, pbkr_success, 'orthogonal complement status')
      call assert_int(size(basis, 2), 1, 'orthogonal complement dimension')
      call assert_close(abs(basis(1, 1)), 1.0_dp, 1.0e-12_dp, 'orthogonal complement first coordinate')
      call assert_close(basis(2, 1), 0.0_dp, 1.0e-12_dp, 'orthogonal complement second coordinate')

      call make_model_matrix(x, l, x_small, status)
      call assert_int(status, pbkr_success, 'make_model_matrix status')
      call assert_int(size(x_small, 2), 1, 'small model columns')
      call assert_close(maxval(abs(x_small(:, 1) - x_small(1, 1))), 0.0_dp, 1.0e-12_dp, &
         'small model intercept span')

      call make_restriction_matrix(x, x_intercept, l_recovered, status)
      call assert_int(status, pbkr_success, 'make_restriction_matrix status')
      call assert_int(size(l_recovered, 1), 1, 'restriction rows')
      call assert_close(l_recovered(1, 1), 0.0_dp, 1.0e-11_dp, 'restriction intercept coefficient')
      call assert_close(abs(l_recovered(1, 2)), 1.0_dp, 1.0e-11_dp, 'restriction slope coefficient')

      l_redundant(1, :) = [0.0_dp, 1.0_dp]
      l_redundant(2, :) = [0.0_dp, 2.0_dp]
      call force_full_rank(l_redundant, l_full, status)
      call assert_int(status, pbkr_success, 'force_full_rank status')
      call assert_int(size(l_full, 1), 1, 'force_full_rank rows')
      call assert_close(abs(l_full(1, 2)), 1.0_dp, 1.0e-11_dp, 'force_full_rank row space')
   end subroutine test_column_spaces

   subroutine test_hypothesis_helpers()
      type(bootstrap_result_t) :: boot
      type(lrt_result_t) :: lrt
      real(dp) :: ref(8)
      integer :: status

      call likelihood_ratio_test(-10.0_dp, -12.0_dp, 5, 3, lrt, status)
      call assert_int(status, pbkr_success, 'LRT status')
      call assert_close(lrt%statistic, 4.0_dp, 1.0e-13_dp, 'LRT statistic')
      call assert_int(lrt%df, 2, 'LRT df')
      call assert_close(lrt%p_value, 0.1353352832366127_dp, 2.0e-11_dp, 'LRT p-value')

      ref = [-1.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      call bootstrap_p_values(3.2_dp, 2, ref, boot, status)
      call assert_int(status, pbkr_success, 'bootstrap status')
      call assert_int(boot%nsim, 8, 'bootstrap nsim')
      call assert_int(boot%npos, 7, 'bootstrap npos')
      call assert_int(boot%n_extreme, 3, 'bootstrap extreme count')
      call assert_close(boot%mean_positive, 3.0714285714285716_dp, 1.0e-13_dp, 'bootstrap mean')
      call assert_close(boot%variance_positive, 4.2023809523809526_dp, 1.0e-13_dp, 'bootstrap variance')
      call assert_close(boot%p_bootstrap, 0.5_dp, 1.0e-13_dp, 'bootstrap p')
      call assert_close(boot%p_bootstrap_all, 0.4444444444444444_dp, 1.0e-13_dp, 'bootstrap all p')
      call assert_close(boot%standard_error, 0.1890_dp, 1.0e-13_dp, 'bootstrap se')
      call assert_close(boot%ci_low, 0.1296_dp, 1.0e-13_dp, 'bootstrap CI low')
      call assert_close(boot%ci_high, 0.8704_dp, 1.0e-13_dp, 'bootstrap CI high')
      call assert_close(boot%bartlett_statistic, 2.083720930232558_dp, 2.0e-12_dp, 'Bartlett statistic')
      call assert_close(boot%gamma_scale, 1.3682170542635659_dp, 2.0e-12_dp, 'gamma scale')
      call assert_close(boot%gamma_shape, 2.244840145690004_dp, 2.0e-12_dp, 'gamma shape')
      call assert_close(boot%f_ddf, 5.733333333333333_dp, 2.0e-12_dp, 'bootstrap F ddf')
   end subroutine test_hypothesis_helpers

   subroutine test_kr_kernels()
      type(kr_result_t) :: kr
      type(vcov_adjustment_t) :: adj
      real(dp) :: beta(2)
      real(dp) :: beta_h(2)
      real(dp) :: ddf
      real(dp) :: g(6, 6, 2)
      real(dp) :: l(1, 2)
      real(dp) :: phi(2, 2)
      real(dp) :: sigma(6, 6)
      real(dp) :: x(6, 2)
      integer :: i
      integer :: status

      x(:, 1) = 1.0_dp
      x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      g = 0.0_dp
      do i = 1, 6
         g(i, i, 2) = 1.0_dp
      end do
      g(1, 1, 1) = 1.0_dp
      g(2, 2, 1) = 1.0_dp
      g(3, 3, 1) = 2.0_dp
      g(4, 4, 1) = 2.0_dp
      g(5, 5, 1) = 3.0_dp
      g(6, 6, 1) = 3.0_dp
      sigma = 0.8_dp * g(:, :, 1) + 1.2_dp * g(:, :, 2)
      phi = reshape([ &
         1.1306557991461641_dp, -0.3273093840063068_dp, &
         -0.3273093840063068_dp, 0.1552412667094590_dp], [2, 2])
      call vcov_adjust_kr(phi, sigma, g, x, adj, status)
      call assert_int(status, pbkr_success, 'vcov_adjust_kr status')
      call assert_close(adj%condition, 0.0621213417583446_dp, 2.0e-10_dp, 'KR information condition')
      call assert_close(adj%w(1, 1), 7.5121849461909385_dp, 2.0e-9_dp, 'KR W11')
      call assert_close(adj%w(1, 2), -12.86672560749181_dp, 3.0e-9_dp, 'KR W12')
      call assert_close(adj%w(2, 2), 25.487867591616922_dp, 3.0e-9_dp, 'KR W22')
      call assert_close(adj%phi_adjusted(1, 1), 1.768154429001076_dp, 3.0e-9_dp, 'KR adjusted covariance 11')
      call assert_close(adj%phi_adjusted(1, 2), -0.5699506224507946_dp, 3.0e-9_dp, 'KR adjusted covariance 12')
      call assert_close(adj%phi_adjusted(2, 2), 0.2573122005246890_dp, 3.0e-9_dp, 'KR adjusted covariance 22')

      l = reshape([0.0_dp, 1.0_dp], [1, 2])
      beta = [1.0_dp, 0.5_dp]
      beta_h = 0.0_dp
      call kr_adjust(adj%phi_adjusted, phi, l, beta, beta_h, adj%p_matrices, adj%w, kr, status)
      call assert_int(status, pbkr_success, 'kr_adjust status')
      call assert_int(kr%ndf, 1, 'KR ndf')
      call assert_close(kr%a1, 0.5036337843832914_dp, 3.0e-9_dp, 'KR A1')
      call assert_close(kr%a2, 0.5036337843832914_dp, 3.0e-9_dp, 'KR A2')
      call assert_close(kr%ddf, 3.971139470814167_dp, 3.0e-9_dp, 'KR ddf')
      call assert_close(kr%f_stat_unscaled, 0.9715823792662043_dp, 3.0e-9_dp, 'KR F')
      call assert_close(kr%wald_unadjusted, 1.610396547896548_dp, 3.0e-9_dp, 'KR unadjusted Wald')
      call assert_close(kr%p_value, 0.3804662483828975_dp, 3.0e-8_dp, 'KR p-value')

      call lb_ddf(l, phi, adj%p_matrices, adj%w, ddf, status)
      call assert_int(status, pbkr_success, 'lb_ddf status')
      call assert_close(ddf, kr%ddf, 3.0e-9_dp, 'Lb ddf equals KR ddf')
      call ddf_lb_scalar(phi, l(1, :), adj%p_matrices, adj%w, ddf, status)
      call assert_int(status, pbkr_success, 'ddf_lb_scalar status')
      call assert_close(ddf, kr%ddf, 3.0e-9_dp, 'scalar Lb ddf equals KR ddf')
   end subroutine test_kr_kernels

   subroutine test_satterthwaite_kernels()
      type(satterthwaite_result_t) :: sat
      real(dp) :: beta(2)
      real(dp) :: beta_h(2)
      real(dp) :: jac(2, 2, 2)
      real(dp) :: l(2, 2)
      real(dp) :: nu(3)
      real(dp) :: vcov_beta(2, 2)
      real(dp) :: vcov_varpar(2, 2)
      integer :: status

      beta = [1.0_dp, 0.5_dp]
      beta_h = 0.0_dp
      l = identity_matrix(2)
      vcov_beta = reshape([0.5_dp, 0.1_dp, 0.1_dp, 0.3_dp], [2, 2])
      vcov_varpar = reshape([0.04_dp, 0.0_dp, 0.0_dp, 0.09_dp], [2, 2])
      jac(:, :, 1) = reshape([0.2_dp, 0.05_dp, 0.05_dp, 0.1_dp], [2, 2])
      jac(:, :, 2) = reshape([-0.1_dp, 0.02_dp, 0.02_dp, 0.15_dp], [2, 2])
      call satterthwaite_test(beta, beta_h, l, vcov_beta, vcov_varpar, jac, sat, status)
      call assert_int(status, pbkr_success, 'Satterthwaite status')
      call assert_int(sat%ndf, 2, 'Satterthwaite ndf')
      call assert_close(sat%f_stat, 1.160714285714286_dp, 2.0e-12_dp, 'Satterthwaite F')
      call assert_close(sat%nu(1), 270.57202937_dp, 3.0e-6_dp, 'Satterthwaite nu1')
      call assert_close(sat%nu(2), 117.51418251_dp, 3.0e-6_dp, 'Satterthwaite nu2')
      call assert_close(sat%ddf, 163.5464312898107_dp, 3.0e-7_dp, 'Satterthwaite ddf')
      call assert_close(sat%p_value, 0.3158292215756831_dp, 3.0e-8_dp, 'Satterthwaite p-value')

      nu = [10.0_dp, 10.0_dp + 1.0e-10_dp, 10.0_dp + 2.0e-10_dp]
      call assert_close(get_fstat_ddf(nu), sum(nu) / 3.0_dp, 1.0e-12_dp, 'equal-nu shortcut')
      nu = [1.9_dp, 10.0_dp, 20.0_dp]
      call assert_close(get_fstat_ddf(nu), 2.0_dp, 1.0e-13_dp, 'small-nu shortcut')
   end subroutine test_satterthwaite_kernels

   subroutine test_auxiliary_derivatives()
      type(auxiliary_callbacks_t) :: callbacks
      type(auxiliary_result_t) :: aux
      real(dp) :: x(2)
      integer :: status

      callbacks%deviance => quadratic_deviance
      callbacks%covbeta_vector => covariance_vector
      x = [1.0_dp, -0.5_dp]
      call compute_auxiliary_numeric(x, 2, callbacks, aux, status)
      call assert_int(status, pbkr_success, 'auxiliary derivative status')
      call assert_matrix_close(aux%hessian, reshape([2.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [2, 2]), &
         3.0e-4_dp, 'auxiliary Hessian')
      call assert_matrix_close(aux%vcov_varpar, reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [2, 2]), &
         3.0e-4_dp, 'auxiliary variance covariance')
      call assert_matrix_close(aux%jacobian(:, :, 1), &
         reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2]), 3.0e-5_dp, 'auxiliary Jacobian 1')
      call assert_matrix_close(aux%jacobian(:, :, 2), &
         reshape([0.0_dp, 0.1_dp, 0.1_dp, 0.0_dp], [2, 2]), 3.0e-5_dp, 'auxiliary Jacobian 2')
      call assert_int(aux%negative_eigenvalues, 0, 'auxiliary negative modes')
      call assert_int(aux%near_zero_eigenvalues, 0, 'auxiliary zero modes')
   end subroutine test_auxiliary_derivatives

   pure real(dp) function quadratic_deviance(x) result(value)
      real(dp), intent(in) :: x(:) !! Two variance parameters used in the deterministic quadratic test objective.

      value = (x(1) - 1.0_dp) ** 2 + 2.0_dp * (x(2) + 0.5_dp) ** 2
   end function quadratic_deviance

   function covariance_vector(x) result(value)
      real(dp), intent(in) :: x(:) !! Two variance parameters defining the deterministic covariance test function.
      real(dp), allocatable :: value(:)

      allocate(value(4))
      value = [1.0_dp + x(1), 0.1_dp * x(2), 0.1_dp * x(2), 2.0_dp + 2.0_dp * x(1)]
   end function covariance_vector

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n !! Order of the square identity matrix.
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed scalar value.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference.
      character(len=*), intent(in) :: label !! Test label printed if the assertion fails.

      if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_int(actual, expected, label)
      integer, intent(in) :: actual !! Computed integer value.
      integer, intent(in) :: expected !! Reference integer value.
      character(len=*), intent(in) :: label !! Test label printed if the assertion fails.

      if (actual /= expected) then
         print '(a,2(1x,i0))', trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_int

   subroutine assert_vector_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Computed vector.
      real(dp), intent(in) :: expected(:) !! Reference vector with the same shape.
      real(dp), intent(in) :: tolerance !! Maximum allowed elementwise absolute difference.
      character(len=*), intent(in) :: label !! Test label printed if the assertion fails.

      if (size(actual) /= size(expected)) then
         print '(a)', trim(label) // ': shape mismatch'
         error stop 1
      end if
      if (any(ieee_is_nan(actual)) .or. maxval(abs(actual - expected)) > tolerance) then
         print '(a,1x,es24.16)', trim(label), maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_vector_close

   subroutine assert_matrix_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:, :) !! Computed matrix.
      real(dp), intent(in) :: expected(:, :) !! Reference matrix with the same shape.
      real(dp), intent(in) :: tolerance !! Maximum allowed elementwise absolute difference.
      character(len=*), intent(in) :: label !! Test label printed if the assertion fails.

      if (any(shape(actual) /= shape(expected))) then
         print '(a)', trim(label) // ': shape mismatch'
         error stop 1
      end if
      if (any(ieee_is_nan(actual)) .or. maxval(abs(actual - expected)) > tolerance) then
         print '(a,1x,es24.16)', trim(label), maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_matrix_close

end program test_pbkrtest
