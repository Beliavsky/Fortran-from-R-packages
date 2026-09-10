! SPDX-License-Identifier: GPL-2.0-only
program test_marss
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use marss_api
   use kfas, only : kfas_model
   implicit none

   integer :: failures

   failures = 0
   call test_utilities(failures)
   call test_initialization(failures)
   call test_model_builder(failures)
   call test_constraint_initialization(failures)
   call test_high_level_workflow(failures)
   call test_kemcheck_degeneracy(failures)
   call test_kalman_reference(failures)
   call test_kfas_bridge(failures)
   call test_kfas_edge_conversion(failures)
   call test_kfas_backend_parity(failures)
   call test_time_varying_deterministic(failures)
   call test_marxss_numerical_form(failures)
   call test_rich_constraints(failures)
   call test_constraint_labels(failures)
   call test_dfa_model(failures)
   call test_dfa_spec(failures)
   call test_vectorization(failures)
   call test_free_vectorization(failures)
   call test_simulation(failures)
   call test_analysis_helpers(failures)
   call test_residual_hierarchy(failures)
   call test_harvey_residuals(failures)
   call test_hatyt_correlated_missing(failures)
   call test_hatyt_full(failures)
   call test_em_fit(failures)
   call test_missing_em(failures)
   call test_time_varying_fixed_em(failures)
   call test_bootstrap(failures)
   call test_constrained_bootstrap(failures)
   call test_bootstrap_aic(failures)
   call test_innovations_bootstrap(failures)
   call test_parameter_cis(failures)
   call test_bootstrap_parameter_cis(failures)
   call test_harvey_information(failures)
   call test_hessian_summary(failures)
   call test_optim_fit(failures)
   call test_linear_constraints(failures)
   call test_constrained_em(failures)
   call test_cross_validation(failures)
   call test_constrained_cross_validation(failures)
   if (failures /= 0) then
      write (*, '(a,i0)') 'MARSS tests failed: ', failures
      error stop 1
   end if
   write (*, '(a)') 'All MARSS deterministic tests passed.'

contains

   subroutine test_utilities(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      real(dp) :: x(3)
      real(dp) :: zs(3)
      real(dp) :: v(2)
      real(dp), allocatable :: d(:, :)

      x = [1.0_dp, 2.0_dp, 3.0_dp]
      zs = zscore(x)
      call assert_vector_close('zscore vector', zs, [-1.0_dp, 0.0_dp, 1.0_dp], 1.0e-12_dp, failures)
      d = ldiag(2.0_dp, 2)
      call assert_matrix_close('ldiag scalar', d, reshape([2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp], [2, 2]), &
         1.0e-12_dp, failures)
      v = [1.0_dp, 3.0_dp]
      d = ldiag(v)
      call assert_matrix_close('ldiag vector', d, reshape([1.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2]), &
         1.0e-12_dp, failures)
   end subroutine test_utilities

   subroutine test_initialization(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      real(dp) :: y(2, 3)
      integer :: info
      logical :: ok

      y = reshape([1.0_dp, 2.0_dp, 1.2_dp, 2.2_dp, 1.4_dp, 2.4_dp], [2, 3])
      call marss_inits(y, 2, model, 0)
      call marss_kemcheck(model, ok, info)
      call assert_true('MARSSinits/MARSSkemcheck', ok .and. info == 0, failures)
      call assert_close('MARSSinits B(1,1)', model%b(1, 1), 1.0_dp, 1.0e-12_dp, failures)
      call assert_close('MARSSinits Q(1,1)', model%q(1, 1), 0.05_dp, 1.0e-12_dp, failures)
   end subroutine test_initialization

   subroutine test_model_builder(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model_spec) :: spec
      type(marss_model) :: model
      type(marss_constraints) :: constraints
      real(dp) :: y(2, 4)
      real(dp) :: cdata(1, 4)
      real(dp) :: ddata(2, 1)
      real(dp), allocatable :: theta(:)
      integer :: info

      y = reshape([1.0_dp, 2.0_dp, 1.1_dp, 2.1_dp, 1.2_dp, 2.2_dp, 1.3_dp, 2.3_dp], [2, 4])
      call marss_build(y, spec, model, constraints, info)
      call assert_true('MARSS builder default status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSS builder default free count', &
            marss_constraints_parameter_count(constraints) == 7, failures)
         call assert_matrix_close('MARSS builder identity Z', model%z, &
            reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 0.0_dp, failures)
         call assert_true('MARSS builder scaling A fixed zero', &
            size(constraints%a%design, 2) == 0 .and. maxval(abs(model%a)) <= 1.0e-15_dp, failures)
         call assert_true('MARSS builder zero V0', maxval(abs(model%v0)) <= 1.0e-15_dp, failures)
         call marss_vectorize_free(model, constraints, theta, info)
         call assert_true('MARSS builder beta projection', info == 0 .and. size(theta) == 7, failures)
      end if

      spec = marss_model_spec()
      spec%z = 'onestate'
      call marss_build(y, spec, model, constraints, info)
      call assert_true('MARSS builder onestate status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSS builder onestate dimension', size(model%b, 1) == 1, failures)
         call assert_true('MARSS builder onestate Z', maxval(abs(model%z - 1.0_dp)) <= 1.0e-15_dp, failures)
         call assert_true('MARSS builder scaling A references', size(constraints%a%design, 2) == 1, failures)
         call assert_true('MARSS builder onestate free count', &
            marss_constraints_parameter_count(constraints) == 5, failures)
      end if

      spec = marss_model_spec()
      spec%form = 'marxss'
      spec%c = 'unconstrained'
      spec%d = 'diagonal and equal'
      cdata(1, :) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      ddata(:, 1) = [1.0_dp, -1.0_dp]
      call marss_build(y, spec, model, constraints, info, cdata, ddata)
      call assert_true('MARSS builder marxss covariate status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSS builder state covariates', allocated(model%c_coef) .and. &
            all(shape(model%c_coef) == [2, 1]), failures)
         call assert_true('MARSS builder observation covariate recycling', allocated(model%obs_covariates) .and. &
            size(model%obs_covariates, 2) == 4 .and. &
            maxval(abs(model%obs_covariates(:, 1) - model%obs_covariates(:, 4))) <= 1.0e-15_dp, failures)
         call assert_true('MARSS builder covariate free count', &
            marss_constraints_parameter_count(constraints) == 10, failures)
      end if

      spec = marss_model_spec()
      spec%form = 'dfa'
      call marss_build(y, spec, model, constraints, info)
      call assert_true('MARSS builder rejects unsupported form', info == 5, failures)
   end subroutine test_model_builder

   subroutine test_constraint_initialization(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model_spec) :: spec
      type(marss_model) :: template
      type(marss_model) :: initialized_model
      type(marss_model) :: time_template
      type(marss_model) :: time_model
      type(marss_constraints) :: constraints
      type(marss_constraints) :: initialized
      type(marss_constraints) :: time_constraints
      type(marss_constraints) :: time_initialized
      real(dp) :: y(2, 4)
      character(len=16) :: init_names(2)
      integer :: info

      y = reshape([1.0_dp, 2.0_dp, 1.1_dp, 2.1_dp, 1.2_dp, 2.2_dp, 1.3_dp, 2.3_dp], [2, 4])
      spec%b = 'unconstrained'
      spec%q = 'unconstrained'
      call marss_build(y, spec, template, constraints, info)
      call assert_true('MARSSinits constrained setup', info == 0, failures)
      if (info /= 0) return
      call marss_inits_linear(template, constraints, initialized, initialized_model, info)
      call assert_true('MARSSinits constrained status', info == 0, failures)
      if (info == 0) then
         call assert_matrix_close('MARSSinits constrained B diagonal target', initialized_model%b, &
            reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 1.0e-10_dp, failures)
         call assert_matrix_close('MARSSinits constrained Q diagonal target', initialized_model%q_noise, &
            0.05_dp * reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 1.0e-10_dp, failures)
         call assert_vector_close('MARSSinits constrained x0 solve', initialized_model%x0, &
            [1.0_dp, 2.0_dp], 1.0e-10_dp, failures)
      end if

      init_names = [character(len=16) :: 'B.1', 'Q.1']
      call marss_inits_named(template, constraints, init_names, [0.7_dp, 0.2_dp], initialized, initialized_model, info)
      call assert_true('MARSSinits partial named override status', info == 0, failures)
      if (info == 0) then
         call assert_close('MARSSinits named B override', initialized_model%b(1, 1), 0.7_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSinits unnamed B default retained', initialized_model%b(2, 2), 1.0_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSinits named Q override', initialized_model%q_noise(1, 1), 0.2_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSinits unnamed Q default retained', initialized_model%q_noise(2, 2), &
            0.05_dp, 1.0e-12_dp, failures)
      end if

      spec = marss_model_spec()
      spec%z = 'onestate'
      call marss_build(y, spec, template, constraints, info)
      if (info /= 0) then
         call assert_true('MARSSinits onestate setup', .false., failures)
         return
      end if
      call marss_inits_linear(template, constraints, initialized, initialized_model, info)
      call assert_true('MARSSinits onestate status', info == 0, failures)
      if (info == 0) then
         call assert_close('MARSSinits onestate x0 least squares', initialized_model%x0(1), &
            1.5_dp, 1.0e-10_dp, failures)
      end if

      allocate(time_template%y(1, 2), time_template%b(1, 1), time_template%u(1), time_template%q(1, 1))
      allocate(time_template%z(1, 1), time_template%a(1), time_template%r(1, 1))
      allocate(time_template%x0(1), time_template%v0(1, 1), time_template%b_t(1, 1, 2))
      time_template%y = 0.0_dp
      time_template%b = 0.5_dp
      time_template%b_t = 0.5_dp
      time_template%u = 0.0_dp
      time_template%q = 0.1_dp
      time_template%z = 1.0_dp
      time_template%a = 0.0_dp
      time_template%r = 0.2_dp
      time_template%x0 = 0.0_dp
      time_template%v0 = 1.0_dp
      time_template%tinitx = 1
      allocate(time_constraints%b%fixed(2), time_constraints%b%design(2, 1), time_constraints%b%start(1))
      time_constraints%b%fixed = 0.0_dp
      time_constraints%b%design(:, 1) = [1.0_dp, 2.0_dp]
      time_constraints%b%start = 0.1_dp
      call marss_inits_linear(time_template, time_constraints, time_initialized, time_model, info)
      call assert_true('MARSSinits time-average projection status', info == 0, failures)
      if (info == 0) then
         call assert_close('MARSSinits time-average scalar projection', &
            time_initialized%b%start(1), 2.0_dp / 3.0_dp, 1.0e-12_dp, failures)
      end if
   end subroutine test_constraint_initialization

   subroutine test_high_level_workflow(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model_spec) :: spec
      type(marss_constraints) :: constraints
      type(marss_fit_result) :: fit
      real(dp) :: y(1, 8)
      real(dp) :: cdata(1, 8)
      integer :: info

      y(1, :) = [0.2_dp, 0.4_dp, 0.35_dp, 0.6_dp, 0.55_dp, 0.7_dp, 0.65_dp, 0.8_dp]
      cdata(1, :) = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
      spec%form = 'marxss'
      spec%z = 'onestate'
      spec%a = 'zero'
      spec%b = 'diagonal and equal'
      spec%u = 'zero'
      spec%c = 'equal'
      spec%v0 = 'zero'
      call marss_from_data(y, spec, 4, 1.0e-5_dp, fit, constraints, 'kem', info, cdata)
      call assert_true('MARSS high-level KEM workflow', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSS high-level beta output', allocated(fit%free_parameters), failures)
         call assert_true('MARSS high-level finite likelihood', fit%loglik > -huge(1.0_dp), failures)
      end if
      call marss_from_data(y, spec, 3, 1.0e-5_dp, fit, constraints, 'bfgs', info, cdata)
      call assert_true('MARSS high-level BFGS workflow', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSS high-level BFGS beta output', allocated(fit%free_parameters), failures)
      end if
   end subroutine test_high_level_workflow

   subroutine test_kemcheck_degeneracy(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: unit_model
      type(marss_model) :: complex_model
      type(marss_model) :: adjacency_model
      type(marss_model) :: classification_model
      type(marss_model) :: indirectly_stochastic_model
      type(marss_model) :: consistency_model
      type(marss_model) :: time_model
      type(marss_model) :: short_model
      type(marss_constraints) :: fixed_constraints
      type(marss_constraints) :: z_free
      type(marss_constraints) :: a_free
      type(marss_constraints) :: b_free
      type(marss_constraints) :: u_free
      type(marss_constraints) :: indirect_u_free
      type(marss_constraints) :: indirect_c_free
      type(marss_constraints) :: x0_free
      type(marss_constraints) :: free_b_constraints
      logical :: ok
      integer :: info
      integer :: t
      real(dp) :: y2(2, 4)

      call make_reference_model(short_model)
      deallocate(short_model%y)
      allocate(short_model%y(1, 2))
      short_model%y = 0.0_dp
      call marss_kemcheck(short_model, ok, info)
      call assert_true('MARSSkemcheck minimum time length', .not. ok .and. info == 2, failures)

      call make_reference_model(unit_model)
      unit_model%b(1, 1) = 1.0_dp + 0.5_dp * sqrt(epsilon(1.0_dp))
      call marss_kemcheck(unit_model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck fixed B unit-circle tolerance', ok .and. info == 0, failures)
      unit_model%b(1, 1) = 1.01_dp
      call marss_kemcheck(unit_model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck fixed B unit circle', .not. ok .and. info == 11, failures)
      call make_free_block(free_b_constraints%b, reshape(unit_model%b, [size(unit_model%b)]), 1)
      call marss_kemcheck(unit_model, ok, info, free_b_constraints)
      call assert_true('MARSSkemcheck estimated B may start outside unit circle', ok .and. info == 0, failures)

      y2 = reshape([0.1_dp, -0.2_dp, 0.2_dp, -0.1_dp, 0.3_dp, 0.0_dp, 0.4_dp, 0.1_dp], [2, 4])
      call marss_inits(y2, 2, complex_model, 0)
      complex_model%b = reshape([0.0_dp, 1.02_dp, -1.02_dp, 0.0_dp], [2, 2])
      call marss_kemcheck(complex_model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck fixed B complex eigenvalue modulus', .not. ok .and. info == 11, failures)

      call make_reference_model(model)
      model%r = 0.0_dp
      call make_free_block(z_free%z, reshape(model%z, [size(model%z)]), 1)
      call marss_kemcheck(model, ok, info, z_free)
      call assert_true('MARSSkemcheck R0 requires fixed Z', .not. ok .and. info == 20, failures)

      call make_free_block(a_free%a, model%a, 1)
      call marss_kemcheck(model, ok, info, a_free)
      call assert_true('MARSSkemcheck R0 requires fixed A', .not. ok .and. info == 21, failures)

      model%r = 0.4_dp
      model%q = 0.0_dp
      call make_free_block(b_free%b, reshape(model%b, [size(model%b)]), 1)
      call marss_kemcheck(model, ok, info, b_free)
      call assert_true('MARSSkemcheck Q0 requires fixed B row', .not. ok .and. info == 31, failures)

      model%r = 0.0_dp
      call make_free_block(u_free%u, model%u, 1)
      call marss_kemcheck(model, ok, info, u_free)
      call assert_true('MARSSkemcheck linked R0 Z Q0 requires fixed U', .not. ok .and. info == 33, failures)

      model%y(1, :) = [0.62_dp, 0.716_dp, 0.7928_dp, 0.85424_dp]
      call marss_kemcheck(model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck fixed degenerate model', ok .and. info == 0, failures)

      call make_reference_model(consistency_model)
      consistency_model%y = 0.0_dp
      consistency_model%b = 0.5_dp
      consistency_model%u = 0.0_dp
      consistency_model%q = 0.0_dp
      consistency_model%z = 1.0_dp
      consistency_model%a = 0.0_dp
      consistency_model%r = 0.0_dp
      consistency_model%x0 = 0.0_dp
      consistency_model%v0 = 0.0_dp
      call marss_kemcheck(consistency_model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck R0 fitted-value consistency', ok .and. info == 0, failures)
      consistency_model%y(1, 2) = 1.0_dp
      call marss_kemcheck(consistency_model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck R0 fitted-value inconsistency', .not. ok .and. info == 22, failures)

      allocate(time_model%y(1, 4), time_model%b(2, 2), time_model%u(2), time_model%q(2, 2))
      allocate(time_model%z(1, 2), time_model%a(1), time_model%r(1, 1), time_model%x0(2), time_model%v0(2, 2))
      allocate(time_model%q_t(2, 2, 4))
      time_model%y = 0.0_dp
      time_model%b = 0.0_dp
      time_model%b(1, 1) = 0.5_dp
      time_model%b(2, 2) = 0.5_dp
      time_model%u = 0.0_dp
      time_model%q = 0.0_dp
      time_model%q(1, 1) = 1.0_dp
      time_model%q(2, 2) = 1.0_dp
      time_model%z = reshape([1.0_dp, 0.0_dp], [1, 2])
      time_model%a = 0.0_dp
      time_model%r = 0.2_dp
      time_model%x0 = 0.0_dp
      time_model%v0 = 0.0_dp
      time_model%v0(1, 1) = 1.0_dp
      time_model%v0(2, 2) = 1.0_dp
      do t = 1, 4
         time_model%q_t(:, :, t) = 0.0_dp
      end do
      time_model%q_t(2, 2, 1) = 1.0_dp
      time_model%q_t(2, 2, 2) = 1.0_dp
      time_model%q_t(1, 1, 3) = 1.0_dp
      time_model%q_t(1, 1, 4) = 1.0_dp
      call marss_kemcheck(time_model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck Q0 placement time constant', .not. ok .and. info == 30, failures)

      call make_kem_graph_model(adjacency_model)
      adjacency_model%b_t(2, 1, 2) = 0.25_dp
      call make_free_block(x0_free%x0, adjacency_model%x0, 1)
      call marss_kemcheck(adjacency_model, ok, info, x0_free)
      call assert_true('MARSSkemcheck estimated deterministic x0 requires fixed B adjacency', &
         .not. ok .and. info == 34, failures)

      call make_kem_graph_model(classification_model)
      classification_model%b_t(1, 2, 2) = 0.25_dp
      call marss_kemcheck(classification_model, ok, info, fixed_constraints)
      call assert_true('MARSSkemcheck deterministic state class time constant', &
         .not. ok .and. info == 35, failures)

      call make_kem_graph_model(indirectly_stochastic_model)
      indirectly_stochastic_model%b(1, 2) = 0.25_dp
      do t = 1, 4
         indirectly_stochastic_model%b_t(:, :, t) = indirectly_stochastic_model%b
      end do
      call make_free_block(indirect_u_free%u, indirectly_stochastic_model%u, 1)
      call marss_kemcheck(indirectly_stochastic_model, ok, info, indirect_u_free)
      call assert_true('MARSSkemcheck indirectly stochastic state requires fixed U', &
         .not. ok .and. info == 37, failures)

      allocate(indirectly_stochastic_model%state_covariates(1, 4))
      allocate(indirectly_stochastic_model%c_coef(2, 1))
      indirectly_stochastic_model%state_covariates = 1.0_dp
      indirectly_stochastic_model%c_coef = 0.0_dp
      call make_free_block(indirect_c_free%c, reshape(indirectly_stochastic_model%c_coef, &
         [size(indirectly_stochastic_model%c_coef)]), 1)
      call marss_kemcheck(indirectly_stochastic_model, ok, info, indirect_c_free)
      call assert_true('MARSSkemcheck indirectly stochastic state requires fixed C row', &
         .not. ok .and. info == 37, failures)
   end subroutine test_kemcheck_degeneracy

   subroutine test_kalman_reference(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      real(dp) :: nan

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      call make_reference_model(model)
      model%y(1, :) = [0.9_dp, 0.2_dp, nan, 1.1_dp]
      call marss_kfss(model, kf)
      call assert_true('MARSSkfss status', kf%ok .and. kf%info == 0, failures)
      call assert_close('MARSSkfss logLik', kf%loglik, -3.376636469758858_dp, 2.0e-8_dp, failures)
      call assert_vector_close('MARSSkfss xpred', kf%x_pred(1, :), &
         [0.6_dp, 0.81611785_dp, 0.57783645_dp, 0.66226916_dp], 2.0e-8_dp, failures)
      call assert_vector_close('MARSSkfss xfilt', kf%x_filt(1, :), &
         [0.77014731_dp, 0.47229557_dp, 0.57783645_dp, 0.88816726_dp], 2.0e-8_dp, failures)
      call assert_vector_close('MARSSkfss xsmooth', kf%x_smooth(1, :), &
         [0.65686347_dp, 0.51576422_dp, 0.70923486_dp, 0.88816726_dp], 2.0e-8_dp, failures)
      call assert_close('MARSSkfss x0smooth', kf%x0_smooth(1), 0.5425715817251424_dp, 2.0e-8_dp, failures)
      call assert_close('MARSSkfss v0smooth', kf%v0_smooth(1, 1), 0.3723675757698826_dp, 2.0e-8_dp, failures)
      call assert_true('MARSSkfss missing innovation', abs(kf%innov(1, 3)) < 1.0e-14_dp, failures)
   end subroutine test_kalman_reference

   subroutine test_kfas_bridge(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: native
      type(marss_kf_result) :: dispatched
      type(kfas_model) :: kmodel
      real(dp) :: nan
      integer :: info

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(model%y(1, 3), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
      model%y(1, :) = [0.9_dp, nan, 0.4_dp]
      model%b = 0.8_dp
      model%u = 0.2_dp
      model%q = 0.3_dp
      model%z = 1.2_dp
      model%a = -0.1_dp
      model%r = 0.4_dp
      model%x0 = 0.5_dp
      model%v0 = 0.7_dp
      model%tinitx = 1
      model%diffuse = .true.

      call assert_true('diffuse model validation', marss_model_valid(model), failures)
      call marss_to_kfas_model(model, kmodel, info)
      call assert_true('MARSSkfas bridge status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSkfas bridge state dimension', size(kmodel%a1) == 4, failures)
         call assert_true('MARSSkfas bridge missing flag', kmodel%missing(2, 1) == 1, failures)
         call assert_close('MARSSkfas bridge Z', kmodel%z(1, 1, 1), 1.2_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSkfas bridge A', kmodel%z(1, 2, 1), -0.1_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSkfas bridge B', kmodel%tmat(1, 1, 1), 0.8_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSkfas bridge U', kmodel%tmat(1, 2, 1), 0.2_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSkfas bridge lag copy', kmodel%tmat(3, 1, 1), 1.0_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSkfas bridge Q', kmodel%q(1, 1, 1), 0.3_dp, 1.0e-12_dp, failures)
         call assert_close('MARSSkfas bridge P1inf', kmodel%p1inf(1, 1), 0.7_dp, 1.0e-12_dp, failures)
         call assert_true('MARSSkfas bridge diffuse rank', kmodel%diffuse_rank == 1, failures)
      end if
      call marss_kfss(model, native)
      call assert_true('native KF rejects diffuse initialization', .not. native%ok .and. native%info == 2, failures)
      call marss_kf(model, dispatched)
      call assert_true('MARSSkf dispatches diffuse model to KFAS', dispatched%ok .and. dispatched%info == 0, failures)
   end subroutine test_kfas_bridge

   subroutine test_kfas_edge_conversion(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: initial_model
      type(marss_model) :: loading_model
      type(kfas_model) :: kmodel
      real(dp) :: nan
      integer :: info
      integer :: t

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(model%y(2, 3), model%b(2, 2), model%u(2), model%q(2, 2))
      allocate(model%z(2, 2), model%a(2), model%r(2, 2), model%x0(2), model%v0(2, 2))
      allocate(model%g(2, 1), model%q_noise(1, 1))
      model%y = reshape([1.0_dp, 0.2_dp, nan, -0.1_dp, 0.7_dp, 0.3_dp], [2, 3])
      model%b = reshape([0.8_dp, 0.0_dp, 0.1_dp, 0.6_dp], [2, 2])
      model%u = [0.1_dp, -0.2_dp]
      model%q = 0.0_dp
      model%z = reshape([1.0_dp, 0.4_dp, -0.3_dp, 0.8_dp], [2, 2])
      model%a = [0.2_dp, -0.1_dp]
      model%r = reshape([0.4_dp, 0.1_dp, 0.1_dp, 0.3_dp], [2, 2])
      model%x0 = [0.5_dp, -0.25_dp]
      model%v0 = 0.0_dp
      model%v0(1, 1) = 1.0_dp
      model%g(:, 1) = [1.0_dp, 2.0_dp]
      model%q_noise(1, 1) = 0.25_dp
      model%tinitx = 1
      model%diffuse = .true.

      call assert_true('KFAS edge static model valid', marss_model_valid(model), failures)
      call marss_to_kfas_model(model, kmodel, info)
      call assert_true('KFAS edge static conversion status', info == 0, failures)
      if (info == 0) then
         call assert_true('KFAS edge static flags', all(kmodel%time_varying == 0), failures)
         call assert_true('KFAS edge static Z slices', size(kmodel%z, 3) == 1, failures)
         call assert_true('KFAS edge static H slices', size(kmodel%h, 3) == 1, failures)
         call assert_true('KFAS edge static T slices', size(kmodel%tmat, 3) == 1, failures)
         call assert_true('KFAS edge static G slices', size(kmodel%rmat, 3) == 1, failures)
         call assert_true('KFAS edge static Q slices', size(kmodel%q, 3) == 1, failures)
         call assert_true('KFAS edge stacked disturbance dimension', size(kmodel%q, 1) == 2, failures)
         call assert_vector_close('KFAS edge raw G', kmodel%rmat(1:2, 1, 1), &
            [1.0_dp, 2.0_dp], 1.0e-12_dp, failures)
         call assert_close('KFAS edge raw Q', kmodel%q(1, 1, 1), 0.25_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge unused disturbance Q', kmodel%q(2, 2, 1), 0.0_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge correlated H', kmodel%h(1, 2, 1), 0.1_dp, 1.0e-12_dp, failures)
         call assert_true('KFAS edge diffuse rank one', kmodel%diffuse_rank == 1, failures)
         call assert_true('KFAS edge missing observation flag', kmodel%missing(2, 1) == 1, failures)
      end if
      call marss_to_kfas_model(model, kmodel, info, return_lag_one=.false.)
      call assert_true('KFAS edge unstacked conversion status', info == 0, failures)
      if (info == 0) then
         call assert_true('KFAS edge unstacked state dimension', size(kmodel%a1) == 3, failures)
         call assert_true('KFAS edge unstacked disturbance dimension', size(kmodel%q, 1) == 1, failures)
         call assert_true('KFAS edge unstacked Z dimension', size(kmodel%z, 2) == 3, failures)
      end if

      allocate(model%b_t(2, 2, 3), model%u_t(2, 3), model%z_t(2, 2, 3))
      allocate(model%a_t(2, 3), model%r_t(2, 2, 3))
      allocate(model%g_t(2, 1, 3), model%q_noise_t(1, 1, 3))
      do t = 1, 3
         model%b_t(:, :, t) = model%b
         model%b_t(1, 1, t) = 0.5_dp + 0.1_dp * real(t, dp)
         model%u_t(:, t) = [0.1_dp * real(t, dp), -0.05_dp * real(t, dp)]
         model%z_t(:, :, t) = model%z
         model%z_t(2, 2, t) = 0.7_dp + 0.1_dp * real(t, dp)
         model%a_t(:, t) = [0.2_dp * real(t, dp), -0.1_dp * real(t, dp)]
         model%r_t(:, :, t) = model%r
         model%r_t(1, 1, t) = 0.3_dp + 0.1_dp * real(t, dp)
         model%g_t(:, 1, t) = [real(t, dp), -0.5_dp * real(t, dp)]
         model%q_noise_t(1, 1, t) = 0.1_dp * real(t, dp)
      end do
      call assert_true('KFAS edge time-varying model valid', marss_model_valid(model), failures)
      call marss_to_kfas_model(model, kmodel, info)
      call assert_true('KFAS edge time-varying conversion status', info == 0, failures)
      if (info == 0) then
         call assert_true('KFAS edge time-varying flags', all(kmodel%time_varying == 1), failures)
         call assert_true('KFAS edge time-varying slices', size(kmodel%tmat, 3) == 3, failures)
         call assert_close('KFAS edge B index shift', kmodel%tmat(1, 1, 1), &
            model%b_t(1, 1, 2), 1.0e-12_dp, failures)
         call assert_close('KFAS edge U index shift', kmodel%tmat(1, 3, 2), &
            model%u_t(1, 3), 1.0e-12_dp, failures)
         call assert_close('KFAS edge terminal T zero', maxval(abs(kmodel%tmat(:, :, 3))), &
            0.0_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge G index shift', kmodel%rmat(1, 1, 1), &
            model%g_t(1, 1, 2), 1.0e-12_dp, failures)
         call assert_close('KFAS edge terminal G zero', maxval(abs(kmodel%rmat(:, :, 3))), &
            0.0_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge Q index shift', kmodel%q(1, 1, 1), &
            model%q_noise_t(1, 1, 2), 1.0e-12_dp, failures)
         call assert_close('KFAS edge terminal Q zero', maxval(abs(kmodel%q(:, :, 3))), &
            0.0_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge Z current time', kmodel%z(2, 2, 3), &
            model%z_t(2, 2, 3), 1.0e-12_dp, failures)
         call assert_close('KFAS edge R current time', kmodel%h(1, 1, 2), &
            model%r_t(1, 1, 2), 1.0e-12_dp, failures)
      end if

      allocate(initial_model%y(1, 2), initial_model%b(1, 1), initial_model%u(1))
      allocate(initial_model%q(1, 1), initial_model%z(1, 1), initial_model%a(1))
      allocate(initial_model%r(1, 1), initial_model%x0(1), initial_model%v0(1, 1))
      initial_model%y = 0.0_dp
      initial_model%b = 0.8_dp
      initial_model%u = 0.2_dp
      initial_model%q = 0.3_dp
      initial_model%z = 1.0_dp
      initial_model%a = 0.0_dp
      initial_model%r = 0.4_dp
      initial_model%x0 = 0.5_dp
      initial_model%v0 = 0.7_dp
      initial_model%tinitx = 0
      initial_model%diffuse = .false.
      call marss_to_kfas_model(initial_model, kmodel, info)
      call assert_true('KFAS edge tinitx0 conversion status', info == 0, failures)
      if (info == 0) then
         call assert_close('KFAS edge tinitx0 a1', kmodel%a1(1), 0.6_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge tinitx0 P1 x1', kmodel%p1(1, 1), 0.748_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge tinitx0 P1 x0', kmodel%p1(3, 3), 0.7_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge tinitx0 cross covariance', kmodel%p1(1, 3), &
            0.56_dp, 1.0e-12_dp, failures)
      end if


      allocate(loading_model%y(2, 2), loading_model%b(1, 1), loading_model%u(1))
      allocate(loading_model%q(1, 1), loading_model%z(2, 1), loading_model%a(2))
      allocate(loading_model%r(2, 2), loading_model%x0(1), loading_model%v0(1, 1))
      allocate(loading_model%h(2, 1), loading_model%r_noise(1, 1))
      allocate(loading_model%l(1, 1), loading_model%v0_noise(1, 1))
      loading_model%y = 0.0_dp
      loading_model%b = 0.6_dp
      loading_model%u = 0.0_dp
      loading_model%q = 0.2_dp
      loading_model%z(:, 1) = [1.0_dp, 0.5_dp]
      loading_model%a = 0.0_dp
      loading_model%r = 0.0_dp
      loading_model%x0 = 0.0_dp
      loading_model%v0 = 0.0_dp
      loading_model%h(:, 1) = [1.0_dp, 2.0_dp]
      loading_model%r_noise(1, 1) = 0.3_dp
      loading_model%l(1, 1) = 2.0_dp
      loading_model%v0_noise(1, 1) = 0.4_dp
      loading_model%tinitx = 1
      loading_model%diffuse = .false.
      call marss_to_kfas_model(loading_model, kmodel, info)
      call assert_true('KFAS edge H/L loading conversion status', info == 0, failures)
      if (info == 0) then
         call assert_close('KFAS edge H loading covariance 11', kmodel%h(1, 1, 1), 0.3_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge H loading covariance 12', kmodel%h(1, 2, 1), 0.6_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge H loading covariance 22', kmodel%h(2, 2, 1), 1.2_dp, 1.0e-12_dp, failures)
         call assert_close('KFAS edge L loading initial covariance', kmodel%p1(1, 1), 1.6_dp, 1.0e-12_dp, failures)
      end if
   end subroutine test_kfas_edge_conversion


   subroutine test_kfas_backend_parity(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: initial_model
      type(marss_model) :: diffuse_model
      type(marss_kf_result) :: native
      type(marss_kf_result) :: kresult
      real(dp) :: loglik_only
      real(dp) :: nan
      integer :: info
      integer :: t

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(model%y(2, 5), model%b(2, 2), model%u(2), model%q(2, 2))
      allocate(model%z(2, 2), model%a(2), model%r(2, 2), model%x0(2), model%v0(2, 2))
      allocate(model%g(2, 1), model%q_noise(1, 1))
      allocate(model%b_t(2, 2, 5), model%u_t(2, 5), model%z_t(2, 2, 5))
      allocate(model%a_t(2, 5), model%r_t(2, 2, 5))
      allocate(model%g_t(2, 1, 5), model%q_noise_t(1, 1, 5))
      model%y = reshape([ &
         0.30_dp, -0.10_dp, &
         nan, nan, &
         0.52_dp, 0.08_dp, &
         0.41_dp, nan, &
         0.63_dp, 0.21_dp], [2, 5])
      model%b = reshape([0.72_dp, 0.05_dp, -0.08_dp, 0.61_dp], [2, 2])
      model%u = [0.04_dp, -0.02_dp]
      model%q = 0.0_dp
      model%z = reshape([1.00_dp, 0.35_dp, -0.20_dp, 0.85_dp], [2, 2])
      model%a = [0.08_dp, -0.03_dp]
      model%r = reshape([0.40_dp, 0.12_dp, 0.12_dp, 0.32_dp], [2, 2])
      model%x0 = [0.10_dp, -0.15_dp]
      model%v0 = reshape([0.70_dp, 0.10_dp, 0.10_dp, 0.45_dp], [2, 2])
      model%g(:, 1) = [1.0_dp, 0.45_dp]
      model%q_noise(1, 1) = 0.18_dp
      do t = 1, 5
         model%b_t(:, :, t) = model%b
         model%b_t(1, 1, t) = 0.68_dp + 0.015_dp * real(t, dp)
         model%b_t(2, 2, t) = 0.58_dp + 0.010_dp * real(t, dp)
         model%u_t(:, t) = [0.015_dp * real(t, dp), -0.008_dp * real(t, dp)]
         model%z_t(:, :, t) = model%z
         model%z_t(1, 2, t) = -0.12_dp - 0.01_dp * real(t, dp)
         model%a_t(:, t) = [0.02_dp * real(t, dp), -0.01_dp * real(t, dp)]
         model%r_t(:, :, t) = model%r
         model%r_t(1, 1, t) = 0.34_dp + 0.015_dp * real(t, dp)
         model%r_t(2, 2, t) = 0.27_dp + 0.010_dp * real(t, dp)
         model%g_t(:, 1, t) = [1.0_dp, 0.30_dp + 0.03_dp * real(t, dp)]
         model%q_noise_t(1, 1, t) = 0.12_dp + 0.015_dp * real(t, dp)
      end do
      model%tinitx = 1
      model%diffuse = .false.
      call marss_kfss(model, native)
      call marss_kfas(model, kresult)
      call marss_kfas_loglik(model, loglik_only, info)
      call assert_true('KFAS runtime logLik-only status', info == 0, failures)
      if (info == 0 .and. native%ok) then
         call assert_close('KFAS runtime logLik-only value', loglik_only, native%loglik, 2.0e-7_dp, failures)
      end if
      call marss_kfas_loglik(model, loglik_only, info, return_lag_one=.true.)
      call assert_true('KFAS runtime stacked logLik-only status', info == 0, failures)
      if (info == 0 .and. native%ok) then
         call assert_close('KFAS runtime stacked logLik-only value', loglik_only, native%loglik, 2.0e-7_dp, failures)
      end if
      call assert_true('KFAS runtime non-diffuse native status', native%ok, failures)
      call assert_true('KFAS runtime non-diffuse status', kresult%ok, failures)
      if (native%ok .and. kresult%ok) then
         call assert_close('KFAS runtime non-diffuse logLik', kresult%loglik, native%loglik, 2.0e-7_dp, failures)
         call assert_close('KFAS runtime non-diffuse xpred', maxval(abs(kresult%x_pred - native%x_pred)), &
            0.0_dp, 2.0e-7_dp, failures)
         call assert_close('KFAS runtime non-diffuse xfilt', maxval(abs(kresult%x_filt - native%x_filt)), &
            0.0_dp, 2.0e-7_dp, failures)
         call assert_close('KFAS runtime non-diffuse xsmooth', maxval(abs(kresult%x_smooth - native%x_smooth)), &
            0.0_dp, 5.0e-7_dp, failures)
         call assert_close('KFAS runtime non-diffuse Ppred', maxval(abs(kresult%p_pred - native%p_pred)), &
            0.0_dp, 5.0e-7_dp, failures)
         call assert_close('KFAS runtime non-diffuse Pfilt', maxval(abs(kresult%p_filt - native%p_filt)), &
            0.0_dp, 5.0e-7_dp, failures)
         call assert_close('KFAS runtime non-diffuse Psmooth', maxval(abs(kresult%p_smooth - native%p_smooth)), &
            0.0_dp, 8.0e-7_dp, failures)
         call assert_close('KFAS runtime non-diffuse lag', &
            maxval(abs(kresult%p_lag(:, :, 2:5) - native%p_lag(:, :, 2:5))), 0.0_dp, 1.0e-6_dp, failures)
         call assert_true('KFAS runtime correlated H transformed', kresult%observation_covariance_transformed, failures)
         call assert_close('KFAS runtime finite Pinf', maxval(abs(kresult%p_pred_inf)), 0.0_dp, 1.0e-14_dp, failures)
         call assert_close('KFAS runtime finite Finf', maxval(abs(kresult%sigma_inf)), 0.0_dp, 1.0e-14_dp, failures)
         call assert_true('KFAS runtime lag available', kresult%lag_one_available, failures)
      end if
      call marss_kfas(model, kresult, return_lag_one=.false.)
      call assert_true('KFAS runtime unstacked status', kresult%ok, failures)
      if (kresult%ok) then
         call assert_close('KFAS runtime unstacked xsmooth', maxval(abs(kresult%x_smooth - native%x_smooth)), &
            0.0_dp, 5.0e-7_dp, failures)
         call assert_true('KFAS runtime unstacked lag unavailable', .not. kresult%lag_one_available, failures)
      end if

      allocate(initial_model%y(1, 4), initial_model%b(1, 1), initial_model%u(1))
      allocate(initial_model%q(1, 1), initial_model%z(1, 1), initial_model%a(1))
      allocate(initial_model%r(1, 1), initial_model%x0(1), initial_model%v0(1, 1))
      initial_model%y(1, :) = [0.4_dp, 0.1_dp, 0.7_dp, 0.5_dp]
      initial_model%b = 0.75_dp
      initial_model%u = 0.12_dp
      initial_model%q = 0.20_dp
      initial_model%z = 1.1_dp
      initial_model%a = -0.05_dp
      initial_model%r = 0.30_dp
      initial_model%x0 = 0.25_dp
      initial_model%v0 = 0.65_dp
      initial_model%tinitx = 0
      initial_model%diffuse = .false.
      call marss_kfss(initial_model, native)
      call marss_kfas(initial_model, kresult)
      call assert_true('KFAS runtime tinitx0 status', native%ok .and. kresult%ok, failures)
      if (native%ok .and. kresult%ok) then
         call assert_close('KFAS runtime tinitx0 x0T', kresult%x0_smooth(1), &
            native%x0_smooth(1), 2.0e-7_dp, failures)
         call assert_close('KFAS runtime tinitx0 V0T', kresult%v0_smooth(1, 1), &
            native%v0_smooth(1, 1), 2.0e-7_dp, failures)
         call assert_close('KFAS runtime tinitx0 lag1', kresult%p_lag(1, 1, 1), &
            native%p_lag(1, 1, 1), 5.0e-7_dp, failures)
      end if

      allocate(diffuse_model%y(2, 5), diffuse_model%b(2, 2), diffuse_model%u(2), diffuse_model%q(2, 2))
      allocate(diffuse_model%z(2, 2), diffuse_model%a(2), diffuse_model%r(2, 2))
      allocate(diffuse_model%x0(2), diffuse_model%v0(2, 2))
      allocate(diffuse_model%g(2, 1), diffuse_model%q_noise(1, 1))
      diffuse_model%y = reshape([ &
         nan, 0.2_dp, &
         0.1_dp, -0.1_dp, &
         0.4_dp, 0.3_dp, &
         0.2_dp, nan, &
         0.5_dp, 0.15_dp], [2, 5])
      diffuse_model%b = reshape([0.8_dp, 0.0_dp, 0.1_dp, 0.65_dp], [2, 2])
      diffuse_model%u = [0.0_dp, 0.02_dp]
      diffuse_model%q = 0.0_dp
      diffuse_model%z = reshape([1.0_dp, 0.4_dp, -0.25_dp, 0.9_dp], [2, 2])
      diffuse_model%a = 0.0_dp
      diffuse_model%r = reshape([0.35_dp, 0.08_dp, 0.08_dp, 0.28_dp], [2, 2])
      diffuse_model%x0 = 0.0_dp
      diffuse_model%v0 = 0.0_dp
      diffuse_model%v0(1, 1) = 1.0_dp
      diffuse_model%g(:, 1) = [1.0_dp, 0.5_dp]
      diffuse_model%q_noise(1, 1) = 0.10_dp
      diffuse_model%tinitx = 1
      diffuse_model%diffuse = .true.
      call marss_kfas(diffuse_model, kresult)
      call marss_kfas_loglik(diffuse_model, loglik_only, info)
      call assert_true('KFAS runtime diffuse logLik-only status', info == 0, failures)
      if (info == 0 .and. kresult%ok) then
         call assert_close('KFAS runtime diffuse logLik-only value', loglik_only, kresult%loglik, 2.0e-7_dp, failures)
      end if
      call assert_true('KFAS runtime diffuse status', kresult%ok, failures)
      if (kresult%ok) then
         call assert_true('KFAS runtime diffuse Pinf starts nonzero', &
            maxval(abs(kresult%p_pred_inf(:, :, 1))) > 0.0_dp, failures)
         call assert_true('KFAS runtime diffuse Finf observed', maxval(abs(kresult%sigma_inf)) > 0.0_dp, failures)
         call assert_true('KFAS runtime diffuse rank resolved', kresult%remaining_diffuse_rank == 0, failures)
         call assert_true('KFAS runtime diffuse end reported', kresult%diffuse_end >= 1, failures)
         call assert_true('KFAS runtime diffuse lag1 NA', all(ieee_is_nan(kresult%p_lag(:, :, 1))), failures)
         call assert_true('KFAS runtime diffuse correlated H transformed', &
            kresult%observation_covariance_transformed, failures)
         call assert_true('KFAS runtime diffuse gain exposed', maxval(abs(kresult%gain_inf)) > 0.0_dp, failures)
      end if

      deallocate(initial_model%y, initial_model%b, initial_model%u, initial_model%q)
      deallocate(initial_model%z, initial_model%a, initial_model%r, initial_model%x0, initial_model%v0)
      allocate(initial_model%y(2, 4), initial_model%b(1, 1), initial_model%u(1))
      allocate(initial_model%q(1, 1), initial_model%z(2, 1), initial_model%a(2))
      allocate(initial_model%r(2, 2), initial_model%x0(1), initial_model%v0(1, 1))
      initial_model%y = reshape([0.3_dp, 0.2_dp, 0.5_dp, 0.4_dp, 0.1_dp, -0.1_dp, 0.6_dp, 0.5_dp], [2, 4])
      initial_model%b = 0.7_dp
      initial_model%u = 0.05_dp
      initial_model%q = 0.16_dp
      initial_model%z(:, 1) = [1.0_dp, 0.7_dp]
      initial_model%a = [0.0_dp, 0.02_dp]
      initial_model%r = 0.0_dp
      initial_model%r(2, 2) = 0.25_dp
      initial_model%x0 = 0.2_dp
      initial_model%v0 = 0.5_dp
      initial_model%tinitx = 1
      initial_model%diffuse = .false.
      call marss_kfss(initial_model, native)
      call marss_kfas(initial_model, kresult)
      call assert_true('KFAS runtime singular R status', native%ok .and. kresult%ok, failures)
      if (native%ok .and. kresult%ok) then
         call assert_close('KFAS runtime singular R logLik', kresult%loglik, native%loglik, 2.0e-7_dp, failures)
         call assert_close('KFAS runtime singular R state', maxval(abs(kresult%x_smooth - native%x_smooth)), &
            0.0_dp, 5.0e-7_dp, failures)
         call assert_close('KFAS runtime singular R covariance', maxval(abs(kresult%p_smooth - native%p_smooth)), &
            0.0_dp, 5.0e-7_dp, failures)
         call assert_close('KFAS runtime deterministic filtered variance', maxval(abs(kresult%p_filt)), &
            0.0_dp, 2.0e-14_dp, failures)
      end if
   end subroutine test_kfas_backend_parity

   subroutine test_time_varying_deterministic(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: restored
      type(marss_kf_result) :: kf
      type(marss_simulation) :: simulation
      real(dp), allocatable :: theta(:)
      real(dp) :: states(4)
      integer :: info

      allocate(model%y(1, 4), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
      allocate(model%b_t(1, 1, 4), model%u_t(1, 4), model%q_t(1, 1, 4))
      allocate(model%z_t(1, 1, 4), model%a_t(1, 4), model%r_t(1, 1, 4))
      model%b = 1.0_dp
      model%u = 0.0_dp
      model%q = 0.0_dp
      model%z = 1.0_dp
      model%a = 0.0_dp
      model%r = 0.0_dp
      model%x0 = 1.0_dp
      model%v0 = 0.0_dp
      model%tinitx = 1
      model%b_t(1, 1, :) = [1.0_dp, 2.0_dp, 0.5_dp, 1.5_dp]
      model%u_t(1, :) = [0.0_dp, 0.5_dp, -0.5_dp, 1.0_dp]
      model%q_t = 0.0_dp
      model%z_t(1, 1, :) = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp]
      model%a_t(1, :) = [0.0_dp, 1.0_dp, 0.25_dp, -0.5_dp]
      model%r_t = 0.0_dp
      states = [1.0_dp, 2.5_dp, 0.75_dp, 2.125_dp]
      model%y(1, :) = [1.0_dp, 6.0_dp, -0.5_dp, 0.5625_dp]

      call assert_true('time-varying model validation', marss_model_valid(model), failures)
      call marss_kfss(model, kf)
      call assert_true('time-varying deterministic KF status', kf%ok .and. kf%info == 0, failures)
      if (kf%ok) then
         call assert_vector_close('time-varying deterministic states', kf%x_smooth(1, :), states, 1.0e-12_dp, failures)
         call assert_close('time-varying deterministic logLik', kf%loglik, 0.0_dp, 1.0e-12_dp, failures)
      end if
      call marss_simulate(model, 4, 1, 991_i8, simulation, info)
      call assert_true('time-varying deterministic simulation status', info == 0, failures)
      if (info == 0) then
         call assert_vector_close('time-varying simulated states', simulation%states(1, :, 1), states, 1.0e-12_dp, failures)
         call assert_vector_close('time-varying simulated data', simulation%data(1, :, 1), model%y(1, :), 1.0e-12_dp, failures)
      end if
      call marss_vectorizeparam(model, theta)
      call assert_true('time-varying vectorization count', size(theta) == 26, failures)
      call marss_unvectorizeparam(model, theta, restored, info)
      call assert_true('time-varying vectorization status', info == 0, failures)
      if (info == 0) then
         call assert_true('time-varying B roundtrip', maxval(abs(restored%b_t - model%b_t)) <= 1.0e-12_dp, failures)
         call assert_true('time-varying Q roundtrip', maxval(abs(restored%q_t - model%q_t)) <= 1.0e-12_dp, failures)
         call assert_true('time-varying Z roundtrip', maxval(abs(restored%z_t - model%z_t)) <= 1.0e-12_dp, failures)
         call assert_true('time-varying R roundtrip', maxval(abs(restored%r_t - model%r_t)) <= 1.0e-12_dp, failures)
      end if
   end subroutine test_time_varying_deterministic


   subroutine test_marxss_numerical_form(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: direct
      type(marss_model) :: extended
      type(marss_kf_result) :: kf_direct
      type(marss_kf_result) :: kf_extended
      type(marss_simulation) :: sim_direct
      type(marss_simulation) :: sim_extended
      integer :: info_direct
      integer :: info_extended

      allocate(direct%y(1, 4), direct%b(1, 1), direct%u(1), direct%q(1, 1))
      allocate(direct%z(1, 1), direct%a(1), direct%r(1, 1), direct%x0(1), direct%v0(1, 1))
      allocate(direct%u_t(1, 4), direct%a_t(1, 4))
      direct%y(1, :) = [0.4_dp, 0.8_dp, 0.7_dp, 1.1_dp]
      direct%b = 0.6_dp
      direct%u = 0.1_dp
      direct%q = 0.25_dp
      direct%z = 1.1_dp
      direct%a = 0.2_dp
      direct%r = 0.36_dp
      direct%x0 = -0.3_dp
      direct%v0 = 0.49_dp
      direct%u_t(1, :) = 0.1_dp + 0.5_dp * [1.0_dp, -1.0_dp, 0.5_dp, 2.0_dp]
      direct%a_t(1, :) = 0.2_dp - 0.2_dp * [0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp]
      direct%tinitx = 0

      extended = direct
      deallocate(extended%u_t, extended%a_t)
      allocate(extended%state_covariates(1, 4), extended%c_coef(1, 1))
      allocate(extended%obs_covariates(1, 4), extended%d_coef(1, 1))
      allocate(extended%g(1, 1), extended%q_noise(1, 1))
      allocate(extended%h(1, 1), extended%r_noise(1, 1))
      allocate(extended%l(1, 1), extended%v0_noise(1, 1))
      extended%state_covariates(1, :) = [1.0_dp, -1.0_dp, 0.5_dp, 2.0_dp]
      extended%c_coef = 0.5_dp
      extended%obs_covariates(1, :) = [0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp]
      extended%d_coef = -0.2_dp
      extended%g = 2.0_dp
      extended%q_noise = 0.0625_dp
      extended%h = 3.0_dp
      extended%r_noise = 0.04_dp
      extended%l = 0.7_dp
      extended%v0_noise = 1.0_dp

      call assert_true('marxss numerical model valid', marss_model_valid(extended), failures)
      call marss_kfss(direct, kf_direct)
      call marss_kfss(extended, kf_extended)
      call assert_true('marxss numerical filter status', kf_direct%ok .and. kf_extended%ok, failures)
      if (kf_direct%ok .and. kf_extended%ok) then
         call assert_close('marxss numerical logLik', kf_extended%loglik, kf_direct%loglik, 1.0e-11_dp, failures)
         call assert_matrix_close('marxss numerical predictions', kf_extended%x_pred, kf_direct%x_pred, 1.0e-11_dp, failures)
      end if
      call marss_simulate(direct, 4, 1, 4242_i8, sim_direct, info_direct)
      call marss_simulate(extended, 4, 1, 4242_i8, sim_extended, info_extended)
      call assert_true('marxss numerical simulation status', info_direct == 0 .and. info_extended == 0, failures)
      if (info_direct == 0 .and. info_extended == 0) then
         call assert_matrix_close('marxss numerical simulated states', sim_extended%states(:, :, 1), &
            sim_direct%states(:, :, 1), 1.0e-11_dp, failures)
         call assert_matrix_close('marxss numerical simulated data', sim_extended%data(:, :, 1), &
            sim_direct%data(:, :, 1), 1.0e-11_dp, failures)
      end if
   end subroutine test_marxss_numerical_form


   subroutine test_rich_constraints(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: rebuilt
      type(marss_constraints) :: constraints
      real(dp), allocatable :: theta(:)
      integer :: info

      allocate(model%y(1, 3), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
      allocate(model%state_covariates(1, 3), model%c_coef(1, 1))
      allocate(model%obs_covariates(1, 3), model%d_coef(1, 1))
      allocate(model%g(1, 1), model%q_noise(1, 1))
      allocate(model%h(1, 1), model%r_noise(1, 1))
      allocate(model%l(1, 1), model%v0_noise(1, 1))
      model%y(1, :) = [0.2_dp, 0.4_dp, 0.7_dp]
      model%b = 0.7_dp
      model%u = 0.1_dp
      model%q = 0.2_dp
      model%z = 1.0_dp
      model%a = -0.1_dp
      model%r = 0.3_dp
      model%x0 = 0.0_dp
      model%v0 = 0.4_dp
      model%state_covariates(1, :) = [1.0_dp, -0.5_dp, 0.25_dp]
      model%c_coef = 0.15_dp
      model%obs_covariates(1, :) = [0.5_dp, 1.0_dp, -1.0_dp]
      model%d_coef = -0.2_dp
      model%g = 1.2_dp
      model%q_noise = 0.2_dp / (1.2_dp**2)
      model%h = 0.8_dp
      model%r_noise = 0.3_dp / (0.8_dp**2)
      model%l = 0.5_dp
      model%v0_noise = 0.4_dp / (0.5_dp**2)
      model%tinitx = 0

      call make_scalar_constraint(constraints%q, model%q_noise(1, 1))
      call make_scalar_constraint(constraints%r, model%r_noise(1, 1))
      call make_scalar_constraint(constraints%v0, model%v0_noise(1, 1))
      call make_scalar_constraint(constraints%c, model%c_coef(1, 1))
      call make_scalar_constraint(constraints%d, model%d_coef(1, 1))
      call make_scalar_constraint(constraints%g, model%g(1, 1))
      call make_scalar_constraint(constraints%h, model%h(1, 1))
      call make_scalar_constraint(constraints%l, model%l(1, 1))
      call assert_true('rich constraints valid', marss_constraints_valid(model, constraints), failures)
      call assert_true('rich constraints parameter count', marss_constraints_parameter_count(constraints) == 8, failures)
      call marss_constraint_start_vector(constraints, theta)
      call marss_apply_constraints(model, constraints, theta, rebuilt, info)
      call assert_true('rich constraints apply status', info == 0, failures)
      if (info == 0) then
         call assert_close('rich constraint C', rebuilt%c_coef(1, 1), model%c_coef(1, 1), 1.0e-12_dp, failures)
         call assert_close('rich constraint D', rebuilt%d_coef(1, 1), model%d_coef(1, 1), 1.0e-12_dp, failures)
         call assert_close('rich constraint G', rebuilt%g(1, 1), model%g(1, 1), 1.0e-12_dp, failures)
         call assert_close('rich constraint H', rebuilt%h(1, 1), model%h(1, 1), 1.0e-12_dp, failures)
         call assert_close('rich constraint L', rebuilt%l(1, 1), model%l(1, 1), 1.0e-12_dp, failures)
         call assert_close('rich constraint Q noise', rebuilt%q_noise(1, 1), model%q_noise(1, 1), 1.0e-12_dp, failures)
         call assert_close('rich constraint R noise', rebuilt%r_noise(1, 1), model%r_noise(1, 1), 1.0e-12_dp, failures)
         call assert_close('rich constraint V0 noise', rebuilt%v0_noise(1, 1), model%v0_noise(1, 1), 1.0e-12_dp, failures)
      end if
   end subroutine test_rich_constraints

   subroutine test_constraint_labels(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_constraint_block) :: block
      type(marss_constraint_block) :: affine_block
      type(marss_constraint_block) :: invalid_block
      type(marss_constraints) :: constraints
      type(marss_constraints) :: updated
      character(len=64), allocatable :: names(:)
      character(len=64), allocatable :: vector_names(:)
      character(len=16) :: custom_names(2)
      character(len=32) :: input_names(2)
      character(len=32) :: one_name(1)
      real(dp), allocatable :: theta(:)
      real(dp) :: affine_design(3, 2)
      real(dp) :: affine_fixed(3)
      real(dp) :: affine_start(2)
      real(dp) :: entry_numeric(3)
      real(dp) :: entry_starts(2)
      real(dp) :: input_values(2)
      real(dp) :: values(4)
      character(len=64) :: entry_expressions(3)
      logical :: entry_is_expression(3)
      integer :: info
      integer :: labels(4)

      values = [0.2_dp, 5.0_dp, 0.4_dp, 0.6_dp]
      labels = [1, 0, 2, 1]
      custom_names = [character(len=16) :: 'shared-loading', 'second-loading']
      call marss_constraint_from_labels(values, labels, block, info, custom_names)
      call assert_true('labeled constraint status', info == 0, failures)
      if (info == 0) then
         call assert_close('labeled constraint fixed value', block%fixed(2), 5.0_dp, 0.0_dp, failures)
         call assert_close('labeled constraint shared start', block%start(1), 0.4_dp, 1.0e-15_dp, failures)
         call assert_close('labeled constraint second start', block%start(2), 0.4_dp, 1.0e-15_dp, failures)
         call assert_true('labeled constraint shared design', &
            abs(block%design(1, 1) - 1.0_dp) <= 1.0e-15_dp .and. &
            abs(block%design(4, 1) - 1.0_dp) <= 1.0e-15_dp, failures)
         constraints%u = block
         call marss_free_parameter_names(constraints, names)
         call assert_true('labeled constraint names preserved', size(names) == 2 .and. &
            trim(names(1)) == 'shared-loading' .and. trim(names(2)) == 'second-loading', failures)
         call marss_vectorized_parameter_names(constraints, vector_names)
         call assert_true('vectorized names carry parameter prefix', size(vector_names) == 2 .and. &
            trim(vector_names(1)) == 'U.shared-loading' .and. &
            trim(vector_names(2)) == 'U.second-loading', failures)
         input_names = [character(len=32) :: 'second-loading', 'shared-loading']
         input_values = [9.0_dp, 8.0_dp]
         call marss_reorder_free_parameters(constraints, input_names, input_values, theta, info)
         call assert_true('named parameter reorder status', info == 0, failures)
         if (info == 0) then
            call assert_vector_close('named parameter reorder values', theta, [8.0_dp, 9.0_dp], 0.0_dp, failures)
         end if
         call marss_constraints_set_start_named(constraints, input_names, input_values, updated, info)
         call assert_true('named constraint start status', info == 0, failures)
         if (info == 0) then
            call assert_vector_close('named constraint start values', updated%u%start, &
               [8.0_dp, 9.0_dp], 0.0_dp, failures)
         end if
         input_names = [character(len=32) :: 'U.second-loading', 'U.shared-loading']
         call marss_reorder_free_parameters(constraints, input_names, input_values, theta, info)
         call assert_true('prefixed parameter reorder status', info == 0, failures)
         if (info == 0) then
            call assert_vector_close('prefixed parameter reorder values', theta, [8.0_dp, 9.0_dp], 0.0_dp, failures)
         end if
      end if

      affine_fixed = [2.0_dp, -1.0_dp, 0.5_dp]
      affine_design = reshape([1.0_dp, 0.0_dp, 2.0_dp, -0.5_dp, 3.0_dp, 0.0_dp], [3, 2])
      affine_start = [4.0_dp, -2.0_dp]
      custom_names = [character(len=16) :: 'trend', 'offset']
      call marss_constraint_from_affine(affine_fixed, affine_design, affine_start, affine_block, info, custom_names)
      call assert_true('affine constraint status', info == 0, failures)
      if (info == 0) then
         call assert_vector_close('affine constraint fixed offsets', affine_block%fixed, affine_fixed, 0.0_dp, failures)
         call assert_matrix_close('affine constraint coefficients', affine_block%design, affine_design, 0.0_dp, failures)
         call assert_vector_close('affine constraint beta starts', affine_block%start, affine_start, 0.0_dp, failures)
         constraints = marss_constraints()
         constraints%a = affine_block
         call marss_free_parameter_names(constraints, names)
         call assert_true('affine constraint names preserved', size(names) == 2 .and. &
            trim(names(1)) == 'trend' .and. trim(names(2)) == 'offset', failures)
      end if

      entry_numeric = [5.0_dp, 0.0_dp, 0.0_dp]
      entry_expressions = [character(len=64) :: '', 'alpha', '2+0.5*alpha+3*beta+1.5*alpha']
      entry_is_expression = [.false., .true., .true.]
      entry_starts = [0.25_dp, -0.5_dp]
      call marss_constraint_from_entries(entry_numeric, entry_expressions, entry_is_expression, &
         affine_block, info, entry_starts)
      call assert_true('list-entry constraint status', info == 0, failures)
      if (info == 0) then
         call assert_vector_close('list-entry fixed offsets', affine_block%fixed, &
            [5.0_dp, 0.0_dp, 2.0_dp], 0.0_dp, failures)
         call assert_matrix_close('list-entry affine coefficients', affine_block%design, &
            reshape([0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [3, 2]), 1.0e-15_dp, failures)
         call assert_vector_close('list-entry starts', affine_block%start, entry_starts, 0.0_dp, failures)
         call assert_true('list-entry parameter names', trim(affine_block%free_names(1)) == 'alpha' .and. &
            trim(affine_block%free_names(2)) == 'beta', failures)
      end if

      entry_numeric = 0.0_dp
      entry_expressions = [character(len=64) :: '1', '', '']
      entry_is_expression = [.true., .false., .false.]
      call marss_constraint_from_entries(entry_numeric, entry_expressions, entry_is_expression, affine_block, info)
      call assert_true('numeric-looking character label remains free', info == 0 .and. &
         size(affine_block%free_names) == 1 .and. trim(affine_block%free_names(1)) == '1', failures)

      entry_expressions = [character(len=64) :: '2+3', '', '']
      call marss_constraint_from_entries(entry_numeric, entry_expressions, entry_is_expression, invalid_block, info)
      call assert_true('list-entry numeric term after plus rejected', info /= 0, failures)

      one_name(1) = 'diag'
      call marss_constraint_from_labels([0.5_dp], [1], block, info, one_name)
      call assert_true('duplicate raw-label setup', info == 0, failures)
      if (info == 0) then
         constraints = marss_constraints()
         constraints%u = block
         constraints%a = block
         input_names = [character(len=32) :: 'diag', 'diag']
         input_values = [7.0_dp, 8.0_dp]
         call marss_reorder_free_parameters(constraints, input_names, input_values, theta, info)
         call assert_true('ambiguous raw parameter labels rejected', info /= 0, failures)
         input_names = [character(len=32) :: 'A.diag', 'U.diag']
         call marss_reorder_free_parameters(constraints, input_names, input_values, theta, info)
         call assert_true('prefixed duplicate-label reorder status', info == 0, failures)
         if (info == 0) then
            call assert_vector_close('prefixed duplicate-label reorder values', theta, &
               [8.0_dp, 7.0_dp], 0.0_dp, failures)
         end if
      end if

      custom_names = [character(len=16) :: 'duplicate', 'duplicate']
      call marss_constraint_from_labels(values, labels, invalid_block, info, custom_names)
      call assert_true('labeled constraint duplicate names rejected', info == 6, failures)
      call marss_constraint_from_affine(affine_fixed, affine_design, affine_start, invalid_block, info, custom_names)
      call assert_true('affine constraint duplicate names rejected', info == 4, failures)
   end subroutine test_constraint_labels

   subroutine test_dfa_model(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_constraints) :: constraints
      type(marss_fit_result) :: fit
      type(marss_constraints) :: fit_constraints
      real(dp) :: data(3, 6)
      real(dp) :: covariates(1, 6)
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: scale(:)
      real(dp) :: row_mean
      real(dp) :: row_sd
      integer :: i
      integer :: info

      data(1, :) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      data(2, :) = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 6.0_dp, 5.0_dp]
      data(3, :) = [-1.0_dp, 0.0_dp, 1.0_dp, 1.5_dp, 2.0_dp, 3.5_dp]
      covariates(1, :) = [-1.0_dp, -0.5_dp, 0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp]
      call marss_dfa_model(data, 2, model, constraints, info, covariates, center=center, scale=scale)
      call assert_true('DFA model status', info == 0, failures)
      if (info /= 0) return
      call assert_true('DFA model valid', marss_model_valid(model), failures)
      call assert_true('DFA constraints valid', marss_constraints_valid(model, constraints), failures)
      call assert_true('DFA parameter count', marss_constraints_parameter_count(constraints) == 9, failures)
      call assert_matrix_close('DFA B identity', model%b, reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
         1.0e-12_dp, failures)
      call assert_matrix_close('DFA Q identity', model%q, reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
         1.0e-12_dp, failures)
      call assert_matrix_close('DFA V0 default', model%v0, reshape([5.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], [2, 2]), &
         1.0e-12_dp, failures)
      call assert_close('DFA triangular fixed loading', model%z(1, 2), 0.0_dp, 1.0e-12_dp, failures)
      call assert_close('DFA first loading anchor', model%z(1, 1), 1.0_dp, 1.0e-12_dp, failures)
      call assert_close('DFA second loading anchor', model%z(2, 2), 1.0_dp, 1.0e-12_dp, failures)
      call assert_close('DFA equal R first', model%r(1, 1), 0.5_dp, 1.0e-12_dp, failures)
      call assert_close('DFA equal R second', model%r(2, 2), 0.5_dp, 1.0e-12_dp, failures)
      call assert_close('DFA zero R off diagonal', model%r(1, 2), 0.0_dp, 1.0e-12_dp, failures)
      call assert_true('DFA covariate coefficient allocated', allocated(model%d_coef), failures)
      do i = 1, 3
         row_mean = sum(model%y(i, :)) / 6.0_dp
         row_sd = sqrt(sum((model%y(i, :) - row_mean)**2) / 5.0_dp)
         call assert_close('DFA transformed row mean', row_mean, 0.0_dp, 1.0e-12_dp, failures)
         call assert_close('DFA transformed row sd', row_sd, 1.0_dp, 1.0e-12_dp, failures)
      end do
      call assert_true('DFA center returned', size(center) == 3, failures)
      call assert_true('DFA scale returned', size(scale) == 3 .and. all(scale > 0.0_dp), failures)

      call marss_dfa_fit(data, 1, 2, 1.0e-5_dp, fit, fit_constraints, info, mstep_iter=6)
      call assert_true('DFA fit status', info == 0, failures)
      if (info == 0) then
         call assert_true('DFA fit finite likelihood', .not. ieee_is_nan(fit%loglik), failures)
         call assert_true('DFA fit retains constraints', marss_constraints_valid(fit%model, fit_constraints), failures)
      end if
   end subroutine test_dfa_model

   subroutine test_dfa_spec(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_dfa_spec) :: spec
      type(marss_model) :: model
      type(marss_constraints) :: constraints
      type(marss_fit_result) :: fit
      real(dp) :: data(3, 6)
      real(dp) :: covariates(1, 6)
      integer :: info

      data = reshape([ &
         1.0_dp, 2.0_dp, 3.0_dp, &
         1.2_dp, 1.9_dp, 3.1_dp, &
         1.4_dp, 2.1_dp, 2.9_dp, &
         1.5_dp, 2.3_dp, 3.2_dp, &
         1.7_dp, 2.2_dp, 3.4_dp, &
         1.8_dp, 2.4_dp, 3.3_dp], [3, 6])
      covariates(1, :) = [-1.0_dp, -0.5_dp, 0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp]
      spec%ntrends = 2
      spec%b = 'diagonal and equal'
      spec%q = 'diagonal and unequal'
      spec%a = 'unequal'
      spec%r = 'equalvarcov'
      spec%x0 = 'unconstrained'
      spec%v0 = 'identity'
      spec%d = 'unconstrained'
      call marss_dfa_build(data, spec, model, constraints, info, covariates)
      call assert_true('DFA spec builder status', info == 0, failures)
      if (info == 0) then
         call assert_true('DFA spec B diagonal equal count', size(constraints%b%design, 2) == 1, failures)
         call assert_true('DFA spec Q diagonal unequal count', size(constraints%q%design, 2) == 2, failures)
         call assert_true('DFA spec A unequal count', size(constraints%a%design, 2) == 3, failures)
         call assert_true('DFA spec R equalvarcov count', size(constraints%r%design, 2) == 2, failures)
         call assert_true('DFA spec x0 unequal count', size(constraints%x0%design, 2) == 2, failures)
         call assert_true('DFA spec D covariate count', size(constraints%d%design, 2) == 3, failures)
         call assert_matrix_close('DFA spec fixed V0 identity', model%v0, &
            reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 0.0_dp, failures)
      end if

      spec = marss_dfa_spec()
      spec%ntrends = 1
      spec%demean = .false.
      spec%z_score = .false.
      call marss_from_data(data(1:2, :), spec, 2, 1.0e-4_dp, fit, constraints, 'bfgs', info)
      call assert_true('DFA high-level BFGS dispatch', info == 0, failures)
      if (info == 0) then
         call assert_true('DFA high-level BFGS beta output', allocated(fit%free_parameters), failures)
      end if

      spec = marss_dfa_spec()
      spec%ntrends = 1
      spec%d = 'unconstrained'
      call marss_from_data(data(1:2, :), spec, 2, 1.0e-4_dp, fit, constraints, 'kem', info, &
         obs_covariates=covariates)
      call assert_true('DFA high-level KEM covariate dispatch', info == 0, failures)
      if (info == 0) then
         call assert_true('DFA high-level KEM D coefficients', allocated(fit%model%d_coef), failures)
      end if
   end subroutine test_dfa_spec

   subroutine make_scalar_constraint(block, start_value)
      type(marss_constraint_block), intent(out) :: block !! One-parameter affine block used by rich-constraint regression tests.
      real(dp), intent(in) :: start_value !! Initial scalar coefficient represented by the one-column design matrix.

      allocate(block%fixed(1), block%design(1, 1), block%start(1))
      block%fixed = 0.0_dp
      block%design = 1.0_dp
      block%start = start_value
   end subroutine make_scalar_constraint

   subroutine test_vectorization(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: restored
      real(dp), allocatable :: theta(:)
      integer :: info

      call make_reference_model(model)
      model%y = 0.0_dp
      call marss_vectorizeparam(model, theta)
      call assert_true('MARSSvectorizeparam count', size(theta) == marss_parameter_count(model), failures)
      call marss_unvectorizeparam(model, theta, restored, info)
      call assert_true('MARSSvectorizeparam roundtrip status', info == 0, failures)
      call assert_matrix_close('vectorize B roundtrip', restored%b, model%b, 1.0e-12_dp, failures)
      call assert_matrix_close('vectorize Q roundtrip', restored%q, model%q, 1.0e-12_dp, failures)
      call assert_matrix_close('vectorize R roundtrip', restored%r, model%r, 1.0e-12_dp, failures)
      call assert_matrix_close('vectorize V0 roundtrip', restored%v0, model%v0, 1.0e-12_dp, failures)
   end subroutine test_vectorization

   subroutine test_free_vectorization(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model_spec) :: spec
      type(marss_model) :: model
      type(marss_model) :: reconstructed
      type(marss_constraints) :: constraints
      real(dp) :: y(2, 4)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: beta_roundtrip(:)
      real(dp), allocatable :: raw(:)
      real(dp), allocatable :: raw_roundtrip(:)
      integer :: info

      y = reshape([1.0_dp, 2.0_dp, 1.1_dp, 2.1_dp, 1.2_dp, 2.2_dp, 1.3_dp, 2.3_dp], [2, 4])
      call marss_build(y, spec, model, constraints, info)
      if (info /= 0) then
         call assert_true('MARSS free vectorization setup', .false., failures)
         return
      end if
      call marss_vectorizeparam(model, beta, constraints, info)
      call assert_true('MARSSvectorizeparam free coordinates', info == 0 .and. &
         size(beta) == marss_parameter_count(model, constraints), failures)
      if (info == 0 .and. size(beta) > 0) beta(1) = beta(1) + 0.125_dp
      call marss_unvectorizeparam(model, beta, reconstructed, info, constraints)
      call assert_true('MARSSvectorizeparam constrained unvectorize', info == 0, failures)
      if (info == 0) then
         call marss_vectorizeparam(reconstructed, beta_roundtrip, constraints, info)
         call assert_true('MARSSvectorizeparam beta round trip status', info == 0, failures)
         if (info == 0) then
            call assert_vector_close('MARSSvectorizeparam beta round trip', beta_roundtrip, beta, &
               1.0e-10_dp, failures)
         end if
      end if

      call marss_vectorizeparam(model, raw)
      call marss_unvectorizeparam(model, raw, reconstructed, info)
      call assert_true('MARSS raw extended vectorization status', info == 0, failures)
      if (info == 0) then
         call marss_vectorizeparam(reconstructed, raw_roundtrip)
         call assert_vector_close('MARSS raw extended vectorization', raw_roundtrip, raw, &
            1.0e-12_dp, failures)
      end if
   end subroutine test_free_vectorization

   subroutine test_simulation(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_simulation) :: sim1
      type(marss_simulation) :: sim2
      logical :: missing2(1, 6)
      logical :: missing3(1, 6, 2)
      integer :: info1
      integer :: info2

      call make_reference_model(model)
      model%y = 0.0_dp
      call marss_simulate(model, 6, 2, 12345_i8, sim1, info1)
      call marss_simulate(model, 6, 2, 12345_i8, sim2, info2)
      call assert_true('MARSSsimulate status', info1 == 0 .and. info2 == 0, failures)
      call assert_true('MARSSsimulate shape', size(sim1%data, 2) == 6 .and. size(sim1%data, 3) == 2, failures)
      call assert_true('MARSSsimulate deterministic data', maxval(abs(sim1%data - sim2%data)) <= 0.0_dp, failures)
      call assert_true('MARSSsimulate deterministic states', maxval(abs(sim1%states - sim2%states)) <= 0.0_dp, failures)

      missing2 = .false.
      missing2(1, 2) = .true.
      call marss_simulate(model, 6, 2, 12345_i8, sim1, info1, missing_mask=missing2)
      call assert_true('MARSSsimulate repeated missing mask status', info1 == 0, failures)
      if (info1 == 0) then
         call assert_true('MARSSsimulate repeated missing mask locations', &
            ieee_is_nan(sim1%data(1, 2, 1)) .and. ieee_is_nan(sim1%data(1, 2, 2)), failures)
         call assert_true('MARSSsimulate repeated missing mask preserves other values', &
            .not. ieee_is_nan(sim1%data(1, 1, 1)), failures)
      end if

      missing3 = .false.
      missing3(1, 3, 2) = .true.
      call marss_simulate(model, 6, 2, 12345_i8, sim2, info2, missing_mask_by_sim=missing3)
      call assert_true('MARSSsimulate replicate-specific missing mask status', info2 == 0, failures)
      if (info2 == 0) then
         call assert_true('MARSSsimulate replicate-specific missing mask locations', &
            .not. ieee_is_nan(sim2%data(1, 3, 1)) .and. ieee_is_nan(sim2%data(1, 3, 2)), failures)
      end if
   end subroutine test_simulation

   subroutine test_analysis_helpers(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      type(marss_residual_result) :: residual
      type(marss_hatyt_result) :: hat
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: fisher(:, :)
      real(dp) :: aic
      real(dp) :: aicc
      integer :: info_h
      integer :: info_f

      call make_reference_model(model)
      model%y(1, :) = [0.9_dp, 0.2_dp, 0.7_dp, 1.1_dp]
      call marss_kfss(model, kf)
      call marss_residuals(model, kf, residual)
      call marss_hatyt(model, kf, hat)
      call assert_true('MARSSresiduals finite', all(.not. ieee_is_nan(residual%innovations)), failures)
      call assert_true('MARSShatyt shape', size(hat%yt, 2) == size(model%y, 2), failures)
      call marss_aic(-10.0_dp, 3, 50, aic, aicc)
      call assert_close('MARSSaic AIC', aic, 26.0_dp, 1.0e-12_dp, failures)
      call assert_close('MARSSaic AICc', aicc, 26.52173913043478_dp, 1.0e-12_dp, failures)
      call marss_hessian(model, h, 1.0e-4_dp, info_h)
      call assert_true('MARSShessian status', info_h == 0, failures)
      if (info_h == 0) then
         call assert_matrix_close('MARSShessian symmetry', h, transpose(h), 2.0e-4_dp, failures)
      end if
      call marss_fisher_i(model, fisher, 1.0e-4_dp, info_f)
      call assert_true('MARSSFisherI status', info_f == 0, failures)
      if (info_h == 0 .and. info_f == 0) then
         call assert_matrix_close('MARSSFisherI relation', fisher, -h, 2.0e-4_dp, failures)
      end if
   end subroutine test_analysis_helpers

   subroutine test_residual_hierarchy(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: correlated
      type(marss_kf_result) :: kf
      type(marss_kf_result) :: kf_correlated
      type(marss_residual_result) :: smoothed
      type(marss_residual_result) :: filtered
      type(marss_residual_result) :: one_step
      real(dp) :: nan_value

      call make_reference_model(model)
      model%y(1, :) = [0.9_dp, 0.2_dp, 0.7_dp, 1.1_dp]
      call marss_kfss(model, kf)
      call marss_residuals_smoothed(model, kf, smoothed)
      call marss_residuals_filtered(model, kf, filtered)
      call marss_residuals_one_step(model, kf, one_step)
      call assert_true('MARSSresiduals tT status', smoothed%info == 0, failures)
      call assert_true('MARSSresiduals tt status', filtered%info == 0, failures)
      call assert_true('MARSSresiduals tt1 status', one_step%info == 0, failures)
      if (smoothed%info == 0) then
         call assert_true('MARSSresiduals tT final state undefined', ieee_is_nan(smoothed%state_residuals(1, 4)), failures)
         call assert_true('MARSSresiduals tT model finite', all(.not. ieee_is_nan(smoothed%model_residuals)), failures)
      end if
      if (filtered%info == 0) then
         call assert_true('MARSSresiduals tt state undefined', all(ieee_is_nan(filtered%state_residuals)), failures)
      end if
      if (one_step%info == 0) then
         call assert_vector_close('MARSSresiduals tt1 innovations', one_step%model_residuals(1, :), &
            kf%innov(1, :), 2.0e-12_dp, failures)
         call assert_close('MARSSresiduals tt1 variance', one_step%var_residuals(1, 1, 2), &
            kf%sigma(1, 1, 2), 2.0e-12_dp, failures)
      end if

      allocate(correlated%y(2, 1), correlated%b(1, 1), correlated%u(1), correlated%q(1, 1))
      allocate(correlated%z(2, 1), correlated%a(2), correlated%r(2, 2))
      allocate(correlated%x0(1), correlated%v0(1, 1))
      correlated%y(:, 1) = [2.0_dp, 1.0_dp]
      correlated%b = 0.0_dp
      correlated%u = 0.0_dp
      correlated%q = 0.0_dp
      correlated%z = 0.0_dp
      correlated%a = 0.0_dp
      correlated%r = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
      correlated%x0 = 0.0_dp
      correlated%v0 = 0.0_dp
      correlated%tinitx = 1
      call marss_kfss(correlated, kf_correlated)
      call marss_residuals_one_step(correlated, kf_correlated, one_step)
      call assert_true('MARSSresiduals correlated Cholesky status', one_step%info == 0, failures)
      if (one_step%info == 0) then
         call assert_vector_close('MARSSresiduals ordered Cholesky', one_step%std_residuals(1:2, 1), &
            [1.0_dp, 0.0_dp], 2.0e-12_dp, failures)
      end if

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      model%y(1, 3) = nan_value
      call marss_kfss(model, kf)
      call marss_residuals_smoothed(model, kf, smoothed, normalize=.true.)
      call assert_true('MARSSresiduals missing normalized status', smoothed%info == 0, failures)
      if (smoothed%info == 0) then
         call assert_true('MARSSresiduals missing response masked', ieee_is_nan(smoothed%model_residuals(1, 3)), failures)
         call assert_true('MARSSresiduals missing expectation finite', &
            .not. ieee_is_nan(smoothed%expected_observed_residuals(1, 3)), failures)
      end if
   end subroutine test_residual_hierarchy

   subroutine test_harvey_residuals(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      type(marss_residual_result) :: harvey
      type(marss_residual_result) :: normalized
      real(dp) :: expected_var(2, 2)

      call make_reference_model(model)
      model%y(1, :) = [0.9_dp, 0.2_dp, 0.5_dp, 1.1_dp]
      call marss_kfss(model, kf)
      call marss_residuals_harvey(model, kf, harvey)
      call assert_true('MARSSresiduals Harvey status', harvey%info == 0, failures)
      if (harvey%info == 0) then
         call assert_vector_close('MARSSresiduals Harvey observation', harvey%model_residuals(1, :), &
            [0.227518008801703_dp, -0.277147439765325_dp, -0.124819304003807_dp, &
             0.182761806152382_dp], 3.0e-12_dp, failures)
         call assert_vector_close('MARSSresiduals Harvey state', harvey%state_residuals(1, 1:3), &
            [-0.234031794327761_dp, 0.0192511268262891_dp, 0.164485625537144_dp], &
            3.0e-12_dp, failures)
         call assert_true('MARSSresiduals Harvey final state undefined', &
            ieee_is_nan(harvey%state_residuals(1, 4)), failures)
         expected_var = reshape([0.167784950657329_dp, -0.0917430923242745_dp, &
                                 -0.0917430923242745_dp, 0.141539476660712_dp], [2, 2])
         call assert_matrix_close('MARSSresiduals Harvey variance', harvey%var_residuals(:, :, 1), &
            expected_var, 3.0e-12_dp, failures)
         call assert_true('MARSSresiduals Harvey final joint standardized undefined', &
            all(ieee_is_nan(harvey%std_residuals(:, 4))), failures)
         call assert_true('MARSSresiduals Harvey final marginal observation finite', &
            .not. ieee_is_nan(harvey%marginal_residuals(1, 4)), failures)
         call assert_true('MARSSresiduals Harvey final block observation finite', &
            .not. ieee_is_nan(harvey%block_cholesky_residuals(1, 4)), failures)
      end if

      call marss_residuals_harvey(model, kf, normalized, normalize=.true.)
      call assert_true('MARSSresiduals Harvey normalized status', normalized%info == 0, failures)
      if (normalized%info == 0) then
         call assert_close('MARSSresiduals Harvey normalized observation', normalized%model_residuals(1, 1), &
            0.227518008801703_dp / sqrt(0.4_dp), 3.0e-12_dp, failures)
         call assert_close('MARSSresiduals Harvey normalized state', normalized%state_residuals(1, 1), &
            -0.234031794327761_dp / sqrt(0.3_dp), 3.0e-12_dp, failures)
      end if
   end subroutine test_harvey_residuals

   subroutine test_hatyt_correlated_missing(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      type(marss_hatyt_result) :: hat
      real(dp) :: nan

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(model%y(2, 1), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(2, 1), model%a(2), model%r(2, 2), model%x0(1), model%v0(1, 1))
      model%y(:, 1) = [2.0_dp, nan]
      model%b = 0.0_dp
      model%u = 0.0_dp
      model%q = 0.0_dp
      model%z = 0.0_dp
      model%a = 0.0_dp
      model%r = reshape([1.0_dp, 0.5_dp, 0.5_dp, 2.0_dp], [2, 2])
      model%x0 = 0.0_dp
      model%v0 = 0.0_dp
      model%tinitx = 1
      allocate(kf%x_smooth(1, 1), kf%p_smooth(1, 1, 1))
      kf%x_smooth = 0.0_dp
      kf%p_smooth = 0.0_dp

      call marss_hatyt(model, kf, hat)
      call assert_vector_close('MARSShatyt correlated conditional mean', hat%yt(:, 1), &
         [2.0_dp, 1.0_dp], 1.0e-12_dp, failures)
      call assert_close('MARSShatyt correlated conditional second moment', &
         hat%ot(2, 2, 1), 2.75_dp, 1.0e-12_dp, failures)
      call assert_close('MARSShatyt observed-missing cross moment', &
         hat%ot(1, 2, 1), 2.0_dp, 1.0e-12_dp, failures)
   end subroutine test_hatyt_correlated_missing


   subroutine test_hatyt_full(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      type(marss_hatyt_result) :: hat
      real(dp) :: expected_mean
      real(dp) :: expected_var

      call make_reference_model(model)
      model%y(1, :) = [0.9_dp, 0.2_dp, 0.5_dp, 1.1_dp]
      call marss_kfss(model, kf)
      call assert_true('MARSShatyt full filter status', kf%ok, failures)
      if (.not. kf%ok) return
      call marss_hatyt(model, kf, hat, only_kem=.false.)
      call assert_true('MARSShatyt full status', hat%ok .and. hat%info == 0, failures)
      if (.not. hat%ok) return
      expected_mean = model%z(1, 1) * kf%x_pred(1, 1) + model%a(1)
      expected_var = model%r(1, 1) + model%z(1, 1)**2 * kf%p_pred(1, 1, 1)
      call assert_close('MARSShatyt ytt1 mean', hat%ytt1(1, 1), expected_mean, 1.0e-12_dp, failures)
      call assert_close('MARSShatyt ytt1 variance', hat%var_ytt1(1, 1, 1), expected_var, 1.0e-12_dp, failures)
      call assert_close('MARSShatyt ytt observed', hat%ytt(1, 1), model%y(1, 1), 1.0e-12_dp, failures)
      call assert_close('MARSShatyt var yt observed', hat%var_yt(1, 1, 1), 0.0_dp, 1.0e-12_dp, failures)
      call assert_close('MARSShatyt yxtt observed', hat%yxtt(1, 1, 1), &
         model%y(1, 1) * kf%x_filt(1, 1), 1.0e-12_dp, failures)
      call assert_close('MARSShatyt yxt previous smoothed', hat%yxt_prev_smooth(1, 1, 1), &
         model%y(1, 1) * kf%x0_smooth(1), 1.0e-12_dp, failures)
      call assert_true('MARSShatyt final yxtp undefined', ieee_is_nan(hat%yxtp(1, 1, size(model%y, 2))), failures)
   end subroutine test_hatyt_full

   subroutine test_em_fit(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: initial_kf
      type(marss_fit_result) :: fit

      call make_em_model(model)
      call marss_kfss(model, initial_kf)
      call marss_kem(model, 30, 1.0e-7_dp, fit, estimate_z=.false., estimate_a=.false., estimate_v0=.false.)
      call assert_true('MARSSkem status', fit%info == 0, failures)
      if (fit%info == 0) then
         call assert_true('MARSSkem likelihood monotonic endpoint', fit%loglik >= initial_kf%loglik - 1.0e-7_dp, failures)
         call assert_true('MARSSkem finite Q', fit%model%q(1, 1) > 0.0_dp, failures)
         call assert_true('MARSSkem finite R', fit%model%r(1, 1) > 0.0_dp, failures)
      end if
   end subroutine test_em_fit

   subroutine test_missing_em(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_fit_result) :: fit
      real(dp) :: nan

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(model%y(1, 5), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
      model%y(1, :) = [1.0_dp, 2.0_dp, nan, 4.0_dp, 5.0_dp]
      model%b = 0.0_dp
      model%u = 0.0_dp
      model%q = 0.0_dp
      model%z = 0.0_dp
      model%a = 0.0_dp
      model%r = 4.0_dp
      model%x0 = 0.0_dp
      model%v0 = 0.0_dp
      model%tinitx = 1

      call marss_kem(model, 100, 1.0e-10_dp, fit, estimate_b=.false., estimate_u=.false., &
         estimate_q=.false., estimate_z=.false., estimate_a=.true., estimate_r=.true., &
         estimate_x0=.false., estimate_v0=.false.)
      call assert_true('MARSSkem missing-data status', fit%info == 0, failures)
      if (fit%info == 0) then
         call assert_close('MARSSkem missing-data A', fit%model%a(1), 3.0_dp, 2.0e-5_dp, failures)
         call assert_close('MARSSkem missing-data R', fit%model%r(1, 1), 2.5_dp, 5.0e-5_dp, failures)
      end if
   end subroutine test_missing_em

   subroutine test_time_varying_fixed_em(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_fit_result) :: fit
      real(dp) :: states(4)

      allocate(model%y(1, 4), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
      allocate(model%b_t(1, 1, 4), model%u_t(1, 4), model%q_t(1, 1, 4))
      allocate(model%z_t(1, 1, 4), model%r_t(1, 1, 4))
      model%b = 1.0_dp
      model%u = 0.0_dp
      model%q = 0.0_dp
      model%z = 1.0_dp
      model%a = 0.0_dp
      model%r = 1.0_dp
      model%x0 = 1.0_dp
      model%v0 = 0.0_dp
      model%tinitx = 1
      model%b_t(1, 1, :) = [1.0_dp, 2.0_dp, 0.5_dp, 1.5_dp]
      model%u_t(1, :) = [0.0_dp, 0.5_dp, -0.5_dp, 1.0_dp]
      model%q_t = 0.0_dp
      model%z_t(1, 1, :) = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp]
      model%r_t = 1.0_dp
      states = [1.0_dp, 2.5_dp, 0.75_dp, 2.125_dp]
      model%y(1, :) = [2.0_dp, 7.0_dp, 2.25_dp, 3.0625_dp]

      call marss_kem(model, 20, 1.0e-10_dp, fit, estimate_b=.false., estimate_u=.false., &
         estimate_q=.false., estimate_z=.false., estimate_a=.true., estimate_r=.false., &
         estimate_x0=.false., estimate_v0=.false.)
      call assert_true('MARSSkem fixed time-varying status', fit%info == 0, failures)
      if (fit%info == 0) then
         call assert_close('MARSSkem fixed time-varying A', fit%model%a(1), 2.0_dp, 1.0e-10_dp, failures)
         call assert_vector_close('MARSSkem preserves B(t)', fit%model%b_t(1, 1, :), &
            model%b_t(1, 1, :), 0.0_dp, failures)
      end if
   end subroutine test_time_varying_fixed_em

   subroutine test_bootstrap(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      real(dp), allocatable :: params(:, :)
      real(dp), allocatable :: loglik(:)
      integer :: info

      call make_em_model(model)
      call marss_boot(model, 2, 4, 1.0e-5_dp, 20260905_i8, params, loglik, info, &
         estimate_z=.false., estimate_a=.false., estimate_v0=.false.)
      call assert_true('MARSSboot status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSboot dimensions', size(params, 2) == 2 .and. size(loglik) == 2, failures)
      end if
   end subroutine test_bootstrap

   subroutine test_constrained_bootstrap(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_constraints) :: constraints
      type(marss_constraints) :: updated
      type(marss_kf_result) :: kf
      real(dp), allocatable :: bias(:)
      real(dp), allocatable :: boot_data(:, :, :)
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: loglik(:)
      real(dp), allocatable :: hess_params(:, :)
      real(dp), allocatable :: hess_params2(:, :)
      real(dp), allocatable :: params(:, :)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: aicbb
      real(dp) :: aicbp
      integer :: info

      call make_equal_b_model(model, constraints)
      call marss_constraints_set_start(constraints, [0.35_dp], updated, info)
      call assert_true('constraint start update status', info == 0, failures)
      if (info == 0) then
         call assert_close('constraint start update value', updated%b%start(1), 0.35_dp, 0.0_dp, failures)
         call assert_close('constraint start source unchanged', constraints%b%start(1), 0.4_dp, 0.0_dp, failures)
      end if

      call marss_boot(model, 2, 5, 1.0e-4_dp, 6101_i8, params, loglik, info, constraints=constraints, &
         free_parameters=[0.4_dp], fit_method='BFGS')
      call assert_true('MARSSboot constrained BFGS status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSboot constrained beta dimensions', &
            size(params, 1) == 1 .and. size(params, 2) == 2, failures)
         call assert_true('MARSSboot constrained finite beta', all(abs(params) < huge(1.0_dp)), failures)
      end if

      call marss_boot(model, 1, 2, 1.0e-4_dp, 6102_i8, params, loglik, info, constraints=constraints, &
         free_parameters=[0.4_dp], fit_method='kem')
      call assert_true('MARSSboot constrained KEM dispatch', info == 0, failures)

      call marss_boot_hessian(model, 4, 6110_i8, hess_params, info, constraints=constraints, &
         free_parameters=[0.4_dp])
      call assert_true('MARSSboot Hessian parameter generation status', info == 0, failures)
      if (info == 0) then
         call marss_boot_hessian(model, 4, 6110_i8, hess_params2, info, constraints=constraints, &
            free_parameters=[0.4_dp])
         call assert_true('MARSSboot Hessian deterministic repeat status', info == 0, failures)
         if (info == 0) then
            call assert_true('MARSSboot Hessian beta dimensions', all(shape(hess_params) == [1, 4]), failures)
            call assert_matrix_close('MARSSboot Hessian deterministic draws', hess_params, hess_params2, 0.0_dp, failures)
            call assert_true('MARSSboot Hessian nondegenerate draws', &
               maxval(abs(hess_params(:, 2:4) - spread(hess_params(:, 1), 2, 3))) > 0.0_dp, failures)
         end if
      end if

      call marss_boot(model, 2, 5, 1.0e-4_dp, 6105_i8, params, loglik, info, constraints=constraints, &
         free_parameters=[0.4_dp], fit_method='bfgs', simulation_method='innovations', boot_data=boot_data)
      call assert_true('MARSSboot constrained innovations status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSboot innovations data dimensions', &
            all(shape(boot_data) == [2, 8, 2]), failures)
         call assert_true('MARSSboot innovations data finite', all(.not. ieee_is_nan(boot_data)), failures)
      end if

      call marss_kfss(model, kf, smoother=.false.)
      call marss_bootstrap_aic(model, kf%loglik, 1, 5, 1.0e-4_dp, 6103_i8, aicbp, aicbb, info, &
         constraints=constraints, free_parameters=[0.4_dp], fit_method='bfgs')
      call assert_true('MARSSaic constrained bootstrap status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSaic constrained AICbp finite', abs(aicbp) < huge(1.0_dp), failures)
         call assert_true('MARSSaic constrained AICbb finite', abs(aicbb) < huge(1.0_dp), failures)
      end if

      call marss_bootstrap_param_cis(model, 'innovations', 0.1_dp, 2, 5, 1.0e-4_dp, 6104_i8, &
         estimate, se, bias, lower, upper, info, params, constraints=constraints, &
         free_parameters=[0.4_dp], fit_method='bfgs')
      call assert_true('MARSSparamCIs constrained innovations status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSparamCIs constrained beta dimensions', &
            size(estimate) == 1 .and. size(params, 1) == 1, failures)
         call assert_close('MARSSparamCIs constrained fitted beta', estimate(1), 0.4_dp, 0.0_dp, failures)
         call assert_true('MARSSparamCIs constrained interval ordering', lower(1) <= upper(1), failures)
      end if
   end subroutine test_constrained_bootstrap


   subroutine test_bootstrap_aic(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      real(dp), allocatable :: innov_star(:)
      real(dp), allocatable :: param_star(:)
      real(dp) :: aicbb
      real(dp) :: aicbp
      integer :: info

      call make_em_model(model)
      call marss_kfss(model, kf, smoother=.false.)
      call marss_bootstrap_aic(model, kf%loglik, 2, 4, 1.0e-5_dp, 8081_i8, aicbp, aicbb, info, &
         param_star, innov_star, estimate_z=.false., estimate_a=.false., estimate_v0=.false.)
      call assert_true('MARSSaic bootstrap status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSaic bootstrap finite AICbp', abs(aicbp) < huge(1.0_dp), failures)
         call assert_true('MARSSaic bootstrap finite AICbb', abs(aicbb) < huge(1.0_dp), failures)
         call assert_close('MARSSaic AICbp formula', aicbp, &
            -4.0_dp * sum(param_star) / real(size(param_star), dp) + 2.0_dp * kf%loglik, 1.0e-12_dp, failures)
         call assert_close('MARSSaic AICbb formula', aicbb, &
            -4.0_dp * sum(innov_star) / real(size(innov_star), dp) + 2.0_dp * kf%loglik, 1.0e-12_dp, failures)
      end if

      model%y(1, 3) = ieee_value(0.0_dp, ieee_quiet_nan)
      call marss_kfss(model, kf, smoother=.false.)
      call marss_bootstrap_aic(model, kf%loglik, 1, 4, 1.0e-5_dp, 8082_i8, aicbp, aicbb, info, &
         estimate_z=.false., estimate_a=.false., estimate_v0=.false., compute_innovations=.false.)
      call assert_true('MARSSaic missing-data AICbp status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSaic missing-data AICbp finite', abs(aicbp) < huge(1.0_dp), failures)
      end if
   end subroutine test_bootstrap_aic

   subroutine test_innovations_bootstrap(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_innov_boot_result) :: boot1
      type(marss_innov_boot_result) :: boot2
      integer :: info1
      integer :: info2

      call make_em_model(model)
      call marss_innovations_boot(model, 2, 2, 7001_i8, boot1, info1)
      call marss_innovations_boot(model, 2, 2, 7001_i8, boot2, info2)
      call assert_true('MARSSinnovationsboot status', info1 == 0 .and. info2 == 0, failures)
      if (info1 == 0 .and. info2 == 0) then
         call assert_true('MARSSinnovationsboot data shape', size(boot1%data, 3) == 2, failures)
         call assert_true('MARSSinnovationsboot state shape', size(boot1%states, 2) == size(model%y, 2), failures)
         call assert_true('MARSSinnovationsboot deterministic data', &
            maxval(abs(boot1%data - boot2%data)) <= 0.0_dp, failures)
         call assert_true('MARSSinnovationsboot deterministic states', &
            maxval(abs(boot1%states - boot2%states)) <= 0.0_dp, failures)
      end if
   end subroutine test_innovations_bootstrap

   subroutine test_parameter_cis(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      integer :: info

      call make_em_model(model)
      model%b = 0.8_dp
      model%u = 0.1_dp
      model%q = 0.2_dp
      model%r = 0.3_dp
      model%x0 = 0.2_dp
      model%v0 = 0.5_dp
      call marss_param_cis(model, 0.05_dp, estimate, se, lower, upper, info, 1.0e-4_dp)
      call assert_true('MARSSparamCIs status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSparamCIs dimensions', size(estimate) == marss_parameter_count(model), failures)
         call assert_true('MARSSparamCIs nonnegative SE', all(se >= 0.0_dp), failures)
         call assert_true('MARSSparamCIs interval ordering', all(lower <= upper), failures)
      end if
   end subroutine test_parameter_cis

   subroutine test_bootstrap_parameter_cis(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      real(dp), allocatable :: bias(:)
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: params(:, :)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: expected_bias
      real(dp) :: expected_lower
      real(dp) :: expected_se
      real(dp) :: expected_upper
      real(dp) :: mean_value
      real(dp) :: middle_value
      integer :: i
      integer :: info

      call make_em_model(model)
      call marss_bootstrap_param_cis(model, 'PARAMETRIC', 0.5_dp, 3, 4, 1.0e-5_dp, 9091_i8, &
         estimate, se, bias, lower, upper, info, params, estimate_z=.false., estimate_a=.false., &
         estimate_v0=.false.)
      call assert_true('MARSSparamCIs parametric bootstrap status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSparamCIs parametric dimensions', size(params, 2) == 3, failures)
         do i = 1, size(estimate)
            mean_value = sum(params(i, :)) / 3.0_dp
            expected_bias = estimate(i) - mean_value
            expected_se = sqrt(sum((params(i, :) - mean_value)**2) / 2.0_dp)
            call assert_close('MARSSparamCIs bootstrap bias', bias(i), expected_bias, 1.0e-12_dp, failures)
            call assert_close('MARSSparamCIs bootstrap SE', se(i), expected_se, 1.0e-12_dp, failures)
         end do
         middle_value = sum(params(1, :)) - minval(params(1, :)) - maxval(params(1, :))
         expected_lower = 0.5_dp * (minval(params(1, :)) + middle_value)
         expected_upper = 0.5_dp * (middle_value + maxval(params(1, :)))
         call assert_close('MARSSparamCIs type-7 lower', lower(1), expected_lower, 1.0e-12_dp, failures)
         call assert_close('MARSSparamCIs type-7 upper', upper(1), expected_upper, 1.0e-12_dp, failures)
      end if

      call marss_bootstrap_param_cis(model, 'innovations', 0.1_dp, 2, 4, 1.0e-5_dp, 9092_i8, &
         estimate, se, bias, lower, upper, info, estimate_z=.false., estimate_a=.false., estimate_v0=.false.)
      call assert_true('MARSSparamCIs innovations bootstrap status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSSparamCIs innovations nonnegative SE', all(se >= 0.0_dp), failures)
         call assert_true('MARSSparamCIs innovations interval ordering', all(lower <= upper), failures)
      end if
   end subroutine test_bootstrap_parameter_cis

   subroutine test_harvey_information(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_constraints) :: constraints
      real(dp), allocatable :: fisher(:, :)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: upper(:)
      integer :: info

      call make_reference_model(model)
      model%y(1, :) = [0.9_dp, 0.2_dp, 0.5_dp, 1.1_dp]
      allocate(constraints%a%fixed(1), constraints%a%design(1, 1), constraints%a%start(1))
      constraints%a%fixed = 0.0_dp
      constraints%a%design = 1.0_dp
      constraints%a%start = -0.1_dp
      call marss_constraint_start_vector(constraints, theta)
      call marss_fisher_i_harvey_linear(model, constraints, theta, fisher, info)
      call assert_true('Harvey Fisher status', info == 0, failures)
      if (info == 0) then
         call assert_close('Harvey Fisher scalar reference', fisher(1, 1), &
            1.0541595869596039_dp, 3.0e-12_dp, failures)
      end if
      call marss_param_cis_harvey_linear(model, constraints, theta, 0.05_dp, se, lower, upper, info)
      call assert_true('Harvey Fisher CI status', info == 0, failures)
      if (info == 0) then
         call assert_close('Harvey Fisher CI standard error', se(1), &
            1.0_dp / sqrt(1.0541595869596039_dp), 3.0e-12_dp, failures)
         call assert_true('Harvey Fisher CI ordering', lower(1) < theta(1) .and. theta(1) < upper(1), failures)
      end if
   end subroutine test_harvey_information


   subroutine test_hessian_summary(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_constraints) :: constraints
      type(marss_hessian_result) :: harvey
      type(marss_hessian_result) :: fd
      type(marss_hessian_result) :: optim_result
      type(marss_hessian_result) :: invalid
      real(dp), allocatable :: theta(:)

      call make_reference_model(model)
      model%y(1, :) = [0.9_dp, 0.2_dp, 0.5_dp, 1.1_dp]
      allocate(constraints%a%fixed(1), constraints%a%design(1, 1), constraints%a%start(1))
      constraints%a%fixed = 0.0_dp
      constraints%a%design = 1.0_dp
      constraints%a%start = -0.1_dp
      call marss_constraint_start_vector(constraints, theta)

      call marss_hessian_summary_linear(model, constraints, theta, harvey)
      call assert_true('MARSShessian summary Harvey status', harvey%info == 0, failures)
      if (harvey%info == 0) then
         call assert_close('MARSShessian summary Harvey information', harvey%hessian(1, 1), &
            1.0541595869596039_dp, 3.0e-12_dp, failures)
         call assert_true('MARSShessian summary covariance available', harvey%covariance_available, failures)
         if (harvey%covariance_available) then
            call assert_close('MARSShessian summary covariance', harvey%par_sigma(1, 1), &
               1.0_dp / 1.0541595869596039_dp, 3.0e-12_dp, failures)
         end if
         call assert_close('MARSShessian summary parameter mean', harvey%par_mean(1), theta(1), &
            1.0e-15_dp, failures)
         call assert_true('MARSShessian summary parameter name', trim(harvey%parameter_names(1)) == 'A.1', failures)
      end if

      call marss_hessian_summary_linear(model, constraints, theta, fd, method='fdHess', rel_step=1.0e-4_dp)
      call assert_true('MARSShessian summary fdHess status', fd%info == 0, failures)
      call marss_hessian_summary_linear(model, constraints, theta, optim_result, method='optim', rel_step=1.0e-4_dp)
      call assert_true('MARSShessian summary optim status', optim_result%info == 0, failures)
      if (fd%info == 0 .and. optim_result%info == 0) then
         call assert_matrix_close('MARSShessian fdHess optim equivalence', fd%hessian, optim_result%hessian, &
            1.0e-12_dp, failures)
      end if

      call marss_hessian_summary_linear(model, constraints, theta, invalid, method='unknown')
      call assert_true('MARSShessian summary invalid method', invalid%info == 2, failures)
   end subroutine test_hessian_summary

   subroutine test_optim_fit(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: initial_kf
      type(marss_fit_result) :: fit

      call make_em_model(model)
      call marss_kfss(model, initial_kf, smoother=.false.)
      call marss_optim(model, 20, 1.0e-4_dp, fit, estimate_z=.false., estimate_a=.false., estimate_v0=.false.)
      call assert_true('MARSSoptim status', fit%info == 0, failures)
      if (fit%info == 0) then
         call assert_true('MARSSoptim likelihood improvement', fit%loglik > initial_kf%loglik, failures)
         call assert_true('MARSSoptim positive Q', fit%model%q(1, 1) > 0.0_dp, failures)
         call assert_true('MARSSoptim positive R', fit%model%r(1, 1) > 0.0_dp, failures)
      end if
   end subroutine test_optim_fit

   subroutine test_linear_constraints(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: constrained_model
      type(marss_fit_result) :: fit
      type(marss_constraints) :: constraints
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: upper(:)
      integer :: info

      allocate(model%y(2, 8), model%b(2, 2), model%u(2), model%q(2, 2))
      allocate(model%z(2, 2), model%a(2), model%r(2, 2), model%x0(2), model%v0(2, 2))
      model%y(1, :) = [0.2_dp, 0.4_dp, 0.5_dp, 0.7_dp, 0.8_dp, 1.0_dp, 1.1_dp, 1.3_dp]
      model%y(2, :) = [-0.1_dp, 0.0_dp, 0.2_dp, 0.1_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.7_dp]
      model%b = 0.0_dp
      model%b(1, 1) = 0.4_dp
      model%b(2, 2) = 0.4_dp
      model%u = 0.0_dp
      model%q = 0.15_dp * reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%z = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%a = 0.0_dp
      model%r = 0.1_dp * reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%x0 = 0.0_dp
      model%v0 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%tinitx = 0

      allocate(constraints%b%fixed(4), constraints%b%design(4, 1), constraints%b%start(1))
      constraints%b%fixed = 0.0_dp
      constraints%b%design(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
      constraints%b%start = 0.4_dp
      call assert_true('linear constraints valid', marss_constraints_valid(model, constraints), failures)
      call marss_constraint_start_vector(constraints, theta)
      call assert_true('linear constraints parameter count', size(theta) == 1, failures)
      theta(1) = 0.65_dp
      call marss_apply_constraints(model, constraints, theta, constrained_model, info)
      call assert_true('linear constraints application status', info == 0, failures)
      if (info == 0) then
         call assert_close('linear constraints equal diagonal 1', constrained_model%b(1, 1), 0.65_dp, 1.0e-12_dp, failures)
         call assert_close('linear constraints equal diagonal 2', constrained_model%b(2, 2), 0.65_dp, 1.0e-12_dp, failures)
         call assert_close('linear constraints fixed off diagonal', constrained_model%b(1, 2), 0.0_dp, 1.0e-12_dp, failures)
      end if
      call marss_optim_linear(model, constraints, 8, 1.0e-4_dp, fit)
      call assert_true('linear constrained optimizer status', fit%info == 0, failures)
      if (fit%info == 0) then
         call assert_close('linear optimizer equality', fit%model%b(1, 1), fit%model%b(2, 2), 1.0e-12_dp, failures)
         call assert_close('linear optimizer fixed B12', fit%model%b(1, 2), 0.0_dp, 1.0e-12_dp, failures)
         call assert_close('linear optimizer fixed B21', fit%model%b(2, 1), 0.0_dp, 1.0e-12_dp, failures)
         call assert_true('linear optimizer beta retained', allocated(fit%free_parameters), failures)
         if (allocated(fit%free_parameters)) then
            call marss_hessian_linear(model, constraints, fit%free_parameters, hessian, 1.0e-4_dp, info)
            call assert_true('linear constrained Hessian status', info == 0, failures)
            if (info == 0) then
               call assert_true('linear constrained Hessian curvature', &
                  size(hessian, 1) == 1 .and. hessian(1, 1) < 0.0_dp, failures)
            end if
            call marss_param_cis_linear(model, constraints, fit%free_parameters, 0.05_dp, &
               se, lower, upper, info, 1.0e-4_dp)
            call assert_true('linear constrained CI status', info == 0, failures)
            if (info == 0) then
               call assert_true('linear constrained CI dimensions', size(se) == 1, failures)
               call assert_true('linear constrained CI ordering', lower(1) <= upper(1), failures)
            end if
         end if
      end if
   end subroutine test_linear_constraints


   subroutine test_constrained_em(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_kf_result) :: initial_kf
      type(marss_fit_result) :: fit
      type(marss_constraints) :: constraints

      allocate(model%y(2, 8), model%b(2, 2), model%u(2), model%q(2, 2))
      allocate(model%z(2, 2), model%a(2), model%r(2, 2), model%x0(2), model%v0(2, 2))
      model%y(1, :) = [0.2_dp, 0.4_dp, 0.5_dp, 0.7_dp, 0.8_dp, 1.0_dp, 1.1_dp, 1.3_dp]
      model%y(2, :) = [-0.1_dp, 0.0_dp, 0.2_dp, 0.1_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.7_dp]
      model%b = 0.0_dp
      model%b(1, 1) = 0.2_dp
      model%b(2, 2) = 0.2_dp
      model%u = 0.0_dp
      model%q = 0.15_dp * reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%z = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%a = 0.0_dp
      model%r = 0.1_dp * reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%x0 = 0.0_dp
      model%v0 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%tinitx = 0
      allocate(constraints%b%fixed(4), constraints%b%design(4, 1), constraints%b%start(1))
      constraints%b%fixed = 0.0_dp
      constraints%b%design(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
      constraints%b%start = 0.2_dp
      call marss_kfss(model, initial_kf)
      call marss_kem_linear(model, constraints, 12, 1.0e-5_dp, fit, mstep_iter=20)
      call assert_true('constrained EM status', fit%info == 0, failures)
      if (fit%info == 0) then
         call assert_true('constrained EM likelihood nondecrease', fit%loglik >= initial_kf%loglik - 1.0e-9_dp, failures)
         call assert_close('constrained EM equality', fit%model%b(1, 1), fit%model%b(2, 2), 1.0e-12_dp, failures)
         call assert_close('constrained EM fixed B12', fit%model%b(1, 2), 0.0_dp, 1.0e-12_dp, failures)
         call assert_close('constrained EM fixed B21', fit%model%b(2, 1), 0.0_dp, 1.0e-12_dp, failures)
         call assert_true('constrained EM beta retained', allocated(fit%free_parameters), failures)
      end if
   end subroutine test_constrained_em

   subroutine test_cross_validation(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_cv_result) :: cv
      type(marss_cv_result) :: cv_default
      integer :: folds(1, 12)
      integer :: info
      integer :: info_default
      integer :: t

      call make_em_model(model)
      do t = 1, size(folds, 2)
         folds(1, t) = 1 + modulo(t - 1, 2)
      end do
      call marss_cv(model, 15, 1.0e-4_dp, cv, info, fold_ids=folds, estimate_z=.false., &
         estimate_a=.false., estimate_v0=.false.)
      call assert_true('MARSScv status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSScv prediction coverage', all(.not. ieee_is_nan(cv%prediction)), failures)
         call assert_true('MARSScv finite SE', all(cv%se >= 0.0_dp), failures)
         call assert_true('MARSScv fold IDs', all(cv%fold_id == folds), failures)
      end if
      call marss_cv(model, 1, 1.0e-4_dp, cv_default, info_default, &
         estimate_b=.false., estimate_u=.false., estimate_q=.false., estimate_z=.false., &
         estimate_a=.false., estimate_r=.false., estimate_x0=.false., estimate_v0=.false.)
      call assert_true('MARSScv default folds status', info_default == 0, failures)
      if (info_default == 0) then
         call assert_true('MARSScv default prediction coverage', &
            all(.not. ieee_is_nan(cv_default%prediction)), failures)
      end if
   end subroutine test_cross_validation

   subroutine test_constrained_cross_validation(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions.
      type(marss_model) :: model
      type(marss_model) :: reconstructed
      type(marss_constraints) :: constraints
      type(marss_cv_result) :: cv
      real(dp), allocatable :: fold_beta(:, :)
      integer :: folds(2, 8)
      integer :: info
      integer :: k
      integer :: t

      call make_equal_b_model(model, constraints)
      do t = 1, 8
         folds(1, t) = 1 + modulo(t - 1, 2)
         folds(2, t) = 1 + modulo(t, 2)
      end do
      call marss_cv(model, 6, 1.0e-4_dp, cv, info, fold_ids=folds, constraints=constraints, &
         free_parameters=[0.4_dp], fit_method='bfgs', fold_free_parameters=fold_beta)
      call assert_true('MARSScv constrained BFGS status', info == 0, failures)
      if (info == 0) then
         call assert_true('MARSScv constrained prediction coverage', &
            all(.not. ieee_is_nan(cv%prediction)), failures)
         call assert_true('MARSScv constrained beta dimensions', &
            size(fold_beta, 1) == 1 .and. size(fold_beta, 2) == 2, failures)
         do k = 1, size(fold_beta, 2)
            call marss_apply_constraints(model, constraints, fold_beta(:, k), reconstructed, info)
            call assert_true('MARSScv constrained fold reconstruction', info == 0, failures)
            if (info == 0) then
               call assert_close('MARSScv preserves equal B diagonal', &
                  reconstructed%b(1, 1), reconstructed%b(2, 2), 0.0_dp, failures)
               call assert_close('MARSScv preserves fixed B12', reconstructed%b(1, 2), 0.0_dp, 0.0_dp, failures)
               call assert_close('MARSScv preserves fixed B21', reconstructed%b(2, 1), 0.0_dp, 0.0_dp, failures)
            end if
         end do
      end if

      call marss_cv(model, 2, 1.0e-4_dp, cv, info, fold_ids=folds, constraints=constraints, &
         free_parameters=[0.4_dp], fit_method='kem')
      call assert_true('MARSScv constrained KEM dispatch', info == 0, failures)
   end subroutine test_constrained_cross_validation

   subroutine make_equal_b_model(model, constraints)
      type(marss_model), intent(out) :: model !! Two-state model with a shared diagonal B coefficient for constraint tests.
      type(marss_constraints), intent(out) :: constraints !! One-beta affine constraint tying the two B diagonal elements.

      allocate(model%y(2, 8), model%b(2, 2), model%u(2), model%q(2, 2))
      allocate(model%z(2, 2), model%a(2), model%r(2, 2), model%x0(2), model%v0(2, 2))
      model%y(1, :) = [0.2_dp, 0.4_dp, 0.5_dp, 0.7_dp, 0.8_dp, 1.0_dp, 1.1_dp, 1.3_dp]
      model%y(2, :) = [-0.1_dp, 0.0_dp, 0.2_dp, 0.1_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.7_dp]
      model%b = 0.0_dp
      model%b(1, 1) = 0.4_dp
      model%b(2, 2) = 0.4_dp
      model%u = 0.0_dp
      model%q = 0.15_dp * reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%z = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%a = 0.0_dp
      model%r = 0.1_dp * reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%x0 = 0.0_dp
      model%v0 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      model%tinitx = 0
      allocate(constraints%b%fixed(4), constraints%b%design(4, 1), constraints%b%start(1))
      constraints%b%fixed = 0.0_dp
      constraints%b%design(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
      constraints%b%start = 0.4_dp
   end subroutine make_equal_b_model

   subroutine make_free_block(block, values, free_index)
      type(marss_constraint_block), intent(out) :: block !! Affine block matching current values with one element free.
      real(dp), intent(in) :: values(:) !! Flattened current block values used as the affine fixed offset.
      integer, intent(in) :: free_index !! One-based flattened element allowed to vary through a single beta coordinate.

      allocate(block%fixed(size(values)), block%design(size(values), 1), block%start(1))
      block%fixed = values
      block%design = 0.0_dp
      block%start = 0.0_dp
      block%design(free_index, 1) = 1.0_dp
   end subroutine make_free_block

   subroutine make_reference_model(model)
      type(marss_model), intent(out) :: model !! One-state reference model used by deterministic numerical tests.

      allocate(model%y(1, 4), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
      model%y = 0.0_dp
      model%b = 0.8_dp
      model%u = 0.2_dp
      model%q = 0.3_dp
      model%z = 1.2_dp
      model%a = -0.1_dp
      model%r = 0.4_dp
      model%x0 = 0.5_dp
      model%v0 = 0.7_dp
      model%tinitx = 0
   end subroutine make_reference_model

   subroutine make_kem_graph_model(model)
      type(marss_model), intent(out) :: model !! Two-state degenerate model used for KEM graph-restriction tests.
      integer :: t

      allocate(model%y(1, 4), model%b(2, 2), model%u(2), model%q(2, 2))
      allocate(model%z(1, 2), model%a(1), model%r(1, 1), model%x0(2), model%v0(2, 2))
      allocate(model%b_t(2, 2, 4))
      model%y = 0.0_dp
      model%b = 0.0_dp
      model%b(1, 1) = 0.5_dp
      model%b(2, 2) = 0.5_dp
      model%u = 0.0_dp
      model%q = 0.0_dp
      model%q(2, 2) = 0.2_dp
      model%z = reshape([1.0_dp, 0.0_dp], [1, 2])
      model%a = 0.0_dp
      model%r = 0.3_dp
      model%x0 = 0.0_dp
      model%v0 = 0.0_dp
      do t = 1, 4
         model%b_t(:, :, t) = model%b
      end do
   end subroutine make_kem_graph_model

   subroutine make_em_model(model)
      type(marss_model), intent(out) :: model !! One-state complete-data model used to test EM and bootstrap fitting.

      allocate(model%y(1, 12), model%b(1, 1), model%u(1), model%q(1, 1))
      allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
      model%y(1, :) = [0.5_dp, 0.7_dp, 0.2_dp, 1.0_dp, 1.1_dp, 0.9_dp, &
         1.3_dp, 1.0_dp, 1.4_dp, 1.6_dp, 1.5_dp, 1.8_dp]
      model%b = 0.5_dp
      model%u = 0.0_dp
      model%q = 0.5_dp
      model%z = 1.0_dp
      model%a = 0.0_dp
      model%r = 0.5_dp
      model%x0 = 0.0_dp
      model%v0 = 1.0_dp
      model%tinitx = 0
   end subroutine make_em_model

   subroutine assert_true(name, condition, failures)
      character(len=*), intent(in) :: name !! Human-readable assertion label.
      logical, intent(in) :: condition !! Condition expected to be true.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (.not. condition) then
         failures = failures + 1
         write (*, '(a)') 'FAIL: ' // trim(name)
      end if
   end subroutine assert_true

   subroutine assert_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name !! Human-readable assertion label.
      real(dp), intent(in) :: actual !! Computed scalar value.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tolerance !! Allowed absolute error.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         write (*, '(a,2(1x,es16.8))') 'FAIL: ' // trim(name), actual, expected
      end if
   end subroutine assert_close

   subroutine assert_vector_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name !! Human-readable assertion label.
      real(dp), intent(in) :: actual(:) !! Computed vector.
      real(dp), intent(in) :: expected(:) !! Reference vector with the same shape.
      real(dp), intent(in) :: tolerance !! Allowed maximum absolute element error.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (size(actual) /= size(expected)) then
         failures = failures + 1
         write (*, '(a)') 'FAIL: ' // trim(name) // ' shape'
      else if (maxval(abs(actual - expected)) > tolerance) then
         failures = failures + 1
         write (*, '(a,1x,es16.8)') 'FAIL: ' // trim(name), maxval(abs(actual - expected))
      end if
   end subroutine assert_vector_close

   subroutine assert_matrix_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name !! Human-readable assertion label.
      real(dp), intent(in) :: actual(:, :) !! Computed matrix.
      real(dp), intent(in) :: expected(:, :) !! Reference matrix with the same shape.
      real(dp), intent(in) :: tolerance !! Allowed maximum absolute element error.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (any(shape(actual) /= shape(expected))) then
         failures = failures + 1
         write (*, '(a)') 'FAIL: ' // trim(name) // ' shape'
      else if (maxval(abs(actual - expected)) > tolerance) then
         failures = failures + 1
         write (*, '(a,1x,es16.8)') 'FAIL: ' // trim(name), maxval(abs(actual - expected))
      end if
   end subroutine assert_matrix_close

end program test_marss
