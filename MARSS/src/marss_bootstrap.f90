! SPDX-License-Identifier: GPL-2.0-only
module marss_bootstrap
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use marss_kinds, only : dp, i8
   use marss_types, only : marss_model, marss_simulation, marss_fit_result, marss_kf_result
   use marss_types, only : marss_innov_boot_result, marss_constraints
   use marss_simulation_mod, only : marss_simulate
   use marss_em, only : marss_kem
   use marss_analysis, only : marss_parameter_count, marss_vectorizeparam, marss_unvectorizeparam, marss_fisher_i
   use marss_innovations, only : marss_innovations_boot
   use marss_kalman, only : marss_kfss
   use marss_optim_mod, only : marss_optim
   use marss_constraints_mod, only : marss_constraints_parameter_count, marss_constraint_start_vector
   use marss_constraints_mod, only : marss_constraints_set_start, marss_apply_constraints, marss_optim_linear
   use marss_constraints_mod, only : marss_fisher_i_linear, marss_vectorize_free
   use marss_constrained_em, only : marss_kem_linear
   use marss_random, only : marss_rng, rng_seed, sample_mvn
   use r_linalg, only : inverse_matrix, symmetrize
   implicit none
   private
   public :: marss_boot
   public :: marss_boot_hessian
   public :: marss_bootstrap_aic
   public :: marss_bootstrap_param_cis

contains

   subroutine marss_boot(model, nboot, max_iter, tol, seed, boot_params, boot_loglik, info, &
      estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0, &
      constraints, free_parameters, fit_method, simulation_method, boot_data)
      type(marss_model), intent(in) :: model !! Fitted model used as the bootstrap generator and refit template.
      integer, intent(in) :: nboot !! Number of parametric bootstrap replicates, which must be positive.
      integer, intent(in) :: max_iter !! Maximum fitting iterations for each bootstrap refit.
      real(dp), intent(in) :: tol !! Positive relative convergence tolerance for each bootstrap refit.
      integer(i8), intent(in) :: seed !! Base deterministic seed; the replicate index is added for independent streams.
      real(dp), allocatable, intent(out) :: boot_params(:, :) !! Refit parameters by replicate; beta coordinates if constrained.
      real(dp), allocatable, intent(out) :: boot_loglik(:) !! Maximized log likelihood for each successful refit.
      integer, intent(out) :: info !! Zero on success, negative for invalid input, or a positive failing replicate code.
      logical, intent(in), optional :: estimate_b !! Estimate B in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_u !! Estimate U in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_q !! Estimate Q in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_z !! Estimate Z in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_a !! Estimate A in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_r !! Estimate R in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_x0 !! Estimate x0 in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_v0 !! Estimate V0 in unconstrained refits, defaulting to false.
      type(marss_constraints), intent(in), optional :: constraints !! Affine constraints preserved in every bootstrap refit.
      real(dp), intent(in), optional :: free_parameters(:) !! Fitted beta coordinates overriding constraint block start values.
      character(len=*), intent(in), optional :: fit_method !! Refit method, "kem" or "bfgs"; defaults to "kem".
      character(len=*), intent(in), optional :: simulation_method !! Resampling method, "parametric" or "innovations".
      real(dp), allocatable, intent(out), optional :: boot_data(:, :, :) !! Generated data by replicate, preserving missingness.
      type(marss_innov_boot_result) :: innov_boot
      type(marss_simulation) :: sim
      type(marss_model) :: generator
      type(marss_model) :: start
      type(marss_fit_result) :: fit
      real(dp), allocatable :: beta0(:)
      real(dp), allocatable :: theta(:)
      character(len=:), allocatable :: sim_mode
      logical :: ea
      logical :: eb
      logical :: eq
      logical :: er
      logical :: eu
      logical :: ev0
      logical :: ex0
      logical :: ez
      integer :: b
      integer :: fit_info
      integer :: p
      integer :: prepare_info
      integer :: sim_info

      info = 0
      if (nboot < 1) then
         info = -1
         return
      end if
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         info = -2
         return
      end if
      call resolve_estimate_flags(eb, eu, eq, ez, ea, er, ex0, ev0, &
         estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      call prepare_generator(model, generator, p, beta0, prepare_info, constraints, free_parameters)
      if (prepare_info /= 0) then
         info = -10 - prepare_info
         return
      end if

      sim_mode = 'parametric'
      if (present(simulation_method)) sim_mode = lowercase_ascii(trim(adjustl(simulation_method)))
      if (sim_mode /= 'parametric' .and. sim_mode /= 'innovations') then
         info = -3
         return
      end if
      if (sim_mode == 'innovations') then
         call marss_innovations_boot(generator, nboot, 3, seed, innov_boot, sim_info)
         if (sim_info /= 0) then
            info = 1000 + abs(sim_info)
            return
         end if
      end if

      allocate(boot_params(p, nboot), boot_loglik(nboot))
      if (present(boot_data)) allocate(boot_data(size(generator%y, 1), size(generator%y, 2), nboot))
      do b = 1, nboot
         start = generator
         if (sim_mode == 'parametric') then
            call marss_simulate(generator, size(generator%y, 2), 1, seed + int(b, i8), sim, sim_info)
            if (sim_info /= 0) then
               info = 1100 + b
               return
            end if
            start%y = sim%data(:, :, 1)
            where (ieee_is_nan(model%y)) start%y = model%y
         else
            start%y = innov_boot%data(:, :, b)
         end if
         if (present(boot_data)) boot_data(:, :, b) = start%y
         call refit_model(start, max_iter, tol, fit, fit_info, eb, eu, eq, ez, ea, er, ex0, ev0, &
            constraints, beta0, fit_method)
         if (fit_info /= 0 .or. fit%info /= 0) then
            info = 2000 + b
            return
         end if
         if (present(constraints)) then
            if (.not. allocated(fit%free_parameters) .or. size(fit%free_parameters) /= p) then
               info = 3000 + b
               return
            end if
            boot_params(:, b) = fit%free_parameters
         else
            call marss_vectorizeparam(fit%model, theta)
            boot_params(:, b) = theta
            deallocate(theta)
         end if
         boot_loglik(b) = fit%loglik
      end do
   end subroutine marss_boot

   subroutine marss_boot_hessian(model, nboot, seed, boot_params, info, constraints, free_parameters, rel_step)
      type(marss_model), intent(in) :: model !! Fitted numerical model defining the Hessian likelihood and parameter mean.
      integer, intent(in) :: nboot !! Number of normal parameter draws; must be positive.
      integer(i8), intent(in) :: seed !! Deterministic seed for the package-local multivariate-normal generator.
      real(dp), allocatable, intent(out) :: boot_params(:, :) !! Hessian-normal parameter draws by replicate.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid input, Hessian inversion, or MVN generation failure.
      type(marss_constraints), intent(in), optional :: constraints !! Optional affine constraints defining beta-coordinate draws.
      real(dp), intent(in), optional :: free_parameters(:) !! Optional fitted beta mean overriding projection from model values.
      real(dp), intent(in), optional :: rel_step !! Relative finite-difference step for the observed information calculation.
      type(marss_rng) :: rng
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fisher(:, :)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: draw(:)
      integer :: i
      integer :: local_info

      info = 0
      if (nboot <= 0) then
         allocate(boot_params(0, 0))
         info = -1
         return
      end if
      if (present(constraints)) then
         if (present(free_parameters)) then
            if (size(free_parameters) /= marss_constraints_parameter_count(constraints)) then
               allocate(boot_params(0, 0))
               info = -2
               return
            end if
            allocate(mean(size(free_parameters)))
            mean = free_parameters
         else
            call marss_vectorize_free(model, constraints, mean, local_info)
            if (local_info /= 0) then
               allocate(boot_params(0, 0))
               info = -3
               return
            end if
         end if
         if (present(rel_step)) then
            call marss_fisher_i_linear(model, constraints, mean, fisher, rel_step, local_info)
         else
            call marss_fisher_i_linear(model, constraints, mean, fisher, info=local_info)
         end if
      else
         call marss_vectorizeparam(model, mean)
         if (present(rel_step)) then
            call marss_fisher_i(model, fisher, rel_step, local_info)
         else
            call marss_fisher_i(model, fisher, info=local_info)
         end if
      end if
      if (local_info /= 0 .or. size(fisher, 1) /= size(mean)) then
         allocate(boot_params(0, 0))
         info = 1
         return
      end if
      call inverse_matrix(fisher, covariance, local_info)
      if (local_info /= 0) then
         allocate(boot_params(0, 0))
         info = 2
         return
      end if
      covariance = symmetrize(covariance)
      allocate(boot_params(size(mean), nboot), draw(size(mean)))
      call rng_seed(rng, seed)
      do i = 1, nboot
         call sample_mvn(mean, covariance, rng, draw, local_info)
         if (local_info /= 0) then
            info = 100 + i
            return
         end if
         boot_params(:, i) = draw
      end do
   end subroutine marss_boot_hessian

   subroutine marss_bootstrap_param_cis(model, method, alpha, nboot, max_iter, tol, seed, estimate, se, bias, &
      lower, upper, info, bootstrap_parameters, estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, &
      estimate_r, estimate_x0, estimate_v0, constraints, free_parameters, fit_method)
      type(marss_model), intent(in) :: model !! Fitted model used as the bootstrap generator and refit template.
      character(len=*), intent(in) :: method !! Bootstrap method: "parametric" or "innovations".
      real(dp), intent(in) :: alpha !! Two-sided significance level in the open interval (0,1).
      integer, intent(in) :: nboot !! Number of bootstrap refits; at least two are required for a sample standard error.
      integer, intent(in) :: max_iter !! Maximum fitting iterations for each bootstrap refit.
      real(dp), intent(in) :: tol !! Positive relative convergence tolerance used for each bootstrap refit.
      integer(i8), intent(in) :: seed !! Deterministic bootstrap seed.
      real(dp), allocatable, intent(out) :: estimate(:) !! Fitted parameters; beta coordinates when constraints are supplied.
      real(dp), allocatable, intent(out) :: se(:) !! Bootstrap sample standard deviations of the fitted coordinates.
      real(dp), allocatable, intent(out) :: bias(:) !! Fitted estimate minus bootstrap mean, matching MARSSparamCIs.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower empirical R type-7 bootstrap confidence limits.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper empirical R type-7 bootstrap confidence limits.
      integer, intent(out) :: info !! Zero on success; negative values are invalid input and positive values identify failed refits.
      real(dp), allocatable, intent(out), optional :: bootstrap_parameters(:, :) !! Refit coordinates by bootstrap replicate.
      logical, intent(in), optional :: estimate_b !! Estimate B in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_u !! Estimate U in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_q !! Estimate Q in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_z !! Estimate Z in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_a !! Estimate A in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_r !! Estimate R in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_x0 !! Estimate x0 in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_v0 !! Estimate V0 in unconstrained refits, defaulting to false.
      type(marss_constraints), intent(in), optional :: constraints !! Affine constraints preserved in all bootstrap refits.
      real(dp), intent(in), optional :: free_parameters(:) !! Fitted beta coordinates overriding constraint start values.
      character(len=*), intent(in), optional :: fit_method !! Refit method, "kem" or "bfgs"; defaults to "kem".
      type(marss_model) :: generator
      real(dp), allocatable :: beta0(:)
      real(dp), allocatable :: boot_loglik(:)
      real(dp), allocatable :: params(:, :)
      real(dp) :: mean_value
      logical :: ea
      logical :: eb
      logical :: eq
      logical :: er
      logical :: eu
      logical :: ev0
      logical :: ex0
      logical :: ez
      character(len=:), allocatable :: mode
      integer :: boot_info
      integer :: i
      integer :: p
      integer :: prepare_info

      info = 0
      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         info = -1
         return
      end if
      if (nboot < 2) then
         info = -2
         return
      end if
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         info = -3
         return
      end if
      mode = lowercase_ascii(trim(adjustl(method)))
      if (mode /= 'parametric' .and. mode /= 'innovations') then
         info = -4
         return
      end if

      call resolve_estimate_flags(eb, eu, eq, ez, ea, er, ex0, ev0, &
         estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      call prepare_generator(model, generator, p, beta0, prepare_info, constraints, free_parameters)
      if (prepare_info /= 0) then
         info = -10 - prepare_info
         return
      end if

      call marss_boot(model, nboot, max_iter, tol, seed, params, boot_loglik, boot_info, &
         eb, eu, eq, ez, ea, er, ex0, ev0, constraints, beta0, fit_method, simulation_method=mode)
      if (boot_info /= 0) then
         info = 1000 + abs(boot_info)
         return
      end if

      if (present(constraints)) then
         allocate(estimate(p))
         estimate = beta0
      else
         call marss_vectorizeparam(generator, estimate)
      end if
      allocate(se(p), bias(p), lower(p), upper(p))
      do i = 1, p
         mean_value = sum(params(i, :)) / real(nboot, dp)
         bias(i) = estimate(i) - mean_value
         se(i) = sqrt(sum((params(i, :) - mean_value)**2) / real(nboot - 1, dp))
         lower(i) = quantile_type7(params(i, :), 0.5_dp * alpha)
         upper(i) = quantile_type7(params(i, :), 1.0_dp - 0.5_dp * alpha)
      end do

      if (present(bootstrap_parameters)) then
         allocate(bootstrap_parameters(p, nboot))
         bootstrap_parameters = params
      end if
   end subroutine marss_bootstrap_param_cis

   subroutine marss_bootstrap_aic(model, fitted_loglik, nboot, max_iter, tol, seed, aicbp, aicbb, info, &
      loglik_star_parametric, loglik_star_innovations, estimate_b, estimate_u, estimate_q, estimate_z, &
      estimate_a, estimate_r, estimate_x0, estimate_v0, compute_parametric, compute_innovations, &
      constraints, free_parameters, fit_method)
      type(marss_model), intent(in) :: model !! Fitted model used for original-data likelihood and bootstrap generation.
      real(dp), intent(in) :: fitted_loglik !! Maximized log likelihood of the fitted model on the original observations.
      integer, intent(in) :: nboot !! Number of bootstrap refits for each requested correction, which must be positive.
      integer, intent(in) :: max_iter !! Maximum fitting iterations for every bootstrap refit.
      real(dp), intent(in) :: tol !! Positive relative convergence tolerance used in bootstrap refits.
      integer(i8), intent(in) :: seed !! Base deterministic seed; innovations use a fixed offset from this stream.
      real(dp), intent(out) :: aicbp !! Parametric-bootstrap AIC correction, -4 mean(logL.star) + 2 fitted log likelihood.
      real(dp), intent(out) :: aicbb !! Innovations-bootstrap AIC correction, -4 mean(logL.star) + 2 fitted log likelihood.
      integer, intent(out) :: info !! Zero on success; other values identify invalid input or failed bootstrap stages.
      real(dp), allocatable, intent(out), optional :: loglik_star_parametric(:) !! Parametric original-data logL.star values.
      real(dp), allocatable, intent(out), optional :: loglik_star_innovations(:) !! Innovations original-data logL.star values.
      logical, intent(in), optional :: estimate_b !! Estimate B in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_u !! Estimate U in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_q !! Estimate Q in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_z !! Estimate Z in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_a !! Estimate A in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_r !! Estimate R in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_x0 !! Estimate x0 in unconstrained refits, defaulting to true.
      logical, intent(in), optional :: estimate_v0 !! Estimate V0 in unconstrained refits, defaulting to false.
      logical, intent(in), optional :: compute_parametric !! Compute AICbp and parametric logL.star, defaulting to true.
      logical, intent(in), optional :: compute_innovations !! Compute AICbb and innovations logL.star, defaulting to true.
      type(marss_constraints), intent(in), optional :: constraints !! Affine constraints preserved in every bootstrap refit.
      real(dp), intent(in), optional :: free_parameters(:) !! Fitted beta coordinates overriding constraint start values.
      character(len=*), intent(in), optional :: fit_method !! Refit method, "kem" or "bfgs"; defaults to "kem".
      type(marss_kf_result) :: kf
      type(marss_model) :: candidate
      type(marss_model) :: generator
      real(dp), allocatable :: beta0(:)
      real(dp), allocatable :: boot_loglik(:)
      real(dp), allocatable :: boot_params(:, :)
      real(dp), allocatable :: star_innov(:)
      real(dp), allocatable :: star_param(:)
      logical :: do_innovations
      logical :: do_parametric
      logical :: ea
      logical :: eb
      logical :: eq
      logical :: er
      logical :: eu
      logical :: ev0
      logical :: ex0
      logical :: ez
      integer :: b
      integer :: boot_info
      integer :: p
      integer :: prepare_info
      integer :: unpack_info
      integer(i8), parameter :: innovations_seed_offset = 104729_i8

      info = 0
      aicbp = huge(1.0_dp)
      aicbb = huge(1.0_dp)
      if (nboot < 1) then
         info = -1
         return
      end if
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         info = -2
         return
      end if

      do_parametric = .true.
      do_innovations = .true.
      if (present(compute_parametric)) do_parametric = compute_parametric
      if (present(compute_innovations)) do_innovations = compute_innovations
      if (.not. do_parametric .and. .not. do_innovations) then
         info = -3
         return
      end if
      call resolve_estimate_flags(eb, eu, eq, ez, ea, er, ex0, ev0, &
         estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      call prepare_generator(model, generator, p, beta0, prepare_info, constraints, free_parameters)
      if (prepare_info /= 0) then
         info = -10 - prepare_info
         return
      end if

      if (do_parametric) then
         call marss_boot(model, nboot, max_iter, tol, seed, boot_params, boot_loglik, boot_info, &
            eb, eu, eq, ez, ea, er, ex0, ev0, constraints, beta0, fit_method, simulation_method='parametric')
         if (boot_info /= 0) then
            info = 1000 + abs(boot_info)
            return
         end if
         allocate(star_param(nboot))
         do b = 1, nboot
            if (present(constraints)) then
               call marss_apply_constraints(model, constraints, boot_params(:, b), candidate, unpack_info)
            else
               call marss_unvectorizeparam(model, boot_params(:, b), candidate, unpack_info)
            end if
            if (unpack_info /= 0) then
               info = 1100 + b
               return
            end if
            call marss_kfss(candidate, kf, smoother=.false.)
            if (.not. kf%ok) then
               info = 1200 + b
               return
            end if
            star_param(b) = kf%loglik
         end do
         aicbp = -4.0_dp * sum(star_param) / real(nboot, dp) + 2.0_dp * fitted_loglik
      end if

      if (do_innovations) then
         call marss_boot(model, nboot, max_iter, tol, seed + innovations_seed_offset, boot_params, boot_loglik, &
            boot_info, eb, eu, eq, ez, ea, er, ex0, ev0, constraints, beta0, fit_method, &
            simulation_method='innovations')
         if (boot_info /= 0) then
            info = 2000 + abs(boot_info)
            return
         end if
         allocate(star_innov(nboot))
         do b = 1, nboot
            if (present(constraints)) then
               call marss_apply_constraints(model, constraints, boot_params(:, b), candidate, unpack_info)
            else
               call marss_unvectorizeparam(model, boot_params(:, b), candidate, unpack_info)
            end if
            if (unpack_info /= 0) then
               info = 2200 + b
               return
            end if
            call marss_kfss(candidate, kf, smoother=.false.)
            if (.not. kf%ok) then
               info = 2300 + b
               return
            end if
            star_innov(b) = kf%loglik
         end do
         aicbb = -4.0_dp * sum(star_innov) / real(nboot, dp) + 2.0_dp * fitted_loglik
      end if

      if (present(loglik_star_parametric) .and. allocated(star_param)) then
         allocate(loglik_star_parametric(nboot))
         loglik_star_parametric = star_param
      end if
      if (present(loglik_star_innovations) .and. allocated(star_innov)) then
         allocate(loglik_star_innovations(nboot))
         loglik_star_innovations = star_innov
      end if
   end subroutine marss_bootstrap_aic

   subroutine prepare_generator(model, generator, p, beta0, info, constraints, free_parameters)
      type(marss_model), intent(in) :: model !! Input fitted model used as the unconstrained generator or constrained template.
      type(marss_model), intent(out) :: generator !! Generator made consistent with the supplied constrained beta coordinates.
      integer, intent(out) :: p !! Number of returned bootstrap coordinates.
      real(dp), allocatable, intent(out) :: beta0(:) !! Constrained fitted beta vector, or an empty vector when unconstrained.
      integer, intent(out) :: info !! Zero on success, or nonzero for inconsistent constrained coordinates.
      type(marss_constraints), intent(in), optional :: constraints !! Optional affine constraint specification.
      real(dp), intent(in), optional :: free_parameters(:) !! Optional fitted beta coordinates overriding block starts.
      integer :: apply_info

      generator = model
      if (.not. present(constraints)) then
         if (present(free_parameters)) then
            if (size(free_parameters) /= 0) then
               info = 1
               return
            end if
         end if
         p = marss_parameter_count(model)
         allocate(beta0(0))
         info = 0
         return
      end if

      p = marss_constraints_parameter_count(constraints)
      if (present(free_parameters)) then
         if (size(free_parameters) /= p) then
            info = 2
            return
         end if
         allocate(beta0(p))
         beta0 = free_parameters
      else
         call marss_constraint_start_vector(constraints, beta0)
      end if
      call marss_apply_constraints(model, constraints, beta0, generator, apply_info)
      if (apply_info /= 0) then
         info = 10 + apply_info
         return
      end if
      info = 0
   end subroutine prepare_generator

   subroutine refit_model(start_model, max_iter, tol, fit, info, eb, eu, eq, ez, ea, er, ex0, ev0, &
      constraints, free_parameters, fit_method)
      type(marss_model), intent(in) :: start_model !! Model containing bootstrap or cross-validation observations to refit.
      integer, intent(in) :: max_iter !! Maximum EM or BFGS iterations for this refit.
      real(dp), intent(in) :: tol !! Positive convergence tolerance for this refit.
      type(marss_fit_result), intent(out) :: fit !! Fitted model, likelihood, filter output, and constrained beta coordinates.
      integer, intent(out) :: info !! Zero on valid dispatch, or nonzero for invalid method/constraint coordinates.
      logical, intent(in) :: eb !! Whether unconstrained B is estimated.
      logical, intent(in) :: eu !! Whether unconstrained U is estimated.
      logical, intent(in) :: eq !! Whether unconstrained Q is estimated.
      logical, intent(in) :: ez !! Whether unconstrained Z is estimated.
      logical, intent(in) :: ea !! Whether unconstrained A is estimated.
      logical, intent(in) :: er !! Whether unconstrained R is estimated.
      logical, intent(in) :: ex0 !! Whether unconstrained x0 is estimated.
      logical, intent(in) :: ev0 !! Whether unconstrained V0 is estimated.
      type(marss_constraints), intent(in), optional :: constraints !! Optional affine constraints replacing block estimate flags.
      real(dp), intent(in), optional :: free_parameters(:) !! Starting constrained beta coordinates.
      character(len=*), intent(in), optional :: fit_method !! Fitting method, "kem" or "bfgs"; default is "kem".
      type(marss_constraints) :: fit_constraints
      character(len=:), allocatable :: mode
      integer :: start_info

      mode = 'kem'
      if (present(fit_method)) mode = lowercase_ascii(trim(adjustl(fit_method)))
      if (mode /= 'kem' .and. mode /= 'bfgs') then
         info = 1
         return
      end if

      if (present(constraints)) then
         if (present(free_parameters)) then
            call marss_constraints_set_start(constraints, free_parameters, fit_constraints, start_info)
         else
            fit_constraints = constraints
            start_info = 0
         end if
         if (start_info /= 0) then
            info = 2
            return
         end if
         if (mode == 'kem') then
            call marss_kem_linear(start_model, fit_constraints, max_iter, tol, fit)
         else
            call marss_optim_linear(start_model, fit_constraints, max_iter, tol, fit)
         end if
      else
         if (present(free_parameters)) then
            if (size(free_parameters) /= 0) then
               info = 3
               return
            end if
         end if
         if (mode == 'kem') then
            call marss_kem(start_model, max_iter, tol, fit, eb, eu, eq, ez, ea, er, ex0, ev0)
         else
            call marss_optim(start_model, max_iter, tol, fit, eb, eu, eq, ez, ea, er, ex0, ev0)
         end if
      end if
      info = 0
   end subroutine refit_model

   pure subroutine resolve_estimate_flags(eb, eu, eq, ez, ea, er, ex0, ev0, &
      estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      logical, intent(out) :: eb !! Resolved unconstrained B estimate flag.
      logical, intent(out) :: eu !! Resolved unconstrained U estimate flag.
      logical, intent(out) :: eq !! Resolved unconstrained Q estimate flag.
      logical, intent(out) :: ez !! Resolved unconstrained Z estimate flag.
      logical, intent(out) :: ea !! Resolved unconstrained A estimate flag.
      logical, intent(out) :: er !! Resolved unconstrained R estimate flag.
      logical, intent(out) :: ex0 !! Resolved unconstrained x0 estimate flag.
      logical, intent(out) :: ev0 !! Resolved unconstrained V0 estimate flag.
      logical, intent(in), optional :: estimate_b !! Optional B estimate override.
      logical, intent(in), optional :: estimate_u !! Optional U estimate override.
      logical, intent(in), optional :: estimate_q !! Optional Q estimate override.
      logical, intent(in), optional :: estimate_z !! Optional Z estimate override.
      logical, intent(in), optional :: estimate_a !! Optional A estimate override.
      logical, intent(in), optional :: estimate_r !! Optional R estimate override.
      logical, intent(in), optional :: estimate_x0 !! Optional x0 estimate override.
      logical, intent(in), optional :: estimate_v0 !! Optional V0 estimate override.

      eb = .true.
      eu = .true.
      eq = .true.
      ez = .true.
      ea = .true.
      er = .true.
      ex0 = .true.
      ev0 = .false.
      if (present(estimate_b)) eb = estimate_b
      if (present(estimate_u)) eu = estimate_u
      if (present(estimate_q)) eq = estimate_q
      if (present(estimate_z)) ez = estimate_z
      if (present(estimate_a)) ea = estimate_a
      if (present(estimate_r)) er = estimate_r
      if (present(estimate_x0)) ex0 = estimate_x0
      if (present(estimate_v0)) ev0 = estimate_v0
   end subroutine resolve_estimate_flags

   pure function quantile_type7(values, probability) result(value)
      real(dp), intent(in) :: values(:) !! Sample values whose empirical R type-7 quantile is requested.
      real(dp), intent(in) :: probability !! Quantile probability, conventionally between zero and one inclusive.
      real(dp) :: value
      real(dp), allocatable :: sorted(:)
      real(dp) :: fraction
      real(dp) :: position
      real(dp) :: p
      integer :: j
      integer :: n

      n = size(values)
      if (n < 1) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      p = min(max(probability, 0.0_dp), 1.0_dp)
      allocate(sorted(n))
      sorted = values
      call insertion_sort(sorted)
      if (n == 1) then
         value = sorted(1)
         return
      end if
      position = 1.0_dp + real(n - 1, dp) * p
      j = int(floor(position))
      fraction = position - real(j, dp)
      if (j >= n) then
         value = sorted(n)
      else
         value = (1.0_dp - fraction) * sorted(j) + fraction * sorted(j + 1)
      end if
   end function quantile_type7

   pure subroutine insertion_sort(values)
      real(dp), intent(inout) :: values(:) !! Values sorted into ascending order in place.
      real(dp) :: key
      integer :: i
      integer :: j

      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine insertion_sort

   pure function lowercase_ascii(text) result(lower)
      character(len=*), intent(in) :: text !! ASCII text converted to lowercase without locale-dependent behavior.
      character(len=len(text)) :: lower
      integer :: code
      integer :: i

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
   end function lowercase_ascii

end module marss_bootstrap
