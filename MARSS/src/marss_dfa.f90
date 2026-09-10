! SPDX-License-Identifier: GPL-2.0-only
module marss_dfa_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_dfa_spec, marss_constraints, marss_fit_result
   use marss_utils, only : identity_matrix
   use marss_constraints_mod, only : marss_constraints_valid, marss_constraint_start_vector
   use marss_constraints_mod, only : marss_apply_constraints
   use marss_constrained_em, only : marss_kem_linear
   use marss_constraints_mod, only : marss_optim_linear
   use marss_builder, only : canonical_spec, make_matrix_block, make_vector_block
   use marss_builder, only : make_fixed_block, make_fixed_vector
   implicit none
   private
   public :: marss_dfa_build
   public :: marss_dfa_model
   public :: marss_dfa_fit
   public :: marss_dfa_fit_spec

contains

   pure subroutine marss_dfa_build(data, spec, model, constraints, info, covariates, center, scale)
      real(dp), intent(in) :: data(:, :) !! Observation matrix with variables by time; NaNs denote missing observations.
      type(marss_dfa_spec), intent(in) :: spec !! DFA structure, preprocessing choices, and supported model shortcuts.
      type(marss_model), intent(out) :: model !! Numerical DFA model after preprocessing and applying starting free coordinates.
      type(marss_constraints), intent(out) :: constraints !! DFA affine constraints for Z and requested B/Q/A/R/x0/D blocks.
      integer, intent(out) :: info !! Zero on success; positive values identify invalid dimensions or unsupported DFA options.
      real(dp), intent(in), optional :: covariates(:, :) !! Complete observation covariates with predictors by time.
      real(dp), allocatable, intent(out), optional :: center(:) !! Applied observation-row centers.
      real(dp), allocatable, intent(out), optional :: scale(:) !! Applied observation-row sample-standard-deviation scales.
      type(marss_model) :: template
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: transformed(:, :)
      real(dp), allocatable :: target(:, :)
      real(dp), allocatable :: target_vector(:)
      real(dp) :: row_center(size(data, 1))
      real(dp) :: row_mean(size(data, 1))
      real(dp) :: row_scale(size(data, 1))
      character(len=32) :: aspec
      character(len=32) :: bspec
      character(len=32) :: dspec
      character(len=32) :: qspec
      character(len=32) :: rspec
      character(len=32) :: v0spec
      character(len=32) :: x0spec
      integer :: i
      integer :: idx
      integer :: j
      integer :: k
      integer :: m
      integer :: n
      integer :: nfree_z
      integer :: nobs
      integer :: p
      integer :: t
      integer :: tt

      info = 0
      n = size(data, 1)
      tt = size(data, 2)
      m = spec%ntrends
      if (n < 1 .or. tt < 1 .or. m < 1 .or. m > n) then
         info = 1
         return
      end if
      if (spec%tinitx /= 0 .and. spec%tinitx /= 1) then
         info = 2
         return
      end if
      if (spec%diffuse .and. spec%tinitx /= 1) then
         info = 3
         return
      end if
      if (present(covariates)) then
         if (size(covariates, 1) < 1 .or. size(covariates, 2) /= tt) then
            info = 4
            return
         end if
         if (any(ieee_is_nan(covariates))) then
            info = 5
            return
         end if
      end if

      allocate(transformed(n, tt))
      transformed = data
      row_center = 0.0_dp
      row_scale = 1.0_dp
      row_mean = 0.0_dp
      do i = 1, n
         nobs = count(.not. ieee_is_nan(data(i, :)))
         if (nobs < 1) then
            info = 6
            return
         end if
         do t = 1, tt
            if (.not. ieee_is_nan(data(i, t))) row_mean(i) = row_mean(i) + data(i, t)
         end do
         row_mean(i) = row_mean(i) / real(nobs, dp)
         if (spec%demean) then
            row_center(i) = row_mean(i)
            do t = 1, tt
               if (.not. ieee_is_nan(transformed(i, t))) then
                  transformed(i, t) = transformed(i, t) - row_mean(i)
               end if
            end do
         end if
         if (spec%z_score) then
            if (nobs < 2) then
               info = 7
               return
            end if
            row_scale(i) = 0.0_dp
            do t = 1, tt
               if (.not. ieee_is_nan(data(i, t))) then
                  row_scale(i) = row_scale(i) + (data(i, t) - row_mean(i))**2
               end if
            end do
            row_scale(i) = sqrt(row_scale(i) / real(nobs - 1, dp))
            if (row_scale(i) <= sqrt(epsilon(1.0_dp))) then
               info = 8
               return
            end if
            do t = 1, tt
               if (.not. ieee_is_nan(transformed(i, t))) then
                  transformed(i, t) = transformed(i, t) / row_scale(i)
               end if
            end do
         end if
      end do
      if (present(center)) then
         allocate(center(n))
         center = row_center
      end if
      if (present(scale)) then
         allocate(scale(n))
         scale = row_scale
      end if

      allocate(template%y(n, tt), template%b(m, m), template%u(m), template%q(m, m))
      allocate(template%z(n, m), template%a(n), template%r(n, n))
      allocate(template%x0(m), template%v0(m, m))
      template%y = transformed
      template%b = identity_matrix(m)
      template%u = 0.0_dp
      template%q = identity_matrix(m)
      template%z = 0.0_dp
      template%a = 0.0_dp
      template%r = 0.5_dp * identity_matrix(n)
      template%x0 = 0.0_dp
      template%v0 = 5.0_dp * identity_matrix(m)
      template%tinitx = spec%tinitx
      template%diffuse = spec%diffuse

      bspec = canonical_spec(spec%b)
      select case (trim(bspec))
      case ('identity', 'diagonal and equal', 'diagonal and unequal')
         target = identity_matrix(m)
         call make_matrix_block(bspec, target, .false., constraints%b, info)
      case default
         info = 20
         return
      end select
      if (info /= 0) then
         info = 21
         return
      end if
      call make_fixed_vector(template%u, constraints%u)

      qspec = canonical_spec(spec%q)
      select case (trim(qspec))
      case ('identity', 'diagonal and equal', 'diagonal and unequal')
         target = identity_matrix(m)
         call make_matrix_block(qspec, target, .true., constraints%q, info)
      case default
         info = 30
         return
      end select
      if (info /= 0) then
         info = 31
         return
      end if

      nfree_z = n * m - m * (m - 1) / 2
      allocate(constraints%z%fixed(n * m), constraints%z%design(n * m, nfree_z))
      allocate(constraints%z%start(nfree_z))
      constraints%z%fixed = 0.0_dp
      constraints%z%design = 0.0_dp
      constraints%z%start = 0.0_dp
      k = 0
      do j = 1, m
         do i = 1, n
            if (i < j) cycle
            k = k + 1
            idx = i + (j - 1) * n
            constraints%z%design(idx, k) = 1.0_dp
            if (i == j) then
               constraints%z%start(k) = 1.0_dp
            else
               constraints%z%start(k) = 0.1_dp
            end if
         end do
      end do

      aspec = canonical_spec(spec%a)
      select case (trim(aspec))
      case ('zero', 'unconstrained', 'unequal')
         allocate(target_vector(n))
         target_vector = 0.0_dp
         call make_vector_block(aspec, target_vector, constraints%a, info)
      case default
         info = 40
         return
      end select
      if (info /= 0) then
         info = 41
         return
      end if

      rspec = canonical_spec(spec%r)
      select case (trim(rspec))
      case ('identity', 'zero', 'unconstrained', 'unequal', 'diagonal and unequal', &
         'diagonal and equal', 'equalvarcov')
         target = 0.5_dp * identity_matrix(n)
         call make_matrix_block(rspec, target, .true., constraints%r, info)
      case default
         info = 50
         return
      end select
      if (info /= 0) then
         info = 51
         return
      end if

      x0spec = canonical_spec(spec%x0)
      select case (trim(x0spec))
      case ('zero', 'unconstrained', 'unequal')
         if (allocated(target_vector)) deallocate(target_vector)
         allocate(target_vector(m))
         target_vector = 0.0_dp
         call make_vector_block(x0spec, target_vector, constraints%x0, info)
      case default
         info = 60
         return
      end select
      if (info /= 0) then
         info = 61
         return
      end if

      v0spec = canonical_spec(spec%v0)
      select case (trim(v0spec))
      case ('fixed5', 'default')
         target = 5.0_dp * identity_matrix(m)
         call make_fixed_block(target, constraints%v0)
      case ('identity')
         target = identity_matrix(m)
         call make_fixed_block(target, constraints%v0)
      case ('zero')
         target = 0.0_dp * identity_matrix(m)
         call make_fixed_block(target, constraints%v0)
      case default
         info = 70
         return
      end select

      dspec = canonical_spec(spec%d)
      if (present(covariates)) then
         p = size(covariates, 1)
         allocate(template%obs_covariates(p, tt), template%d_coef(n, p))
         template%obs_covariates = covariates
         template%d_coef = 0.0_dp
         if (dspec == 'auto') dspec = 'unconstrained'
         select case (trim(dspec))
         case ('identity', 'zero', 'unconstrained', 'unequal', 'equal', 'diagonal and unequal', &
            'diagonal and equal', 'equalvarcov')
            target = template%d_coef
            call make_matrix_block(dspec, target, .false., constraints%d, info)
         case default
            info = 80
            return
         end select
         if (info /= 0) then
            info = 81
            return
         end if
      else if (dspec /= 'auto' .and. dspec /= 'zero') then
         info = 82
         return
      end if

      call marss_constraint_start_vector(constraints, theta)
      call marss_apply_constraints(template, constraints, theta, model, info)
      if (info /= 0) then
         info = 90 + info
         return
      end if
      if (.not. marss_constraints_valid(model, constraints)) then
         info = 100
         return
      end if
   end subroutine marss_dfa_build

   pure subroutine marss_dfa_model(data, ntrends, model, constraints, info, covariates, demean, z_score, center, scale)
      real(dp), intent(in) :: data(:, :) !! Observation matrix with variables by time; NaNs denote missing observations.
      integer, intent(in) :: ntrends !! Number of DFA random-walk trends, constrained to 1 through the observation dimension.
      type(marss_model), intent(out) :: model !! Numerical DFA model using standard MARSS DFA defaults and transformed data.
      type(marss_constraints), intent(out) :: constraints !! Free Z, equal-diagonal R, and optional free D coefficient constraints.
      integer, intent(out) :: info !! Zero on success; nonzero for dimensions, covariates, degenerate scaling, or invalid setup.
      real(dp), intent(in), optional :: covariates(:, :) !! Optional complete covariate matrix with predictors by time.
      logical, intent(in), optional :: demean !! If true, subtract each variable's observed mean; defaults to true.
      logical, intent(in), optional :: z_score !! If true, divide by each variable's sample standard deviation; defaults to true.
      real(dp), allocatable, intent(out), optional :: center(:) !! Applied row centers, zero when demeaning is disabled.
      real(dp), allocatable, intent(out), optional :: scale(:) !! Applied row scales, one when z scoring is disabled.
      type(marss_dfa_spec) :: spec

      spec%ntrends = ntrends
      if (present(demean)) spec%demean = demean
      if (present(z_score)) spec%z_score = z_score
      call marss_dfa_build(data, spec, model, constraints, info, covariates, center, scale)
   end subroutine marss_dfa_model

   subroutine marss_dfa_fit(data, ntrends, max_iter, tol, fit, constraints, info, covariates, demean, z_score, &
      mstep_iter, method)
      real(dp), intent(in) :: data(:, :) !! Observation matrix with variables by time; NaNs denote missing observations.
      integer, intent(in) :: ntrends !! Number of DFA random-walk trends, constrained to 1 through the observation dimension.
      integer, intent(in) :: max_iter !! Maximum generalized-EM or BFGS iterations, which must be positive.
      real(dp), intent(in) :: tol !! Relative convergence tolerance, which must be positive.
      type(marss_fit_result), intent(out) :: fit !! Fitted DFA model, free coefficients, Kalman result, and likelihood.
      type(marss_constraints), intent(out) :: constraints !! DFA affine constraints used for the fit.
      integer, intent(out) :: info !! Zero on success; otherwise model-construction or fitting status.
      real(dp), intent(in), optional :: covariates(:, :) !! Optional complete covariate matrix with predictors by time.
      logical, intent(in), optional :: demean !! If true, subtract row means before fitting; defaults to true.
      logical, intent(in), optional :: z_score !! If true, standardize row sample variances before fitting; defaults to true.
      integer, intent(in), optional :: mstep_iter !! Maximum inner generalized-EM BFGS iterations, defaulting internally.
      character(len=*), intent(in), optional :: method !! Fitting method: kem/em or bfgs/optim; defaults to kem.
      type(marss_dfa_spec) :: spec

      spec%ntrends = ntrends
      if (present(demean)) spec%demean = demean
      if (present(z_score)) spec%z_score = z_score
      call marss_dfa_fit_spec(data, spec, max_iter, tol, fit, constraints, info, method, covariates, mstep_iter)
   end subroutine marss_dfa_fit

   subroutine marss_dfa_fit_spec(data, spec, max_iter, tol, fit, constraints, info, method, covariates, mstep_iter)
      real(dp), intent(in) :: data(:, :) !! Observation matrix with variables by time; NaNs denote missing observations.
      type(marss_dfa_spec), intent(in) :: spec !! DFA model and preprocessing specification used to build the fit.
      integer, intent(in) :: max_iter !! Maximum generalized-EM or BFGS iterations, which must be positive.
      real(dp), intent(in) :: tol !! Relative convergence tolerance, which must be positive.
      type(marss_fit_result), intent(out) :: fit !! Fitted DFA model, free beta coordinates, Kalman result, and likelihood.
      type(marss_constraints), intent(out) :: constraints !! DFA affine constraints used for the fit.
      integer, intent(out) :: info !! Zero on success; otherwise model-construction or fitting status.
      character(len=*), intent(in), optional :: method !! Fitting method: kem/em or bfgs/optim; defaults to kem.
      real(dp), intent(in), optional :: covariates(:, :) !! Optional complete observation covariates with predictors by time.
      integer, intent(in), optional :: mstep_iter !! Maximum inner generalized-EM BFGS iterations, defaulting internally.
      type(marss_model) :: model
      character(len=32) :: fit_method
      integer :: build_info

      if (max_iter < 1 .or. tol <= 0.0_dp) then
         info = 1
         fit%info = 1
         return
      end if
      call marss_dfa_build(data, spec, model, constraints, build_info, covariates)
      if (build_info /= 0) then
         info = 100 + build_info
         fit%info = info
         return
      end if
      fit_method = 'kem'
      if (present(method)) fit_method = canonical_spec(method)
      select case (trim(fit_method))
      case ('kem', 'em')
         if (present(mstep_iter)) then
            call marss_kem_linear(model, constraints, max_iter, tol, fit, mstep_iter)
         else
            call marss_kem_linear(model, constraints, max_iter, tol, fit)
         end if
      case ('bfgs', 'optim')
         call marss_optim_linear(model, constraints, max_iter, tol, fit)
      case default
         info = 2
         fit%info = 2
         return
      end select
      info = fit%info
   end subroutine marss_dfa_fit_spec

end module marss_dfa_mod
