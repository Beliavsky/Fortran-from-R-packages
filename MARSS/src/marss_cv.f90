! SPDX-License-Identifier: GPL-2.0-only
module marss_cv_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_fit_result, marss_cv_result, marss_constraints
   use marss_optim_mod, only : marss_optim
   use marss_em, only : marss_kem
   use marss_constraints_mod, only : marss_constraints_parameter_count, marss_constraint_start_vector
   use marss_constraints_mod, only : marss_constraints_set_start, marss_apply_constraints, marss_optim_linear
   use marss_constrained_em, only : marss_kem_linear
   use marss_parameters, only : marss_z_at, marss_a_at, marss_r_at
   implicit none
   private
   public :: marss_cv

contains

   subroutine marss_cv(model, max_iter, tol, result, info, fold_ids, future_cv, n_future_cv, &
      prediction_interval, estimate_b, estimate_u, estimate_q, estimate_z, estimate_a, estimate_r, &
      estimate_x0, estimate_v0, constraints, free_parameters, fit_method, fold_free_parameters)
      type(marss_model), intent(in) :: model !! Prepared numerical MARSS model whose observation matrix is cross-validated.
      integer, intent(in) :: max_iter !! Maximum fitting iterations used for each fold refit.
      real(dp), intent(in) :: tol !! Positive convergence tolerance used for each fold refit.
      type(marss_cv_result), intent(out) :: result !! Cross-validated predictions, standard errors, and fold identifiers.
      integer, intent(out) :: info !! Zero on success, or a positive code for invalid folds or a failed fold refit.
      integer, intent(in), optional :: fold_ids(:, :) !! User fold IDs matching y, positive integers, ignored for future CV.
      logical, intent(in), optional :: future_cv !! If true, perform expanding-history future cross-validation.
      integer, intent(in), optional :: n_future_cv !! Number of final time slices predicted in future CV, default floor(T/3).
      logical, intent(in), optional :: prediction_interval !! Include observation noise R in SEs when true, default false.
      logical, intent(in), optional :: estimate_b !! Optimize B in unconstrained fold refits, defaulting to true.
      logical, intent(in), optional :: estimate_u !! Optimize U in unconstrained fold refits, defaulting to true.
      logical, intent(in), optional :: estimate_q !! Optimize Q in unconstrained fold refits, defaulting to true.
      logical, intent(in), optional :: estimate_z !! Optimize Z in unconstrained fold refits, defaulting to true.
      logical, intent(in), optional :: estimate_a !! Optimize A in unconstrained fold refits, defaulting to true.
      logical, intent(in), optional :: estimate_r !! Optimize R in unconstrained fold refits, defaulting to true.
      logical, intent(in), optional :: estimate_x0 !! Optimize x0 in unconstrained fold refits, defaulting to true.
      logical, intent(in), optional :: estimate_v0 !! Optimize V0 in unconstrained fold refits, defaulting to false.
      type(marss_constraints), intent(in), optional :: constraints !! Affine constraints preserved in every fold refit.
      real(dp), intent(in), optional :: free_parameters(:) !! Fitted beta coordinates overriding constraint block starts.
      character(len=*), intent(in), optional :: fit_method !! Refit method, "bfgs" or "kem"; defaults to "bfgs".
      real(dp), allocatable, intent(out), optional :: fold_free_parameters(:, :) !! Constrained beta vector for each fold fit.
      type(marss_constraints) :: fit_constraints
      type(marss_model) :: base_model
      type(marss_model) :: training
      type(marss_fit_result) :: fit
      real(dp), allocatable :: beta0(:)
      real(dp) :: nan
      real(dp) :: variance
      real(dp) :: at(size(model%a))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      character(len=:), allocatable :: mode
      integer :: apply_info
      integer :: i
      integer :: k
      integer :: n
      integer :: nfold
      integer :: nf
      integer :: p
      integer :: start_info
      integer :: t
      integer :: tt
      logical :: future
      logical :: include_r

      info = 0
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         info = 1
         return
      end if
      n = size(model%y, 1)
      tt = size(model%y, 2)
      if (n < 1 .or. tt < 1) then
         info = 2
         return
      end if
      mode = 'bfgs'
      if (present(fit_method)) mode = lowercase_ascii(trim(adjustl(fit_method)))
      if (mode /= 'bfgs' .and. mode /= 'kem') then
         info = 6
         return
      end if

      base_model = model
      p = 0
      if (present(constraints)) then
         p = marss_constraints_parameter_count(constraints)
         if (present(free_parameters)) then
            if (size(free_parameters) /= p) then
               info = 7
               return
            end if
            allocate(beta0(p))
            beta0 = free_parameters
         else
            call marss_constraint_start_vector(constraints, beta0)
         end if
         call marss_apply_constraints(model, constraints, beta0, base_model, apply_info)
         if (apply_info /= 0) then
            info = 8
            return
         end if
         call marss_constraints_set_start(constraints, beta0, fit_constraints, start_info)
         if (start_info /= 0) then
            info = 9
            return
         end if
      else if (present(free_parameters)) then
         if (size(free_parameters) /= 0) then
            info = 10
            return
         end if
      end if

      future = .false.
      if (present(future_cv)) future = future_cv
      include_r = .false.
      if (present(prediction_interval)) include_r = prediction_interval
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(result%prediction(n, tt), result%se(n, tt), result%fold_id(n, tt))
      result%prediction = nan
      result%se = nan

      if (future) then
         nf = max(1, tt / 3)
         if (present(n_future_cv)) nf = n_future_cv
         if (nf < 1 .or. nf > tt) then
            info = 3
            return
         end if
         do t = 1, tt
            result%fold_id(:, t) = tt - t + 1
         end do
         nfold = nf
      else if (present(fold_ids)) then
         if (any(shape(fold_ids) /= [n, tt])) then
            info = 4
            return
         end if
         if (any(fold_ids < 1)) then
            info = 5
            return
         end if
         result%fold_id = fold_ids
         nfold = maxval(fold_ids)
      else
         nfold = min(10, n * tt)
         do t = 1, tt
            do i = 1, n
               result%fold_id(i, t) = 1 + modulo((t - 1) * n + i - 1, nfold)
            end do
         end do
      end if

      if (present(fold_free_parameters)) then
         if (present(constraints)) then
            allocate(fold_free_parameters(p, nfold))
            fold_free_parameters = nan
         else
            allocate(fold_free_parameters(0, nfold))
         end if
      end if

      do k = 1, nfold
         training = base_model
         if (future) then
            training%y(:, tt - k + 1:tt) = nan
         else
            do t = 1, tt
               do i = 1, n
                  if (result%fold_id(i, t) == k) training%y(i, t) = nan
               end do
            end do
         end if
         if (present(constraints)) then
            if (mode == 'kem') then
               call marss_kem_linear(training, fit_constraints, max_iter, tol, fit)
            else
               call marss_optim_linear(training, fit_constraints, max_iter, tol, fit)
            end if
         else
            if (mode == 'kem') then
               call marss_kem(training, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, estimate_z, &
                  estimate_a, estimate_r, estimate_x0, estimate_v0)
            else
               call marss_optim(training, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, estimate_z, &
                  estimate_a, estimate_r, estimate_x0, estimate_v0)
            end if
         end if
         if (fit%info /= 0) then
            info = 1000 + k
            return
         end if
         if (present(fold_free_parameters) .and. present(constraints)) then
            if (.not. allocated(fit%free_parameters) .or. size(fit%free_parameters) /= p) then
               info = 2000 + k
               return
            end if
            fold_free_parameters(:, k) = fit%free_parameters
         end if
         do t = 1, tt
            zt = marss_z_at(fit%model, t)
            at = marss_a_at(fit%model, t)
            rt = marss_r_at(fit%model, t)
            do i = 1, n
               if (result%fold_id(i, t) /= k) cycle
               if (ieee_is_nan(model%y(i, t))) cycle
               result%prediction(i, t) = dot_product(zt(i, :), fit%kf%x_smooth(:, t)) + at(i)
               variance = dot_product(zt(i, :), matmul(fit%kf%p_smooth(:, :, t), zt(i, :)))
               if (include_r) variance = variance + rt(i, i)
               result%se(i, t) = sqrt(max(variance, 0.0_dp))
            end do
         end do
      end do
   end subroutine marss_cv

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

end module marss_cv_mod
