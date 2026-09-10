! SPDX-License-Identifier: GPL-2.0-only
module marss_workflow
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_model_spec, marss_dfa_spec, marss_constraints, marss_fit_result
   use marss_builder, only : marss_build
   use marss_constrained_em, only : marss_kem_linear
   use marss_constraints_mod, only : marss_optim_linear
   use marss_dfa_mod, only : marss_dfa_fit_spec
   implicit none
   private
   public :: marss_from_data

   interface marss_from_data
      module procedure marss_from_data_model
      module procedure marss_from_data_dfa
   end interface marss_from_data

contains

   subroutine marss_from_data_model(y, spec, max_iter, tol, fit, constraints, method, info, state_covariates, obs_covariates)
      real(dp), intent(in) :: y(:, :) !! Observation matrix, variables by time; NaNs denote missing observations.
      type(marss_model_spec), intent(in) :: spec !! High-level marss/marxss shortcut specification used to build the model.
      integer, intent(in) :: max_iter !! Maximum generalized-EM or BFGS iterations; must be positive.
      real(dp), intent(in) :: tol !! Relative convergence tolerance; must be positive.
      type(marss_fit_result), intent(out) :: fit !! Fitted model, Kalman output, free beta coordinates, and diagnostics.
      type(marss_constraints), intent(out) :: constraints !! Constructed affine f+D*beta representation used for fitting.
      character(len=*), intent(in), optional :: method !! Fitting method: "kem"/"em" or "bfgs"/"optim"; defaults to "kem".
      integer, intent(out) :: info !! Zero when model construction and fitting start successfully; otherwise a setup error code.
      real(dp), intent(in), optional :: state_covariates(:, :) !! Optional c(t) rows passed to form=marxss-style construction.
      real(dp), intent(in), optional :: obs_covariates(:, :) !! Optional d(t) rows passed to form=marxss-style construction.
      type(marss_model) :: start_model
      character(len=16) :: fit_method

      info = 0
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         info = 1
         return
      end if
      fit_method = "kem"
      if (present(method)) fit_method = canonical_method(method)
      call marss_build(y, spec, start_model, constraints, info, state_covariates, obs_covariates)
      if (info /= 0) then
         info = 1000 + info
         return
      end if
      select case (trim(fit_method))
      case ("kem", "em")
         call marss_kem_linear(start_model, constraints, max_iter, tol, fit)
      case ("bfgs", "optim")
         call marss_optim_linear(start_model, constraints, max_iter, tol, fit)
      case default
         info = 2
         return
      end select
      if (fit%info /= 0) info = 2000 + fit%info
   end subroutine marss_from_data_model

   subroutine marss_from_data_dfa(y, spec, max_iter, tol, fit, constraints, method, info, &
      state_covariates, obs_covariates)
      real(dp), intent(in) :: y(:, :) !! Observation matrix with variables by time; NaNs denote missing observations.
      type(marss_dfa_spec), intent(in) :: spec !! High-level DFA form specification and preprocessing choices.
      integer, intent(in) :: max_iter !! Maximum generalized-EM or BFGS iterations; must be positive.
      real(dp), intent(in) :: tol !! Relative convergence tolerance; must be positive.
      type(marss_fit_result), intent(out) :: fit !! Fitted DFA model, Kalman output, free beta coordinates, and diagnostics.
      type(marss_constraints), intent(out) :: constraints !! Constructed DFA affine constraints used for fitting.
      character(len=*), intent(in), optional :: method !! Fitting method: kem/em or bfgs/optim; defaults to kem.
      integer, intent(out) :: info !! Zero when DFA construction/fitting succeeds; otherwise a setup or fitting status.
      real(dp), intent(in), optional :: state_covariates(:, :) !! Unsupported DFA state covariates; presence returns an error.
      real(dp), intent(in), optional :: obs_covariates(:, :) !! Optional DFA observation covariates d(t), complete by construction.

      if (present(state_covariates)) then
         info = 3
         fit%info = 3
         return
      end if
      if (present(obs_covariates)) then
         if (present(method)) then
            call marss_dfa_fit_spec(y, spec, max_iter, tol, fit, constraints, info, method, obs_covariates)
         else
            call marss_dfa_fit_spec(y, spec, max_iter, tol, fit, constraints, info, covariates=obs_covariates)
         end if
      else
         if (present(method)) then
            call marss_dfa_fit_spec(y, spec, max_iter, tol, fit, constraints, info, method)
         else
            call marss_dfa_fit_spec(y, spec, max_iter, tol, fit, constraints, info)
         end if
      end if
   end subroutine marss_from_data_dfa

   pure function canonical_method(text) result(value)
      character(len=*), intent(in) :: text !! Fitting-method text normalized to lowercase for dispatch.
      character(len=16) :: value
      integer :: code
      integer :: i

      value = ""
      value = adjustl(trim(text))
      do i = 1, len_trim(value)
         code = iachar(value(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) value(i:i) = achar(code + 32)
      end do
   end function canonical_method

end module marss_workflow
