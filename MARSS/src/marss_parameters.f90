! SPDX-License-Identifier: GPL-2.0-only
module marss_parameters
   use marss_kinds, only : dp
   use marss_types, only : marss_model
   implicit none
   private
   public :: marss_b_at
   public :: marss_u_at
   public :: marss_c_effect_at
   public :: marss_q_at
   public :: marss_z_at
   public :: marss_a_at
   public :: marss_d_effect_at
   public :: marss_r_at
   public :: marss_v0_effective
   public :: marss_has_time_varying
   public :: marss_time_shape_valid

contains

   pure function marss_b_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying the state-transition matrix or its time-indexed replacement.
      integer, intent(in) :: t !! One-based observation time whose transition matrix is requested.
      real(dp) :: value(size(model%b, 1), size(model%b, 2))

      if (allocated(model%b_t)) then
         value = model%b_t(:, :, t)
      else
         value = model%b
      end if
   end function marss_b_at

   pure function marss_u_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying U plus any fixed C(t)c(t) state-covariate contribution.
      integer, intent(in) :: t !! One-based observation time whose effective state intercept is requested.
      real(dp) :: value(size(model%u))

      if (allocated(model%u_t)) then
         value = model%u_t(:, t)
      else
         value = model%u
      end if
      value = value + marss_c_effect_at(model, t)
   end function marss_u_at

   pure function marss_c_effect_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying optional C(t) coefficients and c(t) state covariates.
      integer, intent(in) :: t !! One-based observation time whose C(t)c(t) contribution is requested.
      real(dp) :: value(size(model%u))

      value = 0.0_dp
      if (allocated(model%c_coef) .and. allocated(model%state_covariates)) then
         if (allocated(model%c_coef_t)) then
            value = matmul(model%c_coef_t(:, :, t), model%state_covariates(:, t))
         else
            value = matmul(model%c_coef, model%state_covariates(:, t))
         end if
      end if
   end function marss_c_effect_at

   pure function marss_q_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying effective G(t)Q(t)G(t)' or direct process covariance.
      integer, intent(in) :: t !! One-based observation time whose effective process covariance is requested.
      real(dp) :: value(size(model%q, 1), size(model%q, 2))

      if (allocated(model%g) .and. allocated(model%q_noise)) then
         if (allocated(model%g_t) .and. allocated(model%q_noise_t)) then
            value = matmul(matmul(model%g_t(:, :, t), model%q_noise_t(:, :, t)), transpose(model%g_t(:, :, t)))
         else if (allocated(model%g_t)) then
            value = matmul(matmul(model%g_t(:, :, t), model%q_noise), transpose(model%g_t(:, :, t)))
         else if (allocated(model%q_noise_t)) then
            value = matmul(matmul(model%g, model%q_noise_t(:, :, t)), transpose(model%g))
         else
            value = matmul(matmul(model%g, model%q_noise), transpose(model%g))
         end if
      else if (allocated(model%q_t)) then
         value = model%q_t(:, :, t)
      else
         value = model%q
      end if
   end function marss_q_at

   pure function marss_z_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying the observation loading matrix or its time-indexed replacement.
      integer, intent(in) :: t !! One-based observation time whose observation loading matrix is requested.
      real(dp) :: value(size(model%z, 1), size(model%z, 2))

      if (allocated(model%z_t)) then
         value = model%z_t(:, :, t)
      else
         value = model%z
      end if
   end function marss_z_at

   pure function marss_a_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying A plus any fixed D(t)d(t) observation-covariate contribution.
      integer, intent(in) :: t !! One-based observation time whose effective observation intercept is requested.
      real(dp) :: value(size(model%a))

      if (allocated(model%a_t)) then
         value = model%a_t(:, t)
      else
         value = model%a
      end if
      value = value + marss_d_effect_at(model, t)
   end function marss_a_at

   pure function marss_d_effect_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying optional D(t) coefficients and d(t) observation covariates.
      integer, intent(in) :: t !! One-based observation time whose D(t)d(t) contribution is requested.
      real(dp) :: value(size(model%a))

      value = 0.0_dp
      if (allocated(model%d_coef) .and. allocated(model%obs_covariates)) then
         if (allocated(model%d_coef_t)) then
            value = matmul(model%d_coef_t(:, :, t), model%obs_covariates(:, t))
         else
            value = matmul(model%d_coef, model%obs_covariates(:, t))
         end if
      end if
   end function marss_d_effect_at

   pure function marss_r_at(model, t) result(value)
      type(marss_model), intent(in) :: model !! Model supplying effective H(t)R(t)H(t)' or direct observation covariance.
      integer, intent(in) :: t !! One-based observation time whose effective observation covariance is requested.
      real(dp) :: value(size(model%r, 1), size(model%r, 2))

      if (allocated(model%h) .and. allocated(model%r_noise)) then
         if (allocated(model%h_t) .and. allocated(model%r_noise_t)) then
            value = matmul(matmul(model%h_t(:, :, t), model%r_noise_t(:, :, t)), transpose(model%h_t(:, :, t)))
         else if (allocated(model%h_t)) then
            value = matmul(matmul(model%h_t(:, :, t), model%r_noise), transpose(model%h_t(:, :, t)))
         else if (allocated(model%r_noise_t)) then
            value = matmul(matmul(model%h, model%r_noise_t(:, :, t)), transpose(model%h))
         else
            value = matmul(matmul(model%h, model%r_noise), transpose(model%h))
         end if
      else if (allocated(model%r_t)) then
         value = model%r_t(:, :, t)
      else
         value = model%r
      end if
   end function marss_r_at

   pure function marss_v0_effective(model) result(value)
      type(marss_model), intent(in) :: model !! Model supplying effective L V0 L' or the direct initial-state covariance.
      real(dp) :: value(size(model%v0, 1), size(model%v0, 2))

      if (allocated(model%l) .and. allocated(model%v0_noise)) then
         value = matmul(matmul(model%l, model%v0_noise), transpose(model%l))
      else
         value = model%v0
      end if
   end function marss_v0_effective

   pure logical function marss_has_time_varying(model) result(has_time_varying)
      type(marss_model), intent(in) :: model !! Model tested for any time-indexed numerical parameter or loading block.

      has_time_varying = allocated(model%b_t) .or. allocated(model%u_t) .or. allocated(model%q_t) .or. &
         allocated(model%z_t) .or. allocated(model%a_t) .or. allocated(model%r_t) .or. &
         allocated(model%c_coef_t) .or. allocated(model%d_coef_t) .or. allocated(model%g_t) .or. &
         allocated(model%h_t) .or. allocated(model%q_noise_t) .or. allocated(model%r_noise_t)
   end function marss_has_time_varying

   pure logical function marss_time_shape_valid(model, ntime) result(ok)
      type(marss_model), intent(in) :: model !! Model whose optional time-indexed parameter and covariate dimensions are validated.
      integer, intent(in) :: ntime !! Required number of time slices in every allocated time-indexed block.
      integer :: g1
      integer :: h1
      integer :: m
      integer :: n
      integer :: p

      ok = .false.
      if (ntime < 1) return
      if (.not. allocated(model%b)) return
      if (.not. allocated(model%u)) return
      if (.not. allocated(model%q)) return
      if (.not. allocated(model%z)) return
      if (.not. allocated(model%a)) return
      if (.not. allocated(model%r)) return
      m = size(model%b, 1)
      n = size(model%z, 1)
      if (allocated(model%b_t)) then
         if (any(shape(model%b_t) /= [m, m, ntime])) return
      end if
      if (allocated(model%u_t)) then
         if (any(shape(model%u_t) /= [m, ntime])) return
      end if
      if (allocated(model%q_t)) then
         if (any(shape(model%q_t) /= [m, m, ntime])) return
      end if
      if (allocated(model%z_t)) then
         if (any(shape(model%z_t) /= [n, m, ntime])) return
      end if
      if (allocated(model%a_t)) then
         if (any(shape(model%a_t) /= [n, ntime])) return
      end if
      if (allocated(model%r_t)) then
         if (any(shape(model%r_t) /= [n, n, ntime])) return
      end if
      if (allocated(model%state_covariates) .neqv. allocated(model%c_coef)) return
      if (allocated(model%state_covariates)) then
         p = size(model%state_covariates, 1)
         if (size(model%state_covariates, 2) /= ntime) return
         if (any(shape(model%c_coef) /= [m, p])) return
         if (allocated(model%c_coef_t)) then
            if (any(shape(model%c_coef_t) /= [m, p, ntime])) return
         end if
      else if (allocated(model%c_coef_t)) then
         return
      end if
      if (allocated(model%obs_covariates) .neqv. allocated(model%d_coef)) return
      if (allocated(model%obs_covariates)) then
         p = size(model%obs_covariates, 1)
         if (size(model%obs_covariates, 2) /= ntime) return
         if (any(shape(model%d_coef) /= [n, p])) return
         if (allocated(model%d_coef_t)) then
            if (any(shape(model%d_coef_t) /= [n, p, ntime])) return
         end if
      else if (allocated(model%d_coef_t)) then
         return
      end if
      if (allocated(model%g) .neqv. allocated(model%q_noise)) return
      if (allocated(model%g)) then
         g1 = size(model%g, 2)
         if (size(model%g, 1) /= m) return
         if (any(shape(model%q_noise) /= [g1, g1])) return
         if (allocated(model%g_t)) then
            if (any(shape(model%g_t) /= [m, g1, ntime])) return
         end if
         if (allocated(model%q_noise_t)) then
            if (any(shape(model%q_noise_t) /= [g1, g1, ntime])) return
         end if
      else if (allocated(model%g_t) .or. allocated(model%q_noise_t)) then
         return
      end if
      if (allocated(model%h) .neqv. allocated(model%r_noise)) return
      if (allocated(model%h)) then
         h1 = size(model%h, 2)
         if (size(model%h, 1) /= n) return
         if (any(shape(model%r_noise) /= [h1, h1])) return
         if (allocated(model%h_t)) then
            if (any(shape(model%h_t) /= [n, h1, ntime])) return
         end if
         if (allocated(model%r_noise_t)) then
            if (any(shape(model%r_noise_t) /= [h1, h1, ntime])) return
         end if
      else if (allocated(model%h_t) .or. allocated(model%r_noise_t)) then
         return
      end if
      ok = .true.
   end function marss_time_shape_valid

end module marss_parameters
