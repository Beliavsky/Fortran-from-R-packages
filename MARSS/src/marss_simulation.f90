! SPDX-License-Identifier: GPL-2.0-only
module marss_simulation_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use marss_kinds, only : dp, i8
   use marss_types, only : marss_model, marss_simulation
   use marss_model_ops, only : marss_model_valid
   use marss_parameters, only : marss_b_at, marss_u_at, marss_q_at, marss_z_at, marss_a_at, marss_r_at
   use marss_parameters, only : marss_v0_effective
   use marss_parameters, only : marss_has_time_varying
   use marss_random, only : marss_rng, rng_seed, sample_mvn
   implicit none
   private
   public :: marss_simulate

contains

   subroutine marss_simulate(model, tsteps, nsim, seed, simulation, info, missing_mask, missing_mask_by_sim)
      type(marss_model), intent(in) :: model !! MARSS parameterization; optional time-indexed blocks are honored at each time.
      integer, intent(in) :: tsteps !! Number of observation times; must equal model data length when any block is time-varying.
      integer, intent(in) :: nsim !! Number of independent series to generate, must be positive.
      integer(i8), intent(in) :: seed !! Deterministic seed for the package-local pseudo-random generator.
      type(marss_simulation), intent(out) :: simulation !! Generated state and observation arrays with third dimension nsim.
      integer, intent(out) :: info !! Zero on success, nonzero for invalid inputs, mask shapes, or covariance failures.
      logical, intent(in), optional :: missing_mask(:, :) !! 2-D mask repeated across simulations; true becomes NaN.
      logical, intent(in), optional :: missing_mask_by_sim(:, :, :) !! Observation-by-time-by-simulation mask; true becomes NaN.
      type(marss_rng) :: rng
      real(dp), allocatable :: obs_error(:)
      real(dp), allocatable :: proc_error(:)
      real(dp), allocatable :: state(:)
      real(dp), allocatable :: initial(:)
      real(dp) :: at(size(model%a))
      real(dp) :: bt(size(model%b, 1), size(model%b, 2))
      real(dp) :: qt(size(model%q, 1), size(model%q, 2))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: v0eff(size(model%v0, 1), size(model%v0, 2))
      real(dp) :: ut(size(model%u))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      integer :: i
      integer :: m
      integer :: n
      integer :: s
      integer :: t
      real(dp) :: nan_value

      info = 0
      if (.not. marss_model_valid(model) .or. tsteps < 1 .or. nsim < 1) then
         info = 1
         return
      end if
      if (marss_has_time_varying(model) .and. tsteps /= size(model%y, 2)) then
         info = 2
         return
      end if
      n = size(model%z, 1)
      m = size(model%b, 1)
      if (present(missing_mask)) then
         if (size(missing_mask, 1) /= n .or. size(missing_mask, 2) /= tsteps) then
            info = 3
            return
         end if
      end if
      if (present(missing_mask_by_sim)) then
         if (size(missing_mask_by_sim, 1) /= n .or. size(missing_mask_by_sim, 2) /= tsteps .or. &
            size(missing_mask_by_sim, 3) /= nsim) then
            info = 4
            return
         end if
      end if
      if (present(missing_mask) .and. present(missing_mask_by_sim)) then
         info = 5
         return
      end if
      allocate(simulation%states(m, tsteps, nsim), simulation%data(n, tsteps, nsim))
      allocate(state(m), initial(m), proc_error(m), obs_error(n))
      call rng_seed(rng, seed)
      v0eff = marss_v0_effective(model)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)

      do s = 1, nsim
         call sample_mvn(model%x0, v0eff, rng, initial, info)
         if (info /= 0) return
         state = initial
         do t = 1, tsteps
            bt = marss_b_at(model, t)
            ut = marss_u_at(model, t)
            qt = marss_q_at(model, t)
            zt = marss_z_at(model, t)
            at = marss_a_at(model, t)
            rt = marss_r_at(model, t)
            call sample_mvn([(0.0_dp, i = 1, m)], qt, rng, proc_error, info)
            if (info /= 0) return
            call sample_mvn([(0.0_dp, i = 1, n)], rt, rng, obs_error, info)
            if (info /= 0) return
            if (.not. (t == 1 .and. model%tinitx == 1)) then
               state = matmul(bt, state) + ut + proc_error
            end if
            simulation%states(:, t, s) = state
            simulation%data(:, t, s) = matmul(zt, state) + at + obs_error
         end do
         if (present(missing_mask)) then
            where (missing_mask)
               simulation%data(:, :, s) = nan_value
            end where
         else if (present(missing_mask_by_sim)) then
            where (missing_mask_by_sim(:, :, s))
               simulation%data(:, :, s) = nan_value
            end where
         end if
      end do
   end subroutine marss_simulate

end module marss_simulation_mod
