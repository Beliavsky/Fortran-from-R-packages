! SPDX-License-Identifier: MIT
module invgamstochvol_model
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use, intrinsic :: iso_fortran_env, only : int64
   use invgamstochvol_kinds, only : dp
   use invgamstochvol_status, only : invgam_success, invgam_invalid_argument, &
      invgam_nonfinite_input, invgam_numerical_failure
   use invgamstochvol_special, only : ourgeo, build_log_factorials, &
      hypergeo_from_tables, log_sum_exp
   use invgamstochvol_rng, only : rng_state
   implicit none
   private

   real(dp), parameter :: log_two = log(2.0_dp)
   real(dp), parameter :: log_two_pi = log(2.0_dp * acos(-1.0_dp))

   type, public :: invgam_likelihood_result
      real(dp) :: total_loglik = 0.0_dp
      real(dp), allocatable :: loglik(:)
      real(dp), allocatable :: all_st(:)
      real(dp), allocatable :: all_ctil(:, :)
      real(dp), allocatable :: alogfac(:, :)
      real(dp), allocatable :: alogfac2(:)
      real(dp), allocatable :: alfac(:)
      integer :: nit = 0
      integer :: niter = 0
      integer :: status = invgam_success
      character(len=160) :: message = ''
   end type invgam_likelihood_result

   interface draw_k0
      module procedure draw_k0_result
      module procedure draw_k0_arrays
   end interface draw_k0

   public :: lik_clo, draw_k0

contains

   subroutine clear_result(result)
      type(invgam_likelihood_result), intent(inout) :: result

      result%total_loglik = 0.0_dp
      result%nit = 0
      result%niter = 0
      result%status = invgam_success
      result%message = ''
      if (allocated(result%loglik)) deallocate(result%loglik)
      if (allocated(result%all_st)) deallocate(result%all_st)
      if (allocated(result%all_ctil)) deallocate(result%all_ctil)
      if (allocated(result%alogfac)) deallocate(result%alogfac)
      if (allocated(result%alogfac2)) deallocate(result%alogfac2)
      if (allocated(result%alfac)) deallocate(result%alfac)
   end subroutine clear_result

   subroutine fail_result(result, status, message)
      type(invgam_likelihood_result), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message

      result%status = status
      result%message = message
   end subroutine fail_result

   subroutine lik_clo(residuals, b2, nu, rho, result, nit, niter, nproc, nproc2)
      real(dp), intent(in) :: residuals(:)
      real(dp), intent(in) :: b2, nu, rho
      type(invgam_likelihood_result), intent(out) :: result
      integer, intent(in), optional :: nit, niter, nproc, nproc2
      integer :: trunc_outer, trunc_hyper
      integer :: t_count, h, hold, tt, stat
      real(dp) :: v_inv, v_inv_til, et, useme, l0, bc
      real(dp) :: delta_h2, z, c2, normsum, s3, delta_h3, c3
      real(dp) :: z3, ccc, useful_log, p_st, st, delta_ht, zt
      real(dp) :: log_ct, log_lik_sum, log_normalization, geo_value
      real(dp), allocatable :: old_ctil(:), new_ctil(:)
      real(dp), allocatable :: all_vinv(:), all_zt(:), all_delta(:)
      real(dp), allocatable :: all_geo(:, :), work(:), alct(:), allik(:)

      call clear_result(result)
      trunc_outer = 200
      trunc_hyper = 200
      if (present(nit)) trunc_outer = nit
      if (present(niter)) trunc_hyper = niter

      if (present(nproc)) then
         if (nproc < 1) then
            call fail_result(result, invgam_invalid_argument, 'nproc must be positive')
            return
         end if
      end if
      if (present(nproc2)) then
         if (nproc2 < 1) then
            call fail_result(result, invgam_invalid_argument, 'nproc2 must be positive')
            return
         end if
      end if

      t_count = size(residuals)
      if (t_count < 3) then
         call fail_result(result, invgam_invalid_argument, &
            'lik_clo requires at least three residuals')
         return
      end if
      if (trunc_outer < 0 .or. trunc_hyper < 1) then
         call fail_result(result, invgam_invalid_argument, &
            'truncation orders must satisfy nit >= 0 and niter >= 1')
         return
      end if
      if (b2 <= 0.0_dp .or. nu <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         call fail_result(result, invgam_invalid_argument, &
            'require b2 > 0, nu > 0, and abs(rho) < 1')
         return
      end if
      if (.not. ieee_is_finite(b2) .or. .not. ieee_is_finite(nu) .or. &
          .not. ieee_is_finite(rho) .or. .not. all(ieee_is_finite(residuals))) then
         call fail_result(result, invgam_nonfinite_input, &
            'all model inputs must be finite')
         return
      end if

      result%nit = trunc_outer
      result%niter = trunc_hyper
      allocate(result%loglik(0:t_count - 1))
      allocate(result%all_st(0:t_count))
      allocate(result%all_ctil(0:t_count - 1, 0:trunc_outer))
      result%loglik = 0.0_dp
      result%all_st = 0.0_dp
      result%all_ctil = 0.0_dp

      call build_log_factorials(nu, trunc_outer, trunc_hyper, &
         result%alogfac, result%alogfac2, result%alfac, stat)
      if (stat /= invgam_success) then
         call fail_result(result, stat, 'failed to construct rising-factorial tables')
         return
      end if

      allocate(old_ctil(0:trunc_outer), new_ctil(0:trunc_outer))
      allocate(all_vinv(0:t_count), all_zt(0:t_count), all_delta(0:t_count))
      allocate(all_geo(0:t_count - 1, 0:trunc_outer))
      allocate(work(0:trunc_outer), alct(0:trunc_outer), allik(0:trunc_outer))
      old_ctil = 0.0_dp
      new_ctil = 0.0_dp
      all_vinv = 0.0_dp
      all_zt = 0.0_dp
      all_delta = 0.0_dp
      all_geo = 0.0_dp

      v_inv = 1.0_dp - rho * rho
      et = residuals(1)
      useme = b2 * et * et
      l0 = -0.5_dp * (nu + 1.0_dp) * log(0.5_dp * (v_inv + useme))
      l0 = l0 + 0.5_dp * nu * log(0.5_dp * v_inv)
      l0 = l0 + log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu)
      l0 = l0 - 0.5_dp * log_two_pi + 0.5_dp * log(b2)
      result%loglik(0) = l0

      v_inv_til = v_inv + b2 * et * et
      bc = 0.5_dp * rho * rho / (rho * rho + v_inv_til)
      et = residuals(2)
      delta_h2 = rho * rho / (rho * rho + v_inv_til)
      z = delta_h2 / (b2 * et * et + 1.0_dp)
      ccc = log_gamma(0.5_dp * (nu + 1.0_dp))
      geo_value = ourgeo(-0.5_dp, -0.5_dp, 0.5_dp * nu, z, trunc_hyper, stat)
      if (stat /= invgam_success .or. geo_value <= 0.0_dp) then
         call fail_result(result, invgam_numerical_failure, &
            'hypergeometric evaluation failed at the second observation')
         return
      end if
      c2 = (1.0_dp - z)**(-0.5_dp * (nu + 2.0_dp)) * geo_value
      normsum = 2.0_dp**(0.5_dp * nu) * gamma(0.5_dp * nu)
      normsum = normsum * (1.0_dp - delta_h2)**(-0.5_dp * (nu + 1.0_dp))
      l0 = -0.5_dp * log_two_pi + 0.5_dp * log(b2)
      l0 = l0 + 0.5_dp * (nu + 1.0_dp) * log_two + ccc
      l0 = l0 - 0.5_dp * (nu + 1.0_dp) * log(b2 * et * et + 1.0_dp)
      result%loglik(1) = l0 + log(c2) - log(normsum)

      s3 = 1.0_dp / (b2 * et * et + 1.0_dp + rho * rho)
      delta_h3 = s3 * (rho * rho / (v_inv_til + rho * rho))
      delta_h3 = delta_h3 / (1.0_dp - s3 * rho * rho)
      geo_value = ourgeo(-0.5_dp, -0.5_dp, 0.5_dp * nu, delta_h3, &
         trunc_hyper, stat)
      if (stat /= invgam_success .or. geo_value <= 0.0_dp) then
         call fail_result(result, invgam_numerical_failure, &
            'hypergeometric evaluation failed at the third observation')
         return
      end if
      c3 = (1.0_dp - delta_h3)**(-0.5_dp * (nu + 2.0_dp)) * geo_value
      c3 = c3 * gamma(0.5_dp * (nu + 1.0_dp))
      c3 = c3 * (1.0_dp - s3 * rho * rho)**(-0.5_dp * (nu + 1.0_dp))
      c3 = c3 * (2.0_dp * s3)**(0.5_dp * (nu + 1.0_dp))

      et = residuals(3)
      v_inv_til = 1.0_dp + b2 * et * et
      z3 = s3 * rho * rho / v_inv_til
      l0 = -0.5_dp * log_two_pi + 0.5_dp * log(b2) - log(c3)
      useful_log = -0.5_dp * (nu + 1.0_dp) * log(v_inv_til)
      useful_log = useful_log + log_gamma(0.5_dp * (nu + 1.0_dp))
      useful_log = useful_log - log_gamma(0.5_dp * nu) + 0.5_dp * log_two

      do h = 0, trunc_outer
         geo_value = hypergeo_from_tables(h, result%alogfac, result%alogfac2, &
            result%alfac, z3, trunc_hyper, stat)
         if (stat /= invgam_success .or. geo_value <= 0.0_dp) then
            call fail_result(result, invgam_numerical_failure, &
               'table-based hypergeometric evaluation failed')
            return
         end if
         old_ctil(h) = result%alogfac(0, h) - result%alogfac2(h) - result%alfac(h)
         if (h > 0) then
            if (bc <= 0.0_dp) then
               old_ctil(h) = -huge(1.0_dp)
            else
               old_ctil(h) = old_ctil(h) + real(h, dp) * log(bc)
            end if
         end if
         work(h) = old_ctil(h) + ccc + result%alogfac(0, h)
         work(h) = work(h) + 0.5_dp * (nu + 1.0_dp + 2.0_dp * real(h, dp)) &
            * log(2.0_dp * s3) + log(geo_value) + useful_log
      end do
      log_lik_sum = log_sum_exp(work, stat)
      if (stat /= invgam_success) then
         call fail_result(result, invgam_numerical_failure, &
            'failed to combine third-observation likelihood terms')
         return
      end if
      result%loglik(2) = l0 + log_lik_sum
      result%all_ctil(1, :) = old_ctil

      st = s3
      result%all_st(2) = st
      do tt = 3, t_count - 1
         st = 1.0_dp / (b2 * et * et + 1.0_dp + rho * rho)
         delta_ht = st * (rho * rho / (v_inv_til + rho * rho))
         delta_ht = delta_ht / (1.0_dp - st * rho * rho)
         et = residuals(tt + 1)
         v_inv_til = 1.0_dp + b2 * et * et
         zt = st * rho * rho / v_inv_til
         result%all_st(tt) = st
         all_delta(tt) = delta_ht
         all_vinv(tt) = v_inv_til
         all_zt(tt) = zt
      end do

      do tt = 3, t_count - 1
         do h = 0, trunc_outer
            geo_value = hypergeo_from_tables(h, result%alogfac, result%alogfac2, &
               result%alfac, all_zt(tt), trunc_hyper, stat)
            if (stat /= invgam_success .or. geo_value <= 0.0_dp) then
               call fail_result(result, invgam_numerical_failure, &
                  'failed to precompute a hypergeometric likelihood term')
               return
            end if
            all_geo(tt, h) = log(geo_value)
         end do
      end do

      l0 = -0.5_dp * log_two_pi + 0.5_dp * log(b2)
      do tt = 3, t_count
         p_st = result%all_st(tt - 1)
         if (tt < t_count) then
            st = result%all_st(tt)
            v_inv_til = all_vinv(tt)
            zt = all_zt(tt)
         end if

         do h = 0, trunc_outer
            do hold = 0, trunc_outer
               work(hold) = old_ctil(hold) + ccc + result%alogfac(0, hold)
               work(hold) = work(hold) + result%alogfac(hold, h)
               work(hold) = work(hold) + 0.5_dp * &
                  (nu + 1.0_dp + 2.0_dp * real(hold, dp)) * log(2.0_dp * p_st)
            end do
            new_ctil(h) = log_sum_exp(work, stat) - result%alogfac2(h) - result%alfac(h)
            if (stat /= invgam_success) then
               call fail_result(result, invgam_numerical_failure, &
                  'failed to update likelihood recursion coefficients')
               return
            end if
            if (h > 0) then
               if (abs(rho) <= tiny(1.0_dp)) then
                  new_ctil(h) = -huge(1.0_dp)
               else
                  new_ctil(h) = new_ctil(h) + real(h, dp) * &
                     log(0.5_dp * rho * rho * p_st)
               end if
            end if

            if (tt < t_count) then
               alct(h) = new_ctil(h) - 0.5_dp * &
                  (nu + 1.0_dp + 2.0_dp * real(h, dp)) * log(1.0_dp - rho * rho * st)
               alct(h) = alct(h) + ccc + result%alogfac(0, h)
               alct(h) = alct(h) + 0.5_dp * &
                  (nu + 1.0_dp + 2.0_dp * real(h, dp)) * log(2.0_dp * st)

               allik(h) = new_ctil(h) + ccc + result%alogfac(0, h)
               allik(h) = allik(h) + 0.5_dp * &
                  (nu + 1.0_dp + 2.0_dp * real(h, dp)) * log(2.0_dp * st)
               allik(h) = allik(h) + all_geo(tt, h)
            end if
         end do

         if (tt < t_count) then
            log_ct = log_sum_exp(alct, stat)
            if (stat /= invgam_success) then
               call fail_result(result, invgam_numerical_failure, &
                  'failed to normalize recursion coefficients')
               return
            end if
            log_lik_sum = log_sum_exp(allik, stat)
            if (stat /= invgam_success) then
               call fail_result(result, invgam_numerical_failure, &
                  'failed to combine likelihood terms')
               return
            end if
            log_normalization = -0.5_dp * (nu + 1.0_dp) * log(v_inv_til)
            log_normalization = log_normalization + &
               log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu)
            log_normalization = log_normalization + 0.5_dp * log_two
            result%loglik(tt) = l0 + log_lik_sum + log_normalization - log_ct
         end if

         new_ctil = new_ctil - new_ctil(0)
         result%all_ctil(tt - 1, :) = new_ctil
         old_ctil = new_ctil
         new_ctil = 0.0_dp
      end do

      result%all_st(1) = 1.0_dp / (b2 * residuals(1)**2 + 1.0_dp)
      result%all_st(t_count) = 1.0_dp / (b2 * residuals(t_count)**2 + 1.0_dp)
      result%total_loglik = sum(result%loglik)
      if (.not. ieee_is_finite(result%total_loglik)) then
         call fail_result(result, invgam_numerical_failure, &
            'the total likelihood is not finite')
         return
      end if
      result%status = invgam_success
      result%message = 'success'
   end subroutine lik_clo

   subroutine draw_k0_result(result, nu, rho, b2, inverse_volatility, seed, status, nproc2)
      type(invgam_likelihood_result), intent(in) :: result
      real(dp), intent(in) :: nu, rho, b2
      real(dp), allocatable, intent(out) :: inverse_volatility(:)
      integer(int64), intent(in), optional :: seed
      integer, intent(out), optional :: status
      integer, intent(in), optional :: nproc2
      integer :: local_status

      if (result%status /= invgam_success) then
         if (allocated(inverse_volatility)) deallocate(inverse_volatility)
         if (present(status)) status = invgam_invalid_argument
         return
      end if
      call draw_k0_arrays(result%all_st, result%all_ctil, result%alogfac, &
         result%alogfac2, result%alfac, nu, rho, b2, inverse_volatility, &
         seed, local_status, nproc2)
      if (present(status)) status = local_status
   end subroutine draw_k0_result

   subroutine draw_k0_arrays(all_st, all_ctil, alogfac, alogfac2, alfac, &
      nu, rho, b2, inverse_volatility, seed, status, nproc2)
      real(dp), intent(in) :: all_st(0:), all_ctil(0:, 0:)
      real(dp), intent(in) :: alogfac(0:, 0:), alogfac2(0:), alfac(0:)
      real(dp), intent(in) :: nu, rho, b2
      real(dp), allocatable, intent(out) :: inverse_volatility(:)
      integer(int64), intent(in), optional :: seed
      integer, intent(out), optional :: status
      integer, intent(in), optional :: nproc2
      type(rng_state) :: rng
      integer :: t_count, trunc_outer, cualt, h, ht, selected, local_status
      real(dp) :: k_here, freedom, log_rho_term
      real(dp), allocatable :: log_weights(:), weights(:), terms(:)

      local_status = invgam_success
      if (present(nproc2)) then
         if (nproc2 < 1) then
            if (allocated(inverse_volatility)) deallocate(inverse_volatility)
            if (present(status)) status = invgam_invalid_argument
            return
         end if
      end if

      t_count = size(all_ctil, 1)
      trunc_outer = ubound(all_ctil, 2)
      if (t_count < 1 .or. ubound(all_st, 1) < t_count .or. &
          ubound(alogfac, 1) < trunc_outer .or. ubound(alogfac2, 1) < trunc_outer .or. &
          ubound(alfac, 1) < trunc_outer .or. nu <= 0.0_dp .or. b2 <= 0.0_dp .or. &
          abs(rho) >= 1.0_dp) then
         if (allocated(inverse_volatility)) deallocate(inverse_volatility)
         if (present(status)) status = invgam_invalid_argument
         return
      end if
      if (.not. all(ieee_is_finite(all_st)) .or. &
          .not. all(ieee_is_finite(all_ctil)) .or. &
          .not. all(ieee_is_finite(alogfac)) .or. &
          .not. all(ieee_is_finite(alogfac2)) .or. &
          .not. all(ieee_is_finite(alfac))) then
         if (allocated(inverse_volatility)) deallocate(inverse_volatility)
         if (present(status)) status = invgam_nonfinite_input
         return
      end if

      if (present(seed)) then
         call rng%seed(seed)
      else
         call rng%seed(104729_int64)
      end if

      allocate(inverse_volatility(1:t_count))
      allocate(log_weights(0:trunc_outer), weights(0:trunc_outer))
      allocate(terms(0:trunc_outer))
      freedom = 0.5_dp * (nu + 1.0_dp)

      cualt = t_count
      do h = 0, trunc_outer
         log_weights(h) = all_ctil(cualt - 1, h) + alogfac(0, h)
         log_weights(h) = log_weights(h) + 0.5_dp * &
            (nu + 1.0_dp + 2.0_dp * real(h, dp)) * log(2.0_dp * all_st(cualt))
      end do
      call normalize_log_weights(log_weights, weights, local_status)
      if (local_status /= invgam_success) then
         deallocate(inverse_volatility)
         if (present(status)) status = local_status
         return
      end if
      selected = sample_discrete(weights, rng)
      k_here = rng%gamma(freedom + real(selected, dp), 2.0_dp * all_st(cualt))
      inverse_volatility(cualt) = k_here

      if (abs(rho) <= tiny(1.0_dp)) then
         log_rho_term = -huge(1.0_dp)
      else
         log_rho_term = log(0.25_dp * rho * rho)
      end if

      do cualt = t_count - 1, 1, -1
         do h = 0, trunc_outer
            if (cualt > 1) then
               do ht = 0, h
                  terms(ht) = all_ctil(cualt - 1, h - ht) - alogfac2(ht) - alfac(ht)
                  if (ht > 0) then
                     if (abs(rho) <= tiny(1.0_dp)) then
                        terms(ht) = -huge(1.0_dp)
                     else
                        terms(ht) = terms(ht) + real(ht, dp) * &
                           (log(k_here) + log_rho_term)
                     end if
                  end if
               end do
               log_weights(h) = log_sum_exp(terms(0:h), local_status)
            else
               log_weights(h) = -alogfac2(h) - alfac(h)
               if (h > 0) then
                  if (abs(rho) <= tiny(1.0_dp)) then
                     log_weights(h) = -huge(1.0_dp)
                  else
                     log_weights(h) = log_weights(h) + real(h, dp) * &
                        (log(k_here) + log_rho_term)
                  end if
               end if
            end if
            if (local_status /= invgam_success) then
               deallocate(inverse_volatility)
               if (present(status)) status = local_status
               return
            end if
            log_weights(h) = log_weights(h) + alogfac(0, h)
            log_weights(h) = log_weights(h) + 0.5_dp * &
               (nu + 1.0_dp + 2.0_dp * real(h, dp)) * log(2.0_dp * all_st(cualt))
         end do

         call normalize_log_weights(log_weights, weights, local_status)
         if (local_status /= invgam_success) then
            deallocate(inverse_volatility)
            if (present(status)) status = local_status
            return
         end if
         selected = sample_discrete(weights, rng)
         k_here = rng%gamma(freedom + real(selected, dp), 2.0_dp * all_st(cualt))
         inverse_volatility(cualt) = k_here
      end do

      if (present(status)) status = invgam_success
   end subroutine draw_k0_arrays

   subroutine normalize_log_weights(log_weights, weights, status)
      real(dp), intent(in) :: log_weights(0:)
      real(dp), intent(out) :: weights(0:)
      integer, intent(out) :: status
      real(dp) :: normalizer

      normalizer = log_sum_exp(log_weights, status)
      if (status /= invgam_success) return
      weights = exp(log_weights - normalizer)
      if (.not. all(ieee_is_finite(weights)) .or. sum(weights) <= 0.0_dp) then
         status = invgam_numerical_failure
         return
      end if
      weights = weights / sum(weights)
   end subroutine normalize_log_weights

   function sample_discrete(weights, rng) result(selected)
      real(dp), intent(in) :: weights(0:)
      type(rng_state), intent(inout) :: rng
      integer :: selected
      real(dp) :: u, cumulative
      integer :: h

      u = rng%uniform()
      cumulative = 0.0_dp
      selected = ubound(weights, 1)
      do h = 0, ubound(weights, 1)
         cumulative = cumulative + weights(h)
         if (u <= cumulative) then
            selected = h
            exit
         end if
      end do
   end function sample_discrete

end module invgamstochvol_model
