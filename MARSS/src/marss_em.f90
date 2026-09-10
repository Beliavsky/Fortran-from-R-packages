! SPDX-License-Identifier: GPL-2.0-only
module marss_em
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_fit_result, marss_kf_result, marss_hatyt_result
   use marss_kalman, only : marss_kfss
   use marss_analysis, only : marss_hatyt
   use marss_utils, only : identity_matrix
   use marss_parameters, only : marss_b_at, marss_u_at, marss_z_at, marss_a_at
   use marss_parameters, only : marss_c_effect_at, marss_d_effect_at
   use r_linalg, only : inverse_matrix, symmetrize
   implicit none
   private
   public :: marss_kem
   public :: marss_fit
   public :: marss

contains

   pure subroutine marss_kem(start_model, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, &
      estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      type(marss_model), intent(in) :: start_model !! Starting values for a time-invariant unconstrained MARSS model.
      integer, intent(in) :: max_iter !! Maximum number of EM iterations, must be positive.
      real(dp), intent(in) :: tol !! Relative log-likelihood convergence tolerance.
      type(marss_fit_result), intent(out) :: fit !! Fitted model, final smoother output, convergence flag, and iteration count.
      logical, intent(in), optional :: estimate_b !! Estimate the state-transition matrix B, defaults to true.
      logical, intent(in), optional :: estimate_u !! Estimate the state intercept U, defaults to true.
      logical, intent(in), optional :: estimate_q !! Estimate the process covariance Q, defaults to true.
      logical, intent(in), optional :: estimate_z !! Estimate the observation loading matrix Z, defaults to true.
      logical, intent(in), optional :: estimate_a !! Estimate the observation intercept A, defaults to true.
      logical, intent(in), optional :: estimate_r !! Estimate the observation covariance R, defaults to true.
      logical, intent(in), optional :: estimate_x0 !! Estimate the initial-state mean x0, defaults to true.
      logical, intent(in), optional :: estimate_v0 !! Estimate the initial-state covariance V0, defaults to false.
      type(marss_model) :: current
      type(marss_model) :: updated
      type(marss_kf_result) :: kf
      real(dp) :: old_loglik
      logical :: eb
      logical :: eu
      logical :: eq
      logical :: ez
      logical :: ea
      logical :: er
      logical :: ex0
      logical :: ev0
      integer :: iter
      integer :: info

      fit%converged = .false.
      fit%iterations = 0
      fit%info = 0
      current = start_model
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
      if (start_model%diffuse) then
         fit%info = 5
         return
      end if
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         fit%info = 1
         return
      end if
      if ((allocated(current%b_t) .and. eb) .or. (allocated(current%u_t) .and. eu) .or. &
         (allocated(current%q_t) .and. eq) .or. (allocated(current%z_t) .and. ez) .or. &
         (allocated(current%a_t) .and. ea) .or. (allocated(current%r_t) .and. er)) then
         fit%info = 3
         return
      end if
      if ((allocated(current%q_noise) .and. eq) .or. (allocated(current%r_noise) .and. er) .or. &
         (allocated(current%v0_noise) .and. ev0)) then
         fit%info = 4
         return
      end if
      old_loglik = -huge(1.0_dp)
      do iter = 1, max_iter
         call marss_kfss(current, kf)
         if (.not. kf%ok) then
            fit%info = 100 + kf%info
            return
         end if
         if (iter > 1) then
            if (abs(kf%loglik - old_loglik) <= tol * (1.0_dp + abs(old_loglik))) then
               fit%converged = .true.
               exit
            end if
         end if
         old_loglik = kf%loglik
         call em_update(current, kf, updated, info, eb, eu, eq, ez, ea, er, ex0, ev0)
         if (info /= 0) then
            fit%info = 200 + info
            return
         end if
         current = updated
      end do
      fit%iterations = min(iter, max_iter)
      call marss_kfss(current, fit%kf)
      if (.not. fit%kf%ok) then
         fit%info = 300 + fit%kf%info
         return
      end if
      if (.not. fit%converged .and. fit%iterations < max_iter) fit%converged = .true.
      fit%model = current
      fit%loglik = fit%kf%loglik
      fit%info = 0
   end subroutine marss_kem

   pure subroutine marss_fit(start_model, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, &
      estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      type(marss_model), intent(in) :: start_model !! Starting MARSS model for the EM fit method.
      integer, intent(in) :: max_iter !! Maximum EM iterations.
      real(dp), intent(in) :: tol !! Relative convergence tolerance for successive log likelihoods.
      type(marss_fit_result), intent(out) :: fit !! Result corresponding to the computational MARSSfit role.
      logical, intent(in), optional :: estimate_b !! Whether B is free in this simplified unconstrained fit.
      logical, intent(in), optional :: estimate_u !! Whether U is free in this simplified unconstrained fit.
      logical, intent(in), optional :: estimate_q !! Whether Q is free in this simplified unconstrained fit.
      logical, intent(in), optional :: estimate_z !! Whether Z is free in this simplified unconstrained fit.
      logical, intent(in), optional :: estimate_a !! Whether A is free in this simplified unconstrained fit.
      logical, intent(in), optional :: estimate_r !! Whether R is free in this simplified unconstrained fit.
      logical, intent(in), optional :: estimate_x0 !! Whether x0 is free in this simplified unconstrained fit.
      logical, intent(in), optional :: estimate_v0 !! Whether V0 is free in this simplified unconstrained fit.

      call marss_kem(start_model, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, &
         estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
   end subroutine marss_fit

   pure subroutine marss(start_model, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, &
      estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
      type(marss_model), intent(in) :: start_model !! Prepared time-invariant model replacing the R model-list/formula interface.
      integer, intent(in) :: max_iter !! Maximum EM iterations.
      real(dp), intent(in) :: tol !! Relative convergence tolerance for successive log likelihoods.
      type(marss_fit_result), intent(out) :: fit !! Fitted MARSS model and numerical diagnostics.
      logical, intent(in), optional :: estimate_b !! Whether B is estimated.
      logical, intent(in), optional :: estimate_u !! Whether U is estimated.
      logical, intent(in), optional :: estimate_q !! Whether Q is estimated.
      logical, intent(in), optional :: estimate_z !! Whether Z is estimated.
      logical, intent(in), optional :: estimate_a !! Whether A is estimated.
      logical, intent(in), optional :: estimate_r !! Whether R is estimated.
      logical, intent(in), optional :: estimate_x0 !! Whether x0 is estimated.
      logical, intent(in), optional :: estimate_v0 !! Whether V0 is estimated.

      call marss_kem(start_model, max_iter, tol, fit, estimate_b, estimate_u, estimate_q, &
         estimate_z, estimate_a, estimate_r, estimate_x0, estimate_v0)
   end subroutine marss

   pure subroutine em_update(current, kf, updated, info, eb, eu, eq, ez, ea, er, ex0, ev0)
      type(marss_model), intent(in) :: current !! Current EM iterate.
      type(marss_kf_result), intent(in) :: kf !! Smoothed sufficient statistics at the current iterate.
      type(marss_model), intent(out) :: updated !! Model after one unconstrained EM M-step.
      integer, intent(out) :: info !! Zero on success, nonzero if a required normal equation is singular.
      logical, intent(in) :: eb !! Update B when true.
      logical, intent(in) :: eu !! Update U when true.
      logical, intent(in) :: eq !! Update Q when true.
      logical, intent(in) :: ez !! Update Z when true.
      logical, intent(in) :: ea !! Update A when true.
      logical, intent(in) :: er !! Update R when true.
      logical, intent(in) :: ex0 !! Update x0 when true.
      logical, intent(in) :: ev0 !! Update V0 when true.
      type(marss_hatyt_result) :: hat
      real(dp), allocatable :: inv(:, :)
      real(dp), allocatable :: sss(:, :)
      real(dp), allocatable :: xss(:, :)
      real(dp), allocatable :: theta(:, :)
      real(dp), allocatable :: eprev(:, :)
      real(dp), allocatable :: ett(:, :)
      real(dp), allocatable :: etprev(:, :)
      real(dp), allocatable :: ess(:, :)
      real(dp), allocatable :: exs(:, :)
      real(dp), allocatable :: eyx(:, :)
      real(dp), allocatable :: yss(:, :)
      real(dp), allocatable :: ytheta(:, :)
      real(dp), allocatable :: prev_mean(:)
      real(dp), allocatable :: prev_cov(:, :)
      real(dp), allocatable :: cross(:, :)
      real(dp), allocatable :: aug_mean(:)
      real(dp) :: delta(size(current%x0))
      integer :: first_transition
      integer :: m
      integer :: n
      integer :: ntrans
      integer :: t
      integer :: tt

      updated = current
      info = 0
      m = size(current%b, 1)
      n = size(current%z, 1)
      tt = size(current%y, 2)
      first_transition = merge(2, 1, current%tinitx == 1)
      ntrans = tt - first_transition + 1

      if (ntrans > 0 .and. (eb .or. eu .or. eq)) then
         allocate(sss(m + 1, m + 1), xss(m, m + 1))
         sss = 0.0_dp
         xss = 0.0_dp
         do t = first_transition, tt
            call transition_moments(current, kf, t, prev_mean, prev_cov, cross)
            eprev = prev_cov + outer_product(prev_mean, prev_mean)
            ett = kf%p_smooth(:, :, t) + outer_product(kf%x_smooth(:, t), kf%x_smooth(:, t))
            etprev = cross + outer_product(kf%x_smooth(:, t), prev_mean)
            allocate(ess(m + 1, m + 1), exs(m, m + 1), aug_mean(m + 1))
            ess = 0.0_dp
            ess(1:m, 1:m) = eprev
            ess(1:m, m + 1) = prev_mean
            ess(m + 1, 1:m) = prev_mean
            ess(m + 1, m + 1) = 1.0_dp
            exs(:, 1:m) = etprev - outer_product(marss_c_effect_at(current, t), prev_mean)
            exs(:, m + 1) = kf%x_smooth(:, t) - marss_c_effect_at(current, t)
            sss = sss + ess
            xss = xss + exs
            deallocate(prev_mean, prev_cov, cross, eprev, ett, etprev, ess, exs, aug_mean)
         end do
         if (eb .and. eu) then
            call inverse_matrix(sss, inv, info)
            if (info /= 0) return
            theta = matmul(xss, inv)
            updated%b = theta(:, 1:m)
            updated%u = theta(:, m + 1)
            deallocate(inv, theta)
         else if (eb) then
            call update_b_only(current, kf, first_transition, updated%b, info)
            if (info /= 0) return
         else if (eu) then
            call update_u_only(current, kf, first_transition, updated%u)
         end if

         if (eq) then
            updated%q = 0.0_dp
            do t = first_transition, tt
               call transition_moments(current, kf, t, prev_mean, prev_cov, cross)
               eprev = prev_cov + outer_product(prev_mean, prev_mean)
               ett = kf%p_smooth(:, :, t) + outer_product(kf%x_smooth(:, t), kf%x_smooth(:, t))
               etprev = cross + outer_product(kf%x_smooth(:, t), prev_mean)
               allocate(ess(m + 1, m + 1), exs(m, m + 1), theta(m, m + 1))
               ess = 0.0_dp
               ess(1:m, 1:m) = eprev
               ess(1:m, m + 1) = prev_mean
               ess(m + 1, 1:m) = prev_mean
               ess(m + 1, m + 1) = 1.0_dp
               exs(:, 1:m) = etprev
               exs(:, m + 1) = kf%x_smooth(:, t)
               theta(:, 1:m) = marss_b_at(updated, t)
               theta(:, m + 1) = marss_u_at(updated, t)
               updated%q = updated%q + ett - matmul(theta, transpose(exs)) - &
                  matmul(exs, transpose(theta)) + matmul(matmul(theta, ess), transpose(theta))
               deallocate(prev_mean, prev_cov, cross, eprev, ett, etprev, ess, exs, theta)
            end do
            updated%q = symmetrize(updated%q / real(ntrans, dp))
            call stabilize_covariance(updated%q)
         end if
         deallocate(sss, xss)
      end if

      if (ez .or. ea .or. er) then
         call marss_hatyt(current, kf, hat)
         allocate(yss(m + 1, m + 1), eyx(n, m + 1))
         yss = 0.0_dp
         eyx = 0.0_dp
         do t = 1, tt
            allocate(ess(m + 1, m + 1))
            ess = 0.0_dp
            ess(1:m, 1:m) = kf%p_smooth(:, :, t) + &
               outer_product(kf%x_smooth(:, t), kf%x_smooth(:, t))
            ess(1:m, m + 1) = kf%x_smooth(:, t)
            ess(m + 1, 1:m) = kf%x_smooth(:, t)
            ess(m + 1, m + 1) = 1.0_dp
            yss = yss + ess
            eyx(:, 1:m) = eyx(:, 1:m) + hat%yxt(:, :, t) - &
               outer_product(marss_d_effect_at(current, t), kf%x_smooth(:, t))
            eyx(:, m + 1) = eyx(:, m + 1) + hat%yt(:, t) - marss_d_effect_at(current, t)
            deallocate(ess)
         end do
         if (ez .and. ea) then
            call inverse_matrix(yss, inv, info)
            if (info /= 0) return
            ytheta = matmul(eyx, inv)
            updated%z = ytheta(:, 1:m)
            updated%a = ytheta(:, m + 1)
            deallocate(inv, ytheta)
         else if (ez) then
            call update_z_only(current, kf, hat, updated%z, info)
            if (info /= 0) return
         else if (ea) then
            call update_a_only(current, kf, hat, updated%a)
         end if
         if (er) then
            allocate(ytheta(n, m + 1))
            updated%r = 0.0_dp
            do t = 1, tt
               ytheta(:, 1:m) = marss_z_at(updated, t)
               ytheta(:, m + 1) = marss_a_at(updated, t)
               allocate(ess(m + 1, m + 1), exs(n, m + 1))
               ess = 0.0_dp
               ess(1:m, 1:m) = kf%p_smooth(:, :, t) + &
                  outer_product(kf%x_smooth(:, t), kf%x_smooth(:, t))
               ess(1:m, m + 1) = kf%x_smooth(:, t)
               ess(m + 1, 1:m) = kf%x_smooth(:, t)
               ess(m + 1, m + 1) = 1.0_dp
               exs(:, 1:m) = hat%yxt(:, :, t)
               exs(:, m + 1) = hat%yt(:, t)
               updated%r = updated%r + hat%ot(:, :, t) - &
                  matmul(ytheta, transpose(exs)) - matmul(exs, transpose(ytheta)) + &
                  matmul(matmul(ytheta, ess), transpose(ytheta))
               deallocate(ess, exs)
            end do
            updated%r = symmetrize(updated%r / real(tt, dp))
            call stabilize_covariance(updated%r)
            deallocate(ytheta)
         end if
         deallocate(yss, eyx)
      end if

      if (current%tinitx == 0) then
         if (ex0) updated%x0 = kf%x0_smooth
         if (ev0) then
            delta = kf%x0_smooth - updated%x0
            updated%v0 = symmetrize(kf%v0_smooth + outer_product(delta, delta))
            call stabilize_covariance(updated%v0)
         end if
      else
         if (ex0) updated%x0 = kf%x_smooth(:, 1)
         if (ev0) then
            delta = kf%x_smooth(:, 1) - updated%x0
            updated%v0 = symmetrize(kf%p_smooth(:, :, 1) + outer_product(delta, delta))
            call stabilize_covariance(updated%v0)
         end if
      end if
   end subroutine em_update

   pure subroutine transition_moments(model, kf, t, mean_prev, cov_prev, cross)
      type(marss_model), intent(in) :: model !! Model providing the tinitx convention.
      type(marss_kf_result), intent(in) :: kf !! Smoother moments from which transition sufficient statistics are extracted.
      integer, intent(in) :: t !! Current observation time whose predecessor is needed.
      real(dp), allocatable, intent(out) :: mean_prev(:) !! Smoothed mean of x(t-1), or x0 when t=1 and tinitx=0.
      real(dp), allocatable, intent(out) :: cov_prev(:, :) !! Smoothed covariance of the predecessor state.
      real(dp), allocatable, intent(out) :: cross(:, :) !! Smoothed covariance Cov[x(t),x(t-1)|Y].
      integer :: m

      m = size(model%b, 1)
      allocate(mean_prev(m), cov_prev(m, m), cross(m, m))
      if (t == 1) then
         mean_prev = kf%x0_smooth
         cov_prev = kf%v0_smooth
      else
         mean_prev = kf%x_smooth(:, t - 1)
         cov_prev = kf%p_smooth(:, :, t - 1)
      end if
      cross = kf%p_lag(:, :, t)
   end subroutine transition_moments

   pure subroutine update_b_only(model, kf, first_transition, b, info)
      type(marss_model), intent(in) :: model !! Current model whose U is held fixed.
      type(marss_kf_result), intent(in) :: kf !! Smoothed moments used in the B normal equations.
      integer, intent(in) :: first_transition !! First transition time included in the update.
      real(dp), intent(out) :: b(:, :) !! Updated transition matrix B.
      integer, intent(out) :: info !! Zero on success, nonzero for a singular normal equation.
      real(dp), allocatable :: cross(:, :)
      real(dp), allocatable :: cov_prev(:, :)
      real(dp), allocatable :: inv(:, :)
      real(dp), allocatable :: mean_prev(:)
      real(dp) :: denom(size(b, 1), size(b, 2))
      real(dp) :: numer(size(b, 1), size(b, 2))
      integer :: t

      denom = 0.0_dp
      numer = 0.0_dp
      do t = first_transition, size(model%y, 2)
         call transition_moments(model, kf, t, mean_prev, cov_prev, cross)
         denom = denom + cov_prev + outer_product(mean_prev, mean_prev)
         numer = numer + cross + outer_product(kf%x_smooth(:, t) - marss_u_at(model, t), mean_prev)
         deallocate(mean_prev, cov_prev, cross)
      end do
      call inverse_matrix(denom, inv, info)
      if (info == 0) b = matmul(numer, inv)
   end subroutine update_b_only

   pure subroutine update_u_only(model, kf, first_transition, u)
      type(marss_model), intent(in) :: model !! Current model whose B matrix is held fixed.
      type(marss_kf_result), intent(in) :: kf !! Smoothed state means used to update the transition intercept.
      integer, intent(in) :: first_transition !! First transition time included in the update.
      real(dp), intent(out) :: u(:) !! Updated state intercept vector.
      real(dp), allocatable :: cross(:, :)
      real(dp), allocatable :: cov_prev(:, :)
      real(dp), allocatable :: mean_prev(:)
      integer :: ntrans
      integer :: t

      u = 0.0_dp
      ntrans = size(model%y, 2) - first_transition + 1
      do t = first_transition, size(model%y, 2)
         call transition_moments(model, kf, t, mean_prev, cov_prev, cross)
         u = u + kf%x_smooth(:, t) - matmul(marss_b_at(model, t), mean_prev) - marss_c_effect_at(model, t)
         deallocate(mean_prev, cov_prev, cross)
      end do
      u = u / real(ntrans, dp)
   end subroutine update_u_only

   pure subroutine update_z_only(model, kf, hat, z, info)
      type(marss_model), intent(in) :: model !! Current model whose A vector is held fixed.
      type(marss_kf_result), intent(in) :: kf !! Smoothed moments used in the observation loading update.
      type(marss_hatyt_result), intent(in) :: hat !! Conditional observation moments accounting for missing responses.
      real(dp), intent(out) :: z(:, :) !! Updated observation loading matrix Z.
      integer, intent(out) :: info !! Zero on success, nonzero for a singular normal equation.
      real(dp), allocatable :: inv(:, :)
      real(dp) :: denom(size(z, 2), size(z, 2))
      real(dp) :: numer(size(z, 1), size(z, 2))
      integer :: t

      denom = 0.0_dp
      numer = 0.0_dp
      do t = 1, size(model%y, 2)
         denom = denom + kf%p_smooth(:, :, t) + &
            outer_product(kf%x_smooth(:, t), kf%x_smooth(:, t))
         numer = numer + hat%yxt(:, :, t) - outer_product(marss_a_at(model, t), kf%x_smooth(:, t))
      end do
      call inverse_matrix(denom, inv, info)
      if (info == 0) z = matmul(numer, inv)
   end subroutine update_z_only

   pure subroutine update_a_only(model, kf, hat, a)
      type(marss_model), intent(in) :: model !! Current model whose Z matrix is held fixed.
      type(marss_kf_result), intent(in) :: kf !! Smoothed state means used in the observation intercept update.
      type(marss_hatyt_result), intent(in) :: hat !! Conditional observation means accounting for missing responses.
      real(dp), intent(out) :: a(:) !! Updated observation intercept vector A.
      integer :: t

      a = 0.0_dp
      do t = 1, size(model%y, 2)
         a = a + hat%yt(:, t) - matmul(marss_z_at(model, t), kf%x_smooth(:, t)) - marss_d_effect_at(model, t)
      end do
      a = a / real(size(model%y, 2), dp)
   end subroutine update_a_only

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the rank-one product.
      real(dp), intent(in) :: y(:) !! Right vector in the rank-one product.
      real(dp) :: a(size(x), size(y))

      a = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

   pure subroutine stabilize_covariance(a)
      real(dp), intent(inout) :: a(:, :) !! Symmetric covariance estimate whose diagonal is floored at a small positive value.
      integer :: i

      a = 0.5_dp * (a + transpose(a))
      do i = 1, size(a, 1)
         if (a(i, i) < 1.0e-10_dp) a(i, i) = 1.0e-10_dp
      end do
   end subroutine stabilize_covariance

end module marss_em
