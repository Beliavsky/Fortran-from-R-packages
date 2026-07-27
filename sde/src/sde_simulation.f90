! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_simulation
   use sde_kinds, only : dp
   use sde_interfaces, only : sde_coefficient, transition_sampler, state_function
   use sde_random, only : random_normal, random_uniform, random_poisson
   use sde_special, only : safe_expm1, safe_log1p, expm1_over_x
   use sde_models, only : ou_conditional_random, gbm_conditional_random, cir_conditional_random
   implicit none
   private

   public :: brownian_motion
   public :: geometric_brownian_motion
   public :: brownian_bridge
   public :: simulate_euler
   public :: simulate_milstein
   public :: simulate_milstein_second_order
   public :: simulate_kps
   public :: simulate_ozaki
   public :: simulate_shoji
   public :: simulate_conditional
   public :: simulate_ou_exact
   public :: simulate_gbm_exact
   public :: simulate_cir_exact
   public :: simulate_exact_ea
   public :: diffusion_bridge_euler

contains

   subroutine brownian_motion(x0, t0, t_end, n_steps, path, times)
      real(dp), intent(in) :: x0, t0, t_end
      integer, intent(in) :: n_steps
      real(dp), allocatable, intent(out) :: path(:)
      real(dp), allocatable, intent(out), optional :: times(:)
      real(dp) :: dt
      integer :: i

      if (n_steps <= 0 .or. t_end <= t0) error stop "brownian_motion: invalid grid"
      dt = (t_end-t0)/real(n_steps, dp)
      allocate(path(n_steps+1))
      path(1) = x0
      do i = 2, n_steps+1
         path(i) = path(i-1)+random_normal(sd=sqrt(dt))
      end do
      if (present(times)) call make_times(t0, dt, n_steps, times)
   end subroutine brownian_motion

   subroutine geometric_brownian_motion(x0, drift_rate, sigma, t0, t_end, n_steps, path, times)
      real(dp), intent(in) :: x0, drift_rate, sigma, t0, t_end
      integer, intent(in) :: n_steps
      real(dp), allocatable, intent(out) :: path(:)
      real(dp), allocatable, intent(out), optional :: times(:)
      real(dp), allocatable :: bm(:), grid(:)

      if (x0 <= 0.0_dp .or. sigma < 0.0_dp) error stop "geometric_brownian_motion: invalid parameters"
      call brownian_motion(0.0_dp, t0, t_end, n_steps, bm, grid)
      allocate(path(n_steps+1))
      path = x0*exp((drift_rate-0.5_dp*sigma*sigma)*(grid-t0)+sigma*bm)
      if (present(times)) then
         allocate(times(n_steps+1))
         times = grid
      end if
   end subroutine geometric_brownian_motion

   subroutine brownian_bridge(x_start, x_end, t0, t_end, n_steps, path, times)
      real(dp), intent(in) :: x_start, x_end, t0, t_end
      integer, intent(in) :: n_steps
      real(dp), allocatable, intent(out) :: path(:)
      real(dp), allocatable, intent(out), optional :: times(:)
      real(dp), allocatable :: bm(:), grid(:)
      real(dp) :: duration

      if (n_steps <= 0 .or. t_end <= t0) error stop "brownian_bridge: invalid grid"
      duration = t_end-t0
      call brownian_motion(0.0_dp, t0, t_end, n_steps, bm, grid)
      allocate(path(n_steps+1))
      path = x_start+bm+(grid-t0)/duration*(x_end-x_start-bm(n_steps+1))
      path(1) = x_start
      path(n_steps+1) = x_end
      if (present(times)) then
         allocate(times(n_steps+1))
         times = grid
      end if
   end subroutine brownian_bridge

   subroutine simulate_euler(x0, t0, dt, n_steps, drift, diffusion, theta, path, &
         predictor_corrector, alpha, eta, diffusion_x)
      real(dp), intent(in) :: x0(:), t0, dt, theta(:)
      integer, intent(in) :: n_steps
      procedure(sde_coefficient) :: drift
      procedure(sde_coefficient) :: diffusion
      real(dp), allocatable, intent(out) :: path(:, :)
      logical, intent(in), optional :: predictor_corrector
      real(dp), intent(in), optional :: alpha, eta
      procedure(sde_coefficient), optional :: diffusion_x
      logical :: correct
      real(dp) :: alpha_value, eta_value, t, next_t, x, y, dw
      real(dp) :: d_current, d_predicted, s_current, s_predicted
      integer :: i, j, n_paths

      call validate_grid(x0, dt, n_steps, n_paths)
      correct = .true.
      alpha_value = 0.5_dp
      eta_value = 0.5_dp
      if (present(predictor_corrector)) correct = predictor_corrector
      if (present(alpha)) alpha_value = alpha
      if (present(eta)) eta_value = eta
      if (correct .and. .not. present(diffusion_x)) then
         error stop "simulate_euler: diffusion_x is required for predictor-corrector"
      end if
      allocate(path(n_steps+1, n_paths))
      path(1, :) = x0
      t = t0
      do i = 2, n_steps+1
         next_t = t+dt
         do j = 1, n_paths
            x = path(i-1, j)
            dw = random_normal(sd=sqrt(dt))
            s_current = diffusion(t, x, theta)
            if (correct) then
               y = x+drift(t, x, theta)*dt+s_current*dw
               s_predicted = diffusion(next_t, y, theta)
               d_current = drift(t, x, theta)-eta_value*s_current*diffusion_x(t, x, theta)
               d_predicted = drift(next_t, y, theta)-eta_value*s_predicted*diffusion_x(next_t, y, theta)
               path(i, j) = x+(alpha_value*d_predicted+(1.0_dp-alpha_value)*d_current)*dt+ &
                  (eta_value*s_predicted+(1.0_dp-eta_value)*s_current)*dw
            else
               path(i, j) = x+drift(t, x, theta)*dt+s_current*dw
            end if
         end do
         t = next_t
      end do
   end subroutine simulate_euler

   subroutine simulate_milstein(x0, t0, dt, n_steps, drift, diffusion, diffusion_x, theta, path)
      real(dp), intent(in) :: x0(:), t0, dt, theta(:)
      integer, intent(in) :: n_steps
      procedure(sde_coefficient) :: drift, diffusion, diffusion_x
      real(dp), allocatable, intent(out) :: path(:, :)
      real(dp) :: t, x, dw, d, s, sx
      integer :: i, j, n_paths

      call validate_grid(x0, dt, n_steps, n_paths)
      allocate(path(n_steps+1, n_paths))
      path(1, :) = x0
      t = t0
      do i = 2, n_steps+1
         do j = 1, n_paths
            x = path(i-1, j)
            dw = random_normal(sd=sqrt(dt))
            d = drift(t, x, theta)
            s = diffusion(t, x, theta)
            sx = diffusion_x(t, x, theta)
            path(i, j) = x+d*dt+s*dw+0.5_dp*s*sx*(dw*dw-dt)
         end do
         t = t+dt
      end do
   end subroutine simulate_milstein

   subroutine simulate_milstein_second_order(x0, t0, dt, n_steps, drift, drift_x, drift_xx, &
         diffusion, diffusion_x, diffusion_xx, theta, path)
      real(dp), intent(in) :: x0(:), t0, dt, theta(:)
      integer, intent(in) :: n_steps
      procedure(sde_coefficient) :: drift, drift_x, drift_xx
      procedure(sde_coefficient) :: diffusion, diffusion_x, diffusion_xx
      real(dp), allocatable, intent(out) :: path(:, :)
      real(dp) :: t, x, dw, d, dx, dxx, s, sx, sxx
      integer :: i, j, n_paths

      call validate_grid(x0, dt, n_steps, n_paths)
      allocate(path(n_steps+1, n_paths))
      path(1, :) = x0
      t = t0
      do i = 2, n_steps+1
         do j = 1, n_paths
            x = path(i-1, j)
            dw = random_normal(sd=sqrt(dt))
            d = drift(t, x, theta)
            dx = drift_x(t, x, theta)
            dxx = drift_xx(t, x, theta)
            s = diffusion(t, x, theta)
            sx = diffusion_x(t, x, theta)
            sxx = diffusion_xx(t, x, theta)
            path(i, j) = x+d*dt+s*dw+0.5_dp*s*sx*(dw*dw-dt)+ &
               dt**1.5_dp*(0.5_dp*d*sx+0.5_dp*dx*s+0.25_dp*s*sxx)*dw+ &
               dt*dt*(0.5_dp*d*dx+0.25_dp*dxx*s*s)
         end do
         t = t+dt
      end do
   end subroutine simulate_milstein_second_order

   subroutine simulate_kps(x0, t0, dt, n_steps, drift, drift_x, drift_xx, &
         diffusion, diffusion_x, diffusion_xx, theta, path)
      real(dp), intent(in) :: x0(:), t0, dt, theta(:)
      integer, intent(in) :: n_steps
      procedure(sde_coefficient) :: drift, drift_x, drift_xx
      procedure(sde_coefficient) :: diffusion, diffusion_x, diffusion_xx
      real(dp), allocatable, intent(out) :: path(:, :)
      real(dp) :: t, x, z, u, d, dx, dxx, s, sx, sxx, n1, n2
      integer :: i, j, n_paths

      call validate_grid(x0, dt, n_steps, n_paths)
      allocate(path(n_steps+1, n_paths))
      path(1, :) = x0
      t = t0
      do i = 2, n_steps+1
         do j = 1, n_paths
            x = path(i-1, j)
            n1 = random_normal()
            n2 = random_normal()
            z = sqrt(dt)*n1
            u = 0.5_dp*dt*z+sqrt(dt**3/12.0_dp)*n2
            d = drift(t, x, theta)
            dx = drift_x(t, x, theta)
            dxx = drift_xx(t, x, theta)
            s = diffusion(t, x, theta)
            sx = diffusion_x(t, x, theta)
            sxx = diffusion_xx(t, x, theta)
            path(i, j) = x+d*dt+s*z+0.5_dp*s*sx*(z*z-dt)+s*dx*u+ &
               0.5_dp*(d*dx+0.5_dp*s*s*dxx)*dt*dt+ &
               (d*sx+0.5_dp*s*s*sxx)*(z*dt-u)+ &
               0.5_dp*s*(sx*sx+s*sxx)*(z*z/3.0_dp-dt)*z
         end do
         t = t+dt
      end do
   end subroutine simulate_kps

   subroutine simulate_ozaki(x0, t0, dt, n_steps, drift, drift_x, sigma, theta, path)
      real(dp), intent(in) :: x0(:), t0, dt, sigma, theta(:)
      integer, intent(in) :: n_steps
      procedure(sde_coefficient) :: drift, drift_x
      real(dp), allocatable, intent(out) :: path(:, :)
      real(dp) :: t, x, d, dx, mean_value, variance_value, kx, adjustment
      integer :: i, j, n_paths

      call validate_grid(x0, dt, n_steps, n_paths)
      if (sigma < 0.0_dp) error stop "simulate_ozaki: sigma must be nonnegative"
      allocate(path(n_steps+1, n_paths))
      path(1, :) = x0
      t = t0
      do i = 2, n_steps+1
         do j = 1, n_paths
            x = path(i-1, j)
            d = drift(t, x, theta)
            dx = drift_x(t, x, theta)
            mean_value = x+d*dt*expm1_over_x(dx*dt)
            if (abs(x) > sqrt(epsilon(1.0_dp))) then
               adjustment = d*dt*expm1_over_x(dx*dt)/x
               if (adjustment > -1.0_dp) then
                  kx = safe_log1p(adjustment)/dt
               else
                  kx = dx
               end if
            else
               kx = dx
            end if
            variance_value = sigma*sigma*dt*expm1_over_x(2.0_dp*kx*dt)
            path(i, j) = random_normal(mean_value, sqrt(max(0.0_dp, variance_value)))
         end do
         t = t+dt
      end do
   end subroutine simulate_ozaki

   subroutine simulate_shoji(x0, t0, dt, n_steps, drift, drift_x, drift_xx, drift_t, sigma, theta, path)
      real(dp), intent(in) :: x0(:), t0, dt, sigma, theta(:)
      integer, intent(in) :: n_steps
      procedure(sde_coefficient) :: drift, drift_x, drift_xx, drift_t
      real(dp), allocatable, intent(out) :: path(:, :)
      real(dp) :: t, x, d, dx, m, z, mean_value, variance_value
      integer :: i, j, n_paths

      call validate_grid(x0, dt, n_steps, n_paths)
      if (sigma < 0.0_dp) error stop "simulate_shoji: sigma must be nonnegative"
      allocate(path(n_steps+1, n_paths))
      path(1, :) = x0
      t = t0
      do i = 2, n_steps+1
         do j = 1, n_paths
            x = path(i-1, j)
            d = drift(t, x, theta)
            dx = drift_x(t, x, theta)
            m = 0.5_dp*sigma*sigma*drift_xx(t, x, theta)+drift_t(t, x, theta)
            z = dx*dt
            mean_value = x+d*dt*expm1_over_x(z)+m*dt*dt*expm1_minus_x_over_x2(z)
            variance_value = sigma*sigma*dt*expm1_over_x(2.0_dp*z)
            path(i, j) = random_normal(mean_value, sqrt(max(0.0_dp, variance_value)))
         end do
         t = t+dt
      end do
   end subroutine simulate_shoji

   subroutine simulate_conditional(x0, dt, n_steps, sampler, theta, path)
      real(dp), intent(in) :: x0(:), dt, theta(:)
      integer, intent(in) :: n_steps
      procedure(transition_sampler) :: sampler
      real(dp), allocatable, intent(out) :: path(:, :)
      integer :: i, j, n_paths

      call validate_grid(x0, dt, n_steps, n_paths)
      allocate(path(n_steps+1, n_paths))
      path(1, :) = x0
      do i = 2, n_steps+1
         do j = 1, n_paths
            path(i, j) = sampler(dt, path(i-1, j), theta)
         end do
      end do
   end subroutine simulate_conditional

   subroutine simulate_ou_exact(x0, dt, n_steps, theta, path)
      real(dp), intent(in) :: x0(:), dt, theta(:)
      integer, intent(in) :: n_steps
      real(dp), allocatable, intent(out) :: path(:, :)
      call simulate_conditional(x0, dt, n_steps, ou_sampler, theta, path)
   contains
      function ou_sampler(local_dt, local_x0, local_theta) result(value)
         real(dp), intent(in) :: local_dt, local_x0, local_theta(:)
         real(dp) :: value
         value = ou_conditional_random(local_dt, local_x0, local_theta)
      end function ou_sampler
   end subroutine simulate_ou_exact

   subroutine simulate_gbm_exact(x0, dt, n_steps, theta, path)
      real(dp), intent(in) :: x0(:), dt, theta(:)
      integer, intent(in) :: n_steps
      real(dp), allocatable, intent(out) :: path(:, :)
      call simulate_conditional(x0, dt, n_steps, gbm_sampler, theta, path)
   contains
      function gbm_sampler(local_dt, local_x0, local_theta) result(value)
         real(dp), intent(in) :: local_dt, local_x0, local_theta(:)
         real(dp) :: value
         value = gbm_conditional_random(local_dt, local_x0, local_theta)
      end function gbm_sampler
   end subroutine simulate_gbm_exact

   subroutine simulate_cir_exact(x0, dt, n_steps, theta, path)
      real(dp), intent(in) :: x0(:), dt, theta(:)
      integer, intent(in) :: n_steps
      real(dp), allocatable, intent(out) :: path(:, :)
      call simulate_conditional(x0, dt, n_steps, cir_sampler, theta, path)
   contains
      function cir_sampler(local_dt, local_x0, local_theta) result(value)
         real(dp), intent(in) :: local_dt, local_x0, local_theta(:)
         real(dp) :: value
         value = cir_conditional_random(local_dt, local_x0, local_theta)
      end function cir_sampler
   end subroutine simulate_cir_exact

   subroutine simulate_exact_ea(x0, dt, n_steps, endpoint_sampler, psi, theta, k1, k2, path, &
         rejection_rate, max_attempts, status)
      real(dp), intent(in) :: x0, dt, theta(:), k1, k2
      integer, intent(in) :: n_steps
      procedure(transition_sampler) :: endpoint_sampler
      procedure(state_function) :: psi
      real(dp), allocatable, intent(out) :: path(:)
      real(dp), intent(out), optional :: rejection_rate
      integer, intent(in), optional :: max_attempts
      integer, intent(out), optional :: status
      real(dp), allocatable :: points(:), uniforms(:), bridge_values(:), w(:)
      real(dp) :: current, endpoint, intensity, temp
      integer :: accepted, rejected, attempts_limit, attempts, k, i
      logical :: accept

      if (dt <= 0.0_dp .or. n_steps <= 0 .or. k2 <= k1) then
         error stop "simulate_exact_ea: invalid arguments"
      end if
      intensity = k2-k1
      attempts_limit = 100000
      if (present(max_attempts)) attempts_limit = max_attempts
      allocate(path(n_steps+1))
      path(1) = x0
      current = x0
      accepted = 0
      rejected = 0
      attempts = 0
      do while (accepted < n_steps .and. attempts < attempts_limit)
         attempts = attempts+1
         endpoint = endpoint_sampler(dt, current, theta)
         k = random_poisson(intensity*dt)
         accept = .true.
         if (k > 0) then
            allocate(points(k+2), uniforms(k), bridge_values(k+2), w(k+2))
            points(1) = 0.0_dp
            points(k+2) = dt
            do i = 1, k
               points(i+1) = random_uniform()*dt
               uniforms(i) = random_uniform()*intensity
            end do
            ! Uniform heights are exchangeable, so only the event times need
            ! to be sorted before constructing the conditional bridge.
            call sort_real(points(2:k+1))
            w(1) = 0.0_dp
            do i = 2, k+2
               w(i) = w(i-1)+random_normal(sd=sqrt(points(i)-points(i-1)))
            end do
            do i = 1, k+2
               bridge_values(i) = current+w(i)+(endpoint-current-w(k+2))*points(i)/dt
            end do
            do i = 1, k
               temp = psi(bridge_values(i+1), theta)-k1
               if (temp < 0.0_dp .or. temp > intensity .or. temp > uniforms(i)) then
                  accept = .false.
                  exit
               end if
            end do
            deallocate(points, uniforms, bridge_values, w)
         end if
         if (accept) then
            accepted = accepted+1
            path(accepted+1) = endpoint
            current = endpoint
         else
            rejected = rejected+1
         end if
      end do
      if (present(status)) then
         if (accepted == n_steps) then
            status = 0
         else
            status = 1
         end if
      end if
      if (accepted < n_steps) path(accepted+2:) = current
      if (present(rejection_rate)) then
         rejection_rate = real(rejected, dp)/real(max(1, accepted+rejected), dp)
      end if
   end subroutine simulate_exact_ea

   subroutine diffusion_bridge_euler(x_start, x_end, t0, t_end, dt, drift, diffusion, theta, &
         bridge, max_attempts, status)
      real(dp), intent(in) :: x_start, x_end, t0, t_end, dt, theta(:)
      procedure(sde_coefficient) :: drift, diffusion
      real(dp), allocatable, intent(out) :: bridge(:)
      integer, intent(in), optional :: max_attempts
      integer, intent(out), optional :: status
      real(dp), allocatable :: forward(:, :), backward(:, :)
      real(dp) :: initial_forward(1), initial_backward(1)
      integer :: n_steps, attempts_limit, attempt, crossing, i
      logical :: starts_above

      if (t_end <= t0 .or. dt <= 0.0_dp) error stop "diffusion_bridge_euler: invalid grid"
      n_steps = nint((t_end-t0)/dt)
      if (n_steps < 2) error stop "diffusion_bridge_euler: at least two steps are required"
      attempts_limit = 10000
      if (present(max_attempts)) attempts_limit = max_attempts
      initial_forward = x_start
      initial_backward = x_end
      crossing = 0
      do attempt = 1, attempts_limit
         call simulate_euler(initial_forward, t0, (t_end-t0)/real(n_steps, dp), n_steps, &
            drift, diffusion, theta, forward, predictor_corrector=.false.)
         call simulate_euler(initial_backward, t0, (t_end-t0)/real(n_steps, dp), n_steps, &
            drift, diffusion, theta, backward, predictor_corrector=.false.)
         backward(:, 1) = backward(n_steps+1:1:-1, 1)
         starts_above = forward(1, 1) >= backward(1, 1)
         do i = 2, n_steps
            if ((starts_above .and. forward(i, 1) <= backward(i, 1)) .or. &
                (.not. starts_above .and. forward(i, 1) >= backward(i, 1))) then
               crossing = i
               exit
            end if
         end do
         if (crossing > 1 .and. crossing < n_steps+1) exit
      end do
      allocate(bridge(n_steps+1))
      if (crossing > 0) then
         bridge(1:crossing) = forward(1:crossing, 1)
         bridge(crossing+1:) = backward(crossing+1:, 1)
         bridge(1) = x_start
         bridge(n_steps+1) = x_end
         if (present(status)) status = 0
      else
         bridge = 0.0_dp
         if (present(status)) status = 1
      end if
   end subroutine diffusion_bridge_euler

   subroutine validate_grid(x0, dt, n_steps, n_paths)
      real(dp), intent(in) :: x0(:), dt
      integer, intent(in) :: n_steps
      integer, intent(out) :: n_paths
      n_paths = size(x0)
      if (n_paths <= 0 .or. n_steps <= 0 .or. dt <= 0.0_dp) then
         error stop "SDE simulation requires nonempty x0, n_steps>0, and dt>0"
      end if
   end subroutine validate_grid

   subroutine make_times(t0, dt, n_steps, times)
      real(dp), intent(in) :: t0, dt
      integer, intent(in) :: n_steps
      real(dp), allocatable, intent(out) :: times(:)
      integer :: i
      allocate(times(n_steps+1))
      do i = 1, n_steps+1
         times(i) = t0+real(i-1, dp)*dt
      end do
   end subroutine make_times

   pure function expm1_minus_x_over_x2(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (abs(x) < 1.0e-5_dp) then
         value = 0.5_dp+x/6.0_dp+x*x/24.0_dp+x*x*x/120.0_dp
      else
         value = (safe_expm1(x)-x)/(x*x)
      end if
   end function expm1_minus_x_over_x2

   subroutine sort_real(values)
      real(dp), intent(inout) :: values(:)
      real(dp) :: key
      integer :: i, j
      do i = 2, size(values)
         key = values(i)
         j = i-1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j+1) = values(j)
            j = j-1
         end do
         values(j+1) = key
      end do
   end subroutine sort_real

end module sde_simulation
