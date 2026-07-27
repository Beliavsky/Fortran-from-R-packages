! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_simulation_density
   use sde
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none

   real(dp), parameter :: theta(3) = [1.0_dp, 0.7_dp, 0.35_dp]
   real(dp), parameter :: x0(2) = [0.2_dp, 0.5_dp]
   real(dp), allocatable :: path1(:), times(:), path2(:, :), bridge(:)
   real(dp) :: d_exact, d_euler, d_elerian, d_kessler, d_ozaki, d_shoji
   real(dp) :: h_density, p_density, rejection_rate
   integer :: failures, status

   failures = 0
   call seed_rng(314159265_i64)

   call brownian_motion(0.0_dp, 0.0_dp, 1.0_dp, 100, path1, times)
   call assert_true(size(path1) == 101 .and. size(times) == 101, "Brownian motion dimensions", failures)
   call assert_close(path1(1), 0.0_dp, 0.0_dp, "Brownian motion initial value", failures)
   call assert_true(all(ieee_is_finite(path1)), "Brownian motion finite", failures)

   call geometric_brownian_motion(1.0_dp, 0.05_dp, 0.2_dp, 0.0_dp, 1.0_dp, 100, path1, times)
   call assert_true(all(path1 > 0.0_dp), "geometric Brownian motion positive", failures)

   call brownian_bridge(-0.5_dp, 0.75_dp, 0.0_dp, 1.0_dp, 100, path1, times)
   call assert_close(path1(1), -0.5_dp, 0.0_dp, "Brownian bridge start", failures)
   call assert_close(path1(size(path1)), 0.75_dp, 0.0_dp, "Brownian bridge end", failures)

   call simulate_euler(x0, 0.0_dp, 0.01_dp, 50, ou_drift, constant_diffusion, theta, path2, &
      predictor_corrector=.true., diffusion_x=zero_coefficient)
   call check_matrix_path(path2, 51, 2, "Euler", failures)

   call simulate_milstein(x0, 0.0_dp, 0.01_dp, 50, ou_drift, constant_diffusion, &
      zero_coefficient, theta, path2)
   call check_matrix_path(path2, 51, 2, "Milstein", failures)

   call simulate_milstein_second_order(x0, 0.0_dp, 0.01_dp, 50, ou_drift, ou_drift_x, &
      zero_coefficient, constant_diffusion, zero_coefficient, zero_coefficient, theta, path2)
   call check_matrix_path(path2, 51, 2, "second-order Milstein", failures)

   call simulate_kps(x0, 0.0_dp, 0.01_dp, 50, ou_drift, ou_drift_x, zero_coefficient, &
      constant_diffusion, zero_coefficient, zero_coefficient, theta, path2)
   call check_matrix_path(path2, 51, 2, "KPS", failures)

   call simulate_ozaki(x0, 0.0_dp, 0.01_dp, 50, ou_drift, ou_drift_x, theta(3), theta, path2)
   call check_matrix_path(path2, 51, 2, "Ozaki", failures)

   call simulate_shoji(x0, 0.0_dp, 0.01_dp, 50, ou_drift, ou_drift_x, zero_coefficient, &
      zero_coefficient, theta(3), theta, path2)
   call check_matrix_path(path2, 51, 2, "Shoji", failures)

   call simulate_ou_exact(x0, 0.01_dp, 50, theta, path2)
   call check_matrix_path(path2, 51, 2, "exact OU", failures)

   call simulate_gbm_exact([1.0_dp, 2.0_dp], 0.01_dp, 50, [0.05_dp, 0.2_dp], path2)
   call check_matrix_path(path2, 51, 2, "exact GBM", failures)
   call assert_true(all(path2 > 0.0_dp), "exact GBM positive", failures)

   call simulate_cir_exact([0.5_dp, 1.0_dp], 0.01_dp, 50, [0.8_dp, 1.2_dp, 0.3_dp], path2)
   call check_matrix_path(path2, 51, 2, "exact CIR", failures)
   call assert_true(all(path2 >= 0.0_dp), "exact CIR nonnegative", failures)

   d_exact = ou_conditional_pdf(0.3_dp, 0.05_dp, 0.2_dp, theta)
   d_euler = transition_density_euler(0.3_dp, 0.05_dp, 0.2_dp, 0.0_dp, theta, &
      ou_drift, constant_diffusion)
   d_elerian = transition_density_elerian(0.3_dp, 0.05_dp, 0.2_dp, 0.0_dp, theta, &
      ou_drift, constant_diffusion, zero_coefficient)
   d_kessler = transition_density_kessler(0.3_dp, 0.05_dp, 0.2_dp, 0.0_dp, theta, &
      ou_drift, ou_drift_x, zero_coefficient, constant_diffusion, zero_coefficient, zero_coefficient)
   d_ozaki = transition_density_ozaki(0.3_dp, 0.05_dp, 0.2_dp, 0.0_dp, theta, &
      ou_drift, ou_drift_x, constant_diffusion)
   d_shoji = transition_density_shoji(0.3_dp, 0.05_dp, 0.2_dp, 0.0_dp, theta, &
      ou_drift, ou_drift_x, zero_coefficient, zero_coefficient, constant_diffusion)
   call assert_true(all([d_exact, d_euler, d_elerian, d_kessler, d_ozaki, d_shoji] > 0.0_dp), &
      "transition densities positive", failures)
   call assert_close(d_elerian, d_euler, 1.0e-13_dp, &
      "Elerian constant-diffusion fallback", failures)
   call assert_close(d_shoji, d_exact, 5.0e-12_dp, "Shoji OU density", failures)

   h_density = hermite_transition_density(0.1_dp, 0.2_dp, 0.3_dp, theta, zero_state, &
      zero_state, zero_state, zero_state, zero_state, zero_state, zero_state, identity_state, &
      unit_state)
   call assert_close(h_density, normal_pdf(0.3_dp, 0.2_dp, sqrt(0.1_dp)), 2.0e-13_dp, &
      "Hermite Brownian density", failures)

   p_density = pedersen_transition_density(0.2_dp, 0.3_dp, 0.05_dp, theta, ou_drift, &
      constant_diffusion, 4, 2000)
   call assert_true(ieee_is_finite(p_density) .and. p_density > 0.0_dp, &
      "Pedersen density positive", failures)
   call assert_true(abs(p_density-d_exact)/d_exact < 0.25_dp, &
      "Pedersen density near exact OU density", failures)

   call simulate_exact_ea(0.0_dp, 0.05_dp, 20, brownian_endpoint, zero_state, &
      [real(dp) ::], -0.1_dp, 0.1_dp, path1, rejection_rate, max_attempts=100000, status=status)
   call assert_true(status == 0, "exact acceptance simulator status", failures)
   call assert_true(size(path1) == 21 .and. all(ieee_is_finite(path1)), &
      "exact acceptance simulator path", failures)
   call assert_true(rejection_rate >= 0.0_dp .and. rejection_rate <= 1.0_dp, &
      "exact acceptance rejection rate", failures)

   call diffusion_bridge_euler(0.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, 0.01_dp, zero_drift, &
      unit_diffusion, [real(dp) ::], bridge, max_attempts=10000, status=status)
   call assert_true(status == 0, "Euler diffusion bridge status", failures)
   if (status == 0) then
      call assert_close(bridge(1), 0.0_dp, 0.0_dp, "Euler diffusion bridge start", failures)
      call assert_close(bridge(size(bridge)), 0.5_dp, 0.0_dp, "Euler diffusion bridge end", failures)
   end if

   if (failures /= 0) then
      write(*, '(a, i0)') "test_simulation_density: failures = ", failures
      error stop 1
   end if
   write(*, '(a)') "test_simulation_density: PASS"

contains

   pure function ou_drift(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = local_theta(1)-local_theta(2)*x+0.0_dp*t
   end function ou_drift

   pure function ou_drift_x(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = -local_theta(2)+0.0_dp*(t+x)
   end function ou_drift_x

   pure function constant_diffusion(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = local_theta(3)+0.0_dp*(t+x)
   end function constant_diffusion

   pure function zero_coefficient(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = 0.0_dp*(t+x+real(size(local_theta), dp))
   end function zero_coefficient

   pure function zero_drift(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = 0.0_dp*(t+x+real(size(local_theta), dp))
   end function zero_drift

   pure function unit_diffusion(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = 1.0_dp+0.0_dp*(t+x+real(size(local_theta), dp))
   end function unit_diffusion

   pure function zero_state(x, local_theta) result(value)
      real(dp), intent(in) :: x, local_theta(:)
      real(dp) :: value
      value = 0.0_dp*(x+real(size(local_theta), dp))
   end function zero_state

   pure function identity_state(x, local_theta) result(value)
      real(dp), intent(in) :: x, local_theta(:)
      real(dp) :: value
      value = x+0.0_dp*real(size(local_theta), dp)
   end function identity_state

   pure function unit_state(x, local_theta) result(value)
      real(dp), intent(in) :: x, local_theta(:)
      real(dp) :: value
      value = 1.0_dp+0.0_dp*(x+real(size(local_theta), dp))
   end function unit_state

   function brownian_endpoint(dt, current, local_theta) result(value)
      real(dp), intent(in) :: dt, current, local_theta(:)
      real(dp) :: value
      value = current+random_normal(sd=sqrt(dt))+0.0_dp*real(size(local_theta), dp)
   end function brownian_endpoint

   subroutine check_matrix_path(path, expected_rows, expected_columns, label, count)
      real(dp), intent(in) :: path(:, :)
      integer, intent(in) :: expected_rows, expected_columns
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      call assert_true(size(path, 1) == expected_rows .and. size(path, 2) == expected_columns, &
         trim(label)//" dimensions", count)
      call assert_true(all(ieee_is_finite(path)), trim(label)//" finite", count)
   end subroutine check_matrix_path

   subroutine assert_close(actual, target, tolerance, label, count)
      real(dp), intent(in) :: actual, target, tolerance
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not. ieee_is_finite(actual) .or. abs(actual-target) > tolerance) then
         write(*, '(a, 2(1x, es23.15), a, es12.4)') "FAIL "//trim(label)//":", actual, target, &
            " tolerance=", tolerance
         count = count+1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label, count)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not. condition) then
         write(*, '(a)') "FAIL "//trim(label)
         count = count+1
      end if
   end subroutine assert_true

end program test_simulation_density
