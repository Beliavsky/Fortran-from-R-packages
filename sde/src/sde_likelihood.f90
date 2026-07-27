! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_likelihood
   use sde_kinds, only : dp
   use sde_interfaces, only : sde_coefficient, state_function
   use sde_special, only : normal_pdf, normal_logpdf, nan_dp
   use sde_random, only : random_normal
   use sde_models, only : ou_conditional_pdf, gbm_conditional_pdf, cir_conditional_pdf
   implicit none
   private

   public :: euler_log_likelihood
   public :: hermite_transition_density
   public :: hermite_log_likelihood
   public :: pedersen_transition_density
   public :: pedersen_log_likelihood
   public :: ou_log_likelihood
   public :: gbm_log_likelihood
   public :: cir_log_likelihood

contains

   function euler_log_likelihood(x, dt, theta, drift, diffusion) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      procedure(sde_coefficient) :: drift, diffusion
      real(dp) :: value
      real(dp) :: mean_value, sd_value
      integer :: i

      if (size(x) < 2 .or. dt <= 0.0_dp) error stop "euler_log_likelihood: invalid data or dt"
      value = 0.0_dp
      do i = 2, size(x)
         mean_value = x(i-1)+drift(0.0_dp, x(i-1), theta)*dt
         sd_value = abs(diffusion(0.0_dp, x(i-1), theta))*sqrt(dt)
         value = value+normal_logpdf(x(i), mean_value, sd_value)
      end do
   end function euler_log_likelihood

   function hermite_transition_density(dt, x0, x1, theta, m0, m1, m2, m3, m4, m5, m6, &
         transform, diffusion) result(value)
      real(dp), intent(in) :: dt, x0, x1, theta(:)
      procedure(state_function) :: m0, m1, m2, m3, m4, m5, m6
      procedure(state_function) :: transform, diffusion
      real(dp) :: value
      real(dp) :: y0, y1, sd, ssd, z

      if (dt <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      sd = sqrt(dt)
      y0 = transform(x0, theta)
      y1 = transform(x1, theta)
      ssd = abs(diffusion(x1, theta))*sd
      if (ssd <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      z = (y1-y0)/sd
      value = hermite_kernel(dt, m0(y0, theta), m1(y0, theta), m2(y0, theta), &
         m3(y0, theta), m4(y0, theta), m5(y0, theta), m6(y0, theta), z, ssd)
   end function hermite_transition_density

   function hermite_log_likelihood(x, dt, theta, m0, m1, m2, m3, m4, m5, m6, &
         transform, diffusion) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      procedure(state_function) :: m0, m1, m2, m3, m4, m5, m6
      procedure(state_function) :: transform, diffusion
      real(dp) :: value
      real(dp) :: density
      integer :: i

      if (size(x) < 2 .or. dt <= 0.0_dp) error stop "hermite_log_likelihood: invalid data or dt"
      value = 0.0_dp
      do i = 2, size(x)
         density = hermite_transition_density(dt, x(i-1), x(i), theta, m0, m1, m2, m3, m4, m5, m6, &
            transform, diffusion)
         if (density <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
         end if
         value = value+log(density)
      end do
   end function hermite_log_likelihood

   function pedersen_transition_density(x0, x1, dt, theta, drift, diffusion, &
         n_substeps, n_simulations) result(value)
      real(dp), intent(in) :: x0, x1, dt, theta(:)
      procedure(sde_coefficient) :: drift, diffusion
      integer, intent(in) :: n_substeps, n_simulations
      real(dp) :: value
      real(dp) :: sub_dt, sub_sd, state_plus, state_minus, z, density, total
      integer :: pair, step, count, n_pairs

      if (dt <= 0.0_dp .or. n_substeps <= 0 .or. n_simulations <= 0) then
         error stop "pedersen_transition_density: invalid simulation controls"
      end if
      sub_dt = dt/real(n_substeps, dp)
      sub_sd = sqrt(sub_dt)
      total = 0.0_dp
      count = 0
      n_pairs = n_simulations/2
      do pair = 1, n_pairs
         state_plus = x0
         state_minus = x0
         do step = 1, n_substeps-1
            z = random_normal()
            state_plus = state_plus+drift(0.0_dp, state_plus, theta)*sub_dt+ &
               diffusion(0.0_dp, state_plus, theta)*sub_sd*z
            state_minus = state_minus+drift(0.0_dp, state_minus, theta)*sub_dt- &
               diffusion(0.0_dp, state_minus, theta)*sub_sd*z
         end do
         density = normal_pdf(x1, state_plus+drift(0.0_dp, state_plus, theta)*sub_dt, &
            abs(diffusion(0.0_dp, state_plus, theta))*sub_sd)
         if (density >= 0.0_dp) then
            total = total+density
            count = count+1
         end if
         density = normal_pdf(x1, state_minus+drift(0.0_dp, state_minus, theta)*sub_dt, &
            abs(diffusion(0.0_dp, state_minus, theta))*sub_sd)
         if (density >= 0.0_dp) then
            total = total+density
            count = count+1
         end if
      end do
      if (2*n_pairs < n_simulations) then
         state_plus = x0
         do step = 1, n_substeps-1
            state_plus = state_plus+drift(0.0_dp, state_plus, theta)*sub_dt+ &
               diffusion(0.0_dp, state_plus, theta)*random_normal(sd=sub_sd)
         end do
         density = normal_pdf(x1, state_plus+drift(0.0_dp, state_plus, theta)*sub_dt, &
            abs(diffusion(0.0_dp, state_plus, theta))*sub_sd)
         total = total+density
         count = count+1
      end if
      if (count == 0) then
         value = nan_dp()
      else
         value = total/real(count, dp)
      end if
   end function pedersen_transition_density

   function pedersen_log_likelihood(x, dt, theta, drift, diffusion, n_substeps, n_simulations) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      procedure(sde_coefficient) :: drift, diffusion
      integer, intent(in) :: n_substeps, n_simulations
      real(dp) :: value
      real(dp) :: density
      integer :: i

      if (size(x) < 2) error stop "pedersen_log_likelihood: at least two observations are required"
      value = 0.0_dp
      do i = 2, size(x)
         density = pedersen_transition_density(x(i-1), x(i), dt, theta, drift, diffusion, &
            n_substeps, n_simulations)
         if (density <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
         end if
         value = value+log(density)
      end do
   end function pedersen_log_likelihood

   function ou_log_likelihood(x, dt, theta) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 2, size(x)
         value = value+ou_conditional_pdf(x(i), dt, x(i-1), theta, .true.)
      end do
   end function ou_log_likelihood

   function gbm_log_likelihood(x, dt, theta) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 2, size(x)
         value = value+gbm_conditional_pdf(x(i), dt, x(i-1), theta, .true.)
      end do
   end function gbm_log_likelihood

   function cir_log_likelihood(x, dt, theta) result(value)
      real(dp), intent(in) :: x(:), dt, theta(:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 2, size(x)
         value = value+cir_conditional_pdf(x(i), dt, x(i-1), theta, .true.)
      end do
   end function cir_log_likelihood

   pure function hermite_kernel(delta, mu0, mu1, mu2, mu3, mu4, mu5, mu6, z, ssd) result(value)
      real(dp), intent(in) :: delta, mu0, mu1, mu2, mu3, mu4, mu5, mu6, z, ssd
      real(dp) :: value
      real(dp) :: eta1, eta2, eta3, eta4, eta5, eta6
      real(dp) :: mu02, mu03, mu04, mu05, mu06, mu12, mu13, mu22
      real(dp) :: polynomial

      mu02 = mu0**2
      mu03 = mu0**3
      mu04 = mu0**4
      mu05 = mu0**5
      mu06 = mu0**6
      mu12 = mu1**2
      mu13 = mu1**3
      mu22 = mu2**2

      eta1 = -mu0*sqrt(delta)-(2.0_dp*mu0*mu1+mu2)*delta**1.5_dp/4.0_dp- &
         (4.0_dp*mu0*mu12+4.0_dp*mu02*mu2+6.0_dp*mu1*mu2+4.0_dp*mu0*mu3+mu4)* &
         delta**2.5_dp/24.0_dp
      eta2 = (mu02+mu1)*delta/2.0_dp+ &
         (6.0_dp*mu02*mu1+4.0_dp*mu12+7.0_dp*mu0*mu2+2.0_dp*mu3)*delta**2/12.0_dp+ &
         (28.0_dp*mu02*mu12+28.0_dp*mu02*mu3+16.0_dp*mu13+16.0_dp*mu03*mu2+ &
          88.0_dp*mu0*mu1*mu2+21.0_dp*mu22+32.0_dp*mu1*mu3+16.0_dp*mu0*mu4+3.0_dp*mu5)* &
         delta**3/96.0_dp
      eta3 = -(mu03+3.0_dp*mu0*mu1+mu2)*delta**1.5_dp/6.0_dp- &
         (12.0_dp*mu03*mu1+28.0_dp*mu0*mu12+22.0_dp*mu02*mu2+24.0_dp*mu1*mu2+ &
          14.0_dp*mu0*mu3+3.0_dp*mu4)*delta**2.5_dp/48.0_dp
      eta4 = (mu04+6.0_dp*mu02*mu1+3.0_dp*mu12+4.0_dp*mu0*mu2+mu3)*delta**2/24.0_dp+ &
         (20.0_dp*mu04*mu1+50.0_dp*mu03*mu2+100.0_dp*mu02*mu12+50.0_dp*mu02*mu3+ &
          23.0_dp*mu0*mu4+180.0_dp*mu0*mu1*mu2+40.0_dp*mu13+34.0_dp*mu22+ &
          52.0_dp*mu1*mu3+4.0_dp*mu5)*delta**3/240.0_dp
      eta5 = -(mu05+10.0_dp*mu03*mu1+15.0_dp*mu0*mu12+10.0_dp*mu02*mu2+ &
         10.0_dp*mu1*mu2+5.0_dp*mu0*mu3+mu4)*delta**2.5_dp/120.0_dp
      eta6 = (mu06+15.0_dp*mu04*mu1+15.0_dp*mu13+20.0_dp*mu03*mu2+15.0_dp*mu0*mu3+ &
         45.0_dp*mu02*mu12+10.0_dp*mu22+15.0_dp*mu02*mu3+60.0_dp*mu0*mu1*mu2+ &
         6.0_dp*mu0*mu4+mu5+0.0_dp*mu6)*delta**3/720.0_dp

      polynomial = 1.0_dp+eta1*h1(z)+eta2*h2(z)+eta3*h3(z)+eta4*h4(z)+eta5*h5(z)+eta6*h6(z)
      value = normal_pdf(z)*polynomial/ssd
   end function hermite_kernel

   pure function h1(z) result(value)
      real(dp), intent(in) :: z
      real(dp) :: value
      value = -z
   end function h1

   pure function h2(z) result(value)
      real(dp), intent(in) :: z
      real(dp) :: value
      value = -1.0_dp+z*z
   end function h2

   pure function h3(z) result(value)
      real(dp), intent(in) :: z
      real(dp) :: value
      value = 3.0_dp*z-z**3
   end function h3

   pure function h4(z) result(value)
      real(dp), intent(in) :: z
      real(dp) :: value
      value = 3.0_dp-6.0_dp*z*z+z**4
   end function h4

   pure function h5(z) result(value)
      real(dp), intent(in) :: z
      real(dp) :: value
      value = -15.0_dp*z+10.0_dp*z**3-z**5
   end function h5

   pure function h6(z) result(value)
      real(dp), intent(in) :: z
      real(dp) :: value
      value = -15.0_dp+45.0_dp*z*z-15.0_dp*z**4+z**6
   end function h6

end module sde_likelihood
