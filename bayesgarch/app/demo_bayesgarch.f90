! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
program demo_bayesgarch
   use bayesgarch_garch, only : simulate_garch11_student
   use bayesgarch_kinds, only : dp
   use bayesgarch_rng, only : seed_rng
   use bayesgarch_sample, only : form_posterior_sample, posterior_mean, posterior_sd
   use bayesgarch_sampler, only : bayesgarch_control, bayesgarch_result, run_bayesgarch
   implicit none

   integer, parameter :: n = 1000
   real(dp), parameter :: true_parameters(4) = [0.03_dp, 0.08_dp, 0.86_dp, 8.0_dp]
   real(dp) :: y(n)
   real(dp) :: h(n + 1)
   real(dp), allocatable :: sample(:, :)
   real(dp) :: means(4)
   real(dp) :: sds(4)
   type(bayesgarch_control) :: control
   type(bayesgarch_result) :: result
   integer :: j

   call seed_rng(115)
   call simulate_garch11_student(n, true_parameters(:3), true_parameters(4), y, h, burn=500)

   control = bayesgarch_control()
   control%n_chains = 2
   control%n_iter = 3000
   control%start = [0.04_dp, 0.10_dp, 0.80_dp, 10.0_dp]
   control%enforce_stationarity = .true.
   call run_bayesgarch(y, result, control=control)
   call form_posterior_sample(result, burn=800, thin=5, sample=sample)
   means = posterior_mean(sample)
   sds = posterior_sd(sample)

   write(*, '(a)') "Bayesian GARCH(1,1) demonstration"
   write(*, '(a)') "parameter            true       post.mean         post.sd"
   do j = 1, 4
      write(*, '(a12,3(1x,f15.6))') parameter_name(j), true_parameters(j), means(j), sds(j)
   end do
   write(*, '(a,i0)') "posterior draws: ", size(sample, 1)
   do j = 1, control%n_chains
      write(*, '(a,i0,a,f8.3,a,f8.3,a,f8.3)') "chain ", j, &
         " alpha acceptance=", real(result%alpha_accept(j), dp) / real(control%n_iter - 1, dp), &
         " beta acceptance=", real(result%beta_accept(j), dp) / real(control%n_iter - 1, dp), &
         " nu updates=", real(result%nu_updates(j), dp) / real(control%n_iter - 1, dp)
   end do

contains

   pure function parameter_name(index) result(name)
      integer, intent(in) :: index
      character(len=12) :: name

      select case (index)
      case (1)
         name = "alpha0"
      case (2)
         name = "alpha1"
      case (3)
         name = "beta"
      case (4)
         name = "nu"
      case default
         name = "unknown"
      end select
   end function parameter_name

end program demo_bayesgarch
