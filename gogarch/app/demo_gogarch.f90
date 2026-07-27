! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program demo_gogarch
   use gogarch
   implicit none
   integer, parameter :: n = 400, m = 2
   real(dp) :: factors(n,m), factor_variance(n,m), data(n,m), true_mixing(m,m)
   real(dp) :: means(5,m), covariance(m,m,5)
   real(dp) :: factor_means(m), omegas(m), arch(m,1), leverage(m,1), garch(m,1)
   real(dp) :: delta(m), shape(m), skew(m)
   type(univariate_spec) :: spec
   type(gogarch_fit) :: fit
   integer :: j

   call seed_rng(20260723)
   call simulate_aparch(n,0.0_dp,0.04_dp,[0.08_dp],[0.18_dp],[0.84_dp],1.4_dp,'sstd',8.0_dp,1.15_dp, &
      factors(:,1),factor_variance(:,1),burnin=600)
   call simulate_aparch(n,0.0_dp,0.06_dp,[0.11_dp],[-0.10_dp],[0.79_dp],1.6_dp,'sstd',8.0_dp,1.15_dp, &
      factors(:,2),factor_variance(:,2),burnin=600)
   true_mixing = reshape([1.0_dp,-0.3_dp,0.5_dp,0.9_dp],[m,m])
   data = matmul(factors,transpose(true_mixing))

   spec%model = 'aparch'
   spec%distribution = 'sstd'
   spec%p = 1
   spec%o = 1
   spec%q = 1
   spec%delta = 1.5_dp
   spec%shape = 8.0_dp
   spec%skew = 1.15_dp
   spec%fit_delta = .true.
   spec%fit_shape = .false.
   spec%fit_skew = .false.
   fit = fit_gogarch_ica(data,max_ica_iterations=500,max_garch_iterations=450,factor_spec=spec)
   call factor_coefficients_full(fit,factor_means,omegas,arch,leverage,garch,delta,shape,skew)
   call forecast_gogarch(fit,5,means,covariance)

   write(*,'(a,a)') 'method: ',trim(fit%method)
   write(*,'(a,a)') 'factor model: ',trim(fit%factor_spec%model)
   write(*,'(a,a)') 'distribution: ',trim(fit%factor_spec%distribution)
   write(*,'(a,i0)') 'status: ',fit%status
   write(*,'(a,f14.5)') 'log likelihood: ',fit%log_likelihood
   write(*,'(a,es12.4)') 'maximum reconstruction error: ',reconstruction_error(fit)
   write(*,'(a)') 'factor parameters: mean omega alpha gamma beta delta shape skew'
   do j = 1, m
      write(*,'(i3,8(1x,f11.6))') j,factor_means(j),omegas(j),arch(j,1),leverage(j,1),garch(j,1), &
         delta(j),shape(j),skew(j)
   end do
   write(*,'(a)') 'one-step covariance forecast:'
   do j = 1, m
      write(*,'(*(1x,f12.6))') covariance(j,:,1)
   end do
end program demo_gogarch
