! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program mcmc_example
   use msgarch
   implicit none
   type(msgarch_spec) :: spec
   type(simulation_result) :: simulation
   type(mcmc_result) :: fit

   call seed_rng(77123)
   spec = create_spec([character(len=12) :: 'sGARCH'], [character(len=8) :: 'std'])
   spec%regime(1)%omega = 0.05_dp
   spec%regime(1)%alpha = 0.10_dp
   spec%regime(1)%beta = 0.84_dp
   spec%regime(1)%shape = 8.0_dp
   simulation = simulate_msgarch(spec, 300, 1)
   fit = fit_mcmc(spec, simulation%draw(1,:), 900, 300, 5, &
      proposal_sd=[0.007_dp, 0.007_dp, 0.007_dp, 0.06_dp])

   write(*,'(a,f8.3)') 'Acceptance rate: ', fit%acceptance_rate
   write(*,'(a,4f11.5)') 'Posterior means: ', fit%posterior_mean
   write(*,'(a,f12.3)') 'DIC: ', fit%dic
end program mcmc_example
