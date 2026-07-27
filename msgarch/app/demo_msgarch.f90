! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program demo_msgarch
   use msgarch
   implicit none
   type(msgarch_spec) :: spec, single
   type(simulation_result) :: simulation
   type(filter_result) :: filtered
   type(fit_result) :: fit
   type(risk_result) :: risk
   real(dp), allocatable :: volatility(:)

   call seed_rng(20260724)
   spec = create_spec([character(len=12) :: 'sGARCH', 'gjrGARCH'], &
      [character(len=8) :: 'norm', 'std'])
   spec%regime(1)%omega = 0.04_dp
   spec%regime(1)%alpha = 0.08_dp
   spec%regime(1)%beta = 0.88_dp
   spec%regime(2)%omega = 0.12_dp
   spec%regime(2)%alpha = 0.10_dp
   spec%regime(2)%gamma = 0.08_dp
   spec%regime(2)%beta = 0.72_dp
   spec%regime(2)%shape = 9.0_dp
   spec%transition(1,:) = [0.96_dp, 0.04_dp]
   spec%transition(2,:) = [0.08_dp, 0.92_dp]

   simulation = simulate_msgarch(spec, 400, 1)
   filtered = hamilton_filter(spec, simulation%draw(1,:))
   volatility = conditional_volatility(spec, simulation%draw(1,:))
   risk = risk_forecast(spec, simulation%draw(1,:), [0.01_dp, 0.05_dp], 3, 4000)

   single = extract_regime(spec, 1)
   fit = fit_ml(single, simulation%draw(1,:), max_iterations=400)

   write(*,'(a,f14.4)') 'Two-regime log likelihood: ', filtered%loglik
   write(*,'(a,2f10.4)') 'Final filtered probabilities: ', filtered%filtered(size(filtered%filtered,1),:)
   write(*,'(a,f10.4)') 'Final conditional volatility: ', volatility(size(volatility))
   write(*,'(a,2f10.4)') 'One-step VaR (1%, 5%): ', risk%var(1,:)
   write(*,'(a,3f10.4)') 'Single-regime fitted parameters: ', fit%parameters
end program demo_msgarch
