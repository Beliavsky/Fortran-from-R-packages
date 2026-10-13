! SPDX-License-Identifier: GPL-3.0-only
program financial_network
   use network_risk_measures
   implicit none
   real(dp) :: assets(3, 3), buffer(3), susceptibility(3), reconstructed(3, 3), error
   type(contagion_result) :: simulation
   integer :: status, i
   assets = reshape([0.0_dp, 10.0_dp, 3.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp], [3, 3])
   buffer = [2.0_dp, 5.0_dp, 2.0_dp]
   call impact_susceptibility(assets, susceptibility, status, buffer)
   if (status /= nr_ok) error stop 'susceptibility failed'
   print '(a,3f10.5)', 'Impact susceptibility: ', susceptibility
   call contagion(assets, simulation, status, buffer)
   if (status /= nr_ok) error stop 'contagion failed'
   do i = 1, 3
      print '(a,i0,a,3f10.5)', 'Default vertex ', i, ': final stress ', simulation%final(:, i)
   end do
   call matrix_estimation([3.0_dp, 2.0_dp, 1.0_dp], [2.0_dp, 2.0_dp, 2.0_dp], reconstructed, status, &
      abs_tol=1.0e-10_dp, residual=error)
   if (status /= nr_ok) error stop 'reconstruction failed'
   print '(a,es12.4)', 'Maximum-entropy marginal residual: ', error
   call matrix_estimation([3.0_dp, 2.0_dp, 1.0_dp], [2.0_dp, 2.0_dp, 2.0_dp], reconstructed, status, &
      method='md', seed=42, residual=error)
   if (status /= nr_ok) error stop 'stochastic reconstruction failed'
   print '(a,es12.4)', 'Minimum-density marginal residual: ', error
end program financial_network
