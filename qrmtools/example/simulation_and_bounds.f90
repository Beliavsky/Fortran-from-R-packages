! SPDX-License-Identifier: GPL-3.0-or-later
program simulation_and_bounds
  use qrmtools, only : dp, r_brownian, brownian_result, &
    rearrange_matrix, rearrangement_result, bound_worst_var
  implicit none

  type(brownian_result) :: paths
  type(rearrangement_result) :: bound
  real(dp) :: times(5)
  real(dp) :: quantiles(4,3)

  times = [0.0_dp,0.25_dp,0.5_dp,0.75_dp,1.0_dp]
  paths = r_brownian(3,times,d=2,process_type='BM',seed=24680)
  if(.not.paths%ok) error stop trim(paths%message)
  print '(a,2f12.6)', 'Final first path: ', paths%paths(1,5,:)

  quantiles = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                       1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                       1.0_dp,2.0_dp,3.0_dp,4.0_dp],[4,3])
  bound = rearrange_matrix(quantiles,bound_worst_var,tolerance=0.0_dp, &
    n_lookback=3,max_iterations=200,already_sorted=.true.)
  if(.not.bound%ok) error stop trim(bound%message)
  print '(a,f12.6)', 'Rearranged worst-VaR bound: ', bound%bound
end program simulation_and_bounds
