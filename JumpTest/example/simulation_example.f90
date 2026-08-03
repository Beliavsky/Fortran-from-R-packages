! SPDX-License-Identifier: MIT
program simulation_example
  use jumptest, only : dp, i8, simulation_result, sv1fj, sv2f
  implicit none

  type(simulation_result) :: one_factor, two_factor

  call sv1fj(390, 2, one_factor, lambda=0.5_dp, seed=20260731_i8)
  print '(a,i0)', 'SV1FJ observations: ', size(one_factor%price)
  print '(a,i0)', 'Jump events: ', sum(one_factor%jump_count)
  print '(a,f12.6)', 'Final log price: ', one_factor%price(size(one_factor%price))

  call sv2f(390, 2, two_factor, seed=20260731_i8)
  print '(a,f12.6)', 'SV2F final log price: ', two_factor%price(size(two_factor%price))
  print '(a,2(f12.6,1x))', 'Final factors: ', &
    two_factor%variance(size(two_factor%variance, 1), :)
end program simulation_example
