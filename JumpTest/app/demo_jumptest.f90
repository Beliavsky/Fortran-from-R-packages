! SPDX-License-Identifier: MIT
program demo_jumptest
  use jumptest, only : dp, i8, simulation_result, pcombine_result, adjp_result, &
    sv1fj, pcombine, ppool
  implicit none

  integer, parameter :: bins = 78, days = 20
  type(simulation_result) :: simulated
  type(pcombine_result) :: combined
  type(adjp_result) :: pooled
  real(dp) :: returns(bins, days)
  real(dp), allocatable :: full_price(:)
  character(len=4) :: methods(3)
  integer :: day, first, last

  call sv1fj(bins, days, simulated, lambda=0.25_dp, seed=8675309_i8)
  allocate(full_price(size(simulated%price) + 1))
  full_price(1) = 3.0_dp
  full_price(2:) = simulated%price
  do day = 1, days
    first = (day - 1)*bins + 1
    last = day*bins
    returns(:, day) = full_price(first + 1:last + 1) - full_price(first:last)
  end do

  methods = ['BNS ', 'Amed', 'Amin']
  call pcombine(returns, methods, combined)
  call ppool(combined%pvalue, pooled, 'SD')

  print '(a)', 'JumpTest modern Fortran demo'
  print '(a,i0)', 'Simulated jump events: ', sum(simulated%jump_count)
  print '(a,i0)', 'Days with pooled BH p-value below 5%: ', count(pooled%adjp < 0.05_dp)
  print '(a,f10.6)', 'Minimum pooled p-value: ', minval(pooled%pvalue)
end program demo_jumptest
