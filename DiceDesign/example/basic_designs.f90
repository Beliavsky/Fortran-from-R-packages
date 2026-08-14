program basic_designs
  use iso_fortran_env, only : int64
  use dicedesign, only : dp, lhs_design, discrepancy_value, mindist, &
    discrep_ese_lhs, lhs_optimization_result
  implicit none

  real(dp), allocatable :: design(:, :)
  type(lhs_optimization_result) :: optimized
  real(dp) :: before, after

  call lhs_design(12, 3, design, randomized=.false., seed=12345_int64)
  before = discrepancy_value(design, 'C2')

  call discrep_ese_lhs(design, optimized, inner_iterations=20, candidates=12, &
    outer_iterations=2, criterion='C2', seed=12345_int64)
  after = discrepancy_value(optimized%design, 'C2')

  print '(a,f10.6)', 'Initial centered-LHS C2 discrepancy: ', before
  print '(a,f10.6)', 'Optimized C2 discrepancy:            ', after
  print '(a,f10.6)', 'Optimized minimum distance:          ', mindist(optimized%design)
end program basic_designs
