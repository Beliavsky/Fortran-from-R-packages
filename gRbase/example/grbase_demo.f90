program grbase_demo
  use grbase
  implicit none

  integer, allocatable :: cycle(:, :)
  integer, allocatable :: triangulated(:, :)
  type(table_t) :: counts
  type(table_t) :: probability
  real(dp) :: levels(4)

  cycle = adjacency_from_edges(4, reshape([1, 2, 3, 4, 2, 3, 4, 1], [2, 4]))
  levels = 2.0_dp
  triangulated = minimal_triangulation(cycle, levels)

  counts = make_table([1, 2], [2, 2], [4.0_dp, 6.0_dp, 3.0_dp, 7.0_dp])
  probability = table_normalize_first(counts)

  write (*, '(a,i0)') 'fill edges added to 4-cycle: ', (sum(triangulated) - sum(cycle)) / 2
  write (*, '(a,4(1x,f7.4))') 'p(variable 1 | variable 2):', probability%value
end program grbase_demo
