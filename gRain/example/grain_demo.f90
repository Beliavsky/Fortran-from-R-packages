program grain_demo
  use grain, only : dp, table_t, grain_network_t, grain_ok
  use grain, only : make_cpt, compile_network, propagate_ls
  use grain, only : set_hard_evidence, probability_of_evidence, query_marginal
  implicit none

  type(table_t) :: cpts(3)
  type(table_t) :: marginal
  type(grain_network_t) :: network
  integer :: status
  real(dp) :: probability

  cpts(1) = make_cpt([1], [2], [0.6_dp, 0.4_dp])
  cpts(2) = make_cpt([2, 1], [2, 2], [0.7_dp, 0.3_dp, 0.2_dp, 0.8_dp])
  cpts(3) = make_cpt([3, 2], [2, 2], [0.9_dp, 0.1_dp, 0.4_dp, 0.6_dp])

  call compile_network(cpts, network, status)
  if (status /= grain_ok) error stop 'network compilation failed'

  call propagate_ls(network, status)
  if (status /= grain_ok) error stop 'initial propagation failed'

  call set_hard_evidence(network, 3, 1, status)
  if (status /= grain_ok) error stop 'evidence insertion failed'

  probability = probability_of_evidence(network, status)
  if (status /= grain_ok) error stop 'evidence probability failed'

  call query_marginal(network, 1, marginal, status)
  if (status /= grain_ok) error stop 'marginal query failed'

  write(*, '(a,f10.6)') 'P(node 3 = state 1) = ', probability
  write(*, '(a,2f10.6)') 'P(node 1 | evidence) = ', marginal%value
end program grain_demo
