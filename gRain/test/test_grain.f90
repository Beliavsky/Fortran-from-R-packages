program test_grain
  use grain, only : dp, table_t, grain_network_t
  use grain, only : grain_ok, make_cpt, compile_network, propagate_ls
  use grain, only : query_marginal, query_joint, query_conditional
  use grain, only : set_hard_evidence, set_soft_evidence, clear_evidence
  use grain, only : probability_of_evidence, simulate_network, replace_network_cpt
  use grain, only : logical_and_cpt, logical_or_cpt, mendel_probability
  use grain, only : estimate_cpts_from_data, estimate_clique_marginals
  use grain, only : marginals_to_potentials, potentials_to_marginals
  use grbase_tables, only : valid_table, table_equal
  implicit none

  type(table_t) :: cpts(3)
  type(table_t) :: tab
  type(table_t) :: joint
  type(table_t) :: cond
  type(table_t), allocatable :: estimated(:)
  type(table_t), allocatable :: marginals(:)
  type(table_t), allocatable :: potentials(:)
  type(table_t), allocatable :: marginals2(:)
  type(grain_network_t) :: net
  integer, allocatable :: sim1(:, :)
  integer, allocatable :: sim2(:, :)
  integer :: status
  integer :: failures
  integer :: i
  real(dp) :: pe
  integer :: data(10, 3)

  failures = 0

  cpts(1) = make_cpt([1], [2], [0.6_dp, 0.4_dp])
  cpts(2) = make_cpt([2, 1], [2, 2], [0.7_dp, 0.3_dp, 0.2_dp, 0.8_dp])
  cpts(3) = make_cpt([3, 2], [2, 2], [0.9_dp, 0.1_dp, 0.4_dp, 0.6_dp])
  call check(valid_table(cpts(1)), 'make_cpt root', failures)
  call check(valid_table(cpts(2)), 'make_cpt conditional', failures)

  call compile_network(cpts, net, status)
  call check(status == grain_ok, 'compile chain network', failures)
  call check(net%rip%cliques%count == 2, 'chain gives two maximal cliques', failures)
  call propagate_ls(net, status)
  call check(status == grain_ok, 'propagate chain network', failures)

  call query_marginal(net, 1, tab, status)
  call check_close(tab%value, [0.6_dp, 0.4_dp], 1.0e-12_dp, 'prior marginal node 1', failures)
  call query_marginal(net, 2, tab, status)
  call check_close(tab%value, [0.5_dp, 0.5_dp], 1.0e-12_dp, 'prior marginal node 2', failures)
  call query_marginal(net, 3, tab, status)
  call check_close(tab%value, [0.65_dp, 0.35_dp], 1.0e-12_dp, 'prior marginal node 3', failures)

  call query_joint(net, [1, 3], joint, status)
  call check(status == grain_ok, 'nonclique joint query', failures)
  call check_close(joint%value, [0.45_dp, 0.20_dp, 0.15_dp, 0.20_dp], 1.0e-12_dp, &
                   'reconstructed joint P(1,3)', failures)

  call set_hard_evidence(net, 3, 1, status)
  call check(status == grain_ok, 'hard evidence and propagation', failures)
  pe = probability_of_evidence(net, status)
  call check_close_scalar(pe, 0.65_dp, 1.0e-12_dp, 'hard evidence probability', failures)
  call query_marginal(net, 1, tab, status)
  call check_close(tab%value, [0.6923076923076923_dp, 0.3076923076923077_dp], &
                   1.0e-12_dp, 'posterior node 1', failures)
  call query_joint(net, [1, 2], joint, status)
  call check_close(joint%value, [0.5815384615384616_dp, 0.1107692307692308_dp, &
                                 0.1107692307692308_dp, 0.1969230769230769_dp], &
                   1.0e-12_dp, 'posterior joint nodes 1 and 2', failures)
  call query_conditional(net, [1, 2], cond, status)
  call check_close(cond%value, [0.84_dp, 0.16_dp, 0.36_dp, 0.64_dp], &
                   1.0e-12_dp, 'conditional node 1 given node 2 and evidence', failures)

  call simulate_network(net, 100, sim1, status, seed=2468)
  call simulate_network(net, 100, sim2, status, seed=2468)
  call check(all(sim1 == sim2), 'portable seeded simulation reproducibility', failures)
  call check(all(sim1(:, 3) == 1), 'hard evidence respected in simulation', failures)
  call check(all(sim1 >= 1 .and. sim1 <= 2), 'simulation state bounds', failures)

  call clear_evidence(net, status, propagate=.false.)
  call set_soft_evidence(net, 2, [0.25_dp, 1.0_dp], status)
  pe = probability_of_evidence(net, status)
  call check_close_scalar(pe, 0.625_dp, 1.0e-12_dp, 'soft evidence probability', failures)
  call query_marginal(net, 3, tab, status)
  call check_close(tab%value, [0.5_dp, 0.5_dp], 1.0e-12_dp, 'soft-evidence posterior', failures)

  call clear_evidence(net, status)
  call replace_network_cpt(net, 1, [0.5_dp, 0.5_dp], status, propagate=.true.)
  call query_marginal(net, 1, tab, status)
  call check_close(tab%value, [0.5_dp, 0.5_dp], 1.0e-12_dp, 'replace CPT without retriangulation', failures)

  tab = logical_and_cpt(1, 2, 3)
  call check_close(tab%value, [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], &
                   0.0_dp, 'logical AND CPT', failures)
  tab = logical_or_cpt(1, 2, 3)
  call check_close(tab%value, [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
                   0.0_dp, 'logical OR CPT', failures)
  call check_close_scalar(mendel_probability(2, 2, 2, 2), 0.5_dp, 1.0e-12_dp, &
                          'Mendel heterozygous child probability', failures)

  data = reshape([ &
      1, 1, 1,  1, 1, 1,  1, 1, 2,  1, 2, 1,  1, 2, 2, &
      2, 1, 1,  2, 2, 1,  2, 2, 2,  2, 2, 2,  2, 2, 2], shape(data), order=[2, 1])
  call estimate_cpts_from_data(data, net%dag, [2, 2, 2], estimated, status, smooth=1.0_dp)
  call check(status == grain_ok, 'estimate CPTs from categorical data', failures)
  do i = 1, size(estimated)
    call check(valid_table(estimated(i)), 'estimated CPT is valid', failures)
    call check_first_normalized(estimated(i), 'estimated CPT columns normalized', failures)
  end do

  call estimate_clique_marginals(data, net%rip, [2, 2, 2], marginals, status, smooth=1.0_dp)
  call check(status == grain_ok, 'estimate clique marginals', failures)
  call marginals_to_potentials(marginals, net%rip, potentials, status)
  call check(status == grain_ok, 'marginals to potentials', failures)
  call potentials_to_marginals(potentials, net%rip, marginals2, status)
  call check(status == grain_ok, 'potentials to marginals', failures)
  do i = 1, size(marginals)
    call check(table_equal(marginals(i), marginals2(i), 1.0e-12_dp), &
               'marginal/potential round trip', failures)
  end do

  if (failures /= 0) then
    write(*, '(a,i0)') 'FAILED tests: ', failures
    error stop 1
  end if
  write(*, '(a)') 'All gRain deterministic tests passed.'

contains

  subroutine check(condition, label, failures)
    logical, value :: condition !! Boolean condition that must be true for the test to pass.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.
    integer, intent(inout) :: failures !! Running failure count incremented when `condition` is false.

    if (.not. condition) then
      failures = failures + 1
      write(*, '(a)') 'FAIL: ' // trim(label)
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, label, failures)
    real(dp), intent(in) :: actual(:) !! Computed vector result.
    real(dp), intent(in) :: expected(:) !! Independently specified reference vector.
    real(dp), value :: tolerance !! Maximum permitted absolute elementwise error.
    character(len=*), intent(in) :: label !! Human-readable test label printed on mismatch.
    integer, intent(inout) :: failures !! Running failure count incremented on size or tolerance mismatch.

    if (size(actual) /= size(expected)) then
      failures = failures + 1
      write(*, '(a)') 'FAIL(size): ' // trim(label)
    else if (any(abs(actual - expected) > tolerance)) then
      failures = failures + 1
      write(*, '(a)') 'FAIL(value): ' // trim(label)
      write(*, '(a,*(1x,es14.6))') ' actual  ', actual
      write(*, '(a,*(1x,es14.6))') ' expected', expected
    end if
  end subroutine check_close

  subroutine check_close_scalar(actual, expected, tolerance, label, failures)
    real(dp), value :: actual !! Computed scalar result.
    real(dp), value :: expected !! Independent reference scalar.
    real(dp), value :: tolerance !! Maximum permitted absolute error.
    character(len=*), intent(in) :: label !! Human-readable test label printed on mismatch.
    integer, intent(inout) :: failures !! Running failure count incremented on tolerance mismatch.

    if (abs(actual - expected) > tolerance) then
      failures = failures + 1
      write(*, '(a,2(1x,es14.6))') 'FAIL: ' // trim(label), actual, expected
    end if
  end subroutine check_close_scalar

  subroutine check_first_normalized(tab, label, failures)
    type(table_t), intent(in) :: tab !! CPT table checked for unit sums over its first/child dimension.
    character(len=*), intent(in) :: label !! Human-readable test label printed on mismatch.
    integer, intent(inout) :: failures !! Running failure count incremented when any conditional column is not normalized.
    integer :: block
    integer :: first_dim
    integer :: offset

    first_dim = tab%dim(1)
    do block = 0, size(tab%value) / first_dim - 1
      offset = block * first_dim
      if (abs(sum(tab%value(offset + 1:offset + first_dim)) - 1.0_dp) > 1.0e-12_dp) then
        failures = failures + 1
        write(*, '(a)') 'FAIL: ' // trim(label)
        return
      end if
    end do
  end subroutine check_first_normalized

end program test_grain
