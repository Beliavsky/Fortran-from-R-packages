program test_ecp_cp3o_parity
  use ecp_api, only: dp, cp_result
  use ecp_api, only: e_cp3o, e_cp3o_delta, ks_cp3o, ks_cp3o_delta
  implicit none

  integer, parameter :: n = 24
  real(dp) :: z(n, 2)
  type(cp_result) :: pruned, exhaustive
  integer :: case_id, i, total_pruned

  total_pruned = 0
  do case_id = 1, 12
    do i = 1, n
      z(i, 1) = sin(0.19_dp*real(i + 2*case_id, dp))
      z(i, 2) = cos(0.13_dp*real(2*i + case_id, dp))
      if (i > 6 + mod(case_id, 3)) z(i, 1) = z(i, 1) + 1.5_dp + 0.1_dp*real(case_id, dp)
      if (i > 12) z(i, 2) = z(i, 2) - 1.2_dp
      if (i > 18 - mod(case_id, 2)) z(i, 1) = z(i, 1) + 2.0_dp
    end do

    pruned = e_cp3o(z, k=3, minsize=3, alpha=1.0_dp, prune=.true.)
    exhaustive = e_cp3o(z, k=3, minsize=3, alpha=1.0_dp, prune=.false.)
    call assert_same(pruned, exhaustive, 'e.cp3o')
    total_pruned = total_pruned + pruned%candidates_pruned

    pruned = e_cp3o_delta(z, k=3, delta=2, alpha=1.0_dp, prune=.true.)
    exhaustive = e_cp3o_delta(z, k=3, delta=2, alpha=1.0_dp, prune=.false.)
    call assert_same(pruned, exhaustive, 'e.cp3o_delta')
    total_pruned = total_pruned + pruned%candidates_pruned

    pruned = ks_cp3o(z, k=3, minsize=3, prune=.true.)
    exhaustive = ks_cp3o(z, k=3, minsize=3, prune=.false.)
    call assert_same(pruned, exhaustive, 'ks.cp3o')
    total_pruned = total_pruned + pruned%candidates_pruned

    pruned = ks_cp3o_delta(z, k=3, minsize=3, prune=.true.)
    exhaustive = ks_cp3o_delta(z, k=3, minsize=3, prune=.false.)
    call assert_same(pruned, exhaustive, 'ks.cp3o_delta')
    total_pruned = total_pruned + pruned%candidates_pruned
  end do

  call assert_true(total_pruned > 0, 'candidate pruning exercised')
  print '(a,i0)', 'CP3O pruning parity cases passed, candidates pruned: ', total_pruned

contains

  subroutine assert_same(actual, reference, label)
    type(cp_result), intent(in) :: actual !! Result from the pruned CP3O implementation.
    type(cp_result), intent(in) :: reference !! Result from the exhaustive reference implementation.
    character(len=*), intent(in) :: label !! Human-readable CP3O variant name printed on failure.

    call assert_true(actual%number == reference%number, trim(label)//' selected count')
    call assert_true(size(actual%gof) == size(reference%gof), trim(label)//' GOF size')
    call assert_true(maxval(abs(actual%gof - reference%gof)) <= 1.0e-11_dp, trim(label)//' GOF values')
    call assert_true(all(actual%cp_loc == reference%cp_loc), trim(label)//' change-point paths')
    call assert_true(all(actual%estimates == reference%estimates), trim(label)//' selected estimates')
    call assert_true(actual%statistic_evaluations <= reference%statistic_evaluations, &
      trim(label)//' statistic evaluation count')
  end subroutine assert_same

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition !! Logical assertion expected to be true.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (.not. condition) then
      write (*, '(a)') 'FAILED '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_ecp_cp3o_parity
