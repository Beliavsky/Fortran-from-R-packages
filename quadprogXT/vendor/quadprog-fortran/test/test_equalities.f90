! SPDX-License-Identifier: GPL-2.0-or-later
program test_equalities
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, qp_success, solve_qp
  use quadprog_test_support, only: check, check_close
  implicit none

  integer :: nblock, n, q, i
  real(dp), allocatable :: d(:, :), dvec(:), a(:, :), b(:), expected(:)
  type(qp_result) :: res

  do nblock = 1, 10
    n = 2 * nblock
    q = 2 + 2 * nblock
    allocate(d(n, n), dvec(n), a(n, q), b(q), expected(n))
    d = 0.0_dp
    do i = 1, n
      d(i, i) = 1.0_dp
    end do
    dvec = 0.0_dp
    a = 0.0_dp
    do i = 1, nblock
      a(2 * i - 1, 1) = 1.0_dp
      a(2 * i, 1) = -1.0_dp
    end do
    a(:, 2) = 1.0_dp
    do i = 1, nblock
      a(2 * i - 1, 2 + 2 * i - 1) = -1.0_dp
      a(2 * i, 2 + 2 * i) = 1.0_dp
    end do
    b(1:2) = 1.0_dp
    do i = 1, nblock
      b(2 + 2 * i - 1) = -1.0_dp
      b(2 + 2 * i) = 0.0_dp
    end do
    expected = 0.0_dp
    expected(1:n:2) = 1.0_dp / real(nblock, dp)

    res = solve_qp(d, dvec, a, b, meq=2)
    call check(res%status == qp_success, 'Talbot-Katz status')
    call check(maxval(abs(res%solution - expected)) < 2.0e-11_dp, &
      'Talbot-Katz solution')
    call check_close(res%value, 0.5_dp / real(nblock, dp), &
      2.0e-11_dp, 'Talbot-Katz objective')
    call check(all(res%active_set == [1, 2]), 'Talbot-Katz active set')
    deallocate(d, dvec, a, b, expected)
  end do

  write(*, '(a)') 'test_equalities: PASS'
end program test_equalities
