program test_compact_normalize
  use quadprog_kinds, only: dp
  use quadprogxt, only: compact_constraints, normalized_constraints, &
    convert_to_compact, normalize_constraints
  use qpxt_test_support, only: assert_true, assert_close
  implicit none

  real(dp) :: a(3, 3), an(3, 4), b(4)
  type(compact_constraints) :: c
  type(normalized_constraints) :: nrm
  integer :: j

  a = reshape([ -4.0_dp, -3.0_dp, 0.0_dp, &
                 2.0_dp,  1.0_dp, 0.0_dp, &
                 0.0_dp, -2.0_dp, 1.0_dp ], [3, 3])
  c = convert_to_compact(a)
  call assert_true(c%succeeded(), 'convert_to_compact succeeds')
  call assert_true(all(c%aind(1, :) == [2, 2, 2]), 'nonzero counts')
  call assert_true(all(c%aind(2, :) == [1, 1, 2]), 'first indices')
  call assert_true(all(c%aind(3, :) == [2, 2, 3]), 'second indices')
  call assert_close(c%amat(1, 1), -4.0_dp, 0.0_dp, 'compact value')
  call assert_close(c%amat(2, 3), 1.0_dp, 0.0_dp, 'compact value')

  c = convert_to_compact(0.0_dp * a)
  call assert_true(.not. c%succeeded(), 'all-zero column is rejected')

  an = 0.0_dp
  an(1, 1) = 1.0_dp
  an(2, 2) = 2.0_dp
  an(3, 3) = 3.0_dp
  an(:, 4) = 11.0_dp
  b = 3.0_dp
  nrm = normalize_constraints(an, b)
  call assert_true(nrm%succeeded(), 'normalize succeeds')
  do j = 1, 4
    call assert_close(dot_product(nrm%amat(:, j), nrm%amat(:, j)), &
      1.0_dp, 2.0e-15_dp, 'normalized column norm')
  end do

  write(*, '(a)') 'PASS test_compact_normalize'
end program test_compact_normalize
