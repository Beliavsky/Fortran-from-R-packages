program test_ecp_algorithms
  use ecp_api, only: dp, cp_result, divisive_result
  use ecp_api, only: e_cp3o, e_cp3o_delta, ks_cp3o, ks_cp3o_delta
  use ecp_api, only: e_divisive, kcpa
  implicit none

  real(dp) :: x6(6, 1), x8(8, 1), x9(9, 1)
  type(cp_result) :: cp
  type(divisive_result) :: div
  integer, allocatable :: bounds(:), all_rows(:)

  x6(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 10.0_dp]
  cp = e_cp3o(x6, k=1, minsize=2, alpha=1.0_dp)
  call assert_true(cp%number == 1, 'e.cp3o number')
  call assert_true(size(cp%estimates) == 1 .and. cp%estimates(1) == 4, 'e.cp3o one change')

  cp = ks_cp3o(x6, k=1, minsize=2)
  call assert_true(size(cp%estimates) == 1 .and. cp%estimates(1) == 4, 'ks.cp3o one change')

  cp = ks_cp3o_delta(x6, k=1, minsize=2)
  call assert_true(size(cp%estimates) == 1 .and. cp%estimates(1) == 4, 'ks.cp3o_delta one change')

  x8(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp]
  cp = e_cp3o_delta(x8, k=1, delta=2, alpha=1.0_dp)
  call assert_true(size(cp%estimates) == 1 .and. cp%estimates(1) == 5, 'e.cp3o_delta one change')

  x9(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 5.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  cp = e_cp3o(x9, k=2, minsize=2, alpha=1.0_dp)
  call assert_true(all(cp%cp_loc(2, 1:2) == [4, 7]), 'e.cp3o two-change path')
  cp = ks_cp3o(x9, k=2, minsize=2)
  call assert_true(all(cp%cp_loc(2, 1:2) == [4, 7]), 'ks.cp3o two-change path')

  div = e_divisive(x9, k=2, min_size=2, alpha=1.0_dp)
  call assert_true(div%k_hat == 3, 'e.divisive cluster count')
  call assert_true(all(div%estimates == [1, 4, 7, 10]), 'e.divisive boundaries')
  call assert_true(all(div%cluster == [1, 1, 1, 2, 2, 2, 3, 3, 3]), 'e.divisive cluster labels')


  call seed_rng(12345)
  div = e_divisive(x6, sig_level=0.05_dp, r=9, min_size=2, alpha=1.0_dp)
  call assert_true(size(div%p_values) >= 1, 'e.divisive permutation p-values')
  call assert_true(all(div%p_values >= 0.0_dp .and. div%p_values <= 1.0_dp), 'e.divisive p-value range')
  call assert_true(all(div%permutations == 9), 'e.divisive permutation counts')

  bounds = kcpa(x6, l=2, c=1.0_dp)
  call assert_true(all(bounds == [1, 4, 7]), 'kcpa boundaries')
  all_rows = [1, 2, 3, 4, 5, 6]
  bounds = kcpa(x6, l=2, c=1.0_dp, bandwidth_rows=all_rows)
  call assert_true(all(bounds == [1, 4, 7]), 'kcpa explicit bandwidth rows')

  print '(a)', 'All ecp algorithm tests passed.'

contains


  subroutine seed_rng(seed_value)
    integer, intent(in) :: seed_value !! Scalar seed expanded deterministically across the compiler RNG state.
    integer, allocatable :: seed(:)
    integer :: i, nseed

    call random_seed(size=nseed)
    allocate(seed(nseed))
    do i = 1, nseed
      seed(i) = modulo(seed_value + 104729*i, huge(1) - 1)
    end do
    call random_seed(put=seed)
  end subroutine seed_rng

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition !! Logical assertion expected to be true.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (.not. condition) then
      write (*, '(a)') 'FAILED '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_ecp_algorithms
