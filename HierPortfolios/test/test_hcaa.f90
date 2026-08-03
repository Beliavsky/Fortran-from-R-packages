! SPDX-License-Identifier: GPL-2.0-only
program test_hcaa
  use hierportfolios, only: dp, portfolio_result, HCAA_Portfolio
  implicit none

  real(dp) :: covar(6, 6)
  type(portfolio_result) :: fixed, automatic, repeated

  call make_block_covariance(covar)
  call HCAA_Portfolio(covar, fixed, clusters=2)
  call check(fixed%ok(), 'HCAA fixed status')
  call check(fixed%n_clusters == 2, 'HCAA fixed cluster count')
  call check(abs(sum(fixed%weights) - 1.0_dp) < 1.0e-12_dp, 'HCAA sum')
  call check(all(abs(fixed%weights - 1.0_dp / 6.0_dp) < 1.0e-12_dp), &
    'HCAA equal within two balanced clusters')
  call check(all(fixed%clusters(1:3) == fixed%clusters(1)), 'first block cluster')
  call check(all(fixed%clusters(4:6) == fixed%clusters(4)), 'second block cluster')
  call check(fixed%clusters(1) /= fixed%clusters(4), 'separate blocks')

  call HCAA_Portfolio(covar, automatic, gap_references=12, seed=777)
  call HCAA_Portfolio(covar, repeated, gap_references=12, seed=777)
  call check(automatic%ok(), 'HCAA gap status')
  call check(automatic%n_clusters >= 2 .and. automatic%n_clusters <= 3, &
    'HCAA gap cluster range')
  call check(maxval(abs(automatic%weights - repeated%weights)) < 1.0e-14_dp, &
    'HCAA deterministic gap selection')
  call check(size(automatic%gap) == 3, 'HCAA gap history')

  print '(a)', 'test_hcaa: PASS'

contains

  subroutine make_block_covariance(c)
    real(dp), intent(out) :: c(:, :)
    integer :: i, j

    c = 0.0_dp
    do i = 1, 6
      c(i, i) = 0.02_dp + 0.005_dp * real(i, dp)
    end do
    do j = 1, 3
      do i = 1, 3
        if (i /= j) c(i, j) = 0.60_dp * sqrt(c(i, i) * c(j, j))
      end do
    end do
    do j = 4, 6
      do i = 4, 6
        if (i /= j) c(i, j) = 0.50_dp * sqrt(c(i, i) * c(j, j))
      end do
    end do
    do j = 4, 6
      do i = 1, 3
        c(i, j) = 0.05_dp * sqrt(c(i, i) * c(j, j))
        c(j, i) = c(i, j)
      end do
    end do
  end subroutine make_block_covariance

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write (*, '(a)') 'FAIL: '//label
      error stop 1
    end if
  end subroutine check

end program test_hcaa
