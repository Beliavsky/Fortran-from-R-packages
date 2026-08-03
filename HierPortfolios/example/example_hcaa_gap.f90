! SPDX-License-Identifier: GPL-2.0-only
program example_hcaa_gap
  use hierportfolios, only: dp, portfolio_result, HCAA_Portfolio
  implicit none

  real(dp) :: covar(6, 6)
  type(portfolio_result) :: result
  integer :: i

  call example_covariance(covar)
  call HCAA_Portfolio(covar, result, linkage='ward', gap_references=25, seed=2025)

  write (*, '(a,i0)') 'Selected clusters: ', result%n_clusters
  do i = 1, size(result%weights)
    write (*, '(i0,2x,f10.6,2x,a,i0)') i, result%weights(i), 'cluster ', result%clusters(i)
  end do

contains

  subroutine example_covariance(c)
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
  end subroutine example_covariance

end program example_hcaa_gap
