! SPDX-License-Identifier: GPL-2.0-only
program demo_hierportfolios
  use hierportfolios, only: dp, portfolio_result, HRP_Portfolio, HCAA_Portfolio, &
    HERC_Portfolio, DHRP_Portfolio
  implicit none

  real(dp) :: covar(6, 6), lower(6), upper(6)
  type(portfolio_result) :: result

  call example_covariance(covar)

  call HRP_Portfolio(covar, result)
  call show('HRP', result)

  call HCAA_Portfolio(covar, result, clusters=2)
  call show('HCAA', result)

  call HERC_Portfolio(covar, result, clusters=2)
  call show('HERC', result)

  lower = 0.05_dp
  upper = 0.22_dp
  call DHRP_Portfolio(covar, result, tau=0.5_dp, lb=lower, ub=upper)
  call show('DHRP', result)

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

  subroutine show(label, value)
    character(len=*), intent(in) :: label
    type(portfolio_result), intent(in) :: value

    write (*, '(a,1x,6(f9.5,1x))') label, value%weights
  end subroutine show

end program demo_hierportfolios
