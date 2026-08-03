! SPDX-License-Identifier: GPL-2.0-only
program example_hrp_herc
  use hierportfolios, only: dp, portfolio_result, HRP_Portfolio, HERC_Portfolio
  implicit none

  real(dp) :: covar(6, 6)
  type(portfolio_result) :: hrp, herc

  call example_covariance(covar)
  call HRP_Portfolio(covar, hrp, linkage='single')
  call HERC_Portfolio(covar, herc, linkage='ward', clusters=2)

  call print_weights('HRP', hrp%weights)
  call print_weights('HERC', herc%weights)

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

  subroutine print_weights(label, weights)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: weights(:)
    integer :: i

    write (*, '(a)') label//' weights:'
    do i = 1, size(weights)
      write (*, '(2x,i0,2x,f10.6)') i, weights(i)
    end do
  end subroutine print_weights

end program example_hrp_herc
