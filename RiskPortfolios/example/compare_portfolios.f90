! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
program compare_portfolios
   use riskportfolios
   implicit none

   integer, parameter :: t = 160, n = 5
   real(dp) :: rets(t, n)
   real(dp), allocatable :: mu(:), sigma(:, :), semidev(:), w(:)
   type(portfolio_control) :: control
   integer :: i, j, method, info
   character(len=30) :: names(7)

   names = [character(len=30) :: 'mean variance', 'minimum variance', &
      'inverse volatility', 'equal risk contribution', &
      'maximum diversification', 'risk efficient', 'maximum decorrelation']

   do j = 1, n
      do i = 1, t
         rets(i, j) = 0.0003_dp * real(j, dp) + &
            0.010_dp * sin(0.041_dp * real(i * (j + 1), dp)) + &
            0.004_dp * cos(0.027_dp * real(i + 2 * j, dp))
      end do
   end do

   call mean_estimation(rets, mu, MEAN_NAIVE, info=info)
   call covariance_estimation(rets, sigma, COV_COR_SHRINKAGE, info=info)
   call semideviation_estimation(rets, semidev, SEMIDEV_NAIVE, info=info)

   control%constraint = CONSTRAINT_LONG_ONLY
   do method = PORT_MEAN_VARIANCE, PORT_MAXIMUM_DECORRELATION
      call optimal_portfolio(sigma, w, method, mu, semidev, control, info)
      write(*, '(a)') trim(names(method))
      write(*, '(*(f10.6,1x))') w
   end do
end program compare_portfolios
