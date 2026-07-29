! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
program demo_riskportfolios
   use riskportfolios
   implicit none

   integer, parameter :: t = 252, n = 6
   real(dp) :: rets(t, n)
   real(dp), allocatable :: mu(:), sigma(:, :), semidev(:), weights(:)
   type(portfolio_control) :: control
   integer :: i, j, info

   do j = 1, n
      do i = 1, t
         rets(i, j) = 0.0002_dp * real(j, dp) + &
            0.012_dp * sin(0.03_dp * real(i * (j + 1), dp)) + &
            0.006_dp * cos(0.07_dp * real(i + j, dp))
      end do
   end do

   call mean_estimation(rets, mu, MEAN_BAYES_STEIN, info=info)
   call covariance_estimation(rets, sigma, COV_LEDOIT_WOLF, info=info)
   call semideviation_estimation(rets, semidev, SEMIDEV_EWMA, 0.94_dp, info)

   control%constraint = CONSTRAINT_LONG_ONLY
   call optimal_portfolio(sigma, weights, PORT_EQUAL_RISK_CONTRIBUTION, &
      mu, semidev, control, info)

   write(*, '(a)') 'Long-only equal-risk-contribution portfolio'
   do i = 1, n
      write(*, '(a,i0,a,f10.6)') 'asset ', i, ': ', weights(i)
   end do
   write(*, '(a,f10.6)') 'sum: ', sum(weights)
end program demo_riskportfolios
