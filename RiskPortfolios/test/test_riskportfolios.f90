! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_riskportfolios
   use riskportfolios
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none

   integer, parameter :: t = 180, n = 8
   real(dp) :: rets(t, n), small(4, 2)
   real(dp), allocatable :: mu(:), semidev(:), sigma(:, :), w(:)
   real(dp) :: expected_cov(2, 2), rc(n), tol
   type(portfolio_control) :: ctr
   integer :: i, j, method, info, failures

   failures = 0
   tol = 1.0e-8_dp

   do j = 1, n
      do i = 1, t
         rets(i, j) = 0.002_dp * real(j, dp) + &
            0.012_dp * sin(0.071_dp * real(i * (j + 1), dp)) + &
            0.007_dp * cos(0.037_dp * real(i + 3 * j, dp)) + &
            0.003_dp * sin(0.013_dp * real(i * i + j, dp))
      end do
   end do

   small(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   small(:, 2) = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp]
   expected_cov = reshape([5.0_dp / 3.0_dp, 1.0_dp, 1.0_dp, &
      5.0_dp / 3.0_dp], [2, 2])
   sigma = naive_covariance(small)
   call check(maxval(abs(sigma - expected_cov)) < tol, 'sample covariance reference')

   call mean_estimation(rets, mu, MEAN_NAIVE, info=info)
   call check(info == 0 .and. size(mu) == n, 'naive mean')
   call mean_estimation(rets, mu, MEAN_EWMA, lambda=0.90_dp, info=info)
   call check(info == 0 .and. all_finite(mu), 'EWMA mean')
   call mean_estimation(rets, mu, MEAN_BAYES_STEIN, info=info)
   call check(info == 0 .and. all_finite(mu), 'Bayes-Stein mean')
   call mean_estimation(rets, mu, MEAN_MARTELLINI, info=info)
   call check(info == 0 .and. all(mu > 0.0_dp), 'Martellini mean')

   call semideviation_estimation(rets, semidev, SEMIDEV_NAIVE, info=info)
   call check(info == 0 .and. all(semidev >= 0.0_dp), 'naive semideviation')
   call semideviation_estimation(rets, semidev, SEMIDEV_EWMA, 0.91_dp, info)
   call check(info == 0 .and. all(semidev >= 0.0_dp), 'EWMA semideviation')

   do method = COV_NAIVE, COV_BAYES_STEIN
      if (method == COV_FACTOR) then
         call covariance_estimation(rets, sigma, method, n_factors=3, info=info)
      else
         call covariance_estimation(rets, sigma, method, lambda=0.92_dp, info=info)
      end if
      call check(info == 0, 'covariance method status ' // int_string(method))
      call check(size(sigma, 1) == n .and. size(sigma, 2) == n, &
         'covariance method dimensions ' // int_string(method))
      call check(maxval(abs(sigma - transpose(sigma))) < 1.0e-10_dp, &
         'covariance method symmetry ' // int_string(method))
      call check(all(diagonal(sigma) > 0.0_dp), &
         'covariance method positive diagonal ' // int_string(method))
      call check(all_finite_matrix(sigma), &
         'covariance method finite ' // int_string(method))
      call check(is_positive_definite(sigma, 1.0e-14_dp), &
         'covariance method positive definite ' // int_string(method))
   end do

   call mean_estimation(rets, mu, MEAN_NAIVE, info=info)
   call covariance_estimation(rets, sigma, COV_LEDOIT_WOLF, info=info)
   call semideviation_estimation(rets, semidev, SEMIDEV_NAIVE, info=info)

   do method = PORT_MEAN_VARIANCE, PORT_MAXIMUM_DECORRELATION
      call optimal_portfolio(sigma, w, method, mu, semidev, info=info)
      call check(info == 0, 'unconstrained portfolio status ' // int_string(method))
      call check(abs(sum(w) - 1.0_dp) < 1.0e-7_dp, &
         'unconstrained portfolio sum ' // int_string(method))
      call check(all_finite(w), 'unconstrained portfolio finite ' // int_string(method))
   end do

   ctr = portfolio_control()
   ctr%constraint = CONSTRAINT_LONG_ONLY
   ctr%max_iterations = 8000
   ctr%tolerance = 1.0e-9_dp
   do method = PORT_MEAN_VARIANCE, PORT_MAXIMUM_DECORRELATION
      call optimal_portfolio(sigma, w, method, mu, semidev, ctr, info)
      call check(info == 0, 'long-only portfolio status ' // int_string(method))
      call check(abs(sum(w) - 1.0_dp) < 1.0e-7_dp, &
         'long-only portfolio sum ' // int_string(method))
      call check(minval(w) >= -1.0e-9_dp, &
         'long-only portfolio bounds ' // int_string(method))
   end do

   call optimal_portfolio(sigma, w, PORT_EQUAL_RISK_CONTRIBUTION, &
      mu, semidev, ctr, info)
   rc = portfolio_risk_contributions(sigma, w)
   call check(maxval(abs(rc - 1.0_dp / real(n, dp))) < 3.0e-4_dp, &
      'ERC equal contributions')

   ctr = portfolio_control()
   ctr%constraint = CONSTRAINT_USER
   allocate(ctr%lower_bounds(n), ctr%upper_bounds(n))
   ctr%lower_bounds = 0.02_dp
   ctr%upper_bounds = 0.40_dp
   ctr%max_iterations = 8000
   call optimal_portfolio(sigma, w, PORT_MINIMUM_VARIANCE, &
      mu, semidev, ctr, info)
   call check(info == 0, 'user-bounded minimum variance status')
   call check(minval(w) >= 0.02_dp - 1.0e-8_dp .and. &
      maxval(w) <= 0.40_dp + 1.0e-8_dp, 'user bounds respected')

   ctr = portfolio_control()
   ctr%constraint = CONSTRAINT_GROSS
   ctr%gross_limit = 1.20_dp
   ctr%max_iterations = 8000
   call optimal_portfolio(sigma, w, PORT_MEAN_VARIANCE, mu, semidev, ctr, info)
   call check(info == 0, 'gross-constrained mean-variance status')
   call check(abs(sum(w) - 1.0_dp) < 1.0e-7_dp, 'gross portfolio sum')
   call check(sum(abs(w)) <= 1.20_dp + 1.0e-7_dp, 'gross exposure respected')

   if (failures > 0) then
      write(*, '(a,1x,i0)') 'FAILED tests:', failures
      error stop 1
   end if
   write(*, '(a)') 'All RiskPortfolios tests passed.'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // label
      end if
   end subroutine check

   pure function diagonal(a) result(d)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: d(min(size(a, 1), size(a, 2)))
      integer :: k
      do k = 1, size(d)
         d(k) = a(k, k)
      end do
   end function diagonal

   pure logical function all_finite(x) result(ok)
      real(dp), intent(in) :: x(:)
      ok = all(ieee_is_finite(x))
   end function all_finite

   pure logical function all_finite_matrix(x) result(ok)
      real(dp), intent(in) :: x(:, :)
      ok = all(ieee_is_finite(x))
   end function all_finite_matrix

   function int_string(x) result(text)
      integer, intent(in) :: x
      character(len=16) :: text
      write(text, '(i0)') x
   end function int_string

end program test_riskportfolios
