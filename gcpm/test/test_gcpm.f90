! SPDX-License-Identifier: GPL-2.0-only
program test_gcpm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use gcpm
   implicit none

   integer :: failures

   failures = 0
   call test_normal_functions(failures)
   call test_analytical_distribution(failures)
   call test_csv_reader(failures)
   call test_simulation_links(failures)

   if (failures > 0) then
      write(*, '(a,i0)') 'FAILURES: ', failures
      error stop 1
   end if
   print '(a)', 'All GCPM tests passed.'

contains

   subroutine test_normal_functions(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: probabilities(5) = [0.001_dp, 0.025_dp, 0.5_dp, 0.975_dp, 0.999_dp]
      real(dp) :: recovered(5)

      recovered = normal_cdf(normal_quantile(probabilities))
      call check(maxval(abs(recovered - probabilities)) < 1.0e-10_dp, &
                 'normal CDF/quantile round trip', failures)
   end subroutine test_normal_functions

   subroutine test_analytical_distribution(failures)
      integer, intent(inout) :: failures
      type(credit_portfolio) :: portfolio
      type(gcpm_model) :: model
      type(loss_distribution) :: distribution
      type(risk_measures) :: measures
      real(dp), allocatable :: contributions(:,:), sd_contribution(:)
      real(dp), parameter :: alpha(2) = [0.95_dp, 0.99_dp]
      integer :: status
      character(len=:), allocatable :: message

      call allocate_portfolio(portfolio, 4, 2)
      portfolio%ead = [100000.0_dp, 80000.0_dp, 120000.0_dp, 70000.0_dp]
      portfolio%lgd = [0.4_dp, 0.5_dp, 0.3_dp, 0.6_dp]
      portfolio%pd = [0.02_dp, 0.03_dp, 0.015_dp, 0.025_dp]
      portfolio%default_kind = default_poisson
      portfolio%weight(1, :) = [0.7_dp, 0.0_dp]
      portfolio%weight(2, :) = [0.5_dp, 0.2_dp]
      portfolio%weight(3, :) = [0.0_dp, 0.8_dp]
      portfolio%weight(4, :) = [0.3_dp, 0.4_dp]

      model%model_kind = model_analytical_crp
      model%link_kind = link_crp
      model%loss_unit = 1000.0_dp
      model%alpha_max = 0.99999_dp
      allocate(model%sector_variance(2))
      model%sector_variance = [0.8_dp, 1.2_dp]

      call calculate_portfolio_statistics(portfolio, model, status, message)
      call check(status == 0, 'portfolio statistics status', failures)
      if (status /= 0) return
      call analytical_creditrisk_plus(portfolio, model, distribution, status, message)
      call check(status == 0 .or. status == 4, 'analytical recursion status', failures)
      call check(all(distribution%pdf >= -1.0e-12_dp), 'nonnegative PDF', failures)
      call check(all(distribution%cdf(2:) >= distribution%cdf(:size(distribution%cdf)-1)), &
                 'monotone CDF', failures)
      call check(distribution%reached_alpha >= 0.999_dp, 'high CDF coverage', failures)
      call check(abs(distribution%expected_loss / portfolio%analytical_expected_loss - 1.0_dp) < 0.03_dp, &
                 'distribution expected loss near analytical expected loss', failures)

      call calculate_risk_measures(distribution, portfolio, alpha, measures, status, message)
      call check(status == 0, 'risk-measure status', failures)
      call check(measures%var(2) >= measures%var(1), 'VaR monotonic in alpha', failures)
      call check(all(measures%expected_shortfall >= measures%var - 1.0e-8_dp), &
                 'expected shortfall at least VaR', failures)

      call standard_deviation_contributions(portfolio, model, sd_contribution, status, message)
      call check(status == 0, 'standard-deviation contribution status', failures)
      call check(abs(sum(sd_contribution) / portfolio%analytical_sd - 1.0_dp) < 1.0e-10_dp, &
                 'standard-deviation contributions add up', failures)

      call analytical_risk_contributions(portfolio, model, distribution, 0.95_dp, &
                                         contributions, status, message)
      call check(status == 0, 'analytical risk-contribution status', failures)
      if (status == 0) then
         call check(all(ieee_is_finite(contributions)), &
                    'finite analytical risk contributions', failures)
      end if
   end subroutine test_analytical_distribution

   subroutine test_csv_reader(failures)
      integer, intent(inout) :: failures
      type(credit_portfolio) :: portfolio
      integer :: status
      character(len=:), allocatable :: message

      call read_gcpm_portfolio('data/sample_portfolio.csv', portfolio, status, message)
      call check(status == 0, 'CSV reader status', failures)
      if (status == 0) then
         call check(portfolio%n_counterparties == 20, 'CSV row count', failures)
         call check(portfolio%n_sectors == 3, 'CSV sector count', failures)
         call check(abs(portfolio%ead(1) - 364942.6056_dp) < 1.0e-6_dp, &
                    'CSV decimal-comma parsing', failures)
      end if
   end subroutine test_csv_reader

   subroutine test_simulation_links(failures)
      integer, intent(inout) :: failures
      type(credit_portfolio) :: portfolio
      type(gcpm_model) :: model
      type(loss_distribution) :: distribution
      real(dp), allocatable :: factors(:,:), contributions(:,:)
      integer :: i, status
      character(len=:), allocatable :: message

      call allocate_portfolio(portfolio, 4, 2)
      portfolio%ead = [100000.0_dp, 80000.0_dp, 120000.0_dp, 70000.0_dp]
      portfolio%lgd = [0.4_dp, 0.5_dp, 0.3_dp, 0.6_dp]
      portfolio%pd = [0.02_dp, 0.03_dp, 0.015_dp, 0.025_dp]
      portfolio%default_kind = default_bernoulli
      portfolio%weight(1, :) = [0.25_dp, 0.00_dp]
      portfolio%weight(2, :) = [0.20_dp, 0.10_dp]
      portfolio%weight(3, :) = [0.00_dp, 0.30_dp]
      portfolio%weight(4, :) = [0.15_dp, 0.15_dp]

      model%model_kind = model_simulation
      model%link_kind = link_cm
      model%loss_unit = 1000.0_dp
      model%n_simulations = 6000
      model%seed = 24680
      model%loss_threshold = 0.0_dp
      model%max_stored_scenarios = model%n_simulations
      call calculate_portfolio_statistics(portfolio, model, status, message)
      call check(status == 0, 'simulation portfolio statistics status', failures)
      if (status /= 0) return

      allocate(factors(model%n_simulations, 2))
      call seed_random_number(model%seed)
      do i = 1, model%n_simulations
         factors(i, 1) = random_standard_normal()
         factors(i, 2) = random_standard_normal()
      end do
      call simulate_loss_distribution(portfolio, model, factors, distribution, status, message)
      call check(status == 0, 'CreditMetrics simulation status', failures)
      call check(abs(distribution%cdf(size(distribution%cdf)) - 1.0_dp) < 1.0e-12_dp, &
                 'simulation CDF sums to one', failures)
      call check(distribution%contributions_available, &
                 'simulation contribution scenarios stored', failures)
      call simulation_risk_contributions(portfolio, distribution, 0.95_dp, &
                                         contributions, status, message)
      call check(status == 0, 'simulation risk-contribution status', failures)
      if (status == 0) then
         call check(all(ieee_is_finite(contributions)), &
                    'finite simulation risk contributions', failures)
      end if
   end subroutine test_simulation_links

   subroutine check(condition, label, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer, intent(inout) :: failures

      if (condition) then
         write(*, '(a)') 'ok: ' // trim(label)
      else
         write(*, '(a)') 'not ok: ' // trim(label)
         failures = failures + 1
      end if
   end subroutine check

end program test_gcpm
