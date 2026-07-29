! SPDX-License-Identifier: GPL-2.0-only
program demo_gcpm
   use gcpm
   implicit none

   type(credit_portfolio) :: portfolio
   type(gcpm_model) :: model
   type(loss_distribution) :: analytical_dist, simulated_dist
   type(risk_measures) :: measures
   real(dp), allocatable :: factors(:,:)
   real(dp), parameter :: alpha(2) = [0.95_dp, 0.99_dp]
   integer :: i, k, status
   character(len=:), allocatable :: message

   call allocate_portfolio(portfolio, 5, 2)
   portfolio%number = [1, 2, 3, 4, 5]
   portfolio%name = ['Firm A', 'Firm B', 'Firm C', 'Firm D', 'Firm E']
   portfolio%sector_name = ['Industry', 'Property']
   portfolio%ead = [100000.0_dp, 80000.0_dp, 120000.0_dp, 60000.0_dp, 90000.0_dp]
   portfolio%lgd = [0.45_dp, 0.40_dp, 0.50_dp, 0.35_dp, 0.55_dp]
   portfolio%pd = [0.020_dp, 0.030_dp, 0.015_dp, 0.040_dp, 0.025_dp]
   portfolio%default_kind = default_poisson
   portfolio%weight(1, :) = [0.70_dp, 0.00_dp]
   portfolio%weight(2, :) = [0.60_dp, 0.10_dp]
   portfolio%weight(3, :) = [0.00_dp, 0.80_dp]
   portfolio%weight(4, :) = [0.20_dp, 0.50_dp]
   portfolio%weight(5, :) = [0.40_dp, 0.30_dp]

   model%model_kind = model_analytical_crp
   model%link_kind = link_crp
   model%loss_unit = 1000.0_dp
   model%alpha_max = 0.9999_dp
   allocate(model%sector_variance(2))
   model%sector_variance = [0.75_dp, 1.10_dp]

   call calculate_portfolio_statistics(portfolio, model, status, message)
   call stop_on_error(status, message)
   call analytical_creditrisk_plus(portfolio, model, analytical_dist, status, message)
   if (status /= 0 .and. status /= 4) call stop_on_error(status, message)
   call calculate_risk_measures(analytical_dist, portfolio, alpha, measures, status, message)
   call stop_on_error(status, message)

   print '(a)', 'Analytical CreditRisk+ results'
   print '(a,f12.2)', 'Expected loss (portfolio): ', portfolio%analytical_expected_loss
   print '(a,f12.2)', 'Expected loss (distribution): ', analytical_dist%expected_loss
   print '(a,f12.2)', 'Analytical standard deviation: ', portfolio%analytical_sd
   do i = 1, size(alpha)
      print '(a,f7.4,a,f12.2,a,f12.2,a,f12.2)', 'alpha=', alpha(i), &
         '  VaR=', measures%var(i), '  EC=', measures%economic_capital(i), &
         '  ES=', measures%expected_shortfall(i)
   end do

   model%model_kind = model_simulation
   model%n_simulations = 30000
   model%seed = 20260727
   model%loss_threshold = 0.0_dp
   model%max_stored_scenarios = model%n_simulations
   allocate(factors(model%n_simulations, portfolio%n_sectors))
   call seed_random_number(model%seed)
   do k = 1, portfolio%n_sectors
      do i = 1, model%n_simulations
         factors(i, k) = random_gamma(1.0_dp / model%sector_variance(k), &
                                      model%sector_variance(k))
      end do
   end do

   call simulate_loss_distribution(portfolio, model, factors, simulated_dist, &
                                   status, message)
   if (status /= 0 .and. status /= 9 .and. status /= 10) call stop_on_error(status, message)
   call calculate_risk_measures(simulated_dist, portfolio, alpha, measures, status, message)
   call stop_on_error(status, message)

   print '(/,a)', 'Monte Carlo CRP results'
   print '(a,f12.2)', 'Expected loss (simulation): ', simulated_dist%expected_loss
   print '(a,f12.2)', 'Standard deviation (simulation): ', simulated_dist%standard_deviation
   do i = 1, size(alpha)
      print '(a,f7.4,a,f12.2,a,f12.2,a,f12.2)', 'alpha=', alpha(i), &
         '  VaR=', measures%var(i), '  EC=', measures%economic_capital(i), &
         '  ES=', measures%expected_shortfall(i)
   end do

contains

   subroutine stop_on_error(code, text)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text

      if (code /= 0) then
         write(*, '(a,i0,2a)') 'Error ', code, ': ', trim(text)
         error stop 1
      end if
   end subroutine stop_on_error

end program demo_gcpm
