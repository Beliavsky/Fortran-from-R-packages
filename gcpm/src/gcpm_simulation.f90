! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_simulation
   use gcpm_kinds, only: dp
   use gcpm_types, only: credit_portfolio, gcpm_model, loss_distribution, &
      default_bernoulli, link_crp, link_cm, model_simulation
   use gcpm_math, only: seed_random_number, random_poisson, normal_cdf, &
      normal_quantile, correlation_matrix, quadratic_form
   implicit none
   private

   public :: simulate_loss_distribution

contains

   subroutine simulate_loss_distribution(portfolio, model, factors, distribution, &
                                          status, message, likelihood_ratio)
      type(credit_portfolio), intent(in) :: portfolio
      type(gcpm_model), intent(in) :: model
      real(dp), intent(in) :: factors(:,:)
      type(loss_distribution), intent(out) :: distribution
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: likelihood_ratio(:)

      integer :: i, n, n_available, n_sim, row, n_defaults
      integer :: n_stored, max_stored, corr_status
      real(dp) :: conditional_pd, total_loss, u, systematic_variance, denom
      logical :: overflow
      real(dp), allocatable :: sigma(:,:), dd(:), cp_loss(:)
      real(dp), allocatable :: cp_store(:,:), sorted_loss(:), sorted_weight(:)
      integer, allocatable :: scenario_store(:)
      character(len=:), allocatable :: corr_message

      status = 0
      message = ''
      distribution%model_kind = model_simulation

      if (model%model_kind /= model_simulation) then
         status = 1
         message = 'simulate_loss_distribution requires model_simulation'
         return
      end if
      if (size(factors, 2) /= portfolio%n_sectors) then
         status = 2
         message = 'factor matrix has the wrong number of columns'
         return
      end if
      n_available = size(factors, 1)
      if (n_available <= 0) then
         status = 3
         message = 'at least one factor scenario is required'
         return
      end if

      n_sim = model%n_simulations
      if (n_sim <= 0) n_sim = n_available
      allocate(distribution%simulated_losses(n_sim), distribution%likelihood(n_sim))
      allocate(cp_loss(portfolio%n_counterparties))
      distribution%simulated_losses = 0.0_dp

      if (present(likelihood_ratio)) then
         if (size(likelihood_ratio) /= n_available) then
            status = 4
            message = 'likelihood_ratio has the wrong length'
            return
         end if
         if (any(likelihood_ratio < 0.0_dp) .or. sum(likelihood_ratio) <= 0.0_dp) then
            status = 5
            message = 'likelihood ratios must be nonnegative and have positive sum'
            return
         end if
         do n = 1, n_sim
            row = modulo(n - 1, n_available) + 1
            distribution%likelihood(n) = likelihood_ratio(row)
         end do
      else
         distribution%likelihood = 1.0_dp
      end if
      distribution%likelihood = distribution%likelihood / sum(distribution%likelihood)

      if (portfolio%n_sectors == 1) then
         allocate(sigma(1, 1))
         sigma(1, 1) = 1.0_dp
      else
         call correlation_matrix(factors, sigma, corr_status, corr_message)
         if (corr_status == 1) then
            status = 6
            message = corr_message
            return
         end if
      end if

      if (model%link_kind == link_cm) then
         allocate(dd(portfolio%n_counterparties))
         dd = normal_quantile(min(1.0_dp - 1.0e-12_dp, &
              max(1.0e-12_dp, portfolio%discretized_pd)))
         do i = 1, portfolio%n_counterparties
            systematic_variance = quadratic_form(portfolio%weight(i, :), sigma)
            if (systematic_variance >= 1.0_dp - 1.0e-12_dp) then
               status = 7
               message = 'CreditMetrics weights imply systematic variance greater than or equal to one'
               return
            end if
         end do
      else if (model%link_kind /= link_crp) then
         status = 8
         message = 'unknown link function'
         return
      end if

      max_stored = max(0, model%max_stored_scenarios)
      allocate(cp_store(portfolio%n_counterparties, max_stored))
      allocate(scenario_store(max_stored))
      cp_store = 0.0_dp
      scenario_store = 0
      n_stored = 0
      overflow = .false.

      call seed_random_number(model%seed)
      do n = 1, n_sim
         row = modulo(n - 1, n_available) + 1
         cp_loss = 0.0_dp
         total_loss = 0.0_dp

         do i = 1, portfolio%n_counterparties
            if (model%link_kind == link_crp) then
               conditional_pd = portfolio%discretized_pd(i) * &
                  (portfolio%idiosyncratic(i) + dot_product(portfolio%weight(i, :), &
                                                            factors(row, :)))
            else
               systematic_variance = quadratic_form(portfolio%weight(i, :), sigma)
               denom = sqrt(max(tiny(1.0_dp), 1.0_dp - systematic_variance))
               conditional_pd = normal_cdf((dd(i) - dot_product(portfolio%weight(i, :), &
                                                                  factors(row, :))) / denom)
            end if

            if (portfolio%default_kind(i) == default_bernoulli) then
               conditional_pd = min(1.0_dp, max(1.0e-7_dp, conditional_pd))
               call random_number(u)
               if (u <= conditional_pd) cp_loss(i) = portfolio%discretized_loss(i)
            else
               conditional_pd = max(1.0e-7_dp, conditional_pd)
               n_defaults = random_poisson(conditional_pd)
               cp_loss(i) = portfolio%discretized_loss(i) * real(n_defaults, dp)
            end if
            total_loss = total_loss + cp_loss(i)
         end do

         distribution%simulated_losses(n) = total_loss
         if (model%loss_threshold < huge(1.0_dp) .and. total_loss >= model%loss_threshold) then
            if (n_stored < max_stored) then
               n_stored = n_stored + 1
               scenario_store(n_stored) = n
               cp_store(:, n_stored) = cp_loss
            else
               overflow = .true.
            end if
         end if
      end do

      if (abs(model%loss_unit) >= 1.0e-8_dp) then
         distribution%simulated_losses = anint(distribution%simulated_losses / &
                                                model%loss_unit) * model%loss_unit
         if (n_stored > 0) then
            cp_store(:, 1:n_stored) = anint(cp_store(:, 1:n_stored) / model%loss_unit) * &
                                       model%loss_unit
         end if
      end if

      allocate(sorted_loss(n_sim), sorted_weight(n_sim))
      sorted_loss = distribution%simulated_losses
      sorted_weight = distribution%likelihood
      call quicksort_pairs(sorted_loss, sorted_weight, 1, n_sim)
      call aggregate_distribution(sorted_loss, sorted_weight, distribution%loss, &
                                  distribution%pdf, distribution%cdf)

      distribution%reached_alpha = distribution%cdf(size(distribution%cdf))
      distribution%expected_loss = sum(distribution%loss * distribution%pdf)
      distribution%standard_deviation = sqrt(sum(distribution%pdf * &
         (distribution%loss - distribution%expected_loss) ** 2))

      if (n_stored > 0 .and. .not. overflow) then
         allocate(distribution%stored_scenario(n_stored))
         allocate(distribution%stored_counterparty_losses(portfolio%n_counterparties, n_stored))
         distribution%stored_scenario = scenario_store(1:n_stored)
         distribution%stored_counterparty_losses = cp_store(:, 1:n_stored)
         distribution%contributions_available = .true.
      else
         allocate(distribution%stored_scenario(0))
         allocate(distribution%stored_counterparty_losses(portfolio%n_counterparties, 0))
         distribution%contributions_available = .false.
      end if

      if (overflow) then
         status = 9
         message = 'stored loss scenarios exceeded max_stored_scenarios; simulation contributions are unavailable'
      else if (model%loss_threshold < huge(1.0_dp) .and. n_stored == 0) then
         status = 10
         message = 'no simulated loss reached loss_threshold; simulation contributions are unavailable'
      end if
   end subroutine simulate_loss_distribution

   subroutine aggregate_distribution(sorted_loss, sorted_weight, loss, pdf, cdf)
      real(dp), intent(in) :: sorted_loss(:)
      real(dp), intent(in) :: sorted_weight(:)
      real(dp), allocatable, intent(out) :: loss(:)
      real(dp), allocatable, intent(out) :: pdf(:)
      real(dp), allocatable, intent(out) :: cdf(:)
      real(dp), allocatable :: loss_work(:), pdf_work(:)
      integer :: i, n_unique, n

      n = size(sorted_loss)
      allocate(loss_work(n), pdf_work(n))
      loss_work = 0.0_dp
      pdf_work = 0.0_dp
      n_unique = 1
      loss_work(1) = sorted_loss(1)
      pdf_work(1) = sorted_weight(1)
      do i = 2, n
         if (abs(sorted_loss(i) - loss_work(n_unique)) < 1.0e-8_dp) then
            pdf_work(n_unique) = pdf_work(n_unique) + sorted_weight(i)
         else
            n_unique = n_unique + 1
            loss_work(n_unique) = sorted_loss(i)
            pdf_work(n_unique) = sorted_weight(i)
         end if
      end do

      allocate(loss(n_unique), pdf(n_unique), cdf(n_unique))
      loss = loss_work(1:n_unique)
      pdf = pdf_work(1:n_unique)
      cdf = 0.0_dp
      cdf(1) = pdf(1)
      do i = 2, n_unique
         cdf(i) = cdf(i - 1) + pdf(i)
      end do
   end subroutine aggregate_distribution

   recursive subroutine quicksort_pairs(x, w, left, right)
      real(dp), intent(inout) :: x(:)
      real(dp), intent(inout) :: w(:)
      integer, intent(in) :: left, right
      integer :: i, j
      real(dp) :: pivot, temp

      if (left >= right) return
      i = left
      j = right
      pivot = x((left + right) / 2)
      do
         do while (x(i) < pivot)
            i = i + 1
         end do
         do while (x(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            temp = x(i)
            x(i) = x(j)
            x(j) = temp
            temp = w(i)
            w(i) = w(j)
            w(j) = temp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (left < j) call quicksort_pairs(x, w, left, j)
      if (i < right) call quicksort_pairs(x, w, i, right)
   end subroutine quicksort_pairs

end module gcpm_simulation
