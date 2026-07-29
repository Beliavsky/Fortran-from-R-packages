! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_portfolio
   use gcpm_kinds, only: dp, name_len, sector_name_len
   use gcpm_types, only: credit_portfolio, gcpm_model, default_bernoulli, &
      default_poisson, link_crp, model_analytical_crp
   implicit none
   private

   public :: allocate_portfolio
   public :: validate_portfolio
   public :: calculate_portfolio_statistics
   public :: sector_information
   public :: standard_deviation_contributions

contains

   subroutine allocate_portfolio(portfolio, n_counterparties, n_sectors)
      type(credit_portfolio), intent(out) :: portfolio
      integer, intent(in) :: n_counterparties
      integer, intent(in) :: n_sectors
      integer :: i

      portfolio%n_counterparties = n_counterparties
      portfolio%n_sectors = n_sectors
      allocate(portfolio%number(n_counterparties))
      allocate(character(len=name_len) :: portfolio%name(n_counterparties))
      allocate(character(len=sector_name_len) :: portfolio%sector_name(n_sectors))
      allocate(portfolio%default_kind(n_counterparties))
      allocate(portfolio%ead(n_counterparties))
      allocate(portfolio%lgd(n_counterparties))
      allocate(portfolio%pd(n_counterparties))
      allocate(portfolio%weight(n_counterparties, n_sectors))
      allocate(portfolio%idiosyncratic(n_counterparties))
      allocate(portfolio%potential_loss(n_counterparties))
      allocate(portfolio%loss_multiple(n_counterparties))
      allocate(portfolio%discretized_loss(n_counterparties))
      allocate(portfolio%discretized_pd(n_counterparties))

      portfolio%number = [(i, i=1, n_counterparties)]
      portfolio%name = ''
      portfolio%sector_name = ''
      portfolio%default_kind = default_poisson
      portfolio%ead = 0.0_dp
      portfolio%lgd = 0.0_dp
      portfolio%pd = 0.0_dp
      portfolio%weight = 0.0_dp
      portfolio%idiosyncratic = 0.0_dp
      portfolio%potential_loss = 0.0_dp
      portfolio%loss_multiple = 0
      portfolio%discretized_loss = 0.0_dp
      portfolio%discretized_pd = 0.0_dp
   end subroutine allocate_portfolio

   subroutine validate_portfolio(portfolio, model, status, message)
      type(credit_portfolio), intent(in) :: portfolio
      type(gcpm_model), intent(in) :: model
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), allocatable :: row_weight(:)

      status = 0
      message = ''

      if (portfolio%n_counterparties <= 0) then
         status = 1
         message = 'the portfolio must contain at least one counterparty'
         return
      end if
      if (portfolio%n_sectors <= 0) then
         status = 2
         message = 'the portfolio must contain at least one sector'
         return
      end if
      if (model%loss_unit <= 0.0_dp) then
         status = 3
         message = 'loss_unit must be positive'
         return
      end if
      if (any(portfolio%ead < 0.0_dp)) then
         status = 4
         message = 'EAD values cannot be negative'
         return
      end if
      if (any(portfolio%lgd < 0.0_dp) .or. any(portfolio%lgd > 1.0_dp)) then
         status = 5
         message = 'LGD values must lie in [0, 1]'
         return
      end if
      if (any(portfolio%pd <= 0.0_dp)) then
         status = 6
         message = 'PD values must be positive'
         return
      end if
      if (any(portfolio%default_kind /= default_bernoulli .and. &
              portfolio%default_kind /= default_poisson)) then
         status = 7
         message = 'default_kind must be default_bernoulli or default_poisson'
         return
      end if
      if (any(portfolio%weight < 0.0_dp)) then
         status = 8
         message = 'sector weights cannot be negative'
         return
      end if

      if (model%link_kind == link_crp .or. model%model_kind == model_analytical_crp) then
         row_weight = sum(portfolio%weight, dim=2)
         if (any(row_weight > 1.0_dp + 1.0e-12_dp)) then
            status = 9
            message = 'CRP sector weights must sum to at most one for each counterparty'
            return
         end if
      end if

      if (model%model_kind == model_analytical_crp) then
         if (.not. allocated(model%sector_variance)) then
            status = 10
            message = 'analytical CreditRisk+ requires sector_variance'
            return
         end if
         if (size(model%sector_variance) /= portfolio%n_sectors) then
            status = 11
            message = 'sector_variance has the wrong length'
            return
         end if
         if (any(model%sector_variance <= 0.0_dp)) then
            status = 12
            message = 'sector variances must be positive'
            return
         end if
      end if
   end subroutine validate_portfolio

   subroutine calculate_portfolio_statistics(portfolio, model, status, message)
      type(credit_portfolio), intent(inout) :: portfolio
      type(gcpm_model), intent(in) :: model
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), allocatable :: sector_mean(:), sector_loss(:)

      call validate_portfolio(portfolio, model, status, message)
      if (status /= 0) return

      portfolio%potential_loss = portfolio%ead * portfolio%lgd
      if (any(portfolio%potential_loss <= 0.0_dp)) then
         status = 13
         message = 'every retained counterparty must have positive EAD times LGD'
         return
      end if

      portfolio%loss_multiple = max(1, floor(portfolio%potential_loss / &
         model%loss_unit + 0.5_dp))
      portfolio%discretized_loss = real(portfolio%loss_multiple, dp) * model%loss_unit
      portfolio%discretized_pd = portfolio%pd * portfolio%potential_loss / &
         portfolio%discretized_loss
      portfolio%analytical_expected_loss = sum(portfolio%potential_loss * portfolio%pd)

      if (model%link_kind == link_crp .or. model%model_kind == model_analytical_crp) then
         portfolio%idiosyncratic = max(0.0_dp, 1.0_dp - sum(portfolio%weight, dim=2))
      else
         portfolio%idiosyncratic = 0.0_dp
      end if

      if (model%model_kind == model_analytical_crp) then
         call sector_information(portfolio, sector_mean, sector_loss)
         portfolio%analytical_sd_systematic = sqrt(sum(model%sector_variance * sector_loss ** 2))
         portfolio%analytical_sd_diversifiable = sqrt(sum( &
            portfolio%discretized_loss ** 2 * portfolio%discretized_pd))
         portfolio%analytical_sd = sqrt(portfolio%analytical_sd_systematic ** 2 + &
            portfolio%analytical_sd_diversifiable ** 2)
      else
         portfolio%analytical_sd_systematic = 0.0_dp
         portfolio%analytical_sd_diversifiable = 0.0_dp
         portfolio%analytical_sd = 0.0_dp
      end if
   end subroutine calculate_portfolio_statistics

   subroutine sector_information(portfolio, sector_mean, sector_loss)
      type(credit_portfolio), intent(in) :: portfolio
      real(dp), allocatable, intent(out) :: sector_mean(:)
      real(dp), allocatable, intent(out) :: sector_loss(:)
      integer :: k

      allocate(sector_mean(portfolio%n_sectors), sector_loss(portfolio%n_sectors))
      do k = 1, portfolio%n_sectors
         sector_mean(k) = sum(portfolio%weight(:, k) * portfolio%discretized_pd)
         sector_loss(k) = sum(portfolio%weight(:, k) * portfolio%discretized_pd * &
                              portfolio%discretized_loss)
      end do
   end subroutine sector_information

   subroutine standard_deviation_contributions(portfolio, model, contribution, status, message)
      type(credit_portfolio), intent(in) :: portfolio
      type(gcpm_model), intent(in) :: model
      real(dp), allocatable, intent(out) :: contribution(:)
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), allocatable :: sector_mean(:), sector_loss(:)
      real(dp) :: systematic_term
      integer :: i

      status = 0
      message = ''
      allocate(contribution(portfolio%n_counterparties))
      contribution = 0.0_dp
      if (model%model_kind /= model_analytical_crp) then
         status = 1
         message = 'standard-deviation contributions require analytical CreditRisk+'
         return
      end if
      if (portfolio%analytical_sd <= tiny(1.0_dp)) then
         status = 2
         message = 'analytical portfolio standard deviation is zero'
         return
      end if

      call sector_information(portfolio, sector_mean, sector_loss)
      do i = 1, portfolio%n_counterparties
         systematic_term = sum(model%sector_variance * portfolio%weight(i, :) * sector_loss)
         contribution(i) = portfolio%discretized_loss(i) * portfolio%discretized_pd(i) / &
            portfolio%analytical_sd * (portfolio%discretized_loss(i) + systematic_term)
      end do
   end subroutine standard_deviation_contributions

end module gcpm_portfolio
