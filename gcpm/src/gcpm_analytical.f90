! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_analytical
   use gcpm_kinds, only: dp
   use gcpm_types, only: credit_portfolio, gcpm_model, loss_distribution, &
      model_analytical_crp
   use gcpm_portfolio, only: sector_information
   implicit none
   private

   public :: analytical_creditrisk_plus

contains

   subroutine analytical_creditrisk_plus(portfolio, model, distribution, status, message, max_bands)
      type(credit_portfolio), intent(in) :: portfolio
      type(gcpm_model), intent(in) :: model
      type(loss_distribution), intent(out) :: distribution
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      integer, intent(in), optional :: max_bands

      integer :: i, j, k, r, m, n_support
      real(dp) :: coefficient_sum
      real(dp), allocatable :: sector_mean(:), sector_loss(:)
      real(dp), allocatable :: a_work(:), pdf_work(:), cdf_work(:)
      real(dp), allocatable :: aa(:,:), bb(:,:)

      status = 0
      message = ''
      distribution%model_kind = model_analytical_crp

      if (model%model_kind /= model_analytical_crp) then
         status = 1
         message = 'analytical_creditrisk_plus requires model_analytical_crp'
         return
      end if
      if (.not. allocated(model%sector_variance)) then
         status = 2
         message = 'sector_variance is not allocated'
         return
      end if
      if (.not. allocated(portfolio%discretized_pd)) then
         status = 3
         message = 'calculate_portfolio_statistics must be called first'
         return
      end if

      m = 100000
      if (present(max_bands)) m = max(2, max_bands)
      allocate(a_work(m), pdf_work(m), cdf_work(m))
      allocate(aa(portfolio%n_sectors, m), bb(portfolio%n_sectors, m))
      a_work = 0.0_dp
      pdf_work = 0.0_dp
      cdf_work = 0.0_dp
      aa = 0.0_dp
      bb = 0.0_dp

      call sector_information(portfolio, sector_mean, sector_loss)

      if (sum(portfolio%weight) <= tiny(1.0_dp)) then
         a_work(1) = -sum(portfolio%idiosyncratic * portfolio%discretized_pd)
         pdf_work(1) = exp(a_work(1))
         cdf_work(1) = pdf_work(1)
         j = 0
         do while (j < m - 1 .and. cdf_work(j + 1) < model%alpha_max)
            j = j + 1
            a_work(j + 1) = 0.0_dp
            do i = 1, portfolio%n_counterparties
               if (portfolio%loss_multiple(i) == j) then
                  a_work(j + 1) = a_work(j + 1) + &
                     portfolio%idiosyncratic(i) * portfolio%discretized_pd(i)
               end if
            end do
            coefficient_sum = 0.0_dp
            do r = 1, j
               coefficient_sum = coefficient_sum + real(r, dp) / real(j, dp) * &
                  a_work(r + 1) * pdf_work(j - r + 1)
            end do
            pdf_work(j + 1) = coefficient_sum
            if (pdf_work(j + 1) < 0.0_dp .and. pdf_work(j + 1) > -1.0e-14_dp) &
               pdf_work(j + 1) = 0.0_dp
            cdf_work(j + 1) = cdf_work(j) + pdf_work(j + 1)
         end do
      else
         do k = 1, portfolio%n_sectors
            aa(k, 1) = 1.0_dp + model%sector_variance(k) * sector_mean(k)
            bb(k, 1) = -log(aa(k, 1))
         end do
         a_work(1) = -sum(portfolio%idiosyncratic * portfolio%discretized_pd) + &
            sum(bb(:, 1) / model%sector_variance)
         pdf_work(1) = exp(a_work(1))
         cdf_work(1) = pdf_work(1)

         j = 0
         do while (j < m - 1 .and. cdf_work(j + 1) < model%alpha_max)
            j = j + 1
            do k = 1, portfolio%n_sectors
               aa(k, j + 1) = 0.0_dp
               do i = 1, portfolio%n_counterparties
                  if (portfolio%loss_multiple(i) == j) then
                     aa(k, j + 1) = aa(k, j + 1) + model%sector_variance(k) * &
                        portfolio%weight(i, k) * portfolio%discretized_pd(i)
                  end if
               end do

               if (j == 1) then
                  bb(k, 2) = aa(k, 2) / aa(k, 1)
               else
                  coefficient_sum = 0.0_dp
                  do r = 1, j - 1
                     coefficient_sum = coefficient_sum + real(r, dp) * &
                        bb(k, r + 1) * aa(k, j - r + 1)
                  end do
                  bb(k, j + 1) = (aa(k, j + 1) + coefficient_sum / real(j, dp)) / &
                     aa(k, 1)
               end if
            end do

            a_work(j + 1) = 0.0_dp
            do i = 1, portfolio%n_counterparties
               if (portfolio%loss_multiple(i) == j) then
                  a_work(j + 1) = a_work(j + 1) + &
                     portfolio%idiosyncratic(i) * portfolio%discretized_pd(i)
               end if
            end do
            a_work(j + 1) = a_work(j + 1) + &
               sum(bb(:, j + 1) / model%sector_variance)

            coefficient_sum = 0.0_dp
            do r = 1, j
               coefficient_sum = coefficient_sum + real(r, dp) / real(j, dp) * &
                  a_work(r + 1) * pdf_work(j - r + 1)
            end do
            pdf_work(j + 1) = coefficient_sum
            if (pdf_work(j + 1) < 0.0_dp .and. pdf_work(j + 1) > -1.0e-14_dp) &
               pdf_work(j + 1) = 0.0_dp
            cdf_work(j + 1) = cdf_work(j) + pdf_work(j + 1)
         end do
      end if

      n_support = j + 1
      allocate(distribution%loss(n_support), distribution%pdf(n_support), &
               distribution%cdf(n_support), distribution%recursion_a(n_support), &
               distribution%recursion_b(portfolio%n_sectors, n_support))
      distribution%loss = [(real(i - 1, dp) * model%loss_unit, i=1, n_support)]
      distribution%pdf = pdf_work(1:n_support)
      distribution%cdf = cdf_work(1:n_support)
      distribution%recursion_a = a_work(1:n_support)
      distribution%recursion_b = bb(:, 1:n_support)
      distribution%reached_alpha = distribution%cdf(n_support)
      distribution%expected_loss = sum(distribution%loss * distribution%pdf)
      distribution%standard_deviation = sqrt(sum(distribution%pdf * &
         (distribution%loss - distribution%expected_loss) ** 2))

      if (distribution%reached_alpha < model%alpha_max) then
         status = 4
         message = 'maximum exposure-band count reached before alpha_max'
      end if
   end subroutine analytical_creditrisk_plus

end module gcpm_analytical
