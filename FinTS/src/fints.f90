! SPDX-License-Identifier: GPL-2.0-or-later
module fints
   use fints_kinds, only : dp, pi_dp
   use fints_status
   use fints_types
   use fints_finance, only : compoundInterest => compound_interest, &
      simple2logReturns => simple_to_log_returns
   use fints_dates, only : as_yearmon2
   use fints_summary_mod, only : FinTS_stats => fints_summary
   use fints_time_series, only : acf, cross_acf, ArchTest => arch_test, &
      AutocorTest => autocor_test, pacf_from_acf
   use fints_apca, only : apca
   use fints_arma, only : arma_true_acf, findConjugates => find_conjugates, &
      polynomial_roots
   use fints_arima, only : ARIMA => arima_fit
   implicit none
   public

contains

   subroutine plotArmaTrueacf(ar, ma, lag_max, result, pacf, complex_eps)
      real(dp), intent(in) :: ar(:), ma(:)
      integer, intent(in) :: lag_max
      type(arma_acf_result), intent(out) :: result
      logical, intent(in), optional :: pacf
      real(dp), intent(in), optional :: complex_eps

      call arma_true_acf(ar, ma, lag_max, result, partial=pacf, &
         complex_tolerance=complex_eps)
   end subroutine plotArmaTrueacf

end module fints
