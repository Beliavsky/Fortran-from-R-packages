! SPDX-License-Identifier: GPL-2.0-or-later
module fincal_ratios
   use fincal_kinds, only : dp
   implicit none
   private

   public :: sf_ratio, sharpe_ratio, coefficient_variation
   public :: cash_ratio, current_ratio, quick_ratio
   public :: debt_ratio, long_term_debt_to_equity, total_debt_to_equity
   public :: financial_leverage, gross_profit_margin, net_profit_margin
   public :: sfratio, gpm, npm, lt_d2e, total_d2e

   interface sfratio
      module procedure sf_ratio
   end interface sfratio
   interface gpm
      module procedure gross_profit_margin
   end interface gpm
   interface npm
      module procedure net_profit_margin
   end interface npm
   interface lt_d2e
      module procedure long_term_debt_to_equity
   end interface lt_d2e
   interface total_d2e
      module procedure total_debt_to_equity
   end interface total_d2e
contains
   elemental pure function sf_ratio(portfolio_return, threshold_return, standard_deviation) result(value)
      real(dp), intent(in) :: portfolio_return, threshold_return, standard_deviation
      real(dp) :: value
      value = (portfolio_return - threshold_return) / standard_deviation
   end function sf_ratio

   elemental pure function sharpe_ratio(portfolio_return, risk_free_return, standard_deviation) result(value)
      real(dp), intent(in) :: portfolio_return, risk_free_return, standard_deviation
      real(dp) :: value
      value = (portfolio_return - risk_free_return) / standard_deviation
   end function sharpe_ratio

   elemental pure function coefficient_variation(standard_deviation, average) result(value)
      real(dp), intent(in) :: standard_deviation, average
      real(dp) :: value
      value = standard_deviation / average
   end function coefficient_variation

   elemental pure function cash_ratio(cash, marketable_securities, current_liabilities) result(value)
      real(dp), intent(in) :: cash, marketable_securities, current_liabilities
      real(dp) :: value
      value = (cash + marketable_securities) / current_liabilities
   end function cash_ratio

   elemental pure function current_ratio(current_assets, current_liabilities) result(value)
      real(dp), intent(in) :: current_assets, current_liabilities
      real(dp) :: value
      value = current_assets / current_liabilities
   end function current_ratio

   elemental pure function quick_ratio(cash, marketable_securities, receivables, current_liabilities) result(value)
      real(dp), intent(in) :: cash, marketable_securities, receivables, current_liabilities
      real(dp) :: value
      value = (cash + marketable_securities + receivables) / current_liabilities
   end function quick_ratio

   elemental pure function debt_ratio(total_debt, total_assets) result(value)
      real(dp), intent(in) :: total_debt, total_assets
      real(dp) :: value
      value = total_debt / total_assets
   end function debt_ratio

   elemental pure function long_term_debt_to_equity(long_term_debt, total_equity) result(value)
      real(dp), intent(in) :: long_term_debt, total_equity
      real(dp) :: value
      value = long_term_debt / total_equity
   end function long_term_debt_to_equity

   elemental pure function total_debt_to_equity(total_debt, total_equity) result(value)
      real(dp), intent(in) :: total_debt, total_equity
      real(dp) :: value
      value = total_debt / total_equity
   end function total_debt_to_equity

   elemental pure function financial_leverage(total_equity, total_assets) result(value)
      real(dp), intent(in) :: total_equity, total_assets
      real(dp) :: value
      value = total_assets / total_equity
   end function financial_leverage

   elemental pure function gross_profit_margin(gross_profit, revenue) result(value)
      real(dp), intent(in) :: gross_profit, revenue
      real(dp) :: value
      value = gross_profit / revenue
   end function gross_profit_margin

   elemental pure function net_profit_margin(net_income, revenue) result(value)
      real(dp), intent(in) :: net_income, revenue
      real(dp) :: value
      value = net_income / revenue
   end function net_profit_margin
end module fincal_ratios
