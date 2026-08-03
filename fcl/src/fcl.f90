! MIT License. Copyright (c) 2024 fcl authors.
module fcl
   use fcl_kinds, only : dp
   use fcl_dates, only : date_type, make_date, add_months, year_frac, days_between
   use fcl_xirr, only : xnpv, xirr
   use fcl_bond, only : fixed_bond_type, bond_value_type, cashflow_type
   use fcl_returns, only : return_series_type
   implicit none
   public
end module fcl
