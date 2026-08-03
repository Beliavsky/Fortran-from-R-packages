! SPDX-License-Identifier: GPL-2.0-or-later
module fer
   use fer_kinds, only : dp, pi
   use fer_vanilla, only : bachelier_price, bachelier_impvol, black_scholes_price, black_scholes_impvol
   use fer_cev, only : cev_price, cev_mass_zero
   use fer_sabr, only : sabr_hagan_2002, sabr_hagan_price, nsvh1_choi_2019
   use fer_spread, only : switch_margrabe, spread_kirk, spread_bjerksund_2014, spread_bachelier
   implicit none
   public
end module fer
