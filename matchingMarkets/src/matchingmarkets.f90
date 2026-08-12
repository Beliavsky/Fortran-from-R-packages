module matchingmarkets
   use matchingmarkets_kinds, only : dp, i8
   use matchingmarkets_types
   use matchingmarkets_rng, only : rng_t
   use matchingmarkets_mechanisms
   use matchingmarkets_stats
   use matchingmarkets_plp
   use matchingmarkets_sim
   use matchingmarkets_helpers
   use matchingr, only : marriage_market, marriage_market_preferences, &
      college_admissions, college_admissions_preferences, stable_roommates, &
      stable_roommates_preferences, top_trading_cycles, top_trading_cycles_preferences
   implicit none
   public

   interface hri
      procedure :: hri_all
   end interface
   interface hri2
      procedure :: hri2_couples_exact
   end interface
   interface hri3
      procedure :: hri3_eadam
   end interface
   interface sri
      procedure :: sri_all
   end interface
   interface ttc
      procedure :: ttc_tenants
   end interface
   interface ttc2
      procedure :: ttc_school
   end interface
   interface ttcc
      procedure :: ttcc_kidney
   end interface
   interface stabchk
      procedure :: stability_check
   end interface
end module matchingmarkets
