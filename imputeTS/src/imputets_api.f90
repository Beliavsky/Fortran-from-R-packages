module imputets_api
   use imputets_kinds, only : dp
   use imputets_types, only : na_stats_result
   use imputets_basic, only : na_locf, na_ma, na_mean, na_random, na_remove, na_replace
   use imputets_interpolation, only : na_interpolation
   use imputets_kalman, only : na_kalman
   use imputets_seasonal, only : na_seadec, na_seasplit
   use imputets_stats, only : stats_na
   implicit none
   private

   public :: dp
   public :: na_stats_result
   public :: na_interpolation
   public :: na_kalman
   public :: na_locf
   public :: na_ma
   public :: na_mean
   public :: na_random
   public :: na_remove
   public :: na_replace
   public :: na_seadec
   public :: na_seasplit
   public :: stats_na
end module imputets_api
