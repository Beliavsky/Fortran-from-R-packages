! SPDX-License-Identifier: GPL-2.0-or-later
module jfe_kinds
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: i8 = int64

   integer, parameter, public :: jfe_success = 0
   integer, parameter, public :: jfe_invalid_argument = 1
   integer, parameter, public :: jfe_insufficient_data = 2
   integer, parameter, public :: jfe_zero_denominator = 3
   integer, parameter, public :: jfe_nonfinite_result = 4

   integer, parameter, public :: risk_stddev = 1
   integer, parameter, public :: risk_var = 2
   integer, parameter, public :: risk_es = 3

   integer, parameter, public :: downside_full = 1
   integer, parameter, public :: downside_subset = 2

   integer, parameter, public :: appraisal_standard = 1
   integer, parameter, public :: appraisal_modified = 2
   integer, parameter, public :: appraisal_alternative = 3

   integer, parameter, public :: volatility_ratio = 1
   integer, parameter, public :: variability_ratio = 2

   integer, parameter, public :: skew_moment = 1
   integer, parameter, public :: skew_fisher = 2
   integer, parameter, public :: skew_sample = 3

   integer, parameter, public :: kurt_excess = 1
   integer, parameter, public :: kurt_moment = 2
   integer, parameter, public :: kurt_fisher = 3
   integer, parameter, public :: kurt_sample = 4
   integer, parameter, public :: kurt_sample_excess = 5
end module jfe_kinds
