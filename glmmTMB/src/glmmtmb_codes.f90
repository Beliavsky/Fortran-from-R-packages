! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_codes
   implicit none
   private
   integer, parameter, public :: gaussian_family = 0
   integer, parameter, public :: binomial_family = 100
   integer, parameter, public :: betabinomial_family = 101
   integer, parameter, public :: beta_family = 200
   integer, parameter, public :: ordbeta_family = 201
   integer, parameter, public :: gamma_family = 300
   integer, parameter, public :: poisson_family = 400
   integer, parameter, public :: truncated_poisson_family = 401
   integer, parameter, public :: genpois_family = 402
   integer, parameter, public :: compois_family = 403
   integer, parameter, public :: truncated_genpois_family = 404
   integer, parameter, public :: truncated_compois_family = 405
   integer, parameter, public :: nbinom1_family = 500
   integer, parameter, public :: nbinom2_family = 501
   integer, parameter, public :: nbinom12_family = 502
   integer, parameter, public :: truncated_nbinom1_family = 550
   integer, parameter, public :: truncated_nbinom2_family = 551
   integer, parameter, public :: t_family = 600
   integer, parameter, public :: tweedie_family = 700
   integer, parameter, public :: lognormal_family = 800
   integer, parameter, public :: skewnormal_family = 900
   integer, parameter, public :: bell_family = 1000

   integer, parameter, public :: log_link = 0
   integer, parameter, public :: logit_link = 1
   integer, parameter, public :: probit_link = 2
   integer, parameter, public :: inverse_link = 3
   integer, parameter, public :: cloglog_link = 4
   integer, parameter, public :: identity_link = 5
   integer, parameter, public :: sqrt_link = 6
   integer, parameter, public :: lambertw_link = 7

   integer, parameter, public :: diag_covstruct = 0
   integer, parameter, public :: us_covstruct = 1
   integer, parameter, public :: cs_covstruct = 2
   integer, parameter, public :: ar1_covstruct = 3
   integer, parameter, public :: ou_covstruct = 4
   integer, parameter, public :: exp_covstruct = 5
   integer, parameter, public :: gau_covstruct = 6
   integer, parameter, public :: mat_covstruct = 7
   integer, parameter, public :: toep_covstruct = 8
   integer, parameter, public :: rr_covstruct = 9
   integer, parameter, public :: homdiag_covstruct = 10
   integer, parameter, public :: propto_covstruct = 11
   integer, parameter, public :: hetar1_covstruct = 12
   integer, parameter, public :: homcs_covstruct = 13
   integer, parameter, public :: homtoep_covstruct = 14
   integer, parameter, public :: equalto_covstruct = 15

   integer, parameter, public :: normal_prior = 0
   integer, parameter, public :: t_prior = 1
   integer, parameter, public :: cauchy_prior = 2
   integer, parameter, public :: gamma_prior = 10
   integer, parameter, public :: beta_prior = 20
   integer, parameter, public :: lkj_prior = 30
end module glmmtmb_codes
