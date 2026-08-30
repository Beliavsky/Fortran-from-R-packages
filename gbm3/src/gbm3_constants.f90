! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_constants
   implicit none
   private

   integer, parameter, public :: GBM_GAUSSIAN = 1
   integer, parameter, public :: GBM_BERNOULLI = 2
   integer, parameter, public :: GBM_POISSON = 3
   integer, parameter, public :: GBM_GAMMA = 4
   integer, parameter, public :: GBM_LAPLACE = 5
   integer, parameter, public :: GBM_TDIST = 6
   integer, parameter, public :: GBM_QUANTILE = 7
   integer, parameter, public :: GBM_ADABOOST = 8
   integer, parameter, public :: GBM_HUBERIZED = 9
   integer, parameter, public :: GBM_TWEEDIE = 10
   integer, parameter, public :: GBM_COXPH = 11
   integer, parameter, public :: GBM_PAIRWISE = 12

   integer, parameter, public :: GBM_TIES_BRESLOW = 0
   integer, parameter, public :: GBM_TIES_EFRON = 1

   integer, parameter, public :: GBM_METRIC_NDCG = 1
   integer, parameter, public :: GBM_METRIC_CONC = 2
   integer, parameter, public :: GBM_METRIC_MAP = 3
   integer, parameter, public :: GBM_METRIC_MRR = 4
end module gbm3_constants
