! SPDX-License-Identifier: GPL-3.0-only
! Derived from computational code in R package pan 2.0.
! Upstream authorship/maintenance: Joseph L. Schafer and Jing Hua Zhao.
! Public facade for the modern Fortran translation of R package pan 2.0.
module pan
   use pan_ecme, only : ecme_fit
   use pan_kinds, only : dp
   use pan_sampler, only : pan_bd_mcmc, pan_mcmc
   use pan_types, only : PAN_ERR_ARGUMENT, PAN_ERR_DIMENSION, PAN_ERR_LINALG, PAN_ERR_NUMERIC, PAN_OK, &
      ecme_result, pan_bd_prior, pan_bd_result, pan_bd_state, pan_prior, pan_result, pan_state
   implicit none
   private

   public :: dp
   public :: PAN_OK
   public :: PAN_ERR_DIMENSION
   public :: PAN_ERR_ARGUMENT
   public :: PAN_ERR_LINALG
   public :: PAN_ERR_NUMERIC
   public :: pan_prior
   public :: pan_bd_prior
   public :: pan_state
   public :: pan_bd_state
   public :: pan_result
   public :: pan_bd_result
   public :: ecme_result
   public :: pan_mcmc
   public :: pan_bd_mcmc
   public :: ecme_fit

end module pan
