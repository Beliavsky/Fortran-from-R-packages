module pbs_mod
   use pbs_api, only : dp, pbs_basis, pbs, predict_pbs, PBS_OK, PBS_INVALID_ARGUMENT, PBS_OUTSIDE_BOUNDARY
   implicit none
   private

   public :: dp
   public :: pbs_basis
   public :: pbs
   public :: predict_pbs
   public :: PBS_OK
   public :: PBS_INVALID_ARGUMENT
   public :: PBS_OUTSIDE_BOUNDARY
end module pbs_mod
