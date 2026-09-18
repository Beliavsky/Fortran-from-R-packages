module randompack
   use randompack_kinds, only : dp
   use randompack_rng_mod, only : randompack_rng, randompack_rng_type, randompack_snapshot, randompack_engines
   implicit none
   private
   public :: dp, randompack_rng, randompack_rng_type, randompack_snapshot, randompack_engines
end module randompack
