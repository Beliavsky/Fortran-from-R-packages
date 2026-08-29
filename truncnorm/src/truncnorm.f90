module truncnorm
use r_compat, only: dp, set_seed_int
use truncnorm_core, only: dtruncnorm_scalar, ptruncnorm_scalar, qtruncnorm_scalar, &
   etruncnorm_scalar, vtruncnorm_scalar, rtruncnorm_scalar, &
   dtruncnorm_recycle, ptruncnorm_recycle, qtruncnorm_recycle, &
   dtruncnorm_vec, ptruncnorm_vec, qtruncnorm_vec, rtruncnorm_n, &
   etruncnorm_recycle, vtruncnorm_recycle, rtruncnorm_recycle
implicit none
private
public :: dp, set_seed_int
public :: dtruncnorm, ptruncnorm, qtruncnorm, etruncnorm, vtruncnorm, rtruncnorm
public :: dtruncnorm_recycle, ptruncnorm_recycle, qtruncnorm_recycle
public :: etruncnorm_recycle, vtruncnorm_recycle, rtruncnorm_recycle

interface dtruncnorm
   module procedure dtruncnorm_scalar
   module procedure dtruncnorm_vec
   module procedure dtruncnorm_recycle
end interface
interface ptruncnorm
   module procedure ptruncnorm_scalar
   module procedure ptruncnorm_vec
   module procedure ptruncnorm_recycle
end interface
interface qtruncnorm
   module procedure qtruncnorm_scalar
   module procedure qtruncnorm_vec
   module procedure qtruncnorm_recycle
end interface
interface etruncnorm
   module procedure etruncnorm_scalar
end interface
interface vtruncnorm
   module procedure vtruncnorm_scalar
end interface
interface rtruncnorm
   module procedure rtruncnorm_scalar
   module procedure rtruncnorm_n
   module procedure rtruncnorm_recycle
end interface
end module truncnorm
