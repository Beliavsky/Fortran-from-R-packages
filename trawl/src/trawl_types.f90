module trawl_types
  use trawl_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: trawl_ok = 0
  integer, parameter, public :: trawl_invalid_argument = 1
  integer, parameter, public :: trawl_optimization_failed = 2
  integer, parameter, public :: trawl_numerical_failure = 3

  type, public :: trawl_spec
    character(len=8) :: kind = 'Exp'
    real(dp) :: lambda1 = 0.0_dp
    real(dp) :: lambda2 = 0.0_dp
    real(dp) :: w = 0.0_dp
    real(dp) :: delta = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: h = 0.0_dp
  end type trawl_spec

  type, public :: trawl_fit_result
    real(dp) :: w = 0.0_dp
    real(dp) :: lambda1 = 0.0_dp
    real(dp) :: lambda2 = 0.0_dp
    real(dp) :: delta = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: h = 0.0_dp
    real(dp) :: lm = 0.0_dp
    real(dp) :: objective = 0.0_dp
    integer :: status = trawl_ok
  end type trawl_fit_result

  type, public :: poisson_fit_result
    real(dp) :: v = 0.0_dp
    integer :: status = trawl_ok
  end type poisson_fit_result

  type, public :: nb_fit_result
    real(dp) :: m = 0.0_dp
    real(dp) :: theta = 0.0_dp
    real(dp) :: a = 0.0_dp
    integer :: status = trawl_ok
  end type nb_fit_result

end module trawl_types
