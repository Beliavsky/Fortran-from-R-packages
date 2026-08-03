! SPDX-License-Identifier: GPL-3.0-only
module indgenerrors
  use indgen_kinds, only : dp
  use indgen_types
  use indgen_core, only : cvm_2series_core, cvm_3series_core, &
    crosscor_2series_core, crosscor_3series_core, &
    crossdep_2series_core, crossdep_3series_core
  implicit none
  private

  public :: dp
  public :: indgen_success, indgen_invalid_argument, indgen_numerical_error
  public :: lag_test_result, cvm_test_result, four_lag_test_result
  public :: cvm_three_result, dependence_two_result, dependence_three_result
  public :: cvm_2series, cvm_3series, crosscor_2series, crosscor_3series
  public :: crossdep_2series, crossdep_3series

contains

  function cvm_2series(x,y,lag) result(out)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lag
    type(cvm_test_result) :: out
    out = cvm_2series_core(x,y,lag)
  end function cvm_2series

  function cvm_3series(x,y,z,lag2,lag3) result(out)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    type(cvm_three_result) :: out
    out = cvm_3series_core(x,y,z,lag2,lag3)
  end function cvm_3series

  function crosscor_2series(x,y,lag) result(out)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lag
    type(lag_test_result) :: out
    out = crosscor_2series_core(x,y,lag)
  end function crosscor_2series

  function crosscor_3series(x,y,z,lag2,lag3) result(out)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    type(four_lag_test_result) :: out
    out = crosscor_3series_core(x,y,z,lag2,lag3)
  end function crosscor_3series

  function crossdep_2series(x,y,lag) result(out)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lag
    type(dependence_two_result) :: out
    out = crossdep_2series_core(x,y,lag)
  end function crossdep_2series

  function crossdep_3series(x,y,z,lag2,lag3) result(out)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    type(dependence_three_result) :: out
    out = crossdep_3series_core(x,y,z,lag2,lag3)
  end function crossdep_3series

end module indgenerrors
