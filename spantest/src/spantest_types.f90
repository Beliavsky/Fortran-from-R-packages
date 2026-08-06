! SPDX-License-Identifier: GPL-3.0-only
module spantest_types
  use spantest_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: span_ok = 0
  integer, parameter, public :: span_invalid_input = 1
  integer, parameter, public :: span_singular = 2
  integer, parameter, public :: span_insufficient_df = 3
  integer, parameter, public :: span_numerical_failure = 4

  type, public :: span_result
    real(dp) :: pval = 0.0_dp
    real(dp) :: stat = 0.0_dp
    character(len=40) :: h0 = ''
    integer :: status = span_ok
    character(len=160) :: message = ''
  end type span_result

  type, public :: gl_result
    real(dp) :: pval_lmc = 0.0_dp
    real(dp) :: pval_bmc = 0.0_dp
    real(dp) :: stat = 0.0_dp
    integer :: decision = -1
    character(len=12) :: decision_string = 'Inconclusive'
    character(len=40) :: h0 = ''
    integer :: status = span_ok
    character(len=160) :: message = ''
  end type gl_result

  type, public :: as_result
    real(dp), allocatable :: pvalues(:)
    character(len=40), allocatable :: names(:)
    integer :: status = span_ok
    character(len=160) :: message = ''
  end type as_result

  type, public :: simulation_result
    real(dp), allocatable :: r1(:,:)
    real(dp), allocatable :: r2(:,:)
    integer :: status = span_ok
    character(len=160) :: message = ''
  end type simulation_result

end module spantest_types
