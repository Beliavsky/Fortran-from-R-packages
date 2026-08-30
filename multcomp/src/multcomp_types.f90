! SPDX-License-Identifier: GPL-2.0-only
module multcomp_types
  use multcomp_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: alt_two_sided = 1
  integer, parameter, public :: alt_less = 2
  integer, parameter, public :: alt_greater = 3

  type, public :: parm_type
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: vcov(:, :)
    real(dp) :: df = 0.0_dp
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type parm_type

  type, public :: glht_type
    real(dp), allocatable :: linfct(:, :)
    real(dp), allocatable :: rhs(:)
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: standard_error(:)
    real(dp), allocatable :: statistic(:)
    real(dp), allocatable :: correlation(:, :)
    real(dp) :: df = 0.0_dp
    integer :: alternative = alt_two_sided
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type glht_type

  type, public :: mtest_type
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: standard_error(:)
    real(dp), allocatable :: statistic(:)
    real(dp), allocatable :: pvalue(:)
    real(dp) :: error = 0.0_dp
    character(len=32) :: method = ''
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type mtest_type

  type, public :: confidence_interval_type
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    real(dp) :: level = 0.95_dp
    real(dp) :: critical = 0.0_dp
    real(dp) :: error = 0.0_dp
    logical :: adjusted = .true.
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type confidence_interval_type

  type, public :: global_test_type
    real(dp) :: ssh = 0.0_dp
    real(dp) :: statistic = 0.0_dp
    real(dp) :: pvalue = 1.0_dp
    integer :: rank = 0
    real(dp) :: denominator_df = 0.0_dp
    logical :: f_test = .false.
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type global_test_type

  type, public :: contrast_matrix_type
    real(dp), allocatable :: value(:, :)
    character(len=32) :: contrast_type = ''
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type contrast_matrix_type

  type, public :: integer_set
    integer, allocatable :: value(:)
  end type integer_set

  type, public :: set_collection
    type(integer_set), allocatable :: set(:)
  end type set_collection

  type, public :: cld_type
    logical, allocatable :: letter_matrix(:, :)
    character(len=256), allocatable :: letters(:)
    character(len=256), allocatable :: monospaced_letters(:)
    character(len=32), allocatable :: column_labels(:)
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type cld_type

end module multcomp_types
