! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_types
  use copula_kinds, only : dp, i8
  implicit none
  private
  integer, parameter, public :: family_independence = 0
  integer, parameter, public :: family_gaussian = 1
  integer, parameter, public :: family_student = 2
  integer, parameter, public :: family_clayton = 3
  integer, parameter, public :: family_gumbel = 4
  integer, parameter, public :: family_frank = 5
  integer, parameter, public :: family_amh = 6
  integer, parameter, public :: family_joe = 7
  integer, parameter, public :: family_fgm = 8
  integer, parameter, public :: family_plackett = 9
  integer, parameter, public :: family_marshall_olkin = 10
  integer, parameter, public :: family_lower_fh = 11
  integer, parameter, public :: family_upper_fh = 12
  integer, parameter, public :: family_galambos = 13
  integer, parameter, public :: family_husler_reiss = 14
  integer, parameter, public :: family_tawn = 15

  integer, parameter, public :: rotation_none = 0
  integer, parameter, public :: rotation_90 = 90
  integer, parameter, public :: rotation_180 = 180
  integer, parameter, public :: rotation_270 = 270

  type, public :: copula_model
    integer :: family = family_independence
    integer :: dimension = 2
    integer :: rotation = rotation_none
    real(dp) :: theta = 0.0_dp
    real(dp) :: df = 4.0_dp
    real(dp) :: alpha1 = 1.0_dp
    real(dp) :: alpha2 = 1.0_dp
    real(dp), allocatable :: correlation(:,:)
  contains
    procedure :: valid => copula_model_valid
  end type copula_model

  type, public :: probability_control
    integer :: max_evaluations = 50000
    integer :: batches = 10
    real(dp) :: absolute_tolerance = 1.0e-6_dp
    integer(i8) :: seed = 1234567_i8
  end type probability_control

  type, public :: probability_result
    real(dp) :: value = 0.0_dp
    real(dp) :: error = 0.0_dp
    integer :: evaluations = 0
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type probability_result

  type, public :: fit_result
    type(copula_model) :: model
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: standard_error = huge(1.0_dp)
    integer :: iterations = 0
    logical :: converged = .false.
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type fit_result

  type, public :: test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    integer :: replicates = 0
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type test_result
contains
  logical function copula_model_valid(self) result(ok)
    class(copula_model), intent(in) :: self
    integer :: d
    d = self%dimension
    ok = d >= 2
    if (.not. ok) return
    select case (self%family)
    case (family_independence)
      ok = .true.
    case (family_gaussian)
      ok = allocated(self%correlation)
      if (ok) ok = size(self%correlation,1) == d .and. size(self%correlation,2) == d
    case (family_student)
      ok = allocated(self%correlation) .and. self%df > 0.0_dp
      if (ok) ok = size(self%correlation,1) == d .and. size(self%correlation,2) == d
    case (family_clayton)
      ok = self%theta > 0.0_dp
    case (family_gumbel, family_joe)
      ok = self%theta >= 1.0_dp
    case (family_frank)
      ok = abs(self%theta) > sqrt(epsilon(1.0_dp))
      if (d > 2) ok = ok .and. self%theta > 0.0_dp
    case (family_amh, family_fgm)
      ok = d == 2 .and. self%theta >= -1.0_dp .and. self%theta <= 1.0_dp
    case (family_plackett)
      ok = d == 2 .and. self%theta > 0.0_dp
    case (family_marshall_olkin)
      ok = d == 2 .and. self%alpha1 >= 0.0_dp .and. self%alpha1 <= 1.0_dp .and. &
        self%alpha2 >= 0.0_dp .and. self%alpha2 <= 1.0_dp
    case (family_lower_fh, family_upper_fh)
      ok = d == 2
    case (family_galambos, family_husler_reiss)
      ok = d == 2 .and. self%theta > 0.0_dp
    case (family_tawn)
      ok = d == 2 .and. self%theta >= 1.0_dp .and. self%alpha1 >= 0.0_dp .and. &
        self%alpha1 <= 1.0_dp .and. self%alpha2 >= 0.0_dp .and. self%alpha2 <= 1.0_dp
    case default
      ok = .false.
    end select
    if (self%rotation /= rotation_none .and. self%rotation /= rotation_90 .and. &
        self%rotation /= rotation_180 .and. self%rotation /= rotation_270) ok = .false.
    if (self%rotation /= rotation_none .and. d /= 2) ok = .false.
  end function copula_model_valid
end module copula_types
