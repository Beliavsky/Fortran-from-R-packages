! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_types
  use mvtnorm_kinds, only : dp
  implicit none

  integer, parameter :: method_genz_bretz = 1
  integer, parameter :: method_tvpack = 2
  integer, parameter :: method_miwa = 3

  type :: probability_control
    integer :: method = method_genz_bretz
    integer :: maxpts = 25000
    integer :: batches = 12
    integer :: miwa_steps = 32
    real(dp) :: abseps = 1.0e-3_dp
    real(dp) :: releps = 0.0_dp
    integer :: seed = 12345
  end type probability_control

  type :: probability_result
    real(dp) :: value = 0.0_dp
    real(dp) :: error = 1.0_dp
    integer :: inform = 0
    integer :: evaluations = 0
    character(len=256) :: message = ''
  end type probability_result

  type :: quantile_result
    real(dp) :: quantile = 0.0_dp
    real(dp) :: probability = 0.0_dp
    real(dp) :: error = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
    character(len=256) :: message = ''
  end type quantile_result

  type :: likelihood_result
    real(dp), allocatable :: loglik(:)
    real(dp), allocatable :: score(:,:)
    real(dp) :: total = 0.0_dp
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type likelihood_result

  type :: mvnormal_model
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: chol(:,:)
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type mvnormal_model

  type :: conditional_result
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: covariance(:,:)
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type conditional_result

contains

  function genz_bretz(maxpts,abseps,releps,seed,batches) result(control)
    integer,intent(in),optional::maxpts,seed,batches
    real(dp),intent(in),optional::abseps,releps
    type(probability_control)::control
    control%method=method_genz_bretz
    if(present(maxpts)) control%maxpts=maxpts
    if(present(abseps)) control%abseps=abseps
    if(present(releps)) control%releps=releps
    if(present(seed)) control%seed=seed
    if(present(batches)) control%batches=batches
  end function genz_bretz

  function tvpack(abseps) result(control)
    real(dp),intent(in),optional::abseps
    type(probability_control)::control
    control%method=method_tvpack
    control%maxpts=200000
    control%batches=24
    if(present(abseps)) control%abseps=abseps
  end function tvpack

  function miwa(steps,abseps) result(control)
    integer,intent(in),optional::steps
    real(dp),intent(in),optional::abseps
    type(probability_control)::control
    control%method=method_miwa
    if(present(steps)) control%miwa_steps=steps
    if(present(abseps)) control%abseps=abseps
  end function miwa

end module mvtnorm_types
