! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_types
  use nvmix_kinds, only : dp, i8
  implicit none
  private
  integer, parameter, public :: mix_constant=1, mix_inverse_gamma=2, mix_pareto=3, mix_gamma=4
  type, public :: integration_control
    integer :: samples=16384
    integer :: batches=8
    integer(i8) :: seed=5489_i8
  end type
  type, public :: probability_result
    real(dp) :: value=0.0_dp
    real(dp) :: error=0.0_dp
    integer :: evaluations=0
    logical :: ok=.true.
    character(len=160) :: message=''
  end type
  type, public :: vector_result
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: errors(:)
    logical :: ok=.true.
    character(len=160) :: message=''
  end type
  type, public :: nvmix_model
    real(dp), allocatable :: loc(:)
    real(dp), allocatable :: scale(:,:)
    integer, allocatable :: groupings(:)
    integer, allocatable :: mix_family(:)
    real(dp), allocatable :: mix_parameter(:)
  contains
    procedure :: dimension => model_dimension
    procedure :: groups => model_groups
  end type
  type, public :: sample_result
    real(dp), allocatable :: x(:,:)
    logical :: ok=.true.
    character(len=160) :: message=''
  end type
  type, public :: fit_result
    real(dp), allocatable :: loc(:)
    real(dp), allocatable :: scale(:,:)
    real(dp), allocatable :: mixing_parameter(:)
    real(dp) :: log_likelihood=-huge(1.0_dp)
    real(dp) :: aic=huge(1.0_dp)
    real(dp) :: bic=huge(1.0_dp)
    integer :: iterations=0
    logical :: converged=.false.
    logical :: ok=.true.
    character(len=160) :: message=''
  end type
  type, public :: qq_result
    real(dp), allocatable :: observed(:)
    real(dp), allocatable :: theoretical(:)
    logical :: ok=.true.
    character(len=160) :: message=''
  end type
contains
  integer function model_dimension(self) result(d)
    class(nvmix_model), intent(in) :: self
    if(allocated(self%loc)) then; d=size(self%loc); else; d=0; end if
  end function
  integer function model_groups(self) result(g)
    class(nvmix_model), intent(in) :: self
    if(allocated(self%mix_family)) then; g=size(self%mix_family); else; g=0; end if
  end function
end module nvmix_types
