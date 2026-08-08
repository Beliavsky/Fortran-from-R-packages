! SPDX-License-Identifier: GPL-2.0-only
module calibrar_objective
  use calibrar_kinds, only : dp
  use calibrar_fitness, only : fitness_norm2, fitness_lnorm2, fitness_lnorm3, &
    fitness_lnorm4, fitness_lnorm4b, fitness_pois, fitness_normp, fitness_penalty
  implicit none
  private
  public :: calibration_term, calibration_objective_value

  type :: calibration_term
    integer :: first = 1
    integer :: last = 0
    character(len=16) :: fitness_type = "norm2"
    real(dp) :: weight = 1.0_dp
    logical :: calibrate = .true.
  end type calibration_term

contains

  subroutine calibration_objective_value(obs, sim, terms, values, aggregate, total)
    real(dp), intent(in) :: obs(:), sim(:)
    type(calibration_term), intent(in) :: terms(:)
    real(dp), intent(out) :: values(:)
    logical, intent(in), optional :: aggregate
    real(dp), intent(out), optional :: total
    integer :: i,a,b
    logical :: agg
    real(dp) :: s
    if(size(obs)/=size(sim)) error stop "calibration_objective_value: size mismatch"
    if(size(values)/=size(terms)) error stop "calibration_objective_value: values size mismatch"
    values=0.0_dp
    do i=1,size(terms)
      if(.not.terms(i)%calibrate) cycle
      a=terms(i)%first; b=terms(i)%last
      if(a<1 .or. b>size(obs) .or. b<a) error stop "calibration_objective_value: invalid term range"
      select case(trim(terms(i)%fitness_type))
      case("norm2"); values(i)=fitness_norm2(obs(a:b),sim(a:b))
      case("lnorm2"); values(i)=fitness_lnorm2(obs(a:b),sim(a:b))
      case("lnorm3"); values(i)=fitness_lnorm3(obs(a:b),sim(a:b))
      case("lnorm4"); values(i)=fitness_lnorm4(obs(a:b),sim(a:b))
      case("lnorm4b"); values(i)=fitness_lnorm4b(obs(a:b),sim(a:b))
      case("pois"); values(i)=fitness_pois(obs(a:b),sim(a:b))
      case("normp","re"); values(i)=fitness_normp(sim(a:b))
      case("penalty","penalty0","penalty1","penalty2"); values(i)=fitness_penalty(sim(a:b))
      case default; error stop "calibration_objective_value: unsupported fitness type"
      end select
    end do
    agg=.false.;if(present(aggregate))agg=aggregate
    if(present(total))then
      s=0.0_dp
      do i=1,size(terms)
        if(terms(i)%calibrate)s=s+terms(i)%weight*values(i)
      end do
      total=s
    else if(agg)then
      continue
    end if
  end subroutine calibration_objective_value
end module calibrar_objective
