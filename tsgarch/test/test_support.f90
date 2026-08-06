module test_support
  use ghyp_kinds, only : dp
  use tsgarch
  implicit none
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then
      write(*,'(a)')'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  function standard_spec(model,distribution) result(spec)
    character(len=*),intent(in)::model
    character(len=*),intent(in),optional::distribution
    type(garch_spec)::spec
    spec%model=model
    spec%p=1
    spec%q=1
    spec%constant=.true.
    if(present(distribution))spec%distribution=distribution
  end function standard_spec

  function standard_parameters(y,spec) result(par)
    real(dp),intent(in)::y(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters)::par
    par=initialize_parameters(y,spec)
    par%mu=0.0_dp
    par%omega=0.02_dp
    par%alpha=0.06_dp
    par%beta=0.90_dp
    select case(trim(spec%model))
    case('gjrgarch')
    par%alpha=0.05_dp
    par%gamma=0.05_dp
    par%beta=0.86_dp
    case('egarch')
    par%omega=-0.12_dp
    par%alpha=-0.05_dp
    par%gamma=0.12_dp
    par%beta=0.92_dp
    case('aparch')
    par%alpha=0.07_dp
    par%gamma=0.10_dp
    par%beta=0.86_dp
    par%delta=1.5_dp
    case('fgarch')
    par%alpha=0.06_dp
    par%gamma=0.10_dp
    par%eta=0.05_dp
    par%beta=0.86_dp
    par%delta=1.5_dp
    case('cgarch')
    par%omega=0.01_dp
    par%alpha=0.04_dp
    par%beta=0.80_dp
    par%rho=0.95_dp
    par%phi=0.03_dp
    case('igarch')
    par%omega=0.01_dp
    par%alpha=0.06_dp
    par%beta=0.94_dp
    case('ewma')
    par%omega=0.0_dp
    par%alpha=0.06_dp
    par%beta=0.94_dp
    end select
  end function standard_parameters
end module test_support
