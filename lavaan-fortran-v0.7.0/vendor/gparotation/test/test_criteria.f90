program test_criteria
  use gpa_kinds, only: dp
  use gpa_criteria
  implicit none
  real(dp) :: l(6,3)
  character(len=16), parameter :: names(10) = [character(len=16) :: &
       'varimax','quartimin','cf','geomin','entropy','oblimax','binormamin', &
       'tandemI','tandemII','bifactor']
  type(criterion_options) :: o
  type(criterion_result) :: r
  integer :: i, fails
  fails=0
  l(1,:)=[0.80_dp,0.20_dp,0.10_dp]
  l(2,:)=[0.75_dp,0.25_dp,0.05_dp]
  l(3,:)=[0.10_dp,0.82_dp,0.20_dp]
  l(4,:)=[0.15_dp,0.76_dp,0.25_dp]
  l(5,:)=[0.12_dp,0.18_dp,0.78_dp]
  l(6,:)=[0.20_dp,0.10_dp,0.73_dp]
  o%kappa=0.2_dp
  o%delta=0.01_dp
  do i=1,size(names)
    call evaluate_criterion(trim(names(i)),l,r,o)
    if(r%info/=0 .or. .not.(abs(r%f)<huge(1.0_dp))) then
      print *, 'criterion failure: ',trim(names(i)),r%info,r%f
      fails=fails+1
    end if
  end do
  call gradient_check('varimax',l,o,2.0e-5_dp,fails)
  call gradient_check('quartimin',l,o,2.0e-5_dp,fails)
  call gradient_check('cf',l,o,3.0e-5_dp,fails)
  call gradient_check('geomin',l,o,3.0e-5_dp,fails)
  call gradient_check('entropy',l,o,3.0e-5_dp,fails)
  call gradient_check('oblimin',l,o,3.0e-5_dp,fails)
  call gradient_check('infomax',l,o,8.0e-5_dp,fails)
  call gradient_check('mccammon',l,o,2.0e-4_dp,fails)
  call gradient_check('tandemI',l,o,3.0e-4_dp,fails)
  call gradient_check('tandemII',l,o,3.0e-4_dp,fails)
  call gradient_check('oblimax',l,o,8.0e-5_dp,fails)
  call gradient_check('bentler',l,o,2.0e-4_dp,fails)
  call gradient_check('binormamin',l,o,2.0e-4_dp,fails)
  call gradient_check('varimin',l,o,2.0e-5_dp,fails)
  call gradient_check('bifactor',l,o,3.0e-4_dp,fails)
  call gradient_check('bigeomin',l,o,6.0e-5_dp,fails)
  call check_target(l,fails)
  if(fails/=0) error stop 1
  print *, 'test_criteria: PASS'
contains
  subroutine gradient_check(name,x,o,tol,fails)
    character(len=*),intent(in)::name
    real(dp),intent(in)::x(:,:)
    type(criterion_options),intent(in)::o
    real(dp),intent(in)::tol
    integer,intent(inout)::fails
    type(criterion_result)::r,rp,rm
    real(dp)::xp(size(x,1),size(x,2)),h,num,err
    integer::i,j
    h=1.0e-6_dp
    call evaluate_criterion(name,x,r,o)
    err=0.0_dp
    do j=1,size(x,2)
    do i=1,size(x,1)
      xp=x
      xp(i,j)=xp(i,j)+h
      call evaluate_criterion(name,xp,rp,o)
      xp=x
      xp(i,j)=xp(i,j)-h
      call evaluate_criterion(name,xp,rm,o)
      num=(rp%f-rm%f)/(2.0_dp*h)
      err=max(err,abs(num-r%gq(i,j)))
    end do
    end do
    if(err>tol) then
    print *, 'gradient failure ',trim(name),err
    fails=fails+1
    end if
  end subroutine gradient_check
  subroutine check_target(x,fails)
    real(dp),intent(in)::x(:,:)
    integer,intent(inout)::fails
    type(criterion_options)::oo
    type(criterion_result)::rr
    allocate(oo%target(size(x,1),size(x,2)))
    oo%target=0.0_dp
    call evaluate_criterion('target',x,rr,oo)
    if(abs(rr%f-sum(x*x))>1.0e-12_dp .or. maxval(abs(rr%gq-2.0_dp*x))>1.0e-12_dp) then
      print *, 'target criterion failure'
      fails=fails+1
    end if
  end subroutine check_target
end program test_criteria
