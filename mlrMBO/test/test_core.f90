program test_core
  use mlrmbo
  implicit none
  real(dp), allocatable :: p(:,:),ref(:),cand(:,:),front(:,:),v(:),w(:),ys(:)
  real(dp) :: hv,ei
  logical, allocatable :: dom(:)
  type(mbo_space) :: space
  real(dp) :: x(3)
  call init_space(space,[mbo_real,mbo_integer,mbo_categorical], [0.0_dp,1.0_dp,1.0_dp], &
    [1.0_dp,5.0_dp,3.0_dp],nlevels=[0,0,3])
  x=[-1.0_dp,3.7_dp,9.0_dp]; call space%repair(x)
  call check(abs(x(1))<1e-14_dp .and. nint(x(2))==4 .and. nint(x(3))==3,'space repair')

  p=reshape([1.0_dp,2.0_dp,3.0_dp,3.0_dp,2.0_dp,1.0_dp],[3,2])
  ! reshape is column-major -> rows are (1,3),(2,2),(3,1)
  ref=[4.0_dp,4.0_dp]
  hv=dominated_hypervolume(p,ref)
  call check(abs(hv-6.0_dp)<1e-12_dp,'hypervolume')
  dom=dominated_mask(p); call check(.not.any(dom),'nondomination')
  cand=reshape([1.5_dp,3.5_dp],[1,2])
  front=p
  v=eps_indicator_values(cand,front)
  call check(size(v)==1 .and. abs(v(1)-0.5_dp)<1e-12_dp,'epsilon indicator')

  ei=expected_improvement(0.0_dp,1.0_dp,-1.0_dp,1.0e-12_dp)
  call check(abs(ei-0.0833154705876863_dp)<1e-13_dp,'EI formula')

  w=[0.25_dp,0.75_dp]
  p=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2])
  call parego_scalarize(p,[.true.,.true.],w,0.05_dp,ys)
  call check(size(ys)==2 .and. all(ys>=0.0_dp),'parego scalarization')
  print *, 'test_core: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
  end subroutine check
end program test_core
