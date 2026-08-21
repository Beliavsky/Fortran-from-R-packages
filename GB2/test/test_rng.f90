program test_rng
  use gb2_kinds, only : dp
  use gb2_distribution, only : rgb2
  implicit none
  integer, parameter :: n=20000
  integer :: ns,i
  integer, allocatable :: seed(:)
  real(dp), allocatable :: x(:)
  real(dp) :: empirical,target

  call random_seed(size=ns)
  allocate(seed(ns),x(n))
  do i=1,ns
    seed(i)=13579+104729*i
  end do
  call random_seed(put=seed)

  call rgb2(n,2.3_dp,4.2_dp,1.7_dp,3.4_dp,x)
  call assert_true(all(x>0.0_dp),'all draws positive')
  empirical=sum(x)/real(n,dp)
  target=3.1980398622384523_dp
  call assert_true(abs(empirical-target)<0.08_dp,'sample mean')
  print '(a)', 'test_rng: PASS'
contains
  subroutine assert_true(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if(.not.ok) then
      print '(a,1x,a)', 'FAIL:',msg
      error stop 1
    end if
  end subroutine assert_true
end program test_rng
