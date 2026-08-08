program test_custom_bounds
  use ao
  implicit none
  type(ao_result) :: r
  type(ao_options) :: o
  real(dp) :: x0(4), lo(4), hi(4)
  integer :: j
  x0=[2.0_dp,2.0_dp,2.0_dp,2.0_dp]
  lo=-1.0_dp; hi=3.0_dp
  o%partition=AO_PARTITION_CUSTOM
  allocate(o%custom_partition(3))
  allocate(o%custom_partition(1)%index(2)); o%custom_partition(1)%index=[1,2]
  allocate(o%custom_partition(2)%index(1)); o%custom_partition(2)%index=[3]
  allocate(o%custom_partition(3)%index(2)); o%custom_partition(3)%index=[2,4]
  o%iteration_limit=20
  call ao_optimize(quad,x0,r,o,lo,hi)
  if(r%value>1.0e-10_dp) error stop 'custom partition failed'
  if(any(r%estimate<lo) .or. any(r%estimate>hi)) error stop 'bounds violated'
  do j=1,r%details%n
    if(any(r%details%parameter(:,j)<lo) .or. any(r%details%parameter(:,j)>hi)) error stop 'history bounds violated'
  end do
  print *, 'PASS test_custom_bounds',r%value
contains
  function quad(x) result(f)
    real(dp),intent(in)::x(:); real(dp)::f
    f=sum((x-[0.25_dp,-0.5_dp,1.0_dp,0.75_dp])**2)
  end function
end program
