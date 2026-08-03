! SPDX-License-Identifier: GPL-2.0-or-later
program test_distributions
  use mixtools
  implicit none
  real(dp) :: v
  real(dp), allocatable :: root(:,:), points(:,:), perms_r(:,:)
  real(dp) :: a(2,2), q
  integer, allocatable :: perms(:,:)
  integer :: status

  allocate(root(0,0),points(0,0),perms(0,0))

  v = ddirichlet([0.2_dp,0.3_dp,0.5_dp],[1.0_dp,1.0_dp,1.0_dp])
  call check(abs(v-2.0_dp)<1.0e-12_dp,"Dirichlet density")
  a=0.0_dp;a(1,1)=1.0_dp;a(2,2)=1.0_dp
  v=dmvnorm([0.0_dp,0.0_dp],[0.0_dp,0.0_dp],a,status)
  call check(status==MIXTOOLS_SUCCESS.and.abs(v-1.0_dp/(2.0_dp*pi))<1.0e-12_dp,"MVN density")
  q=wquantile([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[1.0_dp,1.0_dp,1.0_dp,1.0_dp],0.5_dp,status)
  call check(status==0.and.abs(q-2.0_dp)<1.0e-12_dp,"weighted quantile")
  a=0.0_dp;a(1,1)=4.0_dp;a(2,2)=9.0_dp
  call matsqrt(a,root,status)
  call check(status==0,"matrix square-root status")
  if(status==0) call check(maxval(abs(matmul(root,root)-a))<1.0e-10_dp,"matrix square root")
  call ellipse([0.0_dp,0.0_dp],a,0.05_dp,16,points,status)
  call check(status==0,"ellipse status")
  if(status==0) call check(size(points,1)==16,"ellipse points")
  call perm(4,2,perms,status)
  call check(status==0,"permutation status")
  if(status==0)then
    call check(size(perms,1)==12.and.size(perms,2)==2,"permutations")
    allocate(perms_r(size(perms,1),size(perms,2)))
    perms_r=real(perms,dp)
    call check(minval(perms_r)>=1.0_dp.and.maxval(perms_r)<=4.0_dp,"permutation range")
  end if
  print '(a)', 'test_distributions: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine check
end program test_distributions
