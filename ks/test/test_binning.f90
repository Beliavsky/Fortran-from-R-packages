! SPDX-License-Identifier: GPL-2.0-only
program test_binning
  use ks, only: dp, linear_binning, grid_interpolate, symconv_1d
  implicit none
  real(dp)::x(5,2),a(2),b(2),w(5),pts(3,2),vals(3)
  real(dp),allocatable::counts(:),fun(:),z(:)
  integer::m(2),i,j,k
  x=reshape([0.10_dp,0.20_dp, 0.40_dp,0.55_dp, 0.80_dp,0.75_dp, 0.33_dp,0.91_dp, 0.62_dp,0.12_dp],[5,2],order=[2,1])
  a=[0.0_dp,0.0_dp];b=[1.0_dp,1.0_dp];m=[4,5];w=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  call linear_binning(x,a,b,m,counts,w)
  if(abs(sum(counts)-sum(w))>2e-13_dp) error stop 'binning mass'
  allocate(fun(product(m)));k=0
  do j=0,m(2)-1
    do i=0,m(1)-1
      k=k+1;fun(k)=(a(1)+(b(1)-a(1))*real(i,dp)/real(m(1)-1,dp)) + &
                    2.0_dp*(a(2)+(b(2)-a(2))*real(j,dp)/real(m(2)-1,dp))
    end do
  end do
  pts=reshape([0.15_dp,0.22_dp, 0.51_dp,0.64_dp, 0.88_dp,0.05_dp],[3,2],order=[2,1])
  call grid_interpolate(pts,a,b,m,fun,vals)
  do i=1,3
    if(abs(vals(i)-(pts(i,1)+2.0_dp*pts(i,2)))>2e-13_dp) error stop 'grid interpolation'
  end do
  call symconv_1d([1.0_dp,2.0_dp,3.0_dp],[4.0_dp,5.0_dp,6.0_dp],z)
  if(any(abs(z)>huge(1.0_dp)/2.0_dp)) error stop 'symconv overflow'
  print *, 'test_binning: PASS'
end program
