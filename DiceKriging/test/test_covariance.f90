! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
program test_covariance
  use dicekriging
  implicit none
  real(dp) :: x(3,2), expected, dx1, dx2
  real(dp), allocatable :: c(:,:), p(:), lo(:), up(:)
  type(covariance_model) :: cov
  type(scaling_axis) :: ax(2)

  x=reshape([0.1_dp,0.4_dp,0.8_dp, 0.2_dp,0.7_dp,0.3_dp],[3,2])
  cov%kind=cov_gauss; cov%sd2=1.5_dp; allocate(cov%range(2)); cov%range=[0.4_dp,0.7_dp]
  call covariance_matrix(cov,x,c)
  dx1=(x(1,1)-x(2,1))/(cov%range(1)/(sqrt(2.0_dp)/2.0_dp))
  dx2=(x(1,2)-x(2,2))/(cov%range(2)/(sqrt(2.0_dp)/2.0_dp))
  expected=cov%sd2*exp(-(dx1*dx1+dx2*dx2))
  if(abs(c(1,2)-expected)>1.0e-14_dp) error stop 'Gaussian covariance mismatch'
  if(maxval(abs(c-transpose(c)))>1.0e-14_dp) error stop 'covariance not symmetric'

  cov%nugget_flag=.true.; cov%nugget=0.05_dp
  call covariance_matrix(cov,x,c)
  if(maxval(abs([c(1,1),c(2,2),c(3,3)]-(cov%sd2+cov%nugget)))>1.0e-14_dp) error stop 'nugget mismatch'

  cov%nugget_flag=.false.; cov%kind=cov_powexp; allocate(cov%shape(2)); cov%shape=[1.3_dp,1.8_dp]
  call covariance_bounds(cov,x,lo,up)
  if(size(lo)/=4 .or. maxval(abs(up(3:4)-2.0_dp))>1.0e-14_dp) error stop 'powexp bounds mismatch'
  call get_cov_params(cov,2,p)
  if(maxval(abs(p-[0.4_dp,0.7_dp,1.3_dp,1.8_dp]))>1.0e-14_dp) error stop 'parameter flatten mismatch'

  cov%kind=cov_matern52; cov%iso=.true.; deallocate(cov%shape)
  call set_cov_params(cov,2,[0.6_dp])
  if(maxval(abs(cov%range-0.6_dp))>1.0e-14_dp) error stop 'isotropic parameter mismatch'

  allocate(ax(1)%knots(2),ax(1)%eta(2),ax(2)%knots(2),ax(2)%eta(2))
  ax(1)%knots=[0.0_dp,1.0_dp]; ax(1)%eta=[2.0_dp,2.0_dp]
  ax(2)%knots=[0.0_dp,1.0_dp]; ax(2)%eta=[3.0_dp,3.0_dp]
  if(abs(scaling_fun1d(0.4_dp,ax(1)%knots,ax(1)%eta)-0.8_dp)>1.0e-14_dp) error stop 'scaling value mismatch'
  if(abs(scaling_fun1d_dx(0.4_dp,ax(2)%knots,ax(2)%eta)-3.0_dp)>1.0e-14_dp) error stop 'scaling derivative mismatch'

  print *, 'test_covariance: PASS'
end program test_covariance
