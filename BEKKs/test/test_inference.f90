! SPDX-License-Identifier: MIT
program test_inference
  use iso_fortran_env, only: int64
  use bekks
  implicit none

  real(dp), parameter :: tol = 2.0e-5_dp
  type(rng_state) :: state
  real(dp), allocatable :: theta(:), score(:,:), hessian(:,:), covariance(:,:), robust(:,:)
  real(dp), allocatable :: random_theta(:), simulated(:,:), h(:,:,:), innovations(:,:)
  real(dp) :: data(160,2), signs(2), likelihood, random_likelihood
  real(dp) :: tp(5), tm(5), gradient_fd(5), gradient_score(5), step
  real(dp) :: a(2,2), ainv(2,2), recovered(2,2)
  integer :: i, j, status, info

  do i=1,size(data,1)
    data(i,1)=0.16_dp*sin(0.17_dp*real(i,dp))+0.04_dp*cos(0.041_dp*real(i*i,dp))
    data(i,2)=0.10_dp*cos(0.13_dp*real(i,dp))-0.03_dp*sin(0.071_dp*real(i*i,dp))
  end do
  signs=[-1.0_dp,1.0_dp]

  call grid_search_sbekk(data,theta,likelihood,status)
  call assert_true(status==bekk_ok,'deterministic starting values')
  call assert_true(size(theta)==5,'scalar parameter count')
  call assert_true(likelihood>-huge(1.0_dp)/2.0_dp,'finite starting likelihood')

  call score_sbekk(theta,data,score,status)
  call assert_true(status==bekk_ok,'score status')
  gradient_score=sum(score,dim=1)
  do j=1,size(theta)
    step=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(theta(j)))
    tp=theta;tm=theta
    tp(j)=tp(j)+step;tm(j)=tm(j)-step
    gradient_fd(j)=(loglike_sbekk(tp,data)-loglike_sbekk(tm,data))/(2.0_dp*step)
  end do
  call assert_true(maxval(abs(gradient_score-gradient_fd))<tol,'score finite differences')

  call hesse_sbekk(theta,data,hessian,status)
  call assert_true(status==bekk_ok,'Hessian status')
  call assert_true(maxval(abs(hessian-transpose(hessian)))<1.0e-10_dp,'Hessian symmetry')

  call qml_covariance(theta,data,bekk_scalar,.false.,covariance=covariance, &
    robust=robust,status=status)
  call assert_true(status==bekk_ok,'QML covariance status')
  call assert_true(maxval(abs(covariance-transpose(covariance)))<1.0e-8_dp, &
    'OPG covariance symmetry')
  call assert_true(maxval(abs(robust-transpose(robust)))<1.0e-8_dp, &
    'sandwich covariance symmetry')

  call rng_seed(state,918273_int64)
  call random_grid_search_asymmetric_sbekk(data,signs,state,random_theta, &
    random_likelihood,status,n_trials=12)
  call assert_true(status==bekk_ok,'random starting values')
  call assert_true(size(random_theta)==6,'asymmetric scalar parameter count')
  call assert_true(random_likelihood>-huge(1.0_dp)/2.0_dp,'random starting likelihood')

  allocate(innovations(7,2))
  innovations(:,1)=[0.2_dp,-0.1_dp,0.4_dp,-0.3_dp,0.0_dp,0.5_dp,-0.2_dp]
  innovations(:,2)=[-0.4_dp,0.3_dp,0.1_dp,-0.2_dp,0.6_dp,-0.1_dp,0.2_dp]
  call rng_seed(state,1_int64)
  call simulate_bekk_model(theta,7,2,bekk_scalar,.false.,state,simulated,h,status, &
    innovations=innovations)
  call assert_true(status==bekk_ok,'fixed-innovation public simulation')
  call assert_true(all(shape(simulated)==[7,2]),'simulation dimensions')
  call assert_true(all(shape(h)==[2,2,7]),'covariance dimensions')

  a=reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp],[2,2])
  call general_inverse(a,ainv,info)
  call assert_true(info==0,'general inverse status')
  recovered=matmul(a,matmul(ainv,a))
  call assert_true(maxval(abs(recovered-a))<1.0e-9_dp,'Moore-Penrose recovery')

  print '(a)', 'test_inference: PASS'

contains

  subroutine assert_true(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition)then
      print '(a)', 'FAILED: '//message
      error stop 1
    end if
  end subroutine assert_true

end program test_inference
