! SPDX-License-Identifier: MIT
program test_core
  use bekks
  use bekks_model, only: simulate_bekk
  implicit none
  real(dp), parameter :: tol=2.0e-8_dp
  real(dp) :: data(8,2), signs(2), ll
  real(dp), allocatable :: theta(:),sim(:,:),h(:,:,:),innov(:,:)
  type(bekk_parameters) :: par
  type(rng_state) :: state
  integer :: st

  data=reshape([ &
    0.12_dp,-0.03_dp,0.07_dp,-0.15_dp,0.04_dp,-0.06_dp,0.10_dp,-0.02_dp, &
   -0.08_dp,0.11_dp,-0.02_dp,-0.05_dp,0.09_dp,-0.12_dp,0.03_dp,0.07_dp],[8,2])
  signs=-1.0_dp

  par%model_type=bekk_full;par%asymmetric=.false.
  allocate(par%c(2,2),par%a(2,2),par%b(2,2),par%g(2,2))
  par%c=reshape([0.2_dp,0.05_dp,0.0_dp,0.15_dp],[2,2])
  par%a=reshape([0.25_dp,0.02_dp,-0.01_dp,0.2_dp],[2,2])
  par%b=0.0_dp
  par%g=reshape([0.8_dp,0.01_dp,0.0_dp,0.75_dp],[2,2])
  theta=pack_parameters(par)
  ll=loglike_bekk(theta,data)
  call assert_close(ll,7.559158365670399_dp,tol,'full likelihood')
  if(.not.valid_bekk(theta,2))error stop 'full validity'

  par%asymmetric=.true.;par%b=reshape([0.1_dp,0.0_dp,0.0_dp,0.08_dp],[2,2])
  theta=pack_parameters(par)
  ll=loglike_asymm_bekk(theta,data,signs)
  call assert_close(ll,7.555407727379082_dp,tol,'asymmetric full likelihood')

  par%model_type=bekk_diagonal;par%asymmetric=.false.;par%a=0.0_dp;par%g=0.0_dp
  par%a(1,1)=0.25_dp;par%a(2,2)=0.2_dp;par%g(1,1)=0.8_dp;par%g(2,2)=0.75_dp
  theta=pack_parameters(par)
  call assert_close(loglike_dbekk(theta,data),7.559188653634845_dp,tol,'diagonal likelihood')

  par%model_type=bekk_scalar;par%asymmetric=.false.;par%a_scalar=0.1_dp;par%g_scalar=0.8_dp
  theta=pack_parameters(par)
  call assert_close(loglike_sbekk(theta,data),5.560859831717117_dp,tol,'scalar likelihood')

  par%asymmetric=.true.;par%a_scalar=0.1_dp;par%b_scalar=0.05_dp;par%g_scalar=0.75_dp
  theta=pack_parameters(par)
  call assert_close(loglike_asymm_sbekk(theta,data,signs),6.051637028692405_dp,tol,'asymmetric scalar likelihood')

  par%model_type=bekk_full;par%asymmetric=.false.
  par%a=reshape([0.25_dp,0.02_dp,-0.01_dp,0.2_dp],[2,2]);par%b=0.0_dp
  par%g=reshape([0.8_dp,0.01_dp,0.0_dp,0.75_dp],[2,2])
  theta=pack_parameters(par)
  allocate(innov(5,2))
  innov=reshape([0.2_dp,1.0_dp,-0.5_dp,0.1_dp,0.8_dp, &
                -0.4_dp,0.3_dp,0.7_dp,-1.2_dp,0.2_dp],[5,2])
  call simulate_bekk(theta,5,2,bekk_full,.false.,state=state,data=sim,h=h,status=st,innovations=innov)
  if(st/=bekk_ok)error stop 'simulation status'
  call assert_close(sim(5,1),0.27864212_dp,2.0e-7_dp,'simulation y1')
  call assert_close(sim(5,2),0.10804287_dp,2.0e-7_dp,'simulation y2')
  call assert_close(h(1,1,5),0.12131474_dp,2.0e-7_dp,'simulation h11')
  call assert_close(h(1,2,5),0.02652982_dp,2.0e-7_dp,'simulation h12')

  print '(a)', 'test_core: PASS'
contains
  subroutine assert_close(x,y,eps,label)
    real(dp),intent(in)::x,y,eps
    character(len=*),intent(in)::label
    if(abs(x-y)>eps*max(1.0_dp,abs(y)))then
      print '(a,2es24.15)',trim(label)//': ',x,y
      error stop 1
    end if
  end subroutine
end program test_core
