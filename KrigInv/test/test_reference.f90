program test_reference
  use kriginv
  use kriginv_math, only : bvn_cdf
  implicit none
  type(krig_model) :: m
  type(krig_prediction) :: p
  real(dp) :: x(5,2), y(5), ell(2), q(4,2), noise(5)
  real(dp), parameter :: mean_ref(4)=[ &
    1.46818996502313226e-1_dp,-3.77905199201018505e-2_dp,4.39308488800953634e-1_dp,6.40556214912274768e-1_dp]
  real(dp), parameter :: sd_ref(4)=[ &
    4.77487983019966700e-1_dp,4.05735450489121185e-1_dp,5.11112621961744784e-1_dp,6.49247937770970673e-1_dp]
  logical :: ok
  real(dp) :: bv
  x=reshape([0.0_dp,0.0_dp, 0.2_dp,0.8_dp, 0.5_dp,0.3_dp, 0.7_dp,0.9_dp, 1.0_dp,0.1_dp],[5,2],order=[2,1])
  y=[0.1_dp,0.7_dp,-0.2_dp,0.4_dp,1.0_dp]
  ell=[0.45_dp,0.7_dp]; noise=[0.01_dp,0.02_dp,0.015_dp,0.01_dp,0.025_dp]
  q=reshape([0.1_dp,0.2_dp,0.4_dp,0.4_dp,0.8_dp,0.2_dp,0.9_dp,0.8_dp],[4,2],order=[2,1])
  call init_krig_model(m,x,y,ell,variance=1.3_dp,nugget=0.05_dp,noise=noise,covariance='matern5_2',trend_order=1,ok=ok)
  call check(ok,'reference model')
  p=predict_nobias_km(m,q,'UK',.true.)
  call check(maxval(abs(p%mean-mean_ref))<2.0e-14_dp,'independent UK mean reference')
  call check(maxval(abs(p%sd-sd_ref))<2.0e-14_dp,'independent UK sd reference')
  bv=bvn_cdf(-1.3258790811200347_dp,-0.813229195158315_dp,0.7150960496290253_dp)
  call check(abs(bv-0.06749963151140081_dp)<2.0e-12_dp,'bivariate normal randomized reference')
  print '(a)', 'test_reference: PASS'
contains
  subroutine check(cond,msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if(.not.cond) then
      print '(a)', 'FAIL: '//trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_reference
