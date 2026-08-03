program test_garch_variants
  use rumidas
  implicit none
  integer, parameter :: n=8, k1=2, k2=3
  real(dp) :: ret(n), x(n), mv1(k1+1,n), mv2(k2+1,n)
  real(dp) :: p_gmx(7), p_gm2(7), p_dagm(8), p_dagm2(12)
  real(dp), allocatable :: ll(:), cv(:), lr(:), sr(:), a(:)
  type(garch_midas_spec) :: spec
  integer :: i, status

  do i=1,n
    ret(i)=0.012_dp*sin(real(i,dp))
    x(i)=0.02_dp+0.001_dp*real(i,dp)
    mv1(:,i)=[-0.2_dp+0.01_dp*i, 0.1_dp+0.02_dp*i, 0.3_dp+0.01_dp*i]
    mv2(:,i)=[0.1_dp*i, -0.05_dp*i, 0.02_dp*i, 0.03_dp*i]
  end do

  p_gmx=[0.05_dp,0.80_dp,0.04_dp,0.10_dp,-3.0_dp,0.2_dp,2.5_dp]
  spec=garch_midas_spec(RUMIDAS_GMX,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,k1,0,.true.)
  call garch_midas_evaluate(p_gmx,spec,ret,mv1,ll,cv,lr,sr,status,x_variable=x)
  call check(status==0,'GMX status')
  call check(all(cv>0.0_dp),'GMX positive')

  p_gm2=[0.05_dp,0.85_dp,-3.0_dp,0.2_dp,2.2_dp,-0.1_dp,3.0_dp]
  spec=garch_midas_spec(RUMIDAS_GM2M,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,k1,k2,.false.)
  call garch_midas_evaluate(p_gm2,spec,ret,mv1,ll,cv,lr,sr,status,mv_m_2=mv2)
  call check(status==0,'GM2M status')
  call gm_2m_long_run_vol_no_skew(p_gm2,ret,mv1,mv2,k1,k2,a,status)
  call check(maxval(abs(a-lr))<1.0e-13_dp,'GM2M wrapper')

  p_dagm=[0.05_dp,0.82_dp,0.04_dp,-3.0_dp,0.3_dp,2.0_dp,-0.2_dp,3.0_dp]
  spec=garch_midas_spec(RUMIDAS_DAGM,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,k1,0,.true.)
  call garch_midas_evaluate(p_dagm,spec,ret,mv1,ll,cv,lr,sr,status)
  call check(status==0,'DAGM status')
  call check(garch_midas_parameter_count(spec)==8,'DAGM count')

  p_dagm2=[0.04_dp,0.82_dp,0.04_dp,-3.0_dp,0.2_dp,2.0_dp,-0.1_dp,2.5_dp, &
           0.15_dp,2.2_dp,-0.05_dp,3.0_dp]
  spec=garch_midas_spec(RUMIDAS_DAGM2M,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,k1,k2,.true.)
  call garch_midas_evaluate(p_dagm2,spec,ret,mv1,ll,cv,lr,sr,status,mv_m_2=mv2)
  call check(status==0,'DAGM2M status')
  call check(all(abs(cv*cv-sr*lr*lr)<1.0e-12_dp),'variance factorization')

  print '(a)', 'test_garch_variants: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) error stop message
  end subroutine check
end program test_garch_variants
