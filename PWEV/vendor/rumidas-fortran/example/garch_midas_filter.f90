program garch_midas_filter
  use rumidas
  implicit none
  integer, parameter :: n=20,k=3
  real(dp) :: ret(n),mv(k+1,n),param(6)
  real(dp),allocatable::ll(:),cv(:),lr(:),sr(:)
  type(garch_midas_spec)::spec
  integer::i,status

  do i=1,n
    ret(i)=0.01_dp*sin(0.7_dp*real(i,dp))
    mv(:,i)=[0.1_dp+0.01_dp*i,0.2_dp+0.01_dp*i,0.3_dp+0.01_dp*i,0.4_dp+0.01_dp*i]
  end do
  param=[0.06_dp,0.86_dp,0.05_dp,-8.0_dp,0.25_dp,3.0_dp]
  spec=garch_midas_spec(RUMIDAS_GM,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,k,0,.true.)
  call garch_midas_evaluate(param,spec,ret,mv,ll,cv,lr,sr,status)
  if(status/=RUMIDAS_SUCCESS) error stop 'GARCH-MIDAS evaluation failed'
  print '(a,f12.4)', 'log likelihood: ',sum(ll)
  print '(a,3es14.5)', 'last long/short/conditional: ',lr(n),sr(n),cv(n)
end program garch_midas_filter
