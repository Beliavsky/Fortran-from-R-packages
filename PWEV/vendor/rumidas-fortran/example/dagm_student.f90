program dagm_student
  use rumidas
  implicit none
  integer,parameter::n=16,k=2
  real(dp)::ret(n),mv(k+1,n),param(9)
  real(dp),allocatable::ll(:),cv(:),lr(:),sr(:)
  type(garch_midas_spec)::spec
  integer::i,status

  do i=1,n
    ret(i)=0.015_dp*cos(0.5_dp*real(i,dp))
    mv(:,i)=[-0.2_dp+0.01_dp*i,0.15_dp-0.02_dp*i,0.25_dp+0.01_dp*i]
  end do
  param=[0.05_dp,0.84_dp,0.06_dp,-8.0_dp,0.25_dp,2.0_dp,-0.15_dp,3.0_dp,7.0_dp]
  spec=garch_midas_spec(RUMIDAS_DAGM,RUMIDAS_STUDENT_T,RUMIDAS_BETA_LAG,k,0,.true.)
  call garch_midas_evaluate(param,spec,ret,mv,ll,cv,lr,sr,status)
  if(status/=0) error stop 'DAGM evaluation failed'
  print '(a,f12.4)', 'Student-t DAGM log likelihood: ',sum(ll)
  print '(a,es14.5)', 'last conditional volatility: ',cv(n)
end program dagm_student
