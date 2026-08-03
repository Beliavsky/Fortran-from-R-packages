program mem_midas
  use rumidas
  implicit none
  integer,parameter::n=20,k=2
  real(dp)::x(n),ret(n),mv(k+1,n),param(6)
  real(dp),allocatable::ll(:),pred(:),lr(:),sr(:)
  type(mem_spec)::spec
  integer::i,status

  do i=1,n
    x(i)=1.0_dp+0.15_dp*sin(0.8_dp*real(i,dp))
    ret(i)=0.01_dp*cos(0.6_dp*real(i,dp))
    mv(:,i)=[0.1_dp+0.01_dp*i,0.2_dp+0.01_dp*i,0.3_dp+0.01_dp*i]
  end do
  param=[0.20_dp,0.65_dp,0.05_dp,0.0_dp,0.15_dp,2.5_dp]
  spec=mem_spec(RUMIDAS_MEM_MIDAS,k,.true.)
  call mem_evaluate(param,spec,x,ll,pred,lr,sr,status,daily_ret=ret,mv_m=mv)
  if(status/=0) error stop 'MEM-MIDAS evaluation failed'
  print '(a,f12.4)', 'MEM-MIDAS log likelihood: ',sum(ll)
  print '(a,3f12.6)', 'last long/short/prediction: ',lr(n),sr(n),pred(n)
end program mem_midas
