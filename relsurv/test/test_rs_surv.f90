program test_rs_surv
  use relsurv
  implicit none
  type(ratetable_type)::tab
  type(rs_surv_result)::res
  integer::dims(1),factor(1),ncuts(1),status(2)
  real(dp)::cuts(1,1),rates(1),x(2,1),y(2),times(2),target(2),lambda
  dims=1; factor=1; ncuts=0; cuts=0.0_dp; lambda=0.02_dp; rates=lambda
  tab=make_ratetable(dims,factor,cuts,ncuts,rates); x=1.0_dp
  y=[1.0_dp,2.0_dp]; status=[1,0]; times=y
  call rs_surv(tab,x,y,status,times,res,method=method_ederer1)
  target=[0.5_dp/exp(-lambda),0.5_dp/exp(-2.0_dp*lambda)]
  if(maxval(abs(res%surv-target))>2.0e-12_dp) then
    print *, 'FAIL Ederer I',res%surv,target; error stop 1
  end if
  call rs_surv(tab,x,y,status,times,res,method=method_ederer2)
  if(any(.not.(res%surv>0.0_dp))) error stop 1
  print *, 'test_rs_surv: PASS'
end program
