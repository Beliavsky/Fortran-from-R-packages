program demo_relsurv
  use relsurv
  implicit none
  integer :: i
  type(ratetable_type)::tab
  type(rs_surv_result)::fit
  integer::dims(1),factor(1),ncuts(1),status(4)
  real(dp)::cuts(1,1),rate(1),x(4,1),time(4),times(4)
  dims=1; factor=1; ncuts=0; cuts=0.0_dp; rate=0.00005_dp
  tab=make_ratetable(dims,factor,cuts,ncuts,rate)
  x=1.0_dp; time=[100.0_dp,200.0_dp,300.0_dp,400.0_dp]; status=[1,0,1,0]; times=time
  call rs_surv(tab,x,time,status,times,fit,method=method_pohar_perme,precision=1.0_dp)
  print '(a)', ' time       net survival'
  print '(f8.1,2x,f12.7)', (fit%time(i),fit%surv(i),i=1,size(fit%time))
end program
