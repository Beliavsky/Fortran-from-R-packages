program demo_compound_cox
  use compound_cox
  implicit none
  real(dp) :: t(8), x(8,2)
  integer :: d(8)
  type(univariate_result) :: u
  type(cg_result) :: cg
  t=[1.0_dp,3.0_dp,5.0_dp,4.0_dp,7.0_dp,8.0_dp,10.0_dp,13.0_dp]
  d=[1,0,0,1,1,0,1,0]
  x(:,1)=[0.2_dp,-0.1_dp,0.4_dp,0.1_dp,0.7_dp,-0.2_dp,0.9_dp,0.0_dp]
  x(:,2)=[1.0_dp,0.4_dp,0.2_dp,0.3_dp,-0.1_dp,-0.5_dp,-0.6_dp,-0.8_dp]
  call uni_score(t,d,x,u)
  call cg_clayton(t,d,2.0_dp,cg)
  print '(a,2f12.5)', 'score beta: ',u%beta
  print '(a,f12.5)', 'Clayton tau: ',cg%tau
end program demo_compound_cox
