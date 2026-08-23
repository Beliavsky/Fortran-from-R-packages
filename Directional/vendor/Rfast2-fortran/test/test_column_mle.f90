program test_column_mle
   use rfast2
   implicit none
   real(dp) :: xb(6,2),xp(6,2),u(5,2),a(2,3),c(2,3),pl(2,2),sp(2,2)

   xb(:,1)=[0.1_dp,0.2_dp,0.3_dp,0.5_dp,0.7_dp,0.8_dp]
   xb(:,2)=[0.2_dp,0.25_dp,0.4_dp,0.55_dp,0.65_dp,0.9_dp]
   a=col_beta_mle(xb)
   if (any(a(:,1:2) <= 0.0_dp)) error stop 1
   c=col_cauchy_mle(2.0_dp*xb-1.0_dp)
   if (any(c(:,3) <= 0.0_dp)) error stop 2
   xp(:,1)=[1.0_dp,1.2_dp,1.5_dp,2.0_dp,3.0_dp,5.0_dp]
   xp(:,2)=[2.0_dp,2.2_dp,2.6_dp,3.1_dp,4.0_dp,6.0_dp]
   pl=col_powerlaw_mle(xp)
   if (any(pl(:,1) <= 1.0_dp)) error stop 3
   u(:,1)=[0.1_dp,0.2_dp,0.4_dp,0.6_dp,0.8_dp]
   u(:,2)=[0.15_dp,0.25_dp,0.35_dp,0.55_dp,0.75_dp]
   sp=col_sp_mle(u)
   if (any(sp(:,1) <= 0.0_dp)) error stop 4
   print '(a)', 'test_column_mle: PASS'
end program test_column_mle
