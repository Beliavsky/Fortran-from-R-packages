program lawson_examples
   use lsei
   implicit none
   real(dp) :: a(4,2),b(4),e(3,2),f(3)
   type(ls_result) :: r
   a(:,1)=[.25_dp,.5_dp,.5_dp,.8_dp]; a(:,2)=1.0_dp
   b=[.5_dp,.6_dp,.7_dp,1.2_dp]
   e(:,1)=[1.0_dp,0.0_dp,-1.0_dp]; e(:,2)=[0.0_dp,1.0_dp,-1.0_dp]
   f=[0.0_dp,0.0_dp,-1.0_dp]
   call lsi_solve(a,b,e,f,r)
   print '(a,2f14.8)', 'LSI solution: ',r%x
end program
