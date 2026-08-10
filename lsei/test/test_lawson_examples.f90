program test_lawson_examples
   use lsei
   implicit none
   real(dp) :: a(4,2), b(4), e(3,2), f(3)
   real(dp) :: al(2,2), bl(2), c(1,2), d(1), el(3,2), fl(3)
   type(ls_result) :: r
   real(dp), parameter :: tol=2.0e-7_dp
   a(:,1)=[.25_dp,.5_dp,.5_dp,.8_dp]; a(:,2)=1.0_dp
   b=[.5_dp,.6_dp,.7_dp,1.2_dp]
   e(:,1)=[1.0_dp,0.0_dp,-1.0_dp]; e(:,2)=[0.0_dp,1.0_dp,-1.0_dp]
   f=[0.0_dp,0.0_dp,-1.0_dp]
   call lsi_solve(a,b,e,f,r)
   if (.not.r%succeeded()) error stop 1
   if (maxval(abs(r%x-[.6213152_dp,.3786848_dp]))>tol) error stop 2
   el(:,1)=[-.207_dp,-.392_dp,.599_dp]; el(:,2)=[2.558_dp,-1.351_dp,-1.206_dp]
   fl=[-1.3_dp,-.084_dp,.384_dp]
   call ldp_solve(el,fl,r)
   if (.not.r%succeeded()) error stop 3
   if (maxval(abs(r%x-[.1268538_dp,-.2554018_dp]))>tol) error stop 4
   al(:,1)=[.4302_dp,.6246_dp]; al(:,2)=[.3516_dp,.3384_dp]
   bl=[.6593_dp,.9666_dp]; c(1,:)=[.4087_dp,.1593_dp]; d=[.1376_dp]
   call lsei_solve(al,bl,c=c,d=d,res=r)
   if (.not.r%succeeded()) error stop 5
   if (maxval(abs(r%x-[-1.177499_dp,3.884770_dp]))>2.0e-6_dp) error stop 6
   print *, 'PASS test_lawson_examples'
end program
