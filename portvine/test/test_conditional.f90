program test_conditional
   use portvine, only : dp, conditional_dvine_sample
   use rvinecopulib, only : dvine_model, make_dvine, make_bicop, &
      bicop_gaussian, bicop_indep
   implicit none
   type(dvine_model) :: model
   real(dp) :: sample(3,250), z(3,250)
   integer :: status
   model=make_dvine(3,[1,2,3])
   model%pair(1,2)=make_bicop(bicop_gaussian,parameters=[0.65_dp])
   model%pair(2,3)=make_bicop(bicop_gaussian,parameters=[0.35_dp])
   model%pair(1,3)=make_bicop(bicop_indep)
   call conditional_dvine_sample(model,250,[0.2_dp],.false.,sample,status)
   if(status/=0 .or. maxval(abs(sample(1,:)-0.2_dp))>1.0e-8_dp)error stop 1
   call conditional_dvine_sample(model,250,[0.2_dp,0.7_dp],.true.,sample,status)
   call model%rosenblatt(sample,z)
   if(status/=0 .or. maxval(abs(sample(1,:)-0.2_dp))>1.0e-8_dp)error stop 2
   if(maxval(abs(z(2,:)-0.7_dp))>2.0e-7_dp)error stop 3
   print '(a)', 'test_conditional: PASS'
end program test_conditional
