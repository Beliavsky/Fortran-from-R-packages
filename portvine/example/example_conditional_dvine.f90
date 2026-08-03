program example_conditional_dvine
   use portvine, only : dp, conditional_dvine_sample
   use rvinecopulib, only : dvine_model,make_dvine,make_bicop,bicop_gaussian,bicop_indep
   implicit none
   type(dvine_model)::vine
   real(dp)::u(3,1000)
   integer::status
   vine=make_dvine(3,[1,2,3])
   vine%pair(1,2)=make_bicop(bicop_gaussian,parameters=[0.7_dp])
   vine%pair(2,3)=make_bicop(bicop_gaussian,parameters=[0.4_dp])
   vine%pair(1,3)=make_bicop(bicop_indep)
   call conditional_dvine_sample(vine,1000,[0.05_dp],.true.,u,status)
   print '(a,f10.6)','fixed first uniform: ',sum(u(1,:))/1000.0_dp
   print '(a,f10.6)','mean second uniform: ',sum(u(2,:))/1000.0_dp
end program example_conditional_dvine
