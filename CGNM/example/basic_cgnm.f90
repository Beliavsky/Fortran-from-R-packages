program basic_cgnm
   use cgnm, only : dp,cgnm_problem,cgnm_options,cgnm_result,cgnm_init_problem,cgnm_fit
   implicit none
   type(cgnm_problem)::prob
   type(cgnm_options)::opt
   type(cgnm_result)::res
   real(dp)::target(3),lo(2),hi(2),tt(3)
   integer::ierr,i
   tt=[0._dp,1._dp,2._dp]; target=2._dp+0.5_dp*tt
   lo=[-3._dp,-2._dp]; hi=[5._dp,3._dp]
   call cgnm_init_problem(prob,model,target,lo,hi,ierr=ierr)
   opt%num_minimizers=40; opt%num_iterations=12; opt%seed=7
   call cgnm_fit(prob,opt,res)
   i=minloc(res%residual_history(:,res%iterations+1),dim=1)
   print '(a,2f14.8)', 'best parameters: ',res%theta(i,:)
   print '(a,es14.6)', 'SSR: ',res%residual_history(i,res%iterations+1)
contains
   subroutine model(x,y,ierr)
      real(dp),intent(in)::x(:); real(dp),intent(out)::y(:); integer,intent(out)::ierr
      y=x(1)+x(2)*tt; ierr=0
   end subroutine model
end program basic_cgnm
