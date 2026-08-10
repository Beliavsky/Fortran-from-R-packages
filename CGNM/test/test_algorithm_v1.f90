program test_algorithm_v1
   use cgnm, only : dp,cgnm_problem,cgnm_options,cgnm_result,cgnm_init_problem,cgnm_fit
   implicit none
   type(cgnm_problem)::prob
   type(cgnm_options)::opt
   type(cgnm_result)::res
   real(dp)::target(2),lo(2),hi(2)
   integer::ierr,i
   target=[3._dp,-1._dp]; lo=-4._dp; hi=4._dp
   call cgnm_init_problem(prob,model,target,lo,hi,ierr=ierr)
   opt%num_minimizers=30; opt%num_iterations=12; opt%seed=101
   call cgnm_fit(prob,opt,res,algorithm_version=1)
   if(res%status/=0) error stop trim(res%message)
   i=minloc(res%residual_history(:,res%iterations+1),dim=1)
   if(maxval(abs(res%theta(i,:)-[1._dp,2._dp]))>1.e-5_dp) error stop 'algorithm v1'
contains
   subroutine model(x,y,ierr)
      real(dp),intent(in)::x(:); real(dp),intent(out)::y(:); integer,intent(out)::ierr
      y=[x(1)+x(2),x(1)-x(2)]; ierr=0
   end subroutine model
end program test_algorithm_v1
