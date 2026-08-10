program test_linear
   use cgnm, only : dp,cgnm_problem,cgnm_options,cgnm_result,cgnm_init_problem,cgnm_fit
   implicit none
   type(cgnm_problem) :: prob
   type(cgnm_options) :: opt
   type(cgnm_result) :: res
   real(dp) :: target(2), lo(2), hi(2), best(2), err
   integer :: ierr,i
   target=[3.0_dp,-1.0_dp] ! x=(1,2)
   lo=[-4.0_dp,-4.0_dp]; hi=[4.0_dp,4.0_dp]
   call cgnm_init_problem(prob,model,target,lo,hi,ierr=ierr)
   if(ierr/=0) error stop 'init'
   opt%num_minimizers=36; opt%num_iterations=12; opt%seed=17
   opt%initial_lambda=1.0_dp; opt%gamma=2.0_dp
   call cgnm_fit(prob,opt,res,algorithm_version=3)
   if(res%status/=0) error stop trim(res%message)
   i=minloc(res%residual_history(:,res%iterations+1),dim=1)
   best=res%theta(i,:); err=maxval(abs(best-[1.0_dp,2.0_dp]))
   if(err>1.0e-5_dp) then
      print *,best,res%residual_history(i,res%iterations+1)
      error stop 'linear convergence'
   end if
contains
   subroutine model(x,y,ierr)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      integer,intent(out)::ierr
      y(1)=x(1)+x(2); y(2)=x(1)-x(2); ierr=0
   end subroutine model
end program test_linear
