program test_bounds_mo
   use cgnm, only : dp,cgnm_problem,cgnm_options,cgnm_result,cgnm_init_problem,cgnm_fit
   implicit none
   type(cgnm_problem) :: prob
   type(cgnm_options) :: opt
   type(cgnm_result) :: res
   real(dp) :: target(1),lo(1),hi(1),lb(1),ub(1),mow(1),mov(1)
   integer :: ierr,i
   target=4.0_dp; lo=0.2_dp; hi=4.8_dp; lb=0.0_dp; ub=5.0_dp
   mow=0.05_dp; mov=2.0_dp
   call cgnm_init_problem(prob,model,target,lo,hi,lower_bound=lb,upper_bound=ub, &
                          mo_weights=mow,mo_values=mov,ierr=ierr)
   if(ierr/=0) error stop 'init bounds'
   opt%num_minimizers=24; opt%num_iterations=15; opt%seed=32
   opt%initial_lambda=0.1_dp; opt%stay_in_initial_range=.true.
   call cgnm_fit(prob,opt,res)
   if(res%status/=0) error stop trim(res%message)
   i=minloc(res%residual_history(:,res%iterations+1),dim=1)
   if(abs(res%theta(i,1)-2.0_dp)>2.0e-3_dp) then
      print *,res%theta(i,1),res%residual_history(i,res%iterations+1)
      error stop 'bounded/MO fit'
   end if
   if(any(res%theta(:,1)<=0.0_dp) .or. any(res%theta(:,1)>=5.0_dp)) error stop 'bounds'
contains
   subroutine model(x,y,ierr)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      integer,intent(out)::ierr
      y(1)=x(1)*x(1); ierr=0
   end subroutine model
end program test_bounds_mo
