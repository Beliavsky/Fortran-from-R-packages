program multiple_minima
   use cgnm, only : dp,cgnm_problem,cgnm_options,cgnm_result,cgnm_init_problem,cgnm_fit
   implicit none
   type(cgnm_problem)::prob
   type(cgnm_options)::opt
   type(cgnm_result)::res
   real(dp)::target(1),lo(1),hi(1)
   integer::ierr,i
   target=1._dp; lo=-2._dp; hi=2._dp
   call cgnm_init_problem(prob,model,target,lo,hi,ierr=ierr)
   opt%num_minimizers=60; opt%num_iterations=15; opt%seed=15
   call cgnm_fit(prob,opt,res)
   do i=1,min(10,size(res%theta,1))
      print '(i3,2x,f12.7,2x,es12.4)',i,res%theta(i,1), &
            res%residual_history(i,res%iterations+1)
   end do
contains
   subroutine model(x,y,ierr)
      real(dp),intent(in)::x(:); real(dp),intent(out)::y(:); integer,intent(out)::ierr
      y(1)=x(1)*x(1); ierr=0
   end subroutine model
end program multiple_minima
