program test_extensions
   use cgnm, only : dp,cgnm_problem,cgnm_options,cgnm_result,cgnm_init_problem,cgnm_fit, &
                    cgnm_bootstrap,cgnm_ebe,top_indices
   implicit none
   type(cgnm_problem) :: prob
   type(cgnm_options) :: opt
   type(cgnm_result) :: base,boot,ebe
   real(dp) :: target(4),lo(2),hi(2),t(4)
   integer, allocatable :: idx(:)
   integer :: groups(4),ierr
   t=[0._dp,1._dp,2._dp,3._dp]; target=1.5_dp+0.7_dp*t
   lo=[-2._dp,-2._dp]; hi=[4._dp,3._dp]
   call cgnm_init_problem(prob,model,target,lo,hi,ierr=ierr)
   if(ierr/=0) error stop 'init ext'
   opt%num_minimizers=30; opt%num_iterations=10; opt%seed=91
   call cgnm_fit(prob,opt,base)
   if(base%status/=0) error stop 'base'
   call cgnm_bootstrap(prob,opt,base,16,1,boot,seed=19)
   if(boot%status/=0 .or. size(boot%x,1)/=16) error stop 'bootstrap'
   groups=[1,1,2,2]
   call cgnm_ebe(prob,opt,base,groups,2,1.0_dp,ebe,seed=23)
   if(ebe%status/=0 .or. size(ebe%x,1)/=4) error stop 'ebe'
   call top_indices(base,3,idx)
   if(size(idx)/=3) error stop 'top indices'
contains
   subroutine model(x,y,ierr)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      integer,intent(out)::ierr
      y=x(1)+x(2)*t; ierr=0
   end subroutine model
end program test_extensions
