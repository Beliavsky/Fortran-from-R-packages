! SPDX-License-Identifier: GPL-3.0-only
program test_optimization
   use nmof
   implicit none
   integer :: failures
   type(optimization_result) :: res, best
   type(binary_optimization_result) :: bres
   real(dp), allocatable :: vals(:)
   real(dp) :: levels(5,2)
   integer :: counts(2)

   failures=0
   call de_opt(sphere,[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],res,n_population=40,n_generations=150, &
      minmax_constraint=.true.,seed=12345_i8)
   call check_true('DE status',res%status==nmof_ok)
   call check_true('DE sphere',res%ofvalue<1.0e-7_dp)

   call ps_opt(sphere,[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],res,n_population=50,n_generations=180, &
      inertia=0.7_dp,c1=1.4_dp,c2=1.4_dp,max_velocity=1.0_dp,minmax_constraint=.true.,seed=23456_i8)
   call check_true('PS status',res%status==nmof_ok)
   call check_true('PS sphere',res%ofvalue<1.0e-6_dp)

   call ga_opt(binary_cost,12,bres,n_population=40,n_generations=120,mutation_probability=0.03_dp,seed=34567_i8)
   call check_true('GA status',bres%status==nmof_ok)
   call check_true('GA optimum',count(bres%xbest)>=11)

   call local_search(sphere,toward_zero,[4.0_dp,-3.0_dp],res,n_steps=80,seed=45678_i8)
   call check_true('LS improves',res%ofvalue<1.0e-8_dp)

   call simulated_annealing(sphere,toward_zero,[4.0_dp,-3.0_dp],res,n_temperatures=4,steps_per_temperature=20, &
      initial_temperature=1.0_dp,seed=56789_i8)
   call check_true('SA improves',res%ofvalue<1.0e-8_dp)

   call threshold_accepting(sphere,toward_zero,[4.0_dp,-3.0_dp],res,n_thresholds=4,steps_per_threshold=20, &
      thresholds=[0.2_dp,0.1_dp,0.01_dp,0.0_dp],seed=67890_i8)
   call check_true('TA improves',res%ofvalue<1.0e-8_dp)

   call greedy_search(sphere,axis_neighbours,[1.0_dp,-1.0_dp],res,max_iterations=20)
   call check_true('greedy optimum',res%ofvalue<1.0e-12_dp)

   levels=0.0_dp
   levels(:,1)=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp]
   levels(:,2)=levels(:,1); counts=[5,5]
   call grid_search(shifted_sphere,levels,counts,res)
   call check_close('grid x1',res%xbest(1),1.0_dp,0.0_dp)
   call check_close('grid x2',res%xbest(2),-1.0_dp,0.0_dp)

   call restart_opt(simple_optimizer,3,best,vals)
   call check_close('restart best',best%ofvalue,1.0_dp,0.0_dp)
   call check_true('restart values',size(vals)==3)

   if(failures>0) then
      write(*,'(a,i0)') 'test_optimization failures: ',failures
      error stop 1
   end if
   write(*,'(a)') 'test_optimization: PASS'
contains
   function sphere(x,context) result(f)
      real(dp),intent(in)::x(:)
      class(*),intent(in),optional::context
      real(dp)::f
      f=dot_product(x,x)
   end function sphere
   function shifted_sphere(x,context) result(f)
      real(dp),intent(in)::x(:)
      class(*),intent(in),optional::context
      real(dp)::f
      f=(x(1)-1.0_dp)**2+(x(2)+1.0_dp)**2
   end function shifted_sphere
   function binary_cost(x,context) result(f)
      logical,intent(in)::x(:)
      class(*),intent(in),optional::context
      real(dp)::f
      f=-real(count(x),dp)
   end function binary_cost
   subroutine toward_zero(x,xn,rng,iteration,total_iterations,context)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::xn(:)
      type(rng_state),intent(inout)::rng
      integer,intent(in)::iteration,total_iterations
      class(*),intent(in),optional::context
      xn=0.5_dp*x
   end subroutine toward_zero
   subroutine axis_neighbours(x,neighbours,context)
      real(dp),intent(in)::x(:)
      real(dp),allocatable,intent(out)::neighbours(:,:)
      class(*),intent(in),optional::context
      integer::i,j
      allocate(neighbours(size(x),2*size(x))); neighbours=spread(x,2,2*size(x)); j=0
      do i=1,size(x)
         j=j+1; neighbours(i,j)=x(i)-0.5_dp
         j=j+1; neighbours(i,j)=x(i)+0.5_dp
      end do
   end subroutine axis_neighbours
   subroutine simple_optimizer(result,context)
      type(optimization_result),intent(out)::result
      class(*),intent(in),optional::context
      integer,save::counter=0
      counter=counter+1
      allocate(result%xbest(1)); result%xbest=real(counter,dp); result%ofvalue=real(4-counter,dp); result%status=nmof_ok
   end subroutine simple_optimizer
   subroutine check_close(name,actual,expected,tol)
      character(len=*),intent(in)::name
      real(dp),intent(in)::actual,expected,tol
      if(abs(actual-expected)>tol) then
         failures=failures+1
         write(*,'(a,2(1x,es16.8))') trim(name)//' failed:',actual,expected
      end if
   end subroutine check_close
   subroutine check_true(name,condition)
      character(len=*),intent(in)::name
      logical,intent(in)::condition
      if(.not.condition) then; failures=failures+1; write(*,'(a)') trim(name)//' failed'; end if
   end subroutine check_true
end program test_optimization
