! SPDX-License-Identifier: GPL-3.0-only
module ao
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ao_kinds, only : dp
  use ao_types
  use ao_random, only : ao_seed, generate_random_partition
  use ao_history_mod
  use ao_base_optimizer, only : solve_block
  implicit none
  private
  public :: dp, ao_options, ao_block, ao_real_vector, ao_result, ao_multi_result
  public :: AO_PARTITION_SEQUENTIAL, AO_PARTITION_RANDOM, AO_PARTITION_NONE, AO_PARTITION_CUSTOM
  public :: AO_BASE_BFGS, AO_BASE_NELDER_MEAD, AO_BASE_NEWTON
  public :: ao_optimize, ao_optimize_multiple, ao_seed, generate_random_partition
  public :: split_estimate, euclidean_parameter_norm
contains
  function objective_dispatch(objective,x) result(value)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value=objective(x)
  end function objective_dispatch

  function norm_dispatch(norm_fn,x,y) result(value)
    procedure(ao_norm_fn) :: norm_fn
    real(dp), intent(in) :: x(:),y(:)
    real(dp) :: value
    value=norm_fn(x,y)
  end function norm_dispatch

  function euclidean_parameter_norm(x,y) result(value)
    real(dp), intent(in) :: x(:),y(:)
    real(dp) :: value
    value=sqrt(sum((x-y)**2))
  end function euclidean_parameter_norm

  subroutine split_estimate(estimate,npar,parts,ok)
    real(dp), intent(in) :: estimate(:)
    integer, intent(in) :: npar(:)
    type(ao_real_vector), allocatable, intent(out) :: parts(:)
    logical, intent(out), optional :: ok
    integer :: i,s,e
    logical :: good
    good=all(npar>=0) .and. sum(npar)==size(estimate)
    if(.not.good) then
      allocate(parts(0)); if(present(ok)) ok=.false.; return
    end if
    allocate(parts(size(npar))); s=1
    do i=1,size(npar)
      e=s+npar(i)-1; allocate(parts(i)%value(npar(i)))
      if(npar(i)>0) parts(i)%value=estimate(s:e)
      s=e+1
    end do
    if(present(ok)) ok=.true.
  end subroutine split_estimate

  subroutine build_partition(options,npar,blocks,ok)
    type(ao_options), intent(in) :: options
    integer, intent(in) :: npar
    type(ao_block), allocatable, intent(out) :: blocks(:)
    logical, intent(out) :: ok
    integer :: i
    ok=.true.
    select case(options%partition)
    case(AO_PARTITION_SEQUENTIAL)
      allocate(blocks(npar))
      do i=1,npar; allocate(blocks(i)%index(1)); blocks(i)%index=i; end do
    case(AO_PARTITION_RANDOM)
      call generate_random_partition(npar,options%new_block_probability, &
           options%minimum_block_number,blocks)
    case(AO_PARTITION_NONE)
      allocate(blocks(1)); allocate(blocks(1)%index(npar)); blocks(1)%index=[(i,i=1,npar)]
    case(AO_PARTITION_CUSTOM)
      if(.not.allocated(options%custom_partition)) then; allocate(blocks(0)); ok=.false.; return; end if
      blocks=options%custom_partition
      do i=1,size(blocks)
        if(.not.allocated(blocks(i)%index)) then; ok=.false.; return; end if
        if(size(blocks(i)%index)==0) then; ok=.false.; return; end if
        if(any(blocks(i)%index<1) .or. any(blocks(i)%index>npar)) then; ok=.false.; return; end if
      end do
    case default
      allocate(blocks(0)); ok=.false.
    end select
  end subroutine build_partition

  subroutine fill_bounds(n,lower,upper,lb,ub,ok)
    integer,intent(in)::n
    real(dp),intent(in),optional::lower(:),upper(:)
    real(dp),allocatable,intent(out)::lb(:),ub(:)
    logical,intent(out)::ok
    allocate(lb(n),ub(n)); lb=-huge(1.0_dp); ub=huge(1.0_dp); ok=.true.
    if(present(lower)) then
      if(size(lower)==1) then; lb=lower(1); else if(size(lower)==n) then; lb=lower; else; ok=.false.; return; end if
    end if
    if(present(upper)) then
      if(size(upper)==1) then; ub=upper(1); else if(size(upper)==n) then; ub=upper; else; ok=.false.; return; end if
    end if
    if(any(lb>ub)) ok=.false.
  end subroutine fill_bounds

  subroutine finalize(result,history,iteration,minimize,reason,converged)
    type(ao_result), intent(out) :: result
    type(ao_history), intent(in) :: history
    integer,intent(in)::iteration
    logical,intent(in)::minimize,converged
    character(len=*),intent(in)::reason
    integer::ib
    ib=history_best_index(history,minimize)
    result%estimate=history%parameter(:,ib); result%value=history%value(ib)
    result%seconds=history_total_seconds(history); result%iterations=iteration
    result%stopping_reason=reason; result%converged=converged; result%details=history
  end subroutine finalize

  subroutine ao_optimize(objective,initial,result,options,lower,upper,gradient,hessian,parameter_norm)
    procedure(ao_objective_fn) :: objective
    real(dp),intent(in)::initial(:)
    type(ao_result),intent(out)::result
    type(ao_options),intent(in),optional::options
    real(dp),intent(in),optional::lower(:),upper(:)
    procedure(ao_gradient_fn),optional::gradient
    procedure(ao_hessian_fn),optional::hessian
    procedure(ao_norm_fn),optional::parameter_norm
    type(ao_options)::opt
    type(ao_history)::history
    type(ao_block),allocatable::blocks(:)
    real(dp),allocatable::theta(:),lb(:),ub(:),pblock(:)
    real(dp)::value,newvalue,cpu0,cpu1,elapsed,change
    integer::iteration,j,past_idx
    logical::ok,block_ok,conv
    character(len=160)::reason

    opt=ao_options(); if(present(options)) opt=options
    if(size(initial)==0) then
      allocate(result%estimate(0)); result%stopping_reason='empty initial parameter'; return
    end if
    call fill_bounds(size(initial),lower,upper,lb,ub,ok)
    if(.not.ok .or. any(initial<lb) .or. any(initial>ub)) then
      result%estimate=initial; result%value=huge(1.0_dp); result%stopping_reason='invalid bounds or initial value'; return
    end if
    if(opt%new_block_probability<0.0_dp .or. opt%new_block_probability>1.0_dp .or. &
       opt%minimum_block_number<1 .or. opt%minimum_block_number>size(initial) .or. &
       opt%tolerance_history<1 .or. opt%base_max_iterations<1) then
      result%estimate=initial; result%stopping_reason='invalid options'; return
    end if
    theta=initial; value=objective_dispatch(objective,theta)
    if(.not.ieee_is_finite(value)) then
      result%estimate=theta; result%value=value; result%stopping_reason='objective invalid at initial value'; return
    end if
    call history_initialize(history,theta,value)
    iteration=0; conv=.false.; reason='not terminated yet'

    do
      if(iteration>=opt%iteration_limit) then
        write(reason,'(a,i0,a)') 'iteration limit of ',opt%iteration_limit,' reached'; exit
      end if
      elapsed=history_total_seconds(history)
      if(elapsed>=opt%seconds_limit) then
        reason='time limit reached'; exit
      end if
      if(iteration>=opt%tolerance_history) then
        past_idx=history_first_index(history,iteration-opt%tolerance_history)
        if(past_idx>0) then
          change=abs(history%value(history%n)-history%value(past_idx))
          if(opt%tolerance_value>0.0_dp .and. change<opt%tolerance_value) then
            reason='function-value tolerance reached'; conv=.true.; exit
          end if
          if(present(parameter_norm)) then
            change=norm_dispatch(parameter_norm,history%parameter(:,history%n),history%parameter(:,past_idx))
          else
            change=euclidean_parameter_norm(history%parameter(:,history%n),history%parameter(:,past_idx))
          end if
          if(opt%tolerance_parameter>0.0_dp .and. change<opt%tolerance_parameter) then
            reason='parameter tolerance reached'; conv=.true.; exit
          end if
        end if
      end if

      iteration=iteration+1
      call build_partition(opt,size(theta),blocks,ok)
      if(.not.ok .or. size(blocks)==0) then; reason='invalid parameter partition'; exit; end if
      do j=1,size(blocks)
        allocate(pblock(size(blocks(j)%index)))
        call cpu_time(cpu0)
        if(present(gradient) .and. present(hessian)) then
          call solve_block(objective,theta,blocks(j)%index,lb,ub,opt%minimize,opt%base_optimizer, &
               opt%base_max_iterations,opt%base_gradient_tolerance,opt%finite_difference_step, &
               pblock,newvalue,block_ok,gradient,hessian)
        else if(present(gradient)) then
          call solve_block(objective,theta,blocks(j)%index,lb,ub,opt%minimize,opt%base_optimizer, &
               opt%base_max_iterations,opt%base_gradient_tolerance,opt%finite_difference_step, &
               pblock,newvalue,block_ok,gradient=gradient)
        else if(present(hessian)) then
          call solve_block(objective,theta,blocks(j)%index,lb,ub,opt%minimize,opt%base_optimizer, &
               opt%base_max_iterations,opt%base_gradient_tolerance,opt%finite_difference_step, &
               pblock,newvalue,block_ok,hessian=hessian)
        else
          call solve_block(objective,theta,blocks(j)%index,lb,ub,opt%minimize,opt%base_optimizer, &
               opt%base_max_iterations,opt%base_gradient_tolerance,opt%finite_difference_step, &
               pblock,newvalue,block_ok)
        end if
        call cpu_time(cpu1)
        if(.not.block_ok .or. .not.ieee_is_finite(newvalue)) then
          reason='solving a parameter block failed'; deallocate(pblock); exit
        end if
        theta(blocks(j)%index)=pblock; value=newvalue
        call history_append(history,iteration,value,theta,blocks(j)%index,max(0.0_dp,cpu1-cpu0))
        if(opt%verbose) then
          write(*,'(a,i0,a,i0,a,es16.8)') 'iteration ',iteration,', block ',j,', value ',value
        end if
        deallocate(pblock)
      end do
      if(trim(reason)/='not terminated yet') exit
    end do
    call finalize(result,history,iteration,opt%minimize,reason,conv)
  end subroutine ao_optimize

  subroutine ao_optimize_multiple(objective,initials,partition_modes,base_methods,result,options, &
       lower,upper,gradient,hessian,parameter_norm)
    procedure(ao_objective_fn) :: objective
    real(dp),intent(in)::initials(:,:)
    integer,intent(in)::partition_modes(:),base_methods(:)
    type(ao_multi_result),intent(out)::result
    type(ao_options),intent(in),optional::options
    real(dp),intent(in),optional::lower(:),upper(:)
    procedure(ao_gradient_fn),optional::gradient
    procedure(ao_hessian_fn),optional::hessian
    procedure(ao_norm_fn),optional::parameter_norm
    type(ao_options)::opt,one
    integer::i,j,k,q,best,np
    real(dp)::bestval
    opt=ao_options(); if(present(options)) opt=options
    np=size(initials,2)*size(partition_modes)*size(base_methods)
    allocate(result%process(np)); q=0
    do k=1,size(base_methods)
      do j=1,size(partition_modes)
        do i=1,size(initials,2)
          q=q+1; one=opt; one%partition=partition_modes(j); one%base_optimizer=base_methods(k)
          if(present(gradient) .and. present(hessian) .and. present(parameter_norm)) then
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper, &
                 gradient,hessian,parameter_norm)
          else if(present(gradient) .and. present(hessian)) then
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper, &
                 gradient,hessian)
          else if(present(gradient) .and. present(parameter_norm)) then
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper, &
                 gradient=gradient,parameter_norm=parameter_norm)
          else if(present(hessian) .and. present(parameter_norm)) then
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper, &
                 hessian=hessian,parameter_norm=parameter_norm)
          else if(present(gradient)) then
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper, &
                 gradient=gradient)
          else if(present(hessian)) then
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper, &
                 hessian=hessian)
          else if(present(parameter_norm)) then
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper, &
                 parameter_norm=parameter_norm)
          else
            call ao_optimize(objective,initials(:,i),result%process(q),one,lower,upper)
          end if
        end do
      end do
    end do
    best=1; bestval=result%process(1)%value
    do q=2,np
      if((opt%minimize .and. result%process(q)%value<bestval) .or. &
         (.not.opt%minimize .and. result%process(q)%value>bestval)) then
        best=q; bestval=result%process(q)%value
      end if
    end do
    result%best_process=best; result%estimate=result%process(best)%estimate; result%value=bestval
    result%seconds=0.0_dp
    do q=1,np; result%seconds=result%seconds+result%process(q)%seconds; end do
  end subroutine ao_optimize_multiple
end module ao
