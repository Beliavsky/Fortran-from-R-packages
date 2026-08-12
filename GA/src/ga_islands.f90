! Coarse-grained island GA translated from GA/R/gaIslands.R.
! This implementation is serial; independent islands can be parallelized by callers.
! License: GPL-2.0-or-later.
module ga_islands
  use ga_kinds, only : dp
  use ga_random, only : ga_seed, shuffle_int
  use ga_utils, only : argsort_desc, garun
  use ga_operators, only : random_real_population, random_binary_population, random_perm_population
  use ga_core, only : ga_control_type, ga_real_result, ga_int_result
  use ga_core, only : ga_real, ga_binary, ga_permutation, real_fitness_fn, int_fitness_fn
  implicit none
  private

  type, public :: island_control_type
    integer :: pop_size = 100
    integer :: num_islands = 4
    real(dp) :: migration_rate = 0.10_dp
    integer :: migration_interval = 10
    integer :: max_iter = 1000
    integer :: run = -1
    integer :: seed = 0
    type(ga_control_type) :: ga
  end type island_control_type

  type, public :: island_real_result
    integer :: epochs = 0
    integer :: iter = 0
    real(dp) :: fitness_value = -huge(1.0_dp)
    real(dp), allocatable :: solution(:)
    real(dp), allocatable :: island_fitness(:)
  end type island_real_result

  type, public :: island_int_result
    integer :: epochs = 0
    integer :: iter = 0
    real(dp) :: fitness_value = -huge(1.0_dp)
    integer, allocatable :: solution(:)
    real(dp), allocatable :: island_fitness(:)
  end type island_int_result

  public :: ga_islands_real, ga_islands_binary, ga_islands_permutation
contains
  subroutine ga_islands_real(fitness,lower,upper,control,result,suggestions)
    procedure(real_fitness_fn)::fitness
    real(dp),intent(in)::lower(:),upper(:)
    type(island_control_type),intent(in),optional::control
    type(island_real_result),intent(out)::result
    real(dp),intent(in),optional::suggestions(:,:)
    type(island_control_type)::ic
    type(ga_control_type)::gc
    type(ga_real_result),allocatable::r(:)
    real(dp),allocatable::pops(:,:,:),nextpops(:,:,:),fit(:,:),migrants(:,:),best_hist(:,:)
    integer,allocatable::ord(:),eligible(:)
    integer::ni,isize,n,mig,nepoch,e,ep,k,to,j,ng,bestk,besti,runlim,istart,istop
    ic=island_control_type();if(present(control))ic=control
    ni=ic%num_islands;n=size(lower);isize=max(10,ic%pop_size/max(1,ni))
    mig=max(1,int(ic%migration_rate*real(isize,dp)))
    nepoch=max(1,ic%max_iter/max(1,ic%migration_interval))
    runlim=ic%run;if(runlim<0)runlim=ic%max_iter
    if(ic%seed/=0)call ga_seed(ic%seed)
    allocate(pops(ni,isize,n),nextpops(ni,isize,n),fit(ni,isize),r(ni),ord(isize),eligible(isize))
    allocate(best_hist(ni,nepoch*ic%migration_interval));best_hist=-huge(1.0_dp)
    do k=1,ni
      call random_real_population(pops(k,:,:),lower,upper)
      if(present(suggestions))then
        ng=min(isize,size(suggestions,1));if(ng>0)pops(k,1:ng,:)=suggestions(1:ng,:)
      end if
    end do
    gc=ic%ga;gc%pop_size=isize;gc%max_iter=ic%migration_interval
    gc%run=ic%migration_interval;gc%seed=0;gc%keep_history=.false.
    e=gc%elitism;if(e<0)e=max(1,nint(0.05_dp*real(isize,dp)));e=min(e,isize-1)
    do ep=1,nepoch
      istart=(ep-1)*ic%migration_interval+1
      istop=ep*ic%migration_interval
      do k=1,ni
        call ga_real(fitness,lower,upper,gc,r(k),pops(k,:,:))
        pops(k,:,:)=r(k)%population;fit(k,:)=r(k)%fitness
        best_hist(k,istart:istart+r(k)%iter-1)=r(k)%summary(1:r(k)%iter,1)
        if(r(k)%iter<ic%migration_interval) &
          best_hist(k,istart+r(k)%iter:istop)=r(k)%fitness_value
      end do
      nextpops=pops
      do k=1,ni
        to=mod(k,ni)+1;call argsort_desc(fit(k,:),ord);allocate(migrants(mig,n));migrants=pops(k,ord(1:mig),:)
        call argsort_desc(fit(to,:),ord);eligible(1:isize-e)=ord(e+1:isize);call shuffle_int(eligible(1:isize-e))
        do j=1,mig;nextpops(to,eligible(j),:)=migrants(j,:);end do;deallocate(migrants)
      end do
      result%epochs=ep;result%iter=ep*ic%migration_interval
      bestk=1;besti=maxloc(fit(1,:),dim=1)
      do k=2,ni;j=maxloc(fit(k,:),dim=1);if(fit(k,j)>fit(bestk,besti))then;bestk=k;besti=j;end if;end do
      if(all_run_converged(best_hist,istop,runlim))exit
      if(all(best_hist(:,1:istop)>=ic%ga%max_fitness))exit
      if(ep<nepoch)pops=nextpops
    end do
    allocate(result%island_fitness(ni));do k=1,ni;result%island_fitness(k)=maxval(fit(k,:));end do
    bestk=maxloc(result%island_fitness,dim=1);besti=maxloc(fit(bestk,:),dim=1)
    result%fitness_value=fit(bestk,besti);allocate(result%solution(n));result%solution=pops(bestk,besti,:)
  end subroutine ga_islands_real

  subroutine ga_islands_binary(fitness,n_bits,control,result,suggestions)
    procedure(int_fitness_fn)::fitness
    integer,intent(in)::n_bits
    type(island_control_type),intent(in),optional::control
    type(island_int_result),intent(out)::result
    integer,intent(in),optional::suggestions(:,:)
    type(island_control_type)::ic;type(ga_control_type)::gc;type(ga_int_result),allocatable::r(:)
    integer,allocatable::pops(:,:,:),nextpops(:,:,:),ord(:),eligible(:),migrants(:,:)
    real(dp),allocatable::fit(:,:),best_hist(:,:)
    integer::ni,isize,mig,nepoch,e,ep,k,to,j,ng,bestk,besti,runlim,istart,istop
    ic=island_control_type();if(present(control))ic=control;ni=ic%num_islands;isize=max(10,ic%pop_size/max(1,ni))
    mig=max(1,int(ic%migration_rate*real(isize,dp)))
    nepoch=max(1,ic%max_iter/max(1,ic%migration_interval))
    runlim=ic%run;if(runlim<0)runlim=ic%max_iter
    if(ic%seed/=0)call ga_seed(ic%seed)
    allocate(pops(ni,isize,n_bits),nextpops(ni,isize,n_bits),fit(ni,isize),r(ni),ord(isize),eligible(isize))
    allocate(best_hist(ni,nepoch*ic%migration_interval));best_hist=-huge(1.0_dp)
    do k=1,ni
    call random_binary_population(pops(k,:,:))
    if(present(suggestions))then
    ng=min(isize,size(suggestions,1))
    if(ng>0)pops(k,1:ng,:)=suggestions(1:ng,:)
    end if
    end do
    gc=ic%ga
    gc%pop_size=isize
    gc%max_iter=ic%migration_interval
    gc%run=ic%migration_interval
    gc%seed=0
    gc%keep_history=.false.
    e=gc%elitism;if(e<0)e=max(1,nint(0.05_dp*real(isize,dp)));e=min(e,isize-1)
    do ep=1,nepoch
      istart=(ep-1)*ic%migration_interval+1
      istop=ep*ic%migration_interval
      do k=1,ni
      call ga_binary(fitness,n_bits,gc,r(k),pops(k,:,:))
      pops(k,:,:)=r(k)%population
      fit(k,:)=r(k)%fitness
      best_hist(k,istart:istart+r(k)%iter-1)=r(k)%summary(1:r(k)%iter,1)
      if(r(k)%iter<ic%migration_interval) &
        best_hist(k,istart+r(k)%iter:istop)=r(k)%fitness_value
      end do
      nextpops=pops
      do k=1,ni
      to=mod(k,ni)+1
      call argsort_desc(fit(k,:),ord)
      allocate(migrants(mig,n_bits))
      migrants=pops(k,ord(1:mig),:)
      call argsort_desc(fit(to,:),ord)
      eligible(1:isize-e)=ord(e+1:isize)
      call shuffle_int(eligible(1:isize-e))
      do j=1,mig
      nextpops(to,eligible(j),:)=migrants(j,:)
      end do
      deallocate(migrants)
      end do
      result%epochs=ep;result%iter=ep*ic%migration_interval
      if(all_run_converged(best_hist,istop,runlim))exit
      if(all(best_hist(:,1:istop)>=ic%ga%max_fitness))exit
      if(ep<nepoch)pops=nextpops
    end do
    allocate(result%island_fitness(ni));do k=1,ni;result%island_fitness(k)=maxval(fit(k,:));end do
    bestk=maxloc(result%island_fitness,dim=1)
    besti=maxloc(fit(bestk,:),dim=1)
    result%fitness_value=fit(bestk,besti)
    allocate(result%solution(n_bits))
    result%solution=pops(bestk,besti,:)
  end subroutine ga_islands_binary

  subroutine ga_islands_permutation(fitness,lower,upper,control,result,suggestions)
    procedure(int_fitness_fn)::fitness
    integer,intent(in)::lower,upper
    type(island_control_type),intent(in),optional::control
    type(island_int_result),intent(out)::result
    integer,intent(in),optional::suggestions(:,:)
    type(island_control_type)::ic;type(ga_control_type)::gc;type(ga_int_result),allocatable::r(:)
    integer,allocatable::pops(:,:,:),nextpops(:,:,:),ord(:),eligible(:),migrants(:,:)
    real(dp),allocatable::fit(:,:),best_hist(:,:)
    integer::n,ni,isize,mig,nepoch,e,ep,k,to,j,ng,bestk,besti,runlim,istart,istop
    n=upper-lower+1
    ic=island_control_type()
    if(present(control))ic=control
    ni=ic%num_islands
    isize=max(10,ic%pop_size/max(1,ni))
    mig=max(1,int(ic%migration_rate*real(isize,dp)))
    nepoch=max(1,ic%max_iter/max(1,ic%migration_interval))
    runlim=ic%run;if(runlim<0)runlim=ic%max_iter
    if(ic%seed/=0)call ga_seed(ic%seed)
    allocate(pops(ni,isize,n),nextpops(ni,isize,n),fit(ni,isize),r(ni),ord(isize),eligible(isize))
    allocate(best_hist(ni,nepoch*ic%migration_interval));best_hist=-huge(1.0_dp)
    do k=1,ni
    call random_perm_population(pops(k,:,:),lower)
    if(present(suggestions))then
    ng=min(isize,size(suggestions,1))
    if(ng>0)pops(k,1:ng,:)=suggestions(1:ng,:)
    end if
    end do
    gc=ic%ga
    gc%pop_size=isize
    gc%max_iter=ic%migration_interval
    gc%run=ic%migration_interval
    gc%seed=0
    gc%keep_history=.false.
    e=gc%elitism;if(e<0)e=max(1,nint(0.05_dp*real(isize,dp)));e=min(e,isize-1)
    do ep=1,nepoch
      istart=(ep-1)*ic%migration_interval+1
      istop=ep*ic%migration_interval
      do k=1,ni
      call ga_permutation(fitness,lower,upper,gc,r(k),pops(k,:,:))
      pops(k,:,:)=r(k)%population
      fit(k,:)=r(k)%fitness
      best_hist(k,istart:istart+r(k)%iter-1)=r(k)%summary(1:r(k)%iter,1)
      if(r(k)%iter<ic%migration_interval) &
        best_hist(k,istart+r(k)%iter:istop)=r(k)%fitness_value
      end do
      nextpops=pops
      do k=1,ni
      to=mod(k,ni)+1
      call argsort_desc(fit(k,:),ord)
      allocate(migrants(mig,n))
      migrants=pops(k,ord(1:mig),:)
      call argsort_desc(fit(to,:),ord)
      eligible(1:isize-e)=ord(e+1:isize)
      call shuffle_int(eligible(1:isize-e))
      do j=1,mig
      nextpops(to,eligible(j),:)=migrants(j,:)
      end do
      deallocate(migrants)
      end do
      result%epochs=ep;result%iter=ep*ic%migration_interval
      if(all_run_converged(best_hist,istop,runlim))exit
      if(all(best_hist(:,1:istop)>=ic%ga%max_fitness))exit
      if(ep<nepoch)pops=nextpops
    end do
    allocate(result%island_fitness(ni));do k=1,ni;result%island_fitness(k)=maxval(fit(k,:));end do
    bestk=maxloc(result%island_fitness,dim=1)
    besti=maxloc(fit(bestk,:),dim=1)
    result%fitness_value=fit(bestk,besti)
    allocate(result%solution(n))
    result%solution=pops(bestk,besti,:)
  end subroutine ga_islands_permutation
  logical function all_run_converged(hist,nused,runlim) result(done)
    real(dp),intent(in)::hist(:,:)
    integer,intent(in)::nused,runlim
    integer::k
    done=.true.
    do k=1,size(hist,1)
      if(garun(hist(k,1:nused))<=runlim)then
        done=.false.
        return
      end if
    end do
  end function all_run_converged

end module ga_islands
