! GA generational engines translated from GA/R/ga.R and GA/R/gade.R.
! License: GPL-2.0-or-later.
module ga_core
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use ga_kinds, only : dp
  use ga_random, only : ga_seed, runif, randint, sample_without_replacement, shuffle_int
  use ga_utils, only : fitness_summary, garun, argsort_desc, argsort_asc, ga_pmutation
  use ga_operators
  implicit none
  private

  type, public :: ga_control_type
    integer :: pop_size = 50
    integer :: max_iter = 100
    integer :: run = -1
    real(dp) :: max_fitness = huge(1.0_dp)
    real(dp) :: pcrossover = 0.8_dp
    real(dp) :: pmutation = 0.1_dp
    integer :: elitism = -1
    integer :: selection = 0
    integer :: crossover = 0
    integer :: mutation = 0
    integer :: tournament_k = 3
    real(dp) :: selection_q = -1.0_dp
    real(dp) :: selection_r = -1.0_dp
    real(dp) :: blx_alpha = 0.5_dp
    real(dp) :: laplace_a = 0.0_dp
    real(dp) :: laplace_b = 0.15_dp
    real(dp) :: power_mutation = 10.0_dp
    logical :: adaptive_mutation = .false.
    real(dp) :: pmutation_initial = 0.5_dp
    real(dp) :: pmutation_final = 0.01_dp
    real(dp) :: pmutation_time = -1.0_dp
    logical :: keep_history = .true.
    integer :: seed = 0
  end type ga_control_type

  type, public :: ga_real_result
    integer :: iter = 0
    integer :: run = 0
    integer :: evaluations = 0
    real(dp) :: fitness_value = -huge(1.0_dp)
    real(dp), allocatable :: solution(:)
    real(dp), allocatable :: population(:,:)
    real(dp), allocatable :: fitness(:)
    real(dp), allocatable :: summary(:,:)
    real(dp), allocatable :: best_history(:,:)
  end type ga_real_result

  type, public :: ga_int_result
    integer :: iter = 0
    integer :: run = 0
    integer :: evaluations = 0
    real(dp) :: fitness_value = -huge(1.0_dp)
    integer, allocatable :: solution(:)
    integer, allocatable :: population(:,:)
    real(dp), allocatable :: fitness(:)
    real(dp), allocatable :: summary(:,:)
    integer, allocatable :: best_history(:,:)
  end type ga_int_result

  type, public :: de_control_type
    integer :: pop_size = -1
    integer :: max_iter = 100
    integer :: run = -1
    real(dp) :: max_fitness = huge(1.0_dp)
    real(dp) :: stepsize = 0.8_dp
    logical :: dither = .false.
    real(dp) :: pcrossover = 0.5_dp
    integer :: seed = 0
    logical :: keep_history = .true.
  end type de_control_type

  abstract interface
    function real_fitness_fn(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function real_fitness_fn
    function int_fitness_fn(x) result(f)
      import dp
      integer, intent(in) :: x(:)
      real(dp) :: f
    end function int_fitness_fn
  end interface

  public :: ga_real, ga_binary, ga_permutation, de_real
  public :: real_fitness_fn, int_fitness_fn
contains

  subroutine ga_real(fitness,lower,upper,control,result,suggestions)
    procedure(real_fitness_fn) :: fitness
    real(dp), intent(in) :: lower(:), upper(:)
    type(ga_control_type), intent(in), optional :: control
    type(ga_real_result), intent(out) :: result
    real(dp), intent(in), optional :: suggestions(:,:)
    type(ga_control_type) :: c
    real(dp), allocatable :: pop(:,:), newpop(:,:), fit(:), newfit(:)
    real(dp), allocatable :: old_sorted(:,:), elite_pop(:,:), child1(:), child2(:), mutant(:)
    real(dp), allocatable :: elite_fit(:), power(:), la(:), lb(:)
    integer, allocatable :: sel(:), ord(:), mating(:), worst(:)
    integer :: n,ps,iter,i,j,e,nm,p1,p2,ng,selstr,crossstr,mutstr,runlim
    real(dp) :: nanv,pm,tv,fsum(6)

    c=ga_control_type(); if(present(control)) c=control
    n=size(lower); if(size(upper)/=n) error stop "ga_real: bounds dimensions"
    ps=c%pop_size; if(ps<2) error stop "ga_real: population size must be >= 2"
    if(any(upper<lower)) error stop "ga_real: upper must be >= lower"
    if(c%seed/=0) call ga_seed(c%seed)
    e=c%elitism; if(e<0)e=max(1,nint(0.05_dp*real(ps,dp))); e=min(max(e,0),ps)
    runlim=c%run; if(runlim<0)runlim=c%max_iter
    selstr=c%selection; if(selstr==0)selstr=SEL_REAL_LINEAR_SCALE
    crossstr=c%crossover; if(crossstr==0)crossstr=CROSS_REAL_LOCAL
    mutstr=c%mutation; if(mutstr==0)mutstr=MUT_REAL_RANDOM
    nanv=ieee_value(1.0_dp,ieee_quiet_nan)

    allocate(pop(ps,n),newpop(ps,n),fit(ps),newfit(ps),sel(ps),ord(ps),worst(ps))
    allocate(old_sorted(ps,n),child1(n),child2(n),mutant(n),mating(2*(ps/2)))
    allocate(elite_pop(max(1,e),n),elite_fit(max(1,e)))
    allocate(power(n),la(n),lb(n)); power=c%power_mutation; la=c%laplace_a; lb=c%laplace_b
    call random_real_population(pop,lower,upper)
    if(present(suggestions)) then
      if(size(suggestions,2)/=n) error stop "ga_real: suggestion dimensions"
      ng=min(ps,size(suggestions,1)); if(ng>0)pop(1:ng,:)=suggestions(1:ng,:)
    end if
    fit=nanv
    allocate(result%summary(c%max_iter,6)); result%summary=nanv
    if(c%keep_history) then
      allocate(result%best_history(c%max_iter,n)); result%best_history=nanv
    end if

    do iter=1,c%max_iter
      do i=1,ps
        if(.not.ieee_is_finite(fit(i))) then
          fit(i)=fitness(pop(i,:)); result%evaluations=result%evaluations+1
        end if
      end do
      call fitness_summary(fit,fsum)
      result%summary(iter,:)=fsum
      call argsort_desc(fit,ord)
      if(c%keep_history) result%best_history(iter,:)=pop(ord(1),:)
      result%iter=iter
      if(iter>1) result%run=garun(result%summary(1:iter,1))
      if(result%run>=runlim .or. maxval(fit)>=c%max_fitness .or. iter==c%max_iter) exit

      old_sorted=pop(ord,:)
      if(e>0) then
        elite_pop(1:e,:)=old_sorted(1:e,:)
        elite_fit(1:e)=fit(ord(1:e))
      end if

      call selection_dispatch(fit,selstr,sel,c)
      newpop=pop(sel,:); newfit=fit(sel)
      pop=newpop; fit=newfit

      nm=ps/2; mating=[(i,i=1,2*nm)]; call shuffle_int(mating)
      if(c%pcrossover>0.0_dp) then
        do j=1,nm
          if(runif()<c%pcrossover) then
            p1=mating(2*j-1); p2=mating(2*j)
            call crossover_real(pop(p1,:),pop(p2,:),lower,upper,crossstr,child1,child2, &
                                c%blx_alpha,la,lb)
            pop(p1,:)=child1; pop(p2,:)=child2; fit(p1)=nanv; fit(p2)=nanv
          end if
        end do
      end if

      pm=c%pmutation
      if(c%adaptive_mutation) then
        tv=c%pmutation_time
        if(tv<=0.0_dp)tv=real(nint(real(c%max_iter,dp)/2.0_dp),dp)
        pm=ga_pmutation(iter,c%max_iter,c%pmutation_initial,c%pmutation_final,tv)
      end if
      if(pm>0.0_dp) then
        do i=1,ps
          if(runif()<pm) then
            call mutation_real(pop(i,:),lower,upper,mutstr,iter,c%max_iter,mutant,power)
            pop(i,:)=mutant; fit(i)=nanv
          end if
        end do
      end if

      if(e>0) then
        call fitness_order_with_nan(fit,.false.,worst)
        do i=1,e
          pop(worst(i),:)=elite_pop(i,:); fit(worst(i))=elite_fit(i)
        end do
      end if
    end do
    call finalize_real(pop,fit,result)
  end subroutine ga_real

  subroutine ga_binary(fitness,n_bits,control,result,suggestions)
    procedure(int_fitness_fn) :: fitness
    integer,intent(in)::n_bits
    type(ga_control_type),intent(in),optional::control
    type(ga_int_result),intent(out)::result
    integer,intent(in),optional::suggestions(:,:)
    type(ga_control_type)::c
    integer,allocatable::pop(:,:),newpop(:,:),elite_pop(:,:),c1(:),c2(:),mut(:)
    integer,allocatable::sel(:),ord(:),worst(:),mating(:)
    real(dp),allocatable::fit(:),newfit(:),elite_fit(:)
    integer::ps,e,runlim,selstr,crossstr,iter,i,j,nm,p1,p2,ng
    real(dp)::nanv,pm,tv,fsum(6)
    c=ga_control_type();if(present(control))c=control
    ps=c%pop_size;if(c%seed/=0)call ga_seed(c%seed);e=c%elitism
    if(e<0)e=max(1,nint(0.05_dp*real(ps,dp)));e=min(max(e,0),ps)
    runlim=c%run;if(runlim<0)runlim=c%max_iter
    selstr=c%selection;if(selstr==0)selstr=SEL_LINEAR_RANK
    crossstr=c%crossover;if(crossstr==0)crossstr=CROSS_SINGLE_POINT
    nanv=ieee_value(1.0_dp,ieee_quiet_nan)
    allocate(pop(ps,n_bits),newpop(ps,n_bits),c1(n_bits),c2(n_bits),mut(n_bits))
    allocate(elite_pop(max(1,e),n_bits),elite_fit(max(1,e)))
    allocate(fit(ps),newfit(ps),sel(ps),ord(ps),worst(ps),mating(2*(ps/2)))
    call random_binary_population(pop)
    if(present(suggestions))then
      if(size(suggestions,2)/=n_bits)error stop "ga_binary: suggestion dimensions"
      ng=min(ps,size(suggestions,1));if(ng>0)pop(1:ng,:)=suggestions(1:ng,:)
    end if
    fit=nanv;allocate(result%summary(c%max_iter,6));result%summary=nanv
    if(c%keep_history)then;allocate(result%best_history(c%max_iter,n_bits));result%best_history=0;end if
    do iter=1,c%max_iter
      do i=1,ps
      if(.not.ieee_is_finite(fit(i)))then
      fit(i)=fitness(pop(i,:))
      result%evaluations=result%evaluations+1
      end if
      end do
      call fitness_summary(fit,fsum)
      result%summary(iter,:)=fsum
      call argsort_desc(fit,ord)
      if(c%keep_history)result%best_history(iter,:)=pop(ord(1),:);result%iter=iter
      if(iter>1)result%run=garun(result%summary(1:iter,1))
      if(result%run>=runlim.or.maxval(fit)>=c%max_fitness.or.iter==c%max_iter)exit
      if(e>0)then
        elite_pop(1:e,:)=pop(ord(1:e),:)
        elite_fit(1:e)=fit(ord(1:e))
      end if
      call selection_dispatch(fit,selstr,sel,c);newpop=pop(sel,:);newfit=fit(sel);pop=newpop;fit=newfit
      nm=ps/2;mating=[(i,i=1,2*nm)];call shuffle_int(mating)
      do j=1,nm
      if(runif()<c%pcrossover)then
      p1=mating(2*j-1)
      p2=mating(2*j)
      call crossover_int(pop(p1,:),pop(p2,:),crossstr,c1,c2)
      pop(p1,:)=c1
      pop(p2,:)=c2
      fit(p1)=nanv
      fit(p2)=nanv
      end if
      end do
      pm=c%pmutation
      if(c%adaptive_mutation)then
      tv=c%pmutation_time
      if(tv<=0)tv=real(nint(real(c%max_iter,dp)/2),dp)
      pm=ga_pmutation(iter,c%max_iter,c%pmutation_initial,c%pmutation_final,tv)
      end if
      do i=1,ps;if(runif()<pm)then;call mutation_binary(pop(i,:),mut);pop(i,:)=mut;fit(i)=nanv;end if;end do
      if(e>0)then
      call fitness_order_with_nan(fit,.false.,worst)
      do i=1,e
      pop(worst(i),:)=elite_pop(i,:)
      fit(worst(i))=elite_fit(i)
      end do
      end if
    end do
    call finalize_int(pop,fit,result)
  end subroutine ga_binary

  subroutine ga_permutation(fitness,lower,upper,control,result,suggestions)
    procedure(int_fitness_fn)::fitness
    integer,intent(in)::lower,upper
    type(ga_control_type),intent(in),optional::control
    type(ga_int_result),intent(out)::result
    integer,intent(in),optional::suggestions(:,:)
    type(ga_control_type)::c
    integer,allocatable::pop(:,:),newpop(:,:),elite_pop(:,:),c1(:),c2(:),mut(:)
    integer,allocatable::sel(:),ord(:),worst(:),mating(:)
    real(dp),allocatable::fit(:),newfit(:),elite_fit(:)
    integer::n,ps,e,runlim,selstr,crossstr,mutstr,iter,i,j,nm,p1,p2,ng
    real(dp)::nanv,pm,tv,fsum(6)
    n=upper-lower+1;if(n<2)error stop "ga_permutation: invalid range"
    c=ga_control_type();if(present(control))c=control
    ps=c%pop_size;if(c%seed/=0)call ga_seed(c%seed);e=c%elitism
    if(e<0)e=max(1,nint(0.05_dp*real(ps,dp)));e=min(max(e,0),ps)
    runlim=c%run;if(runlim<0)runlim=c%max_iter
    selstr=c%selection;if(selstr==0)selstr=SEL_LINEAR_RANK
    crossstr=c%crossover;if(crossstr==0)crossstr=CROSS_PERM_OX
    mutstr=c%mutation;if(mutstr==0)mutstr=MUT_PERM_INVERSION
    nanv=ieee_value(1.0_dp,ieee_quiet_nan)
    allocate(pop(ps,n),newpop(ps,n),c1(n),c2(n),mut(n))
    allocate(elite_pop(max(1,e),n),elite_fit(max(1,e)))
    allocate(fit(ps),newfit(ps),sel(ps),ord(ps),worst(ps),mating(2*(ps/2)))
    call random_perm_population(pop,lower)
    if(present(suggestions))then
    if(size(suggestions,2)/=n)error stop "ga_permutation: suggestion dimensions"
    ng=min(ps,size(suggestions,1))
    if(ng>0)pop(1:ng,:)=suggestions(1:ng,:)
    end if
    fit=nanv;allocate(result%summary(c%max_iter,6));result%summary=nanv
    if(c%keep_history)then;allocate(result%best_history(c%max_iter,n));result%best_history=0;end if
    do iter=1,c%max_iter
      do i=1,ps
      if(.not.ieee_is_finite(fit(i)))then
      fit(i)=fitness(pop(i,:))
      result%evaluations=result%evaluations+1
      end if
      end do
      call fitness_summary(fit,fsum)
      result%summary(iter,:)=fsum
      call argsort_desc(fit,ord)
      if(c%keep_history)result%best_history(iter,:)=pop(ord(1),:)
      result%iter=iter
      if(iter>1)result%run=garun(result%summary(1:iter,1))
      if(result%run>=runlim.or.maxval(fit)>=c%max_fitness.or.iter==c%max_iter)exit
      if(e>0)then
        elite_pop(1:e,:)=pop(ord(1:e),:)
        elite_fit(1:e)=fit(ord(1:e))
      end if
      call selection_dispatch(fit,selstr,sel,c);newpop=pop(sel,:);newfit=fit(sel);pop=newpop;fit=newfit
      nm=ps/2;mating=[(i,i=1,2*nm)];call shuffle_int(mating)
      do j=1,nm
      if(runif()<c%pcrossover)then
      p1=mating(2*j-1)
      p2=mating(2*j)
      call crossover_int(pop(p1,:),pop(p2,:),crossstr,c1,c2)
      pop(p1,:)=c1
      pop(p2,:)=c2
      fit(p1)=nanv
      fit(p2)=nanv
      end if
      end do
      pm=c%pmutation
      if(c%adaptive_mutation)then
      tv=c%pmutation_time
      if(tv<=0)tv=real(nint(real(c%max_iter,dp)/2),dp)
      pm=ga_pmutation(iter,c%max_iter,c%pmutation_initial,c%pmutation_final,tv)
      end if
      do i=1,ps;if(runif()<pm)then;call mutation_permutation(pop(i,:),mutstr,mut);pop(i,:)=mut;fit(i)=nanv;end if;end do
      if(e>0)then
      call fitness_order_with_nan(fit,.false.,worst)
      do i=1,e
      pop(worst(i),:)=elite_pop(i,:)
      fit(worst(i))=elite_fit(i)
      end do
      end if
    end do
    call finalize_int(pop,fit,result)
  end subroutine ga_permutation

  subroutine de_real(fitness,lower,upper,control,result,suggestions)
    procedure(real_fitness_fn)::fitness
    real(dp),intent(in)::lower(:),upper(:)
    type(de_control_type),intent(in),optional::control
    type(ga_real_result),intent(out)::result
    real(dp),intent(in),optional::suggestions(:,:)
    type(de_control_type)::c
    real(dp),allocatable::pop(:,:),fit(:),trial(:),v(:)
    integer::n,ps,iter,i,j,jforce,r(3),ng,runlim,best
    real(dp)::fi,ft,fsum(6)
    c=de_control_type()
    if(present(control))c=control
    n=size(lower)
    if(size(upper)/=n)error stop "de_real: bounds dimensions"
    ps=c%pop_size
    if(ps<0)ps=10*n
    if(ps<4)error stop "de_real: population must be >= 4"
    if(c%seed/=0)call ga_seed(c%seed)
    runlim=c%run;if(runlim<0)runlim=c%max_iter
    allocate(pop(ps,n),fit(ps),trial(n),v(n));call random_real_population(pop,lower,upper)
    if(present(suggestions))then
    if(size(suggestions,2)/=n)error stop "de_real: suggestions"
    ng=min(ps,size(suggestions,1))
    if(ng>0)pop(1:ng,:)=suggestions(1:ng,:)
    end if
    do i=1,ps;fit(i)=fitness(pop(i,:));result%evaluations=result%evaluations+1;end do
    allocate(result%summary(c%max_iter,6));result%summary=0
    if(c%keep_history)allocate(result%best_history(c%max_iter,n))
    do iter=1,c%max_iter
      call fitness_summary(fit,fsum)
      result%summary(iter,:)=fsum
      best=maxloc(fit,dim=1)
      if(c%keep_history)result%best_history(iter,:)=pop(best,:)
      result%iter=iter
      if(iter>1)result%run=garun(result%summary(1:iter,1))
      if(result%run>=runlim.or.maxval(fit)>=c%max_fitness.or.iter==c%max_iter)exit
      do i=1,ps
        call sample_without_replacement(ps,3,r);fi=c%stepsize;if(c%dither)fi=runif(0.5_dp,1.0_dp)
        v=pop(r(1),:)+fi*(pop(r(2),:)-pop(r(3),:));trial=pop(i,:);jforce=randint(1,n)
        do j=1,n
          if(runif()<clamp01(c%pcrossover).or.j==jforce)trial(j)=v(j)
          if(trial(j)<lower(j))trial(j)=lower(j)+runif()*(upper(j)-lower(j))
          if(trial(j)>upper(j))trial(j)=upper(j)-runif()*(upper(j)-lower(j))
        end do
        ft=fitness(trial);result%evaluations=result%evaluations+1;if(ft>fit(i))then;fit(i)=ft;pop(i,:)=trial;end if
      end do
    end do
    call finalize_real(pop,fit,result)
  contains
    pure real(dp) function clamp01(x);real(dp),intent(in)::x;clamp01=min(max(x,0.0_dp),1.0_dp);end function
  end subroutine de_real

  subroutine selection_dispatch(fit,strategy,sel,c)
    real(dp),intent(in)::fit(:)
    integer,intent(in)::strategy
    integer,intent(out)::sel(:)
    type(ga_control_type),intent(in)::c
    if(c%selection_q>=0.0_dp.and.c%selection_r>=0.0_dp)then
      call select_indices(fit,strategy,sel,c%selection_q,c%selection_r,c%tournament_k)
    else if(c%selection_q>=0.0_dp)then
      call select_indices(fit,strategy,sel,q=c%selection_q,tournament_k=c%tournament_k)
    else if(c%selection_r>=0.0_dp)then
      call select_indices(fit,strategy,sel,r=c%selection_r,tournament_k=c%tournament_k)
    else
      call select_indices(fit,strategy,sel,tournament_k=c%tournament_k)
    end if
  end subroutine selection_dispatch

  subroutine fitness_order_with_nan(fit,decreasing,idx)
    real(dp),intent(in)::fit(:)
    logical,intent(in)::decreasing
    integer,intent(out)::idx(size(fit))
    real(dp),allocatable::x(:)
    integer::i
    allocate(x(size(fit)));x=fit
    do i=1,size(x);if(.not.ieee_is_finite(x(i)))x(i)=merge(-huge(1.0_dp),huge(1.0_dp),decreasing);end do
    if(decreasing)then;call argsort_desc(x,idx);else;call argsort_asc(x,idx);end if
  end subroutine fitness_order_with_nan

  subroutine finalize_real(pop,fit,result)
    real(dp),intent(in)::pop(:,:),fit(:);type(ga_real_result),intent(inout)::result;integer::i
    i=maxloc(fit,dim=1);result%fitness_value=fit(i);allocate(result%solution(size(pop,2)));result%solution=pop(i,:)
    allocate(result%population(size(pop,1),size(pop,2)),result%fitness(size(fit)))
    result%population=pop
    result%fitness=fit
  end subroutine finalize_real

  subroutine finalize_int(pop,fit,result)
    integer,intent(in)::pop(:,:);real(dp),intent(in)::fit(:);type(ga_int_result),intent(inout)::result;integer::i
    i=maxloc(fit,dim=1);result%fitness_value=fit(i);allocate(result%solution(size(pop,2)));result%solution=pop(i,:)
    allocate(result%population(size(pop,1),size(pop,2)),result%fitness(size(fit)))
    result%population=pop
    result%fitness=fit
  end subroutine finalize_int
end module ga_core
