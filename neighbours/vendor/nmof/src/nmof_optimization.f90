! SPDX-License-Identifier: GPL-3.0-only
module nmof_optimization
   use nmof_kinds, only: dp, i8
   use nmof_rng, only: rng_state, rng_seed, rng_uniform, rng_normal, rng_integer, rng_shuffle, rng_logical
   use nmof_math, only: quantile_type7
   use nmof_types, only: optimization_result, binary_optimization_result, nmof_ok, nmof_invalid_input
   implicit none
   private
   public :: de_opt, ps_opt, ga_opt, local_search, simulated_annealing
   public :: threshold_accepting, greedy_search, grid_search, restart_opt
   public :: real_objective, real_neighbour, real_repair, real_velocity_change
   public :: binary_objective, binary_repair, all_neighbours, optimizer_callback

   abstract interface
      function real_objective(x, context) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         class(*), intent(in), optional :: context
         real(dp) :: f
      end function real_objective
      subroutine real_neighbour(x, xn, rng, iteration, total_iterations, context)
         import dp, rng_state
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: xn(:)
         type(rng_state), intent(inout) :: rng
         integer, intent(in) :: iteration, total_iterations
         class(*), intent(in), optional :: context
      end subroutine real_neighbour
      subroutine real_repair(x, context)
         import dp
         real(dp), intent(inout) :: x(:)
         class(*), intent(in), optional :: context
      end subroutine real_repair
      subroutine real_velocity_change(v, context)
         import dp
         real(dp), intent(inout) :: v(:)
         class(*), intent(in), optional :: context
      end subroutine real_velocity_change
      function binary_objective(x, context) result(f)
         import dp
         logical, intent(in) :: x(:)
         class(*), intent(in), optional :: context
         real(dp) :: f
      end function binary_objective
      subroutine binary_repair(x, context)
         logical, intent(inout) :: x(:)
         class(*), intent(in), optional :: context
      end subroutine binary_repair
      subroutine all_neighbours(x, neighbours, context)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), allocatable, intent(out) :: neighbours(:, :)
         class(*), intent(in), optional :: context
      end subroutine all_neighbours
      subroutine optimizer_callback(result, context)
         import optimization_result
         type(optimization_result), intent(out) :: result
         class(*), intent(in), optional :: context
      end subroutine optimizer_callback
   end interface
contains
   subroutine call_objective(fun, x, f, context)
      procedure(real_objective) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
      class(*), intent(in), optional :: context
      if (present(context)) then
         f = fun(x, context)
      else
         f = fun(x)
      end if
   end subroutine call_objective

   subroutine call_penalty(fun, x, f, context)
      procedure(real_objective) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
      class(*), intent(in), optional :: context
      call call_objective(fun, x, f, context)
   end subroutine call_penalty

   subroutine call_repair(fun, x, context)
      procedure(real_repair) :: fun
      real(dp), intent(inout) :: x(:)
      class(*), intent(in), optional :: context
      if (present(context)) then
         call fun(x, context)
      else
         call fun(x)
      end if
   end subroutine call_repair

   subroutine call_neighbour(fun, x, xn, rng, iteration, total_iterations, context)
      procedure(real_neighbour) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: xn(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(in) :: iteration, total_iterations
      class(*), intent(in), optional :: context
      if (present(context)) then
         call fun(x, xn, rng, iteration, total_iterations, context)
      else
         call fun(x, xn, rng, iteration, total_iterations)
      end if
   end subroutine call_neighbour

   subroutine evaluate_population(objective, population, values, context, penalty)
      procedure(real_objective) :: objective
      real(dp), intent(in) :: population(:, :)
      real(dp), intent(out) :: values(:)
      class(*), intent(in), optional :: context
      procedure(real_objective), optional :: penalty
      real(dp) :: p
      integer :: j
      do j = 1, size(population, 2)
         call call_objective(objective, population(:, j), values(j), context)
         if (present(penalty)) then
            call call_penalty(penalty, population(:, j), p, context)
            values(j) = values(j) + p
         end if
      end do
   end subroutine evaluate_population

   subroutine bounds_repair(x, lower, upper)
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: lower(:), upper(:)
      x = min(max(x, lower), upper)
   end subroutine bounds_repair

   subroutine de_opt(objective, lower, upper, result, n_population, n_generations, &
                     differential_weight, crossover_probability, seed, initial_population, &
                     minmax_constraint, context, penalty, repair)
      procedure(real_objective) :: objective
      real(dp), intent(in) :: lower(:), upper(:)
      type(optimization_result), intent(out) :: result
      integer, intent(in), optional :: n_population, n_generations
      real(dp), intent(in), optional :: differential_weight(:), crossover_probability
      integer(i8), intent(in), optional :: seed
      real(dp), intent(in), optional :: initial_population(:, :)
      logical, intent(in), optional :: minmax_constraint
      class(*), intent(in), optional :: context
      procedure(real_objective), optional :: penalty
      procedure(real_repair), optional :: repair
      type(rng_state) :: rng
      integer(i8) :: seedv
      real(dp), allocatable :: pop(:, :), trial(:, :), fval(:), ftrial(:), fw(:)
      real(dp) :: cr
      integer, allocatable :: perm(:), r1(:), r2(:), r3(:)
      integer :: d, np, ng, g, j, k, best
      logical :: mm
      d = size(lower)
      np = 50; if (present(n_population)) np = n_population
      ng = 300; if (present(n_generations)) ng = n_generations
      cr = 0.9_dp; if (present(crossover_probability)) cr = crossover_probability
      mm = .false.; if (present(minmax_constraint)) mm = minmax_constraint
      if (d < 1 .or. size(upper) /= d .or. any(lower > upper) .or. np < 4 .or. ng < 1 .or. cr < 0 .or. cr > 1) then
         result%status = nmof_invalid_input
         return
      end if
      seedv=104729_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(pop(d,np), trial(d,np), fval(np), ftrial(np), perm(np), r1(np), r2(np), r3(np), fw(d))
      if (present(initial_population)) then
         if (any(shape(initial_population) /= [d,np])) then
            result%status = nmof_invalid_input; return
         end if
         pop = initial_population
      else
         do j=1,np; do k=1,d
            pop(k,j)=lower(k)+rng_uniform(rng)*(upper(k)-lower(k))
         end do; end do
      end if
      if (mm) then; do j=1,np; call bounds_repair(pop(:,j),lower,upper); end do; end if
      if (present(repair)) then; do j=1,np; call call_repair(repair,pop(:,j),context); end do; end if
      call evaluate_population(objective,pop,fval,context,penalty)
      allocate(result%history(ng,np)); result%history = 0.0_dp
      do g=1,ng
         perm=[(j,j=1,np)]; call rng_shuffle(rng,perm)
         r1(1)=perm(np); r1(2:np)=perm(1:np-1)
         r2(1)=r1(np); r2(2:np)=r1(1:np-1)
         r3(1)=r2(np); r3(2:np)=r2(1:np-1)
         do j=1,np
            if (present(differential_weight)) then
               if (size(differential_weight)==d) then
                  fw=differential_weight
               else
                  fw=differential_weight(1)
               end if
            else
               fw=0.5_dp
            end if
            trial(:,j)=pop(:,r1(j))+fw*(pop(:,r2(j))-pop(:,r3(j)))
            do k=1,d
               if (rng_uniform(rng)>cr) trial(k,j)=pop(k,j)
            end do
            if (mm) call bounds_repair(trial(:,j),lower,upper)
            if (present(repair)) call call_repair(repair,trial(:,j),context)
         end do
         call evaluate_population(objective,trial,ftrial,context,penalty)
         do j=1,np
            if (ftrial(j)<fval(j)) then; pop(:,j)=trial(:,j); fval(j)=ftrial(j); end if
         end do
         result%history(g,:)=fval
      end do
      best=minloc(fval,dim=1)
      allocate(result%xbest(d),result%population_values(np))
      result%xbest=pop(:,best); result%ofvalue=fval(best); result%population_values=fval
      result%iterations=ng; result%status=nmof_ok
   end subroutine de_opt

   subroutine ps_opt(objective, lower, upper, result, n_population, n_generations, c1, c2, &
                     inertia, initial_velocity_scale, max_velocity, seed, initial_population, &
                     minmax_constraint, context, penalty, repair, change_velocity)
      procedure(real_objective) :: objective
      real(dp), intent(in) :: lower(:), upper(:)
      type(optimization_result), intent(out) :: result
      integer, intent(in), optional :: n_population,n_generations
      real(dp), intent(in), optional :: c1,c2,inertia,initial_velocity_scale,max_velocity
      integer(i8), intent(in), optional :: seed
      real(dp), intent(in), optional :: initial_population(:, :)
      logical, intent(in), optional :: minmax_constraint
      class(*), intent(in), optional :: context
      procedure(real_objective), optional :: penalty
      procedure(real_repair), optional :: repair
      procedure(real_velocity_change), optional :: change_velocity
      type(rng_state) :: rng
      integer(i8) :: seedv
      real(dp), allocatable :: pop(:,:), vel(:,:), pbest(:,:), f(:), fbest(:)
      real(dp) :: a1,a2,iner,ivs,mv
      integer :: d,np,ng,g,j,k,best
      logical :: mm
      d=size(lower); np=100; ng=500
      if(present(n_population)) np=n_population; if(present(n_generations)) ng=n_generations
      a1=1.0_dp; a2=1.0_dp; iner=0.9_dp; ivs=1.0_dp; mv=1.0_dp
      if(present(c1)) a1=c1; if(present(c2)) a2=c2; if(present(inertia)) iner=inertia
      if(present(initial_velocity_scale)) ivs=initial_velocity_scale
      if(present(max_velocity)) mv=max_velocity
      mm=.false.; if(present(minmax_constraint)) mm=minmax_constraint
      if(d<1 .or. size(upper)/=d .or. any(lower>upper) .or. np<1 .or. ng<1) then
         result%status=nmof_invalid_input; return
      end if
      seedv=130363_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(pop(d,np),vel(d,np),pbest(d,np),f(np),fbest(np),result%history(ng,np))
      if(present(initial_population)) then
         if(any(shape(initial_population)/=[d,np])) then; result%status=nmof_invalid_input; return; end if
         pop=initial_population
      else
         do j=1,np; do k=1,d; pop(k,j)=lower(k)+rng_uniform(rng)*(upper(k)-lower(k)); end do; end do
      end if
      do j=1,np; do k=1,d; vel(k,j)=ivs*rng_normal(rng); end do; end do
      if(mm) then; do j=1,np; call bounds_repair(pop(:,j),lower,upper); end do; end if
      if(present(repair)) then; do j=1,np; call call_repair(repair,pop(:,j),context); end do; end if
      call evaluate_population(objective,pop,f,context,penalty)
      pbest=pop; fbest=f; best=minloc(fbest,dim=1)
      do g=1,ng
         do j=1,np
            do k=1,d
               vel(k,j)=iner*vel(k,j)+a1*rng_uniform(rng)*(pbest(k,j)-pop(k,j))+ &
                        a2*rng_uniform(rng)*(pbest(k,best)-pop(k,j))
               vel(k,j)=min(max(vel(k,j),-mv),mv)
            end do
            if(present(change_velocity)) then
               if(present(context)) then; call change_velocity(vel(:,j),context); else; call change_velocity(vel(:,j)); end if
            end if
            pop(:,j)=pop(:,j)+vel(:,j)
            if(mm) call bounds_repair(pop(:,j),lower,upper)
            if(present(repair)) call call_repair(repair,pop(:,j),context)
         end do
         call evaluate_population(objective,pop,f,context,penalty)
         do j=1,np
            if(f(j)<fbest(j)) then; pbest(:,j)=pop(:,j); fbest(j)=f(j); end if
         end do
         best=minloc(fbest,dim=1); result%history(g,:)=fbest
      end do
      allocate(result%xbest(d),result%population_values(np))
      result%xbest=pbest(:,best); result%ofvalue=fbest(best); result%population_values=fbest
      result%iterations=ng; result%status=nmof_ok
   end subroutine ps_opt

   subroutine ga_opt(objective, n_bits, result, n_population, n_generations, mutation_probability, &
                     crossover, seed, initial_population, context, penalty, repair)
      procedure(binary_objective) :: objective
      integer, intent(in) :: n_bits
      type(binary_optimization_result), intent(out) :: result
      integer, intent(in), optional :: n_population,n_generations
      real(dp), intent(in), optional :: mutation_probability
      character(len=*), intent(in), optional :: crossover
      integer(i8), intent(in), optional :: seed
      logical, intent(in), optional :: initial_population(:, :)
      class(*), intent(in), optional :: context
      procedure(binary_objective), optional :: penalty
      procedure(binary_repair), optional :: repair
      type(rng_state) :: rng
      integer(i8) :: seedv
      logical, allocatable :: pop(:,:), child(:,:)
      real(dp), allocatable :: f(:),fc(:)
      integer, allocatable :: order(:)
      real(dp) :: prob,pv
      integer :: np,ng,g,j,k,cut,best,next
      character(len=16) :: cross
      np=50; ng=300; prob=0.01_dp; cross='onepoint'
      if(present(n_population)) np=n_population; if(present(n_generations)) ng=n_generations
      if(present(mutation_probability)) prob=mutation_probability
      if(present(crossover)) cross=adjustl(crossover)
      if(n_bits<2 .or. np<2 .or. ng<1 .or. prob<0 .or. prob>1) then
         result%status=nmof_invalid_input; return
      end if
      seedv=15485863_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(pop(n_bits,np),child(n_bits,np),f(np),fc(np),order(np),result%history(ng,np))
      if(present(initial_population)) then
         if(any(shape(initial_population)/=[n_bits,np])) then; result%status=nmof_invalid_input; return; end if
         pop=initial_population
      else
         do j=1,np; do k=1,n_bits; pop(k,j)=rng_logical(rng,0.5_dp); end do; end do
      end if
      if(present(repair)) then
         do j=1,np
            if(present(context)) then; call repair(pop(:,j),context); else; call repair(pop(:,j)); end if
         end do
      end if
      call eval_binary(pop,f)
      do g=1,ng
         order=[(j,j=1,np)]; call rng_shuffle(rng,order)
         do j=1,np
            next=j+1; if(next>np) next=1
            child(:,j)=pop(:,order(j))
            if(index(lowercase(cross),'uniform')==1) then
               do k=1,n_bits; if(rng_uniform(rng)>0.5_dp) child(k,j)=pop(k,order(next)); end do
            else
               cut=1+rng_integer(rng,n_bits-1)
               child(cut:n_bits,j)=pop(cut:n_bits,order(next))
            end if
            do k=1,n_bits; if(rng_uniform(rng)<prob) child(k,j)=.not.child(k,j); end do
            if(present(repair)) then
               if(present(context)) then; call repair(child(:,j),context); else; call repair(child(:,j)); end if
            end if
         end do
         call eval_binary(child,fc)
         do j=1,np; if(fc(j)<f(j)) then; pop(:,j)=child(:,j); f(j)=fc(j); end if; end do
         result%history(g,:)=f
      end do
      best=minloc(f,dim=1); allocate(result%xbest(n_bits),result%population_values(np))
      result%xbest=pop(:,best); result%ofvalue=f(best); result%population_values=f
      result%iterations=ng; result%status=nmof_ok
   contains
      subroutine eval_binary(p,v)
         logical,intent(in)::p(:,:); real(dp),intent(out)::v(:)
         integer::ii
         do ii=1,size(p,2)
            if(present(context)) then; v(ii)=objective(p(:,ii),context); else; v(ii)=objective(p(:,ii)); end if
            if(present(penalty)) then
               if(present(context)) then; pv=penalty(p(:,ii),context); else; pv=penalty(p(:,ii)); end if
               v(ii)=v(ii)+pv
            end if
         end do
      end subroutine eval_binary
      pure function lowercase(s) result(t)
         character(len=*),intent(in)::s; character(len=len(s))::t; integer::ii,c
         do ii=1,len(s); c=iachar(s(ii:ii)); if(c>=65.and.c<=90) then; t(ii:ii)=achar(c+32); else; t(ii:ii)=s(ii:ii); end if; end do
      end function lowercase
   end subroutine ga_opt

   subroutine local_search(objective, neighbour, x0, result, n_steps, target, seed, context)
      procedure(real_objective) :: objective
      procedure(real_neighbour) :: neighbour
      real(dp), intent(in) :: x0(:)
      type(optimization_result), intent(out) :: result
      integer, intent(in), optional :: n_steps
      real(dp), intent(in), optional :: target
      integer(i8), intent(in), optional :: seed
      class(*), intent(in), optional :: context
      type(rng_state) :: rng
      integer(i8) :: seedv
      real(dp), allocatable :: xc(:),xn(:)
      real(dp)::fc,fn
      integer::ns,s
      ns=1000; if(present(n_steps)) ns=n_steps
      if(ns<1) then; result%status=nmof_invalid_input; return; end if
      seedv=32452843_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(xc(size(x0)),xn(size(x0)),result%history(ns,2)); xc=x0
      call call_objective(objective,xc,fc,context)
      do s=1,ns
         call call_neighbour(neighbour,xc,xn,rng,s,ns,context)
         call call_objective(objective,xn,fn,context)
         if(fn<=fc) then; xc=xn; fc=fn; end if
         result%history(s,:)=[fn,fc]
         if(present(target)) then; if(fc<=target) exit; end if
      end do
      allocate(result%xbest(size(x0)),result%population_values(1)); result%xbest=xc
      result%ofvalue=fc; result%population_values=[fc]; result%iterations=s; result%status=nmof_ok
   end subroutine local_search

   subroutine simulated_annealing(objective, neighbour, x0, result, n_temperatures, steps_per_temperature, &
                                  calibration_steps, initial_temperature, final_temperature, initial_probability, &
                                  alpha, step_multiplier, target, seed, context)
      procedure(real_objective)::objective
      procedure(real_neighbour)::neighbour
      real(dp),intent(in)::x0(:)
      type(optimization_result),intent(out)::result
      integer,intent(in),optional::n_temperatures,steps_per_temperature,calibration_steps
      real(dp),intent(in),optional::initial_temperature,final_temperature,initial_probability,alpha,step_multiplier,target
      integer(i8),intent(in),optional::seed
      class(*),intent(in),optional::context
      type(rng_state)::rng
      integer(i8)::seedv
      real(dp),allocatable::xc(:),xn(:),xbest(:),diffs(:)
      real(dp)::fc,fn,fbest,temp,ft,ip,ar,mult,u
      integer::nt,ns,nd,t,s,counter,total,i,ns_current
      nt=10; ns=1000; nd=2000; ft=0.0_dp; ip=0.4_dp; ar=0.9_dp; mult=1.0_dp
      if(present(n_temperatures)) nt=n_temperatures; if(present(steps_per_temperature)) ns=steps_per_temperature
      if(present(calibration_steps)) nd=calibration_steps; if(present(final_temperature)) ft=final_temperature
      if(present(initial_probability)) ip=initial_probability; if(present(alpha)) ar=alpha
      if(present(step_multiplier)) mult=step_multiplier
      if(nt<1.or.ns<1.or.nd<1.or.ip<=0.or.ip>=1) then; result%status=nmof_invalid_input; return; end if
      seedv=49979687_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(xc(size(x0)),xn(size(x0)),xbest(size(x0))); xc=x0
      if(present(initial_temperature)) then
         temp=initial_temperature
      else
         allocate(diffs(nd)); call call_objective(objective,xc,fc,context)
         do i=1,nd
            call call_neighbour(neighbour,xc,xn,rng,i,nd,context); call call_objective(objective,xn,fn,context)
            diffs(i)=abs(fc-fn); xc=xn; fc=fn
         end do
         temp=-sum(diffs)/real(nd,dp)/log(ip)
         if(temp<=0.0_dp) temp=1.0_dp
      end if
      total=nt*ns; allocate(result%history(total,2)); result%history=0.0_dp
      xc=x0; call call_objective(objective,xc,fc,context); xbest=xc; fbest=fc
      counter=0; ns_current=ns
      do t=1,nt
         do s=1,ns_current
            counter=counter+1; if(counter>size(result%history,1)) exit
            call call_neighbour(neighbour,xc,xn,rng,counter,total,context); call call_objective(objective,xn,fn,context)
            u=rng_uniform(rng)
            if(fn<=fc .or. exp((fc-fn)/max(temp,tiny(1.0_dp)))>u) then
               xc=xn; fc=fn; if(fn<=fbest) then; xbest=xn; fbest=fn; end if
            end if
            result%history(counter,:)=[fn,fc]
            if(present(target)) then; if(fbest<=target) exit; end if
         end do
         if(present(target)) then; if(fbest<=target) exit; end if
         if(temp<=ft) exit
         ns_current=max(1,nint(mult*real(ns_current,dp))); temp=temp*ar
      end do
      allocate(result%xbest(size(x0)),result%population_values(1)); result%xbest=xbest
      result%ofvalue=fbest; result%population_values=[fbest]; result%iterations=counter; result%status=nmof_ok
   end subroutine simulated_annealing

   subroutine threshold_accepting(objective, neighbour, x0, result, n_thresholds, steps_per_threshold, &
                                  calibration_steps, quantile_scale, thresholds, scale, drop_zero, step_up, &
                                  target, seed, context)
      procedure(real_objective)::objective
      procedure(real_neighbour)::neighbour
      real(dp),intent(in)::x0(:)
      type(optimization_result),intent(out)::result
      integer,intent(in),optional::n_thresholds,steps_per_threshold,calibration_steps,step_up
      real(dp),intent(in),optional::quantile_scale,thresholds(:),scale,target
      logical,intent(in),optional::drop_zero
      integer(i8),intent(in),optional::seed
      class(*),intent(in),optional::context
      type(rng_state)::rng
      integer(i8)::seedv
      real(dp),allocatable::vt(:),base(:),diffs(:),xc(:),xn(:),xbest(:)
      real(dp)::fc,fn,fbest,q,sc,p
      integer::nt,ns,nd,su,t,s,counter,total,i,nkeep,k
      logical::dz
      nt=10; ns=1000; nd=2000; su=0; q=0.5_dp; sc=1.0_dp; dz=.false.
      if(present(n_thresholds)) nt=n_thresholds; if(present(steps_per_threshold)) ns=steps_per_threshold
      if(present(calibration_steps)) nd=calibration_steps; if(present(step_up)) su=step_up
      if(present(quantile_scale)) q=quantile_scale; if(present(scale)) sc=max(0.0_dp,scale)
      if(present(drop_zero)) dz=drop_zero
      if(nt<1.or.ns<1.or.nd<1.or.su<0) then; result%status=nmof_invalid_input; return; end if
      seedv=67867967_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(xc(size(x0)),xn(size(x0)),xbest(size(x0))); xc=x0
      if(present(thresholds)) then
         allocate(base(size(thresholds))); base=thresholds
      else if(q<sqrt(epsilon(1.0_dp))) then
         allocate(base(nt)); base=0.0_dp
      else
         allocate(diffs(nd)); call call_objective(objective,xc,fc,context)
         do i=1,nd
            call call_neighbour(neighbour,xc,xn,rng,i,nd,context); call call_objective(objective,xn,fn,context)
            diffs(i)=abs(fc-fn); xc=xn; fc=fn
         end do
         if(dz) then
            nkeep=count(abs(diffs)>tiny(1.0_dp)); allocate(base(nt))
            if(nkeep==0) then; base=0.0_dp
            else
               do i=1,nt
                  p=q*real(nt-i,dp)/real(nt,dp)
                  base(i)=quantile_type7(pack(diffs,abs(diffs)>tiny(1.0_dp)),p)
               end do
            end if
         else
            allocate(base(nt)); do i=1,nt
               p=q*real(nt-i,dp)/real(nt,dp); base(i)=quantile_type7(diffs,p)
            end do
         end if
         base(nt)=0.0_dp
      end if
      allocate(vt(size(base)*(su+1))); do k=0,su; vt(k*size(base)+1:(k+1)*size(base))=base*sc; end do
      total=ns*size(vt); allocate(result%history(total,2)); result%history=0.0_dp
      xc=x0; call call_objective(objective,xc,fc,context); xbest=xc; fbest=fc; counter=0
      do t=1,size(vt)
         do s=1,ns
            counter=counter+1; call call_neighbour(neighbour,xc,xn,rng,counter,total,context)
            call call_objective(objective,xn,fn,context)
            if(fn<=fc+vt(t)) then; xc=xn; fc=fn; if(fn<=fbest) then; xbest=xn; fbest=fn; end if; end if
            result%history(counter,:)=[fn,fc]
            if(present(target)) then; if(fbest<=target) exit; end if
         end do
         if(present(target)) then; if(fbest<=target) exit; end if
      end do
      allocate(result%xbest(size(x0)),result%population_values(size(vt))); result%xbest=xbest
      result%ofvalue=fbest; result%population_values=vt; result%iterations=counter; result%status=nmof_ok
   end subroutine threshold_accepting

   subroutine greedy_search(objective, neighbours_fun, x0, result, max_iterations, context)
      procedure(real_objective)::objective
      procedure(all_neighbours)::neighbours_fun
      real(dp),intent(in)::x0(:)
      type(optimization_result),intent(out)::result
      integer,intent(in),optional::max_iterations
      class(*),intent(in),optional::context
      real(dp),allocatable::xbest(:),neigh(:,:)
      real(dp)::fbest,f
      integer::mit,it,j,jbest
      mit=1000; if(present(max_iterations)) mit=max_iterations
      allocate(xbest(size(x0))); xbest=x0; call call_objective(objective,xbest,fbest,context)
      do it=1,mit
         if(present(context)) then; call neighbours_fun(xbest,neigh,context); else; call neighbours_fun(xbest,neigh); end if
         jbest=0
         do j=1,size(neigh,2)
            call call_objective(objective,neigh(:,j),f,context)
            if(f<fbest) then; fbest=f; jbest=j; end if
         end do
         if(jbest==0) exit
         xbest=neigh(:,jbest)
      end do
      allocate(result%xbest(size(x0)),result%population_values(1)); result%xbest=xbest
      result%ofvalue=fbest; result%population_values=[fbest]; result%iterations=it; result%status=nmof_ok
   end subroutine greedy_search

   subroutine grid_search(objective, levels, counts, result, context)
      procedure(real_objective)::objective
      real(dp),intent(in)::levels(:,:)
      integer,intent(in)::counts(:)
      type(optimization_result),intent(out)::result
      class(*),intent(in),optional::context
      integer::np,total,i,j,k,tmp,best
      real(dp),allocatable::x(:)
      np=size(counts)
      if(np<1.or.size(levels,2)/=np.or.any(counts<1).or.any(counts>size(levels,1))) then
         result%status=nmof_invalid_input; return
      end if
      total=product(counts); allocate(result%history(total,np),result%population_values(total),x(np))
      do i=0,total-1
         tmp=i
         do j=1,np
            k=mod(tmp,counts(j))+1; tmp=tmp/counts(j); x(j)=levels(k,j)
         end do
         result%history(i+1,:)=x; call call_objective(objective,x,result%population_values(i+1),context)
      end do
      best=minloc(result%population_values,dim=1); allocate(result%xbest(np)); result%xbest=result%history(best,:)
      result%ofvalue=result%population_values(best); result%iterations=total; result%status=nmof_ok
   end subroutine grid_search

   subroutine restart_opt(optimizer, n, best_result, all_values, context)
      procedure(optimizer_callback)::optimizer
      integer,intent(in)::n
      type(optimization_result),intent(out)::best_result
      real(dp),allocatable,intent(out),optional::all_values(:)
      class(*),intent(in),optional::context
      type(optimization_result)::tmp
      integer::i
      if(n<1) then; best_result%status=nmof_invalid_input; return; end if
      if(present(all_values)) allocate(all_values(n))
      best_result%ofvalue=huge(1.0_dp)
      do i=1,n
         if(present(context)) then; call optimizer(tmp,context); else; call optimizer(tmp); end if
         if(present(all_values)) all_values(i)=tmp%ofvalue
         if(tmp%ofvalue<best_result%ofvalue) best_result=tmp
      end do
   end subroutine restart_opt
end module nmof_optimization
