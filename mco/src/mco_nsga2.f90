! SPDX-License-Identifier: GPL-2.0-only
module mco_nsga2
   use mco_kinds, only : dp
   use mco_random, only : seed_random, random_uniform
   use mco_pareto, only : nondominated_sort, crowding_distance
   implicit none
   private
   public :: nsga2_options, nsga2_result, objective_function, constraint_function
   public :: nsga2_optimize

   abstract interface
      subroutine objective_function(x, f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: f(:)
      end subroutine objective_function
      subroutine constraint_function(x, g)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: g(:)
      end subroutine constraint_function
   end interface

   type :: nsga2_options
      integer :: population_size = 100
      integer :: generations = 100
      real(dp) :: crossover_probability = 0.7_dp
      real(dp) :: crossover_distribution = 5.0_dp
      real(dp) :: mutation_probability = -1.0_dp
      real(dp) :: mutation_distribution = 10.0_dp
      integer :: seed = 12345
   end type nsga2_options

   type :: nsga2_result
      real(dp), allocatable :: par(:,:)
      real(dp), allocatable :: value(:,:)
      real(dp), allocatable :: violation(:)
      integer, allocatable :: rank(:)
      real(dp), allocatable :: crowding(:)
      logical, allocatable :: pareto_optimal(:)
      integer :: generations_completed = 0
      integer :: evaluations = 0
      integer :: status = 0
      character(len=160) :: message = ""
   end type nsga2_result
contains
   subroutine nsga2_optimize(objective, nvar, nobj, lower, upper, result, options, constraints, ncon)
      procedure(objective_function) :: objective
      integer, intent(in) :: nvar, nobj
      real(dp), intent(in) :: lower(nvar), upper(nvar)
      type(nsga2_result), intent(out) :: result
      type(nsga2_options), intent(in), optional :: options
      procedure(constraint_function), optional :: constraints
      integer, intent(in), optional :: ncon
      type(nsga2_options) :: opt
      real(dp), allocatable :: x(:,:), f(:,:), cv(:), child_x(:,:), child_f(:,:), child_cv(:)
      real(dp), allocatable :: mix_x(:,:), mix_f(:,:), mix_cv(:), dist(:), mix_dist(:)
      integer, allocatable :: rank(:), mix_rank(:), order(:)
      integer :: np, nc, i, g, neval
      real(dp) :: pmut

      opt = nsga2_options()
      if (present(options)) opt = options
      np = opt%population_size
      if (nvar < 1 .or. nobj < 1) then
         result%status=1; result%message="nvar and nobj must be positive"; return
      end if
      if (np < 4 .or. modulo(np,4) /= 0) then
         result%status=2; result%message="population_size must be a multiple of 4"; return
      end if
      if (opt%generations < 0) then
         result%status=3; result%message="generations must be nonnegative"; return
      end if
      if (any(upper <= lower)) then
         result%status=4; result%message="each upper bound must exceed its lower bound"; return
      end if
      nc=0
      if (present(constraints)) then
         if (.not.present(ncon)) error stop "nsga2_optimize: ncon required with constraints"
         if (ncon < 1) error stop "nsga2_optimize: ncon must be positive"
         nc=ncon
      end if
      pmut=opt%mutation_probability
      if (pmut<0.0_dp) pmut=1.0_dp/real(nvar,dp)
      call seed_random(opt%seed)
      allocate(x(nvar,np),f(nobj,np),cv(np),rank(np),dist(np))
      do i=1,np
         call random_vector(x(:,i))
         x(:,i)=lower+(upper-lower)*x(:,i)
      end do
      call evaluate_population(x,f,cv,objective,constraints,nc,neval)
      call assign_metrics(f,cv,rank,dist)
      allocate(child_x(nvar,np),child_f(nobj,np),child_cv(np))
      allocate(mix_x(nvar,2*np),mix_f(nobj,2*np),mix_cv(2*np))
      allocate(mix_rank(2*np),mix_dist(2*np),order(2*np))
      do g=1,opt%generations
         call make_offspring(x,rank,dist,lower,upper,opt,pmut,child_x)
         call evaluate_population(child_x,child_f,child_cv,objective,constraints,nc,i)
         neval=neval+i
         mix_x(:,1:np)=x; mix_x(:,np+1:)=child_x
         mix_f(:,1:np)=f; mix_f(:,np+1:)=child_f
         mix_cv(1:np)=cv; mix_cv(np+1:)=child_cv
         call assign_metrics(mix_f,mix_cv,mix_rank,mix_dist)
         do i=1,2*np; order(i)=i; end do
         call sort_selection(order,mix_rank,mix_dist)
         do i=1,np
            x(:,i)=mix_x(:,order(i)); f(:,i)=mix_f(:,order(i)); cv(i)=mix_cv(order(i))
         end do
         call assign_metrics(f,cv,rank,dist)
      end do
      allocate(result%par(nvar,np),result%value(nobj,np),result%violation(np))
      allocate(result%rank(np),result%crowding(np),result%pareto_optimal(np))
      result%par=x; result%value=f; result%violation=cv
      result%rank=rank; result%crowding=dist
      result%pareto_optimal=(rank==1 .and. cv<=100.0_dp*epsilon(1.0_dp))
      result%generations_completed=opt%generations
      result%evaluations=neval
      result%status=0; result%message="success"
   end subroutine nsga2_optimize

   subroutine random_vector(x)
      real(dp),intent(out)::x(:)
      integer::i
      do i=1,size(x); x(i)=random_uniform(); end do
   end subroutine random_vector

   subroutine evaluate_population(x,f,cv,objective,constraints,nc,neval)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::f(:,:),cv(:)
      procedure(objective_function)::objective
      procedure(constraint_function),optional::constraints
      integer,intent(in)::nc
      integer,intent(out)::neval
      real(dp),allocatable::g(:)
      integer::i
      if(present(constraints)) allocate(g(nc))
      do i=1,size(x,2)
         call objective(x(:,i),f(:,i))
         if(any(.not.(f(:,i)<huge(1.0_dp)))) f(:,i)=huge(1.0_dp)/100.0_dp
         if(present(constraints)) then
            call constraints(x(:,i),g)
            cv(i)=sum(max(0.0_dp,-g))
         else
            cv(i)=0.0_dp
         end if
      end do
      neval=size(x,2)
   end subroutine evaluate_population

   subroutine assign_metrics(f,cv,rank,dist)
      real(dp),intent(in)::f(:,:),cv(:)
      integer,intent(out)::rank(:)
      real(dp),intent(out)::dist(:)
      call nondominated_sort(f,cv,rank)
      call crowding_distance(f,rank,dist)
   end subroutine assign_metrics

   subroutine sort_selection(order,rank,dist)
      integer,intent(inout)::order(:)
      integer,intent(in)::rank(:)
      real(dp),intent(in)::dist(:)
      integer::i,j,key
      do i=2,size(order)
         key=order(i); j=i-1
         do while(j>=1)
            if(.not.better(key,order(j),rank,dist)) exit
            order(j+1)=order(j); j=j-1
         end do
         order(j+1)=key
      end do
   end subroutine sort_selection

   logical function better(i,j,rank,dist) result(ans)
      integer,intent(in)::i,j,rank(:)
      real(dp),intent(in)::dist(:)
      if(rank(i)<rank(j)) then
         ans=.true.
      else if(rank(i)>rank(j)) then
         ans=.false.
      else if(dist(i)>dist(j)) then
         ans=.true.
      else if(dist(i)<dist(j)) then
         ans=.false.
      else
         ans=i<j
      end if
   end function better

   subroutine make_offspring(x,rank,dist,lower,upper,opt,pmut,child)
      real(dp),intent(in)::x(:,:),dist(:),lower(:),upper(:),pmut
      integer,intent(in)::rank(:)
      type(nsga2_options),intent(in)::opt
      real(dp),intent(out)::child(:,:)
      integer::i,p1,p2
      do i=1,size(x,2),2
         p1=tournament(rank,dist); p2=tournament(rank,dist)
         call sbx(x(:,p1),x(:,p2),lower,upper,opt%crossover_probability, &
                  opt%crossover_distribution,child(:,i),child(:,i+1))
         call polynomial_mutation(child(:,i),lower,upper,pmut,opt%mutation_distribution)
         call polynomial_mutation(child(:,i+1),lower,upper,pmut,opt%mutation_distribution)
      end do
   end subroutine make_offspring

   integer function tournament(rank,dist) result(winner)
      integer,intent(in)::rank(:)
      real(dp),intent(in)::dist(:)
      integer::a,b,n
      n=size(rank); a=1+int(random_uniform()*real(n,dp)); b=1+int(random_uniform()*real(n,dp))
      a=min(a,n); b=min(b,n)
      if(rank(a)<rank(b)) then; winner=a
      else if(rank(b)<rank(a)) then; winner=b
      else if(dist(a)>dist(b)) then; winner=a
      else if(dist(b)>dist(a)) then; winner=b
      else if(random_uniform()<0.5_dp) then; winner=a
      else; winner=b
      end if
   end function tournament

   subroutine sbx(p1,p2,lower,upper,prob,eta,c1,c2)
      real(dp),intent(in)::p1(:),p2(:),lower(:),upper(:),prob,eta
      real(dp),intent(out)::c1(:),c2(:)
      integer::i
      real(dp)::y1,y2,yl,yu,r,beta,alpha,betaq
      if(random_uniform()>prob) then; c1=p1; c2=p2; return; end if
      do i=1,size(p1)
         if(random_uniform()>0.5_dp .or. abs(p1(i)-p2(i))<=1.0e-14_dp) then
            c1(i)=p1(i); c2(i)=p2(i); cycle
         end if
         y1=min(p1(i),p2(i)); y2=max(p1(i),p2(i)); yl=lower(i); yu=upper(i)
         r=random_uniform(); beta=1.0_dp+2.0_dp*(y1-yl)/(y2-y1)
         alpha=2.0_dp-beta**(-(eta+1.0_dp))
         if(r<=1.0_dp/alpha) then; betaq=(r*alpha)**(1.0_dp/(eta+1.0_dp))
         else; betaq=(1.0_dp/(2.0_dp-r*alpha))**(1.0_dp/(eta+1.0_dp)); end if
         c1(i)=0.5_dp*((y1+y2)-betaq*(y2-y1))
         beta=1.0_dp+2.0_dp*(yu-y2)/(y2-y1); alpha=2.0_dp-beta**(-(eta+1.0_dp))
         if(r<=1.0_dp/alpha) then; betaq=(r*alpha)**(1.0_dp/(eta+1.0_dp))
         else; betaq=(1.0_dp/(2.0_dp-r*alpha))**(1.0_dp/(eta+1.0_dp)); end if
         c2(i)=0.5_dp*((y1+y2)+betaq*(y2-y1))
         c1(i)=max(yl,min(yu,c1(i))); c2(i)=max(yl,min(yu,c2(i)))
         if(random_uniform()<0.5_dp) then; r=c1(i); c1(i)=c2(i); c2(i)=r; end if
      end do
   end subroutine sbx

   subroutine polynomial_mutation(x,lower,upper,prob,eta)
      real(dp),intent(inout)::x(:)
      real(dp),intent(in)::lower(:),upper(:),prob,eta
      integer::i
      real(dp)::delta1,delta2,r,mut_pow,xy,val,deltaq
      mut_pow=1.0_dp/(eta+1.0_dp)
      do i=1,size(x)
         if(random_uniform()>prob) cycle
         delta1=(x(i)-lower(i))/(upper(i)-lower(i)); delta2=(upper(i)-x(i))/(upper(i)-lower(i))
         r=random_uniform()
         if(r<=0.5_dp) then
            xy=1.0_dp-delta1; val=2.0_dp*r+(1.0_dp-2.0_dp*r)*xy**(eta+1.0_dp)
            deltaq=val**mut_pow-1.0_dp
         else
            xy=1.0_dp-delta2; val=2.0_dp*(1.0_dp-r)+2.0_dp*(r-0.5_dp)*xy**(eta+1.0_dp)
            deltaq=1.0_dp-val**mut_pow
         end if
         x(i)=max(lower(i),min(upper(i),x(i)+deltaq*(upper(i)-lower(i))))
      end do
   end subroutine polynomial_mutation
end module mco_nsga2
