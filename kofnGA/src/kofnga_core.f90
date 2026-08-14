module kofnga_core
  use kofnga_kinds, only : dp
  use kofnga_rng, only : rng_state
  use kofnga_types, only : kofnga_control, kofnga_result, kofnga_summary, objective_function
  implicit none
  private
  public :: kofn_ga, summarize_result, mutation_probability

contains

  pure function mutation_probability(k, mutfrac) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: mutfrac
    real(dp) :: p
    if (k <= 0) then
      p = 0.0_dp
    else
      p = 1.0_dp - (1.0_dp-mutfrac)**(1.0_dp/real(k,dp))
    end if
  end function mutation_probability

  subroutine kofn_ga(n, k, objective, result, control, initpop)
    integer, intent(in) :: n, k
    procedure(objective_function) :: objective
    type(kofnga_result), intent(out) :: result
    type(kofnga_control), intent(in), optional :: control
    integer, intent(in), optional :: initpop(:,:)

    type(kofnga_control) :: ctl
    type(rng_state) :: rng
    integer :: popsize, keepbest, ngen, tourneysize, gen, i, j, nchosen
    integer :: p1, p2, idx, alltime, last_best_gen
    integer, allocatable :: pop(:,:), offspring(:,:), combo(:), tourney(:), chosen_pos(:)
    integer, allocatable :: old_order(:), new_order(:), unused(:), final_order(:), tmpk(:), pa(:), pb(:)
    logical, allocatable :: ischosen(:), used(:)
    real(dp), allocatable :: fit_old(:), fit_new(:), tie(:)
    real(dp) :: mutprob, bestval

    ctl = kofnga_control()
    if (present(control)) ctl = control
    popsize = ctl%popsize
    keepbest = ctl%keepbest
    if (keepbest < 0) keepbest = popsize/10
    ngen = ctl%ngen
    tourneysize = ctl%tourneysize
    if (tourneysize < 0) tourneysize = max((popsize+9)/10,2)
    mutprob = ctl%mutprob
    if (ctl%mutfrac >= 0.0_dp) mutprob = mutation_probability(k, ctl%mutfrac)
    call validate_inputs(n,k,popsize,keepbest,ngen,tourneysize,mutprob,initpop)
    call rng%seed(ctl%seed)

    allocate(pop(popsize,k), offspring(popsize,k), fit_old(popsize), fit_new(popsize))
    allocate(result%best_history(ngen+1,k), result%obj_history(ngen+1), result%avg_history(ngen+1))
    allocate(combo(2*k), tourney(tourneysize), ischosen(k), chosen_pos(k), used(n), unused(n))
    allocate(old_order(popsize), new_order(popsize), tie(popsize), final_order(popsize), tmpk(k), pa(k), pb(k))

    if (present(initpop)) then
      pop = initpop
    else
      do i=1,popsize
        call rng%sample_without_replacement(n,k,tmpk)
        pop(i,:)=tmpk
      end do
    end if
    call evaluate_population(pop,objective,fit_old)
    call record_generation(pop,fit_old,rng,result%best_history(1,:),result%obj_history(1),result%avg_history(1))

    do gen=1,ngen
      do i=1,popsize
        call rng%sample_without_replacement(popsize,tourneysize,tourney)
        p1 = tournament_winner(tourney,fit_old,rng)
        call rng%sample_without_replacement(popsize,tourneysize,tourney)
        p2 = tournament_winner(tourney,fit_old,rng)
        pa=pop(p1,:); pb=pop(p2,:)
        call make_offspring(pa,pb,k,rng,tmpk)
        offspring(i,:)=tmpk

        ischosen = .false.
        nchosen = 0
        do j=1,k
          if (rng%uniform() < mutprob) then
            nchosen = nchosen + 1
            ischosen(j) = .true.
            chosen_pos(nchosen) = j
          end if
        end do
        if (nchosen > 0) then
          used = .false.
          do j=1,k
            used(offspring(i,j)) = .true.
          end do
          idx = 0
          do j=1,n
            if (.not.used(j)) then
              idx = idx + 1
              unused(idx) = j
            end if
          end do
          if (nchosen > idx) error stop "kofn_ga: mutation requested more unused indices than available"
          call rng%shuffle(unused(1:idx))
          do j=1,nchosen
            offspring(i,chosen_pos(j)) = unused(j)
          end do
        end if
      end do

      call evaluate_population(offspring,objective,fit_new)
      if (keepbest == 0) then
        pop = offspring
        fit_old = fit_new
      else
        do i=1,popsize
          tie(i)=rng%uniform()
        end do
        call order_by_value_tie(fit_old,tie,old_order)
        do i=1,popsize
          tie(i)=rng%uniform()
        end do
        call order_by_value_tie(fit_new,tie,new_order)
        pop(1:keepbest,:) = pop(old_order(1:keepbest),:)
        fit_old(1:keepbest) = fit_old(old_order(1:keepbest))
        pop(keepbest+1:popsize,:) = offspring(new_order(1:popsize-keepbest),:)
        fit_old(keepbest+1:popsize) = fit_new(new_order(1:popsize-keepbest))
      end if
      call record_generation(pop,fit_old,rng,result%best_history(gen+1,:), &
                             result%obj_history(gen+1),result%avg_history(gen+1))
      if (ctl%verbose > 0) then
        if (mod(gen,ctl%verbose)==0) then
          write(*,'(a,i0,a,es16.8)') 'Finished iteration ',gen,'. Best OF value = ',result%obj_history(gen+1)
        end if
      end if
    end do

    do i=1,popsize
      tie(i)=real(i,dp)
    end do
    call order_by_value_tie(fit_old,tie,final_order)
    allocate(result%pop(popsize,k),result%obj(popsize),result%bestsol(k))
    do i=1,popsize
      tmpk=pop(final_order(i),:)
      call sort_int(tmpk)
      result%pop(i,:)=tmpk
      result%obj(i) = fit_old(final_order(i))
    end do

    bestval = minval(result%obj_history)
    last_best_gen = 1
    do i=1,ngen+1
      if (same_real(result%obj_history(i),bestval)) last_best_gen = i
    end do
    alltime = last_best_gen
    result%bestsol = result%best_history(alltime,:)
    result%bestobj = result%obj_history(alltime)
  end subroutine kofn_ga

  subroutine validate_inputs(n,k,popsize,keepbest,ngen,tourneysize,mutprob,initpop)
    integer,intent(in)::n,k,popsize,keepbest,ngen,tourneysize
    real(dp),intent(in)::mutprob
    integer,intent(in),optional::initpop(:,:)
    integer :: i,j
    logical, allocatable :: seen(:)
    if(n<=0 .or. k<=0 .or. k>n) error stop "kofn_ga: require n >= k > 0"
    if(popsize<2 .or. keepbest<0 .or. keepbest>=popsize) error stop "kofn_ga: invalid population/elitism sizes"
    if(ngen<1) error stop "kofn_ga: ngen must be positive"
    if(tourneysize<2 .or. tourneysize>=popsize) error stop "kofn_ga: invalid tournament size"
    if(mutprob<0.0_dp .or. mutprob>1.0_dp) error stop "kofn_ga: mutprob must be in [0,1]"
    if(present(initpop)) then
      if(size(initpop,1)/=popsize .or. size(initpop,2)/=k) error stop "kofn_ga: initpop has wrong shape"
      allocate(seen(n))
      do i=1,popsize
        seen=.false.
        do j=1,k
          if(initpop(i,j)<1 .or. initpop(i,j)>n) error stop "kofn_ga: initpop index out of range"
          if(seen(initpop(i,j))) error stop "kofn_ga: duplicate index within initpop row"
          seen(initpop(i,j))=.true.
        end do
      end do
    end if
  end subroutine validate_inputs

  subroutine evaluate_population(pop, objective, f)
    integer, intent(in) :: pop(:,:)
    procedure(objective_function) :: objective
    real(dp), intent(out) :: f(size(pop,1))
    integer :: i
    integer, allocatable :: subset(:)
    allocate(subset(size(pop,2)))
    do i=1,size(pop,1)
      subset=pop(i,:)
      f(i)=objective(subset)
    end do
  end subroutine evaluate_population

  function tournament_winner(tourney,fitness,rng) result(winner)
    integer,intent(in)::tourney(:)
    real(dp),intent(in)::fitness(:)
    type(rng_state),intent(inout)::rng
    integer::winner,i,j,m,ntie
    real(dp),allocatable::w(:)
    real(dp)::u,cum,total,fi,fj
    m=size(tourney); allocate(w(m)); w=1.0_dp
    do i=1,m
      fi=-fitness(tourney(i))
      w(i)=1.0_dp
      ntie=1
      do j=1,m
        if(j==i) cycle
        fj=-fitness(tourney(j))
        if(fj<fi) then
          w(i)=w(i)+1.0_dp
        else if(same_real(fj,fi)) then
          ntie=ntie+1
        end if
      end do
      w(i)=w(i)+0.5_dp*real(ntie-1,dp)
    end do
    total=sum(w); u=rng%uniform()*total; cum=0.0_dp; winner=tourney(m)
    do i=1,m
      cum=cum+w(i)
      if(u<=cum) then
        winner=tourney(i); return
      end if
    end do
  end function tournament_winner

  subroutine make_offspring(a,b,k,rng,child)
    integer,intent(in)::a(:),b(:),k
    type(rng_state),intent(inout)::rng
    integer,intent(out)::child(k)
    integer::combo(2*k),nc,i,j
    logical::found
    nc=0
    do i=1,size(a)
      nc=nc+1; combo(nc)=a(i)
    end do
    do i=1,size(b)
      found=.false.
      do j=1,nc
        if(combo(j)==b(i)) then; found=.true.; exit; end if
      end do
      if(.not.found) then; nc=nc+1; combo(nc)=b(i); end if
    end do
    call rng%shuffle(combo(1:nc))
    child=combo(1:k)
  end subroutine make_offspring

  subroutine record_generation(pop,fit,rng,best,bestobj,avgobj)
    integer,intent(in)::pop(:,:)
    real(dp),intent(in)::fit(:)
    type(rng_state),intent(inout)::rng
    integer,intent(out)::best(:)
    real(dp),intent(out)::bestobj,avgobj
    integer,allocatable::candidates(:)
    integer::i,nc,chosen
    bestobj=minval(fit); avgobj=sum(fit)/real(size(fit),dp)
    allocate(candidates(size(fit))); nc=0
    do i=1,size(fit)
      if(same_real(fit(i),bestobj)) then; nc=nc+1; candidates(nc)=i; end if
    end do
    chosen=candidates(rng%randint(1,nc))
    best=pop(chosen,:); call sort_int(best)
  end subroutine record_generation

  subroutine order_by_value_tie(x,tie,ord)
    real(dp),intent(in)::x(:),tie(:)
    integer,intent(out)::ord(size(x))
    integer::i,j,key
    do i=1,size(x); ord(i)=i; end do
    do i=2,size(x)
      key=ord(i); j=i-1
      do while(j>=1)
        if(x(ord(j)) < x(key)) exit
        if(same_real(x(ord(j)),x(key))) then
          if(tie(ord(j)) <= tie(key)) exit
        end if
        ord(j+1)=ord(j); j=j-1
      end do
      ord(j+1)=key
    end do
  end subroutine order_by_value_tie

  subroutine sort_int(x)
    integer,intent(inout)::x(:)
    integer::i,j,key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_int

  pure logical function same_real(a,b) result(eq)
    real(dp), intent(in) :: a,b
    eq = (a <= b .and. a >= b)
  end function same_real

  function summarize_result(result) result(s)
    type(kofnga_result),intent(in)::result
    type(kofnga_summary)::s
    integer::i,j,n,u
    logical::same
    s%generations=size(result%obj_history)-1
    allocate(s%best_solution(size(result%bestsol))); s%best_solution=result%bestsol
    s%initial_average=result%avg_history(1); s%initial_minimum=result%obj_history(1)
    s%final_average=result%avg_history(size(result%avg_history)); s%final_minimum=result%obj_history(size(result%obj_history))
    s%best_generation=0
    do i=1,size(result%obj_history)
      if(same_real(result%obj_history(i),result%bestobj)) then; s%best_generation=i-1; exit; end if
    end do
    n=size(result%pop,1); u=0
    do i=1,n
      same=.false.
      do j=1,i-1
        if(all(result%pop(i,:)==result%pop(j,:))) then; same=.true.; exit; end if
      end do
      if(.not.same) u=u+1
    end do
    s%unique_final=u
  end function summarize_result

end module kofnga_core
