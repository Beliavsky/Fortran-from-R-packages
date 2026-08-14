! High-level multiobjective evolutionary algorithms translated from rmoo.
module rmoo_core
  use ga_kinds, only : dp
  use ga_random, only : ga_seed, runif, randint
  use ga_operators, only : random_real_population
  use ga_operators, only : random_binary_population, random_perm_population
  use rmoo_pareto, only : non_dominated_sort, crowding_distance
  use rmoo_reference, only : generate_reference_points
  use rmoo_operators, only : tournament_nsga1, tournament_nsga2, tournament_rank
  use rmoo_operators, only : sbx_crossover, polynomial_mutation
  use rmoo_operators, only : single_point_crossover_int
  use rmoo_operators, only : hux_crossover, uniform_crossover_int
  use rmoo_operators, only : ox_crossover, inversion_mutation
  use rmoo_operators, only : random_binary_mutation, uniform_integer_mutation
  use rmoo_survival, only : nsga2_survivors, nsga3_survivors
  use rmoo_survival, only : rnsga2_survivors, sharing_dummy_fitness
  implicit none
  private

  integer, parameter, public :: ALG_NSGA1 = 1
  integer, parameter, public :: ALG_NSGA2 = 2
  integer, parameter, public :: ALG_NSGA3 = 3
  integer, parameter, public :: ALG_RNSGA2 = 4

  integer, parameter, public :: REP_BINARY = 1
  integer, parameter, public :: REP_INTEGER = 2
  integer, parameter, public :: REP_PERMUTATION = 3

  type, public :: rmoo_real_result
    real(dp), allocatable :: population(:,:)
    real(dp), allocatable :: fitness(:,:)
    real(dp), allocatable :: crowding(:)
    integer, allocatable :: rank(:)
    real(dp), allocatable :: reference_points(:,:)
    real(dp), allocatable :: ideal_point(:)
    real(dp), allocatable :: nadir_point(:)
    integer :: iterations = 0
    integer :: algorithm = ALG_NSGA2
  end type rmoo_real_result

  type, public :: rmoo_integer_result
    integer, allocatable :: population(:,:)
    real(dp), allocatable :: fitness(:,:)
    real(dp), allocatable :: crowding(:)
    integer, allocatable :: rank(:)
    real(dp), allocatable :: reference_points(:,:)
    real(dp), allocatable :: ideal_point(:)
    real(dp), allocatable :: nadir_point(:)
    integer :: iterations = 0
    integer :: algorithm = ALG_NSGA2
    integer :: representation = REP_INTEGER
  end type rmoo_integer_result

  abstract interface
    subroutine rmoo_real_objective(x, f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
    end subroutine rmoo_real_objective

    subroutine rmoo_integer_objective(x, f)
      import dp
      integer, intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
    end subroutine rmoo_integer_objective
  end interface

  public :: rmoo_optimize_real
  public :: rmoo_optimize_integer
  public :: rmoo_optimize_binary
  public :: rmoo_optimize_permutation

contains

  subroutine rmoo_optimize_real(objective, lower, upper, nobj, popsize, &
      maxiter, result, algorithm, pcrossover, pmutation, reference_dirs, &
      epsilon, weights, seed, eta_c, eta_m, suggestions)
    procedure(rmoo_real_objective) :: objective
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(in) :: nobj, popsize, maxiter
    type(rmoo_real_result), intent(out) :: result
    integer, intent(in), optional :: algorithm, seed
    real(dp), intent(in), optional :: pcrossover, pmutation
    real(dp), intent(in), optional :: reference_dirs(:,:)
    real(dp), intent(in), optional :: epsilon, weights(:), eta_c, eta_m
    real(dp), intent(in), optional :: suggestions(:,:)

    real(dp), allocatable :: pop(:,:), off(:,:), comb(:,:)
    real(dp), allocatable :: fit(:,:), offfit(:,:), combfit(:,:)
    real(dp), allocatable :: crowd(:), ref(:,:), weightv(:), dummy(:)
    real(dp), allocatable :: ideal(:), worst(:), nadir(:), smin(:)
    real(dp), allocatable :: extreme(:,:), mut(:), child1(:), child2(:)
    integer, allocatable :: rank(:), sel(:), keep(:)
    real(dp) :: pc, pm, eps, ec, em
    integer :: alg, nvar, iter, i, p1, p2, ng

    if (size(lower) /= size(upper)) then
      error stop "rmoo_optimize_real: lower/upper mismatch"
    end if
    if (any(upper <= lower)) then
      error stop "rmoo_optimize_real: invalid bounds"
    end if
    if (popsize < 2 .or. nobj < 1 .or. maxiter < 0) then
      error stop "rmoo_optimize_real: invalid size"
    end if

    nvar = size(lower)
    alg = ALG_NSGA2
    if (present(algorithm)) alg = algorithm
    pc = 0.8_dp
    if (present(pcrossover)) pc = pcrossover
    pm = 0.1_dp
    if (present(pmutation)) pm = pmutation
    eps = 1.0e-3_dp
    if (present(epsilon)) eps = epsilon
    ec = 20.0_dp
    if (present(eta_c)) ec = eta_c
    em = 20.0_dp
    if (present(eta_m)) em = eta_m
    if (present(seed)) call ga_seed(seed)

    allocate(pop(popsize,nvar), off(popsize,nvar))
    allocate(fit(popsize,nobj), offfit(popsize,nobj))
    allocate(rank(popsize), crowd(popsize), sel(popsize), mut(nvar))
    allocate(dummy(popsize))
    allocate(child1(nvar),child2(nvar))
    call random_real_population(pop, lower, upper)
    if (present(suggestions)) then
      if (size(suggestions,2) /= nvar) then
        error stop "rmoo_optimize_real: suggestions dimension mismatch"
      end if
      ng = min(popsize,size(suggestions,1))
      if (ng > 0) pop(1:ng,:) = suggestions(1:ng,:)
    end if
    call evaluate_real_population(objective, pop, fit)

    allocate(ideal(nobj), worst(nobj), nadir(nobj), smin(nobj))
    allocate(extreme(nobj,nobj), weightv(nobj))
    ideal = huge(1.0_dp)
    worst = -huge(1.0_dp)
    nadir = 0.0_dp
    smin = huge(1.0_dp)
    extreme = 0.0_dp

    call prepare_reference_dirs(nobj, popsize, alg, reference_dirs, ref)
    if (present(weights)) then
      if (size(weights) /= nobj) then
        error stop "rmoo_optimize_real: weights size mismatch"
      end if
      weightv = weights
    else
      weightv = 1.0_dp/real(nobj,dp)
    end if

    call refresh_rank(fit, rank, crowd)

    do iter = 1, maxiter
      select case (alg)
      case (ALG_NSGA1)
        call sharing_dummy_fitness(fit,rank,0.5_dp, &
          0.1_dp*real(popsize,dp),dummy)
        call tournament_nsga1(rank,dummy,sel)
      case (ALG_NSGA2)
        call tournament_nsga2(rank,crowd,sel)
      case default
        call tournament_rank(rank,sel)
      end select

      i = 1
      do while (i <= popsize)
        p1 = sel(i)
        if (i < popsize) then
          p2 = sel(i+1)
          if (runif() < pc) then
            call sbx_crossover(pop(p1,:), pop(p2,:), lower, upper, &
              child1, child2, ec, 0.5_dp)
            off(i,:) = child1
            off(i+1,:) = child2
          else
            off(i,:) = pop(p1,:)
            off(i+1,:) = pop(p2,:)
          end if
        else
          off(i,:) = pop(p1,:)
        end if

        if (runif() < pm) then
          call polynomial_mutation(off(i,:), lower, upper, mut, em, &
            min(1.0_dp,1.0_dp/real(nvar,dp)))
          off(i,:) = mut
        end if
        if (i < popsize) then
          if (runif() < pm) then
            call polynomial_mutation(off(i+1,:), lower, upper, mut, em, &
              min(1.0_dp,1.0_dp/real(nvar,dp)))
            off(i+1,:) = mut
          end if
        end if
        i = i + 2
      end do

      call evaluate_real_population(objective, off, offfit)
      if (alg == ALG_NSGA1) then
        pop = off
        fit = offfit
      else
        allocate(comb(2*popsize,nvar), combfit(2*popsize,nobj))
        allocate(keep(popsize))
        comb(1:popsize,:) = pop
        comb(popsize+1:,:) = off
        combfit(1:popsize,:) = fit
        combfit(popsize+1:,:) = offfit

        call choose_survivors(combfit, popsize, alg, ref, eps, weightv, &
          ideal, worst, smin, extreme, nadir, keep)
        do i = 1, popsize
          pop(i,:) = comb(keep(i),:)
          fit(i,:) = combfit(keep(i),:)
        end do
        deallocate(comb, combfit, keep)
      end if
      call refresh_rank(fit, rank, crowd)
    end do

    result%population = pop
    result%fitness = fit
    result%rank = rank
    result%crowding = crowd
    result%reference_points = ref
    result%iterations = maxiter
    result%algorithm = alg
    result%ideal_point = minval(fit,dim=1)
    result%nadir_point = maxval(fit,dim=1)
  end subroutine rmoo_optimize_real

  subroutine rmoo_optimize_binary(objective, nbits, nobj, popsize, maxiter, &
      result, algorithm, pcrossover, pmutation, reference_dirs, epsilon, &
      weights, seed, suggestions)
    procedure(rmoo_integer_objective) :: objective
    integer, intent(in) :: nbits, nobj, popsize, maxiter
    type(rmoo_integer_result), intent(out) :: result
    integer, intent(in), optional :: algorithm, seed
    real(dp), intent(in), optional :: pcrossover, pmutation
    real(dp), intent(in), optional :: reference_dirs(:,:), epsilon, weights(:)
    integer, intent(in), optional :: suggestions(:,:)
    integer, allocatable :: lower(:), upper(:)

    allocate(lower(nbits), upper(nbits))
    lower = 0
    upper = 1
    call optimize_integer_common(objective, lower, upper, nobj, popsize, &
      maxiter, result, REP_BINARY, algorithm, pcrossover, pmutation, &
      reference_dirs, epsilon, weights, seed, suggestions)
  end subroutine rmoo_optimize_binary

  subroutine rmoo_optimize_integer(objective, lower, upper, nobj, popsize, &
      maxiter, result, algorithm, pcrossover, pmutation, reference_dirs, &
      epsilon, weights, seed, suggestions)
    procedure(rmoo_integer_objective) :: objective
    integer, intent(in) :: lower(:), upper(:)
    integer, intent(in) :: nobj, popsize, maxiter
    type(rmoo_integer_result), intent(out) :: result
    integer, intent(in), optional :: algorithm, seed
    real(dp), intent(in), optional :: pcrossover, pmutation
    real(dp), intent(in), optional :: reference_dirs(:,:), epsilon, weights(:)
    integer, intent(in), optional :: suggestions(:,:)

    call optimize_integer_common(objective, lower, upper, nobj, popsize, &
      maxiter, result, REP_INTEGER, algorithm, pcrossover, pmutation, &
      reference_dirs, epsilon, weights, seed, suggestions)
  end subroutine rmoo_optimize_integer

  subroutine rmoo_optimize_permutation(objective, lower_value, upper_value, &
      nobj, popsize, maxiter, result, algorithm, pcrossover, pmutation, &
      reference_dirs, epsilon, weights, seed, suggestions)
    procedure(rmoo_integer_objective) :: objective
    integer, intent(in) :: lower_value, upper_value, nobj, popsize, maxiter
    type(rmoo_integer_result), intent(out) :: result
    integer, intent(in), optional :: algorithm, seed
    real(dp), intent(in), optional :: pcrossover, pmutation
    real(dp), intent(in), optional :: reference_dirs(:,:), epsilon, weights(:)
    integer, intent(in), optional :: suggestions(:,:)
    integer, allocatable :: lower(:), upper(:)
    integer :: nvar

    if (upper_value < lower_value) then
      error stop "rmoo_optimize_permutation: invalid range"
    end if
    nvar = upper_value - lower_value + 1
    allocate(lower(nvar), upper(nvar))
    lower = lower_value
    upper = upper_value
    call optimize_integer_common(objective, lower, upper, nobj, popsize, &
      maxiter, result, REP_PERMUTATION, algorithm, pcrossover, pmutation, &
      reference_dirs, epsilon, weights, seed, suggestions)
  end subroutine rmoo_optimize_permutation

  subroutine optimize_integer_common(objective, lower, upper, nobj, popsize, &
      maxiter, result, representation, algorithm, pcrossover, pmutation, &
      reference_dirs, epsilon, weights, seed, suggestions)
    procedure(rmoo_integer_objective) :: objective
    integer, intent(in) :: lower(:), upper(:)
    integer, intent(in) :: nobj, popsize, maxiter, representation
    type(rmoo_integer_result), intent(out) :: result
    integer, intent(in), optional :: algorithm, seed
    real(dp), intent(in), optional :: pcrossover, pmutation
    real(dp), intent(in), optional :: reference_dirs(:,:), epsilon, weights(:)
    integer, intent(in), optional :: suggestions(:,:)

    integer, allocatable :: pop(:,:), off(:,:), comb(:,:), mut(:)
    integer, allocatable :: child1(:), child2(:)
    integer, allocatable :: rank(:), sel(:), keep(:)
    real(dp), allocatable :: fit(:,:), offfit(:,:), combfit(:,:)
    real(dp), allocatable :: crowd(:), ref(:,:), weightv(:), dummy(:)
    real(dp), allocatable :: ideal(:), worst(:), nadir(:), smin(:)
    real(dp), allocatable :: extreme(:,:)
    real(dp) :: pc, pm, eps
    integer :: alg, nvar, iter, i, j, p1, p2, ng

    nvar = size(lower)
    if (size(upper) /= nvar) then
      error stop "optimize_integer_common: bounds mismatch"
    end if
    if (any(upper < lower)) then
      error stop "optimize_integer_common: invalid bounds"
    end if

    alg = ALG_NSGA2
    if (present(algorithm)) alg = algorithm
    pc = 0.8_dp
    if (present(pcrossover)) pc = pcrossover
    pm = 0.1_dp
    if (present(pmutation)) pm = pmutation
    eps = 1.0e-3_dp
    if (present(epsilon)) eps = epsilon
    if (present(seed)) call ga_seed(seed)

    allocate(pop(popsize,nvar), off(popsize,nvar), mut(nvar))
    allocate(child1(nvar),child2(nvar),dummy(popsize))
    allocate(fit(popsize,nobj), offfit(popsize,nobj))
    allocate(rank(popsize), crowd(popsize), sel(popsize))

    select case (representation)
    case (REP_BINARY)
      call random_binary_population(pop)
    case (REP_PERMUTATION)
      call random_perm_population(pop,lower(1))
    case default
      do i = 1, popsize
        do j = 1, nvar
          pop(i,j) = randint(lower(j),upper(j))
        end do
      end do
    end select
    if (present(suggestions)) then
      if (size(suggestions,2) /= nvar) then
        error stop "optimize_integer_common: suggestions dimension mismatch"
      end if
      ng = min(popsize,size(suggestions,1))
      if (ng > 0) pop(1:ng,:) = suggestions(1:ng,:)
    end if
    call evaluate_integer_population(objective,pop,fit)

    allocate(ideal(nobj), worst(nobj), nadir(nobj), smin(nobj))
    allocate(extreme(nobj,nobj), weightv(nobj))
    ideal = huge(1.0_dp)
    worst = -huge(1.0_dp)
    nadir = 0.0_dp
    smin = huge(1.0_dp)
    extreme = 0.0_dp
    call prepare_reference_dirs(nobj,popsize,alg,reference_dirs,ref)

    if (present(weights)) then
      if (size(weights) /= nobj) then
        error stop "optimize_integer_common: weights size mismatch"
      end if
      weightv = weights
    else
      weightv = 1.0_dp/real(nobj,dp)
    end if

    call refresh_rank(fit,rank,crowd)
    do iter = 1, maxiter
      select case (alg)
      case (ALG_NSGA1)
        call sharing_dummy_fitness(fit,rank,0.5_dp, &
          0.1_dp*real(popsize,dp),dummy)
        call tournament_nsga1(rank,dummy,sel)
      case (ALG_NSGA2)
        call tournament_nsga2(rank,crowd,sel)
      case default
        call tournament_rank(rank,sel)
      end select

      i = 1
      do while (i <= popsize)
        p1 = sel(i)
        if (i < popsize) then
          p2 = sel(i+1)
          if (runif() < pc) then
            select case (representation)
            case (REP_PERMUTATION)
              call ox_crossover(pop(p1,:),pop(p2,:),child1,child2)
              off(i,:)=child1
              off(i+1,:)=child2
            case (REP_BINARY)
              call single_point_crossover_int(pop(p1,:),pop(p2,:),child1,child2)
              off(i,:)=child1
              off(i+1,:)=child2
            case default
              call uniform_crossover_int(pop(p1,:),pop(p2,:),child1,child2)
              off(i,:)=child1
              off(i+1,:)=child2
            end select
          else
            off(i,:) = pop(p1,:)
            off(i+1,:) = pop(p2,:)
          end if
        else
          off(i,:) = pop(p1,:)
        end if

        if (runif() < pm) then
          call mutate_integer(representation,off(i,:),lower,upper,mut)
          off(i,:) = mut
        end if
        if (i < popsize) then
          if (runif() < pm) then
            call mutate_integer(representation,off(i+1,:),lower,upper,mut)
            off(i+1,:) = mut
          end if
        end if
        i = i + 2
      end do

      call evaluate_integer_population(objective,off,offfit)
      if (alg == ALG_NSGA1) then
        pop = off
        fit = offfit
      else
        allocate(comb(2*popsize,nvar),combfit(2*popsize,nobj))
        allocate(keep(popsize))
        comb(1:popsize,:) = pop
        comb(popsize+1:,:) = off
        combfit(1:popsize,:) = fit
        combfit(popsize+1:,:) = offfit
        call choose_survivors(combfit,popsize,alg,ref,eps,weightv, &
          ideal,worst,smin,extreme,nadir,keep)
        do i = 1, popsize
          pop(i,:) = comb(keep(i),:)
          fit(i,:) = combfit(keep(i),:)
        end do
        deallocate(comb,combfit,keep)
      end if
      call refresh_rank(fit,rank,crowd)
    end do

    result%population = pop
    result%fitness = fit
    result%rank = rank
    result%crowding = crowd
    result%reference_points = ref
    result%iterations = maxiter
    result%algorithm = alg
    result%representation = representation
    result%ideal_point = minval(fit,dim=1)
    result%nadir_point = maxval(fit,dim=1)
  end subroutine optimize_integer_common

  subroutine mutate_integer(representation,parent,lower,upper,mutant)
    integer, intent(in) :: representation, parent(:), lower(:), upper(:)
    integer, intent(out) :: mutant(size(parent))
    select case (representation)
    case (REP_PERMUTATION)
      call inversion_mutation(parent,mutant)
    case (REP_BINARY)
      call random_binary_mutation(parent,mutant)
    case default
      call uniform_integer_mutation(parent,lower,upper,mutant, &
        min(1.0_dp,1.0_dp/real(size(parent),dp)))
    end select
  end subroutine mutate_integer

  subroutine prepare_reference_dirs(nobj,popsize,algorithm,input_dirs,ref)
    integer, intent(in) :: nobj, popsize, algorithm
    real(dp), intent(in), optional :: input_dirs(:,:)
    real(dp), allocatable, intent(out) :: ref(:,:)
    integer :: h

    if (algorithm /= ALG_NSGA3 .and. algorithm /= ALG_RNSGA2) then
      allocate(ref(0,nobj))
      return
    end if
    if (present(input_dirs)) then
      if (size(input_dirs,2) /= nobj) then
        error stop "prepare_reference_dirs: objective dimension mismatch"
      end if
      allocate(ref(size(input_dirs,1),nobj))
      ref = input_dirs
      return
    end if
    h = 1
    do while (reference_count(nobj,h) < popsize)
      h = h + 1
    end do
    call generate_reference_points(nobj,h,ref)
  end subroutine prepare_reference_dirs

  integer function reference_count(m,h) result(n)
    integer, intent(in) :: m,h
    integer :: i,k
    integer(kind=8) :: v
    k = min(m-1,h)
    v = 1_8
    do i = 1, k
      v = v*int(h+m-1-k+i,8)/int(i,8)
    end do
    if (v > huge(n)) error stop "reference_count: overflow"
    n = int(v)
  end function reference_count

  subroutine choose_survivors(fitness,target,algorithm,ref,epsilon,weights, &
      ideal,worst,smin,extreme,nadir,keep)
    real(dp), intent(in) :: fitness(:,:),ref(:,:),epsilon,weights(:)
    integer, intent(in) :: target,algorithm
    real(dp), intent(inout) :: ideal(:),worst(:),smin(:),extreme(:,:),nadir(:)
    integer, intent(out) :: keep(target)

    select case (algorithm)
    case (ALG_NSGA2)
      call nsga2_survivors(fitness,target,keep)
    case (ALG_NSGA3)
      call nsga3_survivors(fitness,target,ref,ideal,worst,smin, &
        extreme,nadir,keep)
    case (ALG_RNSGA2)
      call rnsga2_survivors(fitness,target,ref,epsilon,keep,weights)
    case default
      error stop "choose_survivors: unsupported algorithm"
    end select
  end subroutine choose_survivors

  subroutine evaluate_real_population(objective,pop,fitness)
    procedure(rmoo_real_objective) :: objective
    real(dp), intent(in) :: pop(:,:)
    real(dp), intent(out) :: fitness(:,:)
    integer :: i
    do i = 1, size(pop,1)
      call objective(pop(i,:),fitness(i,:))
    end do
  end subroutine evaluate_real_population

  subroutine evaluate_integer_population(objective,pop,fitness)
    procedure(rmoo_integer_objective) :: objective
    integer, intent(in) :: pop(:,:)
    real(dp), intent(out) :: fitness(:,:)
    integer :: i
    do i = 1, size(pop,1)
      call objective(pop(i,:),fitness(i,:))
    end do
  end subroutine evaluate_integer_population

  subroutine refresh_rank(fitness,rank,crowding)
    real(dp), intent(in) :: fitness(:,:)
    integer, intent(out) :: rank(:)
    real(dp), intent(out) :: crowding(:)
    call non_dominated_sort(fitness,rank)
    call crowding_distance(fitness,rank,crowding)
  end subroutine refresh_rank

end module rmoo_core
