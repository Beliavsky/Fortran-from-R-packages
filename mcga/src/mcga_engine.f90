! SPDX-License-Identifier: GPL-2.0-or-later
module mcga_engine
  use mcga_kinds, only : dp
  use mcga_random, only : set_random_seed, runif_scalar, randint
  use mcga_bytes, only : uniform_crossover_doubles, byte_code_mutation_doubles
  implicit none
  private

  abstract interface
    function mcga_objective(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function mcga_objective

    subroutine mcga_multi_objective(x, f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
    end subroutine mcga_multi_objective
  end interface

  type, public :: mcga_result
    real(dp), allocatable :: population(:,:) ! (popsize, chsize), sorted by increasing cost
    real(dp), allocatable :: costs(:)
    integer :: generations = 0
  contains
    procedure :: best => mcga_result_best
  end type mcga_result

  type, public :: multi_mcga_result
    real(dp), allocatable :: population(:,:) ! (popsize, chsize), sorted by decreasing upstream rank score
    real(dp), allocatable :: costs(:,:)      ! (popsize, numfunc)
    real(dp), allocatable :: ranks(:)
    integer :: generations = 0
  end type multi_mcga_result

  public :: mcga_optimize, multi_mcga_optimize, calculate_rank_scores
  public :: mcga_objective, mcga_multi_objective

contains
  function mcga_result_best(self) result(x)
    class(mcga_result), intent(in) :: self
    real(dp), allocatable :: x(:)
    if (.not. allocated(self%population)) then
      allocate(x(0))
    else
      allocate(x(size(self%population, 2)))
      x = self%population(1,:)
    end if
  end function mcga_result_best

  subroutine mcga_optimize(popsize, chsize, minval, maxval, objective, result, crossprob, mutateprob, elitism, maxiter, seed)
    integer, intent(in) :: popsize, chsize
    real(dp), intent(in) :: minval, maxval
    procedure(mcga_objective) :: objective
    type(mcga_result), intent(out) :: result
    real(dp), intent(in), optional :: crossprob, mutateprob
    integer, intent(in), optional :: elitism, maxiter, seed

    real(dp), allocatable :: pop(:,:), nextpop(:,:), costs(:), nextcosts(:)
    real(dp) :: cp, mp
    integer :: elite, ngen, gen, i, j

    call validate_common(popsize, chsize, minval, maxval)
    cp = 1.0_dp
    mp = 0.01_dp
    elite = 1
    ngen = 10
    if (present(crossprob)) cp = crossprob
    if (present(mutateprob)) mp = mutateprob
    if (present(elitism)) elite = elitism
    if (present(maxiter)) ngen = maxiter
    if (cp < 0.0_dp .or. cp > 1.0_dp) error stop "mcga_optimize: crossprob must be in [0,1]"
    if (mp < 0.0_dp .or. mp > 1.0_dp) error stop "mcga_optimize: mutateprob must be in [0,1]"
    if (elite < 0 .or. elite >= popsize) error stop "mcga_optimize: invalid elitism"
    if (ngen < 0) error stop "mcga_optimize: negative maxiter"
    if (present(seed)) call set_random_seed(seed)

    allocate(pop(chsize,popsize), nextpop(chsize,popsize), costs(popsize), nextcosts(popsize))
    do i = 1, popsize
      do j = 1, chsize
        pop(j,i) = runif_scalar(minval, maxval)
      end do
    end do
    ! Source-compatible with exported R mcga(): initial costs are zero and
    ! the randomized population is not evaluated before the first tournament.
    costs = 0.0_dp
    nextpop = 0.0_dp
    nextcosts = 0.0_dp

    do gen = 1, ngen
      call tournament_generation(pop, costs, nextpop, nextcosts, cp, mp, elite)
      pop = nextpop
      costs = nextcosts
      call evaluate_population(pop, objective, costs)
    end do

    call sort_single(pop, costs)
    allocate(result%population(popsize,chsize), result%costs(popsize))
    result%population = transpose(pop)
    result%costs = costs
    result%generations = ngen
  end subroutine mcga_optimize

  subroutine multi_mcga_optimize(popsize, chsize, numfunc, minval, maxval, objective, result, &
                                 crossprob, mutateprob, elitism, maxiter, seed)
    integer, intent(in) :: popsize, chsize, numfunc
    real(dp), intent(in) :: minval, maxval
    procedure(mcga_multi_objective) :: objective
    type(multi_mcga_result), intent(out) :: result
    real(dp), intent(in), optional :: crossprob, mutateprob
    integer, intent(in), optional :: elitism, maxiter, seed

    real(dp), allocatable :: pop(:,:), nextpop(:,:), costs(:,:), nextcosts(:,:), ranks(:), nextranks(:)
    real(dp) :: cp, mp
    integer :: elite, ngen, gen, i, j

    call validate_common(popsize, chsize, minval, maxval)
    if (numfunc < 1) error stop "multi_mcga_optimize: numfunc must be positive"
    cp = 1.0_dp
    mp = 0.01_dp
    elite = 1
    ngen = 10
    if (present(crossprob)) cp = crossprob
    if (present(mutateprob)) mp = mutateprob
    if (present(elitism)) elite = elitism
    if (present(maxiter)) ngen = maxiter
    if (cp < 0.0_dp .or. cp > 1.0_dp) error stop "multi_mcga_optimize: crossprob must be in [0,1]"
    if (mp < 0.0_dp .or. mp > 1.0_dp) error stop "multi_mcga_optimize: mutateprob must be in [0,1]"
    if (elite < 0 .or. elite >= popsize) error stop "multi_mcga_optimize: invalid elitism"
    if (ngen < 0) error stop "multi_mcga_optimize: negative maxiter"
    if (present(seed)) call set_random_seed(seed)

    allocate(pop(chsize,popsize), nextpop(chsize,popsize))
    allocate(costs(numfunc,popsize), nextcosts(numfunc,popsize), ranks(popsize), nextranks(popsize))
    do i = 1, popsize
      do j = 1, chsize
        pop(j,i) = runif_scalar(minval, maxval)
      end do
    end do
    costs = 0.0_dp
    ranks = 0.0_dp
    nextpop = 0.0_dp
    nextcosts = 0.0_dp
    nextranks = 0.0_dp

    do gen = 1, ngen
      call multi_tournament_generation(pop, costs, ranks, nextpop, nextcosts, nextranks, cp, mp, elite)
      pop = nextpop
      costs = nextcosts
      ranks = nextranks
      call evaluate_multi_population(pop, objective, costs)
      call calculate_rank_scores(costs, ranks)
    end do

    call calculate_rank_scores(costs, ranks)
    call sort_multi(pop, costs, ranks)
    allocate(result%population(popsize,chsize), result%costs(popsize,numfunc), result%ranks(popsize))
    result%population = transpose(pop)
    result%costs = transpose(costs)
    result%ranks = ranks
    result%generations = ngen
  end subroutine multi_mcga_optimize

  subroutine validate_common(popsize, chsize, minval, maxval)
    integer, intent(in) :: popsize, chsize
    real(dp), intent(in) :: minval, maxval
    if (popsize < 2) error stop "mcga: popsize must be at least 2"
    if (chsize < 1) error stop "mcga: chsize must be positive"
    if (maxval < minval) error stop "mcga: maxval must be >= minval"
  end subroutine validate_common

  subroutine evaluate_population(pop, objective, costs)
    real(dp), intent(in) :: pop(:,:)
    procedure(mcga_objective) :: objective
    real(dp), intent(out) :: costs(:)
    integer :: i
    do i = 1, size(pop,2)
      costs(i) = objective(pop(:,i))
    end do
  end subroutine evaluate_population

  subroutine evaluate_multi_population(pop, objective, costs)
    real(dp), intent(in) :: pop(:,:)
    procedure(mcga_multi_objective) :: objective
    real(dp), intent(out) :: costs(:,:)
    integer :: i
    do i = 1, size(pop,2)
      call objective(pop(:,i), costs(:,i))
    end do
  end subroutine evaluate_multi_population

  subroutine tournament_generation(pop, costs, nextpop, nextcosts, crossprob, mutateprob, elitism)
    real(dp), intent(inout) :: pop(:,:), costs(:)
    real(dp), intent(out) :: nextpop(:,:), nextcosts(:)
    real(dp), intent(in) :: crossprob, mutateprob
    integer, intent(in) :: elitism
    integer :: n, selected, i1, i2, i3, i4, w1, w2, p1, p2, i
    real(dp), allocatable :: c1(:), c2(:)

    n = size(pop,2)
    nextpop = 0.0_dp
    nextcosts = 0.0_dp
    selected = 0
    if (elitism > 0) then
      call sort_single(pop, costs)
      do i = 1, elitism
        nextpop(:,i) = pop(:,i)
        nextcosts(i) = costs(i)
      end do
      selected = elitism
    end if

    do while (selected < n)
      call draw_distinct_pair(n, i1, i2)
      call draw_distinct_pair(n, i3, i4)
      if (costs(i1) < costs(i2)) then
        w1 = i1
      else
        w1 = i2
      end if
      if (costs(i3) < costs(i4)) then
        w2 = i3
      else
        w2 = i4
      end if

      selected = selected + 1
      p1 = selected
      nextpop(:,p1) = pop(:,w1)
      nextcosts(p1) = costs(w1)
      if (selected >= n) exit
      selected = selected + 1
      p2 = selected
      nextpop(:,p2) = pop(:,w2)
      nextcosts(p2) = costs(w2)

      if (runif_scalar(0.0_dp, 1.0_dp) < crossprob) then
        call uniform_crossover_doubles(nextpop(:,p1), nextpop(:,p2), c1, c2)
        nextpop(:,p1) = c1
        nextpop(:,p2) = c2
      end if
      call byte_code_mutation_doubles(nextpop(:,p1), mutateprob)
      call byte_code_mutation_doubles(nextpop(:,p2), mutateprob)
    end do
  end subroutine tournament_generation

  subroutine multi_tournament_generation(pop, costs, ranks, nextpop, nextcosts, nextranks, crossprob, mutateprob, elitism)
    real(dp), intent(inout) :: pop(:,:), costs(:,:), ranks(:)
    real(dp), intent(out) :: nextpop(:,:), nextcosts(:,:), nextranks(:)
    real(dp), intent(in) :: crossprob, mutateprob
    integer, intent(in) :: elitism
    integer :: n, selected, i1, i2, i3, i4, w1, w2, p1, p2, i
    real(dp), allocatable :: c1(:), c2(:)

    n = size(pop,2)
    nextpop = 0.0_dp
    nextcosts = 0.0_dp
    nextranks = 0.0_dp
    selected = 0
    if (elitism > 0) then
      call sort_multi(pop, costs, ranks)
      do i = 1, elitism
        nextpop(:,i) = pop(:,i)
        nextcosts(:,i) = costs(:,i)
        nextranks(i) = ranks(i)
      end do
      selected = elitism
    end if

    do while (selected < n)
      call draw_distinct_pair(n, i1, i2)
      call draw_distinct_pair(n, i3, i4)
      if (ranks(i1) > ranks(i2)) then
        w1 = i1
      else
        w1 = i2
      end if
      ! Preserve upstream multi_mcga.c exactly: the second tournament uses
      ! the opposite comparison, selecting the lower rank score.
      if (ranks(i3) < ranks(i4)) then
        w2 = i3
      else
        w2 = i4
      end if

      selected = selected + 1
      p1 = selected
      nextpop(:,p1) = pop(:,w1)
      nextcosts(:,p1) = costs(:,w1)
      nextranks(p1) = ranks(w1)
      if (selected >= n) exit
      selected = selected + 1
      p2 = selected
      nextpop(:,p2) = pop(:,w2)
      nextcosts(:,p2) = costs(:,w2)
      nextranks(p2) = ranks(w2)

      if (runif_scalar(0.0_dp, 1.0_dp) < crossprob) then
        call uniform_crossover_doubles(nextpop(:,p1), nextpop(:,p2), c1, c2)
        nextpop(:,p1) = c1
        nextpop(:,p2) = c2
      end if
      call byte_code_mutation_doubles(nextpop(:,p1), mutateprob)
      call byte_code_mutation_doubles(nextpop(:,p2), mutateprob)
    end do
  end subroutine multi_tournament_generation

  subroutine draw_distinct_pair(n, i1, i2)
    integer, intent(in) :: n
    integer, intent(out) :: i1, i2
    i1 = randint(1, n)
    do
      i2 = randint(1, n)
      if (i2 /= i1) exit
    end do
  end subroutine draw_distinct_pair

  subroutine calculate_rank_scores(costs, ranks)
    real(dp), intent(in) :: costs(:,:) ! (numfunc,popsize)
    real(dp), intent(out) :: ranks(:)
    integer :: i, j, h, n

    n = size(costs,2)
    if (size(ranks) /= n) error stop "calculate_rank_scores: shape mismatch"
    ranks = 0.0_dp
    do i = 1, n
      do j = 1, n
        do h = 1, size(costs,1)
          if (costs(h,i) < costs(h,j)) then
            ranks(i) = ranks(i) + 1.0_dp
            exit
          end if
        end do
      end do
    end do
  end subroutine calculate_rank_scores

  subroutine sort_single(pop, costs)
    real(dp), intent(inout) :: pop(:,:), costs(:)
    integer :: i, j, k, n
    real(dp) :: tc
    real(dp), allocatable :: tx(:)

    n = size(costs)
    allocate(tx(size(pop,1)))
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if (costs(j) < costs(k)) k = j
      end do
      if (k /= i) then
        tx = pop(:,i)
        pop(:,i) = pop(:,k)
        pop(:,k) = tx
        tc = costs(i)
        costs(i) = costs(k)
        costs(k) = tc
      end if
    end do
  end subroutine sort_single

  subroutine sort_multi(pop, costs, ranks)
    real(dp), intent(inout) :: pop(:,:), costs(:,:), ranks(:)
    integer :: i, j, k, n
    real(dp) :: tr
    real(dp), allocatable :: tx(:), tc(:)

    n = size(ranks)
    allocate(tx(size(pop,1)), tc(size(costs,1)))
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if (ranks(j) > ranks(k)) k = j
      end do
      if (k /= i) then
        tx = pop(:,i)
        pop(:,i) = pop(:,k)
        pop(:,k) = tx
        tc = costs(:,i)
        costs(:,i) = costs(:,k)
        costs(:,k) = tc
        tr = ranks(i)
        ranks(i) = ranks(k)
        ranks(k) = tr
      end if
    end do
  end subroutine sort_multi
end module mcga_engine
