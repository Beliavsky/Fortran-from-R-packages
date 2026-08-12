module abcoptim
  use iso_fortran_env, only : int64
  use ieee_arithmetic, only : ieee_is_finite
  use abc_rng, only : abc_rng_state
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

  type, public :: abc_control
    integer :: food_number = 20
    integer :: limit = 100
    integer :: max_cycle = 1000
    integer :: criter = 50
    logical :: optiinteger = .false.
    integer(int64) :: seed = 123456789_int64
    logical :: legacy_r_scaling = .true.
  end type abc_control

  type, public :: abc_result
    real(dp), allocatable :: foods(:,:)
    real(dp), allocatable :: f(:)
    real(dp), allocatable :: fitness(:)
    integer, allocatable :: trial(:)
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: hist(:,:)
    real(dp) :: value = huge(1.0_dp)
    integer :: counts = 0
    integer :: objective_evaluations = 0
    integer :: scouts = 0
    logical :: converged = .false.
  end type abc_result

  abstract interface
    function abc_objective(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function abc_objective
  end interface

  public :: abc_optim, abc_cpp, calculate_fitness

contains

  function calculate_fitness(value) result(fit)
    real(dp), intent(in) :: value
    real(dp) :: fit

    if (value >= 0.0_dp) then
      fit = 1.0_dp / (value + 1.0_dp)
    else
      fit = 1.0_dp + abs(value)
    end if
  end function calculate_fitness

  subroutine abc_optim(par, fn, lb, ub, result, control, parscale, fnscale)
    real(dp), intent(in) :: par(:), lb(:), ub(:)
    procedure(abc_objective) :: fn
    type(abc_result), intent(out) :: result
    type(abc_control), intent(in), optional :: control
    real(dp), intent(in), optional :: parscale(:)
    real(dp), intent(in), optional :: fnscale

    type(abc_control) :: ctl

    ctl = abc_control()
    if (present(control)) ctl = control
    call run_abc(par, fn, lb, ub, result, ctl, .true., parscale, fnscale)
  end subroutine abc_optim

  subroutine abc_cpp(par, fn, lb, ub, result, control, parscale, fnscale)
    real(dp), intent(in) :: par(:), lb(:), ub(:)
    procedure(abc_objective) :: fn
    type(abc_result), intent(out) :: result
    type(abc_control), intent(in), optional :: control
    real(dp), intent(in), optional :: parscale(:)
    real(dp), intent(in), optional :: fnscale

    type(abc_control) :: ctl

    ctl = abc_control()
    if (present(control)) ctl = control
    ctl%optiinteger = .false.
    call run_abc(par, fn, lb, ub, result, ctl, .false., parscale, fnscale)
  end subroutine abc_cpp

  subroutine run_abc(par, fn, lb_in, ub_in, result, ctl, r_semantics, parscale_in, fnscale_in)
    real(dp), intent(in) :: par(:), lb_in(:), ub_in(:)
    procedure(abc_objective) :: fn
    type(abc_result), intent(out) :: result
    type(abc_control), intent(in) :: ctl
    logical, intent(in) :: r_semantics
    real(dp), intent(in), optional :: parscale_in(:)
    real(dp), intent(in), optional :: fnscale_in

    type(abc_rng_state) :: rng
    real(dp), allocatable :: foods(:,:), f(:), fitness(:), prob(:)
    real(dp), allocatable :: lb(:), ub(:), scale(:), global_params(:), hist_work(:,:)
    real(dp) :: global_min, fnscale
    integer, allocatable :: trial(:)
    integer :: d, nf, iter, persistence, n_hist

    d = size(par)
    nf = ctl%food_number
    call validate_inputs(d, nf, ctl, lb_in, ub_in, parscale_in, fnscale_in)

    allocate(lb(d), ub(d), scale(d))
    call expand_bounds(lb_in, d, lb)
    call expand_bounds(ub_in, d, ub)
    call sanitize_bounds(lb, ub)
    scale = 1.0_dp
    if (present(parscale_in)) then
      if (size(parscale_in) == 1) then
        scale = parscale_in(1)
      else
        scale = parscale_in
      end if
    end if
    fnscale = 1.0_dp
    if (present(fnscale_in)) fnscale = fnscale_in

    allocate(foods(d,nf), f(nf), fitness(nf), prob(nf), trial(nf))
    allocate(global_params(d), hist_work(d, max(1, ctl%max_cycle)))
    call rng%seed(ctl%seed)

    result%objective_evaluations = 0
    result%scouts = 0
    global_params = par
    if (r_semantics .and. ctl%legacy_r_scaling) then
      global_min = eval_raw(par)
    else
      global_min = eval_scaled(par)
    end if

    call initialize_foods()
    persistence = 0
    call memorize_best()

    if (.not. r_semantics) persistence = 0
    n_hist = 0
    if (.not. r_semantics) then
      n_hist = 1
      hist_work(:,1) = global_params
    end if

    iter = 0
    do while (iter + 1 < ctl%max_cycle)
      iter = iter + 1
      call employed_bees()
      call calculate_probabilities()
      if (r_semantics) then
        call onlooker_bees_r()
      else
        call onlooker_bees_cpp()
      end if
      call memorize_best()

      if (r_semantics) then
        n_hist = n_hist + 1
        hist_work(:,n_hist) = global_params
        if (persistence > ctl%criter) exit
      else
        n_hist = n_hist + 1
        hist_work(:,n_hist) = global_params
        if (persistence >= ctl%criter) exit
      end if
      call scout_bee()
    end do

    result%counts = iter
    result%converged = merge(persistence > ctl%criter, persistence >= ctl%criter, r_semantics)
    allocate(result%foods(d,nf), result%f(nf), result%fitness(nf), result%trial(nf), result%par(d))
    result%foods = foods
    result%f = f
    result%fitness = fitness
    result%trial = trial
    result%par = global_params
    allocate(result%hist(d,n_hist))
    if (n_hist > 0) result%hist = hist_work(:,1:n_hist)

    if (r_semantics) then
      result%value = eval_scaled(global_params)
    else
      result%value = global_min
    end if

  contains

    function eval_raw(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y

      y = fn(x)
      result%objective_evaluations = result%objective_evaluations + 1
      if (.not. ieee_is_finite(y)) error stop "ABCoptim: objective returned a non-finite value"
    end function eval_raw

    function eval_scaled(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y

      y = fn(x / scale) / fnscale
      result%objective_evaluations = result%objective_evaluations + 1
      if (.not. ieee_is_finite(y)) error stop "ABCoptim: objective returned a non-finite value"
    end function eval_scaled

    subroutine initialize_foods()
      integer :: ii, j
      real(dp) :: alpha

      do ii = 1, nf
        if (nf == 1) then
          alpha = 0.0_dp
        else
          alpha = real(ii - 1, dp) / real(nf - 1, dp)
        end if
        do j = 1, d
          foods(j,ii) = lb(j) + (ub(j) - lb(j)) * alpha
        end do
        f(ii) = eval_scaled(foods(:,ii))
        fitness(ii) = calculate_fitness(f(ii))
        trial(ii) = 0
      end do
    end subroutine initialize_foods

    subroutine initialize_source(index)
      integer, intent(in) :: index
      integer :: j

      if (ctl%optiinteger .and. r_semantics) then
        do j = 1, d
          foods(j,index) = merge(1.0_dp, 0.0_dp, rng%uniform() > 0.5_dp)
        end do
      else
        do j = 1, d
          foods(j,index) = lb(j) + (ub(j) - lb(j)) * rng%uniform()
        end do
      end if
      f(index) = eval_scaled(foods(:,index))
      fitness(index) = calculate_fitness(f(index))
      trial(index) = 0
      result%scouts = result%scouts + 1
    end subroutine initialize_source

    subroutine propose_and_select(index)
      integer, intent(in) :: index
      real(dp) :: solution(d), obj, fit_new
      integer :: param, neighbour

      param = rng%randint(d)
      neighbour = index
      do while (neighbour == index)
        neighbour = rng%randint(nf)
      end do
      solution = foods(:,index)
      if (ctl%optiinteger .and. r_semantics) then
        solution(param) = merge(1.0_dp, 0.0_dp, rng%uniform() > 0.5_dp)
      else
        solution(param) = foods(param,index) + &
          (foods(param,index) - foods(param,neighbour)) * (rng%uniform() - 0.5_dp) * 2.0_dp
        solution(param) = max(lb(param), min(ub(param), solution(param)))
      end if
      obj = eval_scaled(solution)
      fit_new = calculate_fitness(obj)
      if (fit_new > fitness(index)) then
        foods(:,index) = solution
        f(index) = obj
        fitness(index) = fit_new
        trial(index) = 0
      else
        trial(index) = trial(index) + 1
      end if
    end subroutine propose_and_select

    subroutine employed_bees()
      integer :: ii
      do ii = 1, nf
        call propose_and_select(ii)
      end do
    end subroutine employed_bees

    subroutine calculate_probabilities()
      real(dp) :: maxfit
      maxfit = maxval(fitness)
      if (r_semantics) then
        prob = 0.9_dp * fitness / (maxfit + 1.0e-20_dp) + 0.1_dp
      else
        prob = 0.9_dp * (fitness + 1.0e-40_dp) / (maxfit + 1.0e-40_dp) + 0.1_dp
      end if
    end subroutine calculate_probabilities

    subroutine onlooker_bees_cpp()
      integer :: accepted, index
      accepted = 0
      index = 1
      do while (accepted < nf)
        if (rng%uniform() < prob(index)) then
          accepted = accepted + 1
          call propose_and_select(index)
        else
          index = index + 1
          if (index > nf) index = 1
        end if
      end do
    end subroutine onlooker_bees_cpp

    subroutine onlooker_bees_r()
      integer :: accepted, index
      accepted = 0
      index = 1
      do while (accepted < nf)
        if (rng%uniform() < prob(index)) then
          accepted = accepted + 1
          call propose_and_select(index)
        end if
        index = index + 1
        if (index == nf) index = 1
      end do
    end subroutine onlooker_bees_r

    subroutine memorize_best()
      integer :: ii
      logical :: did_improve

      did_improve = .false.
      do ii = 1, nf
        if (f(ii) < global_min) then
          global_min = f(ii)
          global_params = foods(:,ii)
          did_improve = .true.
        end if
      end do
      if (did_improve) then
        persistence = 0
      else
        persistence = persistence + 1
      end if
    end subroutine memorize_best

    subroutine scout_bee()
      integer :: max_index, ii
      max_index = 1
      do ii = 2, nf
        if (trial(ii) > trial(max_index)) max_index = ii
      end do
      if (trial(max_index) >= ctl%limit) call initialize_source(max_index)
    end subroutine scout_bee

  end subroutine run_abc

  subroutine validate_inputs(d, nf, ctl, lb, ub, parscale, fnscale)
    integer, intent(in) :: d, nf
    type(abc_control), intent(in) :: ctl
    real(dp), intent(in) :: lb(:), ub(:)
    real(dp), intent(in), optional :: parscale(:), fnscale

    if (d < 1) error stop "ABCoptim: par must be nonempty"
    if (nf < 2) error stop "ABCoptim: food_number must be at least 2"
    if (ctl%max_cycle < 2) error stop "ABCoptim: max_cycle must be at least 2"
    if (ctl%limit < 0) error stop "ABCoptim: limit must be nonnegative"
    if (ctl%criter < 0) error stop "ABCoptim: criter must be nonnegative"
    if (.not. (size(lb) == 1 .or. size(lb) == d)) error stop "ABCoptim: invalid lower-bound size"
    if (.not. (size(ub) == 1 .or. size(ub) == d)) error stop "ABCoptim: invalid upper-bound size"
    if (present(parscale)) then
      if (.not. (size(parscale) == 1 .or. size(parscale) == d)) error stop "ABCoptim: invalid parscale size"
      if (any(abs(parscale) <= tiny(1.0_dp))) error stop "ABCoptim: parscale values must be nonzero"
    end if
    if (present(fnscale)) then
      if (abs(fnscale) <= tiny(1.0_dp)) error stop "ABCoptim: fnscale must be nonzero"
    end if
  end subroutine validate_inputs

  subroutine expand_bounds(input, d, output)
    real(dp), intent(in) :: input(:)
    integer, intent(in) :: d
    real(dp), intent(out) :: output(d)

    if (size(input) == 1) then
      output = input(1)
    else
      output = input
    end if
  end subroutine expand_bounds

  subroutine sanitize_bounds(lb, ub)
    real(dp), intent(inout) :: lb(:), ub(:)
    integer :: i
    real(dp), parameter :: finite_big = huge(1.0_dp) * 1.0e-10_dp

    do i = 1, size(lb)
      if (.not. ieee_is_finite(lb(i))) lb(i) = -finite_big
      if (.not. ieee_is_finite(ub(i))) ub(i) = finite_big
      if (lb(i) > ub(i)) error stop "ABCoptim: lower bound exceeds upper bound"
    end do
  end subroutine sanitize_bounds

end module abcoptim
