! cmaes-fortran - modern Fortran translation of cmaes 1.0-12
! Original R package authors: Heike Trautmann, Olaf Mersmann, David Arnu.
! Source algorithm based by the R package on N. Hansen's purecmaes.m.
! License: GPL-2.0-only. See COPYING and original/.
module cmaes
  use iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf
  use cmaes_kinds, only : dp
  use cmaes_rng, only : rng_state, rng_seed, rng_normal
  use cmaes_linalg, only : symmetric_eigen_descending
  implicit none
  private

  public :: dp, cma_control, cma_result, cma_es, extract_population
  public :: cma_objective, cma_vector_objective

  abstract interface
    function cma_objective(x) result(value)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function cma_objective

    subroutine cma_vector_objective(x, value)
      import :: dp
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: value(:)
    end subroutine cma_vector_objective
  end interface

  type :: cma_control
    logical :: trace = .false.
    real(dp) :: fnscale = 1.0_dp
    real(dp) :: stopfitness = -huge(1.0_dp)
    integer :: maxit = 0
    real(dp) :: sigma = 0.5_dp
    real(dp) :: stop_tolx = -1.0_dp
    logical :: keep_best = .true.
    logical :: vectorized = .false.
    integer :: lambda = 0
    integer :: mu = 0
    real(dp), allocatable :: weights(:)
    real(dp) :: mueff = -1.0_dp
    real(dp) :: ccum = -1.0_dp
    real(dp) :: cs = -1.0_dp
    real(dp) :: ccov_mu = -1.0_dp
    real(dp) :: ccov_1 = -1.0_dp
    real(dp) :: damps = -1.0_dp
    logical :: diag = .false.
    logical :: diag_sigma = .false.
    logical :: diag_eigen = .false.
    logical :: diag_value = .false.
    logical :: diag_pop = .false.
    integer(int64) :: seed = 1_int64
  end type cma_control

  type :: cma_result
    real(dp), allocatable :: par(:)
    real(dp) :: value = huge(1.0_dp)
    integer :: function_evaluations = 0
    integer :: gradient_evaluations = 0
    integer :: convergence = 1
    character(len=:), allocatable :: message
    integer :: constraint_violations = 0
    integer :: iterations = 0
    real(dp), allocatable :: sigma_history(:)
    real(dp), allocatable :: eigen_history(:, :)
    real(dp), allocatable :: value_history(:, :)
    real(dp), allocatable :: population_history(:, :, :)
  end type cma_result
contains
  function cma_es(par, fn, lower, upper, control, fn_vector) result(res)
    real(dp), intent(in) :: par(:)
    procedure(cma_objective) :: fn
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(cma_control), intent(in), optional :: control
    procedure(cma_vector_objective), optional :: fn_vector
    type(cma_result) :: res

    type(cma_control) :: ctrl
    type(rng_state) :: rng
    integer :: n, lambda, mu, maxiter, iter, counteval, cviol
    integer :: j, k, flat_idx, info
    real(dp) :: sigma, stop_tolx, fnscale, stopfitness
    real(dp) :: mueff, cc, cs, mucov, ccov, damps, chin, psnorm
    real(dp) :: infv
    logical :: hsig, log_sigma, log_eigen, log_value, log_pop
    logical :: has_best
    real(dp), allocatable :: lo(:), hi(:), weights(:)
    real(dp), allocatable :: xmean(:), pc(:), ps(:), b(:, :), d(:), bd(:, :), c(:, :)
    real(dp), allocatable :: arz(:, :), arx(:, :), vx(:, :), y(:), pen(:), arfitness(:)
    integer, allocatable :: arindex(:), aripop(:)
    real(dp), allocatable :: selx(:, :), selz(:, :), zmean(:), bdz(:, :), rankmu(:, :)
    real(dp), allocatable :: evals(:), evecs(:, :)
    real(dp), allocatable :: sigma_log(:), eigen_log(:, :), value_log(:, :), pop_log(:, :, :)

    n = size(par)
    if (n < 1) error stop "cma_es: par must be nonempty"
    ctrl = cma_control()
    if (present(control)) ctrl = control
    if (ctrl%fnscale <= 0.0_dp .and. ctrl%fnscale >= 0.0_dp) &
      error stop "cma_es: fnscale must be nonzero"
    fnscale = ctrl%fnscale
    stopfitness = ctrl%stopfitness
    maxiter = ctrl%maxit
    if (maxiter <= 0) maxiter = 100 * n * n
    sigma = ctrl%sigma
    if (sigma <= 0.0_dp) error stop "cma_es: sigma must be positive"
    stop_tolx = ctrl%stop_tolx
    if (stop_tolx < 0.0_dp) stop_tolx = 1.0e-12_dp * sigma

    allocate(lo(n), hi(n))
    infv = ieee_value(1.0_dp, ieee_positive_inf)
    lo = -infv
    hi = infv
    if (present(lower)) then
      if (size(lower) == 1) then
        lo = lower(1)
      else if (size(lower) == n) then
        lo = lower
      else
        error stop "cma_es: lower must have size 1 or size(par)"
      end if
    end if
    if (present(upper)) then
      if (size(upper) == 1) then
        hi = upper(1)
      else if (size(upper) == n) then
        hi = upper
      else
        error stop "cma_es: upper must have size 1 or size(par)"
      end if
    end if
    if (any(lo >= hi)) error stop "cma_es: each lower bound must be below upper bound"

    lambda = ctrl%lambda
    if (lambda <= 0) lambda = 4 + floor(3.0_dp * log(real(n, dp)))
    mu = ctrl%mu
    if (mu <= 0) mu = lambda / 2
    if (mu < 1 .or. lambda < mu) error stop "cma_es: require lambda >= mu >= 1"

    allocate(weights(mu))
    if (allocated(ctrl%weights)) then
      if (size(ctrl%weights) /= mu) error stop "cma_es: weights must have length mu"
      weights = ctrl%weights
    else
      do j = 1, mu
        weights(j) = log(real(mu + 1, dp)) - log(real(j, dp))
      end do
    end if
    if (sum(weights) <= 0.0_dp .and. sum(weights) >= 0.0_dp) &
      error stop "cma_es: weights sum to zero"
    weights = weights / sum(weights)

    mueff = ctrl%mueff
    if (mueff <= 0.0_dp) mueff = sum(weights)**2 / sum(weights**2)
    cc = ctrl%ccum
    if (cc <= 0.0_dp) cc = 4.0_dp / real(n + 4, dp)
    cs = ctrl%cs
    if (cs <= 0.0_dp) cs = (mueff + 2.0_dp) / (real(n, dp) + mueff + 3.0_dp)
    mucov = ctrl%ccov_mu
    if (mucov <= 0.0_dp) mucov = mueff
    ccov = ctrl%ccov_1
    if (ccov <= 0.0_dp) then
      ccov = (1.0_dp / mucov) * 2.0_dp / (real(n, dp) + 1.4_dp)**2 + &
        (1.0_dp - 1.0_dp / mucov) * ((2.0_dp * mucov - 1.0_dp) / &
        ((real(n + 2, dp))**2 + 2.0_dp * mucov))
    end if
    damps = ctrl%damps
    if (damps <= 0.0_dp) damps = 1.0_dp + 2.0_dp * max(0.0_dp, &
      sqrt((mueff - 1.0_dp) / real(n + 1, dp)) - 1.0_dp) + cs

    log_sigma = ctrl%diag .or. ctrl%diag_sigma
    log_eigen = ctrl%diag .or. ctrl%diag_eigen
    log_value = ctrl%diag .or. ctrl%diag_value
    log_pop = ctrl%diag .or. ctrl%diag_pop

    allocate(sigma_log(0), eigen_log(0, 0), value_log(0, 0), pop_log(0, 0, 0))
    if (log_sigma) then
      deallocate(sigma_log)
      allocate(sigma_log(maxiter))
    end if
    if (log_eigen) then
      deallocate(eigen_log)
      allocate(eigen_log(maxiter, n))
    end if
    if (log_value) then
      deallocate(value_log)
      allocate(value_log(maxiter, mu))
    end if
    if (log_pop) then
      deallocate(pop_log)
      allocate(pop_log(n, mu, maxiter))
    end if

    allocate(xmean(n), pc(n), ps(n), b(n, n), d(n), bd(n, n), c(n, n))
    allocate(arz(n, lambda), arx(n, lambda), vx(n, lambda), y(lambda), pen(lambda))
    allocate(arfitness(lambda), arindex(lambda), aripop(mu))
    allocate(selx(n, mu), selz(n, mu), zmean(n), bdz(n, mu), rankmu(n, n))

    xmean = par
    pc = 0.0_dp
    ps = 0.0_dp
    b = 0.0_dp
    do j = 1, n
      b(j, j) = 1.0_dp
    end do
    d = 1.0_dp
    bd = b
    c = matmul(bd, transpose(bd))
    chin = sqrt(real(n, dp)) * (1.0_dp - 1.0_dp / (4.0_dp * real(n, dp)) + &
      1.0_dp / (21.0_dp * real(n * n, dp)))

    call rng_seed(rng, ctrl%seed)
    res%value = huge(1.0_dp)
    allocate(res%par(n))
    res%par = par
    has_best = .false.
    counteval = 0
    cviol = 0
    res%message = ""

    do iter = 1, maxiter
      if (.not. ctrl%keep_best) then
        res%value = huge(1.0_dp)
        has_best = .false.
      end if
      if (log_sigma) sigma_log(iter) = sigma

      do k = 1, lambda
        do j = 1, n
          arz(j, k) = rng_normal(rng)
        end do
      end do
      arx = matmul(bd, arz)
      do k = 1, lambda
        arx(:, k) = xmean + sigma * arx(:, k)
      end do

      do k = 1, lambda
        pen(k) = 1.0_dp
        do j = 1, n
          vx(j, k) = min(max(arx(j, k), lo(j)), hi(j))
          pen(k) = pen(k) + (arx(j, k) - vx(j, k))**2
        end do
        if (.not. ieee_is_finite(pen(k))) pen(k) = huge(1.0_dp) / 2.0_dp
        if (pen(k) > 1.0_dp) cviol = cviol + 1
      end do

      if (ctrl%vectorized) then
        if (.not. present(fn_vector)) error stop "cma_es: vectorized=.true. requires fn_vector"
        call fn_vector(vx, y)
        if (size(y) /= lambda) error stop "cma_es: vector objective returned wrong size"
        y = y * fnscale
      else
        do k = 1, lambda
          y(k) = fn(vx(:, k)) * fnscale
        end do
      end if
      counteval = counteval + lambda
      arfitness = y * pen

      do k = 1, lambda
        if (pen(k) <= 1.0_dp) then
          if (.not. has_best .or. y(k) < res%value) then
            res%value = y(k)
            res%par = arx(:, k)
            has_best = .true.
          end if
        end if
      end do

      call argsort_ascending(arfitness, arindex)
      arfitness = arfitness(arindex)
      aripop = arindex(1:mu)
      do k = 1, mu
        selx(:, k) = arx(:, aripop(k))
        selz(:, k) = arz(:, aripop(k))
      end do
      xmean = matmul(selx, weights)
      zmean = matmul(selz, weights)

      if (log_pop) pop_log(:, :, iter) = selx
      if (log_value) value_log(iter, :) = arfitness(aripop)

      ps = (1.0_dp - cs) * ps + sqrt(cs * (2.0_dp - cs) * mueff) * matmul(b, zmean)
      psnorm = sqrt(dot_product(ps, ps))
      hsig = (psnorm / sqrt(1.0_dp - (1.0_dp - cs)**(2.0_dp * &
        real(counteval, dp) / real(lambda, dp))) / chin) < (1.4_dp + 2.0_dp / real(n + 1, dp))
      pc = (1.0_dp - cc) * pc
      if (hsig) pc = pc + sqrt(cc * (2.0_dp - cc) * mueff) * matmul(bd, zmean)

      bdz = matmul(bd, selz)
      rankmu = 0.0_dp
      do k = 1, mu
        do j = 1, n
          rankmu(:, j) = rankmu(:, j) + weights(k) * bdz(:, k) * bdz(j, k)
        end do
      end do
      c = (1.0_dp - ccov) * c + ccov * (1.0_dp / mucov) * &
        (outer(pc, pc) + merge(0.0_dp, 1.0_dp, hsig) * cc * (2.0_dp - cc) * c) + &
        ccov * (1.0_dp - 1.0_dp / mucov) * rankmu

      sigma = sigma * exp((psnorm / chin - 1.0_dp) * cs / damps)
      call symmetric_eigen_descending(c, evals, evecs, info)
      if (info /= 0) then
        res%message = "Symmetric eigensolver failed."
        exit
      end if
      if (log_eigen) eigen_log(iter, :) = evals
      if (any(evals < sqrt(epsilon(1.0_dp)) * abs(evals(1)))) then
        res%message = "Covariance matrix C is numerically not positive definite."
        exit
      end if
      b = evecs
      d = sqrt(evals)
      do j = 1, n
        bd(:, j) = b(:, j) * d(j)
      end do

      if (arfitness(1) <= stopfitness * fnscale) then
        res%message = "Stop fitness reached."
        exit
      end if

      if (all(d < stop_tolx) .and. all(sigma * pc < stop_tolx)) then
        res%message = "All standard deviations smaller than tolerance."
        exit
      end if

      flat_idx = min(1 + lambda / 2, 2 + ceiling(real(lambda, dp) / 4.0_dp))
      if (arfitness(1) <= arfitness(flat_idx) .and. arfitness(1) >= arfitness(flat_idx)) then
        sigma = sigma * exp(0.2_dp + cs / damps)
        if (ctrl%trace) write(*, '(a)') "Flat fitness function. Increasing sigma."
      end if
      if (ctrl%trace) then
        write(*, '(a,i0,a,i0,a,es16.8)') "Iteration ", iter, " of ", maxiter, &
          ": current fitness ", arfitness(1) * fnscale
      end if
    end do

    res%iterations = min(iter, maxiter)
    res%function_evaluations = counteval
    res%constraint_violations = cviol
    res%gradient_evaluations = 0
    if (has_best) then
      res%value = res%value / fnscale
    else
      res%value = huge(1.0_dp)
      if (allocated(res%par)) deallocate(res%par)
      allocate(res%par(0))
    end if
    if (res%iterations >= maxiter) then
      res%convergence = 1
    else
      res%convergence = 0
    end if

    if (log_sigma) then
      allocate(res%sigma_history(res%iterations))
      res%sigma_history = sigma_log(1:res%iterations)
    end if
    if (log_eigen) then
      allocate(res%eigen_history(res%iterations, n))
      res%eigen_history = eigen_log(1:res%iterations, :)
    end if
    if (log_value) then
      allocate(res%value_history(res%iterations, mu))
      res%value_history = value_log(1:res%iterations, :)
    end if
    if (log_pop) then
      allocate(res%population_history(n, mu, res%iterations))
      res%population_history = pop_log(:, :, 1:res%iterations)
    end if
  end function cma_es

  subroutine extract_population(res, iter, population, values)
    type(cma_result), intent(in) :: res
    integer, intent(in) :: iter
    real(dp), allocatable, intent(out) :: population(:, :)
    real(dp), allocatable, intent(out), optional :: values(:)

    if (.not. allocated(res%population_history)) &
      error stop "extract_population: result contains no population history"
    if (iter < 1 .or. iter > size(res%population_history, 3)) &
      error stop "extract_population: iter out of range"
    allocate(population(size(res%population_history, 1), size(res%population_history, 2)))
    population = res%population_history(:, :, iter)
    if (present(values)) then
      if (allocated(res%value_history)) then
        allocate(values(size(res%value_history, 2)))
        values = res%value_history(iter, :)
      else
        allocate(values(0))
      end if
    end if
  end subroutine extract_population

  pure function outer(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: j
    do j = 1, size(y)
      a(:, j) = x * y(j)
    end do
  end function outer

  subroutine argsort_ascending(x, idx)
    real(dp), intent(in) :: x(:)
    integer, intent(out) :: idx(:)
    integer :: i, j, key
    if (size(idx) /= size(x)) error stop "argsort_ascending: size mismatch"
    do i = 1, size(x)
      idx(i) = i
    end do
    do i = 2, size(x)
      key = idx(i)
      j = i - 1
      do while (j >= 1)
        if (x(idx(j)) <= x(key)) exit
        idx(j + 1) = idx(j)
        j = j - 1
      end do
      idx(j + 1) = key
    end do
  end subroutine argsort_ascending
end module cmaes
