module rmalschains_local_search
  use rmalschains_kinds, only : dp
  use rmalschains_rng, only : rng_state, rng_uniform, rng_normal, rng_int
  use rmalschains_types, only : mals_control, mals_objective, ls_state
  use rmalschains_operators, only : nearest_distance, nearest_coordinate_distances
  use rmalschains_linalg, only : symmetric_eigen_descending
  implicit none
  private
  public :: apply_local_search, reset_ls_state
contains
  subroutine reset_ls_state(state)
    type(ls_state), intent(inout) :: state
    state%initialized = .false.
    state%kind = 0
    state%delta = 0.0_dp
    state%num_failed = 0
    state%num_success = 0
    state%sigma = 1.0_dp
    state%lambda = 0
    state%mu = 0
    state%cma_evals = 0
    if (allocated(state%delta_vec)) deallocate(state%delta_vec)
    if (allocated(state%delta_init)) deallocate(state%delta_init)
    if (allocated(state%bias)) deallocate(state%bias)
    if (allocated(state%simplex)) deallocate(state%simplex)
    if (allocated(state%simplex_fit)) deallocate(state%simplex_fit)
    if (allocated(state%xmean)) deallocate(state%xmean)
    if (allocated(state%pc)) deallocate(state%pc)
    if (allocated(state%ps)) deallocate(state%ps)
    if (allocated(state%c)) deallocate(state%c)
    if (allocated(state%b)) deallocate(state%b)
    if (allocated(state%d)) deallocate(state%d)
    if (allocated(state%bd)) deallocate(state%bd)
    if (allocated(state%weights)) deallocate(state%weights)
  end subroutine reset_ls_state

  subroutine apply_local_search(name, state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    character(len=*), intent(in) :: name
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:)
    real(dp), intent(inout) :: fx
    integer, intent(in) :: maxeval, idx
    real(dp), intent(in) :: pop(:, :), lower(:), upper(:)
    type(rng_state), intent(inout) :: rng
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    character(len=:), allocatable :: key

    key = trim(adjustl(name))
    select case (key)
    case ('sw')
      call apply_sw(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    case ('ssw')
      call apply_ssw(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    case ('simplex')
      call apply_simplex(state, x, fx, maxeval, lower, upper, fn, ctrl, used)
    case ('cmaes', 'cmaesmyrandom', 'cmaesalways')
      call apply_cmaes(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    case ('mts1', 'mtsls1')
      call apply_mts1(state, x, fx, maxeval, pop, idx, lower, upper, fn, ctrl, used)
    case ('mts2')
      call apply_mts2(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    case default
      error stop "apply_local_search: unsupported local search"
    end select
  end subroutine apply_local_search

  pure logical function target_hit(fx, ctrl)
    real(dp), intent(in) :: fx
    type(mals_control), intent(in) :: ctrl
    target_hit = fx <= ctrl%optimum + ctrl%threshold
  end function target_hit

  subroutine apply_sw(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval, idx
    real(dp), intent(in) :: pop(:, :), lower(:), upper(:)
    type(rng_state), intent(inout) :: rng
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    real(dp), allocatable :: dif(:), trial(:)
    real(dp) :: fnew, dist
    integer :: i

    if (.not. state%initialized .or. state%kind /= 1) then
      call reset_ls_state(state)
      state%kind = 1
      state%initialized = .true.
      allocate(state%bias(size(x)))
      state%bias = 0.0_dp
      dist = nearest_distance(x, pop, idx)
      state%delta = min(0.2_dp, 0.5_dp * dist)
      if (state%delta <= 0.0_dp) state%delta = 0.2_dp
    end if
    allocate(dif(size(x)), trial(size(x)))
    used = 0
    do while (used < maxeval .and. .not. target_hit(fx, ctrl))
      do i = 1, size(x)
        dif(i) = rng_normal(rng, state%delta)
      end do
      trial = min(max(x + state%bias + dif, lower), upper)
      fnew = fn(trial)
      used = used + 1
      if (fnew < fx) then
        x = trial
        fx = fnew
        state%bias = 0.2_dp * state%bias + 0.4_dp * (dif + state%bias)
        state%num_success = state%num_success + 1
        state%num_failed = 0
      else if (used < maxeval) then
        trial = min(max(x - state%bias - dif, lower), upper)
        fnew = fn(trial)
        used = used + 1
        if (fnew < fx) then
          x = trial
          fx = fnew
          state%bias = state%bias - 0.4_dp * (dif + state%bias)
          state%num_success = state%num_success + 1
          state%num_failed = 0
        else
          state%bias = 0.5_dp * state%bias
          state%num_failed = state%num_failed + 1
          state%num_success = 0
        end if
      end if
      if (state%num_success >= 5) then
        state%num_success = 0
        state%delta = 2.0_dp * state%delta
      else if (state%num_failed >= 3) then
        state%num_failed = 0
        state%delta = 0.5_dp * state%delta
      end if
    end do
  end subroutine apply_sw

  subroutine apply_ssw(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval, idx
    real(dp), intent(in) :: pop(:, :), lower(:), upper(:)
    type(rng_state), intent(inout) :: rng
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    real(dp), allocatable :: dif(:), trial(:), cdist(:)
    logical, allocatable :: changed(:)
    real(dp) :: fnew
    integer :: i, start

    if (.not. state%initialized .or. state%kind /= 2) then
      call reset_ls_state(state)
      state%kind = 2
      state%initialized = .true.
      allocate(state%bias(size(x)), state%delta_vec(size(x)), state%delta_init(size(x)), cdist(size(x)))
      call nearest_coordinate_distances(x, pop, cdist, idx)
      do i = 1, size(x)
        state%delta_init(i) = 0.5_dp * cdist(i)
        if (state%delta_init(i) <= 0.0_dp) state%delta_init(i) = 0.4_dp
        state%delta_init(i) = min(max(state%delta_init(i), 1.0e-15_dp), 0.4_dp)
      end do
      state%delta_vec = state%delta_init
      state%bias = 0.0_dp
    end if
    allocate(dif(size(x)), trial(size(x)), changed(size(x)))
    changed = .false.
    start = rng_int(rng, 1, size(x))
    do i = 0, size(x) - 1
      if (rng_uniform(rng) < 0.1_dp) changed(1 + modulo(start - 1 + i, size(x))) = .true.
    end do
    used = 0
    do while (used < maxeval .and. .not. target_hit(fx, ctrl))
      do i = 1, size(x)
        if (changed(i)) then
          dif(i) = rng_normal(rng, state%delta_vec(i))
        else
          dif(i) = 0.0_dp
        end if
      end do
      trial = min(max(x + state%bias + dif, lower), upper)
      fnew = fn(trial)
      used = used + 1
      if (fnew < fx) then
        x = trial
        fx = fnew
        state%bias = 0.2_dp * state%bias + 0.4_dp * (dif + state%bias)
        state%num_success = state%num_success + 1
        state%num_failed = 0
      else if (used < maxeval) then
        trial = min(max(x - state%bias - dif, lower), upper)
        fnew = fn(trial)
        used = used + 1
        if (fnew < fx) then
          x = trial
          fx = fnew
          do i = 1, size(x)
            if (changed(i)) state%bias(i) = state%bias(i) - 0.4_dp * (dif(i) + state%bias(i))
          end do
          state%num_success = state%num_success + 1
          state%num_failed = 0
        else
          do i = 1, size(x)
            if (changed(i)) state%bias(i) = 0.5_dp * state%bias(i)
          end do
          state%num_failed = state%num_failed + 1
          state%num_success = 0
        end if
      end if
      if (state%num_success >= 5) then
        state%num_success = 0
        do i = 1, size(x)
          if (changed(i)) state%delta_vec(i) = min(0.4_dp, 2.0_dp * state%delta_vec(i))
        end do
      else if (state%num_failed >= 3) then
        state%num_failed = 0
        do i = 1, size(x)
          if (changed(i)) then
            state%delta_vec(i) = 0.5_dp * state%delta_vec(i)
            if (state%delta_vec(i) < 1.0e-15_dp) state%delta_vec(i) = state%delta_init(i)
          end if
        end do
      end if
    end do
  end subroutine apply_ssw

  subroutine apply_mts1(state, x, fx, maxeval, pop, idx, lower, upper, fn, ctrl, used)
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval, idx
    real(dp), intent(in) :: pop(:, :), lower(:), upper(:)
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    real(dp) :: old, fnew, dist, initial
    integer :: i
    logical :: improved
    if (.not. state%initialized .or. state%kind /= 5) then
      call reset_ls_state(state)
      state%kind = 5
      state%initialized = .true.
      dist = nearest_distance(x, pop, idx)
      state%delta = min(0.4_dp, 0.5_dp * dist)
      if (state%delta <= 0.0_dp) state%delta = 0.2_dp
      allocate(state%delta_init(1))
      state%delta_init(1) = state%delta
      state%num_success = 1
    end if
    initial = state%delta_init(1)
    used = 0
    do while (used < maxeval .and. .not. target_hit(fx, ctrl))
      if (state%num_success == 0) then
        state%delta = 0.5_dp * state%delta
        if (state%delta < 1.0e-15_dp) state%delta = initial
      end if
      improved = .false.
      do i = 1, size(x)
        if (used >= maxeval) exit
        old = x(i)
        x(i) = min(max(old - state%delta, lower(i)), upper(i))
        fnew = fn(x); used = used + 1
        if (fnew < fx) then
          fx = fnew; improved = .true.
        else
          x(i) = old
          if (used < maxeval) then
            x(i) = min(max(old + 0.5_dp * state%delta, lower(i)), upper(i))
            fnew = fn(x); used = used + 1
            if (fnew < fx) then
              fx = fnew; improved = .true.
            else
              x(i) = old
            end if
          end if
        end if
        if (target_hit(fx, ctrl)) exit
      end do
      state%num_success = merge(1, 0, improved)
    end do
  end subroutine apply_mts1

  subroutine apply_mts2(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval, idx
    real(dp), intent(in) :: pop(:, :), lower(:), upper(:)
    type(rng_state), intent(inout) :: rng
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    real(dp), allocatable :: oldx(:), trial(:)
    integer, allocatable :: dsgn(:), r(:)
    real(dp) :: fnew, dist, initial
    integer :: i
    logical :: improved
    if (.not. state%initialized .or. state%kind /= 6) then
      call reset_ls_state(state)
      state%kind = 6
      state%initialized = .true.
      dist = nearest_distance(x, pop, idx)
      state%delta = min(0.4_dp, 0.5_dp * dist)
      if (state%delta <= 0.0_dp) state%delta = 0.2_dp
      allocate(state%delta_init(1))
      state%delta_init(1) = state%delta
      state%num_success = 1
    end if
    initial = state%delta_init(1)
    if (state%num_success == 0) then
      state%delta = 0.5_dp * state%delta
      if (state%delta < 1.0e-15_dp) state%delta = initial
    end if
    state%num_success = 0
    allocate(oldx(size(x)), trial(size(x)), dsgn(size(x)), r(size(x)))
    used = 0
    do while (used < maxeval .and. .not. target_hit(fx, ctrl))
      oldx = x
      do i = 1, size(x)
        dsgn(i) = merge(1, -1, rng_uniform(rng) >= 0.5_dp)
        r(i) = rng_int(rng, 0, 3)
      end do
      trial = x
      do i = 1, size(x)
        if (r(i) == 0) trial(i) = min(max(trial(i) - real(dsgn(i), dp) * state%delta, lower(i)), upper(i))
      end do
      fnew = fn(trial); used = used + 1
      improved = fnew < fx
      if (improved) then
        x = trial; fx = fnew; state%num_success = 1
      else if (used < maxeval) then
        trial = oldx
        do i = 1, size(x)
          if (r(i) == 0) trial(i) = min(max(trial(i) + 0.5_dp * real(dsgn(i), dp) * state%delta, lower(i)), upper(i))
        end do
        fnew = fn(trial); used = used + 1
        if (fnew < fx) then
          x = trial; fx = fnew; state%num_success = 1
        end if
      end if
    end do
  end subroutine apply_mts2

  subroutine apply_simplex(state, x, fx, maxeval, lower, upper, fn, ctrl, used)
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval
    real(dp), intent(in) :: lower(:), upper(:)
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    integer :: n, i, best, worst, second, k
    real(dp), allocatable :: centroid_sum(:), trial(:)
    real(dp) :: fnew, rfac, factor
    n = size(x)
    used = 0
    if (.not. state%initialized .or. state%kind /= 3) then
      call reset_ls_state(state)
      state%kind = 3
      state%initialized = .true.
      if (maxeval < n) then
        call coordinate_fallback(x, fx, maxeval, lower, upper, fn, ctrl, used)
        state%initialized = .false.
        return
      end if
      allocate(state%simplex(n, n + 1), state%simplex_fit(n + 1))
      state%simplex(:, 1) = x
      state%simplex_fit(1) = fx
      do i = 1, n
        state%simplex(:, i + 1) = x
        state%simplex(i, i + 1) = min(max(x(i) + 0.1_dp * (upper(i) - lower(i)), lower(i)), upper(i))
        state%simplex_fit(i + 1) = fn(state%simplex(:, i + 1))
        used = used + 1
      end do
    end if
    allocate(centroid_sum(n), trial(n))
    factor = 1.0_dp
    do while (used < maxeval .and. .not. target_hit(fx, ctrl))
      call simplex_extremes(state%simplex_fit, best, second, worst)
      centroid_sum = sum(state%simplex, dim=2)
      rfac = (1.0_dp + factor) / real(n, dp)
      trial = centroid_sum * rfac - state%simplex(:, worst) * (rfac + factor)
      trial = min(max(trial, lower), upper)
      fnew = fn(trial); used = used + 1
      if (fnew < state%simplex_fit(worst)) then
        state%simplex(:, worst) = trial; state%simplex_fit(worst) = fnew
      end if
      if (fnew <= state%simplex_fit(best) .and. used < maxeval) then
        centroid_sum = sum(state%simplex, dim=2)
        rfac = (1.0_dp - 2.0_dp) / real(n, dp)
        trial = centroid_sum * rfac - state%simplex(:, worst) * (rfac - 2.0_dp)
        trial = min(max(trial, lower), upper)
        fnew = fn(trial); used = used + 1
        if (fnew < state%simplex_fit(worst)) then
          state%simplex(:, worst) = trial; state%simplex_fit(worst) = fnew
        end if
      else if (fnew >= state%simplex_fit(second) .and. used < maxeval) then
        centroid_sum = sum(state%simplex, dim=2)
        rfac = (1.0_dp - 0.5_dp) / real(n, dp)
        trial = centroid_sum * rfac - state%simplex(:, worst) * (rfac - 0.5_dp)
        trial = min(max(trial, lower), upper)
        fnew = fn(trial); used = used + 1
        if (fnew < state%simplex_fit(worst)) then
          state%simplex(:, worst) = trial; state%simplex_fit(worst) = fnew
        end if
        do k = 1, n + 1
          if (k == best .or. used >= maxeval) cycle
          state%simplex(:, k) = 0.5_dp * (state%simplex(:, k) + state%simplex(:, best))
          state%simplex_fit(k) = fn(state%simplex(:, k)); used = used + 1
        end do
      end if
      best = minloc(state%simplex_fit, dim=1)
      if (state%simplex_fit(best) < fx) then
        fx = state%simplex_fit(best)
        x = state%simplex(:, best)
      end if
    end do
    best = minloc(state%simplex_fit, dim=1)
    if (state%simplex_fit(best) < fx) then
      fx = state%simplex_fit(best)
      x = state%simplex(:, best)
    end if
  end subroutine apply_simplex

  subroutine coordinate_fallback(x, fx, maxeval, lower, upper, fn, ctrl, used)
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval
    real(dp), intent(in) :: lower(:), upper(:)
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    integer :: i
    real(dp) :: old, fnew
    used = 0
    do i = 1, size(x)
      if (used >= maxeval .or. target_hit(fx, ctrl)) exit
      old = x(i)
      x(i) = min(max(old + 0.1_dp * (upper(i) - lower(i)), lower(i)), upper(i))
      fnew = fn(x); used = used + 1
      if (fnew < fx) then
        fx = fnew
      else
        x(i) = old
      end if
    end do
  end subroutine coordinate_fallback

  subroutine simplex_extremes(f, best, second, worst)
    real(dp), intent(in) :: f(:)
    integer, intent(out) :: best, second, worst
    integer :: i
    best = minloc(f, dim=1)
    worst = maxloc(f, dim=1)
    second = best
    do i = 1, size(f)
      if (i == worst) cycle
      if (second == worst .or. f(i) > f(second)) second = i
    end do
  end subroutine simplex_extremes

  subroutine apply_cmaes(state, x, fx, maxeval, pop, idx, lower, upper, rng, fn, ctrl, used)
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval, idx
    real(dp), intent(in) :: pop(:, :), lower(:), upper(:)
    type(rng_state), intent(inout) :: rng
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    integer :: n, lam, mu, i, j, k, info
    real(dp) :: mueff, cc, cs, mucov, ccov, damps, chin, psnorm
    logical :: hsig
    real(dp), allocatable :: stddev(:), arz(:, :), arx(:, :), fit(:), selz(:, :), selx(:, :)
    real(dp), allocatable :: zmean(:), bdz(:, :), rankmu(:, :), evals(:), evecs(:, :)
    integer, allocatable :: order(:)
    n = size(x)
    if (.not. state%initialized .or. state%kind /= 4) then
      call reset_ls_state(state)
      state%kind = 4
      state%initialized = .true.
      lam = ctrl%ls_param1
      if (lam <= 0) lam = 4 + floor(3.0_dp * log(real(n, dp)))
      mu = ctrl%ls_param2
      if (mu <= 0 .or. mu >= lam) mu = lam / 2
      state%lambda = lam; state%mu = mu
      allocate(state%xmean(n), state%pc(n), state%ps(n), state%c(n, n), state%b(n, n))
      allocate(state%d(n), state%bd(n, n), state%weights(mu), stddev(n))
      state%xmean = x; state%pc = 0.0_dp; state%ps = 0.0_dp
      call nearest_coordinate_distances(x, pop, stddev, idx)
      do i = 1, n
        stddev(i) = 0.5_dp * stddev(i) + 0.001_dp
        if (stddev(i) <= 0.001_dp) stddev(i) = 0.25_dp * (upper(i) - lower(i)) + 0.001_dp
      end do
      state%sigma = 1.0_dp
      state%c = 0.0_dp; state%b = 0.0_dp; state%bd = 0.0_dp
      do i = 1, n
        state%c(i, i) = stddev(i)**2
        state%b(i, i) = 1.0_dp
        state%d(i) = stddev(i)
        state%bd(i, i) = stddev(i)
      end do
      do i = 1, mu
        state%weights(i) = log(real(mu + 1, dp)) - log(real(i, dp))
      end do
      state%weights = state%weights / sum(state%weights)
    end if
    lam = state%lambda; mu = state%mu
    used = 0
    if (maxeval < lam) then
      call cma_small_step(state, x, fx, maxeval, lower, upper, rng, fn, ctrl, used)
      return
    end if
    allocate(arz(n, lam), arx(n, lam), fit(lam), order(lam))
    allocate(selz(n, mu), selx(n, mu), zmean(n), bdz(n, mu), rankmu(n, n))
    mueff = sum(state%weights)**2 / sum(state%weights**2)
    cc = 4.0_dp / real(n + 4, dp)
    cs = (mueff + 2.0_dp) / (real(n, dp) + mueff + 3.0_dp)
    mucov = mueff
    ccov = (1.0_dp / mucov) * 2.0_dp / (real(n, dp) + 1.4_dp)**2 + &
      (1.0_dp - 1.0_dp / mucov) * (2.0_dp * mucov - 1.0_dp) / ((real(n + 2, dp))**2 + 2.0_dp * mucov)
    damps = 1.0_dp + 2.0_dp * max(0.0_dp, sqrt((mueff - 1.0_dp) / real(n + 1, dp)) - 1.0_dp) + cs
    chin = sqrt(real(n, dp)) * (1.0_dp - 1.0_dp / (4.0_dp * real(n, dp)) + 1.0_dp / (21.0_dp * real(n*n, dp)))
    do while (used + lam <= maxeval .and. .not. target_hit(fx, ctrl))
      do k = 1, lam
        do i = 1, n
          arz(i, k) = rng_normal(rng)
        end do
      end do
      arx = matmul(state%bd, arz)
      do k = 1, lam
        arx(:, k) = min(max(state%xmean + state%sigma * arx(:, k), lower), upper)
        fit(k) = fn(arx(:, k))
        used = used + 1
        if (fit(k) < fx) then
          fx = fit(k); x = arx(:, k)
        end if
      end do
      call argsort(fit, order)
      do k = 1, mu
        selx(:, k) = arx(:, order(k)); selz(:, k) = arz(:, order(k))
      end do
      state%xmean = matmul(selx, state%weights)
      zmean = matmul(selz, state%weights)
      state%ps = (1.0_dp - cs) * state%ps + sqrt(cs * (2.0_dp - cs) * mueff) * matmul(state%b, zmean)
      psnorm = sqrt(dot_product(state%ps, state%ps))
      hsig = psnorm / sqrt(max(tiny(1.0_dp), 1.0_dp - (1.0_dp - cs)**(2.0_dp * &
        real(state%cma_evals + used, dp) / real(lam, dp)))) / chin < 1.4_dp + 2.0_dp / real(n + 1, dp)
      state%pc = (1.0_dp - cc) * state%pc
      if (hsig) state%pc = state%pc + sqrt(cc * (2.0_dp - cc) * mueff) * matmul(state%bd, zmean)
      bdz = matmul(state%bd, selz)
      rankmu = 0.0_dp
      do k = 1, mu
        do j = 1, n
          rankmu(:, j) = rankmu(:, j) + state%weights(k) * bdz(:, k) * bdz(j, k)
        end do
      end do
      state%c = (1.0_dp - ccov) * state%c + ccov / mucov * outer(state%pc, state%pc) + &
        ccov * (1.0_dp - 1.0_dp / mucov) * rankmu
      if (.not. hsig) state%c = state%c + ccov / mucov * cc * (2.0_dp - cc) * state%c
      state%sigma = state%sigma * exp((psnorm / chin - 1.0_dp) * cs / damps)
      call symmetric_eigen_descending(state%c, evals, evecs, info)
      if (info /= 0 .or. any(evals <= 0.0_dp)) exit
      state%b = evecs; state%d = sqrt(max(evals, tiny(1.0_dp)))
      do j = 1, n
        state%bd(:, j) = state%b(:, j) * state%d(j)
      end do
    end do
    state%cma_evals = state%cma_evals + used
    if (used < maxeval .and. .not. target_hit(fx, ctrl)) then
      i = maxeval - used
      call cma_small_step(state, x, fx, i, lower, upper, rng, fn, ctrl, k)
      used = used + k
    end if
  end subroutine apply_cmaes

  subroutine cma_small_step(state, x, fx, maxeval, lower, upper, rng, fn, ctrl, used)
    type(ls_state), intent(inout) :: state
    real(dp), intent(inout) :: x(:), fx
    integer, intent(in) :: maxeval
    real(dp), intent(in) :: lower(:), upper(:)
    type(rng_state), intent(inout) :: rng
    procedure(mals_objective) :: fn
    type(mals_control), intent(in) :: ctrl
    integer, intent(out) :: used
    real(dp), allocatable :: z(:), trial(:)
    real(dp) :: fnew
    integer :: i
    allocate(z(size(x)), trial(size(x)))
    used = 0
    do while (used < maxeval .and. .not. target_hit(fx, ctrl))
      do i = 1, size(x)
        z(i) = rng_normal(rng)
      end do
      trial = min(max(state%xmean + state%sigma * matmul(state%bd, z), lower), upper)
      fnew = fn(trial); used = used + 1
      if (fnew < fx) then
        fx = fnew; x = trial; state%xmean = trial
      end if
    end do
  end subroutine cma_small_step

  pure function outer(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: j
    do j = 1, size(y)
      a(:, j) = x * y(j)
    end do
  end function outer

  subroutine argsort(x, idx)
    real(dp), intent(in) :: x(:)
    integer, intent(out) :: idx(:)
    integer :: i, j, key
    do i = 1, size(x)
      idx(i) = i
    end do
    do i = 2, size(x)
      key = idx(i); j = i - 1
      do while (j >= 1)
        if (x(idx(j)) <= x(key)) exit
        idx(j + 1) = idx(j); j = j - 1
      end do
      idx(j + 1) = key
    end do
  end subroutine argsort
end module rmalschains_local_search
