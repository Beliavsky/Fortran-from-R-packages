module pcapp_pca
  use pcapp_kinds, only: dp
  use pcapp_types, only: scale_result, median_result, pca_result, covariance_result, tuning_result
  use pcapp_stats, only: column_center, column_scale, robust_scale, orient_columns
  use pcapp_stats, only: random_normal_matrix
  use pcapp_l1median, only: l1median_nlm
  implicit none
  private

  public :: scale_adv, pca_grid, spca_grid, pca_proj
  public :: cov_pc, cov_pca_grid, cov_pca_proj
  public :: opt_tpo, opt_bic, data_zou

contains

  pure function scale_adv(x, center_method, scale_method, center_values, scale_values) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix to center and/or scale.
    integer, intent(in), optional :: center_method !! Center selector: -1=none, 0=mean, 1=coordinate median, 2=spatial median.
    integer, intent(in), optional :: scale_method !! Scale selector: -1=none, 0=SD, 1=MAD, 2=Qn.
    real(dp), intent(in), optional :: center_values(:) !! Optional explicit center vector overriding center_method.
    real(dp), intent(in), optional :: scale_values(:) !! Optional explicit scale vector overriding scale_method.
    type(scale_result) :: res
    type(median_result) :: medres
    integer :: cm, sm, j

    allocate(res%x(size(x, 1), size(x, 2)))
    allocate(res%center(size(x, 2)), res%scale(size(x, 2)))
    cm = 0
    if (present(center_method)) cm = center_method
    sm = 0
    if (present(scale_method)) sm = scale_method

    if (present(center_values)) then
      if (size(center_values) == size(x, 2)) then
        res%center = center_values
      else
        res%center = 0.0_dp
      end if
    else if (cm == 2) then
      medres = l1median_nlm(x)
      res%center = medres%par
    else if (cm < 0) then
      res%center = 0.0_dp
    else
      call column_center(x, cm, res%center)
    end if

    res%x = x
    do j = 1, size(x, 2)
      res%x(:, j) = res%x(:, j) - res%center(j)
    end do

    if (present(scale_values)) then
      if (size(scale_values) == size(x, 2)) then
        res%scale = scale_values
      else
        res%scale = 1.0_dp
      end if
    else
      call column_scale(res%x, sm, res%scale)
    end if
    do j = 1, size(x, 2)
      if (abs(res%scale(j)) <= tiny(1.0_dp)) res%scale(j) = 1.0_dp
      res%x(:, j) = res%x(:, j) / res%scale(j)
    end do
  end function scale_adv

  pure subroutine descending_scale_order(y, method, order)
    real(dp), intent(in) :: y(:,:) !! Current residual coordinates in the active PCA subspace.
    integer, intent(in) :: method !! Robust scale selector used to rank coordinate directions.
    integer, intent(out) :: order(size(y, 2)) !! Column indices ordered from largest to smallest robust scale.
    real(dp), allocatable :: scales(:)
    real(dp) :: ts
    integer :: i, j, best, ti

    allocate(scales(size(y, 2)))
    do j = 1, size(y, 2)
      scales(j) = robust_scale(y(:, j), method)
      order(j) = j
    end do
    do i = 1, size(order) - 1
      best = i
      do j = i + 1, size(order)
        if (scales(j) > scales(best)) best = j
      end do
      if (best /= i) then
        ts = scales(i)
        scales(i) = scales(best)
        scales(best) = ts
        ti = order(i)
        order(i) = order(best)
        order(best) = ti
      end if
    end do
  end subroutine descending_scale_order

  pure subroutine householder_first(direction, h)
    real(dp), intent(in) :: direction(:) !! Unit vector to become the first column of an orthogonal basis.
    real(dp), intent(out) :: h(size(direction), size(direction)) !! Orthogonal Householder matrix with first column direction.
    real(dp), allocatable :: v(:)
    real(dp) :: norm2
    integer :: i, p

    p = size(direction)
    h = 0.0_dp
    do i = 1, p
      h(i, i) = 1.0_dp
    end do
    allocate(v(p))
    v = -direction
    v(1) = v(1) + 1.0_dp
    norm2 = dot_product(v, v)
    if (norm2 <= epsilon(1.0_dp)) return
    h = h - 2.0_dp * spread(v, 2, p) * spread(v, 1, p) / norm2
  end subroutine householder_first

  pure real(dp) function loading_penalty(load, norm_q, norm_s) result(value)
    real(dp), intent(in) :: load(:) !! Candidate loading vector in the original variable coordinates.
    real(dp), intent(in) :: norm_q !! Absolute-power exponent q in the sparse loading penalty.
    real(dp), intent(in) :: norm_s !! Outer exponent s in the sparse loading penalty.
    real(dp) :: base

    base = sum(abs(load)**norm_q)
    value = base**norm_s
  end function loading_penalty

  pure subroutine grid_component(y, basis, method, splitcircle, maxiter, zero_tol, lambda, norm_q, norm_s, &
      global_scale, direction, scatter, objective)
    real(dp), intent(in) :: y(:,:) !! Data represented in the current orthogonal complement.
    real(dp), intent(in) :: basis(:,:) !! Original-coordinate orthonormal basis for the current subspace.
    integer, intent(in) :: method !! Robust scale selector: 0=SD, 1=MAD, 2=Qn.
    integer, intent(in) :: splitcircle !! Number of candidate angles evaluated in each coordinate plane.
    integer, intent(in) :: maxiter !! Number of progressive grid refinements after the initial grid.
    real(dp), intent(in) :: zero_tol !! Numerical threshold for degenerate coordinate rotations.
    real(dp), intent(in) :: lambda !! Sparse loading penalty multiplier; zero gives ordinary PCAgrid.
    real(dp), intent(in) :: norm_q !! Absolute-power exponent q used by the sparse penalty.
    real(dp), intent(in) :: norm_s !! Outer exponent s used by the sparse penalty.
    real(dp), intent(in) :: global_scale !! Scatter normalization multiplying the sparse loading penalty.
    real(dp), intent(out) :: direction(size(y, 2)) !! Unit loading direction within the current subspace.
    real(dp), intent(out) :: scatter !! Robust scale of the selected projected scores.
    real(dp), intent(out) :: objective !! Selected grid objective value including any sparse penalty.
    integer, allocatable :: order(:)
    real(dp), allocatable :: a(:), best_a(:), base_a(:), proj(:), base_proj(:), cand_proj(:)
    real(dp), allocatable :: cand_a(:), original_load(:)
    real(dp) :: aj, comp, angle0, angle, c, s, cur_scatter, cur_obj
    real(dp) :: best_obj, best_scatter, pi, denom
    integer :: i, j, m, idx, psub, nang

    psub = size(y, 2)
    nang = max(3, splitcircle)
    allocate(order(psub), a(psub), best_a(psub), base_a(psub), cand_a(psub))
    allocate(proj(size(y, 1)), base_proj(size(y, 1)), cand_proj(size(y, 1)))
    allocate(original_load(size(basis, 1)))
    call descending_scale_order(y, method, order)
    a = 0.0_dp
    a(order(1)) = 1.0_dp
    proj = y(:, order(1))
    best_a = a
    best_scatter = robust_scale(proj, method)
    original_load = matmul(basis, a)
    best_obj = best_scatter**2 - lambda * loading_penalty(original_load, norm_q, norm_s) * global_scale**2
    pi = acos(-1.0_dp)

    do i = 0, maxiter
      do m = 1, psub
        idx = order(m)
        aj = max(-1.0_dp, min(1.0_dp, a(idx)))
        if (abs(abs(aj) - 1.0_dp) <= zero_tol) cycle
        comp = sqrt(max(0.0_dp, 1.0_dp - aj * aj))
        if (comp <= zero_tol) cycle
        base_proj = (proj - aj * y(:, idx)) / comp
        base_a = a / comp
        base_a(idx) = 0.0_dp
        angle0 = asin(aj)
        cur_obj = -huge(1.0_dp)
        cur_scatter = 0.0_dp
        cand_a = a
        do j = 1, nang
          denom = real(max(1, 2**i), dp)
          angle = angle0 + (-0.5_dp * pi + pi * real(j - 1, dp) / real(nang - 1, dp)) / denom
          c = cos(angle)
          s = sin(angle)
          cand_proj = c * base_proj + s * y(:, idx)
          original_load = matmul(basis, c * base_a)
          original_load = original_load + s * basis(:, idx)
          if (robust_scale(cand_proj, method)**2 - &
              lambda * loading_penalty(original_load, norm_q, norm_s) * global_scale**2 > cur_obj) then
            cur_scatter = robust_scale(cand_proj, method)
            cur_obj = cur_scatter**2 - &
              lambda * loading_penalty(original_load, norm_q, norm_s) * global_scale**2
            cand_a = c * base_a
            cand_a(idx) = s
          end if
        end do
        a = cand_a
        a = a / max(sqrt(sum(a**2)), tiny(1.0_dp))
        proj = matmul(y, a)
        if (cur_obj > best_obj) then
          best_obj = cur_obj
          best_scatter = cur_scatter
          best_a = a
        end if
      end do
    end do

    direction = best_a / max(sqrt(sum(best_a**2)), tiny(1.0_dp))
    scatter = best_scatter
    objective = best_obj
  end subroutine grid_component

  pure function grid_core(x, k, method, splitcircle, maxiter, zero_tol, lambda, norm_q, norm_s, &
      center_method, scale_method, sparse_mode) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix for grid projection pursuit.
    integer, intent(in) :: k !! Number of robust principal components to compute.
    integer, intent(in) :: method !! Robust projected-scale selector: 0=SD, 1=MAD, 2=Qn.
    integer, intent(in) :: splitcircle !! Number of candidate angles per coordinate-plane grid.
    integer, intent(in) :: maxiter !! Number of progressive grid refinements after the initial grid.
    real(dp), intent(in) :: zero_tol !! Numerical threshold for rotations and zero loadings.
    real(dp), intent(in) :: lambda(:) !! Sparse penalty multipliers, length one or at least k.
    real(dp), intent(in) :: norm_q !! Absolute-power exponent q in the sparse penalty.
    real(dp), intent(in) :: norm_s !! Outer exponent s in the sparse penalty.
    integer, intent(in) :: center_method !! Center selector passed to scale_adv.
    integer, intent(in) :: scale_method !! Scale selector passed to scale_adv; -1 disables scaling.
    logical, intent(in) :: sparse_mode !! True for sPCAgrid penalty behavior; false for PCAgrid.
    type(pca_result) :: res
    type(scale_result) :: scaled
    real(dp), allocatable :: basis(:,:), next_basis(:,:), y(:,:), h(:,:), a(:), load(:)
    real(dp), allocatable :: col_scales(:), lam(:)
    real(dp) :: global_scale, scat, obj
    integer :: comp, j, kk, p, psub

    p = size(x, 2)
    kk = min(max(k, 0), min(size(x, 1), p))
    scaled = scale_adv(x, center_method=center_method, scale_method=scale_method)
    allocate(res%loadings(p, kk), res%sdev(kk), res%objective(kk), res%lambda(kk))
    allocate(res%scores(size(x, 1), kk), res%center(p), res%scale(p))
    res%center = scaled%center
    res%scale = scaled%scale
    res%k = kk
    res%n_obs = size(x, 1)
    if (kk == 0) return

    allocate(basis(p, p), col_scales(p), lam(kk))
    basis = 0.0_dp
    do j = 1, p
      basis(j, j) = 1.0_dp
      col_scales(j) = robust_scale(scaled%x(:, j), method)
    end do
    global_scale = sqrt(sum(col_scales**2) / real(max(1, p), dp))
    if (.not. sparse_mode) global_scale = 1.0_dp
    if (size(lambda) == 1) then
      lam = lambda(1)
    else
      lam = lambda(1:kk)
    end if

    do comp = 1, kk
      psub = size(basis, 2)
      allocate(y(size(x, 1), psub), a(psub), load(p), h(psub, psub))
      y = matmul(scaled%x, basis)
      call grid_component(y, basis, method, splitcircle, maxiter, zero_tol, lam(comp), norm_q, norm_s, &
        global_scale, a, scat, obj)
      load = matmul(basis, a)
      res%loadings(:, comp) = load / max(sqrt(sum(load**2)), tiny(1.0_dp))
      res%sdev(comp) = scat
      res%objective(comp) = obj
      res%lambda(comp) = lam(comp)
      if (comp < kk) then
        call householder_first(a, h)
        allocate(next_basis(p, psub - 1))
        next_basis = matmul(basis, h(:, 2:psub))
        call move_alloc(next_basis, basis)
      end if
      deallocate(y, a, load, h)
    end do
    call orient_columns(res%loadings)
    res%scores = matmul(scaled%x, res%loadings)
    call sort_pca_desc(res)
  end function grid_core

  pure function pca_grid(x, k, method, maxiter, splitcircle, zero_tol, center_method, scale_method) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix for robust PCAgrid fitting.
    integer, intent(in), optional :: k !! Number of components; default min(2,n,p).
    integer, intent(in), optional :: method !! Robust scale selector: 0=SD, 1=MAD, 2=Qn; default MAD.
    integer, intent(in), optional :: maxiter !! Grid refinement count after the initial grid; default 10.
    integer, intent(in), optional :: splitcircle !! Candidate angles per grid; default 25.
    real(dp), intent(in), optional :: zero_tol !! Numerical zero threshold; default 1e-16.
    integer, intent(in), optional :: center_method !! Center selector; default spatial median (2).
    integer, intent(in), optional :: scale_method !! Scale selector; default -1 (no scaling), matching PCAgrid.
    type(pca_result) :: res
    integer :: kk, meth, mit, split, cm, sm
    real(dp) :: ztol

    kk = min(2, min(size(x, 1), size(x, 2)))
    if (present(k)) kk = k
    meth = 1
    if (present(method)) meth = method
    mit = 10
    if (present(maxiter)) mit = maxiter
    split = 25
    if (present(splitcircle)) split = splitcircle
    ztol = 1.0e-16_dp
    if (present(zero_tol)) ztol = zero_tol
    cm = 2
    if (present(center_method)) cm = center_method
    sm = -1
    if (present(scale_method)) sm = scale_method
    res = grid_core(x, kk, meth, split, mit, ztol, [0.0_dp], 1.0_dp, 1.0_dp, cm, sm, .false.)
  end function pca_grid

  pure function spca_grid(x, lambda, k, method, maxiter, splitcircle, zero_tol, norm_q, norm_s, &
      center_method, scale_method) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix for sparse robust PCAgrid fitting.
    real(dp), intent(in) :: lambda(:) !! Sparse penalty multiplier, supplied once or once per requested component.
    integer, intent(in), optional :: k !! Number of components; default min(2,n,p).
    integer, intent(in), optional :: method !! Robust scale selector: 0=SD, 1=MAD, 2=Qn; default MAD.
    integer, intent(in), optional :: maxiter !! Grid refinement count after the initial grid; default 10.
    integer, intent(in), optional :: splitcircle !! Candidate angles per grid; default 25.
    real(dp), intent(in), optional :: zero_tol !! Numerical zero threshold; default 1e-16.
    real(dp), intent(in), optional :: norm_q !! Absolute-power exponent q for the loading penalty; default 1.
    real(dp), intent(in), optional :: norm_s !! Outer exponent s for the loading penalty; default 1.
    integer, intent(in), optional :: center_method !! Center selector; default spatial median (2).
    integer, intent(in), optional :: scale_method !! Scale selector; default -1 (no scaling).
    type(pca_result) :: res
    integer :: kk, meth, mit, split, cm, sm
    real(dp) :: ztol, nq, ns

    kk = min(2, min(size(x, 1), size(x, 2)))
    if (present(k)) kk = k
    meth = 1
    if (present(method)) meth = method
    mit = 10
    if (present(maxiter)) mit = maxiter
    split = 25
    if (present(splitcircle)) split = splitcircle
    ztol = 1.0e-16_dp
    if (present(zero_tol)) ztol = zero_tol
    nq = 1.0_dp
    if (present(norm_q)) nq = norm_q
    ns = 1.0_dp
    if (present(norm_s)) ns = norm_s
    cm = 2
    if (present(center_method)) cm = center_method
    sm = -1
    if (present(scale_method)) sm = scale_method
    res = grid_core(x, kk, meth, split, mit, ztol, lambda, nq, ns, cm, sm, .true.)
  end function spca_grid

  subroutine append_projection_candidates(work, calc_method, nmax, candidates)
    real(dp), intent(in) :: work(:,:) !! Current deflated observations used to construct projection candidates.
    integer, intent(in) :: calc_method !! Candidate mode: 0=observations, 1=random linear combinations, 2=random sphere.
    integer, intent(in) :: nmax !! Maximum candidate count for random augmentation modes.
    real(dp), allocatable, intent(out) :: candidates(:,:) !! Candidate vectors before row normalization.
    real(dp), allocatable :: extra(:,:), weights(:,:)
    integer :: n, p, target

    n = size(work, 1)
    p = size(work, 2)
    target = n
    if (calc_method /= 0) target = max(n, nmax)
    allocate(candidates(target, p))
    candidates(1:n, :) = work
    if (target == n) return
    allocate(extra(target - n, p))
    if (calc_method == 1) then
      allocate(weights(target - n, n))
      call random_number(weights)
      extra = matmul(weights, work)
    else
      call random_normal_matrix(extra)
    end if
    candidates(n + 1:target, :) = extra
  end subroutine append_projection_candidates

  subroutine normalize_rows(a, zero_tol, keep)
    real(dp), intent(inout) :: a(:,:) !! Candidate rows normalized in place when their norm exceeds zero_tol.
    real(dp), intent(in) :: zero_tol !! Row-norm threshold used to mark degenerate candidates.
    logical, intent(out) :: keep(size(a, 1)) !! True for candidate rows with usable nonzero norms.
    real(dp) :: normv
    integer :: i

    do i = 1, size(a, 1)
      normv = sqrt(sum(a(i, :)**2))
      keep(i) = normv > zero_tol
      if (keep(i)) a(i, :) = a(i, :) / normv
    end do
  end subroutine normalize_rows

  pure subroutine orthogonalize_direction(v, previous, zero_tol)
    real(dp), intent(inout) :: v(:) !! Candidate loading normalized after removing previous loading components.
    real(dp), intent(in) :: previous(:,:) !! Previously selected loading columns; may have zero columns.
    real(dp), intent(in) :: zero_tol !! Norm threshold below which normalization is skipped.
    real(dp) :: normv
    integer :: j

    do j = 1, size(previous, 2)
      v = v - dot_product(v, previous(:, j)) * previous(:, j)
    end do
    normv = sqrt(sum(v**2))
    if (normv > zero_tol) v = v / normv
  end subroutine orthogonalize_direction

  subroutine refine_projection(work, direction, method, maxit, maxhalf, zero_tol, value)
    real(dp), intent(in) :: work(:,:) !! Current deflated observation matrix.
    real(dp), intent(inout) :: direction(:) !! Candidate loading direction refined in place.
    integer, intent(in) :: method !! Robust projected-scale selector.
    integer, intent(in) :: maxit !! Maximum sign-update iterations.
    integer, intent(in) :: maxhalf !! Maximum step halvings when a sign update reduces the objective.
    real(dp), intent(in) :: zero_tol !! Numerical threshold for candidate normalization.
    real(dp), intent(inout) :: value !! Current projected scale, updated with accepted refinements.
    real(dp), allocatable :: score(:), new_dir(:), old_dir(:), trial(:)
    real(dp) :: normv, new_value
    integer :: i, j, half

    allocate(score(size(work, 1)), new_dir(size(direction)), old_dir(size(direction)), trial(size(direction)))
    score = matmul(work, direction)
    do i = 1, maxit
      new_dir = 0.0_dp
      do j = 1, size(work, 1)
        if (abs(score(j)) <= zero_tol) cycle
        new_dir = new_dir + sign(1.0_dp, score(j)) * work(j, :) / &
          max(sqrt(sum(work(j, :)**2)), zero_tol)
      end do
      normv = sqrt(sum(new_dir**2))
      if (normv <= zero_tol) exit
      new_dir = new_dir / normv
      old_dir = direction
      new_value = value
      do half = 0, maxhalf
        if (half == 0) then
          trial = new_dir
        else
          trial = 0.5_dp * (trial + old_dir)
          normv = sqrt(sum(trial**2))
          if (normv > zero_tol) trial = trial / normv
        end if
        new_value = robust_scale(matmul(work, trial), method)
        if (new_value >= value) exit
      end do
      if (new_value < value) exit
      direction = trial
      value = new_value
      score = matmul(work, direction)
    end do
  end subroutine refine_projection

  function pca_proj(x, k, method, calc_method, nmax, update, maxit, maxhalf, zero_tol, &
      center_method, scale_method) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for projection-pursuit robust PCA.
    integer, intent(in), optional :: k !! Number of components; default min(2,n,p).
    integer, intent(in), optional :: method !! Robust projected-scale selector: 0=SD, 1=MAD, 2=Qn; default MAD.
    integer, intent(in), optional :: calc_method !! Candidate mode: 0=each observation, 1=lincomb, 2=sphere; default 0.
    integer, intent(in), optional :: nmax !! Maximum random candidate count for modes 1 and 2; default 1000.
    logical, intent(in), optional :: update !! Whether to use iterative sign refinement; default true.
    integer, intent(in), optional :: maxit !! Maximum sign-refinement iterations; default 5.
    integer, intent(in), optional :: maxhalf !! Maximum step halvings per update; default 5.
    real(dp), intent(in), optional :: zero_tol !! Numerical zero threshold; default 1e-16.
    integer, intent(in), optional :: center_method !! Center selector; default spatial median (2).
    integer, intent(in), optional :: scale_method !! Scale selector; default -1 (no scaling).
    type(pca_result) :: res
    type(scale_result) :: scaled
    real(dp), allocatable :: work(:,:), candidates(:,:), direction(:), scores(:)
    logical, allocatable :: keep(:)
    real(dp) :: value, cur, ztol, normv
    integer :: kk, meth, mode, max_candidates, mit, mhalf, cm, sm
    integer :: comp, i, best, p
    logical :: do_update

    p = size(x, 2)
    kk = min(2, min(size(x, 1), p))
    if (present(k)) kk = min(k, min(size(x, 1), p))
    meth = 1
    if (present(method)) meth = method
    mode = 0
    if (present(calc_method)) mode = calc_method
    max_candidates = 1000
    if (present(nmax)) max_candidates = nmax
    do_update = .true.
    if (present(update)) do_update = update
    mit = 5
    if (present(maxit)) mit = maxit
    mhalf = 5
    if (present(maxhalf)) mhalf = maxhalf
    ztol = 1.0e-16_dp
    if (present(zero_tol)) ztol = zero_tol
    cm = 2
    if (present(center_method)) cm = center_method
    sm = -1
    if (present(scale_method)) sm = scale_method

    scaled = scale_adv(x, center_method=cm, scale_method=sm)
    work = scaled%x
    allocate(res%loadings(p, kk), res%sdev(kk), res%objective(kk), res%lambda(kk))
    allocate(res%scores(size(x, 1), kk), res%center(p), res%scale(p))
    res%loadings = 0.0_dp
    res%sdev = 0.0_dp
    res%objective = 0.0_dp
    res%lambda = 0.0_dp
    res%center = scaled%center
    res%scale = scaled%scale
    res%k = kk
    res%n_obs = size(x, 1)

    do comp = 1, kk
      call append_projection_candidates(work, mode, max_candidates, candidates)
      allocate(keep(size(candidates, 1)))
      call normalize_rows(candidates, ztol, keep)
      allocate(direction(p), scores(size(x, 1)))
      best = 0
      value = -huge(1.0_dp)
      do i = 1, size(candidates, 1)
        if (.not. keep(i)) cycle
        direction = candidates(i, :)
        if (comp > 1) call orthogonalize_direction(direction, res%loadings(:, 1:comp - 1), ztol)
        normv = sqrt(sum(direction**2))
        if (normv <= ztol) cycle
        direction = direction / normv
        scores = matmul(work, direction)
        cur = robust_scale(scores, meth)
        if (cur > value) then
          value = cur
          best = i
          res%loadings(:, comp) = direction
        end if
      end do
      if (best == 0) exit
      direction = res%loadings(:, comp)
      if (do_update) call refine_projection(work, direction, meth, mit, mhalf, ztol, value)
      if (comp > 1) call orthogonalize_direction(direction, res%loadings(:, 1:comp - 1), ztol)
      res%loadings(:, comp) = direction
      res%sdev(comp) = robust_scale(matmul(work, direction), meth)
      res%objective(comp) = res%sdev(comp)
      scores = matmul(work, direction)
      work = work - spread(scores, 2, p) * spread(direction, 1, size(x, 1))
      deallocate(candidates, keep, direction, scores)
    end do
    call orient_columns(res%loadings)
    res%scores = matmul(scaled%x, res%loadings)
    call sort_pca_desc(res)
  end function pca_proj

  pure subroutine sort_pca_desc(res)
    type(pca_result), intent(inout) :: res !! PCA result whose component fields are reordered by decreasing sdev.
    real(dp), allocatable :: temp_load(:), temp_score(:)
    real(dp) :: ts
    integer :: i, j, best

    if (.not. allocated(res%sdev)) return
    allocate(temp_load(size(res%loadings, 1)), temp_score(size(res%scores, 1)))
    do i = 1, size(res%sdev) - 1
      best = i
      do j = i + 1, size(res%sdev)
        if (res%sdev(j) > res%sdev(best)) best = j
      end do
      if (best /= i) then
        ts = res%sdev(i)
        res%sdev(i) = res%sdev(best)
        res%sdev(best) = ts
        temp_load = res%loadings(:, i)
        res%loadings(:, i) = res%loadings(:, best)
        res%loadings(:, best) = temp_load
        temp_score = res%scores(:, i)
        res%scores(:, i) = res%scores(:, best)
        res%scores(:, best) = temp_score
        ts = res%objective(i)
        res%objective(i) = res%objective(best)
        res%objective(best) = ts
        ts = res%lambda(i)
        res%lambda(i) = res%lambda(best)
        res%lambda(best) = ts
      end if
    end do
  end subroutine sort_pca_desc

  pure function cov_pc(pc, k) result(res)
    type(pca_result), intent(in) :: pc !! PCA result providing loadings, component scales, and center.
    integer, intent(in), optional :: k !! Number of leading components used in covariance reconstruction.
    type(covariance_result) :: res
    real(dp), allocatable :: weighted(:,:)
    integer :: i, kk, p

    p = size(pc%loadings, 1)
    kk = size(pc%loadings, 2)
    if (present(k)) kk = min(k, kk)
    allocate(weighted(p, kk), res%covariance(p, p), res%center(p))
    weighted = pc%loadings(:, 1:kk)
    do i = 1, kk
      weighted(:, i) = weighted(:, i) * pc%sdev(i)**2
    end do
    res%covariance = matmul(weighted, transpose(pc%loadings(:, 1:kk)))
    res%center = pc%center
  end function cov_pc

  pure function cov_pca_grid(x, method) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for full-rank PCAgrid covariance estimation.
    integer, intent(in), optional :: method !! Robust projected-scale selector; default MAD.
    type(covariance_result) :: res
    type(pca_result) :: pc
    integer :: meth

    meth = 1
    if (present(method)) meth = method
    pc = pca_grid(x, k=min(size(x, 1), size(x, 2)), method=meth)
    res = cov_pc(pc)
  end function cov_pca_grid

  function cov_pca_proj(x, method) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for full-rank PCAproj covariance estimation.
    integer, intent(in), optional :: method !! Robust projected-scale selector; default MAD.
    type(covariance_result) :: res
    type(pca_result) :: pc
    integer :: meth

    meth = 1
    if (present(method)) meth = method
    pc = pca_proj(x, k=min(size(x, 1), size(x, 2)), method=meth)
    res = cov_pc(pc)
  end function cov_pca_proj

  pure real(dp) function residual_scale_sum(x, pc, k, method) result(value)
    real(dp), intent(in) :: x(:,:) !! Original data matrix used for sparse-PCA tuning.
    type(pca_result), intent(in) :: pc !! Candidate sparse PCA fit.
    integer, intent(in) :: k !! Number of leading components included in reconstruction.
    integer, intent(in) :: method !! Robust column-scale selector for residual scatter.
    real(dp), allocatable :: centered(:,:), residual(:,:), proj(:,:)
    integer :: j, kk

    kk = min(k, size(pc%loadings, 2))
    centered = x - spread(pc%center, 1, size(x, 1))
    proj = matmul(centered, pc%loadings(:, 1:kk))
    residual = centered - matmul(proj, transpose(pc%loadings(:, 1:kk)))
    value = 0.0_dp
    do j = 1, size(residual, 2)
      value = value + robust_scale(residual(:, j), method)**2
    end do
  end function residual_scale_sum

  pure function opt_tpo(x, k_max, n_lambda, lambda_max, method) result(res)
    real(dp), intent(in) :: x(:,:) !! Data matrix for sparse-PCA tradeoff tuning.
    integer, intent(in), optional :: k_max !! Number of sparse components optimized; default ncol(x) capped by nrow(x).
    integer, intent(in), optional :: n_lambda !! Number of equally spaced penalty values; default 30.
    real(dp), intent(in), optional :: lambda_max !! Upper end of penalty grid; default 2.
    integer, intent(in), optional :: method !! Robust projected-scale selector; default MAD.
    type(tuning_result) :: res
    type(pca_result), allocatable :: fits(:)
    real(dp), allocatable :: sparse(:), variance(:)
    real(dp) :: lmax, sxden, vyden
    integer :: i, j, kk, nl, meth, zeros

    kk = min(size(x, 2), size(x, 1))
    if (present(k_max)) kk = min(k_max, kk)
    nl = 30
    if (present(n_lambda)) nl = max(2, n_lambda)
    lmax = 2.0_dp
    if (present(lambda_max)) lmax = max(0.0_dp, lambda_max)
    meth = 1
    if (present(method)) meth = method
    allocate(res%lambda_grid(nl), res%criterion(nl), fits(nl), sparse(nl), variance(nl))
    do i = 1, nl
      res%lambda_grid(i) = lmax * real(i - 1, dp) / real(nl - 1, dp)
      fits(i) = spca_grid(x, [res%lambda_grid(i)], k=kk, method=meth)
      zeros = 0
      do j = 1, size(fits(i)%loadings, 2)
        zeros = zeros + count(abs(fits(i)%loadings(:, j)) <= 1.0e-12_dp)
      end do
      sparse(i) = real(zeros, dp)
      variance(i) = sum(fits(i)%sdev**2)
    end do
    sxden = max(abs(sparse(nl) - sparse(1)), 1.0_dp)
    vyden = max(abs(variance(1) - variance(nl)), epsilon(1.0_dp))
    do i = 1, nl
      res%criterion(i) = -((sparse(i) - sparse(1)) / sxden) * ((variance(i) - variance(nl)) / vyden)
    end do
    res%best_index = minloc(res%criterion, dim=1)
    res%pc = fits(res%best_index)
  end function opt_tpo

  pure function opt_bic(x, k_max, n_lambda, lambda_max, method) result(res)
    real(dp), intent(in) :: x(:,:) !! Data matrix for BIC-style sparse-PCA tuning.
    integer, intent(in), optional :: k_max !! Number of sparse components optimized; default ncol(x) capped by nrow(x).
    integer, intent(in), optional :: n_lambda !! Number of equally spaced penalty values; default 30.
    real(dp), intent(in), optional :: lambda_max !! Upper end of penalty grid; default 2.
    integer, intent(in), optional :: method !! Robust projected-scale selector; default MAD.
    type(tuning_result) :: res
    type(pca_result), allocatable :: fits(:)
    real(dp) :: lmax, rss
    integer :: i, j, kk, nl, meth, df

    kk = min(size(x, 2), size(x, 1))
    if (present(k_max)) kk = min(k_max, kk)
    nl = 30
    if (present(n_lambda)) nl = max(2, n_lambda)
    lmax = 2.0_dp
    if (present(lambda_max)) lmax = max(0.0_dp, lambda_max)
    meth = 1
    if (present(method)) meth = method
    allocate(res%lambda_grid(nl), res%criterion(nl), fits(nl))
    do i = 1, nl
      res%lambda_grid(i) = lmax * real(i - 1, dp) / real(nl - 1, dp)
      fits(i) = spca_grid(x, [res%lambda_grid(i)], k=kk, method=meth)
      rss = residual_scale_sum(x, fits(i), kk, meth)
      df = 0
      do j = 1, size(fits(i)%loadings, 2)
        df = df + count(abs(fits(i)%loadings(:, j)) > 1.0e-12_dp)
      end do
      res%criterion(i) = real(size(x, 1), dp) * log(max(rss, tiny(1.0_dp))) + &
        real(df, dp) * log(real(max(2, size(x, 1)), dp))
    end do
    res%best_index = minloc(res%criterion, dim=1)
    res%pc = fits(res%best_index)
  end function opt_bic

  subroutine data_zou(n, p, x)
    integer, intent(in) :: n !! Number of simulated observations; pcaPP defaults to 250.
    integer, intent(in) :: p(3) !! Numbers of replicated variables from the three latent factors.
    real(dp), allocatable, intent(out) :: x(:,:) !! Simulated Zou-style data matrix with n rows and sum(p) columns.
    real(dp), allocatable :: v(:,:), err(:,:)
    integer :: j, col

    allocate(v(n, 3), x(n, sum(p)), err(n, sum(p)))
    call random_normal_matrix(v)
    v(:, 1) = sqrt(290.0_dp) * v(:, 1)
    v(:, 2) = sqrt(300.0_dp) * v(:, 2)
    v(:, 3) = -0.3_dp * v(:, 1) + 0.925_dp * v(:, 2) + v(:, 3)
    col = 0
    do j = 1, p(1)
      col = col + 1
      x(:, col) = v(:, 1)
    end do
    do j = 1, p(2)
      col = col + 1
      x(:, col) = v(:, 2)
    end do
    do j = 1, p(3)
      col = col + 1
      x(:, col) = v(:, 3)
    end do
    call random_normal_matrix(err)
    x = x + err
  end subroutine data_zou

end module pcapp_pca
