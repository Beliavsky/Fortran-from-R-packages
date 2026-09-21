module ecp_cp3o
  use ecp_kinds, only: dp
  use ecp_types, only: cp_result
  use ecp_energy, only: distance_matrix_alpha
  use ecp_energy, only: ks_split_stat
  implicit none
  private

  public :: e_cp3o, e_cp3o_delta, ks_cp3o, ks_cp3o_delta

contains

  pure function e_cp3o(z, k, minsize, alpha, prune) result(res)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    integer, intent(in), optional :: k !! Maximum number of change points; default is one.
    integer, intent(in), optional :: minsize !! Minimum segment size; default is 30.
    real(dp), intent(in), optional :: alpha !! Energy-distance exponent in (0, 2]; default is one.
    logical, intent(in), optional :: prune !! Use CP3O candidate pruning; default is true.
    type(cp_result) :: res
    integer :: kk, mm
    real(dp) :: aa
    logical :: pp

    kk = 1
    if (present(k)) kk = k
    mm = 30
    if (present(minsize)) mm = minsize
    aa = 1.0_dp
    if (present(alpha)) aa = alpha
    pp = .true.
    if (present(prune)) pp = prune
    res = cp3o_engine(z, kk, mm, aa, 1, mm - 1, pp)
  end function e_cp3o

  pure function e_cp3o_delta(z, k, delta, alpha, prune) result(res)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    integer, intent(in), optional :: k !! Maximum number of change points; default is one.
    integer, intent(in), optional :: delta !! Window size for the complete statistic portion; default is 29.
    real(dp), intent(in), optional :: alpha !! Energy-distance exponent in (0, 2]; default is one.
    logical, intent(in), optional :: prune !! Use CP3O candidate pruning; default is true.
    type(cp_result) :: res
    integer :: kk, dd
    real(dp) :: aa
    logical :: pp

    kk = 1
    if (present(k)) kk = k
    dd = 29
    if (present(delta)) dd = delta
    aa = 1.0_dp
    if (present(alpha)) aa = alpha
    pp = .true.
    if (present(prune)) pp = prune
    res = cp3o_engine(z, kk, dd + 1, aa, 2, dd, pp)
  end function e_cp3o_delta

  pure function ks_cp3o(z, k, minsize, prune) result(res)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    integer, intent(in), optional :: k !! Maximum number of change points; default is one.
    integer, intent(in), optional :: minsize !! Minimum segment size; default is 30.
    logical, intent(in), optional :: prune !! Use CP3O candidate pruning; default is true.
    type(cp_result) :: res
    integer :: kk, mm
    logical :: pp

    kk = 1
    if (present(k)) kk = k
    mm = 30
    if (present(minsize)) mm = minsize
    pp = .true.
    if (present(prune)) pp = prune
    res = cp3o_engine(z, kk, mm, 1.0_dp, 3, mm, pp)
  end function ks_cp3o

  pure function ks_cp3o_delta(z, k, minsize, prune) result(res)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    integer, intent(in), optional :: k !! Maximum number of change points; default is one.
    integer, intent(in), optional :: minsize !! Minimum segment and comparison-window size; default is 30.
    logical, intent(in), optional :: prune !! Use CP3O candidate pruning; default is true.
    type(cp_result) :: res
    integer :: kk, mm
    logical :: pp

    kk = 1
    if (present(k)) kk = k
    mm = 30
    if (present(minsize)) mm = minsize
    pp = .true.
    if (present(prune)) pp = prune
    res = cp3o_engine(z, kk, mm, 1.0_dp, 4, mm, pp)
  end function ks_cp3o_delta

  pure function cp3o_engine(z, k, minsize, alpha, mode, window, prune) result(res)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    integer, intent(in) :: k !! Requested maximum number of change points.
    integer, intent(in) :: minsize !! Minimum number of observations in every segment.
    real(dp), intent(in) :: alpha !! Energy-distance exponent, used by energy modes.
    integer, intent(in) :: mode !! Statistic selector: 1 energy, 2 windowed energy, 3 KS, 4 windowed KS.
    integer, intent(in) :: window !! Window size used by the approximate modes.
    logical, intent(in) :: prune !! Apply upstream-style candidate pruning when true.
    type(cp_result) :: res
    real(dp), allocatable :: d(:,:), prefix(:,:), chain(:), f(:,:), stat_cache(:,:), stat_cache_new(:,:)
    integer, allocatable :: a(:,:)
    logical, allocatable :: active(:,:), active_new(:,:)
    integer :: n, k_eff, level, s, selected, evals, pruned
    real(dp) :: neg, best_sse, sse

    n = size(z, 1)
    if (n < 2*minsize .or. k < 1) then
      allocate(res%estimates(0), res%gof(0), res%cp_loc(0, 0))
      return
    end if
    k_eff = min(k, max(1, n/minsize - 1))
    neg = -huge(1.0_dp)
    allocate(f(k_eff, n), source=neg)
    allocate(a(k_eff, n), source=0)
    evals = 0
    pruned = 0

    if (mode <= 2) then
      d = distance_matrix_alpha(z, alpha)
      call build_distance_prefix(d, prefix, chain)
    else
      allocate(d(0, 0), prefix(0, 0), chain(0))
    end if

    if (.not. prune .or. k_eff == 1) then
      call exhaustive_dp(z, prefix, chain, f, a, minsize, mode, window, evals)
    else
      allocate(active(n, n), source=.false.)
      allocate(active_new(n, n), source=.false.)
      allocate(stat_cache(n, n), source=0.0_dp)
      allocate(stat_cache_new(n, n), source=0.0_dp)

      call first_level(z, prefix, chain, f, a, minsize, mode, window, evals)
      call prune_level(z, prefix, chain, f, a, 1, minsize, mode, window, active, stat_cache, evals, pruned)

      do level = 2, k_eff
        if (level == k_eff) then
          call final_pruned_level(f, a, level, minsize, active, stat_cache)
        else
          call next_pruned_level(f, a, level, minsize, active, stat_cache)
          active_new = .false.
          stat_cache_new = 0.0_dp
          call prune_level(z, prefix, chain, f, a, level, minsize, mode, window, active_new, stat_cache_new, &
                           evals, pruned, active)
          active = active_new
          stat_cache = stat_cache_new
        end if
      end do
    end if

    res%statistic_evaluations = evals
    res%candidates_pruned = pruned
    allocate(res%gof(k_eff))
    res%gof = f(:, n)
    allocate(res%cp_loc(k_eff, k_eff), source=0)
    do level = 1, k_eff
      call reconstruct_row(a, level, n, res%cp_loc(level, :))
    end do

    if (k_eff == 1) then
      selected = 1
    else if (k_eff == 2) then
      selected = 1
    else
      selected = 2
      best_sse = huge(1.0_dp)
      do s = 2, k_eff - 1
        sse = piecewise_sse(res%gof, s)
        if (sse < best_sse) then
          best_sse = sse
          selected = s
        end if
      end do
    end if
    if (a(selected, n) <= 0) then
      selected = 1
      do level = 1, k_eff
        if (a(level, n) > 0) selected = level
      end do
    end if
    res%number = selected
    allocate(res%estimates(selected))
    res%estimates = res%cp_loc(selected, 1:selected)
  end function cp3o_engine

  pure subroutine exhaustive_dp(z, prefix, chain, f, a, minsize, mode, window, evals)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix sums of pairwise distances, or zero-size for KS.
    real(dp), intent(in) :: chain(0:) !! Prefix sums of adjacent distances, or zero-size for KS.
    real(dp), intent(inout) :: f(:,:) !! Dynamic-programming goodness-of-fit table.
    integer, intent(inout) :: a(:,:) !! Dynamic-programming last-change table.
    integer, intent(in) :: minsize !! Minimum number of observations in every segment.
    integer, intent(in) :: mode !! Statistic selector used by the CP3O engine.
    integer, intent(in) :: window !! Window size used by approximate modes.
    integer, intent(inout) :: evals !! Number of local statistic evaluations performed.
    integer :: level, t, s, u
    real(dp) :: stat, candidate

    do t = 2*minsize, size(z, 1)
      do s = minsize, t - minsize
        stat = local_stat_fast(z, prefix, chain, 1, s, t, mode, window)
        evals = evals + 1
        if (stat > f(1, t)) then
          f(1, t) = stat
          a(1, t) = s
        end if
      end do
    end do

    do level = 2, size(f, 1)
      do t = (level + 1)*minsize, size(z, 1)
        do s = level*minsize, t - minsize
          if (a(level - 1, s) <= 0) cycle
          if (.not. finite_score(f(level - 1, s))) cycle
          u = a(level - 1, s)
          stat = local_stat_fast(z, prefix, chain, u + 1, s, t, mode, window)
          evals = evals + 1
          candidate = f(level - 1, s) + stat
          if (candidate > f(level, t)) then
            f(level, t) = candidate
            a(level, t) = s
          end if
        end do
      end do
    end do
  end subroutine exhaustive_dp

  pure subroutine first_level(z, prefix, chain, f, a, minsize, mode, window, evals)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix sums of pairwise distances, or zero-size for KS.
    real(dp), intent(in) :: chain(0:) !! Prefix sums of adjacent distances, or zero-size for KS.
    real(dp), intent(inout) :: f(:,:) !! Dynamic-programming goodness-of-fit table.
    integer, intent(inout) :: a(:,:) !! Dynamic-programming last-change table.
    integer, intent(in) :: minsize !! Minimum number of observations in every segment.
    integer, intent(in) :: mode !! Statistic selector used by the CP3O engine.
    integer, intent(in) :: window !! Window size used by approximate modes.
    integer, intent(inout) :: evals !! Number of local statistic evaluations performed.
    integer :: t, s
    real(dp) :: stat

    do t = 2*minsize, size(z, 1)
      do s = minsize, t - minsize
        stat = local_stat_fast(z, prefix, chain, 1, s, t, mode, window)
        evals = evals + 1
        if (stat > f(1, t)) then
          f(1, t) = stat
          a(1, t) = s
        end if
      end do
    end do
  end subroutine first_level

  pure subroutine next_pruned_level(f, a, level, minsize, active, stat_cache)
    real(dp), intent(inout) :: f(:,:) !! Dynamic-programming goodness-of-fit table.
    integer, intent(inout) :: a(:,:) !! Dynamic-programming last-change table.
    integer, intent(in) :: level !! Number of change points represented by the row being computed.
    integer, intent(in) :: minsize !! Minimum number of observations in every segment.
    logical, intent(in) :: active(:,:) !! Candidate mask retained by the preceding pruning level.
    real(dp), intent(in) :: stat_cache(:,:) !! Cached local statistics corresponding to active candidate states.
    integer :: t, s
    real(dp) :: candidate

    do t = (level + 1)*minsize, size(f, 2)
      do s = level*minsize, t - minsize
        if (.not. active(t, s)) cycle
        if (.not. finite_score(f(level - 1, s))) cycle
        candidate = f(level - 1, s) + stat_cache(t, s)
        if (candidate > f(level, t)) then
          f(level, t) = candidate
          a(level, t) = s
        end if
      end do
    end do
  end subroutine next_pruned_level

  pure subroutine final_pruned_level(f, a, level, minsize, active, stat_cache)
    real(dp), intent(inout) :: f(:,:) !! Dynamic-programming goodness-of-fit table.
    integer, intent(inout) :: a(:,:) !! Dynamic-programming last-change table.
    integer, intent(in) :: level !! Number of change points represented by the final row.
    integer, intent(in) :: minsize !! Minimum number of observations in every segment.
    logical, intent(in) :: active(:,:) !! Candidate mask retained by the preceding pruning level.
    real(dp), intent(in) :: stat_cache(:,:) !! Cached local statistics corresponding to active candidate states.
    integer :: n, s
    real(dp) :: candidate

    n = size(f, 2)
    do s = level*minsize, n - minsize
      if (.not. active(n, s)) cycle
      if (.not. finite_score(f(level - 1, s))) cycle
      candidate = f(level - 1, s) + stat_cache(n, s)
      if (candidate > f(level, n)) then
        f(level, n) = candidate
        a(level, n) = s
      end if
    end do
  end subroutine final_pruned_level

  pure subroutine prune_level(z, prefix, chain, f, a, level, minsize, mode, window, active_out, stat_cache_out, &
                               evals, pruned, active_in)
    real(dp), intent(in) :: z(:,:) !! Time series matrix with observations in rows.
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix sums of pairwise distances, or zero-size for KS.
    real(dp), intent(in) :: chain(0:) !! Prefix sums of adjacent distances, or zero-size for KS.
    real(dp), intent(in) :: f(:,:) !! Dynamic-programming goodness-of-fit table for completed levels.
    integer, intent(in) :: a(:,:) !! Dynamic-programming last-change table for completed levels.
    integer, intent(in) :: level !! Number of change points represented by the row being pruned.
    integer, intent(in) :: minsize !! Minimum number of observations in every segment.
    integer, intent(in) :: mode !! Statistic selector used by the CP3O engine.
    integer, intent(in) :: window !! Window size used by approximate modes.
    logical, intent(out) :: active_out(:,:) !! Candidate mask retained for the next level.
    real(dp), intent(out) :: stat_cache_out(:,:) !! Local statistics cached for retained candidate states.
    integer, intent(inout) :: evals !! Number of local statistic evaluations performed.
    integer, intent(inout) :: pruned !! Number of candidate states removed by pruning.
    logical, intent(in), optional :: active_in(:,:) !! Candidate mask from the preceding pruning level.
    integer :: t, s, b, a2, b2, start_t
    real(dp) :: ref_score, score, stat
    logical :: eligible

    active_out = .false.
    stat_cache_out = 0.0_dp
    start_t = (level + 2)*minsize
    do t = start_t, size(z, 1)
      a2 = t - minsize
      if (a(level, a2) <= 0 .or. .not. finite_score(f(level, a2))) cycle
      b2 = a(level, a2)
      stat = local_stat_fast(z, prefix, chain, b2 + 1, a2, t, mode, window)
      evals = evals + 1
      ref_score = f(level, a2) + stat

      do s = level*minsize, t - minsize
        eligible = .true.
        if (present(active_in)) eligible = active_in(t, s)
        if (.not. eligible) cycle
        if (a(level, s) <= 0 .or. .not. finite_score(f(level, s))) then
          pruned = pruned + 1
          cycle
        end if
        if (s == a2) then
          stat_cache_out(t, s) = stat
          active_out(t, s) = .true.
          cycle
        end if
        b = a(level, s)
        stat_cache_out(t, s) = local_stat_fast(z, prefix, chain, b + 1, s, t, mode, window)
        evals = evals + 1
        score = f(level, s) + stat_cache_out(t, s)
        if (score >= ref_score) then
          active_out(t, s) = .true.
        else
          stat_cache_out(t, s) = 0.0_dp
          pruned = pruned + 1
        end if
      end do
      active_out(t, a2) = .true.
      stat_cache_out(t, a2) = stat
    end do
  end subroutine prune_level

  pure subroutine build_distance_prefix(d, prefix, chain)
    real(dp), intent(in) :: d(:,:) !! Symmetric pairwise distance matrix with a zero diagonal.
    real(dp), allocatable, intent(out) :: prefix(:,:) !! Two-dimensional block-sum prefix table with zero lower border.
    real(dp), allocatable, intent(out) :: chain(:) !! Prefix sums of adjacent distances d(i,i+1).
    integer :: i, j, n

    n = size(d, 1)
    allocate(prefix(0:n, 0:n), source=0.0_dp)
    do j = 1, n
      do i = 1, n
        prefix(i, j) = d(i, j) + prefix(i - 1, j) + prefix(i, j - 1) - prefix(i - 1, j - 1)
      end do
    end do
    allocate(chain(0:n), source=0.0_dp)
    do i = 1, n - 1
      chain(i) = chain(i - 1) + d(i, i + 1)
    end do
    if (n >= 1) chain(n) = chain(max(0, n - 1))
  end subroutine build_distance_prefix

  pure real(dp) function local_stat_fast(z, prefix, chain, first, split, last, mode, window) result(value)
    real(dp), intent(in) :: z(:,:) !! Time series matrix used by KS modes.
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix sums of pairwise distances, or zero-size for KS.
    real(dp), intent(in) :: chain(0:) !! Prefix sums of adjacent distances, or zero-size for KS.
    integer, intent(in) :: first !! First observation index of the left segment.
    integer, intent(in) :: split !! Last observation index of the left segment.
    integer, intent(in) :: last !! Last observation index of the right segment.
    integer, intent(in) :: mode !! Statistic selector used by the CP3O engine.
    integer, intent(in) :: window !! Window size for approximate modes.

    select case (mode)
    case (1)
      value = energy_stat_prefix(prefix, first, split, last)
    case (2)
      value = windowed_energy_stat_prefix(prefix, chain, first, split, last, window)
    case (3)
      value = ks_split_stat(z, first, split, last)
    case (4)
      value = ks_split_stat(z, first, split, last, window)
    case default
      value = -huge(1.0_dp)
    end select
  end function local_stat_fast

  pure real(dp) function energy_stat_prefix(prefix, first, split, last) result(value)
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix sums of pairwise distances.
    integer, intent(in) :: first !! First observation index of the combined interval.
    integer, intent(in) :: split !! Last observation index of the left segment.
    integer, intent(in) :: last !! Last observation index of the right segment.
    integer :: n, m
    real(dp) :: ll, rr, lr

    n = split - first + 1
    m = last - split
    if (n < 2 .or. m < 2) then
      value = -huge(1.0_dp)
      return
    end if
    ll = 0.5_dp*block_sum(prefix, first, split, first, split)
    rr = 0.5_dp*block_sum(prefix, split + 1, last, split + 1, last)
    lr = block_sum(prefix, first, split, split + 1, last)
    value = 2.0_dp*lr/real(n*m, dp)
    value = value - 2.0_dp*ll/real(n*(n - 1), dp)
    value = value - 2.0_dp*rr/real(m*(m - 1), dp)
    value = value*real(n*m, dp)/real((n + m)*(n + m), dp)
  end function energy_stat_prefix

  pure real(dp) function windowed_energy_stat_prefix(prefix, chain, first, split, last, delta) result(value)
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix sums of pairwise distances.
    real(dp), intent(in) :: chain(0:) !! Prefix sums of adjacent pairwise distances.
    integer, intent(in) :: first !! First observation index of the combined interval.
    integer, intent(in) :: split !! Last observation index of the left segment.
    integer, intent(in) :: last !! Last observation index of the right segment.
    integer, intent(in) :: delta !! Complete window size used on each side of the split.
    integer :: n, m, lf, rr0, cll, crr, clr
    real(dp) :: dll, drr, dlr

    n = split - first + 1
    m = last - split
    if (n < delta + 1 .or. m < delta + 1) then
      value = -huge(1.0_dp)
      return
    end if
    lf = split - delta + 1
    rr0 = split + 1
    dll = 0.5_dp*block_sum(prefix, lf, split, lf, split)
    drr = 0.5_dp*block_sum(prefix, rr0, split + delta, rr0, split + delta)
    dlr = block_sum(prefix, lf, split, rr0, split + delta)
    if (first < lf) dll = dll + chain(lf - 1) - chain(first - 1)
    if (split + delta < last) drr = drr + chain(last - 1) - chain(split + delta - 1)
    cll = delta*(delta - 1)/2 + (n - delta)
    crr = delta*(delta - 1)/2 + (m - delta)
    clr = delta*delta
    value = 2.0_dp*dlr/real(clr, dp) - dll/real(cll, dp) - drr/real(crr, dp)
    value = value*real(n*m, dp)/real((n + m)*(n + m), dp)
  end function windowed_energy_stat_prefix

  pure real(dp) function block_sum(prefix, row_first, row_last, col_first, col_last) result(value)
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix table with a zero lower border.
    integer, intent(in) :: row_first !! First row index of the requested block.
    integer, intent(in) :: row_last !! Last row index of the requested block.
    integer, intent(in) :: col_first !! First column index of the requested block.
    integer, intent(in) :: col_last !! Last column index of the requested block.

    value = prefix(row_last, col_last) - prefix(row_first - 1, col_last)
    value = value - prefix(row_last, col_first - 1) + prefix(row_first - 1, col_first - 1)
  end function block_sum

  pure logical function finite_score(value) result(ok)
    real(dp), intent(in) :: value !! Dynamic-programming score tested against the negative-infinity sentinel.

    ok = value > -0.5_dp*huge(1.0_dp)
  end function finite_score

  pure subroutine reconstruct_row(a, level, n, row)
    integer, intent(in) :: a(:,:) !! Dynamic-programming last-change matrix.
    integer, intent(in) :: level !! Number of change points represented by this row.
    integer, intent(in) :: n !! Final observation index used for reconstruction.
    integer, intent(out) :: row(:) !! Output change-point starts, padded with zeros after level entries.
    integer :: lev, current, cp

    row = 0
    current = n
    do lev = level, 1, -1
      cp = a(lev, current)
      if (cp <= 0) return
      row(lev) = cp + 1
      current = cp
    end do
  end subroutine reconstruct_row

  pure real(dp) function piecewise_sse(y, split) result(sse)
    real(dp), intent(in) :: y(:) !! Goodness-of-fit sequence indexed by number of change points.
    integer, intent(in) :: split !! Shared index where the two fitted lines meet.
    real(dp) :: a1, b1, a2, b2
    integer :: i

    call fit_line(y, 1, split, a1, b1)
    call fit_line(y, split, size(y), a2, b2)
    sse = 0.0_dp
    do i = 1, split
      sse = sse + (a1 + b1*real(i - 1, dp) - y(i))**2
    end do
    do i = split, size(y)
      sse = sse + (a2 + b2*real(i - 1, dp) - y(i))**2
    end do
  end function piecewise_sse

  pure subroutine fit_line(y, first, last, intercept, slope)
    real(dp), intent(in) :: y(:) !! Response sequence fitted by a line over the selected index range.
    integer, intent(in) :: first !! First included response index.
    integer, intent(in) :: last !! Last included response index.
    real(dp), intent(out) :: intercept !! Least-squares intercept using zero-based x coordinates.
    real(dp), intent(out) :: slope !! Least-squares slope using zero-based x coordinates.
    real(dp) :: xm, ym, num, den, xx
    integer :: i, count

    count = last - first + 1
    xm = 0.0_dp
    ym = 0.0_dp
    do i = first, last
      xm = xm + real(i - 1, dp)
      ym = ym + y(i)
    end do
    xm = xm/real(count, dp)
    ym = ym/real(count, dp)
    num = 0.0_dp
    den = 0.0_dp
    do i = first, last
      xx = real(i - 1, dp)
      num = num + (xx - xm)*(y(i) - ym)
      den = den + (xx - xm)**2
    end do
    if (den <= tiny(1.0_dp)) then
      slope = 0.0_dp
    else
      slope = num/den
    end if
    intercept = ym - slope*xm
  end subroutine fit_line

end module ecp_cp3o
