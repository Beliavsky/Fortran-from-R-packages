module ecp_divisive
  use ecp_kinds, only: dp
  use ecp_types, only: divisive_result
  use ecp_energy, only: distance_matrix_alpha, split_point
  implicit none
  private

  public :: e_divisive

contains

  function e_divisive(x, sig_level, r, k, min_size, alpha) result(res)
    real(dp), intent(in) :: x(:,:) !! Time series matrix with observations in rows.
    real(dp), intent(in), optional :: sig_level !! Permutation-test significance level; default is 0.05.
    integer, intent(in), optional :: r !! Number of within-segment permutations per proposed split; default is 199.
    integer, intent(in), optional :: k !! Number of change points to force; omission enables significance testing.
    integer, intent(in), optional :: min_size !! Minimum observations per segment; default is 30.
    real(dp), intent(in), optional :: alpha !! Energy-distance exponent in (0, 2]; default is one.
    type(divisive_result) :: res
    real(dp), allocatable :: d(:,:), pvals(:)
    integer, allocatable :: changes(:), sorted(:), perms(:)
    real(dp) :: level, aa, estat, pval
    integer :: rr, mm, n, max_changes, nchanges, split, interval_first, interval_last
    integer :: iter, target, over, f, proposed, ntests
    logical :: fixed_k

    n = size(x, 1)
    level = 0.05_dp
    if (present(sig_level)) level = sig_level
    rr = 199
    if (present(r)) rr = r
    mm = 30
    if (present(min_size)) mm = min_size
    aa = 1.0_dp
    if (present(alpha)) aa = alpha
    fixed_k = present(k)
    if (fixed_k) then
      target = max(0, k)
      rr = 0
    else
      target = max(0, n/mm - 1)
    end if
    max_changes = min(target, max(0, n/mm - 1))

    if (n == 0) then
      allocate(res%order_found(0), res%estimates(0), res%cluster(0), res%p_values(0), res%permutations(0))
      return
    end if
    d = distance_matrix_alpha(x, aa)
    allocate(changes(max_changes + 2), source=0)
    allocate(pvals(max_changes), source=0.0_dp)
    allocate(perms(max_changes), source=0)
    changes(1) = 1
    changes(2) = n + 1
    nchanges = 2
    res%k_hat = 1
    proposed = -1
    ntests = 0

    do iter = 1, max_changes
      sorted = sorted_prefix(changes, nchanges)
      call best_split(d, sorted, mm, split, estat, interval_first, interval_last)
      proposed = split
      if (split < 0) exit
      if (fixed_k) then
        pval = 0.0_dp
        f = 0
      else
        over = 0
        do f = 1, rr
          if (permuted_best_at_least(d, sorted, mm, estat)) over = over + 1
        end do
        pval = real(1 + over, dp)/real(rr + 1, dp)
      end if
      ntests = ntests + 1
      pvals(ntests) = pval
      if (fixed_k) then
        perms(ntests) = 0
      else
        perms(ntests) = rr
      end if
      if (.not. fixed_k) then
        if (pval > level) exit
      end if
      nchanges = nchanges + 1
      changes(nchanges) = split
      res%k_hat = res%k_hat + 1
    end do

    res%considered_last = proposed
    allocate(res%order_found(nchanges))
    res%order_found = changes(:nchanges)
    res%estimates = sorted_prefix(changes, nchanges)
    res%cluster = cluster_labels(res%estimates, n)
    allocate(res%p_values(ntests))
    allocate(res%permutations(ntests))
    if (ntests > 0) then
      res%p_values = pvals(:ntests)
      res%permutations = perms(:ntests)
    end if
  end function e_divisive

  pure subroutine best_split(d, changes, min_size, split, statistic, interval_first, interval_last)
    real(dp), intent(in) :: d(:,:) !! Pairwise alpha-powered distance matrix.
    integer, intent(in) :: changes(:) !! Sorted segment boundaries including 1 and n+1.
    integer, intent(in) :: min_size !! Minimum observations allowed in each child segment.
    integer, intent(out) :: split !! Best proposed change point as the first index of the right segment.
    real(dp), intent(out) :: statistic !! Largest split statistic over all current segments.
    integer, intent(out) :: interval_first !! First index of the segment containing the best split.
    integer, intent(out) :: interval_last !! Last index of the segment containing the best split.
    integer :: i, candidate_split
    real(dp) :: candidate_stat

    split = -1
    statistic = -huge(1.0_dp)
    interval_first = -1
    interval_last = -1
    do i = 1, size(changes) - 1
      call split_point(d, changes(i), changes(i + 1) - 1, min_size, candidate_split, candidate_stat)
      if (candidate_stat > statistic) then
        statistic = candidate_stat
        split = candidate_split
        interval_first = changes(i)
        interval_last = changes(i + 1) - 1
      end if
    end do
  end subroutine best_split

  logical function permuted_best_at_least(d, changes, min_size, observed) result(over)
    real(dp), intent(in) :: d(:,:) !! Pairwise distance matrix to permute within the current segments.
    integer, intent(in) :: changes(:) !! Sorted segment boundaries including 1 and n+1.
    integer, intent(in) :: min_size !! Minimum segment size used by split search.
    real(dp), intent(in) :: observed !! Observed best split statistic used as the permutation threshold.
    real(dp), allocatable :: dperm(:,:)
    integer, allocatable :: index(:)
    integer :: i, j, n, split, first, last
    real(dp) :: stat

    n = size(d, 1)
    allocate(index(n))
    index = [(i, i=1, n)]
    do i = 1, size(changes) - 1
      call shuffle_range(index, changes(i), changes(i + 1) - 1)
    end do
    allocate(dperm(n, n))
    do j = 1, n
      do i = 1, n
        dperm(i, j) = d(index(i), index(j))
      end do
    end do
    call best_split(dperm, changes, min_size, split, stat, first, last)
    over = stat >= observed
  end function permuted_best_at_least

  subroutine shuffle_range(index, first, last)
    integer, intent(inout) :: index(:) !! Permutation index vector modified within the selected range.
    integer, intent(in) :: first !! First position in index eligible for shuffling.
    integer, intent(in) :: last !! Last position in index eligible for shuffling.
    integer :: i, j, tmp
    real(dp) :: u

    do i = last, first + 1, -1
      call random_number(u)
      j = first + int(u*real(i - first + 1, dp))
      j = min(j, i)
      tmp = index(i)
      index(i) = index(j)
      index(j) = tmp
    end do
  end subroutine shuffle_range

  pure function sorted_prefix(values, nuse) result(sorted)
    integer, intent(in) :: values(:) !! Integer storage containing active values in its leading positions.
    integer, intent(in) :: nuse !! Number of leading values to copy and sort.
    integer, allocatable :: sorted(:)
    integer :: i, j, key

    allocate(sorted(nuse))
    sorted = values(:nuse)
    do i = 2, nuse
      key = sorted(i)
      j = i - 1
      do while (j >= 1)
        if (sorted(j) <= key) exit
        sorted(j + 1) = sorted(j)
        j = j - 1
      end do
      sorted(j + 1) = key
    end do
  end function sorted_prefix

  pure function cluster_labels(boundaries, n) result(cluster)
    integer, intent(in) :: boundaries(:) !! Sorted segment boundaries including 1 and n+1.
    integer, intent(in) :: n !! Number of observations to label.
    integer, allocatable :: cluster(:)
    integer :: i, j

    allocate(cluster(n), source=0)
    do i = 1, size(boundaries) - 1
      do j = boundaries(i), boundaries(i + 1) - 1
        if (j >= 1 .and. j <= n) cluster(j) = i
      end do
    end do
  end function cluster_labels

end module ecp_divisive
