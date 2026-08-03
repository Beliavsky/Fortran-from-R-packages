! SPDX-License-Identifier: GPL-2.0-only
module hierportfolios_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use hierportfolios_kinds, only: dp
  use hierportfolios_types, only: portfolio_result, hierarchy_result, hp_success, &
    hp_invalid_argument, hp_numerical_failure
  use hierportfolios_hierarchy, only: covariance_distance, hclust_distance, &
    diana_distance, cutree, select_gap_clusters
  implicit none
  private

  public :: HRP_Portfolio
  public :: HCAA_Portfolio
  public :: HERC_Portfolio
  public :: DHRP_Portfolio

contains

  subroutine HRP_Portfolio(covar, result, linkage)
    real(dp), intent(in) :: covar(:, :)
    type(portfolio_result), intent(out) :: result
    character(len=*), intent(in), optional :: linkage

    type(hierarchy_result) :: hierarchy
    real(dp), allocatable :: features(:, :), distances(:, :)
    character(len=:), allocatable :: link, message
    integer :: status

    call init_result(result, 'HRP')
    link = 'single'
    if (present(linkage)) link = canonical_linkage(linkage)
    if (len(link) == 0) then
      call fail_result(result, hp_invalid_argument, 'invalid linkage method')
      return
    end if
    call covariance_distance(covar, features, distances, status, message)
    if (status /= hp_success) then
      call fail_result(result, status, message)
      return
    end if
    call hclust_distance(distances, link, hierarchy)
    if (.not. hierarchy%ok()) then
      call fail_result(result, hierarchy%status, hierarchy%message)
      return
    end if
    call recursive_bisection(covar, hierarchy%order, result%weights, result%iterations)
    result%order = hierarchy%order
    allocate(result%clusters(size(covar, 1)))
    result%clusters = 1
    result%n_clusters = 1
    call finalize_weights(result)
  end subroutine HRP_Portfolio

  subroutine HCAA_Portfolio(covar, result, linkage, clusters, gap_references, seed)
    real(dp), intent(in) :: covar(:, :)
    type(portfolio_result), intent(out) :: result
    character(len=*), intent(in), optional :: linkage
    integer, intent(in), optional :: clusters, gap_references, seed

    type(hierarchy_result) :: hierarchy
    real(dp), allocatable :: features(:, :), distances(:, :)
    integer, allocatable :: current(:), previous(:)
    character(len=:), allocatable :: link, message
    integer :: b, c, i, k, kmax, n, selected, status

    call init_result(result, 'HCAA')
    link = 'ward'
    if (present(linkage)) link = canonical_linkage(linkage)
    if (len(link) == 0) then
      call fail_result(result, hp_invalid_argument, 'invalid linkage method')
      return
    end if
    call covariance_distance(covar, features, distances, status, message)
    if (status /= hp_success) then
      call fail_result(result, status, message)
      return
    end if
    call hclust_distance(distances, link, hierarchy)
    if (.not. hierarchy%ok()) then
      call fail_result(result, hierarchy%status, hierarchy%message)
      return
    end if
    n = size(covar, 1)
    if (present(clusters)) then
      selected = clusters
      if (selected < 2 .or. selected > n) then
        call fail_result(result, hp_invalid_argument, &
          'clusters must be between two and the number of assets')
        return
      end if
    else
      kmax = max(2, n / 2)
      b = 100
      if (present(gap_references)) b = gap_references
      if (b < 1) then
        call fail_result(result, hp_invalid_argument, 'gap_references must be positive')
        return
      end if
      if (present(seed)) then
        call select_gap_clusters(features, link, kmax, b, selected, &
          result%gap, result%gap_se, seed=seed, status=status)
      else
        call select_gap_clusters(features, link, kmax, b, selected, &
          result%gap, result%gap_se, status=status)
      end if
      if (status /= hp_success) then
        call fail_result(result, status, 'gap-statistic cluster selection failed')
        return
      end if
      selected = max(2, selected)
    end if

    allocate(result%weights(n))
    result%weights = 1.0_dp
    call cutree(hierarchy, 2, current, status)
    if (status /= hp_success) then
      call fail_result(result, status, 'could not cut hierarchy')
      return
    end if
    do c = 1, 2
      where (current == c) result%weights = result%weights / 2.0_dp
    end do
    previous = current
    do k = 3, selected
      call cutree(hierarchy, k, current, status)
      if (status /= hp_success) then
        call fail_result(result, status, 'could not cut hierarchy')
        return
      end if
      do c = 1, maxval(previous)
        if (count_distinct(pack(current, previous == c)) > 1) then
          where (previous == c) result%weights = result%weights / 2.0_dp
          exit
        end if
      end do
      previous = current
    end do
    result%clusters = current
    do c = 1, selected
      i = count(current == c)
      if (i > 0) then
        where (current == c) result%weights = result%weights / real(i, dp)
      end if
    end do
    result%n_clusters = selected
    result%order = hierarchy%order
    result%iterations = max(1, selected - 1)
    call finalize_weights(result)
  end subroutine HCAA_Portfolio

  subroutine HERC_Portfolio(covar, result, linkage, clusters, gap_references, seed)
    real(dp), intent(in) :: covar(:, :)
    type(portfolio_result), intent(out) :: result
    character(len=*), intent(in), optional :: linkage
    integer, intent(in), optional :: clusters, gap_references, seed

    type(hierarchy_result) :: hierarchy
    real(dp), allocatable :: features(:, :), distances(:, :), local_w(:)
    real(dp), allocatable :: cluster_risk(:), alpha(:)
    integer, allocatable :: current(:), final_labels(:), previous(:)
    character(len=:), allocatable :: link, message
    real(dp) :: risk_a, risk_b, split_alpha
    integer :: b, c, child_a, child_b, k, kmax, n, selected, status

    call init_result(result, 'HERC')
    link = 'ward'
    if (present(linkage)) link = canonical_linkage(linkage)
    if (len(link) == 0) then
      call fail_result(result, hp_invalid_argument, 'invalid linkage method')
      return
    end if
    call covariance_distance(covar, features, distances, status, message)
    if (status /= hp_success) then
      call fail_result(result, status, message)
      return
    end if
    call hclust_distance(distances, link, hierarchy)
    if (.not. hierarchy%ok()) then
      call fail_result(result, hierarchy%status, hierarchy%message)
      return
    end if
    n = size(covar, 1)
    if (present(clusters)) then
      selected = clusters
      if (selected < 2 .or. selected > n) then
        call fail_result(result, hp_invalid_argument, &
          'clusters must be between two and the number of assets')
        return
      end if
    else
      kmax = max(2, n / 2)
      b = 100
      if (present(gap_references)) b = gap_references
      if (b < 1) then
        call fail_result(result, hp_invalid_argument, 'gap_references must be positive')
        return
      end if
      if (present(seed)) then
        call select_gap_clusters(features, link, kmax, b, selected, &
          result%gap, result%gap_se, seed=seed, status=status)
      else
        call select_gap_clusters(features, link, kmax, b, selected, &
          result%gap, result%gap_se, status=status)
      end if
      if (status /= hp_success) then
        call fail_result(result, status, 'gap-statistic cluster selection failed')
        return
      end if
      selected = max(2, selected)
    end if

    call cutree(hierarchy, selected, final_labels, status)
    if (status /= hp_success) then
      call fail_result(result, status, 'could not cut hierarchy')
      return
    end if
    allocate(local_w(n), cluster_risk(selected), alpha(n))
    local_w = 0.0_dp
    cluster_risk = 0.0_dp
    do c = 1, selected
      call inverse_variance_cluster(covar, pack([(k, k=1,n)], final_labels == c), &
        local_w, cluster_risk(c))
    end do
    alpha = 1.0_dp
    allocate(previous(n))
    previous = 1
    do k = 2, selected
      call cutree(hierarchy, k, current, status)
      if (status /= hp_success) then
        call fail_result(result, status, 'could not cut hierarchy')
        return
      end if
      do c = 1, maxval(previous)
        if (count_distinct(pack(current, previous == c)) > 1) then
          call first_two_groups(pack(current, previous == c), child_a, child_b)
          risk_a = terminal_group_risk(current == child_a, final_labels, cluster_risk)
          risk_b = terminal_group_risk(current == child_b, final_labels, cluster_risk)
          if (risk_a + risk_b <= tiny(1.0_dp)) then
            split_alpha = 0.5_dp
          else
            split_alpha = risk_b / (risk_a + risk_b)
          end if
          where (current == child_a) alpha = alpha * split_alpha
          where (current == child_b) alpha = alpha * (1.0_dp - split_alpha)
          exit
        end if
      end do
      previous = current
    end do
    allocate(result%weights(n))
    result%weights = local_w * alpha
    result%clusters = final_labels
    result%n_clusters = selected
    result%order = hierarchy%order
    result%iterations = max(1, selected - 1)
    call finalize_weights(result)
  end subroutine HERC_Portfolio

  subroutine DHRP_Portfolio(covar, result, tau, ub, lb)
    real(dp), intent(in) :: covar(:, :)
    type(portfolio_result), intent(out) :: result
    real(dp), intent(in), optional :: tau
    real(dp), intent(in), optional :: ub(:), lb(:)

    type(hierarchy_result) :: hierarchy
    real(dp), allocatable :: features(:, :), distances(:, :), lower(:), upper(:)
    real(dp) :: alpha, score, tau_value, va, vb
    integer, allocatable :: end_pos(:), indices_a(:), indices_b(:), start_pos(:)
    character(len=:), allocatable :: message
    integer :: best, candidate, e, head, high, low, m, n, s, split, status, tail

    call init_result(result, 'DHRP')
    n = size(covar, 1)
    tau_value = 1.0_dp
    if (present(tau)) tau_value = tau
    if (tau_value < 0.0_dp .or. tau_value > 1.0_dp) then
      call fail_result(result, hp_invalid_argument, 'tau must lie between zero and one')
      return
    end if
    allocate(lower(n), upper(n))
    lower = 0.0_dp
    upper = 1.0_dp
    if (present(lb)) then
      if (size(lb) /= n) then
        call fail_result(result, hp_invalid_argument, 'LB has the wrong size')
        return
      end if
      lower = lb
    end if
    if (present(ub)) then
      if (size(ub) /= n) then
        call fail_result(result, hp_invalid_argument, 'UB has the wrong size')
        return
      end if
      upper = ub
    end if
    if (any(.not. ieee_is_finite(lower)) .or. any(.not. ieee_is_finite(upper)) .or. &
        any(lower < 0.0_dp) .or. any(upper > 1.0_dp) .or. any(lower > upper) .or. &
        sum(lower) > 1.0_dp + 100.0_dp * epsilon(1.0_dp) .or. &
        sum(upper) < 1.0_dp - 100.0_dp * epsilon(1.0_dp)) then
      call fail_result(result, hp_invalid_argument, 'infeasible portfolio bounds')
      return
    end if

    call covariance_distance(covar, features, distances, status, message)
    if (status /= hp_success) then
      call fail_result(result, status, message)
      return
    end if
    call diana_distance(distances, hierarchy)
    if (.not. hierarchy%ok()) then
      call fail_result(result, hierarchy%status, hierarchy%message)
      return
    end if
    allocate(result%weights(n), start_pos(2*n), end_pos(2*n))
    result%weights = 1.0_dp
    head = 1
    tail = 1
    start_pos(1) = 1
    end_pos(1) = n
    result%iterations = 0
    do while (head <= tail)
      s = start_pos(head)
      e = end_pos(head)
      head = head + 1
      m = e - s + 1
      if (m <= 1) cycle
      low = max(1, floor(real(m, dp) / 2.0_dp - &
        (real(m, dp) / 2.0_dp - 1.0_dp) * tau_value))
      high = min(m - 1, floor(real(m, dp) / 2.0_dp + &
        (real(m, dp) / 2.0_dp - 1.0_dp) * tau_value))
      if (low > high) then
        low = max(1, m / 2)
        high = low
      end if
      best = low
      score = -huge(1.0_dp)
      do candidate = low, high
        if (hierarchy%height(min(n - 1, s + candidate - 1)) > score) then
          score = hierarchy%height(min(n - 1, s + candidate - 1))
          best = candidate
        end if
      end do
      split = s + best - 1
      indices_a = hierarchy%order(s:split)
      indices_b = hierarchy%order(split+1:e)
      va = cluster_variance(covar, indices_a)
      vb = cluster_variance(covar, indices_b)
      if (va + vb <= tiny(1.0_dp)) then
        alpha = 0.5_dp
      else
        alpha = vb / (va + vb)
      end if
      result%weights(indices_a) = result%weights(indices_a) * alpha
      result%weights(indices_b) = result%weights(indices_b) * (1.0_dp - alpha)
      if (size(indices_a) > 1) then
        tail = tail + 1
        start_pos(tail) = s
        end_pos(tail) = split
      end if
      if (size(indices_b) > 1) then
        tail = tail + 1
        start_pos(tail) = split + 1
        end_pos(tail) = e
      end if
      result%iterations = result%iterations + 1
    end do
    call project_box_simplex(result%weights, lower, upper, status)
    if (status /= hp_success) then
      call fail_result(result, status, 'could not enforce portfolio bounds')
      return
    end if
    result%order = hierarchy%order
    allocate(result%clusters(n))
    result%clusters = 1
    result%n_clusters = 1
    call finalize_weights(result)
  end subroutine DHRP_Portfolio

  subroutine recursive_bisection(covar, order, weights, iterations)
    real(dp), intent(in) :: covar(:, :)
    integer, intent(in) :: order(:)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: iterations

    integer, allocatable :: end_pos(:), indices_a(:), indices_b(:), start_pos(:)
    real(dp) :: alpha, va, vb
    integer :: e, head, m, n, s, split, tail

    n = size(order)
    allocate(weights(n), start_pos(2*n), end_pos(2*n))
    weights = 1.0_dp
    head = 1
    tail = 1
    start_pos(1) = 1
    end_pos(1) = n
    iterations = 0
    do while (head <= tail)
      s = start_pos(head)
      e = end_pos(head)
      head = head + 1
      m = e - s + 1
      if (m <= 1) cycle
      split = s + m / 2 - 1
      indices_a = order(s:split)
      indices_b = order(split+1:e)
      va = cluster_variance(covar, indices_a)
      vb = cluster_variance(covar, indices_b)
      if (va + vb <= tiny(1.0_dp)) then
        alpha = 0.5_dp
      else
        alpha = vb / (va + vb)
      end if
      weights(indices_a) = weights(indices_a) * alpha
      weights(indices_b) = weights(indices_b) * (1.0_dp - alpha)
      if (size(indices_a) > 1) then
        tail = tail + 1
        start_pos(tail) = s
        end_pos(tail) = split
      end if
      if (size(indices_b) > 1) then
        tail = tail + 1
        start_pos(tail) = split + 1
        end_pos(tail) = e
      end if
      iterations = iterations + 1
    end do
  end subroutine recursive_bisection

  real(dp) function cluster_variance(covar, indices) result(value)
    real(dp), intent(in) :: covar(:, :)
    integer, intent(in) :: indices(:)

    real(dp), allocatable :: w(:)
    integer :: i, j

    allocate(w(size(indices)))
    do i = 1, size(indices)
      w(i) = 1.0_dp / covar(indices(i), indices(i))
    end do
    w = w / sum(w)
    value = 0.0_dp
    do j = 1, size(indices)
      do i = 1, size(indices)
        value = value + w(i) * covar(indices(i), indices(j)) * w(j)
      end do
    end do
  end function cluster_variance

  subroutine inverse_variance_cluster(covar, indices, all_weights, variance)
    real(dp), intent(in) :: covar(:, :)
    integer, intent(in) :: indices(:)
    real(dp), intent(inout) :: all_weights(:)
    real(dp), intent(out) :: variance

    real(dp), allocatable :: w(:)
    integer :: i, j

    allocate(w(size(indices)))
    do i = 1, size(indices)
      w(i) = 1.0_dp / covar(indices(i), indices(i))
    end do
    w = w / sum(w)
    all_weights(indices) = w
    variance = 0.0_dp
    do j = 1, size(indices)
      do i = 1, size(indices)
        variance = variance + w(i) * covar(indices(i), indices(j)) * w(j)
      end do
    end do
  end subroutine inverse_variance_cluster

  real(dp) function terminal_group_risk(mask, final_labels, cluster_risk) result(value)
    logical, intent(in) :: mask(:)
    integer, intent(in) :: final_labels(:)
    real(dp), intent(in) :: cluster_risk(:)

    logical, allocatable :: seen(:)
    integer :: i, label

    allocate(seen(size(cluster_risk)))
    seen = .false.
    value = 0.0_dp
    do i = 1, size(mask)
      if (.not. mask(i)) cycle
      label = final_labels(i)
      if (.not. seen(label)) then
        value = value + cluster_risk(label)
        seen(label) = .true.
      end if
    end do
  end function terminal_group_risk

  subroutine first_two_groups(values, first, second)
    integer, intent(in) :: values(:)
    integer, intent(out) :: first, second
    integer :: i

    first = values(1)
    second = first
    do i = 2, size(values)
      if (values(i) /= first) then
        second = values(i)
        return
      end if
    end do
  end subroutine first_two_groups

  integer function count_distinct(values) result(count_value)
    integer, intent(in) :: values(:)
    integer :: i, j
    logical :: found

    count_value = 0
    do i = 1, size(values)
      found = .false.
      do j = 1, i - 1
        if (values(j) == values(i)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) count_value = count_value + 1
    end do
  end function count_distinct

  subroutine project_box_simplex(values, lower, upper, status)
    real(dp), intent(inout) :: values(:)
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(out) :: status

    real(dp), allocatable :: capacity(:), original(:)
    real(dp) :: lambda, left, residual, right, total
    integer :: iter

    if (size(values) /= size(lower) .or. size(values) /= size(upper)) then
      status = hp_invalid_argument
      return
    end if
    if (all(values >= lower) .and. all(values <= upper) .and. &
        abs(sum(values) - 1.0_dp) < 1.0e-13_dp) then
      status = hp_success
      return
    end if
    allocate(original(size(values)), capacity(size(values)))
    original = values
    left = minval(original - upper) - 1.0_dp
    right = maxval(original - lower) + 1.0_dp
    do iter = 1, 250
      lambda = 0.5_dp * (left + right)
      total = sum(max(lower, min(upper, original - lambda)))
      if (total > 1.0_dp) then
        left = lambda
      else
        right = lambda
      end if
      if (abs(total - 1.0_dp) < 1.0e-14_dp) exit
    end do
    values = max(lower, min(upper, original - 0.5_dp * (left + right)))
    residual = 1.0_dp - sum(values)
    if (residual > 0.0_dp) then
      capacity = max(0.0_dp, upper - values)
      if (sum(capacity) > 0.0_dp) values = values + residual * capacity / sum(capacity)
    else if (residual < 0.0_dp) then
      capacity = max(0.0_dp, values - lower)
      if (sum(capacity) > 0.0_dp) values = values + residual * capacity / sum(capacity)
    end if
    values = max(lower, min(upper, values))
    if (abs(sum(values) - 1.0_dp) > 1.0e-10_dp) then
      status = hp_numerical_failure
    else
      status = hp_success
    end if
  end subroutine project_box_simplex

  function canonical_linkage(name) result(link)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: link
    character(len=:), allocatable :: value

    value = lower_ascii(trim(adjustl(name)))
    select case (value)
    case ('single')
      link = 'single'
    case ('complete')
      link = 'complete'
    case ('average')
      link = 'average'
    case ('ward', 'ward.d2')
      link = 'ward'
    case default
      link = ''
    end select
  end function canonical_linkage

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: code, i

    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      else
        lower(i:i) = text(i:i)
      end if
    end do
  end function lower_ascii

  subroutine init_result(result, method)
    type(portfolio_result), intent(out) :: result
    character(len=*), intent(in) :: method

    result%status = hp_success
    result%message = 'ok'
    result%method = method
    result%n_clusters = 1
    result%iterations = 0
  end subroutine init_result

  subroutine finalize_weights(result)
    type(portfolio_result), intent(inout) :: result

    if (.not. allocated(result%weights) .or. any(.not. ieee_is_finite(result%weights)) .or. &
        any(result%weights < -100.0_dp * epsilon(1.0_dp)) .or. &
        sum(result%weights) <= tiny(1.0_dp)) then
      call fail_result(result, hp_numerical_failure, 'invalid portfolio weights')
      return
    end if
    result%weights = max(0.0_dp, result%weights)
    result%weights = result%weights / sum(result%weights)
    result%status = hp_success
    result%message = 'ok'
  end subroutine finalize_weights

  subroutine fail_result(result, status, message)
    type(portfolio_result), intent(inout) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message

    result%status = status
    result%message = message
  end subroutine fail_result

end module hierportfolios_core
