! SPDX-License-Identifier: GPL-2.0-or-later
module cluster_hierarchy
  use fastcluster_kinds, only: dp
  use fastcluster_types, only: hclust_result
  use fastcluster_core, only: hclust_matrix
  use cluster_types, only: hierarchy_result, mona_result, cluster_success, &
    cluster_invalid_argument, cluster_numerical_failure
  use cluster_daisy, only: daisy
  implicit none
  private

  public :: agnes
  public :: agnes_distance
  public :: diana
  public :: diana_distance
  public :: mona
  public :: coef_hier

contains

  subroutine agnes(x, result, method, metric)
    real(dp), intent(in) :: x(:, :)
    type(hierarchy_result), intent(out) :: result
    character(len=*), intent(in), optional :: method, metric

    real(dp), allocatable :: distances(:, :)
    character(len=:), allocatable :: message
    integer :: status

    if (present(metric)) then
      call daisy(x, distances, metric=metric, status=status, message=message)
    else
      call daisy(x, distances, status=status, message=message)
    end if
    if (status /= cluster_success) then
      call fail_hierarchy(result, status, message)
      return
    end if
    if (present(method)) then
      call agnes_distance(distances, result, method)
    else
      call agnes_distance(distances, result, 'average')
    end if
  end subroutine agnes

  subroutine agnes_distance(distances, result, method)
    real(dp), intent(in) :: distances(:, :)
    type(hierarchy_result), intent(out) :: result
    character(len=*), intent(in), optional :: method

    type(hclust_result) :: hresult
    character(len=:), allocatable :: link

    link = 'average'
    if (present(method)) link = canonical_agnes_method(method)
    if (len(link) == 0) then
      call fail_hierarchy(result, cluster_invalid_argument, &
        'supported AGNES methods are average, single, complete, weighted, and ward')
      return
    end if
    call hclust_matrix(distances, link, hresult)
    if (.not. hresult%ok()) then
      call fail_hierarchy(result, hresult%status, hresult%message)
      return
    end if
    result%n = hresult%n
    result%method = 'agnes-'//link
    result%merge = hresult%merge
    result%height = hresult%height
    result%order = hresult%order
    result%coefficient = agglomerative_coefficient(result%merge, result%height)
    result%status = cluster_success
    result%message = 'ok'
  end subroutine agnes_distance

  subroutine diana(x, result, metric)
    real(dp), intent(in) :: x(:, :)
    type(hierarchy_result), intent(out) :: result
    character(len=*), intent(in), optional :: metric

    real(dp), allocatable :: distances(:, :)
    character(len=:), allocatable :: message
    integer :: status

    if (present(metric)) then
      call daisy(x, distances, metric=metric, status=status, message=message)
    else
      call daisy(x, distances, status=status, message=message)
    end if
    if (status /= cluster_success) then
      call fail_hierarchy(result, status, message)
      return
    end if
    call diana_distance(distances, result)
  end subroutine diana

  subroutine diana_distance(distances, result)
    real(dp), intent(in) :: distances(:, :)
    type(hierarchy_result), intent(out) :: result

    integer, allocatable :: child1(:), child2(:), items(:), order(:)
    real(dp), allocatable :: node_height(:), leaf_parent_height(:)
    integer :: counter, n, root

    n = size(distances, 1)
    if (n < 2 .or. size(distances, 2) /= n .or. any(distances < 0.0_dp)) then
      call fail_hierarchy(result, cluster_invalid_argument, &
        'distances must be a nonnegative square matrix of order at least two')
      return
    end if
    allocate(child1(2*n-1), child2(2*n-1), node_height(2*n-1), &
      leaf_parent_height(n), items(n), order(n))
    child1 = 0
    child2 = 0
    node_height = 0.0_dp
    leaf_parent_height = 0.0_dp
    items = [(counter, counter=1,n)]
    counter = n
    call build_divisive_tree(distances, items, counter, root, child1, child2, node_height)
    if (root /= 2*n-1) then
      call fail_hierarchy(result, cluster_numerical_failure, 'DIANA tree construction failed')
      return
    end if
    allocate(result%merge(n-1, 2), result%height(n-1), result%order(n))
    call tree_to_merge(n, child1, child2, node_height, result%merge, result%height)
    counter = 0
    call tree_order(root, n, child1, child2, order, counter)
    result%order = order
    call leaf_parent_heights(root, 0.0_dp, n, child1, child2, node_height, leaf_parent_height)
    if (node_height(root) > 0.0_dp) then
      result%coefficient = max(0.0_dp, min(1.0_dp, &
        1.0_dp - sum(leaf_parent_height) / (real(n, dp) * node_height(root))))
    else
      result%coefficient = 0.0_dp
    end if
    result%n = n
    result%method = 'diana'
    result%status = cluster_success
    result%message = 'ok'
  end subroutine diana_distance

  recursive subroutine build_divisive_tree(distances, items, counter, node, &
      child1, child2, node_height)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(in) :: items(:)
    integer, intent(inout) :: counter
    integer, intent(out) :: node
    integer, intent(inout) :: child1(:), child2(:)
    real(dp), intent(inout) :: node_height(:)

    integer, allocatable :: left(:), right(:)
    integer :: a, b, i, j
    real(dp) :: diameter

    if (size(items) == 1) then
      node = items(1)
      return
    end if
    call diana_split(distances, items, left, right)
    call build_divisive_tree(distances, left, counter, a, child1, child2, node_height)
    call build_divisive_tree(distances, right, counter, b, child1, child2, node_height)
    counter = counter + 1
    node = counter
    child1(node) = a
    child2(node) = b
    diameter = 0.0_dp
    do i = 1, size(items) - 1
      do j = i + 1, size(items)
        diameter = max(diameter, distances(items(i), items(j)))
      end do
    end do
    node_height(node) = diameter
  end subroutine build_divisive_tree

  subroutine diana_split(distances, items, left, right)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(in) :: items(:)
    integer, allocatable, intent(out) :: left(:), right(:)

    logical, allocatable :: in_splinter(:)
    real(dp) :: average, best_average, difference, best_difference, to_remainder, to_splinter
    integer :: best, i, j, m, nrem, nspl

    m = size(items)
    allocate(in_splinter(m))
    in_splinter = .false.
    best = 1
    best_average = -1.0_dp
    do i = 1, m
      average = 0.0_dp
      do j = 1, m
        if (j /= i) average = average + distances(items(i), items(j))
      end do
      average = average / real(m - 1, dp)
      if (average > best_average) then
        best_average = average
        best = i
      end if
    end do
    in_splinter(best) = .true.
    do
      nspl = count(in_splinter)
      nrem = m - nspl
      if (nrem <= 1) exit
      best = 0
      best_difference = 0.0_dp
      do i = 1, m
        if (in_splinter(i)) cycle
        to_splinter = 0.0_dp
        do j = 1, m
          if (in_splinter(j)) to_splinter = to_splinter + distances(items(i), items(j))
        end do
        to_splinter = to_splinter / real(nspl, dp)
        to_remainder = 0.0_dp
        do j = 1, m
          if (.not. in_splinter(j) .and. j /= i) then
            to_remainder = to_remainder + distances(items(i), items(j))
          end if
        end do
        to_remainder = to_remainder / real(nrem - 1, dp)
        difference = to_remainder - to_splinter
        if (difference > best_difference) then
          best_difference = difference
          best = i
        end if
      end do
      if (best == 0) exit
      in_splinter(best) = .true.
    end do
    allocate(left(count(in_splinter)), right(m-count(in_splinter)))
    left = pack(items, in_splinter)
    right = pack(items, .not. in_splinter)
  end subroutine diana_split

  subroutine mona(x, result, max_clusters)
    integer, intent(in) :: x(:, :)
    type(mona_result), intent(out) :: result
    integer, intent(in), optional :: max_clusters

    integer, allocatable :: labels(:), order(:), split_vars(:)
    integer :: c, count0, count1, i, j, kmax, n, next_cluster, p, selected, best_balance
    logical :: changed

    n = size(x, 1)
    p = size(x, 2)
    if (n < 2 .or. p < 1 .or. any((x /= 0) .and. (x /= 1))) then
      result%status = cluster_invalid_argument
      result%message = 'MONA requires a nonempty binary 0/1 matrix'
      return
    end if
    kmax = n
    if (present(max_clusters)) kmax = min(n, max(1, max_clusters))
    allocate(labels(n), split_vars(max(1, kmax-1)), order(n))
    labels = 1
    split_vars = 0
    next_cluster = 1
    do while (next_cluster < kmax)
      changed = .false.
      do c = 1, next_cluster
        selected = 0
        best_balance = 0
        do j = 1, p
          count0 = 0
          count1 = 0
          do i = 1, n
            if (labels(i) /= c) cycle
            if (x(i, j) == 0) count0 = count0 + 1
            if (x(i, j) == 1) count1 = count1 + 1
          end do
          if (min(count0, count1) > best_balance) then
            best_balance = min(count0, count1)
            selected = j
          end if
        end do
        if (selected > 0) then
          next_cluster = next_cluster + 1
          split_vars(next_cluster-1) = selected
          do i = 1, n
            if (labels(i) == c .and. x(i, selected) == 1) labels(i) = next_cluster
          end do
          changed = .true.
          if (next_cluster >= kmax) exit
        end if
      end do
      if (.not. changed) exit
    end do
    call order_by_labels(labels, order)
    result%n = n
    result%n_clusters = next_cluster
    result%clustering = labels
    result%order = order
    result%split_variable = split_vars(1:max(0,next_cluster-1))
    result%coefficient = 1.0_dp - real(next_cluster, dp) / real(max(2, n), dp)
    result%status = cluster_success
    result%message = 'ok'
  end subroutine mona

  real(dp) function coef_hier(result) result(value)
    type(hierarchy_result), intent(in) :: result
    value = result%coefficient
  end function coef_hier

  real(dp) function agglomerative_coefficient(merge, height) result(coefficient)
    integer, intent(in) :: merge(:, :)
    real(dp), intent(in) :: height(:)

    real(dp), allocatable :: first_height(:)
    integer, allocatable :: members(:, :), counts(:)
    integer :: a, b, i, j, n, row

    n = size(merge, 1) + 1
    allocate(first_height(n), members(n-1, n), counts(n-1))
    first_height = -1.0_dp
    members = 0
    counts = 0
    do row = 1, n - 1
      call node_members(merge(row, 1), counts, a)
      call node_members(merge(row, 2), counts, b)
      counts(row) = a + b
      call copy_node_members(merge(row, 1), members, counts, members(row, 1:a))
      call copy_node_members(merge(row, 2), members, counts, members(row, a+1:a+b))
      do i = 1, counts(row)
        j = members(row, i)
        if (first_height(j) < 0.0_dp) first_height(j) = height(row)
      end do
    end do
    if (height(n-1) > 0.0_dp) then
      coefficient = max(0.0_dp, min(1.0_dp, &
        1.0_dp - sum(first_height) / (real(n, dp) * height(n-1))))
    else
      coefficient = 0.0_dp
    end if
  end function agglomerative_coefficient

  subroutine node_members(code, counts, count_value)
    integer, intent(in) :: code
    integer, intent(in) :: counts(:)
    integer, intent(out) :: count_value
    if (code < 0) then
      count_value = 1
    else
      count_value = counts(code)
    end if
  end subroutine node_members

  subroutine copy_node_members(code, members, counts, output)
    integer, intent(in) :: code
    integer, intent(in) :: members(:, :), counts(:)
    integer, intent(out) :: output(:)
    if (code < 0) then
      output(1) = -code
    else
      output = members(code, 1:counts(code))
    end if
  end subroutine copy_node_members

  subroutine tree_to_merge(n, child1, child2, node_height, merge, height)
    integer, intent(in) :: n, child1(:), child2(:)
    real(dp), intent(in) :: node_height(:)
    integer, intent(out) :: merge(:, :)
    real(dp), intent(out) :: height(:)
    integer :: node, row
    do node = n + 1, 2*n - 1
      row = node - n
      merge(row, 1) = merge_code(child1(node), n)
      merge(row, 2) = merge_code(child2(node), n)
      height(row) = node_height(node)
    end do
  end subroutine tree_to_merge

  integer function merge_code(node, n) result(code)
    integer, intent(in) :: node, n
    if (node <= n) then
      code = -node
    else
      code = node - n
    end if
  end function merge_code

  recursive subroutine tree_order(node, n, child1, child2, order, count_value)
    integer, intent(in) :: node, n, child1(:), child2(:)
    integer, intent(inout) :: order(:), count_value
    if (node <= n) then
      count_value = count_value + 1
      order(count_value) = node
    else
      call tree_order(child1(node), n, child1, child2, order, count_value)
      call tree_order(child2(node), n, child1, child2, order, count_value)
    end if
  end subroutine tree_order

  recursive subroutine leaf_parent_heights(node, parent_height, n, child1, child2, &
      node_height, output)
    integer, intent(in) :: node, n, child1(:), child2(:)
    real(dp), intent(in) :: parent_height, node_height(:)
    real(dp), intent(inout) :: output(:)
    if (node <= n) then
      output(node) = parent_height
    else
      call leaf_parent_heights(child1(node), node_height(node), n, child1, child2, node_height, output)
      call leaf_parent_heights(child2(node), node_height(node), n, child1, child2, node_height, output)
    end if
  end subroutine leaf_parent_heights

  subroutine order_by_labels(labels, order)
    integer, intent(in) :: labels(:)
    integer, intent(out) :: order(:)
    integer :: c, count_value, i
    count_value = 0
    do c = 1, maxval(labels)
      do i = 1, size(labels)
        if (labels(i) == c) then
          count_value = count_value + 1
          order(count_value) = i
        end if
      end do
    end do
  end subroutine order_by_labels

  function canonical_agnes_method(method) result(name)
    character(len=*), intent(in) :: method
    character(len=:), allocatable :: name
    character(len=:), allocatable :: lower
    lower = lower_ascii(trim(adjustl(method)))
    select case (lower)
    case ('average', 'single', 'complete', 'weighted')
      name = lower
    case ('ward', 'ward.d2', 'ward.d')
      name = 'ward.d2'
    case default
      name = ''
    end select
  end function canonical_agnes_method

  subroutine fail_hierarchy(result, status, message)
    type(hierarchy_result), intent(out) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    result%status = status
    result%message = message
  end subroutine fail_hierarchy

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    end do
  end function lower_ascii

end module cluster_hierarchy
