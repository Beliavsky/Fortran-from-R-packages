! SPDX-License-Identifier: GPL-2.0-only
module hierportfolios_hierarchy
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use hierportfolios_kinds, only: dp
  use hierportfolios_types, only: hierarchy_result, hp_success, hp_invalid_argument, &
    hp_numerical_failure
  implicit none
  private

  integer, parameter :: link_single = 1
  integer, parameter :: link_complete = 2
  integer, parameter :: link_average = 3
  integer, parameter :: link_ward = 4

  public :: covariance_distance
  public :: euclidean_distances
  public :: hclust_distance
  public :: diana_distance
  public :: cutree
  public :: select_gap_clusters

contains

  subroutine covariance_distance(covar, features, distances, status, message)
    real(dp), intent(in) :: covar(:, :)
    real(dp), allocatable, intent(out) :: features(:, :), distances(:, :)
    integer, intent(out) :: status
    character(len=:), allocatable, intent(out) :: message

    real(dp), allocatable :: corr(:, :)
    real(dp) :: scale, tol
    integer :: i, j, n

    n = size(covar, 1)
    if (n < 2 .or. size(covar, 2) /= n) then
      call fail(status, message, hp_invalid_argument, &
        'covariance matrix must be square with order at least two')
      allocate(features(0, 0), distances(0, 0))
      return
    end if
    if (any(.not. ieee_is_finite(covar)) .or. any([(covar(i, i) <= 0.0_dp, i=1,n)])) then
      call fail(status, message, hp_invalid_argument, &
        'covariance matrix must be finite with positive diagonal')
      allocate(features(0, 0), distances(0, 0))
      return
    end if
    tol = 1000.0_dp * epsilon(1.0_dp)
    do j = 1, n - 1
      do i = j + 1, n
        scale = max(1.0_dp, abs(covar(i, j)), abs(covar(j, i)))
        if (abs(covar(i, j) - covar(j, i)) > tol * scale) then
          call fail(status, message, hp_invalid_argument, &
            'covariance matrix must be symmetric')
          allocate(features(0, 0), distances(0, 0))
          return
        end if
      end do
    end do

    allocate(corr(n, n), features(n, n))
    do j = 1, n
      do i = 1, n
        corr(i, j) = covar(i, j) / sqrt(covar(i, i) * covar(j, j))
        corr(i, j) = max(-1.0_dp, min(1.0_dp, corr(i, j)))
        features(i, j) = sqrt(max(0.0_dp, 0.5_dp * (1.0_dp - corr(i, j))))
      end do
    end do
    call euclidean_distances(features, distances)
    status = hp_success
    message = 'ok'
  end subroutine covariance_distance

  subroutine euclidean_distances(x, distances)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: distances(:, :)

    integer :: i, j, n

    n = size(x, 1)
    allocate(distances(n, n))
    distances = 0.0_dp
    do j = 1, n - 1
      do i = j + 1, n
        distances(i, j) = sqrt(sum((x(i, :) - x(j, :)) ** 2))
        distances(j, i) = distances(i, j)
      end do
    end do
  end subroutine euclidean_distances

  subroutine hclust_distance(initial, method, result)
    real(dp), intent(in) :: initial(:, :)
    character(len=*), intent(in) :: method
    type(hierarchy_result), intent(out) :: result

    real(dp), allocatable :: distances(:, :), sizes(:)
    logical, allocatable :: active(:)
    real(dp) :: da, dab, db, minimum, new_distance, sa, sb, sc
    integer :: a, b, c, code, i, j, k, max_nodes, n, new_node
    logical :: found, square_input, sqrt_output
    character(len=:), allocatable :: canonical

    call init_hierarchy(result)
    n = size(initial, 1)
    if (n < 2 .or. size(initial, 2) /= n .or. any(.not. ieee_is_finite(initial)) .or. &
        any(initial < 0.0_dp)) then
      call fail_hierarchy(result, hp_invalid_argument, &
        'distance matrix must be finite, nonnegative, and square')
      return
    end if
    call parse_method(method, code, canonical, square_input, sqrt_output)
    if (code == 0) then
      call fail_hierarchy(result, hp_invalid_argument, &
        'linkage must be single, complete, average, or ward')
      return
    end if

    max_nodes = 2 * n - 1
    allocate(distances(max_nodes, max_nodes), sizes(max_nodes), active(max_nodes))
    distances = huge(1.0_dp)
    distances(1:n, 1:n) = initial
    if (square_input) distances(1:n, 1:n) = distances(1:n, 1:n) ** 2
    do i = 1, n
      distances(i, i) = 0.0_dp
    end do
    sizes = 0.0_dp
    sizes(1:n) = 1.0_dp
    active = .false.
    active(1:n) = .true.
    allocate(result%merge(n - 1, 2), result%height(n - 1), result%order(n))
    result%n = n
    result%method = canonical

    do k = 1, n - 1
      found = .false.
      minimum = huge(1.0_dp)
      a = 0
      b = 0
      do i = 1, n + k - 2
        if (.not. active(i)) cycle
        do j = i + 1, n + k - 1
          if (.not. active(j)) cycle
          if (.not. found .or. distances(i, j) < minimum) then
            found = .true.
            minimum = distances(i, j)
            a = i
            b = j
          end if
        end do
      end do
      if (.not. found) then
        call fail_hierarchy(result, hp_numerical_failure, 'clustering merge failed')
        return
      end if
      result%merge(k, 1) = output_node(a, n)
      result%merge(k, 2) = output_node(b, n)
      if (sqrt_output) then
        result%height(k) = sqrt(max(0.0_dp, minimum))
      else
        result%height(k) = minimum
      end if

      new_node = n + k
      sa = sizes(a)
      sb = sizes(b)
      dab = minimum
      do c = 1, new_node - 1
        if (.not. active(c) .or. c == a .or. c == b) cycle
        da = distances(min(a, c), max(a, c))
        db = distances(min(b, c), max(b, c))
        sc = sizes(c)
        select case (code)
        case (link_single)
          new_distance = min(da, db)
        case (link_complete)
          new_distance = max(da, db)
        case (link_average)
          new_distance = (sa * da + sb * db) / (sa + sb)
        case (link_ward)
          new_distance = ((sc + sa) * da - sc * dab + (sc + sb) * db) / &
            (sa + sb + sc)
        end select
        distances(c, new_node) = max(0.0_dp, new_distance)
        distances(new_node, c) = distances(c, new_node)
      end do
      sizes(new_node) = sa + sb
      active(a) = .false.
      active(b) = .false.
      active(new_node) = .true.
    end do
    call make_order(result%merge, result%order)
    result%status = hp_success
    result%message = 'ok'
  end subroutine hclust_distance

  subroutine diana_distance(distances, result)
    real(dp), intent(in) :: distances(:, :)
    type(hierarchy_result), intent(out) :: result

    integer, allocatable :: child1(:), child2(:), items(:)
    real(dp), allocatable :: node_height(:)
    integer :: counter, n, root

    call init_hierarchy(result)
    n = size(distances, 1)
    if (n < 2 .or. size(distances, 2) /= n .or. any(.not. ieee_is_finite(distances)) .or. &
        any(distances < 0.0_dp)) then
      call fail_hierarchy(result, hp_invalid_argument, &
        'distance matrix must be finite, nonnegative, and square')
      return
    end if
    allocate(child1(2*n-1), child2(2*n-1), node_height(2*n-1), items(n))
    child1 = 0
    child2 = 0
    node_height = 0.0_dp
    items = [(counter, counter=1,n)]
    counter = n
    call build_divisive_tree(distances, items, counter, root, child1, child2, node_height)
    if (root /= 2*n-1) then
      call fail_hierarchy(result, hp_numerical_failure, 'DIANA tree construction failed')
      return
    end if
    allocate(result%merge(n-1, 2), result%height(n-1), result%order(n))
    call tree_to_merge(n, child1, child2, node_height, result%merge, result%height)
    counter = 0
    call tree_order(root, n, child1, child2, result%order, counter)
    result%n = n
    result%method = 'diana'
    result%status = hp_success
    result%message = 'ok'
  end subroutine diana_distance

  subroutine cutree(hierarchy, k_clusters, labels, status)
    type(hierarchy_result), intent(in) :: hierarchy
    integer, intent(in) :: k_clusters
    integer, allocatable, intent(out) :: labels(:)
    integer, intent(out), optional :: status

    integer, allocatable :: parent(:), node_rep(:), root_label(:)
    integer :: a, b, i, label, n, s, ra, rb, root

    n = hierarchy%n
    if (.not. hierarchy%ok() .or. k_clusters < 1 .or. k_clusters > n) then
      allocate(labels(0))
      if (present(status)) status = hp_invalid_argument
      return
    end if
    allocate(labels(n), parent(n), node_rep(2*n-1), root_label(n))
    parent = [(i, i=1,n)]
    node_rep = 0
    node_rep(1:n) = [(i, i=1,n)]
    do s = 1, n - k_clusters
      a = hierarchy%merge(s, 1)
      b = hierarchy%merge(s, 2)
      if (a < 0) then
        ra = -a
      else
        ra = node_rep(n + a)
      end if
      if (b < 0) then
        rb = -b
      else
        rb = node_rep(n + b)
      end if
      ra = find_root(parent, ra)
      rb = find_root(parent, rb)
      if (ra /= rb) parent(rb) = ra
      node_rep(n + s) = ra
    end do
    root_label = 0
    label = 0
    do i = 1, n
      root = find_root(parent, i)
      if (root_label(root) == 0) then
        label = label + 1
        root_label(root) = label
      end if
      labels(i) = root_label(root)
    end do
    if (present(status)) status = hp_success
  end subroutine cutree

  subroutine select_gap_clusters(x, method, k_max, b_references, selected, gap, se, seed, status)
    real(dp), intent(in) :: x(:, :)
    character(len=*), intent(in) :: method
    integer, intent(in) :: k_max, b_references
    integer, intent(out) :: selected
    real(dp), allocatable, intent(out) :: gap(:), se(:)
    integer, intent(in), optional :: seed
    integer, intent(out), optional :: status

    type(hierarchy_result) :: hierarchy
    real(dp), allocatable :: d(:, :), reference(:, :), log_ref(:, :), log_w(:)
    real(dp), allocatable :: xmin(:), xmax(:)
    integer, allocatable :: labels(:)
    integer(int64) :: state
    real(dp) :: value, u
    integer :: b, i, j, k, km, n, p

    n = size(x, 1)
    p = size(x, 2)
    km = min(k_max, n)
    if (n < 2 .or. p < 1 .or. km < 1 .or. b_references < 1) then
      allocate(gap(0), se(0))
      selected = 0
      if (present(status)) status = hp_invalid_argument
      return
    end if
    allocate(gap(km), se(km), log_w(km), log_ref(b_references, km))
    allocate(reference(n, p), xmin(p), xmax(p))
    do j = 1, p
      xmin(j) = minval(x(:, j))
      xmax(j) = maxval(x(:, j))
    end do
    call euclidean_distances(x, d)
    call hclust_distance(d, method, hierarchy)
    if (.not. hierarchy%ok()) then
      selected = 0
      gap = 0.0_dp
      se = 0.0_dp
      if (present(status)) status = hierarchy%status
      return
    end if
    do k = 1, km
      call cutree(hierarchy, k, labels)
      value = within_dispersion(d, labels)
      log_w(k) = log(max(value, tiny(1.0_dp)))
    end do

    state = 104729_int64
    if (present(seed)) state = int(max(1, seed), int64)
    do b = 1, b_references
      do j = 1, p
        do i = 1, n
          u = uniform01(state)
          reference(i, j) = xmin(j) + u * (xmax(j) - xmin(j))
        end do
      end do
      call euclidean_distances(reference, d)
      call hclust_distance(d, method, hierarchy)
      if (.not. hierarchy%ok()) then
        selected = 0
        gap = 0.0_dp
        se = 0.0_dp
        if (present(status)) status = hierarchy%status
        return
      end if
      do k = 1, km
        call cutree(hierarchy, k, labels)
        value = within_dispersion(d, labels)
        log_ref(b, k) = log(max(value, tiny(1.0_dp)))
      end do
    end do
    do k = 1, km
      gap(k) = sum(log_ref(:, k)) / real(b_references, dp) - log_w(k)
      if (b_references > 1) then
        value = sum((log_ref(:, k) - sum(log_ref(:, k)) / real(b_references, dp)) ** 2)
        se(k) = sqrt(value / real(b_references - 1, dp)) * &
          sqrt(1.0_dp + 1.0_dp / real(b_references, dp))
      else
        se(k) = 0.0_dp
      end if
    end do
    selected = km
    do k = 1, km - 1
      if (gap(k) >= gap(k + 1) - se(k + 1)) then
        selected = k
        exit
      end if
    end do
    if (present(status)) status = hp_success
  end subroutine select_gap_clusters

  real(dp) function within_dispersion(distances, labels) result(value)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(in) :: labels(:)

    real(dp) :: cluster_sum
    integer :: c, i, j, nc

    value = 0.0_dp
    do c = 1, maxval(labels)
      nc = count(labels == c)
      if (nc <= 1) cycle
      cluster_sum = 0.0_dp
      do j = 1, size(labels) - 1
        if (labels(j) /= c) cycle
        do i = j + 1, size(labels)
          if (labels(i) == c) cluster_sum = cluster_sum + distances(i, j)
        end do
      end do
      value = value + cluster_sum / real(nc, dp)
    end do
  end function within_dispersion

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
    real(dp) :: average, best_average, difference, best_difference
    real(dp) :: to_remainder, to_splinter
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
    allocate(left(count(in_splinter)), right(m - count(in_splinter)))
    left = pack(items, in_splinter)
    right = pack(items, .not. in_splinter)
  end subroutine diana_split

  subroutine tree_to_merge(n, child1, child2, node_height, merge, height)
    integer, intent(in) :: n
    integer, intent(in) :: child1(:), child2(:)
    real(dp), intent(in) :: node_height(:)
    integer, intent(out) :: merge(:, :)
    real(dp), intent(out) :: height(:)

    integer :: node, s

    do s = 1, n - 1
      node = n + s
      merge(s, 1) = output_node(child1(node), n)
      merge(s, 2) = output_node(child2(node), n)
      height(s) = node_height(node)
    end do
  end subroutine tree_to_merge

  recursive subroutine tree_order(node, n, child1, child2, order, counter)
    integer, intent(in) :: node, n
    integer, intent(in) :: child1(:), child2(:)
    integer, intent(inout) :: order(:), counter

    if (node <= n) then
      counter = counter + 1
      order(counter) = node
    else
      call tree_order(child1(node), n, child1, child2, order, counter)
      call tree_order(child2(node), n, child1, child2, order, counter)
    end if
  end subroutine tree_order

  subroutine make_order(merge, order)
    integer, intent(in) :: merge(:, :)
    integer, intent(out) :: order(:)

    integer, allocatable :: stack(:)
    integer :: code, count, n, top

    n = size(order)
    allocate(stack(2*n))
    top = 1
    stack(top) = n - 1
    count = 0
    do while (top > 0)
      code = stack(top)
      top = top - 1
      if (code < 0) then
        count = count + 1
        order(count) = -code
      else
        top = top + 1
        stack(top) = merge(code, 2)
        top = top + 1
        stack(top) = merge(code, 1)
      end if
    end do
  end subroutine make_order

  integer function find_root(parent, item) result(root)
    integer, intent(inout) :: parent(:)
    integer, intent(in) :: item
    integer :: node, next

    node = item
    do while (parent(node) /= node)
      node = parent(node)
    end do
    root = node
    node = item
    do while (parent(node) /= node)
      next = parent(node)
      parent(node) = root
      node = next
    end do
  end function find_root

  integer pure function output_node(node, n) result(code)
    integer, intent(in) :: node, n
    if (node <= n) then
      code = -node
    else
      code = node - n
    end if
  end function output_node

  subroutine parse_method(name, code, canonical, square_input, sqrt_output)
    character(len=*), intent(in) :: name
    integer, intent(out) :: code
    character(len=:), allocatable, intent(out) :: canonical
    logical, intent(out) :: square_input, sqrt_output

    character(len=:), allocatable :: value

    value = lower_ascii(trim(adjustl(name)))
    square_input = .false.
    sqrt_output = .false.
    select case (value)
    case ('single')
      code = link_single
      canonical = 'single'
    case ('complete')
      code = link_complete
      canonical = 'complete'
    case ('average')
      code = link_average
      canonical = 'average'
    case ('ward', 'ward.d2')
      code = link_ward
      canonical = 'ward.D2'
      square_input = .true.
      sqrt_output = .true.
    case default
      code = 0
      canonical = ''
    end select
  end subroutine parse_method

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

  real(dp) function uniform01(state) result(u)
    integer(int64), intent(inout) :: state
    integer(int64), parameter :: a = 48271_int64
    integer(int64), parameter :: m = 2147483647_int64

    state = modulo(a * state, m)
    if (state <= 0_int64) state = state + m - 1_int64
    u = real(state, dp) / real(m, dp)
  end function uniform01

  subroutine init_hierarchy(result)
    type(hierarchy_result), intent(out) :: result
    result%status = hp_success
    result%message = 'ok'
    result%method = ''
    result%n = 0
  end subroutine init_hierarchy

  subroutine fail_hierarchy(result, status, message)
    type(hierarchy_result), intent(inout) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    result%status = status
    result%message = message
  end subroutine fail_hierarchy

  subroutine fail(status, message, code, text)
    integer, intent(out) :: status
    character(len=:), allocatable, intent(out) :: message
    integer, intent(in) :: code
    character(len=*), intent(in) :: text
    status = code
    message = text
  end subroutine fail

end module hierportfolios_hierarchy
