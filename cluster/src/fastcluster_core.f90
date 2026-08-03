! SPDX-License-Identifier: BSD-2-Clause
module fastcluster_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_positive_inf
  use fastcluster_kinds, only: dp
  use fastcluster_types, only: hclust_result, fc_success, fc_invalid_argument, &
    fc_nan_distance, fc_numerical_failure, fc_allocation_failure
  use fastcluster_distances, only: pairwise_distances, condensed_to_matrix
  implicit none
  private

  integer, parameter :: method_single = 1
  integer, parameter :: method_complete = 2
  integer, parameter :: method_average = 3
  integer, parameter :: method_weighted = 4
  integer, parameter :: method_ward = 5
  integer, parameter :: method_centroid = 6
  integer, parameter :: method_median = 7

  public :: hclust
  public :: hclust_matrix
  public :: hclust_condensed
  public :: hclust_vector

  interface hclust
    module procedure hclust_matrix
    module procedure hclust_condensed
  end interface hclust

contains

  subroutine hclust_matrix(distances, method, result, members)
    real(dp), intent(in) :: distances(:, :)
    character(len=*), intent(in) :: method
    type(hclust_result), intent(out) :: result
    real(dp), intent(in), optional :: members(:)

    real(dp), allocatable :: work(:, :)
    character(len=:), allocatable :: canonical
    integer :: code, n, stat
    logical :: square_input, sqrt_output

    call initialize_result(result)
    n = size(distances, 1)
    if (n < 2 .or. size(distances, 2) /= n) then
      call fail_result(result, fc_invalid_argument, &
        'distances must be a square matrix of order at least two')
      return
    end if

    call parse_matrix_method(method, code, canonical, square_input, sqrt_output, result)
    if (.not. result%ok()) return

    allocate(work(n, n), stat=stat)
    if (stat /= 0) then
      call fail_result(result, fc_allocation_failure, 'could not allocate working distances')
      return
    end if
    call copy_and_validate_distances(distances, work, result)
    if (.not. result%ok()) return
    if (square_input) work = work * work

    call cluster_engine(work, code, canonical, result, members, sqrt_output)
  end subroutine hclust_matrix

  subroutine hclust_condensed(d, n, method, result, members)
    real(dp), intent(in) :: d(:)
    integer, intent(in) :: n
    character(len=*), intent(in) :: method
    type(hclust_result), intent(out) :: result
    real(dp), intent(in), optional :: members(:)

    real(dp), allocatable :: distances(:, :)
    character(len=:), allocatable :: message
    integer :: status

    call condensed_to_matrix(d, n, distances, status, message)
    if (status /= fc_success) then
      call initialize_result(result)
      call fail_result(result, status, message)
      return
    end if
    if (present(members)) then
      call hclust_matrix(distances, method, result, members)
    else
      call hclust_matrix(distances, method, result)
    end if
  end subroutine hclust_condensed

  subroutine hclust_vector(x, method, result, members, metric, p)
    real(dp), intent(in) :: x(:, :)
    character(len=*), intent(in) :: method
    type(hclust_result), intent(out) :: result
    real(dp), intent(in), optional :: members(:)
    character(len=*), intent(in), optional :: metric
    real(dp), intent(in), optional :: p

    real(dp), allocatable :: distances(:, :)
    character(len=:), allocatable :: canonical, message, metric_name
    integer :: code, status
    logical :: squared, sqrt_output

    call initialize_result(result)
    metric_name = 'euclidean'
    if (present(metric)) metric_name = lower_ascii(trim(adjustl(metric)))
    call parse_vector_method(method, code, canonical, squared, sqrt_output, result)
    if (.not. result%ok()) return
    if (code /= method_single .and. metric_name /= 'euclidean') then
      call fail_result(result, fc_invalid_argument, &
        'ward, centroid, and median vector methods require Euclidean distances')
      return
    end if

    if (present(p)) then
      call pairwise_distances(x, metric_name, distances, p=p, status=status, &
        message=message, squared_euclidean=squared)
    else
      call pairwise_distances(x, metric_name, distances, status=status, &
        message=message, squared_euclidean=squared)
    end if
    if (status /= fc_success) then
      call fail_result(result, status, message)
      return
    end if

    if (present(members)) then
      call cluster_engine(distances, code, canonical, result, members, sqrt_output)
    else
      call cluster_engine(distances, code, canonical, result, sqrt_output=sqrt_output)
    end if
    if (result%ok()) result%metric = metric_name
  end subroutine hclust_vector

  subroutine cluster_engine(initial_distances, code, canonical, result, members, sqrt_output)
    real(dp), intent(in) :: initial_distances(:, :)
    integer, intent(in) :: code
    character(len=*), intent(in) :: canonical
    type(hclust_result), intent(inout) :: result
    real(dp), intent(in), optional :: members(:)
    logical, intent(in) :: sqrt_output

    real(dp), allocatable :: distances(:, :), sizes(:)
    logical, allocatable :: active(:)
    real(dp) :: da, dab, db, minimum, new_distance, sa, sb, sc
    integer :: a, b, c, i, j, k, max_nodes, n, new_node, stat
    logical :: found

    n = size(initial_distances, 1)
    max_nodes = 2 * n - 1
    allocate(distances(max_nodes, max_nodes), sizes(max_nodes), active(max_nodes), stat=stat)
    if (stat /= 0) then
      call fail_result(result, fc_allocation_failure, 'could not allocate clustering workspace')
      return
    end if
    distances = ieee_value(0.0_dp, ieee_positive_inf)
    sizes = 0.0_dp
    active = .false.
    distances(1:n, 1:n) = initial_distances
    do i = 1, n
      distances(i, i) = 0.0_dp
    end do
    active(1:n) = .true.
    sizes(1:n) = 1.0_dp
    if (present(members)) then
      if (size(members) /= n .or. any(members <= 0.0_dp) .or. &
          any(.not. ieee_is_finite(members))) then
        call fail_result(result, fc_invalid_argument, &
          'members must contain one positive finite value per observation')
        return
      end if
      sizes(1:n) = members
    end if

    allocate(result%merge(n - 1, 2), result%height(n - 1), result%order(n), stat=stat)
    if (stat /= 0) then
      call fail_result(result, fc_allocation_failure, 'could not allocate clustering result')
      return
    end if
    result%n = n
    result%method = canonical
    result%metric = 'precomputed'

    do k = 1, n - 1
      minimum = ieee_value(0.0_dp, ieee_positive_inf)
      found = .false.
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
      if (.not. found .or. ieee_is_nan(minimum)) then
        call fail_result(result, fc_nan_distance, 'NaN dissimilarity value')
        return
      end if

      result%merge(k, 1) = output_node(a, n)
      result%merge(k, 2) = output_node(b, n)
      if (sqrt_output) then
        if (minimum < -100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(minimum))) then
          call fail_result(result, fc_numerical_failure, &
            'a squared linkage distance became negative')
          return
        end if
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
        case (method_single)
          new_distance = min(da, db)
        case (method_complete)
          new_distance = max(da, db)
        case (method_average)
          new_distance = (sa * da + sb * db) / (sa + sb)
        case (method_weighted)
          new_distance = 0.5_dp * (da + db)
        case (method_ward, method_centroid, method_median)
          if (.not. ieee_is_finite(da) .or. .not. ieee_is_finite(db) .or. &
              .not. ieee_is_finite(dab)) then
            new_distance = ieee_value(0.0_dp, ieee_positive_inf)
          else if (code == method_ward) then
            new_distance = ((sc + sa) * da - sc * dab + (sc + sb) * db) / &
              (sa + sb + sc)
          else if (code == method_centroid) then
            new_distance = sa * da / (sa + sb) + sb * db / (sa + sb) - &
              sa * sb * dab / (sa + sb) ** 2
          else
            new_distance = 0.5_dp * (da + db) - 0.25_dp * dab
          end if
        case default
          call fail_result(result, fc_invalid_argument, 'invalid linkage method code')
          return
        end select
        if (ieee_is_nan(new_distance)) then
          call fail_result(result, fc_numerical_failure, &
            'NaN dissimilarity value in an intermediate result')
          return
        end if
        distances(c, new_node) = new_distance
        distances(new_node, c) = new_distance
      end do
      sizes(new_node) = sa + sb
      active(a) = .false.
      active(b) = .false.
      active(new_node) = .true.
    end do

    call make_order(result%merge, result%order)
    result%status = fc_success
    result%message = 'ok'
  end subroutine cluster_engine

  subroutine copy_and_validate_distances(input, output, result)
    real(dp), intent(in) :: input(:, :)
    real(dp), intent(out) :: output(:, :)
    type(hclust_result), intent(inout) :: result

    real(dp) :: scale
    integer :: i, j, n

    n = size(input, 1)
    output = 0.0_dp
    do j = 1, n - 1
      do i = j + 1, n
        if (ieee_is_nan(input(i, j)) .or. ieee_is_nan(input(j, i))) then
          call fail_result(result, fc_nan_distance, 'NaN dissimilarity value')
          return
        end if
        if (ieee_is_finite(input(i, j)) .and. ieee_is_finite(input(j, i))) then
          scale = max(1.0_dp, abs(input(i, j)), abs(input(j, i)))
          if (abs(input(i, j) - input(j, i)) > 100.0_dp * epsilon(1.0_dp) * scale) then
            call fail_result(result, fc_invalid_argument, 'distance matrix must be symmetric')
            return
          end if
        else if (ieee_is_finite(input(i, j)) .neqv. ieee_is_finite(input(j, i))) then
          call fail_result(result, fc_invalid_argument, 'distance matrix must be symmetric')
          return
        end if
        if (input(i, j) < 0.0_dp .or. input(j, i) < 0.0_dp) then
          call fail_result(result, fc_invalid_argument, 'dissimilarities must be nonnegative')
          return
        end if
        output(i, j) = input(i, j)
        output(j, i) = input(i, j)
      end do
    end do
  end subroutine copy_and_validate_distances

  subroutine make_order(merge, order)
    integer, intent(in) :: merge(:, :)
    integer, intent(out) :: order(:)

    integer, allocatable :: stack(:)
    integer :: code, count, n, top

    n = size(order)
    allocate(stack(2 * n))
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

  integer pure function output_node(node, n) result(code)
    integer, intent(in) :: node, n

    if (node <= n) then
      code = -node
    else
      code = node - n
    end if
  end function output_node

  subroutine parse_matrix_method(name, code, canonical, square_input, sqrt_output, result)
    character(len=*), intent(in) :: name
    integer, intent(out) :: code
    character(len=:), allocatable, intent(out) :: canonical
    logical, intent(out) :: square_input, sqrt_output
    type(hclust_result), intent(inout) :: result

    character(len=:), allocatable :: value

    value = lower_ascii(trim(adjustl(name)))
    square_input = .false.
    sqrt_output = .false.
    select case (value)
    case ('single')
      code = method_single
      canonical = 'single'
    case ('complete')
      code = method_complete
      canonical = 'complete'
    case ('average')
      code = method_average
      canonical = 'average'
    case ('mcquitty', 'weighted')
      code = method_weighted
      canonical = 'mcquitty'
    case ('ward', 'ward.d')
      code = method_ward
      canonical = 'ward.D'
    case ('ward.d2')
      code = method_ward
      canonical = 'ward.D2'
      square_input = .true.
      sqrt_output = .true.
    case ('centroid')
      code = method_centroid
      canonical = 'centroid'
    case ('median')
      code = method_median
      canonical = 'median'
    case default
      code = 0
      canonical = ''
      call fail_result(result, fc_invalid_argument, 'invalid clustering method: '//value)
    end select
  end subroutine parse_matrix_method

  subroutine parse_vector_method(name, code, canonical, squared, sqrt_output, result)
    character(len=*), intent(in) :: name
    integer, intent(out) :: code
    character(len=:), allocatable, intent(out) :: canonical
    logical, intent(out) :: squared, sqrt_output
    type(hclust_result), intent(inout) :: result

    character(len=:), allocatable :: value

    value = lower_ascii(trim(adjustl(name)))
    squared = .false.
    sqrt_output = .false.
    select case (value)
    case ('single')
      code = method_single
      canonical = 'single'
    case ('ward')
      code = method_ward
      canonical = 'ward'
      squared = .true.
      sqrt_output = .true.
    case ('centroid')
      code = method_centroid
      canonical = 'centroid'
      squared = .true.
      sqrt_output = .true.
    case ('median')
      code = method_median
      canonical = 'median'
      squared = .true.
      sqrt_output = .true.
    case default
      code = 0
      canonical = ''
      call fail_result(result, fc_invalid_argument, &
        'invalid vector clustering method: '//value)
    end select
  end subroutine parse_vector_method

  subroutine initialize_result(result)
    type(hclust_result), intent(out) :: result

    result%n = 0
    result%status = fc_success
    result%message = 'ok'
    result%method = ''
    result%metric = ''
  end subroutine initialize_result

  subroutine fail_result(result, status, message)
    type(hclust_result), intent(inout) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message

    result%status = status
    result%message = message
  end subroutine fail_result

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

end module fastcluster_core
