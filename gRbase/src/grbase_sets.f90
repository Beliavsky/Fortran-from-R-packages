module grbase_sets
  use grbase_arrays, only : combinations
  use grbase_types, only : integer_set_t, set_list_t
  implicit none
  private

  public :: unique_sorted
  public :: set_union
  public :: set_intersection
  public :: set_difference
  public :: is_subset_of
  public :: contains_set
  public :: all_subsets
  public :: maximal_sets
  public :: minimal_sets
  public :: all_pairs
  public :: append_set

contains

  pure function unique_sorted(x) result(out)
    integer, intent(in) :: x(:) !! Integer values to sort and deduplicate.
    integer, allocatable :: out(:)
    integer, allocatable :: work(:)
    integer :: i
    integer :: j
    integer :: n
    integer :: tmp

    if (size(x) == 0) then
      allocate(out(0))
      return
    end if
    allocate(work(size(x)))
    work = x
    do i = 2, size(work)
      tmp = work(i)
      j = i - 1
      do while (j >= 1)
        if (work(j) <= tmp) exit
        work(j + 1) = work(j)
        j = j - 1
      end do
      work(j + 1) = tmp
    end do
    n = 1
    do i = 2, size(work)
      if (work(i) /= work(n)) then
        n = n + 1
        work(n) = work(i)
      end if
    end do
    allocate(out(n))
    out = work(:n)
  end function unique_sorted

  pure function set_union(a, b) result(out)
    integer, intent(in) :: a(:) !! First integer set; duplicates and ordering are ignored.
    integer, intent(in) :: b(:) !! Second integer set; duplicates and ordering are ignored.
    integer, allocatable :: out(:)
    integer, allocatable :: both(:)

    allocate(both(size(a) + size(b)))
    if (size(a) > 0) both(:size(a)) = a
    if (size(b) > 0) both(size(a) + 1:) = b
    out = unique_sorted(both)
  end function set_union

  pure function set_intersection(a, b) result(out)
    integer, intent(in) :: a(:) !! First integer set.
    integer, intent(in) :: b(:) !! Second integer set.
    integer, allocatable :: out(:)
    integer, allocatable :: aa(:)
    integer, allocatable :: bb(:)
    integer, allocatable :: tmp(:)
    integer :: i
    integer :: j
    integer :: k

    aa = unique_sorted(a)
    bb = unique_sorted(b)
    allocate(tmp(min(size(aa), size(bb))))
    i = 1
    j = 1
    k = 0
    do while (i <= size(aa) .and. j <= size(bb))
      if (aa(i) == bb(j)) then
        k = k + 1
        tmp(k) = aa(i)
        i = i + 1
        j = j + 1
      else if (aa(i) < bb(j)) then
        i = i + 1
      else
        j = j + 1
      end if
    end do
    allocate(out(k))
    if (k > 0) out = tmp(:k)
  end function set_intersection

  pure function set_difference(a, b) result(out)
    integer, intent(in) :: a(:) !! Source integer set.
    integer, intent(in) :: b(:) !! Values removed from the source set.
    integer, allocatable :: out(:)
    integer, allocatable :: aa(:)
    integer, allocatable :: bb(:)
    integer, allocatable :: tmp(:)
    integer :: i
    integer :: j
    integer :: k

    aa = unique_sorted(a)
    bb = unique_sorted(b)
    allocate(tmp(size(aa)))
    j = 1
    k = 0
    do i = 1, size(aa)
      do while (j <= size(bb))
        if (bb(j) >= aa(i)) exit
        j = j + 1
      end do
      if (j <= size(bb)) then
        if (bb(j) == aa(i)) cycle
      end if
      k = k + 1
      tmp(k) = aa(i)
    end do
    allocate(out(k))
    if (k > 0) out = tmp(:k)
  end function set_difference

  pure logical function is_subset_of(a, b) result(ok)
    integer, intent(in) :: a(:) !! Candidate subset.
    integer, intent(in) :: b(:) !! Candidate superset.
    integer, allocatable :: missing(:)

    missing = set_difference(a, b)
    ok = size(missing) == 0
  end function is_subset_of

  pure logical function contains_set(list, values) result(ok)
    type(set_list_t), intent(in) :: list !! List of integer sets to search.
    integer, intent(in) :: values(:) !! Set that must be contained in at least one list element.
    integer :: i

    ok = .false.
    do i = 1, list%count
      if (is_subset_of(values, list%set(i)%value)) then
        ok = .true.
        return
      end if
    end do
  end function contains_set

  pure subroutine append_set(list, values)
    type(set_list_t), intent(inout) :: list !! Set list receiving one appended set.
    integer, intent(in) :: values(:) !! Integer values stored as a sorted unique set.
    type(integer_set_t), allocatable :: tmp(:)
    integer :: n

    n = list%count
    allocate(tmp(n + 1))
    if (n > 0) tmp(:n) = list%set
    tmp(n + 1)%value = unique_sorted(values)
    call move_alloc(tmp, list%set)
    list%count = n + 1
  end subroutine append_set

  pure function all_subsets(values, min_size, max_size) result(out)
    integer, intent(in) :: values(:) !! Values from which subsets are formed.
    integer, intent(in), optional :: min_size !! Minimum subset cardinality; defaults to zero.
    integer, intent(in), optional :: max_size !! Maximum subset cardinality; defaults to all values.
    type(set_list_t) :: out
    integer, allocatable :: vals(:)
    integer, allocatable :: cmb(:, :)
    integer :: i
    integer :: k
    integer :: lo
    integer :: hi

    vals = unique_sorted(values)
    lo = 0
    hi = size(vals)
    if (present(min_size)) lo = max(0, min_size)
    if (present(max_size)) hi = min(size(vals), max_size)
    allocate(out%set(0))
    if (hi < lo) return
    do k = lo, hi
      cmb = combinations(size(vals), k)
      do i = 1, size(cmb, 2)
        call append_set(out, vals(cmb(:, i)))
      end do
    end do
  end function all_subsets

  pure function maximal_sets(list) result(out)
    type(set_list_t), intent(in) :: list !! Input collection whose inclusion-maximal sets are requested.
    type(set_list_t) :: out
    logical, allocatable :: keep(:)
    integer :: i
    integer :: j

    allocate(out%set(0), keep(list%count))
    keep = .true.
    do i = 1, list%count
      do j = 1, list%count
        if (i == j) cycle
        if (is_subset_of(list%set(i)%value, list%set(j)%value)) then
          if (size(list%set(i)%value) < size(list%set(j)%value) .or. j < i) then
            keep(i) = .false.
            exit
          end if
        end if
      end do
    end do
    do i = 1, list%count
      if (keep(i)) call append_set(out, list%set(i)%value)
    end do
  end function maximal_sets

  pure function minimal_sets(list) result(out)
    type(set_list_t), intent(in) :: list !! Input collection whose inclusion-minimal sets are requested.
    type(set_list_t) :: out
    logical, allocatable :: keep(:)
    integer :: i
    integer :: j

    allocate(out%set(0), keep(list%count))
    keep = .true.
    do i = 1, list%count
      do j = 1, list%count
        if (i == j) cycle
        if (is_subset_of(list%set(j)%value, list%set(i)%value)) then
          if (size(list%set(j)%value) < size(list%set(i)%value) .or. j < i) then
            keep(i) = .false.
            exit
          end if
        end if
      end do
    end do
    do i = 1, list%count
      if (keep(i)) call append_set(out, list%set(i)%value)
    end do
  end function minimal_sets

  pure function all_pairs(x, y) result(pairs)
    integer, intent(in) :: x(:) !! First collection of labels; when `y` is absent, unordered pairs within `x` are returned.
    integer, intent(in), optional :: y(:) !! Optional second collection; when present, the Cartesian pair set is returned.
    integer, allocatable :: pairs(:, :)
    integer :: i
    integer :: j
    integer :: k
    integer :: n

    if (present(y)) then
      n = size(x) * size(y)
      allocate(pairs(2, n))
      k = 0
      do j = 1, size(y)
        do i = 1, size(x)
          k = k + 1
          pairs(:, k) = [x(i), y(j)]
        end do
      end do
    else
      n = size(x) * max(0, size(x) - 1) / 2
      allocate(pairs(2, n))
      k = 0
      do i = 1, size(x) - 1
        do j = i + 1, size(x)
          k = k + 1
          pairs(:, k) = [x(i), x(j)]
        end do
      end do
    end if
  end function all_pairs

end module grbase_sets
