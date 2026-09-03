module grbase_arrays
  use grbase_types, only : set_list_t
  implicit none
  private

  public :: cell_to_entry
  public :: entry_to_cell
  public :: make_plevels
  public :: next_cell
  public :: next_cell_slice
  public :: slice_to_entries
  public :: cell_to_entry_perm
  public :: perm_cell_entries
  public :: choose_integer
  public :: combinations
  public :: next_combination_mask

contains

  pure function make_plevels(dim) result(plevels)
    integer, intent(in) :: dim(:) !! Dimension sizes; every element must be positive.
    integer :: plevels(size(dim))
    integer :: i

    if (size(dim) == 0) return
    plevels(1) = 1
    do i = 2, size(dim)
      plevels(i) = plevels(i - 1) * dim(i - 1)
    end do
  end function make_plevels

  pure integer function cell_to_entry(cell, dim) result(entry)
    integer, intent(in) :: cell(:) !! One-based cell coordinates, one coordinate per dimension.
    integer, intent(in) :: dim(:) !! Dimension sizes matching `cell`.
    integer :: i
    integer :: stride

    if (size(cell) /= size(dim)) then
      entry = 0
      return
    end if
    entry = 1
    stride = 1
    do i = 1, size(dim)
      if (cell(i) < 1 .or. cell(i) > dim(i)) then
        entry = 0
        return
      end if
      entry = entry + (cell(i) - 1) * stride
      stride = stride * dim(i)
    end do
  end function cell_to_entry

  pure function entry_to_cell(entry, dim) result(cell)
    integer, intent(in) :: entry !! One-based flattened entry index in Fortran/R column-major order.
    integer, intent(in) :: dim(:) !! Dimension sizes of the table or array.
    integer :: cell(size(dim))
    integer :: i
    integer :: rem
    integer :: total

    total = product(dim)
    if (size(dim) == 0) total = 1
    if (entry < 1 .or. entry > total) then
      cell = 0
      return
    end if
    rem = entry - 1
    do i = 1, size(dim)
      cell(i) = modulo(rem, dim(i)) + 1
      rem = rem / dim(i)
    end do
  end function entry_to_cell

  pure subroutine next_cell(cell, dim, wrapped)
    integer, intent(inout) :: cell(:) !! Current one-based cell; overwritten by the next cell.
    integer, intent(in) :: dim(:) !! Dimension sizes matching `cell`.
    logical, intent(out), optional :: wrapped !! True when increment wrapped from the last cell to the first.
    integer :: j
    logical :: did_wrap

    did_wrap = .true.
    if (size(cell) /= size(dim)) then
      if (present(wrapped)) wrapped = did_wrap
      return
    end if
    do j = 1, size(dim)
      if (cell(j) < dim(j)) then
        cell(j) = cell(j) + 1
        did_wrap = .false.
        exit
      else
        cell(j) = 1
      end if
    end do
    if (present(wrapped)) wrapped = did_wrap
  end subroutine next_cell

  pure subroutine next_cell_slice(cell, dim, fixed_dim, exhausted)
    integer, intent(inout) :: cell(:) !! Current cell; fixed coordinates are preserved and free coordinates are incremented.
    integer, intent(in) :: dim(:) !! Full array dimensions matching `cell`.
    integer, intent(in) :: fixed_dim(:) !! One-based dimension indices held fixed while advancing the slice.
    logical, intent(out) :: exhausted !! True when no further cell remains in the requested slice.
    logical :: fixed(size(dim))
    integer :: i
    integer :: j

    exhausted = .false.
    if (size(cell) /= size(dim)) then
      exhausted = .true.
      return
    end if
    fixed = .false.
    do i = 1, size(fixed_dim)
      if (fixed_dim(i) >= 1 .and. fixed_dim(i) <= size(dim)) fixed(fixed_dim(i)) = .true.
    end do
    do j = 1, size(dim)
      if (fixed(j)) cycle
      if (cell(j) < dim(j)) then
        cell(j) = cell(j) + 1
        return
      end if
      cell(j) = 1
    end do
    exhausted = .true.
  end subroutine next_cell_slice

  pure function slice_to_entries(slice_cell, slice_dim, dim) result(entries)
    integer, intent(in) :: slice_cell(:) !! Fixed one-based coordinate values, corresponding to `slice_dim`.
    integer, intent(in) :: slice_dim(:) !! One-based dimensions fixed to `slice_cell`.
    integer, intent(in) :: dim(:) !! Full array dimensions.
    integer, allocatable :: entries(:)
    integer, allocatable :: cell(:)
    integer :: i
    integer :: nfree
    integer :: nout
    logical :: exhausted

    if (size(slice_cell) /= size(slice_dim)) then
      allocate(entries(0))
      return
    end if
    nfree = product(dim)
    do i = 1, size(slice_dim)
      if (slice_dim(i) < 1 .or. slice_dim(i) > size(dim)) then
        allocate(entries(0))
        return
      end if
      if (slice_cell(i) < 1 .or. slice_cell(i) > dim(slice_dim(i))) then
        allocate(entries(0))
        return
      end if
      nfree = nfree / dim(slice_dim(i))
    end do
    allocate(entries(nfree), cell(size(dim)))
    cell = 1
    do i = 1, size(slice_dim)
      cell(slice_dim(i)) = slice_cell(i)
    end do
    nout = 0
    do
      nout = nout + 1
      entries(nout) = cell_to_entry(cell, dim)
      call next_cell_slice(cell, dim, slice_dim, exhausted)
      if (exhausted) exit
    end do
  end function slice_to_entries

  pure integer function cell_to_entry_perm(cell, dim, perm) result(entry)
    integer, intent(in) :: cell(:) !! One-based coordinates in the original dimension order.
    integer, intent(in) :: dim(:) !! Original dimension sizes.
    integer, intent(in) :: perm(:) !! Permutation giving the new dimension order.
    integer, allocatable :: pcell(:)
    integer, allocatable :: pdim(:)

    if (size(perm) /= size(dim) .or. size(cell) /= size(dim)) then
      entry = 0
      return
    end if
    allocate(pcell(size(cell)), pdim(size(dim)))
    pcell = cell(perm)
    pdim = dim(perm)
    entry = cell_to_entry(pcell, pdim)
  end function cell_to_entry_perm

  pure function perm_cell_entries(dim, perm) result(map)
    integer, intent(in) :: dim(:) !! Original dimension sizes.
    integer, intent(in) :: perm(:) !! Permutation defining the new dimension order.
    integer, allocatable :: map(:)
    integer, allocatable :: new_cell(:)
    integer, allocatable :: old_cell(:)
    integer, allocatable :: new_dim(:)
    integer :: e
    integer :: i
    integer :: n

    if (size(perm) /= size(dim)) then
      allocate(map(0))
      return
    end if
    if (any([(count(perm == i) /= 1, i = 1, size(dim))])) then
      allocate(map(0))
      return
    end if
    n = product(dim)
    if (size(dim) == 0) n = 1
    allocate(map(n), new_cell(size(dim)), old_cell(size(dim)), new_dim(size(dim)))
    new_dim = dim(perm)
    do e = 1, n
      new_cell = entry_to_cell(e, new_dim)
      do i = 1, size(dim)
        old_cell(perm(i)) = new_cell(i)
      end do
      map(e) = cell_to_entry(old_cell, dim)
    end do
  end function perm_cell_entries

  pure integer function choose_integer(n, k) result(value)
    integer, value :: n !! Number of available items; must be nonnegative for a nonzero result.
    integer, value :: k !! Number of items selected; values outside 0..n return zero.
    integer :: i
    integer :: kk

    if (n < 0 .or. k < 0 .or. k > n) then
      value = 0
      return
    end if
    kk = min(k, n - k)
    value = 1
    do i = 1, kk
      value = (value * (n - kk + i)) / i
    end do
  end function choose_integer

  pure function combinations(n, k) result(out)
    integer, value :: n !! Number of labeled items, represented by integers 1..n.
    integer, value :: k !! Combination size in the inclusive range 0..n.
    integer, allocatable :: out(:, :)
    integer, allocatable :: comb(:)
    integer :: col
    integer :: i
    integer :: nc

    nc = choose_integer(n, k)
    if (nc == 0 .and. k /= 0) then
      allocate(out(max(0, k), 0))
      return
    end if
    allocate(out(k, nc))
    if (k == 0) return
    allocate(comb(k))
    comb = [(i, i = 1, k)]
    col = 1
    do
      out(:, col) = comb
      if (comb(1) == n - k + 1 .and. all(comb == [(n - k + i, i = 1, k)])) exit
      i = k
      do while (i >= 1)
        if (comb(i) < n - k + i) exit
        i = i - 1
      end do
      if (i < 1) exit
      comb(i) = comb(i) + 1
      do while (i < k)
        i = i + 1
        comb(i) = comb(i - 1) + 1
      end do
      col = col + 1
    end do
  end function combinations

  pure subroutine next_combination_mask(mask, has_next)
    integer, intent(inout) :: mask(:) !! Zero-one selection mask; overwritten with the next same-cardinality mask.
    logical, intent(out) :: has_next !! True when a distinct next mask was produced.
    integer :: i
    integer :: j
    integer :: ones

    has_next = .false.
    if (size(mask) < 2) return
    do i = size(mask) - 1, 1, -1
      if (mask(i) == 1 .and. mask(i + 1) == 0) then
        mask(i) = 0
        mask(i + 1) = 1
        if (i + 1 < size(mask)) then
          ones = sum(mask(i + 2:))
          mask(i + 2:) = 0
          do j = i + 2, min(size(mask), i + 1 + ones)
            mask(j) = 1
          end do
        end if
        has_next = .true.
        return
      end if
    end do
  end subroutine next_combination_mask

end module grbase_arrays
