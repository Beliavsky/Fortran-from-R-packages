! SPDX-License-Identifier: GPL-2.0-only
module multcomp_cld
  use multcomp_types, only : cld_type
  implicit none
  private

  public :: compact_letter_display

contains

  subroutine compact_letter_display(significant, group_a, group_b, nlevels, result, decreasing)
    logical, intent(in) :: significant(:) !! Pairwise significance flags; true means the two groups must not share a letter.
    integer, intent(in) :: group_a(:) !! One-based first group index for each pairwise comparison.
    integer, intent(in) :: group_b(:) !! One-based second group index for each pairwise comparison.
    integer, intent(in) :: nlevels !! Number of treatment or factor levels represented by the comparisons.
    type(cld_type), intent(out) :: result !! Piepho insert-absorb compact-letter display and logical letter matrix.
    logical, intent(in), optional :: decreasing !! If true, reverse the final letter-column ordering convention.

    logical, allocatable :: matrix(:, :)
    logical, allocatable :: duplicate(:, :)
    integer, allocatable :: asserted(:)
    integer :: a
    integer :: b
    integer :: i
    integer :: j
    integer :: old_columns
    logical :: desc

    result%ok = .false.
    result%message = ''
    if (nlevels < 1) then
      result%message = 'nlevels must be positive'
      return
    end if
    if (size(significant) /= size(group_a) .or. size(significant) /= size(group_b)) then
      result%message = 'comparison arrays must have identical lengths'
      return
    end if
    if (any(group_a < 1) .or. any(group_a > nlevels) .or. &
        any(group_b < 1) .or. any(group_b > nlevels)) then
      result%message = 'comparison group index is outside 1:nlevels'
      return
    end if

    desc = .false.
    if (present(decreasing)) desc = decreasing
    allocate(matrix(nlevels, 1))
    matrix = .true.

    if (any(significant)) then
      do i = 1, size(significant)
        if (.not. significant(i)) cycle
        a = group_a(i)
        b = group_b(i)
        call columns_asserting_both(matrix, a, b, asserted)
        if (size(asserted) == 0) cycle
        old_columns = size(matrix, 2)
        allocate(duplicate(nlevels, size(asserted)))
        duplicate = matrix(:, asserted)
        duplicate(b, :) = .false.
        do j = 1, size(asserted)
          matrix(a, asserted(j)) = .false.
        end do
        call append_logical_columns(matrix, duplicate)
        deallocate(duplicate)
        call absorb_columns(matrix)
        if (size(matrix, 2) > old_columns + size(asserted)) then
          result%message = 'internal compact-letter column accounting failed'
          return
        end if
      end do

      call sort_columns_by_count(matrix)
      call sweep_letters(matrix)
      call sort_columns_by_first(matrix, desc)
      call sort_columns_by_count(matrix)
      call sweep_letters(matrix)
      call sort_columns_by_first(matrix, desc)
    end if

    call populate_result(matrix, result)
    result%ok = .true.
  end subroutine compact_letter_display

  subroutine columns_asserting_both(matrix, a, b, indices)
    logical, intent(in) :: matrix(:, :) !! Current letter membership matrix.
    integer, intent(in) :: a !! First group index whose shared letters are sought.
    integer, intent(in) :: b !! Second group index whose shared letters are sought.
    integer, allocatable, intent(out) :: indices(:) !! Columns currently containing both groups.

    integer :: count_both
    integer :: i
    integer :: j

    count_both = 0
    do i = 1, size(matrix, 2)
      if (matrix(a, i) .and. matrix(b, i)) count_both = count_both + 1
    end do
    allocate(indices(count_both))
    j = 0
    do i = 1, size(matrix, 2)
      if (matrix(a, i) .and. matrix(b, i)) then
        j = j + 1
        indices(j) = i
      end if
    end do
  end subroutine columns_asserting_both

  subroutine append_logical_columns(matrix, extra)
    logical, allocatable, intent(inout) :: matrix(:, :) !! Existing membership matrix, extended in place.
    logical, intent(in) :: extra(:, :) !! New membership columns to append to matrix.

    logical, allocatable :: combined(:, :)
    integer :: old_columns

    old_columns = size(matrix, 2)
    allocate(combined(size(matrix, 1), old_columns + size(extra, 2)))
    combined(:, 1:old_columns) = matrix
    combined(:, old_columns + 1:) = extra
    call move_alloc(combined, matrix)
  end subroutine append_logical_columns

  recursive subroutine absorb_columns(matrix)
    logical, allocatable, intent(inout) :: matrix(:, :) !! Membership matrix from which contained columns are recursively removed.

    integer :: j
    integer :: k

    if (size(matrix, 2) <= 1) return
    do j = 1, size(matrix, 2) - 1
      do k = j + 1, size(matrix, 2)
        if (column_subset(matrix(:, k), matrix(:, j))) then
          call remove_column(matrix, k)
          call absorb_columns(matrix)
          return
        else if (column_subset(matrix(:, j), matrix(:, k))) then
          call remove_column(matrix, j)
          call absorb_columns(matrix)
          return
        end if
      end do
    end do
  end subroutine absorb_columns

  logical function column_subset(candidate, container) result(answer)
    logical, intent(in) :: candidate(:) !! Column whose true entries must all be present in container.
    logical, intent(in) :: container(:) !! Column tested as a superset of candidate.

    integer :: i

    answer = .true.
    do i = 1, size(candidate)
      if (candidate(i) .and. .not. container(i)) then
        answer = .false.
        return
      end if
    end do
  end function column_subset

  subroutine remove_column(matrix, column)
    logical, allocatable, intent(inout) :: matrix(:, :) !! Membership matrix to shrink by one column.
    integer, intent(in) :: column !! One-based column index to remove.

    logical, allocatable :: reduced(:, :)
    integer :: j
    integer :: out_col

    if (size(matrix, 2) == 1) then
      allocate(reduced(size(matrix, 1), 0))
      call move_alloc(reduced, matrix)
      return
    end if
    allocate(reduced(size(matrix, 1), size(matrix, 2) - 1))
    out_col = 0
    do j = 1, size(matrix, 2)
      if (j == column) cycle
      out_col = out_col + 1
      reduced(:, out_col) = matrix(:, j)
    end do
    call move_alloc(reduced, matrix)
  end subroutine remove_column

  subroutine sort_columns_by_count(matrix)
    logical, allocatable, intent(inout) :: matrix(:, :) !! Membership matrix sorted by increasing number of true entries.

    integer :: i
    integer :: j
    integer :: key_count
    logical, allocatable :: key(:)

    if (size(matrix, 2) <= 1) return
    allocate(key(size(matrix, 1)))
    do i = 2, size(matrix, 2)
      key = matrix(:, i)
      key_count = count(key)
      j = i - 1
      do while (j >= 1)
        if (count(matrix(:, j)) <= key_count) exit
        matrix(:, j + 1) = matrix(:, j)
        j = j - 1
      end do
      matrix(:, j + 1) = key
    end do
  end subroutine sort_columns_by_count

  subroutine sort_columns_by_first(matrix, decreasing)
    logical, allocatable, intent(inout) :: matrix(:, :) !! Membership matrix sorted by first true row position.
    logical, intent(in) :: decreasing !! If true, sort first-row positions from largest to smallest.

    integer :: i
    integer :: j
    integer :: key_first
    logical, allocatable :: key(:)

    if (size(matrix, 2) <= 1) return
    allocate(key(size(matrix, 1)))
    do i = 2, size(matrix, 2)
      key = matrix(:, i)
      key_first = first_true(key)
      j = i - 1
      do while (j >= 1)
        if (decreasing) then
          if (first_true(matrix(:, j)) >= key_first) exit
        else
          if (first_true(matrix(:, j)) <= key_first) exit
        end if
        matrix(:, j + 1) = matrix(:, j)
        j = j - 1
      end do
      matrix(:, j + 1) = key
    end do
  end subroutine sort_columns_by_first

  integer function first_true(column) result(index)
    logical, intent(in) :: column(:) !! Membership column whose first true row is required.

    integer :: i

    index = size(column) + 1
    do i = 1, size(column)
      if (column(i)) then
        index = i
        return
      end if
    end do
  end function first_true

  recursive subroutine sweep_letters(matrix)
    logical, allocatable, intent(inout) :: matrix(:, :) !! Letter matrix from which redundant true entries are swept.

    logical, allocatable :: locked(:, :)
    logical, allocatable :: tmp(:, :)
    integer, allocatable :: one(:)
    integer :: check_count
    integer :: hit
    integer :: i
    integer :: j
    integer :: k
    integer :: n_one
    logical :: redundant

    if (size(matrix, 2) <= 1) return
    allocate(locked(size(matrix, 1), size(matrix, 2)))
    locked = .false.

    do i = 1, size(matrix, 2)
      allocate(tmp(size(matrix, 1), size(matrix, 2)))
      tmp = .false.
      do j = 1, size(matrix, 1)
        if (matrix(j, i)) tmp(j, :) = matrix(j, :)
      end do
      call true_indices(matrix(:, i), one)
      n_one = size(one)
      if (n_one == 0) then
        deallocate(tmp, one)
        call remove_column(matrix, i)
        call sweep_letters(matrix)
        return
      end if

      if (all_rows_have_other(tmp, i)) then
        deallocate(tmp, one)
        cycle
      end if

      do j = 1, n_one
        if (locked(one(j), i)) cycle
        check_count = 0
        redundant = .true.
        do k = 1, n_one
          if (j == k) cycle
          hit = last_shared_other_column(tmp, one(j), one(k), i)
          if (hit > 0) then
            check_count = check_count + 1
          else
            redundant = .false.
            exit
          end if
        end do
        if (redundant .and. check_count == n_one - 1 .and. check_count /= 0) then
          do k = 1, n_one
            if (j == k) cycle
            hit = last_shared_other_column(tmp, one(j), one(k), i)
            if (hit > 0) then
              locked(one(j), hit) = .true.
              locked(one(k), hit) = .true.
            end if
          end do
          matrix(one(j), i) = .false.
        end if
      end do
      deallocate(tmp, one)
      if (.not. any(matrix(:, i))) then
        call remove_column(matrix, i)
        call sweep_letters(matrix)
        return
      end if
    end do
  end subroutine sweep_letters

  subroutine true_indices(column, indices)
    logical, intent(in) :: column(:) !! Logical vector whose true element positions are requested.
    integer, allocatable, intent(out) :: indices(:) !! One-based positions of true elements.

    integer :: i
    integer :: j

    allocate(indices(count(column)))
    j = 0
    do i = 1, size(column)
      if (column(i)) then
        j = j + 1
        indices(j) = i
      end if
    end do
  end subroutine true_indices

  logical function all_rows_have_other(matrix, excluded_column) result(answer)
    logical, intent(in) :: matrix(:, :) !! Temporary matrix containing active rows from one letter column.
    integer, intent(in) :: excluded_column !! Column excluded when checking alternative letter coverage.

    integer :: i
    integer :: j
    logical :: any_other

    answer = .true.
    do i = 1, size(matrix, 1)
      any_other = .false.
      do j = 1, size(matrix, 2)
        if (j == excluded_column) cycle
        if (matrix(i, j)) then
          any_other = .true.
          exit
        end if
      end do
      if (.not. any_other) then
        answer = .false.
        return
      end if
    end do
  end function all_rows_have_other

  integer function last_shared_other_column(matrix, row_a, row_b, excluded_column) result(column)
    logical, intent(in) :: matrix(:, :) !! Temporary membership matrix for pairwise shared-letter checks.
    integer, intent(in) :: row_a !! First row in the pair whose common alternative letter is sought.
    integer, intent(in) :: row_b !! Second row in the pair whose common alternative letter is sought.
    integer, intent(in) :: excluded_column !! Current letter column that cannot satisfy the pair.

    integer :: j

    column = 0
    do j = 1, size(matrix, 2)
      if (j == excluded_column) cycle
      if (matrix(row_a, j) .and. matrix(row_b, j)) column = j
    end do
  end function last_shared_other_column

  subroutine populate_result(matrix, result)
    logical, intent(in) :: matrix(:, :) !! Final compact-letter membership matrix.
    type(cld_type), intent(inout) :: result !! Result object receiving labels and per-level strings.

    character(len=32) :: label
    integer :: i
    integer :: j
    integer :: pos

    allocate(result%letter_matrix(size(matrix, 1), size(matrix, 2)))
    allocate(result%letters(size(matrix, 1)))
    allocate(result%monospaced_letters(size(matrix, 1)))
    allocate(result%column_labels(size(matrix, 2)))
    result%letter_matrix = matrix
    result%letters = ''
    result%monospaced_letters = ''

    do j = 1, size(matrix, 2)
      call letter_label(j, label)
      result%column_labels(j) = trim(label)
    end do

    do i = 1, size(matrix, 1)
      pos = 1
      do j = 1, size(matrix, 2)
        label = trim(result%column_labels(j))
        if (matrix(i, j)) then
          result%letters(i) = trim(result%letters(i)) // trim(label)
          call place_text(result%monospaced_letters(i), pos, trim(label))
        end if
        pos = pos + len_trim(label)
      end do
    end do
  end subroutine populate_result

  subroutine letter_label(index, label)
    integer, intent(in) :: index !! One-based letter-column number to encode as a, ..., Z, .a, and so on.
    character(len=*), intent(out) :: label !! Recycled letter label with dot prefixes after each set of 52.

    character(len=52), parameter :: alphabet = &
      'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: base_index
    integer :: prefix_count

    prefix_count = (index - 1) / 52
    base_index = modulo(index - 1, 52) + 1
    label = repeat('.', prefix_count) // alphabet(base_index:base_index)
  end subroutine letter_label

  subroutine place_text(destination, start, text)
    character(len=*), intent(inout) :: destination !! Fixed-width string receiving a monospaced letter at a known offset.
    integer, intent(in) :: start !! One-based starting character position in destination.
    character(len=*), intent(in) :: text !! Letter label to place without collapsing absent columns.

    integer :: last

    last = min(len(destination), start + len_trim(text) - 1)
    if (start <= len(destination) .and. last >= start) then
      destination(start:last) = text(1:last - start + 1)
    end if
  end subroutine place_text

end module multcomp_cld
