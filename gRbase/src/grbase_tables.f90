module grbase_tables
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value
  use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
  use r_kinds, only : dp
  use grbase_arrays, only : cell_to_entry, entry_to_cell
  use grbase_types, only : table_t
  implicit none
  private

  public :: make_table
  public :: valid_table
  public :: table_permute
  public :: table_expand
  public :: table_align
  public :: table_margin
  public :: table_add
  public :: table_subtract
  public :: table_multiply
  public :: table_divide
  public :: table_divide_zero
  public :: table_equal
  public :: table_normalize_all
  public :: table_normalize_first
  public :: table_slice
  public :: table_sample_from_uniforms
  public :: table_list_add
  public :: table_list_multiply

contains

  pure function make_table(var, dim, value) result(tab)
    integer, intent(in) :: var(:) !! Integer variable labels; labels must be unique and correspond positionally to `dim`.
    integer, intent(in) :: dim(:) !! Positive cardinalities for the variables in `var`.
    real(dp), intent(in) :: value(:) !! Flattened table entries in R/Fortran column-major order.
    type(table_t) :: tab

    if (size(var) /= size(dim)) return
    if (any(dim <= 0)) return
    if (has_duplicates(var)) return
    if (size(value) /= product_or_one(dim)) return
    allocate(tab%var(size(var)))
    allocate(tab%dim(size(dim)))
    allocate(tab%value(size(value)))
    tab%var = var
    tab%dim = dim
    tab%value = value
  end function make_table

  pure logical function valid_table(tab) result(ok)
    type(table_t), intent(in) :: tab !! Table object whose metadata and value-array consistency are checked.

    ok = allocated(tab%var) .and. allocated(tab%dim) .and. allocated(tab%value)
    if (.not. ok) return
    ok = size(tab%var) == size(tab%dim)
    if (.not. ok) return
    ok = all(tab%dim > 0)
    if (.not. ok) return
    ok = .not. has_duplicates(tab%var)
    if (.not. ok) return
    ok = size(tab%value) == product_or_one(tab%dim)
  end function valid_table

  pure function table_permute(tab, perm) result(out)
    type(table_t), intent(in) :: tab !! Source table to reorder.
    integer, intent(in) :: perm(:) !! Permutation of 1..rank(tab) giving the new variable order.
    type(table_t) :: out
    integer, allocatable :: old_cell(:)
    integer, allocatable :: new_cell(:)
    integer :: e
    integer :: new_entry

    if (.not. valid_table(tab)) return
    if (.not. valid_permutation(perm, size(tab%var))) return
    allocate(out%var(size(tab%var)))
    allocate(out%dim(size(tab%dim)))
    allocate(out%value(size(tab%value)))
    out%var = tab%var(perm)
    out%dim = tab%dim(perm)
    allocate(old_cell(size(tab%dim)), new_cell(size(tab%dim)))
    do e = 1, size(tab%value)
      old_cell = entry_to_cell(e, tab%dim)
      new_cell = old_cell(perm)
      new_entry = cell_to_entry(new_cell, out%dim)
      out%value(new_entry) = tab%value(e)
    end do
  end function table_permute

  pure function table_expand(tab, aux_var, aux_dim, expand_type) result(out)
    type(table_t), intent(in) :: tab !! Source table whose domain is enlarged or reordered.
    integer, intent(in) :: aux_var(:) !! Auxiliary variable labels, whose order is retained at the slow-varying end of the result.
    integer, intent(in) :: aux_dim(:) !! Cardinalities for `aux_var`; shared variables must match the source.
    integer, intent(in), optional :: expand_type !! Rule: 0 replicate, 1 divide replicated values, 2 zero-fill new cells.
    type(table_t) :: out
    integer, allocatable :: only_tab(:)
    integer, allocatable :: out_cell(:)
    integer, allocatable :: tab_cell(:)
    integer :: e
    integer :: i
    integer :: idx
    integer :: mode
    integer :: new_factor
    integer :: src_entry
    logical :: added_at_first

    if (.not. valid_table(tab)) return
    if (size(aux_var) /= size(aux_dim)) return
    if (any(aux_dim <= 0)) return
    if (has_duplicates(aux_var)) return
    do i = 1, size(aux_var)
      idx = find_position(tab%var, aux_var(i))
      if (idx > 0) then
        if (tab%dim(idx) /= aux_dim(i)) return
      end if
    end do
    mode = 0
    if (present(expand_type)) mode = expand_type
    if (mode < 0 .or. mode > 2) return

    only_tab = set_difference_preserve(tab%var, aux_var)
    allocate(out%var(size(only_tab) + size(aux_var)))
    if (size(only_tab) > 0) out%var(:size(only_tab)) = only_tab
    if (size(aux_var) > 0) out%var(size(only_tab) + 1:) = aux_var
    allocate(out%dim(size(out%var)))
    do i = 1, size(out%var)
      idx = find_position(tab%var, out%var(i))
      if (idx > 0) then
        out%dim(i) = tab%dim(idx)
      else
        idx = find_position(aux_var, out%var(i))
        out%dim(i) = aux_dim(idx)
      end if
    end do

    new_factor = 1
    do i = 1, size(aux_var)
      if (find_position(tab%var, aux_var(i)) == 0) new_factor = new_factor * aux_dim(i)
    end do
    allocate(out%value(product_or_one(out%dim)))
    allocate(out_cell(size(out%dim)), tab_cell(size(tab%dim)))
    do e = 1, size(out%value)
      out_cell = entry_to_cell(e, out%dim)
      added_at_first = .true.
      do i = 1, size(tab%var)
        idx = find_position(out%var, tab%var(i))
        tab_cell(i) = out_cell(idx)
      end do
      if (mode == 2) then
        do i = 1, size(aux_var)
          if (find_position(tab%var, aux_var(i)) == 0) then
            idx = find_position(out%var, aux_var(i))
            if (out_cell(idx) /= 1) added_at_first = .false.
          end if
        end do
      end if
      src_entry = cell_to_entry(tab_cell, tab%dim)
      if (mode == 2 .and. .not. added_at_first) then
        out%value(e) = 0.0_dp
      else if (mode == 1) then
        out%value(e) = tab%value(src_entry) / real(new_factor, dp)
      else
        out%value(e) = tab%value(src_entry)
      end if
    end do
  end function table_expand

  pure function table_align(tab, reference) result(out)
    type(table_t), intent(in) :: tab !! Source table to align.
    type(table_t), intent(in) :: reference !! Reference table whose variable order determines the result.
    type(table_t) :: out
    integer, allocatable :: perm(:)
    integer :: i
    integer :: idx

    if (.not. valid_table(tab) .or. .not. valid_table(reference)) return
    if (size(tab%var) /= size(reference%var)) return
    allocate(perm(size(tab%var)))
    do i = 1, size(reference%var)
      idx = find_position(tab%var, reference%var(i))
      if (idx == 0) return
      if (tab%dim(idx) /= reference%dim(i)) return
      perm(i) = idx
    end do
    out = table_permute(tab, perm)
  end function table_align

  pure function table_margin(tab, margin_var) result(out)
    type(table_t), intent(in) :: tab !! Source table to sum over variables outside the requested margin.
    integer, intent(in) :: margin_var(:) !! Variable labels retained in the result, in the requested output order.
    type(table_t) :: out
    integer, allocatable :: in_cell(:)
    integer, allocatable :: out_cell(:)
    integer :: e
    integer :: i
    integer :: idx
    integer :: oe

    if (.not. valid_table(tab)) return
    if (has_duplicates(margin_var)) return
    do i = 1, size(margin_var)
      if (find_position(tab%var, margin_var(i)) == 0) return
    end do
    allocate(out%var(size(margin_var)))
    allocate(out%dim(size(margin_var)))
    out%var = margin_var
    do i = 1, size(margin_var)
      idx = find_position(tab%var, margin_var(i))
      out%dim(i) = tab%dim(idx)
    end do
    allocate(out%value(product_or_one(out%dim)))
    out%value = 0.0_dp
    allocate(in_cell(size(tab%dim)), out_cell(size(out%dim)))
    do e = 1, size(tab%value)
      in_cell = entry_to_cell(e, tab%dim)
      do i = 1, size(margin_var)
        idx = find_position(tab%var, margin_var(i))
        out_cell(i) = in_cell(idx)
      end do
      oe = cell_to_entry(out_cell, out%dim)
      out%value(oe) = out%value(oe) + tab%value(e)
    end do
  end function table_margin

  pure function table_add(tab1, tab2) result(out)
    type(table_t), intent(in) :: tab1 !! Left table operand.
    type(table_t), intent(in) :: tab2 !! Right table operand.
    type(table_t) :: out

    out = table_binary(tab1, tab2, 1)
  end function table_add

  pure function table_subtract(tab1, tab2) result(out)
    type(table_t), intent(in) :: tab1 !! Left table operand.
    type(table_t), intent(in) :: tab2 !! Right table operand.
    type(table_t) :: out

    out = table_binary(tab1, tab2, 2)
  end function table_subtract

  pure function table_multiply(tab1, tab2) result(out)
    type(table_t), intent(in) :: tab1 !! Left table operand.
    type(table_t), intent(in) :: tab2 !! Right table operand.
    type(table_t) :: out

    out = table_binary(tab1, tab2, 3)
  end function table_multiply

  pure function table_divide(tab1, tab2) result(out)
    type(table_t), intent(in) :: tab1 !! Numerator table.
    type(table_t), intent(in) :: tab2 !! Denominator table.
    type(table_t) :: out

    out = table_binary(tab1, tab2, 4)
  end function table_divide

  pure function table_divide_zero(tab1, tab2) result(out)
    type(table_t), intent(in) :: tab1 !! Numerator table.
    type(table_t), intent(in) :: tab2 !! Denominator table; zero or nonfinite quotients are mapped to zero.
    type(table_t) :: out
    integer :: i

    out = table_binary(tab1, tab2, 5)
    if (.not. valid_table(out)) return
    do i = 1, size(out%value)
      if (ieee_is_nan(out%value(i)) .or. .not. ieee_is_finite(out%value(i))) out%value(i) = 0.0_dp
    end do
  end function table_divide_zero

  pure logical function table_equal(tab1, tab2, eps) result(ok)
    type(table_t), intent(in) :: tab1 !! First table compared after domain alignment.
    type(table_t), intent(in) :: tab2 !! Second table compared after domain alignment.
    real(dp), intent(in), optional :: eps !! Absolute L1 tolerance; defaults to 1e-12.
    type(table_t) :: diff
    real(dp) :: tol

    tol = 1.0e-12_dp
    if (present(eps)) tol = eps
    ok = .false.
    if (.not. same_domain(tab1, tab2)) return
    diff = table_subtract(tab1, tab2)
    if (.not. valid_table(diff)) return
    ok = sum(abs(diff%value)) < tol
  end function table_equal

  pure function table_normalize_all(tab) result(out)
    type(table_t), intent(in) :: tab !! Table whose entire value array is normalized to unit sum.
    type(table_t) :: out
    real(dp) :: total

    if (.not. valid_table(tab)) return
    out = tab
    total = sum(out%value)
    out%value = r_divide(out%value, total)
  end function table_normalize_all

  pure function table_normalize_first(tab) result(out)
    type(table_t), intent(in) :: tab !! Table normalized over its first, fastest-varying dimension for each remaining configuration.
    type(table_t) :: out
    integer :: block
    integer :: first_dim
    integer :: offset
    real(dp) :: total

    if (.not. valid_table(tab)) return
    out = tab
    if (size(tab%dim) == 0) then
      out%value(1) = r_divide(out%value(1), out%value(1))
      return
    end if
    first_dim = tab%dim(1)
    do block = 0, size(tab%value) / first_dim - 1
      offset = block * first_dim
      total = sum(tab%value(offset + 1:offset + first_dim))
      out%value(offset + 1:offset + first_dim) = &
        r_divide(tab%value(offset + 1:offset + first_dim), total)
    end do
  end function table_normalize_first

  pure function table_slice(tab, fixed_var, fixed_level) result(out)
    type(table_t), intent(in) :: tab !! Source table from which a level slice is extracted.
    integer, intent(in) :: fixed_var(:) !! Variable labels fixed to specified levels and removed from the result domain.
    integer, intent(in) :: fixed_level(:) !! One-based levels corresponding to `fixed_var`.
    type(table_t) :: out
    integer, allocatable :: remain(:)
    integer, allocatable :: in_cell(:)
    integer, allocatable :: out_cell(:)
    integer :: e
    integer :: i
    integer :: idx
    integer :: oe
    logical :: match

    if (.not. valid_table(tab)) return
    if (size(fixed_var) /= size(fixed_level)) return
    if (has_duplicates(fixed_var)) return
    do i = 1, size(fixed_var)
      idx = find_position(tab%var, fixed_var(i))
      if (idx == 0) return
      if (fixed_level(i) < 1 .or. fixed_level(i) > tab%dim(idx)) return
    end do
    remain = set_difference_preserve(tab%var, fixed_var)
    out%var = remain
    allocate(out%dim(size(remain)))
    do i = 1, size(remain)
      idx = find_position(tab%var, remain(i))
      out%dim(i) = tab%dim(idx)
    end do
    allocate(out%value(product_or_one(out%dim)))
    out%value = 0.0_dp
    allocate(in_cell(size(tab%dim)), out_cell(size(out%dim)))
    do e = 1, size(tab%value)
      in_cell = entry_to_cell(e, tab%dim)
      match = .true.
      do i = 1, size(fixed_var)
        idx = find_position(tab%var, fixed_var(i))
        if (in_cell(idx) /= fixed_level(i)) match = .false.
      end do
      if (.not. match) cycle
      do i = 1, size(remain)
        idx = find_position(tab%var, remain(i))
        out_cell(i) = in_cell(idx)
      end do
      oe = cell_to_entry(out_cell, out%dim)
      out%value(oe) = tab%value(e)
    end do
  end function table_slice

  pure function table_list_add(tables) result(out)
    type(table_t), intent(in) :: tables(:) !! Nonempty sequence of tables combined by gRbase domain-aware addition.
    type(table_t) :: out
    integer :: i

    if (size(tables) == 0) return
    out = tables(1)
    do i = 2, size(tables)
      out = table_add(out, tables(i))
    end do
  end function table_list_add

  pure function table_list_multiply(tables) result(out)
    type(table_t), intent(in) :: tables(:) !! Nonempty sequence of tables combined by gRbase domain-aware multiplication.
    type(table_t) :: out
    integer :: i

    if (size(tables) == 0) return
    out = tables(1)
    do i = 2, size(tables)
      out = table_multiply(out, tables(i))
    end do
  end function table_list_multiply

  pure function table_sample_from_uniforms(tab, uniforms) result(cells)
    type(table_t), intent(in) :: tab !! Table whose nonnegative entries define unnormalized sampling probabilities.
    real(dp), intent(in) :: uniforms(:) !! Deterministic uniforms in [0,1]; one sample is generated per value.
    integer, allocatable :: cells(:, :)
    real(dp), allocatable :: cumulative(:)
    real(dp) :: threshold
    real(dp) :: total
    integer :: entry
    integer :: i

    allocate(cells(size(uniforms), size(tab%dim)))
    cells = 0
    if (.not. valid_table(tab)) return
    if (size(tab%value) == 0) return
    if (any(tab%value < 0.0_dp)) return

    total = sum(tab%value)
    if (.not. ieee_is_finite(total)) return
    if (total <= 0.0_dp) return

    allocate(cumulative(size(tab%value)))
    cumulative(1) = tab%value(1)
    do i = 2, size(tab%value)
      cumulative(i) = cumulative(i - 1) + tab%value(i)
    end do

    do i = 1, size(uniforms)
      threshold = max(0.0_dp, min(1.0_dp, uniforms(i))) * total
      entry = 1
      do while (entry < size(cumulative) .and. cumulative(entry) <= threshold)
        entry = entry + 1
      end do
      cells(i, :) = entry_to_cell(entry, tab%dim)
    end do
  end function table_sample_from_uniforms

  pure function table_binary(tab1, tab2, operation) result(out)
    type(table_t), intent(in) :: tab1 !! Left operand table.
    type(table_t), intent(in) :: tab2 !! Right operand table.
    integer, value :: operation !! Operation code: 1 add, 2 subtract, 3 multiply, 4 divide, 5 divide-with-zero replacement.
    type(table_t) :: out
    integer, allocatable :: only1(:)
    integer, allocatable :: cell(:)
    integer, allocatable :: cell1(:)
    integer, allocatable :: cell2(:)
    integer :: e
    integer :: i
    integer :: idx
    integer :: e1
    integer :: e2
    real(dp) :: denom

    if (.not. same_shared_cardinalities(tab1, tab2)) return
    only1 = set_difference_preserve(tab1%var, tab2%var)
    allocate(out%var(size(only1) + size(tab2%var)))
    if (size(only1) > 0) out%var(:size(only1)) = only1
    if (size(tab2%var) > 0) out%var(size(only1) + 1:) = tab2%var
    allocate(out%dim(size(out%var)))
    do i = 1, size(out%var)
      idx = find_position(tab1%var, out%var(i))
      if (idx > 0) then
        out%dim(i) = tab1%dim(idx)
      else
        idx = find_position(tab2%var, out%var(i))
        out%dim(i) = tab2%dim(idx)
      end if
    end do
    allocate(out%value(product_or_one(out%dim)))
    allocate(cell(size(out%dim)), cell1(size(tab1%dim)), cell2(size(tab2%dim)))
    do e = 1, size(out%value)
      cell = entry_to_cell(e, out%dim)
      do i = 1, size(tab1%var)
        idx = find_position(out%var, tab1%var(i))
        cell1(i) = cell(idx)
      end do
      do i = 1, size(tab2%var)
        idx = find_position(out%var, tab2%var(i))
        cell2(i) = cell(idx)
      end do
      e1 = cell_to_entry(cell1, tab1%dim)
      e2 = cell_to_entry(cell2, tab2%dim)
      select case (operation)
      case (1)
        out%value(e) = tab1%value(e1) + tab2%value(e2)
      case (2)
        out%value(e) = tab1%value(e1) - tab2%value(e2)
      case (3)
        out%value(e) = tab1%value(e1) * tab2%value(e2)
      case (4)
        out%value(e) = r_divide(tab1%value(e1), tab2%value(e2))
      case (5)
        denom = tab2%value(e2)
        out%value(e) = r_divide(tab1%value(e1), denom)
        if (ieee_is_nan(out%value(e)) .or. .not. ieee_is_finite(out%value(e))) out%value(e) = 0.0_dp
      case default
        deallocate(out%var, out%dim, out%value)
        return
      end select
    end do
  end function table_binary

  pure elemental real(dp) function r_divide(numerator, denominator) result(value)
    real(dp), intent(in) :: numerator !! Numerator whose IEEE quotient is requested.
    real(dp), intent(in) :: denominator !! Denominator; zero produces IEEE infinities or NaN as in R arithmetic.

    if (ieee_is_nan(numerator) .or. ieee_is_nan(denominator)) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (denominator > 0.0_dp .or. denominator < 0.0_dp) then
      value = numerator / denominator
    else if (numerator > 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_positive_inf)
    else if (numerator < 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_negative_inf)
    else
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    end if
  end function r_divide

  pure logical function same_domain(tab1, tab2) result(ok)
    type(table_t), intent(in) :: tab1 !! First table domain.
    type(table_t), intent(in) :: tab2 !! Second table domain.
    integer :: i
    integer :: idx

    ok = valid_table(tab1) .and. valid_table(tab2)
    if (.not. ok) return
    if (size(tab1%var) /= size(tab2%var)) then
      ok = .false.
      return
    end if
    do i = 1, size(tab1%var)
      idx = find_position(tab2%var, tab1%var(i))
      if (idx == 0 .or. tab1%dim(i) /= tab2%dim(idx)) then
        ok = .false.
        return
      end if
    end do
  end function same_domain

  pure logical function same_shared_cardinalities(tab1, tab2) result(ok)
    type(table_t), intent(in) :: tab1 !! First table participating in a binary operation.
    type(table_t), intent(in) :: tab2 !! Second table participating in a binary operation.
    integer :: i
    integer :: idx

    ok = valid_table(tab1) .and. valid_table(tab2)
    if (.not. ok) return
    do i = 1, size(tab1%var)
      idx = find_position(tab2%var, tab1%var(i))
      if (idx > 0) then
        if (tab1%dim(i) /= tab2%dim(idx)) then
          ok = .false.
          return
        end if
      end if
    end do
  end function same_shared_cardinalities

  pure function set_difference_preserve(a, b) result(out)
    integer, intent(in) :: a(:) !! Ordered source sequence whose labels are unique.
    integer, intent(in) :: b(:) !! Labels removed from the source sequence.
    integer, allocatable :: out(:)
    integer, allocatable :: tmp(:)
    integer :: i
    integer :: n

    allocate(tmp(size(a)))
    n = 0
    do i = 1, size(a)
      if (find_position(b, a(i)) == 0) then
        n = n + 1
        tmp(n) = a(i)
      end if
    end do
    allocate(out(n))
    if (n > 0) out = tmp(:n)
  end function set_difference_preserve

  pure integer function find_position(x, value) result(pos)
    integer, intent(in) :: x(:) !! Integer sequence searched from its first element.
    integer, value :: value !! Label whose one-based position is requested.
    integer :: i

    pos = 0
    do i = 1, size(x)
      if (x(i) == value) then
        pos = i
        return
      end if
    end do
  end function find_position

  pure logical function has_duplicates(x) result(found)
    integer, intent(in) :: x(:) !! Integer labels checked for duplicate values.
    integer :: i
    integer :: j

    found = .false.
    do i = 1, size(x) - 1
      do j = i + 1, size(x)
        if (x(i) == x(j)) then
          found = .true.
          return
        end if
      end do
    end do
  end function has_duplicates

  pure logical function valid_permutation(perm, n) result(ok)
    integer, intent(in) :: perm(:) !! Candidate permutation entries.
    integer, value :: n !! Required permutation length and maximum entry.
    integer :: i

    ok = size(perm) == n
    if (.not. ok) return
    do i = 1, n
      if (count(perm == i) /= 1) then
        ok = .false.
        return
      end if
    end do
  end function valid_permutation

  pure integer function product_or_one(dim) result(value)
    integer, intent(in) :: dim(:) !! Dimension sizes whose product is requested, with the empty product defined as one.

    if (size(dim) == 0) then
      value = 1
    else
      value = product(dim)
    end if
  end function product_or_one

end module grbase_tables
