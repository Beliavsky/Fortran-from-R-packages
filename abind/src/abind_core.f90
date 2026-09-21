module abind_core
   use abind_types, only : dp, array_value, int_vector, init_array, set_dimnames, has_dimnames, &
      array_rank, linear_index, unravel_index
   implicit none
   private

   public :: abind_arrays
   public :: asub_array
   public :: afill_array
   public :: adrop_array
   public :: acorn_array
   public :: indices_from_names
   public :: indices_from_mask

contains

   subroutine abind_arrays(inputs, output, status, along, rev_along, arg_names, use_first_dimnames, hier_names, use_dnns)
      type(array_value), intent(in) :: inputs(:) !! Arrays to combine; ranks must differ by at most one after join-rank resolution.
      type(array_value), intent(out) :: output !! Combined array with the join dimension restored to its requested position.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid ranks, dimensions, or join specification.
      real(dp), intent(in), optional :: along !! R-style join location; fractional/out-of-range values insert a new dimension.
      real(dp), intent(in), optional :: rev_along !! Reverse join location measured from the last existing dimension.
      character(len=*), intent(in), optional :: arg_names(:) !! Optional names of the input arguments used for join labels.
      logical, intent(in), optional :: use_first_dimnames !! If true, prefer the first labeled input on non-join dimensions.
      character(len=*), intent(in), optional :: hier_names !! Join-label mode: 'before', 'after', or 'none'.
      logical, intent(in), optional :: use_dnns !! If true, propagate names of dimensions when available.
      integer, allocatable :: expanded_dims(:, :)
      integer, allocatable :: source_coords(:)
      integer, allocatable :: expanded_coords(:)
      integer, allocatable :: out_coords(:)
      integer, allocatable :: out_dims(:)
      character(len=:), allocatable :: local_arg_names(:)
      character(len=:), allocatable :: join_labels(:)
      character(len=:), allocatable :: names_buffer(:)
      character(len=:), allocatable :: mode
      integer :: n_base
      integer :: n_out
      integer :: join_dim
      integer :: n_inputs
      integer :: i
      integer :: j
      integer :: d
      integer :: r
      integer :: source_dim
      integer :: offset
      integer :: src_idx
      integer :: dst_idx
      integer :: name_width
      integer :: pos
      integer :: preferred_start
      integer :: preferred_end
      integer :: preferred_step
      real(dp) :: requested_along
      logical :: take_first
      logical :: preserve_dnn
      logical :: have_join_names

      status = 0
      n_inputs = size(inputs)
      if (n_inputs == 0) then
         status = 1
         return
      end if

      n_base = 1
      do i = 1, n_inputs
         n_base = max(n_base, array_rank(inputs(i)))
      end do
      requested_along = real(n_base, dp)
      if (present(along)) requested_along = along
      if (present(rev_along)) requested_along = real(n_base + 1, dp) - rev_along

      if (requested_along < 1.0_dp .or. requested_along > real(n_base, dp) .or. &
         abs(requested_along - real(nint(requested_along), dp)) > 16.0_dp * epsilon(1.0_dp)) then
         n_out = n_base + 1
         join_dim = ceiling(requested_along)
         join_dim = max(1, min(n_out, join_dim))
      else
         n_out = n_base
         join_dim = nint(requested_along)
      end if
      if (join_dim < 1 .or. join_dim > n_out) then
         status = 2
         return
      end if

      allocate(expanded_dims(n_out, n_inputs))
      do i = 1, n_inputs
         r = array_rank(inputs(i))
         if (r == n_out) then
            expanded_dims(:, i) = inputs(i)%dims
         else if (r == n_out - 1) then
            if (join_dim > 1) expanded_dims(1:join_dim - 1, i) = inputs(i)%dims(1:join_dim - 1)
            expanded_dims(join_dim, i) = 1
            if (join_dim <= r) expanded_dims(join_dim + 1:n_out, i) = inputs(i)%dims(join_dim:r)
         else
            status = 3
            return
         end if
      end do

      do i = 2, n_inputs
         do d = 1, n_out
            if (d == join_dim) cycle
            if (expanded_dims(d, i) /= expanded_dims(d, 1)) then
               status = 4
               return
            end if
         end do
      end do

      allocate(out_dims(n_out))
      out_dims = expanded_dims(:, 1)
      out_dims(join_dim) = sum(expanded_dims(join_dim, :))
      call init_array(output, [(0.0_dp, i = 1, product_or_one(out_dims))], out_dims, status)
      if (status /= 0) return

      offset = 0
      do i = 1, n_inputs
         r = array_rank(inputs(i))
         allocate(source_coords(r), expanded_coords(n_out), out_coords(n_out))
         do src_idx = 1, size(inputs(i)%data)
            call unravel_index(src_idx, inputs(i)%dims, source_coords)
            if (r == n_out) then
               expanded_coords = source_coords
            else
               if (join_dim > 1) expanded_coords(1:join_dim - 1) = source_coords(1:join_dim - 1)
               expanded_coords(join_dim) = 1
               if (join_dim <= r) expanded_coords(join_dim + 1:n_out) = source_coords(join_dim:r)
            end if
            out_coords = expanded_coords
            out_coords(join_dim) = out_coords(join_dim) + offset
            dst_idx = linear_index(out_coords, out_dims)
            output%data(dst_idx) = inputs(i)%data(src_idx)
         end do
         offset = offset + expanded_dims(join_dim, i)
         deallocate(source_coords, expanded_coords, out_coords)
      end do

      take_first = .false.
      if (present(use_first_dimnames)) take_first = use_first_dimnames
      preserve_dnn = .false.
      if (present(use_dnns)) preserve_dnn = use_dnns
      mode = 'none'
      if (present(hier_names)) mode = trim(adjustl(hier_names))

      name_width = 1
      if (present(arg_names)) name_width = max(name_width, len(arg_names))
      do i = 1, n_inputs
         do d = 1, array_rank(inputs(i))
            if (has_dimnames(inputs(i), d)) name_width = max(name_width, len(inputs(i)%dimnames(d)%values))
         end do
      end do
      name_width = max(name_width + 24, 32)
      allocate(character(len=name_width) :: local_arg_names(n_inputs))
      local_arg_names = ''
      if (present(arg_names)) then
         if (size(arg_names) /= n_inputs) then
            status = 5
            return
         end if
         local_arg_names = arg_names
      end if

      if (take_first) then
         preferred_start = 1
         preferred_end = n_inputs
         preferred_step = 1
      else
         preferred_start = n_inputs
         preferred_end = 1
         preferred_step = -1
      end if

      do d = 1, n_out
         if (d == join_dim) cycle
         i = preferred_start
         do
            source_dim = expanded_to_source_dim(array_rank(inputs(i)), n_out, join_dim, d)
            if (source_dim > 0) then
               if (has_dimnames(inputs(i), source_dim)) then
                  allocate(character(len=name_width) :: names_buffer(expanded_dims(d, i)))
                  names_buffer = inputs(i)%dimnames(source_dim)%values
                  call set_dimnames(output, d, names_buffer, status)
                  deallocate(names_buffer)
                  if (status /= 0) return
                  if (preserve_dnn) call copy_dimension_name(inputs(i), source_dim, output, d)
                  exit
               end if
            end if
            if (i == preferred_end) exit
            i = i + preferred_step
         end do
      end do

      have_join_names = any(len_trim(local_arg_names) > 0)
      do i = 1, n_inputs
         source_dim = expanded_to_source_dim(array_rank(inputs(i)), n_out, join_dim, join_dim)
         if (source_dim > 0) have_join_names = have_join_names .or. has_dimnames(inputs(i), source_dim)
      end do
      if (have_join_names .and. out_dims(join_dim) > 0) then
         allocate(character(len=name_width) :: join_labels(out_dims(join_dim)))
         join_labels = ''
         pos = 0
         do i = 1, n_inputs
            source_dim = expanded_to_source_dim(array_rank(inputs(i)), n_out, join_dim, join_dim)
            do j = 1, expanded_dims(join_dim, i)
               pos = pos + 1
               if (source_dim > 0 .and. has_dimnames(inputs(i), source_dim)) then
                  join_labels(pos) = trim(inputs(i)%dimnames(source_dim)%values(j))
                  if (len_trim(local_arg_names(i)) > 0) then
                     if (mode == 'before') join_labels(pos) = trim(local_arg_names(i)) // '.' // trim(join_labels(pos))
                     if (mode == 'after') join_labels(pos) = trim(join_labels(pos)) // '.' // trim(local_arg_names(i))
                  end if
               else if (len_trim(local_arg_names(i)) > 0) then
                  if (expanded_dims(join_dim, i) == 1) then
                     join_labels(pos) = trim(local_arg_names(i))
                  else
                     write(join_labels(pos), '(a,i0)') trim(local_arg_names(i)), j
                  end if
               end if
            end do
         end do
         call set_dimnames(output, join_dim, join_labels, status)
         if (status /= 0) return
      end if
      if (preserve_dnn) then
         i = preferred_start
         do
            source_dim = expanded_to_source_dim(array_rank(inputs(i)), n_out, join_dim, join_dim)
            if (source_dim > 0 .and. allocated(inputs(i)%dimname_names)) then
               if (source_dim <= size(inputs(i)%dimname_names)) then
                  if (len_trim(inputs(i)%dimname_names(source_dim)) > 0) then
                     call copy_dimension_name(inputs(i), source_dim, output, join_dim)
                     exit
                  end if
               end if
            end if
            if (i == preferred_end) exit
            i = i + preferred_step
         end do
      end if
   end subroutine abind_arrays

   subroutine asub_array(x, indices, dims, output, status, drop)
      type(array_value), intent(in) :: x !! Source array to subset.
      type(int_vector), intent(in) :: indices(:) !! One-based selectors for DIMS; empty means all positions.
      integer, intent(in) :: dims(:) !! One-based dimensions on which the corresponding selectors act.
      type(array_value), intent(out) :: output !! Selected array, optionally with singleton dimensions dropped.
      integer, intent(out) :: status !! Zero on success; nonzero for inconsistent selectors or out-of-range indices.
      logical, intent(in), optional :: drop !! If true (default), remove selected dimensions whose resulting length is one.
      type(int_vector), allocatable :: selected(:)
      integer, allocatable :: full_dims(:)
      integer, allocatable :: keep_dims(:)
      integer, allocatable :: full_coords(:)
      integer, allocatable :: src_coords(:)
      integer, allocatable :: out_dims(:)
      integer :: rank_x
      integer :: d
      integer :: i
      integer :: j
      integer :: src_idx
      integer :: n_keep
      integer :: out_idx
      logical :: do_drop

      status = 0
      rank_x = array_rank(x)
      if (size(indices) /= size(dims)) then
         status = 1
         return
      end if
      if (any(dims < 1) .or. any(dims > rank_x)) then
         status = 2
         return
      end if
      do i = 1, size(dims)
         if (count(dims == dims(i)) /= 1) then
            status = 3
            return
         end if
      end do

      allocate(selected(rank_x))
      do d = 1, rank_x
         allocate(selected(d)%values(x%dims(d)))
         selected(d)%values = [(i, i = 1, x%dims(d))]
      end do
      do i = 1, size(dims)
         d = dims(i)
         if (allocated(indices(i)%values)) then
            if (size(indices(i)%values) > 0) then
               if (any(indices(i)%values < 1) .or. any(indices(i)%values > x%dims(d))) then
                  status = 4
                  return
               end if
               selected(d)%values = indices(i)%values
            end if
         end if
      end do

      allocate(full_dims(rank_x))
      do d = 1, rank_x
         full_dims(d) = size(selected(d)%values)
      end do
      do_drop = .true.
      if (present(drop)) do_drop = drop
      if (do_drop) then
         n_keep = count(full_dims /= 1)
      else
         n_keep = rank_x
      end if
      allocate(keep_dims(n_keep))
      j = 0
      do d = 1, rank_x
         if (.not. do_drop .or. full_dims(d) /= 1) then
            j = j + 1
            keep_dims(j) = d
         end if
      end do
      allocate(out_dims(n_keep))
      do i = 1, n_keep
         out_dims(i) = full_dims(keep_dims(i))
      end do

      call init_array(output, [(0.0_dp, i = 1, product_or_one(full_dims))], out_dims, status)
      if (status /= 0) return
      allocate(full_coords(rank_x), src_coords(rank_x))
      do out_idx = 1, size(output%data)
         call unravel_index(out_idx, full_dims, full_coords)
         do d = 1, rank_x
            src_coords(d) = selected(d)%values(full_coords(d))
         end do
         src_idx = linear_index(src_coords, x%dims)
         output%data(out_idx) = x%data(src_idx)
      end do

      do i = 1, n_keep
         d = keep_dims(i)
         if (has_dimnames(x, d)) then
            call subset_names(x, d, selected(d)%values, output, i, status)
            if (status /= 0) return
         end if
         if (allocated(x%dimname_names) .and. allocated(output%dimname_names)) then
            if (d <= size(x%dimname_names)) call copy_dimension_name(x, d, output, i)
         end if
      end do
   end subroutine asub_array

   subroutine afill_array(x, selectors, value, status, excess_ok)
      type(array_value), intent(inout) :: x !! Destination array modified in place.
      type(int_vector), intent(in) :: selectors(:) !! Destination selectors; empty entries match VALUE dimension labels.
      type(array_value), intent(in) :: value !! Source subarray replicated over explicitly selected destination dimensions.
      integer, intent(out) :: status !! Zero on success; nonzero for incompatible ranks, names, or selectors.
      logical, intent(in), optional :: excess_ok !! If true, discard VALUE labels absent from the destination dimension.
      type(int_vector), allocatable :: resolved(:)
      type(int_vector), allocatable :: value_positions(:)
      integer, allocatable :: sel_coords(:)
      integer, allocatable :: x_coords(:)
      integer, allocatable :: value_coords(:)
      integer, allocatable :: selection_dims(:)
      integer, allocatable :: missing_dims(:)
      integer :: rank_x
      integer :: rank_value
      integer :: d
      integer :: j
      integer :: k
      integer :: n_missing
      integer :: flat
      integer :: x_idx
      integer :: value_idx
      logical :: allow_excess

      status = 0
      rank_x = array_rank(x)
      rank_value = array_rank(value)
      if (size(selectors) /= rank_x) then
         status = 1
         return
      end if
      if (rank_value > rank_x) then
         status = 2
         return
      end if
      allow_excess = .false.
      if (present(excess_ok)) allow_excess = excess_ok

      n_missing = 0
      do d = 1, rank_x
         if (.not. allocated(selectors(d)%values)) then
            n_missing = n_missing + 1
         else if (size(selectors(d)%values) == 0) then
            n_missing = n_missing + 1
         end if
      end do
      if (n_missing /= rank_value) then
         status = 3
         return
      end if

      allocate(resolved(rank_x), missing_dims(n_missing), value_positions(rank_value))
      j = 0
      do d = 1, rank_x
         if (.not. allocated(selectors(d)%values) .or. size(selectors(d)%values) == 0) then
            j = j + 1
            missing_dims(j) = d
            call match_value_dimension(x, d, value, j, allow_excess, resolved(d), value_positions(j), status)
            if (status /= 0) return
         else
            if (any(selectors(d)%values < 1) .or. any(selectors(d)%values > x%dims(d))) then
               status = 4
               return
            end if
            resolved(d)%values = selectors(d)%values
         end if
      end do

      allocate(selection_dims(rank_x))
      do d = 1, rank_x
         selection_dims(d) = size(resolved(d)%values)
      end do
      allocate(sel_coords(rank_x), x_coords(rank_x), value_coords(rank_value))
      do flat = 1, product_or_one(selection_dims)
         call unravel_index(flat, selection_dims, sel_coords)
         do d = 1, rank_x
            x_coords(d) = resolved(d)%values(sel_coords(d))
         end do
         do j = 1, rank_value
            d = missing_dims(j)
            k = sel_coords(d)
            value_coords(j) = value_positions(j)%values(k)
         end do
         x_idx = linear_index(x_coords, x%dims)
         value_idx = linear_index(value_coords, value%dims)
         x%data(x_idx) = value%data(value_idx)
      end do
   end subroutine afill_array

   subroutine adrop_array(x, drop_dims, output, status, named_vector, one_d_array)
      type(array_value), intent(in) :: x !! Source array whose selected singleton dimensions are removed.
      integer, intent(in) :: drop_dims(:) !! One-based dimensions to remove; every selected dimension must have length one.
      type(array_value), intent(out) :: output !! Array with the selected dimensions removed.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid or non-singleton drop dimensions.
      logical, intent(in), optional :: named_vector !! If true, retain element labels for a rank-one result.
      logical, intent(in), optional :: one_d_array !! If true, preserve one-dimensional array metadata when possible.
      integer, allocatable :: keep(:)
      integer, allocatable :: out_dims(:)
      integer :: d
      integer :: i
      integer :: j
      integer :: rank_x
      logical :: keep_vector_names
      logical :: preserve_one_d

      status = 0
      rank_x = array_rank(x)
      keep_vector_names = .true.
      if (present(named_vector)) keep_vector_names = named_vector
      preserve_one_d = .false.
      if (present(one_d_array)) preserve_one_d = one_d_array
      if (any(drop_dims < 1) .or. any(drop_dims > rank_x)) then
         status = 1
         return
      end if
      do i = 1, size(drop_dims)
         if (count(drop_dims == drop_dims(i)) /= 1) then
            status = 2
            return
         end if
         if (x%dims(drop_dims(i)) /= 1) then
            status = 3
            return
         end if
      end do
      allocate(keep(rank_x - size(drop_dims)))
      j = 0
      do d = 1, rank_x
         if (.not. any(drop_dims == d)) then
            j = j + 1
            keep(j) = d
         end if
      end do
      allocate(out_dims(size(keep)))
      do i = 1, size(keep)
         out_dims(i) = x%dims(keep(i))
      end do
      call init_array(output, x%data, out_dims, status)
      if (status /= 0) return
      do i = 1, size(keep)
         d = keep(i)
         if (size(keep) == 1 .and. .not. preserve_one_d) then
            if (keep_vector_names .and. has_dimnames(x, d)) then
               call set_dimnames(output, i, x%dimnames(d)%values, status)
               if (status /= 0) return
            end if
         else
            if (has_dimnames(x, d)) then
               call set_dimnames(output, i, x%dimnames(d)%values, status)
               if (status /= 0) return
            end if
            if (allocated(x%dimname_names) .and. allocated(output%dimname_names)) then
               if (d <= size(x%dimname_names)) call copy_dimension_name(x, d, output, i)
            end if
         end if
      end do
   end subroutine adrop_array

   subroutine acorn_array(x, counts, output, status, addrownums)
      type(array_value), intent(in) :: x !! Source array from which a leading/trailing corner is selected.
      integer, intent(in) :: counts(:) !! Signed number of slices per dimension; positive takes the front, negative the end.
      type(array_value), intent(out) :: output !! Corner array with rank preserved.
      integer, intent(out) :: status !! Zero on success; nonzero if no count is supplied for the first dimension.
      logical, intent(in), optional :: addrownums !! If true (default), synthesize '[i]' labels for previously unnamed dimensions.
      type(int_vector), allocatable :: indices(:)
      integer, allocatable :: dims(:)
      integer :: rank_x
      integer :: d
      integer :: i
      integer :: n_take
      integer :: start_at
      logical :: add_labels
      character(len=32), allocatable :: labels(:)

      status = 0
      add_labels = .true.
      if (present(addrownums)) add_labels = addrownums
      rank_x = array_rank(x)
      if (rank_x == 0) then
         output = x
         return
      end if
      if (size(counts) == 0) then
         status = 1
         return
      end if
      allocate(indices(rank_x), dims(rank_x))
      dims = [(d, d = 1, rank_x)]
      do d = 1, rank_x
         if (d <= size(counts)) then
            n_take = sign(1, counts(d)) * min(abs(counts(d)), x%dims(d))
         else
            n_take = min(1, x%dims(d))
         end if
         if (n_take >= 0) then
            allocate(indices(d)%values(n_take))
            if (n_take > 0) indices(d)%values = [(start_at, start_at = 1, n_take)]
         else
            allocate(indices(d)%values(-n_take))
            start_at = x%dims(d) + n_take + 1
            indices(d)%values = [(i, i = start_at, x%dims(d))]
         end if
      end do
      call asub_array(x, indices, dims, output, status, drop=.false.)
      if (status /= 0 .or. .not. add_labels) return
      do d = 1, rank_x
         if (.not. has_dimnames(output, d)) then
            allocate(labels(size(indices(d)%values)))
            do i = 1, size(indices(d)%values)
               write(labels(i), '(a,i0,a)') '[', indices(d)%values(i), ']'
            end do
            call set_dimnames(output, d, labels, status)
            deallocate(labels)
            if (status /= 0) return
         end if
      end do
   end subroutine acorn_array

   subroutine indices_from_names(x, dim_number, names, indices, status)
      type(array_value), intent(in) :: x !! Array whose labels are searched.
      integer, intent(in) :: dim_number !! One-based dimension containing the labels.
      character(len=*), intent(in) :: names(:) !! Labels to resolve in the supplied order.
      type(int_vector), intent(out) :: indices !! One-based positions corresponding to NAMES.
      integer, intent(out) :: status !! Zero on success; nonzero when labels are absent or the dimension is invalid.
      integer :: i
      integer :: j

      status = 0
      if (.not. has_dimnames(x, dim_number)) then
         status = 1
         return
      end if
      allocate(indices%values(size(names)))
      do i = 1, size(names)
         indices%values(i) = 0
         do j = 1, x%dims(dim_number)
            if (trim(names(i)) == trim(x%dimnames(dim_number)%values(j))) then
               indices%values(i) = j
               exit
            end if
         end do
         if (indices%values(i) == 0) then
            status = 2
            return
         end if
      end do
   end subroutine indices_from_names

   subroutine indices_from_mask(mask, indices)
      logical, intent(in) :: mask(:) !! Logical selector converted to one-based positions of true entries.
      type(int_vector), intent(out) :: indices !! One-based positions at which MASK is true.
      integer :: i
      integer :: j

      allocate(indices%values(count(mask)))
      j = 0
      do i = 1, size(mask)
         if (mask(i)) then
            j = j + 1
            indices%values(j) = i
         end if
      end do
   end subroutine indices_from_mask

   subroutine subset_names(x, source_dim, positions, output, output_dim, status)
      type(array_value), intent(in) :: x !! Source array providing dimension labels.
      integer, intent(in) :: source_dim !! One-based source dimension containing the labels.
      integer, intent(in) :: positions(:) !! Positions selected from the source dimension.
      type(array_value), intent(inout) :: output !! Destination array receiving selected labels.
      integer, intent(in) :: output_dim !! One-based destination dimension receiving labels.
      integer, intent(out) :: status !! Zero on success; propagated SET_DIMNAMES error otherwise.
      character(len=:), allocatable :: labels(:)
      integer :: width
      integer :: i

      width = max(1, len(x%dimnames(source_dim)%values))
      allocate(character(len=width) :: labels(size(positions)))
      do i = 1, size(positions)
         labels(i) = x%dimnames(source_dim)%values(positions(i))
      end do
      call set_dimnames(output, output_dim, labels, status)
   end subroutine subset_names

   subroutine match_value_dimension(x, x_dim, value, value_dim, allow_excess, x_positions, value_positions, status)
      type(array_value), intent(in) :: x !! Destination array whose labels determine assignment positions.
      integer, intent(in) :: x_dim !! One-based destination dimension corresponding to VALUE_DIM.
      type(array_value), intent(in) :: value !! RHS array supplying labels and data.
      integer, intent(in) :: value_dim !! One-based RHS dimension whose labels are matched.
      logical, intent(in) :: allow_excess !! If true, unmatched RHS labels are discarded instead of causing an error.
      type(int_vector), intent(out) :: x_positions !! Matched one-based positions in the destination dimension.
      type(int_vector), intent(out) :: value_positions !! RHS positions retained after matching.
      integer, intent(out) :: status !! Zero on success; nonzero for absent labels or unmatched names when disallowed.
      integer, allocatable :: x_tmp(:)
      integer, allocatable :: v_tmp(:)
      integer :: i
      integer :: j
      integer :: n
      integer :: found

      status = 0
      if (.not. has_dimnames(x, x_dim) .or. .not. has_dimnames(value, value_dim)) then
         status = 5
         return
      end if
      allocate(x_tmp(value%dims(value_dim)), v_tmp(value%dims(value_dim)))
      n = 0
      do i = 1, value%dims(value_dim)
         found = 0
         do j = 1, x%dims(x_dim)
            if (trim(value%dimnames(value_dim)%values(i)) == trim(x%dimnames(x_dim)%values(j))) then
               found = j
               exit
            end if
         end do
         if (found == 0) then
            if (.not. allow_excess) then
               status = 6
               return
            end if
         else
            n = n + 1
            x_tmp(n) = found
            v_tmp(n) = i
         end if
      end do
      allocate(x_positions%values(n), value_positions%values(n))
      if (n > 0) then
         x_positions%values = x_tmp(1:n)
         value_positions%values = v_tmp(1:n)
      end if
   end subroutine match_value_dimension

   pure integer function expanded_to_source_dim(source_rank, expanded_rank, join_dim, expanded_dim) result(source_dim)
      integer, intent(in) :: source_rank !! Original rank of one input array.
      integer, intent(in) :: expanded_rank !! Common rank after optional singleton insertion.
      integer, intent(in) :: join_dim !! One-based join dimension in expanded coordinates.
      integer, intent(in) :: expanded_dim !! Expanded dimension to map back to the original input.

      if (source_rank == expanded_rank) then
         source_dim = expanded_dim
      else if (expanded_dim < join_dim) then
         source_dim = expanded_dim
      else if (expanded_dim == join_dim) then
         source_dim = 0
      else
         source_dim = expanded_dim - 1
      end if
   end function expanded_to_source_dim

   subroutine copy_dimension_name(source, source_dim, destination, destination_dim)
      type(array_value), intent(in) :: source !! Array providing a name for one dimension.
      integer, intent(in) :: source_dim !! One-based dimension whose name is copied.
      type(array_value), intent(inout) :: destination !! Array receiving the dimension name.
      integer, intent(in) :: destination_dim !! One-based destination dimension to name.
      integer :: width

      if (.not. allocated(source%dimname_names)) return
      if (source_dim > size(source%dimname_names)) return
      if (.not. allocated(destination%dimname_names)) return
      if (destination_dim > size(destination%dimname_names)) return
      width = max(len(destination%dimname_names), len(source%dimname_names))
      if (len(destination%dimname_names) < width) call widen_dimname_names(destination, width)
      destination%dimname_names(destination_dim) = source%dimname_names(source_dim)
   end subroutine copy_dimension_name

   subroutine widen_dimname_names(x, width)
      type(array_value), intent(inout) :: x !! Array whose dimension-name character width is increased.
      integer, intent(in) :: width !! New fixed character width for all dimension names.
      character(len=:), allocatable :: tmp(:)
      integer :: n

      if (.not. allocated(x%dimname_names)) return
      n = size(x%dimname_names)
      allocate(character(len=width) :: tmp(n))
      tmp = x%dimname_names
      call move_alloc(tmp, x%dimname_names)
   end subroutine widen_dimname_names

   pure integer function product_or_one(dims) result(value)
      integer, intent(in) :: dims(:) !! Dimension lengths whose product is needed; rank-zero gives one.

      if (size(dims) == 0) then
         value = 1
      else
         value = product(dims)
      end if
   end function product_or_one

end module abind_core
