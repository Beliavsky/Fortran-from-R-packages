module magic_tensor
   use magic_kinds, only : ik
   use magic_status, only : magic_error, set_error, MAGIC_INVALID_ARGUMENT
   implicit none
   private

   integer, parameter, public :: PAD_EXTEND = 1
   integer, parameter, public :: PAD_MIRROR = 2
   integer, parameter, public :: PAD_REPEAT = 3

   type, public :: integer_tensor
      integer, allocatable :: shape(:)
      integer(ik), allocatable :: values(:)
   contains
      procedure :: rank => tensor_rank
      procedure :: size => tensor_size
      procedure :: get => tensor_get
      procedure :: set => tensor_set
      procedure :: valid => tensor_valid
   end type integer_tensor

   public :: make_tensor, sequence_tensor, tensor_from_matrix, tensor_to_matrix
   public :: tensor_offset, unravel_index, tensor_permute, tensor_reverse
   public :: tensor_shift, tensor_rotate, tensor_block_diag, tensor_overlay_add
   public :: tensor_axis_coordinates, tensor_subsums, tensor_pad
   public :: tensor_take, tensor_drop, first_nonsingleton_dimensions
   public :: tensor_equal, tensor_is_circulant

contains

   pure integer function product_shape(shape) result(n)
      integer, intent(in) :: shape(:)
      integer :: k
      n = 1
      do k = 1, size(shape)
         if (shape(k) < 0) then
            n = -1
            return
         end if
         n = n * shape(k)
      end do
   end function product_shape

   function make_tensor(shape, fill, err) result(a)
      integer, intent(in) :: shape(:)
      integer(ik), intent(in), optional :: fill
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: a
      integer :: n
      integer(ik) :: value

      n = product_shape(shape)
      if (n < 0) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "tensor dimensions must be nonnegative")
         return
      end if
      value = 0_ik
      if (present(fill)) value = fill
      allocate(a%shape(size(shape)))
      a%shape = shape
      allocate(a%values(n))
      a%values = value
   end function make_tensor

   function sequence_tensor(shape, start, step, err) result(a)
      integer, intent(in) :: shape(:)
      integer(ik), intent(in), optional :: start, step
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: a
      integer(ik) :: first, increment
      integer :: i

      first = 1_ik
      increment = 1_ik
      if (present(start)) first = start
      if (present(step)) increment = step
      a = make_tensor(shape, err=err)
      if (.not. allocated(a%values)) return
      do i = 1, size(a%values)
         a%values(i) = first + int(i - 1, ik) * increment
      end do
   end function sequence_tensor

   function tensor_from_matrix(m) result(a)
      integer(ik), intent(in) :: m(:, :)
      type(integer_tensor) :: a
      a = make_tensor([size(m, 1), size(m, 2)])
      a%values = reshape(m, [size(m)])
   end function tensor_from_matrix

   function tensor_to_matrix(a, err) result(m)
      type(integer_tensor), intent(in) :: a
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: m(:, :)
      if (a%rank() /= 2) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "tensor_to_matrix requires rank two")
         allocate(m(0, 0))
         return
      end if
      allocate(m(a%shape(1), a%shape(2)))
      m = reshape(a%values, [a%shape(1), a%shape(2)])
   end function tensor_to_matrix

   pure integer function tensor_rank(self) result(r)
      class(integer_tensor), intent(in) :: self
      if (allocated(self%shape)) then
         r = size(self%shape)
      else
         r = 0
      end if
   end function tensor_rank

   pure integer function tensor_size(self) result(n)
      class(integer_tensor), intent(in) :: self
      if (allocated(self%values)) then
         n = size(self%values)
      else
         n = 0
      end if
   end function tensor_size

   pure logical function tensor_valid(self) result(ok)
      class(integer_tensor), intent(in) :: self
      ok = allocated(self%shape) .and. allocated(self%values)
      if (ok) ok = product_shape(self%shape) == size(self%values)
   end function tensor_valid

   pure integer function tensor_offset(shape, index) result(offset)
      integer, intent(in) :: shape(:), index(:)
      integer :: k, stride
      if (size(shape) /= size(index)) then
         offset = 0
         return
      end if
      offset = 1
      stride = 1
      do k = 1, size(shape)
         if (index(k) < 1 .or. index(k) > shape(k)) then
            offset = 0
            return
         end if
         offset = offset + (index(k) - 1) * stride
         stride = stride * shape(k)
      end do
   end function tensor_offset

   pure subroutine unravel_index(shape, offset, index)
      integer, intent(in) :: shape(:), offset
      integer, intent(out) :: index(:)
      integer :: k, q
      q = offset - 1
      do k = 1, size(shape)
         if (shape(k) > 0) then
            index(k) = modulo(q, shape(k)) + 1
            q = q / shape(k)
         else
            index(k) = 0
         end if
      end do
   end subroutine unravel_index

   pure integer(ik) function tensor_get(self, index) result(value)
      class(integer_tensor), intent(in) :: self
      integer, intent(in) :: index(:)
      integer :: offset
      offset = tensor_offset(self%shape, index)
      if (offset > 0) then
         value = self%values(offset)
      else
         value = 0_ik
      end if
   end function tensor_get

   subroutine tensor_set(self, index, value)
      class(integer_tensor), intent(inout) :: self
      integer, intent(in) :: index(:)
      integer(ik), intent(in) :: value
      integer :: offset
      offset = tensor_offset(self%shape, index)
      if (offset > 0) self%values(offset) = value
   end subroutine tensor_set

   function tensor_permute(a, order, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: order(:)
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer, allocatable :: old_index(:), new_index(:), seen(:)
      integer :: i, k, old_offset

      if (size(order) /= a%rank()) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "permutation rank mismatch")
         return
      end if
      allocate(seen(a%rank()), source=0)
      do k = 1, size(order)
         if (order(k) < 1 .or. order(k) > a%rank()) then
            call set_error(err, MAGIC_INVALID_ARGUMENT, "invalid permutation")
            return
         end if
         seen(order(k)) = seen(order(k)) + 1
      end do
      if (any(seen /= 1)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "invalid permutation")
         return
      end if
      out = make_tensor(a%shape(order))
      allocate(old_index(a%rank()), new_index(a%rank()))
      do i = 1, out%size()
         call unravel_index(out%shape, i, new_index)
         do k = 1, a%rank()
            old_index(order(k)) = new_index(k)
         end do
         old_offset = tensor_offset(a%shape, old_index)
         out%values(i) = a%values(old_offset)
      end do
   end function tensor_permute

   function tensor_reverse(a, axes, err) result(out)
      type(integer_tensor), intent(in) :: a
      logical, intent(in), optional :: axes(:)
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      logical, allocatable :: reverse_axis(:)
      integer, allocatable :: src(:), dst(:)
      integer :: i, k

      allocate(reverse_axis(a%rank()), source=.true.)
      if (present(axes)) then
         if (size(axes) /= a%rank()) then
            call set_error(err, MAGIC_INVALID_ARGUMENT, "reverse axes rank mismatch")
            return
         end if
         reverse_axis = axes
      end if
      out = make_tensor(a%shape)
      allocate(src(a%rank()), dst(a%rank()))
      do i = 1, a%size()
         call unravel_index(a%shape, i, dst)
         src = dst
         do k = 1, a%rank()
            if (reverse_axis(k)) src(k) = a%shape(k) - dst(k) + 1
         end do
         out%values(i) = a%get(src)
      end do
   end function tensor_reverse

   function tensor_shift(a, shifts, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: shifts(:)
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer, allocatable :: src(:), dst(:), use_shifts(:)
      integer :: i, k

      if (size(shifts) > a%rank()) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "too many tensor shifts")
         out = make_tensor([0])
         return
      end if
      allocate(use_shifts(a%rank()), source=0)
      use_shifts(1:size(shifts)) = shifts
      out = make_tensor(a%shape)
      allocate(src(a%rank()), dst(a%rank()))
      do i = 1, a%size()
         call unravel_index(a%shape, i, dst)
         do k = 1, a%rank()
            if (a%shape(k) > 0) then
               src(k) = modulo(dst(k) - 1 - use_shifts(k), a%shape(k)) + 1
            else
               src(k) = 0
            end if
         end do
         if (a%size() > 0) out%values(i) = a%get(src)
      end do
   end function tensor_shift

   function tensor_rotate(a, rights, pair, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in), optional :: rights
      integer, intent(in), optional :: pair(2)
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out, tmp
      integer :: r, p(2), i
      integer, allocatable :: order(:)
      logical, allocatable :: axes(:)

      r = 1
      if (present(rights)) r = modulo(rights, 4)
      p = [1, 2]
      if (present(pair)) p = pair
      if (a%rank() < 2 .or. any(p < 1) .or. any(p > a%rank()) .or. p(1) == p(2)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "invalid rotation plane")
         return
      end if
      if (r == 0) then
         out = a
         return
      end if
      allocate(order(a%rank()))
      order = [(i, i=1,a%rank())]
      allocate(axes(a%rank()), source=.false.)
      select case (r)
      case (1)
         axes(p(2)) = .true.
         tmp = tensor_reverse(a, axes)
         order(p(1)) = p(2)
         order(p(2)) = p(1)
         out = tensor_permute(tmp, order)
      case (2)
         axes(p(1)) = .true.
         axes(p(2)) = .true.
         out = tensor_reverse(a, axes)
      case (3)
         axes(p(1)) = .true.
         tmp = tensor_reverse(a, axes)
         order(p(1)) = p(2)
         order(p(2)) = p(1)
         out = tensor_permute(tmp, order)
      end select
   end function tensor_rotate

   function tensor_block_diag(a, b, pad, err) result(out)
      type(integer_tensor), intent(in) :: a, b
      integer(ik), intent(in), optional :: pad
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer(ik) :: fill
      integer, allocatable :: idx(:), out_idx(:)
      integer :: i

      if (a%rank() /= b%rank()) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "block-diagonal tensors need equal rank")
         return
      end if
      fill = 0_ik
      if (present(pad)) fill = pad
      out = make_tensor(a%shape + b%shape, fill)
      allocate(idx(a%rank()), out_idx(a%rank()))
      do i = 1, a%size()
         call unravel_index(a%shape, i, idx)
         call out%set(idx, a%values(i))
      end do
      do i = 1, b%size()
         call unravel_index(b%shape, i, idx)
         out_idx = idx + a%shape
         call out%set(out_idx, b%values(i))
      end do
   end function tensor_block_diag

   function tensor_overlay_add(a, b, err) result(out)
      type(integer_tensor), intent(in) :: a, b
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer, allocatable :: idx(:)
      integer :: i, offset

      if (a%rank() /= b%rank()) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "overlay tensors need equal rank")
         return
      end if
      out = make_tensor(max(a%shape, b%shape))
      allocate(idx(a%rank()))
      do i = 1, a%size()
         call unravel_index(a%shape, i, idx)
         offset = tensor_offset(out%shape, idx)
         out%values(offset) = out%values(offset) + a%values(i)
      end do
      do i = 1, b%size()
         call unravel_index(b%shape, i, idx)
         offset = tensor_offset(out%shape, idx)
         out%values(offset) = out%values(offset) + b%values(i)
      end do
   end function tensor_overlay_add

   function tensor_axis_coordinates(shape, axis, err) result(out)
      integer, intent(in) :: shape(:), axis
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer, allocatable :: idx(:)
      integer :: i
      if (axis < 1 .or. axis > size(shape)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "invalid coordinate axis")
         return
      end if
      out = make_tensor(shape)
      allocate(idx(size(shape)))
      do i = 1, out%size()
         call unravel_index(shape, i, idx)
         out%values(i) = int(idx(axis), ik)
      end do
   end function tensor_axis_coordinates

   function tensor_subsums(a, window, wrap, pad, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: window(:)
      logical, intent(in), optional :: wrap
      integer(ik), intent(in), optional :: pad
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      logical :: do_wrap
      integer(ik) :: pad_value, total
      integer, allocatable :: center(:), offset_idx(:), src(:), offset_shape(:)
      integer :: i, j, k, n_offsets
      logical :: inside

      if (size(window) /= a%rank() .or. any(window < 1)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "invalid subarray window")
         return
      end if
      do_wrap = .true.
      if (present(wrap)) do_wrap = wrap
      pad_value = 0_ik
      if (present(pad)) pad_value = pad
      out = make_tensor(a%shape)
      allocate(center(a%rank()), offset_idx(a%rank()), src(a%rank()))
      allocate(offset_shape(a%rank()))
      offset_shape = window
      n_offsets = product_shape(offset_shape)
      do i = 1, a%size()
         call unravel_index(a%shape, i, center)
         total = 0_ik
         do j = 1, n_offsets
            call unravel_index(offset_shape, j, offset_idx)
            src = center - (offset_idx - 1)
            inside = .true.
            do k = 1, a%rank()
               if (do_wrap) then
                  src(k) = modulo(src(k) - 1, a%shape(k)) + 1
               else if (src(k) < 1 .or. src(k) > a%shape(k)) then
                  inside = .false.
               end if
            end do
            if (inside) then
               total = total + a%get(src)
            else
               total = total + pad_value
            end if
         end do
         out%values(i) = total
      end do
   end function tensor_subsums

   pure integer function mapped_pad_index(i, n, amount, method, post) result(index)
      integer, intent(in) :: i, n, amount, method
      logical, intent(in) :: post
      integer :: position, period
      if (post) then
         position = i
      else
         position = i - amount
      end if
      select case (method)
      case (PAD_EXTEND)
         index = min(max(position, 1), n)
      case (PAD_REPEAT)
         index = modulo(position - 1, n) + 1
      case (PAD_MIRROR)
         if (n <= 1) then
            index = 1
         else
            period = 2 * n
            position = modulo(position - 1, period) + 1
            if (position <= n) then
               index = position
            else
               index = period - position + 1
            end if
         end if
      case default
         index = min(max(position, 1), n)
      end select
   end function mapped_pad_index

   function tensor_pad(a, amounts, method, post, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: amounts(:)
      integer, intent(in), optional :: method
      logical, intent(in), optional :: post
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer :: pad_method, i, k
      logical :: at_end
      integer, allocatable :: dst(:), src(:)

      if (size(amounts) /= a%rank() .or. any(amounts < 0)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "invalid padding amounts")
         return
      end if
      pad_method = PAD_EXTEND
      if (present(method)) pad_method = method
      at_end = .true.
      if (present(post)) at_end = post
      out = make_tensor(a%shape + amounts)
      allocate(dst(a%rank()), src(a%rank()))
      do i = 1, out%size()
         call unravel_index(out%shape, i, dst)
         do k = 1, a%rank()
            src(k) = mapped_pad_index(dst(k), a%shape(k), amounts(k), pad_method, at_end)
         end do
         out%values(i) = a%get(src)
      end do
   end function tensor_pad

   function tensor_take(a, counts, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: counts(:)
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer, allocatable :: use_counts(:), src(:), dst(:)
      integer :: i, k

      if (size(counts) > a%rank()) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "too many take counts")
         return
      end if
      allocate(use_counts(a%rank()))
      use_counts = a%shape
      use_counts(1:size(counts)) = counts
      if (any(abs(use_counts) > a%shape)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "take count exceeds dimension")
         return
      end if
      out = make_tensor(abs(use_counts))
      allocate(src(a%rank()), dst(a%rank()))
      do i = 1, out%size()
         call unravel_index(out%shape, i, dst)
         do k = 1, a%rank()
            if (use_counts(k) >= 0) then
               src(k) = dst(k)
            else
               src(k) = a%shape(k) - out%shape(k) + dst(k)
            end if
         end do
         out%values(i) = a%get(src)
      end do
   end function tensor_take

   function tensor_drop(a, counts, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: counts(:)
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      integer, allocatable :: use_counts(:), new_shape(:), src(:), dst(:)
      integer :: i, k

      if (size(counts) > a%rank()) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "too many drop counts")
         return
      end if
      allocate(use_counts(a%rank()), source=0)
      use_counts(1:size(counts)) = counts
      if (any(abs(use_counts) > a%shape)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "drop count exceeds dimension")
         return
      end if
      allocate(new_shape(a%rank()))
      new_shape = a%shape - abs(use_counts)
      out = make_tensor(new_shape)
      allocate(src(a%rank()), dst(a%rank()))
      do i = 1, out%size()
         call unravel_index(out%shape, i, dst)
         do k = 1, a%rank()
            if (use_counts(k) > 0) then
               src(k) = dst(k) + use_counts(k)
            else
               src(k) = dst(k)
            end if
         end do
         out%values(i) = a%get(src)
      end do
   end function tensor_drop

   function first_nonsingleton_dimensions(a, n) result(dimensions)
      type(integer_tensor), intent(in) :: a
      integer, intent(in), optional :: n
      integer, allocatable :: dimensions(:)
      integer :: wanted, nfound, k
      wanted = 1
      if (present(n)) wanted = max(0, n)
      nfound = min(wanted, count(a%shape > 1))
      allocate(dimensions(nfound))
      nfound = 0
      do k = 1, a%rank()
         if (a%shape(k) > 1 .and. nfound < size(dimensions)) then
            nfound = nfound + 1
            dimensions(nfound) = k
         end if
      end do
   end function first_nonsingleton_dimensions

   pure logical function tensor_equal(a, b) result(equal)
      type(integer_tensor), intent(in) :: a, b
      equal = a%rank() == b%rank()
      if (equal) equal = all(a%shape == b%shape)
      if (equal) equal = all(a%values == b%values)
   end function tensor_equal

   function tensor_is_circulant(a, directions) result(answer)
      type(integer_tensor), intent(in) :: a
      integer, intent(in), optional :: directions(:)
      logical :: answer
      integer, allocatable :: shifts(:)
      type(integer_tensor) :: moved
      allocate(shifts(a%rank()), source=1)
      if (present(directions)) then
         shifts = 1
         shifts(1:min(size(directions), a%rank())) = directions(1:min(size(directions), a%rank()))
      end if
      moved = tensor_shift(a, shifts)
      answer = tensor_equal(a, moved)
   end function tensor_is_circulant

end module magic_tensor
