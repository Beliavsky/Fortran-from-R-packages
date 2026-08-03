module magic_hypercube
   use magic_kinds, only : ik
   use magic_status, only : magic_error, set_error, MAGIC_INVALID_ARGUMENT
   use magic_tensor, only : integer_tensor, make_tensor, tensor_axis_coordinates, &
                            tensor_offset, unravel_index
   implicit none
   private

   public :: magiccube_2np1, magichypercube_4n
   public :: is_semimagichypercube, is_diagonally_correct
   public :: is_magichypercube, is_latinhypercube, is_perfect_hypercube
   public :: diagonal_subhypercubes, is_alicehypercube

contains

   function magiccube_2np1(m, err) result(cube)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: cube
      integer :: n, i
      integer, allocatable :: index(:)
      integer(ik) :: x, y, z
      if (m < 0) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "m must be nonnegative")
         return
      end if
      n = 2 * m + 1
      cube = make_tensor([n, n, n])
      allocate(index(3))
      do i = 1, cube%size()
         call unravel_index(cube%shape, i, index)
         x = int(index(1), ik)
         y = int(index(2), ik)
         z = int(index(3), ik)
         cube%values(i) = modulo(x - y + z - 1_ik, int(n, ik)) * int(n * n, ik) + &
                          modulo(x - y - z, int(n, ik)) * int(n, ik) + &
                          modulo(x + y + z - 2_ik, int(n, ik)) + 1_ik
      end do
   end function magiccube_2np1

   function magichypercube_4n(m, dimension, err) result(cube)
      integer, intent(in) :: m
      integer, intent(in), optional :: dimension
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: cube
      integer :: n, d, i, k, selected_count, position
      integer, allocatable :: index(:), selected(:)
      integer :: parity_count

      d = 3
      if (present(dimension)) d = dimension
      if (m < 1 .or. d < 2) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "m must be positive and dimension at least two")
         return
      end if
      n = 4 * m
      cube = make_tensor([(n, k=1,d)])
      do i = 1, cube%size()
         cube%values(i) = int(i, ik)
      end do
      allocate(index(d), selected(cube%size()))
      selected_count = 0
      do i = 1, cube%size()
         call unravel_index(cube%shape, i, index)
         parity_count = d
         do k = 1, d
            position = modulo(index(k) - 1, 4)
            if (position == 1 .or. position == 2) parity_count = parity_count + 1
         end do
         if (modulo(parity_count, 2) == 1) then
            selected_count = selected_count + 1
            selected(selected_count) = i
         end if
      end do
      do i = 1, selected_count / 2
         position = selected(selected_count - i + 1)
         k = selected(i)
         call swap_values(cube%values(k), cube%values(position))
      end do
   end function magichypercube_4n

   pure subroutine swap_values(a, b)
      integer(ik), intent(inout) :: a, b
      integer(ik) :: temp
      temp = a
      a = b
      b = temp
   end subroutine swap_values

   function first_line_sum(a) result(total)
      type(integer_tensor), intent(in) :: a
      integer(ik) :: total
      integer, allocatable :: index(:)
      integer :: k
      allocate(index(a%rank()), source=1)
      total = 0_ik
      do k = 1, a%shape(1)
         index(1) = k
         total = total + a%get(index)
      end do
   end function first_line_sum

   function semimagic_with_target(a, target) result(answer)
      type(integer_tensor), intent(in) :: a
      integer(ik), intent(in) :: target
      logical :: answer
      integer, allocatable :: base_index(:), index(:)
      integer :: axis, i, k
      integer(ik) :: total

      if (a%rank() < 1 .or. any(a%shape /= a%shape(1))) then
         answer = .false.
         return
      end if
      allocate(base_index(a%rank()), index(a%rank()))
      answer = .true.
      do axis = 1, a%rank()
         do i = 1, a%size()
            call unravel_index(a%shape, i, base_index)
            if (base_index(axis) /= 1) cycle
            index = base_index
            total = 0_ik
            do k = 1, a%shape(axis)
               index(axis) = k
               total = total + a%get(index)
            end do
            if (total /= target) then
               answer = .false.
               return
            end if
         end do
      end do
   end function semimagic_with_target

   function is_semimagichypercube(a) result(answer)
      type(integer_tensor), intent(in) :: a
      logical :: answer
      integer(ik) :: target
      if (.not. a%valid() .or. a%rank() < 1 .or. any(a%shape /= a%shape(1))) then
         answer = .false.
         return
      end if
      target = first_line_sum(a)
      answer = semimagic_with_target(a, target)
   end function is_semimagichypercube

   function is_diagonally_correct(a) result(answer)
      type(integer_tensor), intent(in) :: a
      logical :: answer
      integer, allocatable :: index(:)
      integer :: pattern, patterns, k, t, n
      integer(ik) :: total, target
      logical :: have_target

      if (.not. a%valid() .or. a%rank() < 1 .or. any(a%shape /= a%shape(1))) then
         answer = .false.
         return
      end if
      n = a%shape(1)
      patterns = 2 ** a%rank()
      allocate(index(a%rank()))
      have_target = .false.
      target = 0_ik
      answer = .true.
      do pattern = 0, patterns - 1
         total = 0_ik
         do t = 1, n
            do k = 1, a%rank()
               if (btest(pattern, k - 1)) then
                  index(k) = n - t + 1
               else
                  index(k) = t
               end if
            end do
            total = total + a%get(index)
         end do
         if (.not. have_target) then
            target = total
            have_target = .true.
         else if (total /= target) then
            answer = .false.
            return
         end if
      end do
   end function is_diagonally_correct

   function is_magichypercube(a) result(answer)
      type(integer_tensor), intent(in) :: a
      logical :: answer
      answer = is_semimagichypercube(a) .and. is_diagonally_correct(a)
   end function is_magichypercube

   subroutine sort_integer(x)
      integer(ik), intent(inout) :: x(:)
      integer :: i, j
      integer(ik) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_integer

   function is_latinhypercube(a) result(answer)
      type(integer_tensor), intent(in) :: a
      logical :: answer
      integer, allocatable :: base_index(:), index(:)
      integer(ik), allocatable :: line(:), expected(:)
      integer :: axis, i, k, n
      if (.not. a%valid() .or. a%rank() < 1 .or. any(a%shape /= a%shape(1))) then
         answer = .false.
         return
      end if
      n = a%shape(1)
      allocate(base_index(a%rank()), index(a%rank()), line(n), expected(n))
      expected = [(int(k, ik), k=1,n)]
      answer = .true.
      do axis = 1, a%rank()
         do i = 1, a%size()
            call unravel_index(a%shape, i, base_index)
            if (base_index(axis) /= 1) cycle
            index = base_index
            do k = 1, n
               index(axis) = k
               line(k) = a%get(index)
            end do
            call sort_integer(line)
            if (any(line /= expected)) then
               answer = .false.
               return
            end if
         end do
      end do
   end function is_latinhypercube

   function diagonal_subhypercubes(a, err) result(subcubes)
      type(integer_tensor), intent(in) :: a
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor), allocatable :: subcubes(:)
      integer :: total_patterns, pattern_number, d, n, constrained, free_count
      integer :: k, count_patterns, output_rank, i, t, free_position
      integer, allocatable :: pattern(:), free_axes(:), out_index(:), src_index(:), out_shape(:)

      if (.not. a%valid() .or. a%rank() < 2 .or. any(a%shape /= a%shape(1))) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "equal-sided tensor of rank at least two required")
         allocate(subcubes(0))
         return
      end if
      d = a%rank()
      n = a%shape(1)
      total_patterns = 3 ** d
      count_patterns = 0
      allocate(pattern(d))
      do pattern_number = 0, total_patterns - 1
         call decode_pattern(pattern_number, pattern)
         constrained = count(pattern /= 0)
         if (constrained >= 2) count_patterns = count_patterns + 1
      end do
      allocate(subcubes(count_patterns))
      count_patterns = 0
      do pattern_number = 0, total_patterns - 1
         call decode_pattern(pattern_number, pattern)
         constrained = count(pattern /= 0)
         if (constrained < 2) cycle
         count_patterns = count_patterns + 1
         free_count = d - constrained
         output_rank = free_count + 1
         allocate(out_shape(output_rank), source=n)
         subcubes(count_patterns) = make_tensor(out_shape)
         allocate(free_axes(free_count), out_index(output_rank), src_index(d))
         free_position = 0
         do k = 1, d
            if (pattern(k) == 0) then
               free_position = free_position + 1
               free_axes(free_position) = k
            end if
         end do
         do i = 1, subcubes(count_patterns)%size()
            call unravel_index(subcubes(count_patterns)%shape, i, out_index)
            t = out_index(output_rank)
            free_position = 0
            do k = 1, d
               select case (pattern(k))
               case (0)
                  free_position = free_position + 1
                  src_index(k) = out_index(free_position)
               case (1)
                  src_index(k) = t
               case (-1)
                  src_index(k) = n - t + 1
               end select
            end do
            subcubes(count_patterns)%values(i) = a%get(src_index)
         end do
         deallocate(out_shape, free_axes, out_index, src_index)
      end do
   contains
      subroutine decode_pattern(number, decoded)
         integer, intent(in) :: number
         integer, intent(out) :: decoded(:)
         integer :: q, kk, dd
         q = number
         do kk = 1, size(decoded)
            dd = modulo(q, 3)
            q = q / 3
            select case (dd)
            case (0)
               decoded(kk) = 0
            case (1)
               decoded(kk) = 1
            case (2)
               decoded(kk) = -1
            end select
         end do
      end subroutine decode_pattern
   end function diagonal_subhypercubes

   function is_perfect_hypercube(a) result(answer)
      type(integer_tensor), intent(in) :: a
      logical :: answer
      type(integer_tensor), allocatable :: subcubes(:)
      integer(ik) :: target
      integer :: i
      if (.not. is_semimagichypercube(a)) then
         answer = .false.
         return
      end if
      target = first_line_sum(a)
      subcubes = diagonal_subhypercubes(a)
      answer = .true.
      do i = 1, size(subcubes)
         if (subcubes(i)%rank() == 1) then
            if (sum(subcubes(i)%values) /= target) then
               answer = .false.
               return
            end if
         else if (.not. semimagic_with_target(subcubes(i), target)) then
            answer = .false.
            return
         end if
      end do
   end function is_perfect_hypercube

   function is_alicehypercube(a, ndim) result(answer)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: ndim
      logical :: answer
      integer :: d, keep_count, mask, max_mask, i, j, k, retained_size
      integer, allocatable :: full_index(:), retained_shape(:), retained_index(:), retained_axes(:)
      integer(ik) :: total, target
      logical :: have_target, matches

      d = a%rank()
      keep_count = d - ndim
      if (.not. a%valid() .or. keep_count < 0 .or. keep_count > d) then
         answer = .false.
         return
      end if
      max_mask = 2 ** d
      allocate(full_index(d))
      have_target = .false.
      answer = .true.
      do mask = 0, max_mask - 1
         if (popcnt(mask) /= keep_count) cycle
         allocate(retained_axes(keep_count), retained_shape(keep_count), retained_index(keep_count))
         j = 0
         do k = 1, d
            if (btest(mask, k - 1)) then
               j = j + 1
               retained_axes(j) = k
               retained_shape(j) = a%shape(k)
            end if
         end do
         retained_size = 1
         if (keep_count > 0) retained_size = product(retained_shape)
         do i = 1, retained_size
            if (keep_count > 0) call unravel_index(retained_shape, i, retained_index)
            total = 0_ik
            do j = 1, a%size()
               call unravel_index(a%shape, j, full_index)
               matches = .true.
               do k = 1, keep_count
                  matches = matches .and. full_index(retained_axes(k)) == retained_index(k)
               end do
               if (matches) total = total + a%values(j)
            end do
            if (.not. have_target) then
               target = total
               have_target = .true.
            else if (total /= target) then
               answer = .false.
               return
            end if
         end do
         deallocate(retained_axes, retained_shape, retained_index)
      end do
   end function is_alicehypercube

end module magic_hypercube
