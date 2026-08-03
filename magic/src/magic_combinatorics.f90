module magic_combinatorics
   use magic_kinds, only : ik, dp
   use magic_status, only : magic_error, set_error, MAGIC_INVALID_ARGUMENT
   use magic_tensor, only : integer_tensor, make_tensor, tensor_equal
   use magic_square, only : latin_square, is_latin_square
   implicit none
   private

   public :: sylvester_hadamard, is_hadamard, cilleruelo_square
   public :: bernhardsson_a, bernhardsson_b, bernhardsson_matrix
   public :: incidence, unincidence, is_incidence, incidence_move
   public :: another_incidence, another_latin, random_latin_squares
   public :: sam_square

contains

   function kronecker_integer(a, b) result(out)
      integer(ik), intent(in) :: a(:, :), b(:, :)
      integer(ik), allocatable :: out(:, :)
      integer :: i, j, p, q
      allocate(out(size(a, 1) * size(b, 1), size(a, 2) * size(b, 2)))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            do q = 1, size(b, 2)
               do p = 1, size(b, 1)
                  out((i - 1) * size(b, 1) + p, (j - 1) * size(b, 2) + q) = a(i, j) * b(p, q)
               end do
            end do
         end do
      end do
   end function kronecker_integer

   recursive function sylvester_hadamard(k, err) result(h)
      integer, intent(in) :: k
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: h(:, :)
      integer(ik), allocatable :: previous(:, :), base(:, :)
      if (k < 0) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "k must be nonnegative")
         allocate(h(0, 0))
      else if (k == 0) then
         allocate(h(1, 1), source=1_ik)
      else
         previous = sylvester_hadamard(k - 1, err)
         allocate(base(2, 2))
         base = reshape([1_ik, 1_ik, 1_ik, -1_ik], [2, 2])
         h = kronecker_integer(previous, base)
      end if
   end function sylvester_hadamard

   function is_hadamard(h) result(answer)
      integer(ik), intent(in) :: h(:, :)
      logical :: answer
      integer(ik), allocatable :: gram(:, :)
      integer :: i, n
      if (size(h, 1) /= size(h, 2) .or. any(abs(h) /= 1_ik)) then
         answer = .false.
         return
      end if
      n = size(h, 1)
      gram = matmul(transpose(h), h)
      answer = .true.
      do i = 1, n
         if (gram(i, i) /= int(n, ik)) answer = .false.
         gram(i, i) = 0_ik
      end do
      answer = answer .and. all(gram == 0_ik)
   end function is_hadamard

   function cilleruelo_square(n, m) result(a)
      integer, intent(in) :: n, m
      integer(ik), allocatable :: a(:, :)
      allocate(a(4, 4))
      a(1, :) = [int((n + 2) * m, ik), int((n + 3) * (m + 3), ik), &
                 int((n + 1) * (m + 2), ik), int(n * (m + 1), ik)]
      a(2, :) = [int((n + 1) * (m + 1), ik), int(n * (m + 2), ik), &
                 int((n + 2) * (m + 3), ik), int((n + 3) * m, ik)]
      a(3, :) = [int(n * (m + 3), ik), int((n + 1) * m, ik), &
                 int((n + 3) * (m + 1), ik), int((n + 2) * (m + 2), ik)]
      a(4, :) = [int((n + 3) * (m + 2), ik), int((n + 2) * (m + 1), ik), &
                 int(n * m, ik), int((n + 1) * (m + 3), ik)]
   end function cilleruelo_square

   recursive function bernhardsson_a(n, err) result(a)
      integer, intent(in) :: n
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      integer(ik), allocatable :: sub(:, :)
      integer :: j, half
      if (n < 1) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "n must be positive")
         allocate(a(0, 0))
      else if (modulo(n, 2) == 1) then
         sub = bernhardsson_a(n - 1, err)
         allocate(a(n, n), source=0_ik)
         a(1, 1) = 1_ik
         if (n > 1) a(2:n, 2:n) = sub
      else
         allocate(a(n, n), source=0_ik)
         half = n / 2
         do j = 1, half
            a(j, 2 * j) = 1_ik
            a(half + j, 2 * j - 1) = 1_ik
         end do
      end if
   end function bernhardsson_a

   recursive function bernhardsson_b(n, err) result(a)
      integer, intent(in) :: n
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      integer(ik), allocatable :: sub(:, :)
      integer :: j, half, column
      if (n < 1) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "n must be positive")
         allocate(a(0, 0))
      else if (modulo(n, 2) == 1) then
         sub = bernhardsson_b(n - 1, err)
         allocate(a(n, n), source=0_ik)
         a(1, 1) = 1_ik
         if (n > 1) a(2:n, 2:n) = sub
      else
         allocate(a(n, n), source=0_ik)
         half = n / 2
         do j = 1, half
            column = 1 + modulo(2 * (j - 1) + half - 1, n)
            a(j, column) = 1_ik
            column = n - modulo(2 * (j - 1) + half - 1, n)
            a(n + 1 - j, column) = 1_ik
         end do
      end if
   end function bernhardsson_b

   function bernhardsson_matrix(n, err) result(a)
      integer, intent(in) :: n
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      if (modulo(n, 6) == 0 .or. modulo(n, 6) == 1) then
         a = bernhardsson_a(n, err)
      else
         a = bernhardsson_b(n, err)
      end if
   end function bernhardsson_matrix

   function incidence(a, err) result(tensor)
      integer(ik), intent(in) :: a(:, :)
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: tensor
      integer :: nrow, ncol, symbols, i, j, symbol
      if (size(a) == 0 .or. minval(a) < 1_ik) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "incidence input must contain positive symbols")
         return
      end if
      nrow = size(a, 1)
      ncol = size(a, 2)
      symbols = int(maxval(a))
      tensor = make_tensor([nrow, ncol, symbols])
      do j = 1, ncol
         do i = 1, nrow
            symbol = int(a(i, j))
            call tensor%set([i, j, symbol], 1_ik)
         end do
      end do
   end function incidence

   function is_incidence(a, include_improper) result(answer)
      type(integer_tensor), intent(in) :: a
      logical, intent(in), optional :: include_improper
      logical :: answer, allow_improper, line_ok
      integer :: axis, i, k, minus_count
      integer, allocatable :: base_index(:), index(:)
      integer(ik), allocatable :: line(:)

      allow_improper = .false.
      if (present(include_improper)) allow_improper = include_improper
      if (.not. a%valid() .or. a%rank() /= 3 .or. any(a%shape /= a%shape(1))) then
         answer = .false.
         return
      end if
      minus_count = count(a%values == -1_ik)
      if (minus_count > 0 .and. (.not. allow_improper .or. minus_count /= 1)) then
         answer = .false.
         return
      end if
      if (any(a%values < -1_ik) .or. any(a%values > 1_ik)) then
         answer = .false.
         return
      end if
      allocate(base_index(3), index(3), line(a%shape(1)))
      answer = .true.
      do axis = 1, 3
         do i = 1, a%size()
            call unravel_local(a%shape, i, base_index)
            if (base_index(axis) /= 1) cycle
            index = base_index
            do k = 1, a%shape(axis)
               index(axis) = k
               line(k) = a%get(index)
            end do
            line_ok = sum(line) == 1_ik .and. all(line >= -1_ik) .and. all(line <= 1_ik)
            if (.not. allow_improper) line_ok = line_ok .and. all(line >= 0_ik)
            if (.not. line_ok) then
               answer = .false.
               return
            end if
         end do
      end do
   contains
      pure subroutine unravel_local(shape, offset, index)
         integer, intent(in) :: shape(:), offset
         integer, intent(out) :: index(:)
         integer :: q, kk
         q = offset - 1
         do kk = 1, size(shape)
            index(kk) = modulo(q, shape(kk)) + 1
            q = q / shape(kk)
         end do
      end subroutine unravel_local
   end function is_incidence

   function unincidence(tensor, err) result(a)
      type(integer_tensor), intent(in) :: tensor
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      integer :: i, j, k, found
      if (.not. is_incidence(tensor, .false.)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "proper incidence tensor required")
         allocate(a(0, 0))
         return
      end if
      allocate(a(tensor%shape(1), tensor%shape(2)))
      do j = 1, tensor%shape(2)
         do i = 1, tensor%shape(1)
            found = 0
            do k = 1, tensor%shape(3)
               if (tensor%get([i, j, k]) == 1_ik) then
                  found = k
                  exit
               end if
            end do
            a(i, j) = int(found, ik)
         end do
      end do
   end function unincidence

   integer function random_choice(indices) result(choice)
      integer, intent(in) :: indices(:)
      real(dp) :: u
      integer :: position
      call random_number(u)
      position = min(size(indices), int(u * real(size(indices), dp)) + 1)
      choice = indices(position)
   end function random_choice

   function positions_equal_value(a, axis, pivot, value) result(positions)
      type(integer_tensor), intent(in) :: a
      integer, intent(in) :: axis, pivot(3)
      integer(ik), intent(in) :: value
      integer, allocatable :: positions(:)
      integer, allocatable :: work(:)
      integer :: k, count_found
      allocate(work(a%shape(axis)))
      count_found = 0
      do k = 1, a%shape(axis)
         select case (axis)
         case (1)
            if (a%get([k, pivot(2), pivot(3)]) == value) then
               count_found = count_found + 1
               work(count_found) = k
            end if
         case (2)
            if (a%get([pivot(1), k, pivot(3)]) == value) then
               count_found = count_found + 1
               work(count_found) = k
            end if
         case (3)
            if (a%get([pivot(1), pivot(2), k]) == value) then
               count_found = count_found + 1
               work(count_found) = k
            end if
         end select
      end do
      allocate(positions(count_found))
      positions = work(1:count_found)
   end function positions_equal_value

   function incidence_move(a, err) result(out)
      type(integer_tensor), intent(in) :: a
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      logical :: proper
      integer, allocatable :: zero_positions(:), p1s(:), p2s(:), p3s(:)
      integer :: pivot(3), p1, p2, p3, i, count_zero

      if (.not. is_incidence(a, .true.)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "proper or one-minus-one incidence tensor required")
         return
      end if
      out = a
      proper = is_incidence(a, .false.)
      if (proper) then
         count_zero = count(a%values == 0_ik)
         allocate(zero_positions(count_zero))
         zero_positions = pack([(i, i=1,a%size())], a%values == 0_ik)
         i = random_choice(zero_positions)
         call unravel_local(a%shape, i, pivot)
      else
         i = findloc(a%values, -1_ik, dim=1)
         call unravel_local(a%shape, i, pivot)
      end if
      p1s = positions_equal_value(a, 1, pivot, 1_ik)
      p2s = positions_equal_value(a, 2, pivot, 1_ik)
      p3s = positions_equal_value(a, 3, pivot, 1_ik)
      if (size(p1s) == 0 .or. size(p2s) == 0 .or. size(p3s) == 0) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "incidence move could not find pivots")
         return
      end if
      if (proper) then
         p1 = p1s(1)
         p2 = p2s(1)
         p3 = p3s(1)
      else
         p1 = random_choice(p1s)
         p2 = random_choice(p2s)
         p3 = random_choice(p3s)
      end if
      call add_at(out, pivot, 1_ik)
      call add_at(out, [pivot(1), p2, p3], 1_ik)
      call add_at(out, [p1, pivot(2), p3], 1_ik)
      call add_at(out, [p1, p2, pivot(3)], 1_ik)
      call add_at(out, [p1, pivot(2), pivot(3)], -1_ik)
      call add_at(out, [pivot(1), p2, pivot(3)], -1_ik)
      call add_at(out, [pivot(1), pivot(2), p3], -1_ik)
      call add_at(out, [p1, p2, p3], -1_ik)
   contains
      pure subroutine unravel_local(shape, offset, index)
         integer, intent(in) :: shape(:), offset
         integer, intent(out) :: index(:)
         integer :: q, kk
         q = offset - 1
         do kk = 1, size(shape)
            index(kk) = modulo(q, shape(kk)) + 1
            q = q / shape(kk)
         end do
      end subroutine unravel_local
      subroutine add_at(tensor, index, increment)
         type(integer_tensor), intent(inout) :: tensor
         integer, intent(in) :: index(3)
         integer(ik), intent(in) :: increment
         call tensor%set(index, tensor%get(index) + increment)
      end subroutine add_at
   end function incidence_move

   function another_incidence(a, max_steps, err) result(out)
      type(integer_tensor), intent(in) :: a
      integer, intent(in), optional :: max_steps
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out, current
      integer :: limit, step
      if (.not. is_incidence(a, .false.)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "proper incidence tensor required")
         return
      end if
      limit = 100000
      if (present(max_steps)) limit = max_steps
      current = a
      do step = 1, limit
         current = incidence_move(current, err)
         if (is_incidence(current, .false.) .and. .not. tensor_equal(current, a)) then
            out = current
            return
         end if
      end do
      call set_error(err, MAGIC_INVALID_ARGUMENT, "no different proper incidence tensor found")
   end function another_incidence

   function another_latin(a, max_steps, err) result(out)
      integer(ik), intent(in) :: a(:, :)
      integer, intent(in), optional :: max_steps
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :)
      type(integer_tensor) :: inc, changed
      if (.not. is_latin_square(a)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "Latin square required")
         allocate(out(0, 0))
         return
      end if
      inc = incidence(a)
      changed = another_incidence(inc, max_steps, err)
      if (changed%valid()) then
         out = unincidence(changed, err)
      else
         allocate(out(0, 0))
      end if
   end function another_latin

   function random_latin_squares(n, count_squares, burnin, err) result(out)
      integer, intent(in) :: n
      integer, intent(in), optional :: count_squares, burnin
      type(magic_error), intent(inout), optional :: err
      type(integer_tensor) :: out
      type(integer_tensor) :: inc
      integer(ik), allocatable :: square(:, :)
      integer :: count_use, burn, i, j, k
      if (n < 1) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "n must be positive")
         return
      end if
      count_use = 1
      if (present(count_squares)) count_use = count_squares
      burn = n * n
      if (present(burnin)) burn = burnin
      square = latin_square(n)
      inc = incidence(square)
      do i = 1, burn
         inc = another_incidence(inc, err=err)
      end do
      out = make_tensor([n, n, count_use])
      do k = 1, count_use
         square = unincidence(inc, err)
         do j = 1, n
            do i = 1, n
               call out%set([i, j, k], square(i, j))
            end do
         end do
         inc = another_incidence(inc, err=err)
      end do
   end function random_latin_squares

   function sam_square(m, u, square_a, square_b, err) result(out)
      integer, intent(in) :: m, u
      integer(ik), intent(in), optional :: square_a(:, :), square_b(:, :)
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :), a(:, :), b(:, :), c(:, :), s(:, :), t(:, :), dmat(:, :)
      integer, allocatable :: jc(:), jd(:), js(:), jt(:)
      integer :: i, j, r, q, pos

      if (m < 2 .or. u < 1 .or. u >= m) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "SAM requires m>=2 and 1<=u<m")
         allocate(out(0, 0))
         return
      end if
      if (present(square_a)) then
         allocate(a(m, m), source=square_a)
      else
         allocate(a(m, m))
         a = latin_square(m)
      end if
      if (present(square_b)) then
         allocate(b(m, m), source=square_b)
      else
         b = a
      end if
      if (.not. is_latin_square(a)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "SAM inputs must be Latin squares")
         allocate(out(0, 0))
         return
      end if
      if (.not. is_latin_square(b)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "SAM inputs must be Latin squares")
         allocate(out(0, 0))
         return
      end if
      allocate(jc(u), jd(u), js(u + 1), jt(u + 1))
      if (modulo(u, 2) == 1) then
         pos = 1
         jc(pos) = 0
         jd(pos) = 1
         js(1:2) = [2, 4]
         jt(1:2) = [3, 5]
         pos = 2
         do q = 0, (u - 3) / 2
            jc(pos) = 6 + 8 * q
            jd(pos) = 7 + 8 * q
            js(pos + 1) = 8 + 8 * q
            jt(pos + 1) = 9 + 8 * q
            pos = pos + 1
         end do
         do q = 0, (u - 3) / 2
            jc(pos) = 13 + 8 * q
            jd(pos) = 12 + 8 * q
            js(pos + 1) = 11 + 8 * q
            jt(pos + 1) = 10 + 8 * q
            pos = pos + 1
         end do
      else
         jc(1:2) = [2, 3]
         jd(1:2) = [0, 4]
         js(1:3) = [1, 7, 9]
         jt(1:3) = [5, 6, 8]
         pos = 3
         do q = 0, (u - 4) / 2
            jc(pos) = 10 + 8 * q
            jd(pos) = 11 + 8 * q
            js(pos + 1) = 12 + 8 * q
            jt(pos + 1) = 13 + 8 * q
            pos = pos + 1
         end do
         do q = 0, (u - 4) / 2
            jc(pos) = 17 + 8 * q
            jd(pos) = 16 + 8 * q
            js(pos + 1) = 15 + 8 * q
            jt(pos + 1) = 14 + 8 * q
            pos = pos + 1
         end do
      end if
      allocate(c(m, m), s(m, m), t(m, m), dmat(m, m), source=0_ik)
      do r = 1, u
         do j = 1, m
            do i = 1, m
               if (b(i, j) == int(r, ik)) s(i, j) = int(i + m * js(r), ik)
               if (a(i, j) == int(r, ik)) c(i, j) = int((m + 1) - i + m * jc(r), ik)
               if (a(i, j) == int(r, ik)) t(i, j) = int(i + m * jt(r), ik)
               if (b(i, j) == int(r, ik)) dmat(i, j) = int((m + 1) - i + m * jd(r), ik)
            end do
         end do
      end do
      do j = 1, m
         do i = 1, m
            if (b(i, j) == int(u + 1, ik)) s(i, j) = int(i + m * js(u + 1), ik)
            if (a(i, j) == int(u + 1, ik)) t(i, j) = int(i + m * jt(u + 1), ik)
         end do
      end do
      allocate(out(2 * m, 2 * m), source=0_ik)
      out(1:m, 1:m) = c
      out(1:m, m + 1:2 * m) = s
      out(m + 1:2 * m, 1:m) = t
      out(m + 1:2 * m, m + 1:2 * m) = dmat
   end function sam_square

end module magic_combinatorics
