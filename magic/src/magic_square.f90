module magic_square
   use magic_kinds, only : ik, dp
   use magic_status, only : magic_error, set_error, MAGIC_INVALID_ARGUMENT
   implicit none
   private

   type, public :: square_sums
      integer(ik), allocatable :: rows(:)
      integer(ik), allocatable :: columns(:)
      integer(ik), allocatable :: major_diagonals(:)
      integer(ik), allocatable :: minor_diagonals(:)
   end type square_sums

   public :: process_index, shift_vector, recurse_permutation
   public :: circulant, latin_square, diag_off, all_square_sums
   public :: magic_constant, magic_square_of_order, magic_2np1, magic_4n
   public :: magic_4np2, strachey_square, lozenge_square, hudson_square
   public :: magic_prime, magic_product_fast, panmagic_4n, panmagic_6npm1
   public :: panmagic_6np1, panmagic_6nm1, transform_square
   public :: reverse_square, rotate_square, standardize_square
   public :: is_semimagic, is_magic, is_panmagic, is_multiplicative_magic, is_normal_square
   public :: is_associative, is_centrosymmetric, is_persymmetric
   public :: is_circulant_square, is_latin_square, is_mostperfect
   public :: is_antimagic, is_totally_antimagic, is_heterosquare
   public :: is_totally_heterosquare, is_sparse_square, is_sam, is_stam
   public :: lex_equal, lex_less, lex_greater

contains

   elemental integer function process_index(x, n) result(y)
      integer, intent(in) :: x, n
      if (n <= 0) then
         y = 0
      else
         y = modulo(x - 1, n) + 1
      end if
   end function process_index

   function shift_vector(x, amount) result(out)
      integer(ik), intent(in) :: x(:)
      integer, intent(in), optional :: amount
      integer(ik), allocatable :: out(:)
      integer :: n, k, i
      n = size(x)
      allocate(out(n))
      if (n == 0) return
      k = 1
      if (present(amount)) k = amount
      k = modulo(k, n)
      do i = 1, n
         out(i) = x(modulo(i - 1 - k, n) + 1)
      end do
   end function shift_vector

   function invert_permutation(perm) result(inv)
      integer, intent(in) :: perm(:)
      integer, allocatable :: inv(:)
      integer :: i
      allocate(inv(size(perm)))
      do i = 1, size(perm)
         inv(perm(i)) = i
      end do
   end function invert_permutation

   function recurse_permutation(perm, power, start, err) result(out)
      integer, intent(in) :: perm(:)
      integer, intent(in) :: power
      integer, intent(in), optional :: start(:)
      type(magic_error), intent(inout), optional :: err
      integer, allocatable :: out(:)
      integer, allocatable :: current(:), use_perm(:), source(:)
      integer :: i, p

      if (any(perm < 1) .or. any(perm > size(perm))) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "invalid permutation")
         allocate(out(0))
         return
      end if
      allocate(source(size(perm)))
      source = [(i, i=1,size(perm))]
      if (present(start)) then
         if (size(start) /= size(perm)) then
            call set_error(err, MAGIC_INVALID_ARGUMENT, "start vector size mismatch")
            allocate(out(0))
            return
         end if
         source = start
      end if
      if (power < 0) then
         use_perm = invert_permutation(perm)
         p = -power
      else
         allocate(use_perm(size(perm)))
         use_perm = perm
         p = power
      end if
      allocate(current(size(perm)))
      current = [(i, i=1,size(perm))]
      do i = 1, p
         current = use_perm(current)
      end do
      allocate(out(size(perm)))
      out = source(current)
   end function recurse_permutation

   function circulant(values) result(a)
      integer(ik), intent(in) :: values(:)
      integer(ik), allocatable :: a(:, :)
      integer :: i, j, n
      n = size(values)
      allocate(a(n, n))
      do j = 1, n
         do i = 1, n
            a(i, j) = values(process_index(1 - i + j, n))
         end do
      end do
   end function circulant

   function latin_square(n) result(a)
      integer, intent(in) :: n
      integer(ik), allocatable :: a(:, :), values(:)
      integer :: i
      if (n < 0) then
         allocate(a(0, 0))
         return
      end if
      allocate(values(n))
      values = [(int(i, ik), i=1,n)]
      a = circulant(values)
   end function latin_square

   function diag_off(a, offset, northwest_southeast, err) result(values)
      integer(ik), intent(in) :: a(:, :)
      integer, intent(in), optional :: offset
      logical, intent(in), optional :: northwest_southeast
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: values(:)
      integer :: n, k, i, j
      logical :: major
      if (size(a, 1) /= size(a, 2)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "diag_off requires a square matrix")
         allocate(values(0))
         return
      end if
      n = size(a, 1)
      k = 0
      if (present(offset)) k = offset
      major = .true.
      if (present(northwest_southeast)) major = northwest_southeast
      allocate(values(n))
      do i = 1, n
         if (major) then
            j = process_index(i + k, n)
         else
            j = process_index(n - i + 1 + k, n)
         end if
         values(i) = a(i, j)
      end do
   end function diag_off

   function all_square_sums(a, err) result(sums)
      integer(ik), intent(in) :: a(:, :)
      type(magic_error), intent(inout), optional :: err
      type(square_sums) :: sums
      integer :: n, i
      if (size(a, 1) /= size(a, 2)) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "all_square_sums requires a square matrix")
         return
      end if
      n = size(a, 1)
      allocate(sums%rows(n), sums%columns(n), sums%major_diagonals(n), sums%minor_diagonals(n))
      sums%rows = sum(a, dim=2)
      sums%columns = sum(a, dim=1)
      do i = 0, n - 1
         sums%major_diagonals(i + 1) = sum(diag_off(a, i, .true.))
         sums%minor_diagonals(i + 1) = sum(diag_off(a, i, .false.))
      end do
   end function all_square_sums

   pure integer(ik) function magic_constant(n, dimension, start) result(value)
      integer, intent(in) :: n
      integer, intent(in), optional :: dimension
      integer(ik), intent(in), optional :: start
      integer :: d
      integer(ik) :: first, count_values
      d = 2
      if (present(dimension)) d = dimension
      first = 1_ik
      if (present(start)) first = start
      count_values = int(n, ik) ** d
      value = int(n, ik) * (2_ik * first + count_values - 1_ik) / 2_ik
   end function magic_constant

   function magic_2np1(m, err) result(a)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      integer :: n, i, j, next_i, next_j, number
      if (m < 0) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "m must be nonnegative")
         allocate(a(0, 0))
         return
      end if
      n = 2 * m + 1
      allocate(a(n, n), source=0_ik)
      i = 1
      j = (n + 1) / 2
      do number = 1, n * n
         a(i, j) = int(number, ik)
         next_i = process_index(i - 1, n)
         next_j = process_index(j + 1, n)
         if (a(next_i, next_j) /= 0_ik) then
            i = process_index(i + 1, n)
         else
            i = next_i
            j = next_j
         end if
      end do
   end function magic_2np1

   function magic_4n(m, err) result(a)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      integer :: n, i, j
      integer(ik) :: value, complement
      if (m < 1) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "m must be positive")
         allocate(a(0, 0))
         return
      end if
      n = 4 * m
      allocate(a(n, n))
      complement = int(n * n + 1, ik)
      do i = 1, n
         do j = 1, n
            value = int((i - 1) * n + j, ik)
            if (modulo(i - 1, 4) == modulo(j - 1, 4) .or. &
                modulo(i - 1, 4) + modulo(j - 1, 4) == 3) then
               a(i, j) = complement - value
            else
               a(i, j) = value
            end if
         end do
      end do
   end function magic_4n

   function strachey_square(m, base_square, err) result(out)
      integer, intent(in) :: m
      integer(ik), intent(in), optional :: base_square(:, :)
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :), base(:, :)
      integer :: p, n, k, i, j
      integer(ik) :: temp, p2

      if (m < 0) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "m must be nonnegative")
         allocate(out(0, 0))
         return
      end if
      p = 2 * m + 1
      n = 2 * p
      if (present(base_square)) then
         if (size(base_square, 1) /= p .or. size(base_square, 2) /= p) then
            call set_error(err, MAGIC_INVALID_ARGUMENT, "base square has wrong order")
            allocate(out(0, 0))
            return
         end if
         allocate(base(p, p), source=base_square)
      else
         allocate(base(p, p))
         base = magic_2np1(m)
      end if
      p2 = int(p * p, ik)
      allocate(out(n, n))
      out(1:p, 1:p) = base
      out(1:p, p + 1:n) = base + 2_ik * p2
      out(p + 1:n, 1:p) = base + 3_ik * p2
      out(p + 1:n, p + 1:n) = base + p2
      k = (p - 1) / 2
      do j = 1, k
         do i = 1, p
            temp = out(i, j)
            out(i, j) = out(i + p, j)
            out(i + p, j) = temp
         end do
      end do
      temp = out(k + 1, 1)
      out(k + 1, 1) = out(k + 1 + p, 1)
      out(k + 1 + p, 1) = temp
      j = k + 1
      temp = out(k + 1, j)
      out(k + 1, j) = out(k + 1 + p, j)
      out(k + 1 + p, j) = temp
      do j = p - k + 2, p
         do i = 1, p
            temp = out(i, j + p)
            out(i, j + p) = out(i + p, j + p)
            out(i + p, j + p) = temp
         end do
      end do
   end function strachey_square

   function magic_4np2(m, err) result(a)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      a = strachey_square(m, err=err)
   end function magic_4np2

   function magic_square_of_order(n, err) result(a)
      integer, intent(in) :: n
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      if (n < 1 .or. n == 2) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "normal magic squares require n=1 or n>=3")
         allocate(a(0, 0))
      else if (modulo(n, 2) == 1) then
         a = standardize_square(magic_2np1((n - 1) / 2))
      else if (modulo(n, 4) == 0) then
         a = standardize_square(magic_4n(n / 4))
      else
         a = standardize_square(magic_4np2((n - 2) / 4))
      end if
   end function magic_square_of_order

   function lozenge_square(m, err) result(a)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      ! The lozenge construction is equivalent to an odd-order normal square.
      a = magic_2np1(m, err)
   end function lozenge_square

   function hudson_square(n, start_a, start_b, err) result(out)
      integer, intent(in) :: n
      integer, intent(in), optional :: start_a(:), start_b(:)
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :)
      integer, allocatable :: a(:), b(:), perm(:), aa(:), bb(:)
      integer :: i, j
      if (n < 3) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "Hudson construction requires n>=3")
         allocate(out(0, 0))
         return
      end if
      allocate(a(n), b(n), perm(n))
      a(1) = n - 1
      do i = 2, n
         a(i) = i - 2
      end do
      do i = 1, n - 2
         b(i) = i + 1
      end do
      b(n - 1) = n
      b(n) = 1
      if (present(start_a)) then
         if (size(start_a) /= n) then
            call set_error(err, MAGIC_INVALID_ARGUMENT, "start_a size mismatch")
            allocate(out(0, 0))
            return
         end if
         a = start_a
      end if
      if (present(start_b)) then
         if (size(start_b) /= n) then
            call set_error(err, MAGIC_INVALID_ARGUMENT, "start_b size mismatch")
            allocate(out(0, 0))
            return
         end if
         b = start_b
      end if
      perm(1:2) = [n - 1, n]
      perm(3:n) = [(i, i=1,n-2)]
      allocate(out(n, n))
      do i = 0, n - 1
         aa = recurse_permutation(perm, i, a)
         bb = recurse_permutation(perm, -i, b)
         do j = 1, n
            out(i + 1, j) = int(n * aa(j) + bb(j), ik)
         end do
      end do
   end function hudson_square

   function magic_prime(n, multiplier_i, multiplier_j, err) result(a)
      integer, intent(in) :: n
      integer, intent(in), optional :: multiplier_i, multiplier_j
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: a(:, :)
      integer :: i, j, p, q
      if (n < 1) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "n must be positive")
         allocate(a(0, 0))
         return
      end if
      p = 2
      q = 3
      if (present(multiplier_i)) p = multiplier_i
      if (present(multiplier_j)) q = multiplier_j
      allocate(a(n, n))
      do j = 1, n
         do i = 1, n
            a(i, j) = int(n * modulo(j - p * i + p - 1, n) + &
                              modulo(j - q * i + q - 1, n) + 1, ik)
         end do
      end do
   end function magic_prime

   function reverse_square(a, reverse_rows, reverse_columns) result(out)
      integer(ik), intent(in) :: a(:, :)
      logical, intent(in), optional :: reverse_rows, reverse_columns
      integer(ik), allocatable :: out(:, :)
      logical :: rr, rc
      integer :: i, j, ii, jj
      rr = .true.
      rc = .true.
      if (present(reverse_rows)) rr = reverse_rows
      if (present(reverse_columns)) rc = reverse_columns
      allocate(out(size(a, 1), size(a, 2)))
      do j = 1, size(a, 2)
         jj = merge(size(a, 2) - j + 1, j, rc)
         do i = 1, size(a, 1)
            ii = merge(size(a, 1) - i + 1, i, rr)
            out(i, j) = a(ii, jj)
         end do
      end do
   end function reverse_square

   function rotate_square(a, rights) result(out)
      integer(ik), intent(in) :: a(:, :)
      integer, intent(in), optional :: rights
      integer(ik), allocatable :: out(:, :)
      integer :: r, i, j
      r = 1
      if (present(rights)) r = modulo(rights, 4)
      select case (r)
      case (0)
         allocate(out(size(a, 1), size(a, 2)), source=a)
      case (1)
         allocate(out(size(a, 2), size(a, 1)))
         do j = 1, size(out, 2)
            do i = 1, size(out, 1)
               out(i, j) = a(size(a, 1) - j + 1, i)
            end do
         end do
      case (2)
         out = reverse_square(a)
      case (3)
         allocate(out(size(a, 2), size(a, 1)))
         do j = 1, size(out, 2)
            do i = 1, size(out, 1)
               out(i, j) = a(j, size(a, 2) - i + 1)
            end do
         end do
      end select
   end function rotate_square

   function transform_square(a, transformation) result(out)
      integer(ik), intent(in) :: a(:, :)
      integer, intent(in) :: transformation
      integer(ik), allocatable :: out(:, :), tmp(:, :)
      integer :: code
      code = modulo(transformation, 8)
      allocate(tmp(size(a, 1), size(a, 2)), source=a)
      if (modulo(code, 2) == 1) tmp = transpose(tmp)
      if (modulo(code / 2, 2) == 1) tmp = reverse_square(tmp, .true., .false.)
      if (modulo(code / 4, 2) == 1) tmp = reverse_square(tmp, .false., .true.)
      out = tmp
   end function transform_square

   pure logical function lex_equal(a, b) result(answer)
      integer(ik), intent(in) :: a(:, :), b(:, :)
      answer = all(shape(a) == shape(b))
      if (answer) answer = all(a == b)
   end function lex_equal

   pure logical function lex_less(a, b) result(answer)
      integer(ik), intent(in) :: a(:, :), b(:, :)
      integer(ik) :: av(size(a)), bv(size(b))
      integer :: i
      answer = .false.
      if (any(shape(a) /= shape(b))) return
      av = reshape(a, [size(a)])
      bv = reshape(b, [size(b)])
      do i = 1, size(a)
         if (av(i) < bv(i)) then
            answer = .true.
            return
         else if (av(i) > bv(i)) then
            return
         end if
      end do
   end function lex_less

   pure logical function lex_greater(a, b) result(answer)
      integer(ik), intent(in) :: a(:, :), b(:, :)
      answer = .not. lex_equal(a, b) .and. .not. lex_less(a, b)
   end function lex_greater

   recursive function standardize_square(a, one_minus) result(out)
      integer(ik), intent(in) :: a(:, :)
      logical, intent(in), optional :: one_minus
      integer(ik), allocatable :: out(:, :), candidate(:, :), complement(:, :)
      integer :: k
      logical :: compare_complement
      out = transform_square(a, 0)
      do k = 1, 7
         candidate = transform_square(a, k)
         if (lex_less(candidate, out)) out = candidate
      end do
      compare_complement = .false.
      if (present(one_minus)) compare_complement = one_minus
      if (compare_complement) then
         allocate(complement(size(a, 1), size(a, 2)))
         complement = maxval(a) + 1_ik - a
         candidate = standardize_square(complement, .false.)
         if (lex_less(candidate, out)) out = candidate
      end if
   end function standardize_square

   pure logical function all_equal_integer(x) result(answer)
      integer(ik), intent(in) :: x(:)
      if (size(x) <= 1) then
         answer = .true.
      else
         answer = all(x == x(1))
      end if
   end function all_equal_integer

   function is_semimagic(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      type(square_sums) :: sums
      if (size(a, 1) /= size(a, 2)) then
         answer = .false.
         return
      end if
      sums = all_square_sums(a)
      answer = all_equal_integer([sums%rows, sums%columns])
   end function is_semimagic

   function is_magic(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      type(square_sums) :: sums
      if (size(a, 1) /= size(a, 2)) then
         answer = .false.
         return
      end if
      sums = all_square_sums(a)
      answer = all_equal_integer([sums%rows, sums%columns, &
                                  sums%major_diagonals(1:1), sums%minor_diagonals(1:1)])
   end function is_magic

   function is_panmagic(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      type(square_sums) :: sums
      if (size(a, 1) /= size(a, 2)) then
         answer = .false.
         return
      end if
      sums = all_square_sums(a)
      answer = all_equal_integer([sums%rows, sums%columns, sums%major_diagonals, sums%minor_diagonals])
   end function is_panmagic

   function is_multiplicative_magic(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      integer :: n, i, j
      integer(ik) :: target, value
      if (size(a, 1) /= size(a, 2) .or. size(a, 1) == 0) then
         answer = .false.
         return
      end if
      n = size(a, 1)
      target = product(a(1, :))
      answer = .true.
      do i = 1, n
         if (product(a(i, :)) /= target .or. product(a(:, i)) /= target) then
            answer = .false.
            return
         end if
      end do
      value = 1_ik
      do i = 1, n
         value = value * a(i, i)
      end do
      if (value /= target) answer = .false.
      value = 1_ik
      do i = 1, n
         j = n - i + 1
         value = value * a(i, j)
      end do
      if (value /= target) answer = .false.
   end function is_multiplicative_magic

   function is_normal_square(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      integer(ik), allocatable :: sorted(:)
      integer :: i
      if (size(a, 1) /= size(a, 2)) then
         answer = .false.
         return
      end if
      allocate(sorted(size(a)))
      sorted = reshape(a, [size(a)])
      call sort_integer(sorted)
      answer = sorted(1) == 1_ik
      do i = 2, size(sorted)
         answer = answer .and. sorted(i) == sorted(i - 1) + 1_ik
      end do
   end function is_normal_square

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

   function is_associative(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      integer(ik), allocatable :: paired(:, :)
      if (.not. is_magic(a)) then
         answer = .false.
         return
      end if
      paired = a + reverse_square(a)
      answer = all(paired == paired(1, 1))
   end function is_associative

   pure logical function is_centrosymmetric(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      integer :: i, j
      answer = .true.
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            if (a(i, j) /= a(size(a, 1) - i + 1, size(a, 2) - j + 1)) then
               answer = .false.
               return
            end if
         end do
      end do
   end function is_centrosymmetric

   pure logical function is_persymmetric(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      integer :: i, j, n
      if (size(a, 1) /= size(a, 2)) then
         answer = .false.
         return
      end if
      n = size(a, 1)
      answer = .true.
      do j = 1, n
         do i = 1, n
            if (a(i, n - j + 1) /= a(j, n - i + 1)) then
               answer = .false.
               return
            end if
         end do
      end do
   end function is_persymmetric

   function is_circulant_square(a, row_shift, column_shift) result(answer)
      integer(ik), intent(in) :: a(:, :)
      integer, intent(in), optional :: row_shift, column_shift
      logical :: answer
      integer :: rs, cs, i, j, ii, jj
      rs = 1
      cs = 1
      if (present(row_shift)) rs = row_shift
      if (present(column_shift)) cs = column_shift
      answer = .true.
      do j = 1, size(a, 2)
         jj = process_index(j - cs, size(a, 2))
         do i = 1, size(a, 1)
            ii = process_index(i - rs, size(a, 1))
            if (a(i, j) /= a(ii, jj)) then
               answer = .false.
               return
            end if
         end do
      end do
   end function is_circulant_square

   function is_latin_square(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      integer(ik), allocatable :: work(:)
      integer :: n, i, j
      if (size(a, 1) /= size(a, 2)) then
         answer = .false.
         return
      end if
      n = size(a, 1)
      allocate(work(n))
      answer = .true.
      do i = 1, n
         work = a(i, :)
         call sort_integer(work)
         answer = answer .and. all(work == [(int(j, ik), j=1,n)])
         work = a(:, i)
         call sort_integer(work)
         answer = answer .and. all(work == [(int(j, ik), j=1,n)])
         if (.not. answer) return
      end do
   end function is_latin_square

   function is_mostperfect(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      integer :: n, i, j, half
      integer(ik) :: two_sum, pair_sum
      if (size(a, 1) /= size(a, 2) .or. modulo(size(a, 1), 4) /= 0) then
         answer = .false.
         return
      end if
      n = size(a, 1)
      half = n / 2
      two_sum = a(1, 1) + a(1, 2) + a(2, 1) + a(2, 2)
      pair_sum = a(1, 1) + a(process_index(1 + half, n), process_index(1 + half, n))
      answer = .true.
      do j = 1, n
         do i = 1, n
            if (a(i, j) + a(i, process_index(j + 1, n)) + &
                a(process_index(i + 1, n), j) + &
                a(process_index(i + 1, n), process_index(j + 1, n)) /= two_sum) then
               answer = .false.
               return
            end if
            if (a(i, j) + a(process_index(i + half, n), process_index(j + half, n)) /= pair_sum) then
               answer = .false.
               return
            end if
         end do
      end do
   end function is_mostperfect

   function distinct_consecutive(x) result(answer)
      integer(ik), intent(in) :: x(:)
      logical :: answer
      integer(ik), allocatable :: work(:)
      integer :: i
      allocate(work(size(x)), source=x)
      call sort_integer(work)
      answer = .true.
      do i = 2, size(work)
         if (work(i) /= work(i - 1) + 1_ik) then
            answer = .false.
            return
         end if
      end do
   end function distinct_consecutive

   function all_distinct(x) result(answer)
      integer(ik), intent(in) :: x(:)
      logical :: answer
      integer(ik), allocatable :: work(:)
      integer :: i
      allocate(work(size(x)), source=x)
      call sort_integer(work)
      answer = .true.
      do i = 2, size(work)
         if (work(i) == work(i - 1)) then
            answer = .false.
            return
         end if
      end do
   end function all_distinct

   function is_antimagic(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      type(square_sums) :: sums
      sums = all_square_sums(a)
      answer = distinct_consecutive([sums%rows, sums%columns])
   end function is_antimagic

   function is_totally_antimagic(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      type(square_sums) :: sums
      sums = all_square_sums(a)
      answer = distinct_consecutive([sums%rows, sums%columns, &
                                     sums%major_diagonals(1:1), sums%minor_diagonals(1:1)])
   end function is_totally_antimagic

   function is_heterosquare(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      type(square_sums) :: sums
      sums = all_square_sums(a)
      answer = all_distinct([sums%rows, sums%columns])
   end function is_heterosquare

   function is_totally_heterosquare(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      type(square_sums) :: sums
      sums = all_square_sums(a)
      answer = all_distinct([sums%rows, sums%columns, &
                             sums%major_diagonals(1:1), sums%minor_diagonals(1:1)])
   end function is_totally_heterosquare

   function is_sparse_square(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      integer(ik), allocatable :: nonzero(:)
      integer :: i, n
      n = count(a /= 0_ik)
      if (n == 0) then
         answer = .false.
         return
      end if
      allocate(nonzero(n))
      nonzero = pack(a, a /= 0_ik)
      call sort_integer(nonzero)
      answer = nonzero(1) == 1_ik
      do i = 2, n
         answer = answer .and. nonzero(i) == nonzero(i - 1) + 1_ik
      end do
   end function is_sparse_square

   function is_sam(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      answer = is_antimagic(a)
      if (answer) answer = is_sparse_square(a)
   end function is_sam

   function is_stam(a) result(answer)
      integer(ik), intent(in) :: a(:, :)
      logical :: answer
      answer = is_totally_antimagic(a)
      if (answer) answer = is_sparse_square(a)
   end function is_stam

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

   function magic_product_fast(a, b) result(out)
      integer(ik), intent(in) :: a(:, :), b(:, :)
      integer(ik), allocatable :: out(:, :)
      integer(ik), allocatable :: ones_a(:, :), ones_b(:, :)
      integer :: nb
      nb = size(b, 1)
      allocate(ones_a(size(a, 1), size(a, 2)), source=1_ik)
      allocate(ones_b(size(b, 1), size(b, 2)), source=1_ik)
      out = int(nb * nb, ik) * (kronecker_integer(a, ones_b) - 1_ik) + &
            kronecker_integer(ones_a, b)
   end function magic_product_fast

   function panmagic_6npm1(n, err) result(out)
      integer, intent(in) :: n
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :)
      integer, allocatable :: jj(:, :)
      integer :: i, j
      if (n < 5 .or. modulo(n, 6) == 3) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "n must be 6m-1 or 6m+1")
         allocate(out(0, 0))
         return
      end if
      allocate(jj(n, n), out(n, n))
      do j = 1, n
         do i = 1, n
            jj(i, j) = process_index((j - 1) * (n - 2) + i, n)
         end do
      end do
      do j = 1, n
         do i = 1, n
            out(i, j) = int(jj(i, j) + n * jj(j, i) - n, ik)
         end do
      end do
   end function panmagic_6npm1

   function panmagic_6np1(m, err) result(out)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :)
      out = panmagic_6npm1(6 * m + 1, err)
   end function panmagic_6np1

   function panmagic_6nm1(m, err) result(out)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :)
      out = panmagic_6npm1(6 * m - 1, err)
   end function panmagic_6nm1

   function panmagic_4n(m, err) result(out)
      integer, intent(in) :: m
      type(magic_error), intent(inout), optional :: err
      integer(ik), allocatable :: out(:, :), half(:, :), joined(:, :), rotated(:, :)
      integer :: n, i, j, block
      if (m < 1) then
         call set_error(err, MAGIC_INVALID_ARGUMENT, "m must be positive")
         allocate(out(0, 0))
         return
      end if
      n = 4 * m
      allocate(half(n, 2 * m))
      do block = 1, 2 * m
         do j = 1, 2 * m
            half(2 * block - 1, j) = int(j, ik)
            half(2 * block, j) = int(n - j + 1, ik)
         end do
      end do
      allocate(joined(n, n))
      joined(:, 1:2 * m) = half
      do i = 1, n
         joined(i, 2 * m + 1:n) = half(process_index(i - 1, n), :)
      end do
      rotated = rotate_square(joined)
      out = joined + int(n, ik) * (rotated - 1_ik)
   end function panmagic_4n

end module magic_square
