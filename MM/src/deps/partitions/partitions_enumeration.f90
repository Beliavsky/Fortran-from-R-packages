! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

module partitions_enumeration
    use partitions_kinds, only : i8
    use partitions_counts, only : partition_count, distinct_partition_count, restricted_partition_count, &
                                  block_partition_count
    implicit none
    private

    public :: triangular_index
    public :: first_part, next_part, is_last_part, parts
    public :: first_distinct_part, next_distinct_part, is_last_distinct_part, distinct_parts
    public :: first_restricted_part, next_restricted_part, is_last_restricted_part, restricted_parts
    public :: first_block_part, next_block_part, is_last_block_part, block_parts
    public :: first_composition, next_composition, is_last_composition, compositions
    public :: to_binary, to_decimal, composition_to_binary, binary_to_composition
    public :: conjugate, durfee

contains

    pure integer function triangular_index(n) result(value)
        integer, intent(in) :: n
        real :: x

        if (n < 0) then
            value = -1
            return
        end if
        x = real(n)
        value = int(floor((sqrt(1.0 + 8.0 * x) - 1.0) / 2.0))
    end function triangular_index

    function first_part(n) result(part)
        integer, intent(in) :: n
        integer, allocatable :: part(:)

        if (n < 1) error stop "first_part: n must be positive"
        allocate(part(n))
        part = 0
        part(1) = n
    end function first_part

    pure logical function is_last_part(part) result(last)
        integer, intent(in) :: part(:)
        last = size(part) > 0 .and. all(part == 1)
    end function is_last_part

    subroutine next_part(part)
        integer, intent(inout) :: part(:)
        integer :: a, b, yy, nleft, i

        if (size(part) < 2) error stop "next_part: no successor"
        if (is_last_part(part)) error stop "next_part: already at final partition"

        a = 1
        do while (a < size(part))
            if (part(a + 1) <= 0) exit
            a = a + 1
        end do

        b = a
        do while (b >= 1)
            if (part(b) /= 1) exit
            b = b - 1
        end do
        if (b < 1) error stop "next_part: invalid partition"

        if (part(a) > 1) then
            part(a) = part(a) - 1
            part(a + 1) = 1
            return
        end if

        nleft = a - b
        part(b) = part(b) - 1
        yy = part(b)
        nleft = nleft + 1
        i = b
        do while (nleft >= yy)
            i = i + 1
            part(i) = yy
            nleft = nleft - yy
        end do
        if (nleft > 0) then
            i = i + 1
            part(i) = nleft
        end if
        if (i < a) part(i + 1:a) = 0
    end subroutine next_part

    function parts(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        integer, allocatable :: part(:)
        integer(i8) :: nc8
        integer :: nc, j

        if (n < 0) error stop "parts: n must be nonnegative"
        if (n == 0) then
            allocate(out(0,0))
            return
        end if
        nc8 = partition_count(n)
        call count_to_default(nc8, nc, "parts")
        allocate(out(n, nc))
        part = first_part(n)
        out(:, 1) = part
        do j = 2, nc
            call next_part(part)
            out(:, j) = part
        end do
    end function parts

    function first_distinct_part(n) result(part)
        integer, intent(in) :: n
        integer, allocatable :: part(:)
        integer :: nr

        if (n < 1) error stop "first_distinct_part: n must be positive"
        nr = max(1, triangular_index(n))
        allocate(part(nr))
        part = 0
        part(1) = n
    end function first_distinct_part

    pure logical function is_last_distinct_part(part) result(last)
        integer, intent(in) :: part(:)
        integer, allocatable :: revp(:), d(:)
        integer :: i

        if (size(part) <= 1) then
            last = .true.
            return
        end if
        allocate(revp(size(part)), d(size(part) - 1))
        revp = part(size(part):1:-1)
        do i = 1, size(d)
            d(i) = revp(i + 1) - revp(i)
        end do
        last = all(d == 1 .or. d == 2) .and. count(d == 2) == 1
    end function is_last_distinct_part

    subroutine next_distinct_part(part)
        integer, intent(inout) :: part(:)
        integer :: a, aa, d, nleft, yy

        if (is_last_distinct_part(part)) error stop "next_distinct_part: already at final partition"

        a = size(part)
        do while (a >= 1)
            if (part(a) /= 0) exit
            a = a - 1
        end do
        if (a < 1) error stop "next_distinct_part: invalid partition"
        aa = a

        d = 1
        nleft = 0
        do while (part(a) - d < 2)
            nleft = nleft + part(a)
            a = a - 1
            d = d + 1
            if (a < 1) error stop "next_distinct_part: invalid state"
        end do

        part(a) = part(a) - 1
        yy = part(a)
        a = a + 1
        nleft = nleft + 1
        do while (nleft >= yy)
            yy = yy - 1
            part(a) = yy
            a = a + 1
            nleft = nleft - yy
        end do
        part(a) = nleft
        if (a < aa) part(a + 1:aa) = 0
    end subroutine next_distinct_part

    function distinct_parts(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        integer, allocatable :: part(:)
        integer(i8) :: nc8
        integer :: nc, j, nr

        if (n < 1) then
            allocate(out(0,0))
            return
        end if
        nc8 = distinct_partition_count(n)
        call count_to_default(nc8, nc, "distinct_parts")
        nr = max(1, triangular_index(n))
        allocate(out(nr, nc))
        part = first_distinct_part(n)
        out(:,1) = part
        do j = 2, nc
            call next_distinct_part(part)
            out(:,j) = part
        end do
    end function distinct_parts

    function first_restricted_part(n, m, include_zero) result(part)
        integer, intent(in) :: n, m
        logical, intent(in), optional :: include_zero
        integer, allocatable :: part(:)
        logical :: inc0

        if (n < 0 .or. m < 1) error stop "first_restricted_part: invalid n or m"
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        if (.not. inc0 .and. m > n) error stop "first_restricted_part: m > n without zeros"

        allocate(part(m))
        part = 0
        if (inc0) then
            part(1) = n
        else
            part = 1
            part(1) = n - m + 1
        end if
    end function first_restricted_part

    pure logical function is_last_restricted_part(part) result(last)
        integer, intent(in) :: part(:)
        if (size(part) == 0) then
            last = .true.
        else
            last = maxval(part) - minval(part) <= 1
        end if
    end function is_last_restricted_part

    subroutine next_restricted_part(part)
        integer, intent(inout) :: part(:)
        integer, allocatable :: x(:)
        logical :: done

        if (is_last_restricted_part(part)) error stop "next_restricted_part: already at final partition"
        x = part(size(part):1:-1)
        call next_restricted_ascending(x, done)
        if (done) error stop "next_restricted_part: no successor"
        part = x(size(x):1:-1)
    end subroutine next_restricted_part

    subroutine next_restricted_ascending(x, done)
        integer, intent(inout) :: x(:)
        logical, intent(out) :: done
        integer :: a, j, m, r

        done = .false.
        a = size(x)
        m = x(a)
        do while (m - x(a) < 2)
            a = a - 1
            if (a < 1) then
                done = .true.
                return
            end if
        end do

        x(a) = x(a) + 1
        j = x(a)
        r = -1
        do
            r = r + x(a) - j
            x(a) = j
            a = a + 1
            if (a >= size(x)) exit
        end do
        x(size(x)) = x(size(x)) + r
    end subroutine next_restricted_ascending

    function restricted_parts(n, m, include_zero, decreasing) result(out)
        integer, intent(in) :: n, m
        logical, intent(in), optional :: include_zero, decreasing
        integer, allocatable :: out(:,:)
        integer, allocatable :: part(:)
        logical :: inc0, decr
        integer(i8) :: nc8
        integer :: nc, j

        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        decr = .true.
        if (present(decreasing)) decr = decreasing
        if (m < 1 .or. n < 0) error stop "restricted_parts: invalid n or m"
        if (.not. inc0 .and. m > n) error stop "restricted_parts: m > n without zeros"

        nc8 = restricted_partition_count(m, n, inc0)
        call count_to_default(nc8, nc, "restricted_parts")
        allocate(out(m, nc))
        if (nc == 0) return
        part = first_restricted_part(n, m, inc0)
        out(:,1) = part
        do j = 2, nc
            call next_restricted_part(part)
            out(:,j) = part
        end do
        if (.not. decr) out = out(m:1:-1,:)
    end function restricted_parts

    function first_block_part(f, n, include_fewer) result(part)
        integer, intent(in) :: f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        integer, allocatable :: part(:)
        logical :: fewer
        integer :: remaining, i

        if (any(f < 0)) error stop "first_block_part: capacities must be nonnegative"
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        allocate(part(size(f)))
        part = 0
        if (.not. present(n) .or. fewer) return
        if (n < 0 .or. n > sum(f)) error stop "first_block_part: invalid n"

        remaining = n
        do i = 1, size(f)
            part(i) = min(f(i), remaining)
            remaining = remaining - part(i)
        end do
    end function first_block_part

    subroutine next_block_part(part, f, n, include_fewer)
        integer, intent(inout) :: part(:)
        integer, intent(in) :: f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        integer, allocatable :: x(:), cap(:)
        logical :: fewer, done
        integer :: total, s

        if (size(part) /= size(f)) error stop "next_block_part: size mismatch"
        if (any(f < 0) .or. any(part < 0) .or. any(part > f)) error stop "next_block_part: invalid vector"
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        s = sum(f)

        if (.not. present(n)) then
            total = s
            allocate(x(size(part) + 1), cap(size(f) + 1))
            x = [s - sum(part), part]
            cap = [s, f]
        else if (fewer) then
            if (n < 0 .or. n > s) error stop "next_block_part: invalid n"
            total = n
            allocate(x(size(part) + 1), cap(size(f) + 1))
            x = [n - sum(part), part]
            cap = [s, f]
        else
            if (sum(part) /= n) error stop "next_block_part: part does not sum to n"
            total = n
            allocate(x(size(part)), cap(size(f)))
            x = part
            cap = f
        end if

        if (sum(x) /= total) error stop "next_block_part: invalid total"
        call next_bounded_fixed(x, cap, done)
        if (done) error stop "next_block_part: already at final block partition"
        if (size(x) == size(part) + 1) then
            part = x(2:)
        else
            part = x
        end if
    end subroutine next_block_part

    logical function is_last_block_part(part, f, n, include_fewer) result(last)
        integer, intent(in) :: part(:), f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        integer, allocatable :: tmp(:)

        tmp = part
        last = .false.
        if (.not. valid_block_state(tmp, f, n, include_fewer)) then
            last = .false.
            return
        end if
        call try_next_block(tmp, f, n, include_fewer, last)
    end function is_last_block_part

    subroutine try_next_block(part, f, n, include_fewer, done)
        integer, intent(inout) :: part(:)
        integer, intent(in) :: f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        logical, intent(out) :: done
        integer, allocatable :: x(:), cap(:)
        logical :: fewer
        integer :: s, total

        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        s = sum(f)
        if (.not. present(n)) then
            total = s
            allocate(x(size(part) + 1), cap(size(f) + 1))
            x(1) = s - sum(part)
            x(2:) = part
            cap(1) = s
            cap(2:) = f
        else if (fewer) then
            total = n
            allocate(x(size(part) + 1), cap(size(f) + 1))
            x(1) = n - sum(part)
            x(2:) = part
            cap(1) = s
            cap(2:) = f
        else
            total = n
            allocate(x(size(part)), cap(size(f)))
            x = part
            cap = f
        end if
        if (sum(x) /= total) then
            done = .true.
            return
        end if
        call next_bounded_fixed(x, cap, done)
        if (.not. done) then
            if (size(x) == size(part) + 1) then
                part = x(2:)
            else
                part = x
            end if
        end if
    end subroutine try_next_block

    pure logical function valid_block_state(part, f, n, include_fewer) result(ok)
        integer, intent(in) :: part(:), f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        logical :: fewer

        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        ok = size(part) == size(f)
        if (.not. ok) return
        ok = all(part >= 0) .and. all(part <= f)
        if (.not. ok) return
        if (present(n)) then
            if (fewer) then
                ok = sum(part) <= n
            else
                ok = sum(part) == n
            end if
        end if
    end function valid_block_state

    subroutine next_bounded_fixed(x, cap, done)
        integer, intent(inout) :: x(:)
        integer, intent(in) :: cap(:)
        logical, intent(out) :: done
        integer :: a, i, j

        if (size(x) /= size(cap)) error stop "next_bounded_fixed: size mismatch"
        if (size(x) == 0) then
            done = .true.
            return
        end if

        i = 1
        a = x(1)
        do
            if (x(i) == 0) then
                i = i + 1
                if (i > size(x)) then
                    done = .true.
                    return
                end if
                a = a + x(i)
            else
                i = i + 1
                if (i > size(x)) then
                    done = .true.
                    return
                end if
                if (x(i) == cap(i)) then
                    a = a + x(i)
                else
                    exit
                end if
            end if
        end do

        a = a - 1
        x(i) = x(i) + 1
        i = i - 1
        x(i) = x(i) - 1
        i = i + 1

        do j = 1, i - 1
            if (a < cap(j)) then
                x(j) = a
                a = 0
            else
                x(j) = cap(j)
                a = a - cap(j)
            end if
        end do
        done = .false.
    end subroutine next_bounded_fixed

    function block_parts(f, n, include_fewer) result(out)
        integer, intent(in) :: f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        integer, allocatable :: out(:,:)
        integer, allocatable :: part(:)
        logical :: fewer
        integer(i8) :: nc8
        integer :: nc, j, target

        if (any(f < 0)) error stop "block_parts: capacities must be nonnegative"
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer

        if (present(n)) then
            target = n
            if (target < 0 .or. target > sum(f)) error stop "block_parts: invalid n"
            nc8 = block_partition_count(f, target, fewer)
        else
            target = sum(f)
            nc8 = block_partition_count(f, target, .true.)
        end if
        call count_to_default(nc8, nc, "block_parts")
        allocate(out(size(f), nc))
        if (nc == 0) return
        part = first_block_part(f, n, fewer)
        out(:,1) = part
        do j = 2, nc
            call next_block_part(part, f, n, fewer)
            out(:,j) = part
        end do
    end function block_parts

    function first_composition(n, m, include_zero) result(comp)
        integer, intent(in) :: n
        integer, intent(in), optional :: m
        logical, intent(in), optional :: include_zero
        integer, allocatable :: comp(:)
        logical :: inc0

        if (n < 1) error stop "first_composition: n must be positive"
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        if (present(m)) then
            if (m < 1) error stop "first_composition: m must be positive"
            allocate(comp(m))
            if (inc0) then
                comp = 0
                comp(1) = n
            else
                if (m > n) error stop "first_composition: m > n without zeros"
                comp = 1
                comp(1) = n - m + 1
            end if
        else
            allocate(comp(n))
            comp = 0
            comp(1) = n
        end if
    end function first_composition

    logical function is_last_composition(comp, restricted, include_zero) result(last)
        integer, intent(in) :: comp(:)
        logical, intent(in) :: restricted
        logical, intent(in), optional :: include_zero
        logical :: inc0
        integer :: m, n

        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        n = sum(comp)
        if (restricted) then
            m = size(comp)
            if (inc0) then
                if (m == 1) then
                    last = comp(1) == n
                else
                    last = all(comp(1:m - 1) == 0) .and. comp(m) == n
                end if
            else
                if (m == 1) then
                    last = comp(1) == n
                else
                    last = all(comp(1:m - 1) == 1) .and. comp(m) == n - m + 1
                end if
            end if
        else
            last = count(comp > 0) == n .and. all(pack(comp, comp > 0) == 1)
        end if
    end function is_last_composition

    function next_composition(comp, restricted, include_zero) result(next)
        integer, intent(in) :: comp(:)
        logical, intent(in) :: restricted
        logical, intent(in), optional :: include_zero
        integer, allocatable :: next(:), bits(:), rbits(:), decoded(:), work(:)
        logical :: inc0
        integer(i8) :: code
        integer :: n, m

        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        if (is_last_composition(comp, restricted, inc0)) error stop "next_composition: already final"
        n = sum(comp)

        if (restricted) then
            m = size(comp)
            allocate(work(m))
            if (inc0) then
                work = comp
                call next_block_part(work, spread(n, 1, m), n)
                next = work
            else
                work = comp - 1
                call next_block_part(work, spread(n - 1, 1, m), n - m)
                next = work + 1
            end if
        else
            bits = composition_to_binary(comp)
            rbits = bits(size(bits):1:-1)
            code = to_decimal(rbits) + 1_i8
            rbits = to_binary(code, n - 1)
            decoded = binary_to_composition(rbits)
            next = decoded(size(decoded):1:-1)
        end if
    end function next_composition

    function compositions(n, m, include_zero) result(out)
        integer, intent(in) :: n
        integer, intent(in), optional :: m
        logical, intent(in), optional :: include_zero
        integer, allocatable :: out(:,:), comp(:), tmp(:)
        logical :: inc0
        integer(i8) :: nc8
        integer :: rows, nc, j

        if (n < 1) error stop "compositions: n must be positive"
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        if (present(m)) then
            rows = m
            if (inc0) then
                nc8 = block_partition_count(spread(n, 1, m), n, .false.)
            else
                if (m > n) error stop "compositions: m > n without zeros"
                nc8 = block_partition_count(spread(n - 1, 1, m), n - m, .false.)
            end if
        else
            rows = n
            if (n - 1 >= bit_size(0_i8) - 1) error stop "compositions: count exceeds int64"
            nc8 = shiftl(1_i8, n - 1)
        end if
        call count_to_default(nc8, nc, "compositions")
        allocate(out(rows, nc))
        comp = first_composition(n, m, inc0)
        out = 0
        out(1:size(comp),1) = comp
        do j = 2, nc
            tmp = next_composition(comp, present(m), inc0)
            if (present(m)) then
                comp = tmp
            else
                deallocate(comp)
                allocate(comp(rows))
                comp = 0
                comp(1:size(tmp)) = tmp
            end if
            out(:,j) = comp
        end do
    end function compositions

    function to_binary(number, len) result(bits)
        integer(i8), intent(in) :: number
        integer, intent(in) :: len
        integer, allocatable :: bits(:)
        integer(i8) :: n
        integer :: i

        if (number < 0_i8 .or. len < 0) error stop "to_binary: invalid argument"
        allocate(bits(len))
        bits = 0
        n = number
        do i = len, 1, -1
            bits(i) = int(iand(n, 1_i8))
            n = shiftr(n, 1)
        end do
        if (n /= 0_i8) error stop "to_binary: number does not fit requested length"
    end function to_binary

    function to_decimal(bits) result(number)
        integer, intent(in) :: bits(:)
        integer(i8) :: number
        integer :: i

        if (any(bits /= 0 .and. bits /= 1)) error stop "to_decimal: bits must be 0 or 1"
        number = 0_i8
        do i = 1, size(bits)
            if (number > shiftr(huge(number), 1)) error stop "to_decimal: integer overflow"
            number = shiftl(number, 1) + int(bits(i), i8)
        end do
    end function to_decimal

    function composition_to_binary(comp) result(bits)
        integer, intent(in) :: comp(:)
        integer, allocatable :: bits(:), pos(:)
        integer :: n, i, j, p, nz

        if (any(comp < 0)) error stop "composition_to_binary: negative component"
        nz = count(comp > 0)
        if (nz == 0) then
            allocate(bits(0))
            return
        end if
        allocate(pos(nz))
        p = 0
        do i = 1, size(comp)
            if (comp(i) > 0) then
                p = p + 1
                pos(p) = comp(i)
            end if
        end do
        n = sum(pos)
        allocate(bits(max(0, n - 1)))
        if (n <= 1) return
        bits = 0
        p = 1
        do i = 1, nz
            do j = 1, pos(i) - 1
                if (p <= n - 1) bits(p) = 0
                p = p + 1
            end do
            if (i < nz .and. p <= n - 1) then
                bits(p) = 1
                p = p + 1
            end if
        end do
    end function composition_to_binary

    function binary_to_composition(bits) result(comp)
        integer, intent(in) :: bits(:)
        integer, allocatable :: comp(:)
        integer :: i, p

        if (any(bits /= 0 .and. bits /= 1)) error stop "binary_to_composition: bits must be 0 or 1"
        allocate(comp(1 + count(bits /= 0)))
        comp = 1
        p = 1
        do i = 1, size(bits)
            if (bits(i) /= 0) then
                p = p + 1
            else
                comp(p) = comp(p) + 1
            end if
        end do
    end function binary_to_composition

    function conjugate(x, sorted) result(y)
        integer, intent(in) :: x(:)
        logical, intent(in), optional :: sorted
        integer, allocatable :: y(:), work(:)
        logical :: already_sorted
        integer :: mx, j

        if (size(x) == 0) then
            allocate(y(0))
            return
        end if
        if (any(x < 0)) error stop "conjugate: entries must be nonnegative"
        already_sorted = .true.
        if (present(sorted)) already_sorted = sorted
        work = x
        if (.not. already_sorted) call sort_descending(work)
        mx = maxval(work)
        allocate(y(mx))
        y = 0
        do j = 1, mx
            y(j) = count(work >= j)
        end do
    end function conjugate

    function durfee(x, sorted) result(value)
        integer, intent(in) :: x(:)
        logical, intent(in), optional :: sorted
        integer :: value
        integer, allocatable :: work(:)
        logical :: already_sorted
        integer :: i

        if (size(x) == 0) then
            value = 0
            return
        end if
        already_sorted = .true.
        if (present(sorted)) already_sorted = sorted
        work = x
        if (.not. already_sorted) call sort_descending(work)
        value = 0
        do i = 1, size(work)
            if (work(i) >= i) then
                value = i
            else
                exit
            end if
        end do
    end function durfee

    pure subroutine sort_descending(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) >= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_descending

    subroutine count_to_default(n8, n, where)
        integer(i8), intent(in) :: n8
        integer, intent(out) :: n
        character(*), intent(in) :: where

        if (n8 < 0_i8 .or. n8 > int(huge(n), i8)) then
            error stop trim(where) // ": result has too many columns for default integer indexing"
        end if
        n = int(n8)
    end subroutine count_to_default

end module partitions_enumeration
