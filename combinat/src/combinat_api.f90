! SPDX-License-Identifier: GPL-2.0-only
module combinat_api
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
    use combinat_kinds, only : dp
    use combinat_rng, only : combinat_rng_state, rng_uniform
    implicit none
    private

    integer, parameter, public :: combinat_success = 0
    integer, parameter, public :: combinat_invalid_argument = 1
    integer, parameter, public :: combinat_size_overflow = 2

    interface combn
        module procedure combn_indices
        module procedure combn_integer
        module procedure combn_real
    end interface combn

    interface combn2
        module procedure combn2_indices
        module procedure combn2_integer
        module procedure combn2_real
    end interface combn2

    interface permn
        module procedure permn_indices
        module procedure permn_integer
        module procedure permn_real
    end interface permn

    interface x2u
        module procedure x2u_indices
        module procedure x2u_integer
        module procedure x2u_real
    end interface x2u

    public :: dp
    public :: combn, combn_indices, combn_integer, combn_real
    public :: combn2, combn2_indices, combn2_integer, combn2_real
    public :: dmnom, fact, logfact, hcube
    public :: ncm, ncm_vec, nsimplex, nsimplex_vec
    public :: permn, permn_indices, permn_integer, permn_real
    public :: rmultinomial, rmultz2
    public :: x2u, x2u_indices, x2u_integer, x2u_real
    public :: xsimplex

contains

    pure elemental real(dp) function fact(x) result(value)
        real(dp), intent(in) :: x !! Real argument whose generalized factorial gamma(x + 1) is returned.

        value = gamma(x + 1.0_dp)
    end function fact

    pure elemental real(dp) function logfact(x) result(value)
        real(dp), intent(in) :: x !! Real argument whose log generalized factorial log-gamma(x + 1) is returned.

        value = log_gamma(x + 1.0_dp)
    end function logfact

    pure elemental real(dp) function ncm(n, m, tol) result(value)
        real(dp), intent(in) :: n !! Upper generalized-binomial argument, including real or negative values.
        real(dp), intent(in) :: m !! Lower argument, with NaN for nonintegers and zero for negative integers.
        real(dp), intent(in), optional :: tol !! Absolute tolerance for integer rounding, defaulting to 1e-8.

        real(dp) :: use_tol

        use_tol = 1.0e-8_dp
        if (present(tol)) use_tol = tol
        value = ncm_raw(n, m, use_tol)
    end function ncm

    pure subroutine ncm_vec(n, m, out, tol, info)
        real(dp), intent(in) :: n(:) !! Vector of upper generalized-binomial arguments, recycled to the output length.
        real(dp), intent(in) :: m(:) !! Recycled lower arguments, normally integer-valued for finite results.
        real(dp), allocatable, intent(out) :: out(:) !! Coefficients with length max(size(n), size(m)).
        real(dp), intent(in), optional :: tol !! Absolute integer-rounding tolerance passed to ncm, defaults to 1e-8.
        integer, intent(out), optional :: info !! Zero on success, invalid argument for an empty input.

        integer :: i
        integer :: nout
        real(dp) :: use_tol

        if (present(info)) info = combinat_success
        if (size(n) == 0 .or. size(m) == 0) then
            allocate(out(0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if

        use_tol = 1.0e-8_dp
        if (present(tol)) use_tol = tol
        nout = max(size(n), size(m))
        allocate(out(nout))
        do i = 1, nout
            out(i) = ncm_raw(n(modulo(i - 1, size(n)) + 1), m(modulo(i - 1, size(m)) + 1), use_tol)
        end do
    end subroutine ncm_vec

    pure real(dp) function dmnom(x, prob, size_arg) result(value)
        real(dp), intent(in) :: x(:) !! Multinomial cell counts, shorter input is recycled when prob has greater length.
        real(dp), intent(in) :: prob(:) !! Cell weights, recycled to max(size(x), size(prob)) and normalized to sum to one.
        real(dp), intent(in), optional :: size_arg !! Total count, defaulting to sum of the original unrecycled x.

        integer :: i
        integer :: p
        real(dp), allocatable :: prob_work(:)
        real(dp) :: size_value
        real(dp) :: sum_x

        if (size(x) == 0 .or. size(prob) == 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        p = max(size(x), size(prob))
        allocate(prob_work(p))
        do i = 1, p
            prob_work(i) = prob(modulo(i - 1, size(prob)) + 1)
        end do
        prob_work = prob_work/sum(prob_work)

        size_value = sum(x)
        if (present(size_arg)) size_value = size_arg
        sum_x = 0.0_dp
        do i = 1, p
            sum_x = sum_x + x(modulo(i - 1, size(x)) + 1)
        end do

        if (sum_x /= size_value) then
            value = 0.0_dp
            return
        end if

        value = logfact(size_value)
        do i = 1, p
            value = value + x(modulo(i - 1, size(x)) + 1)*log(prob_work(i)) &
                - logfact(x(modulo(i - 1, size(x)) + 1))
        end do
        value = exp(value)
    end function dmnom

    pure subroutine combn_indices(n, m, out, info)
        integer, intent(in) :: n !! Number of source elements, valid values are nonnegative integers.
        integer, intent(in) :: m !! Number of elements selected per combination, must satisfy 0 <= m <= n.
        integer, allocatable, intent(out) :: out(:, :) !! Indices as m rows by choose(n,m) columns in upstream order.
        integer, intent(out), optional :: info !! Status code: zero, invalid argument, or output-size overflow.

        integer, allocatable :: a(:)
        integer :: count
        integer :: e
        integer :: h
        integer :: i
        integer :: j
        integer :: local_info
        integer :: nmmp1

        if (present(info)) info = combinat_success
        if (n < 0 .or. m < 0 .or. n < m) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        if (m == 0) then
            allocate(out(0, 0))
            return
        end if

        call combination_count_int(n, m, count, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 0))
            if (present(info)) info = local_info
            return
        end if

        allocate(out(m, count))
        allocate(a(m))
        do j = 1, m
            a(j) = j
        end do
        out(:, 1) = a

        e = 0
        h = m
        i = 2
        nmmp1 = n - m + 1
        do while (a(1) /= nmmp1)
            if (e < n - h) then
                h = 1
                e = a(m)
            else
                h = h + 1
                e = a(m + 1 - h)
            end if
            do j = 1, h
                a(m - h + j) = e + j
            end do
            out(:, i) = a
            i = i + 1
        end do
    end subroutine combn_indices

    pure subroutine combn_integer(x, m, out, info)
        integer, intent(in) :: x(:) !! Integer source vector from which combinations are selected without replacement.
        integer, intent(in) :: m !! Number of elements selected per combination, must satisfy 0 <= m <= size(x).
        integer, allocatable, intent(out) :: out(:, :) !! Selected values, with m rows and choose(size(x),m) columns.
        integer, intent(out), optional :: info !! Status code propagated from combination-index generation.

        integer, allocatable :: indices(:, :)
        integer :: i
        integer :: j
        integer :: local_info

        call combn_indices(size(x), m, indices, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(size(indices, 1), size(indices, 2)))
        do j = 1, size(indices, 2)
            do i = 1, size(indices, 1)
                out(i, j) = x(indices(i, j))
            end do
        end do
        if (present(info)) info = combinat_success
    end subroutine combn_integer

    pure subroutine combn_real(x, m, out, info)
        real(dp), intent(in) :: x(:) !! Real source vector from which combinations are selected without replacement.
        integer, intent(in) :: m !! Number of elements selected per combination, must satisfy 0 <= m <= size(x).
        real(dp), allocatable, intent(out) :: out(:, :) !! Selected values, with m rows and choose(size(x),m) columns.
        integer, intent(out), optional :: info !! Status code propagated from combination-index generation.

        integer, allocatable :: indices(:, :)
        integer :: i
        integer :: j
        integer :: local_info

        call combn_indices(size(x), m, indices, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(size(indices, 1), size(indices, 2)))
        do j = 1, size(indices, 2)
            do i = 1, size(indices, 1)
                out(i, j) = x(indices(i, j))
            end do
        end do
        if (present(info)) info = combinat_success
    end subroutine combn_real

    pure subroutine combn2_indices(n, out, info)
        integer, intent(in) :: n !! Number of elements whose unordered index pairs are requested, must be nonnegative.
        integer, allocatable, intent(out) :: out(:, :) !! Pair matrix with choose(n,2) rows in upstream order.
        integer, intent(out), optional :: info !! Status code: zero, invalid argument, or output-size overflow.

        integer :: count
        integer :: i
        integer :: j
        integer :: k
        integer :: local_info

        if (present(info)) info = combinat_success
        if (n < 0) then
            allocate(out(0, 2))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        if (n < 2) then
            allocate(out(0, 2))
            return
        end if

        call combination_count_int(n, 2, count, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 2))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(count, 2))
        k = 0
        do i = 1, n - 1
            do j = i + 1, n
                k = k + 1
                out(k, 1) = i
                out(k, 2) = j
            end do
        end do
    end subroutine combn2_indices

    pure subroutine combn2_integer(x, out, info)
        integer, intent(in) :: x(:) !! Integer source vector whose unordered value pairs are requested.
        integer, allocatable, intent(out) :: out(:, :) !! Pair matrix with choose(size(x),2) rows and two columns.
        integer, intent(out), optional :: info !! Status code propagated from pair-index generation.

        integer, allocatable :: indices(:, :)
        integer :: i
        integer :: local_info

        call combn2_indices(size(x), indices, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 2))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(size(indices, 1), 2))
        do i = 1, size(indices, 1)
            out(i, 1) = x(indices(i, 1))
            out(i, 2) = x(indices(i, 2))
        end do
        if (present(info)) info = combinat_success
    end subroutine combn2_integer

    pure subroutine combn2_real(x, out, info)
        real(dp), intent(in) :: x(:) !! Real source vector whose unordered value pairs are requested.
        real(dp), allocatable, intent(out) :: out(:, :) !! Pair matrix with choose(size(x),2) rows and two columns.
        integer, intent(out), optional :: info !! Status code propagated from pair-index generation.

        integer, allocatable :: indices(:, :)
        integer :: i
        integer :: local_info

        call combn2_indices(size(x), indices, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 2))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(size(indices, 1), 2))
        do i = 1, size(indices, 1)
            out(i, 1) = x(indices(i, 1))
            out(i, 2) = x(indices(i, 2))
        end do
        if (present(info)) info = combinat_success
    end subroutine combn2_real

    pure subroutine hcube(x, out, scale, translation, info)
        integer, intent(in) :: x(:) !! Positive or zero lattice extents for each dimension, output column one varies fastest.
        real(dp), allocatable, intent(out) :: out(:, :) !! Hypercuboid lattice with product(x) rows and size(x) columns.
        real(dp), intent(in), optional :: scale(:) !! Recycled dimension scale factors applied before translation.
        real(dp), intent(in), optional :: translation(:) !! Recycled dimension translations applied after scaling.
        integer, intent(out), optional :: info !! Status code: zero, invalid argument, or output-size overflow.

        integer(int64) :: block
        integer(int64) :: nrows64
        integer :: i
        integer :: j
        integer :: nrows
        integer :: local_info
        real(dp) :: coordinate

        if (present(info)) info = combinat_success
        if (size(x) == 0 .or. any(x < 0)) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        if (present(scale)) then
            if (size(scale) == 0) then
                allocate(out(0, 0))
                if (present(info)) info = combinat_invalid_argument
                return
            end if
        end if
        if (present(translation)) then
            if (size(translation) == 0) then
                allocate(out(0, 0))
                if (present(info)) info = combinat_invalid_argument
                return
            end if
        end if

        call extent_product(x, nrows64, local_info)
        if (local_info /= combinat_success .or. nrows64 > int(huge(0), int64)) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_size_overflow
            return
        end if
        nrows = int(nrows64)
        allocate(out(nrows, size(x)))
        if (nrows == 0) return

        block = 1_int64
        do j = 1, size(x)
            do i = 1, nrows
                coordinate = real(modulo((int(i - 1, int64)/block), int(x(j), int64)) + 1_int64, dp)
                if (present(scale)) coordinate = coordinate*scale(modulo(j - 1, size(scale)) + 1)
                if (present(translation)) coordinate = coordinate + translation(modulo(j - 1, size(translation)) + 1)
                out(i, j) = coordinate
            end do
            block = block*int(x(j), int64)
        end do
    end subroutine hcube

    pure elemental real(dp) function nsimplex(p, n) result(value)
        real(dp), intent(in) :: p !! Number of composition parts, negative values are assigned zero as in the upstream R function.
        real(dp), intent(in) :: n !! Composition total, passed as the lower argument of nCm and therefore normally integer-valued.

        value = ncm(n + p - 1.0_dp, n)
        if (p < 0.0_dp) value = 0.0_dp
    end function nsimplex

    pure subroutine nsimplex_vec(p, n, out, info)
        real(dp), intent(in) :: p(:) !! Vector of part counts, recycled to the output length.
        real(dp), intent(in) :: n(:) !! Vector of composition totals, recycled to the output length.
        real(dp), allocatable, intent(out) :: out(:) !! Counts of simplex-lattice points with R-style recycling.
        integer, intent(out), optional :: info !! Status code: zero on success or invalid argument for an empty input vector.

        integer :: i
        integer :: nout

        if (present(info)) info = combinat_success
        if (size(p) == 0 .or. size(n) == 0) then
            allocate(out(0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        nout = max(size(p), size(n))
        allocate(out(nout))
        do i = 1, nout
            out(i) = nsimplex(p(modulo(i - 1, size(p)) + 1), n(modulo(i - 1, size(n)) + 1))
        end do
    end subroutine nsimplex_vec

    pure subroutine xsimplex(p, n, out, info)
        integer, intent(in) :: p !! Number of nonnegative composition parts, values below one return an empty result.
        integer, intent(in) :: n !! Nonnegative integer total distributed among p parts, negative values return an empty result.
        integer, allocatable, intent(out) :: out(:, :) !! Simplex points as p rows by choose(n+p-1,n) columns in upstream order.
        integer, intent(out), optional :: info !! Zero for valid or empty output, or size overflow for impractical allocation.

        integer :: count
        integer :: i
        integer :: local_info
        integer :: p1
        integer :: target
        integer, allocatable :: x(:)

        if (present(info)) info = combinat_success
        if (p < 1 .or. n < 0) then
            allocate(out(0, 0))
            return
        end if

        call combination_count_int(n + p - 1, n, count, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(p, count))
        allocate(x(p))
        x = 0
        x(1) = n
        if (p == 1 .or. n == 0) then
            out(:, 1) = x
            return
        end if

        p1 = p - 1
        target = 1
        i = 0
        do
            i = i + 1
            out(:, i) = x
            x(target) = x(target) - 1
            if (target < p1) then
                target = target + 1
                x(target) = 1 + x(p)
                x(p) = 0
            else
                x(p) = x(p) + 1
                do while (x(target) == 0)
                    target = target - 1
                    if (target == 0) then
                        i = i + 1
                        out(:, i) = x
                        return
                    end if
                end do
            end if
        end do
    end subroutine xsimplex

    pure subroutine permn_indices(n, out, info)
        integer, intent(in) :: n !! Number of elements in the sequence 1:n, must be nonnegative.
        integer, allocatable, intent(out) :: out(:, :) !! Permutations as n rows by n! columns in upstream minimal-change order.
        integer, intent(out), optional :: info !! Status code: zero, invalid argument, or output-size overflow.

        integer, allocatable :: d(:)
        integer, allocatable :: ip(:)
        integer, allocatable :: p(:)
        integer :: count
        integer :: index1
        integer :: index2
        integer :: i
        integer :: j
        integer :: local_info
        integer :: m
        integer :: temp

        if (present(info)) info = combinat_success
        if (n < 0) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        call factorial_count_int(n, count, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(n, count))
        if (n == 0) return

        allocate(d(n), ip(n), p(0:n + 1))
        do j = 1, n
            ip(j) = j
            p(j) = j
            d(j) = -1
        end do
        d(1) = 0
        p(0) = n + 1
        p(n + 1) = n + 1

        m = n + 1
        i = 1
        do while (m /= 1)
            out(:, i) = p(1:n)
            i = i + 1
            m = 1
            do j = n, 1, -1
                if (p(ip(j) + d(j)) <= j) then
                    m = j
                    exit
                end if
            end do
            if (m < n) d(m + 1:n) = -d(m + 1:n)
            index1 = ip(m)
            index2 = p(index1 + d(m))
            p(index1) = index2
            p(index1 + d(m)) = m
            temp = ip(index2)
            ip(index2) = ip(m)
            ip(m) = temp
        end do
    end subroutine permn_indices

    pure subroutine permn_integer(x, out, info)
        integer, intent(in) :: x(:) !! Integer source vector whose complete set of permutations is requested.
        integer, allocatable, intent(out) :: out(:, :) !! Permuted values as size(x) rows by size(x)! columns.
        integer, intent(out), optional :: info !! Status code propagated from permutation-index generation.

        integer, allocatable :: indices(:, :)
        integer :: i
        integer :: j
        integer :: local_info

        call permn_indices(size(x), indices, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(size(indices, 1), size(indices, 2)))
        do j = 1, size(indices, 2)
            do i = 1, size(indices, 1)
                out(i, j) = x(indices(i, j))
            end do
        end do
        if (present(info)) info = combinat_success
    end subroutine permn_integer

    pure subroutine permn_real(x, out, info)
        real(dp), intent(in) :: x(:) !! Real source vector whose complete set of permutations is requested.
        real(dp), allocatable, intent(out) :: out(:, :) !! Permuted values as size(x) rows by size(x)! columns.
        integer, intent(out), optional :: info !! Status code propagated from permutation-index generation.

        integer, allocatable :: indices(:, :)
        integer :: i
        integer :: j
        integer :: local_info

        call permn_indices(size(x), indices, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0, 0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(size(indices, 1), size(indices, 2)))
        do j = 1, size(indices, 2)
            do i = 1, size(indices, 1)
                out(i, j) = x(indices(i, j))
            end do
        end do
        if (present(info)) info = combinat_success
    end subroutine permn_real

    subroutine rmultinomial(n, p, rng, out, rows, info)
        integer, intent(in) :: n(:) !! Recycled multinomial trial counts, all nonnegative.
        real(dp), intent(in) :: p(:, :) !! Probability weights with distributions in rows, recycled as needed.
        type(combinat_rng_state), intent(inout) :: rng !! Explicit reproducible RNG state advanced for every categorical draw.
        integer, allocatable, intent(out) :: out(:, :) !! Simulated count matrix with output rows by size(p,2) categories.
        integer, intent(in), optional :: rows !! Number of distributions to draw, defaults to max(size(n), size(p,1)).
        integer, intent(out), optional :: info !! Zero on success or invalid argument for malformed inputs.

        integer :: category
        integer :: draw
        integer :: i
        integer :: nrows
        integer :: source_row
        real(dp), allocatable :: probability(:)

        if (present(info)) info = combinat_success
        if (size(n) == 0 .or. size(p, 1) == 0 .or. size(p, 2) == 0) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        nrows = max(size(n), size(p, 1))
        if (present(rows)) nrows = rows
        if (nrows < 0) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        allocate(out(nrows, size(p, 2)))
        out = 0
        if (nrows == 0) return
        allocate(probability(size(p, 2)))

        do i = 1, nrows
            if (n(modulo(i - 1, size(n)) + 1) < 0) then
                out = 0
                if (present(info)) info = combinat_invalid_argument
                return
            end if
            source_row = modulo(i - 1, size(p, 1)) + 1
            probability = p(source_row, :)
            if (.not. valid_probability_weights(probability)) then
                out = 0
                if (present(info)) info = combinat_invalid_argument
                return
            end if
            probability = probability/sum(probability)
            do draw = 1, n(modulo(i - 1, size(n)) + 1)
                category = categorical_draw(probability, rng)
                out(i, category) = out(i, category) + 1
            end do
        end do
    end subroutine rmultinomial

    subroutine rmultz2(n, p, rng, out, draws, info)
        integer, intent(in) :: n(:) !! Recycled multinomial trial counts for requested draws, all nonnegative.
        real(dp), intent(in) :: p(:) !! Fixed nonnegative probability weights with positive finite sum.
        type(combinat_rng_state), intent(inout) :: rng !! Explicit reproducible RNG state advanced for every categorical draw.
        integer, allocatable, intent(out) :: out(:, :) !! Category rows by draw columns, matching upstream orientation.
        integer, intent(in), optional :: draws !! Number of multinomial draws, defaults to size(n).
        integer, intent(out), optional :: info !! Zero on success or invalid argument for malformed inputs.

        integer :: category
        integer :: draw
        integer :: j
        integer :: ndraws
        real(dp), allocatable :: probability(:)

        if (present(info)) info = combinat_success
        if (size(n) == 0 .or. size(p) == 0) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        ndraws = size(n)
        if (present(draws)) ndraws = draws
        if (ndraws < 0 .or. .not. valid_probability_weights(p)) then
            allocate(out(0, 0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        allocate(out(size(p), ndraws))
        out = 0
        if (ndraws == 0) return
        probability = p/sum(p)

        do j = 1, ndraws
            if (n(modulo(j - 1, size(n)) + 1) < 0) then
                out = 0
                if (present(info)) info = combinat_invalid_argument
                return
            end if
            do draw = 1, n(modulo(j - 1, size(n)) + 1)
                category = categorical_draw(probability, rng)
                out(category, j) = out(category, j) + 1
            end do
        end do
    end subroutine rmultz2

    pure subroutine x2u_indices(x, out, info)
        integer, intent(in) :: x(:) !! Nonnegative bin counts whose expanded bin numbers are requested.
        integer, allocatable, intent(out) :: out(:) !! Expanded 1-based bin indices, of length sum(x), in increasing bin order.
        integer, intent(out), optional :: info !! Status code: zero, invalid argument, or output-size overflow.

        integer :: i
        integer :: j
        integer :: k
        integer :: total
        integer :: local_info

        call count_sum_int(x, total, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(total))
        k = 0
        do i = 1, size(x)
            do j = 1, x(i)
                k = k + 1
                out(k) = i
            end do
        end do
        if (present(info)) info = combinat_success
    end subroutine x2u_indices

    pure subroutine x2u_integer(x, labels, out, info)
        integer, intent(in) :: x(:) !! Nonnegative bin counts whose labels are expanded according to their multiplicities.
        integer, intent(in) :: labels(:) !! Integer bin labels, length must equal size(x).
        integer, allocatable, intent(out) :: out(:) !! Expanded integer labels, of length sum(x), in bin order.
        integer, intent(out), optional :: info !! Status code: zero, invalid argument, or output-size overflow.

        integer :: i
        integer :: j
        integer :: k
        integer :: total
        integer :: local_info

        if (size(labels) /= size(x)) then
            allocate(out(0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        call count_sum_int(x, total, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(total))
        k = 0
        do i = 1, size(x)
            do j = 1, x(i)
                k = k + 1
                out(k) = labels(i)
            end do
        end do
        if (present(info)) info = combinat_success
    end subroutine x2u_integer

    pure subroutine x2u_real(x, labels, out, info)
        integer, intent(in) :: x(:) !! Nonnegative bin counts whose labels are expanded according to their multiplicities.
        real(dp), intent(in) :: labels(:) !! Real bin labels, length must equal size(x).
        real(dp), allocatable, intent(out) :: out(:) !! Expanded real labels, of length sum(x), in bin order.
        integer, intent(out), optional :: info !! Status code: zero, invalid argument, or output-size overflow.

        integer :: i
        integer :: j
        integer :: k
        integer :: total
        integer :: local_info

        if (size(labels) /= size(x)) then
            allocate(out(0))
            if (present(info)) info = combinat_invalid_argument
            return
        end if
        call count_sum_int(x, total, local_info)
        if (local_info /= combinat_success) then
            allocate(out(0))
            if (present(info)) info = local_info
            return
        end if
        allocate(out(total))
        k = 0
        do i = 1, size(x)
            do j = 1, x(i)
                k = k + 1
                out(k) = labels(i)
            end do
        end do
        if (present(info)) info = combinat_success
    end subroutine x2u_real

    pure recursive real(dp) function ncm_raw(n, m, tol) result(value)
        real(dp), intent(in) :: n !! Scalar upper generalized-binomial argument used internally.
        real(dp), intent(in) :: m !! Scalar lower generalized-binomial argument used internally and required to be integer-valued.
        real(dp), intent(in) :: tol !! Absolute tolerance for optional final rounding to an integer.

        integer :: j
        integer :: mi
        integer :: negatives
        real(dp) :: nearest
        real(dp) :: term

        if (ieee_is_nan(n) .or. ieee_is_nan(m)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (m /= aint(m)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (m == 0.0_dp) then
            value = 1.0_dp
            return
        end if
        if (m < 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(m) .or. m > real(huge(0), dp)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        mi = int(m)

        if (n < 0.0_dp) then
            value = ncm_raw(m - n - 1.0_dp, m, tol)
            if (modulo(mi, 2) == 1) value = -value
        else if (n == 0.0_dp) then
            value = 0.0_dp
        else if (.not. ieee_is_finite(n)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
        else if (n /= aint(n) .and. n < m) then
            value = 0.0_dp
            negatives = 0
            do j = 1, mi
                term = n - m + real(j, dp)
                if (term < 0.0_dp) negatives = negatives + 1
                if (term == 0.0_dp) then
                    value = -huge(1.0_dp)
                else if (value > -huge(1.0_dp)) then
                    value = value + log(abs(term))
                end if
            end do
            value = exp(value - log_gamma(m + 1.0_dp))
            if (modulo(negatives, 2) == 1) value = -value
        else if (n >= m) then
            value = exp(log_gamma(n + 1.0_dp) - log_gamma(m + 1.0_dp) - log_gamma(n - m + 1.0_dp))
        else
            value = 0.0_dp
        end if

        if (.not. ieee_is_nan(value) .and. ieee_is_finite(value)) then
            nearest = anint(value)
            if (abs(nearest - value) < tol) value = nearest
        end if
    end function ncm_raw

    pure subroutine combination_count_int(n, m, count, info)
        integer, intent(in) :: n !! Nonnegative source size used for an exact-size combination allocation.
        integer, intent(in) :: m !! Selection size satisfying 0 <= m <= n.
        integer, intent(out) :: count !! choose(n,m) converted to default integer when representable.
        integer, intent(out) :: info !! Zero on success, invalid argument, or size overflow.

        real(dp) :: value

        count = 0
        if (n < 0 .or. m < 0 .or. m > n) then
            info = combinat_invalid_argument
            return
        end if
        value = ncm(real(n, dp), real(m, dp), 0.0_dp)
        if (.not. ieee_is_finite(value) .or. value > real(huge(0), dp)) then
            info = combinat_size_overflow
            return
        end if
        count = int(anint(value))
        info = combinat_success
    end subroutine combination_count_int

    pure subroutine factorial_count_int(n, count, info)
        integer, intent(in) :: n !! Nonnegative integer whose factorial determines the number of permutations.
        integer, intent(out) :: count !! n! as a default integer when representable.
        integer, intent(out) :: info !! Zero on success, invalid argument, or size overflow.

        integer :: i

        count = 1
        if (n < 0) then
            info = combinat_invalid_argument
            return
        end if
        do i = 2, n
            if (count > huge(count)/i) then
                count = 0
                info = combinat_size_overflow
                return
            end if
            count = count*i
        end do
        info = combinat_success
    end subroutine factorial_count_int

    pure subroutine extent_product(x, total, info)
        integer, intent(in) :: x(:) !! Nonnegative dimension extents whose product is required for a lattice allocation.
        integer(int64), intent(out) :: total !! Product of x in 64-bit integer arithmetic when representable.
        integer, intent(out) :: info !! Zero on success, invalid argument, or 64-bit size overflow.

        integer :: i

        total = 1_int64
        if (any(x < 0)) then
            info = combinat_invalid_argument
            return
        end if
        do i = 1, size(x)
            if (x(i) == 0) then
                total = 0_int64
                info = combinat_success
                return
            end if
            if (total > huge(total)/int(x(i), int64)) then
                total = 0_int64
                info = combinat_size_overflow
                return
            end if
            total = total*int(x(i), int64)
        end do
        info = combinat_success
    end subroutine extent_product

    pure subroutine count_sum_int(x, total, info)
        integer, intent(in) :: x(:) !! Nonnegative multiplicities whose sum determines an expanded-vector length.
        integer, intent(out) :: total !! Sum of multiplicities as a default integer when representable.
        integer, intent(out) :: info !! Zero on success, invalid argument, or size overflow.

        integer :: i

        total = 0
        if (any(x < 0)) then
            info = combinat_invalid_argument
            return
        end if
        do i = 1, size(x)
            if (total > huge(total) - x(i)) then
                total = 0
                info = combinat_size_overflow
                return
            end if
            total = total + x(i)
        end do
        info = combinat_success
    end subroutine count_sum_int

    pure logical function valid_probability_weights(probability) result(valid)
        real(dp), intent(in) :: probability(:) !! Finite nonnegative weights with positive finite sum.

        real(dp) :: total

        valid = size(probability) > 0
        if (.not. valid) return
        valid = all(ieee_is_finite(probability)) .and. all(probability >= 0.0_dp)
        if (.not. valid) return
        total = sum(probability)
        valid = ieee_is_finite(total) .and. total > 0.0_dp
    end function valid_probability_weights

    integer function categorical_draw(probability, rng) result(category)
        real(dp), intent(in) :: probability(:) !! Normalized categorical probabilities with nonnegative entries summing to one.
        type(combinat_rng_state), intent(inout) :: rng !! RNG state advanced once to choose a category.

        integer :: i
        real(dp) :: cumulative
        real(dp) :: u

        u = rng_uniform(rng)
        cumulative = 0.0_dp
        category = size(probability)
        do i = 1, size(probability)
            cumulative = cumulative + probability(i)
            if (u <= cumulative) then
                category = i
                return
            end if
        end do
    end function categorical_draw

end module combinat_api
