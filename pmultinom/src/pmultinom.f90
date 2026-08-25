! Modern Fortran translation of the computational code in the R package
! pmultinom 1.0.0 by Alexander Davis.
!
! Upstream algorithm: Levin (1981), conditioning independent Poisson counts on
! their sum to obtain a multinomial distribution.
!
! Upstream license: GNU Affero General Public License v3 (AGPL-3).
!
! The R implementation normalizes each restricted Poisson distribution,
! convolves the normalized distributions, and later multiplies the normalizing
! constants back in. This translation performs the algebraically equivalent
! computation on scaled, *unnormalized* restricted Poisson masses. It avoids a
! cancellation that is especially delicate for very small truncation events.
module pmultinom_module
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
    use pmultinom_kinds, only : dp
    use pmultinom_fft, only : convolve_prefix
    use pmultinom_math, only : integer_support, log_poisson_pmf
    implicit none
    private

    public :: dp
    public :: pmultinom
    public :: pmultinom_many
    public :: invert_pmultinom
    public :: invert_pmultinom_many

contains

    function pmultinom(size, probs, lower, upper) result(probability)
        integer, intent(in) :: size
        real(dp), intent(in) :: probs(:)
        real(dp), intent(in), optional :: lower(:), upper(:)
        real(dp) :: probability

        real(dp), allocatable :: lo(:), hi(:), p(:)
        integer :: stat
        character(len=256) :: msg

        if (size < 0) error stop "pmultinom: size must be nonnegative"
        call prepare_inputs(probs, lower, upper, lo, hi, p, stat, msg)
        if (stat /= 0) error stop trim(msg)

        probability = probability_at_size(lo, hi, p, size)
    end function pmultinom

    function pmultinom_many(sizes, probs, lower, upper) result(probabilities)
        integer, intent(in) :: sizes(:)
        real(dp), intent(in) :: probs(:)
        real(dp), intent(in), optional :: lower(:), upper(:)
        real(dp), allocatable :: probabilities(:)

        real(dp), allocatable :: lo(:), hi(:), p(:)
        integer :: i, stat
        character(len=256) :: msg

        call prepare_inputs(probs, lower, upper, lo, hi, p, stat, msg)
        if (stat /= 0) error stop trim(msg)

        allocate(probabilities(size(sizes)))
        do i = 1, size(sizes)
            if (sizes(i) < 0) error stop "pmultinom_many: sizes must be nonnegative"
            probabilities(i) = probability_at_size(lo, hi, p, sizes(i))
        end do
    end function pmultinom_many

    function invert_pmultinom(probs, target_prob, lower, upper) result(sample_size)
        real(dp), intent(in) :: probs(:)
        real(dp), intent(in) :: target_prob
        real(dp), intent(in), optional :: lower(:), upper(:)
        integer :: sample_size

        integer, allocatable :: tmp(:)

        tmp = invert_pmultinom_many(probs, [target_prob], lower, upper)
        sample_size = tmp(1)
    end function invert_pmultinom

    function invert_pmultinom_many(probs, target_probs, lower, upper) result(sample_sizes)
        real(dp), intent(in) :: probs(:)
        real(dp), intent(in) :: target_probs(:)
        real(dp), intent(in), optional :: lower(:), upper(:)
        integer, allocatable :: sample_sizes(:)

        real(dp), allocatable :: lo(:), hi(:), p(:), path(:), block(:)
        logical :: lower_absent, upper_absent, increasing
        logical, allocatable :: done(:)
        integer :: i, n, old_n, v, a, stat, max_iter
        character(len=256) :: msg

        call prepare_inputs(probs, lower, upper, lo, hi, p, stat, msg)
        if (stat /= 0) error stop trim(msg)

        if (any(target_probs < 0.0_dp .or. target_probs > 1.0_dp)) then
            error stop "invert_pmultinom: target probabilities must lie in [0, 1]"
        end if

        lower_absent = all(lo < 0.0_dp .and. .not. is_finite_bound(lo))
        upper_absent = all(hi > 0.0_dp .and. .not. is_finite_bound(hi))

        if (.not. lower_absent .and. .not. upper_absent) then
            error stop "invert_pmultinom: lower and upper cannot both be given"
        end if

        allocate(sample_sizes(size(target_probs)))
        if (lower_absent .and. upper_absent) then
            sample_sizes = 0
            return
        end if

        increasing = upper_absent
        if (increasing) then
            do i = 1, size(p)
                if (p(i) == 0.0_dp .and. lo(i) >= 0.0_dp) then
                    error stop "invert_pmultinom: positive count required in a zero-probability category; inverse does not exist"
                end if
            end do
        end if

        allocate(done(size(target_probs)))
        done = .false.
        sample_sizes = -1

        ! Match the upstream stepping scheme: 9, 35, 80, ... in the tuning
        ! parameter. Each block uses one Poisson scale and computes the whole
        ! multinomial path over the newly covered sample sizes.
        a = 3
        v = a
        n = v * v
        old_n = -1
        allocate(path(0:n))
        call probability_block(lo, hi, p, n, 0, n, path)
        call record_crossings(path, 0, target_probs, increasing, done, sample_sizes)

        max_iter = 100000
        do i = 1, max_iter
            if (all(done)) exit
            old_n = n
            v = v + a
            if (v > 1000000) error stop "invert_pmultinom: search range became too large"
            n = v * v - 1
            allocate(block(old_n + 1:n))
            call probability_block(lo, hi, p, n, old_n + 1, n, block)
            call record_crossings(block, old_n + 1, target_probs, increasing, done, sample_sizes)
            deallocate(block)
        end do
        if (.not. all(done)) error stop "invert_pmultinom: inverse search failed to converge"
    end function invert_pmultinom_many

    subroutine record_crossings(values, first_size, targets, increasing, done, answers)
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: first_size
        real(dp), intent(in) :: targets(:)
        logical, intent(in) :: increasing
        logical, intent(inout) :: done(:)
        integer, intent(inout) :: answers(:)

        integer :: i, j, nsize
        logical :: crossed

        do j = 1, size(values)
            nsize = first_size + j - 1
            do i = 1, size(targets)
                if (done(i)) cycle
                if (increasing) then
                    crossed = values(j) >= targets(i)
                else
                    crossed = values(j) <= targets(i)
                end if
                if (crossed) then
                    answers(i) = nsize
                    done(i) = .true.
                end if
            end do
            if (all(done)) return
        end do
    end subroutine record_crossings

    real(dp) function probability_at_size(lo, hi, p, size_n) result(ans)
        real(dp), intent(in) :: lo(:), hi(:), p(:)
        integer, intent(in) :: size_n

        real(dp), allocatable :: value(:)
        integer :: v, tune_n

        if (event_impossible_at_size(lo, hi, size_n)) then
            ans = 0.0_dp
            return
        end if
        if (event_certain_at_size(lo, hi, size_n)) then
            ans = 1.0_dp
            return
        end if
        if (size(p) == 1) then
            if (lo(1) < real(size_n, dp) .and. real(size_n, dp) <= hi(1)) then
                ans = 1.0_dp
            else
                ans = 0.0_dp
            end if
            return
        end if

        v = 3
        tune_n = v * v - 1
        do while (tune_n < size_n)
            v = v + 3
            tune_n = v * v - 1
        end do
        allocate(value(1))
        call probability_block(lo, hi, p, tune_n, size_n, size_n, value)
        ans = value(1)
    end function probability_at_size

    subroutine probability_block(lo, hi, p, tune_n, first_size, last_size, values)
        real(dp), intent(in) :: lo(:), hi(:), p(:)
        integer, intent(in) :: tune_n, first_size, last_size
        real(dp), intent(out) :: values(:)

        real(dp), allocatable :: conv(:), next(:), q(:)
        real(dp) :: lambda, log_peak, lp, max_conv, log_scale, log_ans
        integer :: j, x, ilo, ihi, m, out_index
        logical :: initialized

        if (size(values) /= last_size - first_size + 1) then
            error stop "probability_block: output has wrong size"
        end if

        if (tune_n <= 0) error stop "probability_block: tuning size must be positive"

        allocate(conv(0:tune_n), next(0:tune_n), q(0:tune_n))
        conv = 0.0_dp
        log_scale = 0.0_dp
        initialized = .false.

        do j = 1, size(p)
            lambda = real(tune_n, dp) * p(j)
            call integer_support(lo(j), hi(j), tune_n, ilo, ihi)
            if (ilo > ihi) then
                values = 0.0_dp
                return
            end if

            q = 0.0_dp
            log_peak = -huge(1.0_dp)
            do x = ilo, ihi
                lp = log_poisson_pmf(x, lambda)
                if (lp > log_peak) log_peak = lp
            end do
            if (log_peak <= -0.5_dp * huge(1.0_dp)) then
                values = 0.0_dp
                return
            end if
            do x = ilo, ihi
                lp = log_poisson_pmf(x, lambda)
                if (lp > log_peak - 745.0_dp) q(x) = exp(lp - log_peak)
            end do
            log_scale = log_scale + log_peak

            if (.not. initialized) then
                conv = q
                initialized = .true.
            else
                call convolve_prefix(conv, q, tune_n, next)
                max_conv = maxval(next)
                if (max_conv <= 0.0_dp) then
                    values = 0.0_dp
                    return
                end if
                next = max(next, 0.0_dp)
                next = next / max_conv
                log_scale = log_scale + log(max_conv)
                conv = next
            end if
        end do

        out_index = 0
        do m = first_size, last_size
            out_index = out_index + 1
            if (event_impossible_at_size(lo, hi, m)) then
                values(out_index) = 0.0_dp
            else if (event_certain_at_size(lo, hi, m)) then
                values(out_index) = 1.0_dp
            else if (conv(m) <= 0.0_dp) then
                values(out_index) = 0.0_dp
            else
                log_ans = log(conv(m)) + log_scale - log_poisson_pmf(m, real(tune_n, dp))
                if (log_ans <= log(tiny(1.0_dp))) then
                    values(out_index) = 0.0_dp
                else
                    values(out_index) = exp(min(log_ans, 0.0_dp))
                    values(out_index) = min(1.0_dp, max(0.0_dp, values(out_index)))
                end if
            end if
        end do
    end subroutine probability_block

    logical function event_impossible_at_size(lo, hi, n) result(impossible)
        real(dp), intent(in) :: lo(:), hi(:)
        integer, intent(in) :: n

        integer :: j, ilo, ihi
        integer(kind=8) :: min_sum, max_sum

        min_sum = 0_8
        max_sum = 0_8
        do j = 1, size(lo)
            call integer_support(lo(j), hi(j), n, ilo, ihi)
            if (ilo > ihi) then
                impossible = .true.
                return
            end if
            min_sum = min_sum + int(ilo, kind=8)
            max_sum = max_sum + int(ihi, kind=8)
        end do
        impossible = int(n, kind=8) < min_sum .or. int(n, kind=8) > max_sum
    end function event_impossible_at_size

    logical function event_certain_at_size(lo, hi, n) result(certain)
        real(dp), intent(in) :: lo(:), hi(:)
        integer, intent(in) :: n

        integer :: j, ilo, ihi

        certain = .true.
        do j = 1, size(lo)
            call integer_support(lo(j), hi(j), n, ilo, ihi)
            if (ilo /= 0 .or. ihi /= n) then
                certain = .false.
                return
            end if
        end do
    end function event_certain_at_size

    subroutine prepare_inputs(probs, lower, upper, lo, hi, p, stat, msg)
        real(dp), intent(in) :: probs(:)
        real(dp), intent(in), optional :: lower(:), upper(:)
        real(dp), allocatable, intent(out) :: lo(:), hi(:), p(:)
        integer, intent(out) :: stat
        character(len=*), intent(out) :: msg

        integer :: np, nl, nu, ncat, i
        real(dp) :: pos_inf, neg_inf, psum, tol

        stat = 0
        msg = ""
        np = size(probs)
        if (np == 0) then
            stat = 1
            msg = "pmultinom: probs must not be empty"
            return
        end if
        nl = 1
        if (present(lower)) nl = size(lower)
        nu = 1
        if (present(upper)) nu = size(upper)
        if (nl == 0 .or. nu == 0) then
            stat = 1
            msg = "pmultinom: bounds must not be empty"
            return
        end if

        ncat = max(np, max(nl, nu))
        if (mod(ncat, np) /= 0 .or. mod(ncat, nl) /= 0 .or. mod(ncat, nu) /= 0) then
            stat = 1
            msg = "pmultinom: longer argument is not a multiple of shorter argument length"
            return
        end if

        allocate(lo(ncat), hi(ncat), p(ncat))
        pos_inf = ieee_value(1.0_dp, ieee_positive_inf)
        neg_inf = ieee_value(1.0_dp, ieee_negative_inf)
        do i = 1, ncat
            p(i) = probs(1 + mod(i - 1, np))
            if (present(lower)) then
                lo(i) = lower(1 + mod(i - 1, nl))
            else
                lo(i) = neg_inf
            end if
            if (present(upper)) then
                hi(i) = upper(1 + mod(i - 1, nu))
            else
                hi(i) = pos_inf
            end if
        end do

        if (any(p < 0.0_dp .or. p > 1.0_dp)) then
            stat = 1
            msg = "pmultinom: probs must lie in [0, 1]"
            return
        end if
        if (any(lo >= hi)) then
            ! Retain this information for the caller; a mathematically empty
            ! event is represented naturally as probability zero. No error.
        end if

        psum = sum(p)
        tol = 128.0_dp * epsilon(1.0_dp) * max(1.0_dp, real(ncat, dp))
        if (abs(psum - 1.0_dp) > tol) then
            stat = 1
            write(msg, '(a,es24.16)') "pmultinom: probabilities do not sum to one; sum = ", psum
            return
        end if
    end subroutine prepare_inputs

    elemental logical function is_finite_bound(x) result(finite)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
        real(dp), intent(in) :: x
        finite = ieee_is_finite(x)
    end function is_finite_bound

end module pmultinom_module
