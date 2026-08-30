! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_core
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_negative_inf, ieee_positive_inf
    use classint_kinds, only: dp
    use classint_types, only: class_intervals, classint_options
    use classint_utils, only: finite_values, unique_sorted, mean_dp, sample_sd, quantile_r, sturges_classes, &
                              sample_without_replacement, sort_real, lower_text
    use classint_pretty, only: pretty_breaks
    use classint_dpih, only: dpih_bandwidth
    use classint_fisher, only: fisher_exact, jenks_breaks
    use classint_cluster, only: kmeans_breaks, hclust_fit_1d, hclust_cut_breaks, bclust_fit_breaks, bclust_cut_breaks
    use e1071, only: rng_state, rng_seed
    implicit none
    private

    interface classint_fit
        module procedure classint_fit_n
        module procedure classint_fit_auto
    end interface classint_fit

    public :: classint_fit, get_hclust_class_intervals, get_bclust_class_intervals

contains

    subroutine classint_fit_n(var, n, style, result, options, interval_closure)
        real(dp), intent(in) :: var(:) !! Numeric variable; non-finite values are retained but omitted when finding breaks.
        integer, intent(in) :: n !! Requested class count; values above the unique finite count are reduced automatically.
        character(len=*), intent(in) :: style !! Classification style name, such as quantile, fisher, jenks, equal, kmeans, or box.
        type(class_intervals), intent(out) :: result !! Fitted intervals, original values, metadata, and optional cluster model.
        type(classint_options), intent(in), optional :: options !! Style-specific controls; omitted values use package defaults.
        character(len=*), intent(in), optional :: interval_closure !! Interval closure: "left" or "right"; Jenks forces right.

        call classint_fit_impl(var, n, style, result, options, interval_closure)
    end subroutine classint_fit_n

    subroutine classint_fit_auto(var, style, result, options, interval_closure)
        real(dp), intent(in) :: var(:) !! Numeric variable to partition; non-finite values are ignored for break estimation.
        character(len=*), intent(in) :: style !! Classification style; n-free styles ignore the internally computed Sturges count.
        type(class_intervals), intent(out) :: result !! Fitted intervals with original values and retained clustering state.
        type(classint_options), intent(in), optional :: options !! Style-specific controls; omitted values use package defaults.
        character(len=*), intent(in), optional :: interval_closure !! Requested closure, "left" or "right"; Jenks forces right.
        real(dp), allocatable :: finite(:)
        integer :: n

        finite = finite_values(var)
        if (size(finite) < 1) error stop "classint_fit: no finite observations"
        n = sturges_classes(size(finite))
        call classint_fit_impl(var, n, style, result, options, interval_closure)
    end subroutine classint_fit_auto

    subroutine classint_fit_impl(var, n_requested, style, result, options, interval_closure)
        real(dp), intent(in) :: var(:) !! Original numeric variable, including non-finite entries preserved for classification.
        integer, intent(in) :: n_requested !! Requested or Sturges-derived class count before unique-value adjustment.
        character(len=*), intent(in) :: style !! Classification style dispatched by this computational core.
        type(class_intervals), intent(out) :: result !! Fully populated native classInt interval object.
        type(classint_options), intent(in), optional :: options !! Optional algorithm controls copied before fitting.
        character(len=*), intent(in), optional :: interval_closure !! Caller-selected closure semantics for all styles except Jenks.
        type(classint_options) :: opts
        type(rng_state) :: rng
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: ux(:)
        real(dp), allocatable :: work_x(:)
        real(dp), allocatable :: sampled(:)
        real(dp), allocatable :: sbrks(:)
        real(dp), allocatable :: tmp(:)
        real(dp) :: mu
        real(dp) :: sd
        real(dp) :: h
        real(dp) :: range_x(2)
        character(len=:), allocatable :: use_style
        character(len=:), allocatable :: closure
        integer :: n
        integer :: i
        integer :: nsamp
        logical :: needs_n

        opts = classint_options()
        if (present(options)) opts = options
        call rng_seed(rng, opts%seed)
        x = finite_values(var)
        if (size(x) < 1) error stop "classint_fit: no finite values"
        ux = unique_sorted(x)
        if (size(ux) == 1) error stop "classint_fit: single unique finite value"
        use_style = trim(adjustl(lower_text(style)))
        closure = "left"
        if (present(interval_closure)) closure = trim(adjustl(lower_text(interval_closure)))
        if (closure /= "left" .and. closure /= "right") error stop "classint_fit: closure must be left or right"
        needs_n = use_style /= "dpih" .and. use_style /= "headtails" .and. use_style /= "box"
        n = n_requested
        if (needs_n .and. n < 2) error stop "classint_fit: n must be at least two"
        if (needs_n .and. n > size(ux)) n = size(ux)

        result%values = var
        result%nobs = size(ux)
        result%sampled = .false.
        result%has_hclust = .false.
        result%has_bclust = .false.
        if (use_style == "jenks") closure = "right"
        result%interval_closure = closure

        if (needs_n .and. n == size(ux)) then
            call unique_breaks(ux, result%breaks)
            result%style = "unique"
            return
        end if

        select case (use_style)
        case ("fixed")
            if (.not. allocated(opts%fixed_breaks)) error stop "classint_fit: fixed style requires fixed_breaks"
            result%breaks = opts%fixed_breaks
            call sort_real(result%breaks)
        case ("sd")
            mu = mean_dp(x)
            sd = sample_sd(x)
            if (sd <= 0.0_dp) error stop "classint_fit: zero scale in sd style"
            if (allocated(opts%sd_m)) then
                sbrks = opts%sd_m
            else
                allocate(sbrks, source=pretty_breaks(minval((x - mu) / sd), maxval((x - mu) / sd), n))
            end if
            result%breaks = sbrks * sd + mu
            call sort_real(result%breaks)
            call unique_adjacent(result%breaks, tmp)
            result%breaks = tmp
        case ("equal")
            allocate(result%breaks(n + 1))
            do i = 0, n
                result%breaks(i + 1) = minval(x) + real(i, dp) * (maxval(x) - minval(x)) / real(n, dp)
            end do
        case ("pretty")
            result%breaks = pretty_breaks(minval(x), maxval(x), n)
        case ("quantile")
            allocate(result%breaks(n + 1))
            do i = 0, n
                result%breaks(i + 1) = quantile_r(x, real(i, dp) / real(n, dp), opts%quantile_type)
            end do
        case ("kmeans")
            call kmeans_breaks(x, n, rng, opts%kmeans_iter_max, opts%kmeans_nstart, result%breaks)
        case ("hclust")
            call hclust_fit_1d(x, trim(opts%hclust_method), result%hclust)
            call hclust_cut_breaks(x, result%hclust, n, result%breaks)
            result%has_hclust = .true.
        case ("bclust")
            call bclust_fit_breaks(x, n, rng, opts, result%bclust, result%breaks)
            result%has_bclust = .true.
        case ("fisher", "jenks")
            work_x = x
            if (opts%sample_large_fisher_jenks .and. size(ux) > opts%large_n) then
                nsamp = ceiling(opts%sample_proportion * real(size(ux), dp))
                nsamp = min(opts%large_n, max(n, nsamp))
                nsamp = min(nsamp, size(x))
                call sample_without_replacement(x, nsamp, rng, sampled)
                allocate(work_x(nsamp + 2))
                work_x(1) = minval(x)
                work_x(2) = maxval(x)
                work_x(3:) = sampled
                result%sampled = .true.
            end if
            if (use_style == "fisher") then
                call fisher_exact(work_x, n, result%fisher_stats, result%breaks)
            else
                result%breaks = jenks_breaks(work_x, n)
            end if
        case ("dpih")
            if (opts%dpih_has_range) then
                range_x = opts%dpih_range
            else
                range_x = [minval(x), maxval(x)]
            end if
            h = dpih_bandwidth(x, trim(opts%dpih_scale), opts%dpih_level, opts%dpih_gridsize, range_x, opts%dpih_truncate)
            call sequence_by(range_x(1), range_x(2), h, result%breaks)
        case ("headtails")
            call headtails_breaks(x, opts%headtails_threshold, result%breaks)
        case ("maximum")
            call maximum_breaks(x, n, result%breaks)
        case ("box")
            call box_breaks(x, opts%box_iqr_mult, opts%box_quantile_type, opts%box_legacy, result%breaks)
        case default
            error stop "classint_fit: unknown style"
        end select
        result%style = use_style
    end subroutine classint_fit_impl

    subroutine unique_breaks(unique_x, breaks)
        real(dp), intent(in) :: unique_x(:) !! Sorted or unsorted distinct finite values, each assigned its own class.
        real(dp), allocatable, intent(out) :: breaks(:) !! Midpoint breaks extending half a mean gap beyond the data range.
        real(dp), allocatable :: s(:)
        real(dp) :: mean_gap
        integer :: i

        s = unique_x
        call sort_real(s)
        mean_gap = sum(s(2:) - s(:size(s) - 1)) / real(size(s) - 1, dp)
        allocate(breaks(size(s) + 1))
        breaks(1) = s(1) - 0.5_dp * mean_gap
        do i = 1, size(s) - 1
            breaks(i + 1) = 0.5_dp * (s(i) + s(i + 1))
        end do
        breaks(size(s) + 1) = s(size(s)) + 0.5_dp * mean_gap
    end subroutine unique_breaks

    subroutine unique_adjacent(x, y)
        real(dp), intent(in) :: x(:) !! Sorted real vector from which exactly repeated adjacent values are removed.
        real(dp), allocatable, intent(out) :: y(:) !! Sorted copy containing one representative of each exact value.
        integer :: i
        integer :: n

        if (size(x) == 0) then
            allocate(y(0))
            return
        end if
        n = 1
        do i = 2, size(x)
            if (x(i) < x(i - 1) .or. x(i) > x(i - 1)) n = n + 1
        end do
        allocate(y(n))
        y(1) = x(1)
        n = 1
        do i = 2, size(x)
            if (x(i) < x(i - 1) .or. x(i) > x(i - 1)) then
                n = n + 1
                y(n) = x(i)
            end if
        end do
    end subroutine unique_adjacent

    subroutine sequence_by(lo, hi, step, values)
        real(dp), intent(in) :: lo !! First value of an R seq-style increasing sequence.
        real(dp), intent(in) :: hi !! Upper bound not intentionally exceeded by the generated sequence.
        real(dp), intent(in) :: step !! Positive increment between consecutive values.
        real(dp), allocatable, intent(out) :: values(:) !! Sequence from lo by step through the last value not above hi.
        integer :: n
        integer :: i

        if (step <= 0.0_dp) error stop "sequence_by: step must be positive"
        n = floor((hi - lo) / step + 1.0e-12_dp) + 1
        n = max(1, n)
        allocate(values(n))
        do i = 1, n
            values(i) = lo + real(i - 1, dp) * step
        end do
    end subroutine sequence_by

    subroutine headtails_breaks(x, threshold, breaks)
        real(dp), intent(in) :: x(:) !! Finite observations partitioned recursively above successive head means.
        real(dp), intent(in) :: threshold !! Maximum head proportion for another iteration; clipped to the range [0,1].
        real(dp), allocatable, intent(out) :: breaks(:) !! Sorted head/tail breaks including the global minimum and maximum.
        real(dp), allocatable :: head(:)
        real(dp), allocatable :: next(:)
        real(dp), allocatable :: raw(:)
        real(dp), allocatable :: sorted(:)
        real(dp) :: thr
        real(dp) :: mu
        real(dp) :: prop
        integer :: nraw
        integer :: iter

        thr = min(1.0_dp, max(0.0_dp, threshold))
        head = x
        allocate(raw(102))
        nraw = 1
        raw(1) = minval(x)
        do iter = 1, 100
            mu = mean_dp(head)
            nraw = nraw + 1
            raw(nraw) = mu
            next = pack(head, head > mu)
            prop = real(size(next), dp) / real(size(head), dp)
            head = next
            if (.not. (prop <= thr .and. size(head) > 1)) exit
        end do
        nraw = nraw + 1
        raw(nraw) = maxval(x)
        sorted = raw(:nraw)
        call sort_real(sorted)
        call unique_adjacent(sorted, breaks)
    end subroutine headtails_breaks

    subroutine maximum_breaks(x, n, breaks)
        real(dp), intent(in) :: x(:) !! Finite observations whose largest adjacent sorted gaps determine interval boundaries.
        integer, intent(in) :: n !! Nominal class count; ties among the largest gaps can create extra classes.
        real(dp), allocatable, intent(out) :: breaks(:) !! Sorted endpoints and midpoints at every selected maximum gap.
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: diffs(:)
        real(dp), allocatable :: ranked(:)
        real(dp) :: cutoff
        integer :: i
        integer :: count_selected
        integer :: pos

        s = x
        call sort_real(s)
        allocate(diffs(size(s) - 1))
        diffs = s(2:) - s(:size(s) - 1)
        ranked = diffs
        call sort_real(ranked)
        cutoff = ranked(max(1, size(ranked) - min(n - 1, size(ranked)) + 1))
        count_selected = count(diffs >= cutoff)
        allocate(breaks(count_selected + 2))
        breaks(1) = s(1)
        pos = 1
        do i = 1, size(diffs)
            if (diffs(i) >= cutoff) then
                pos = pos + 1
                breaks(pos) = 0.5_dp * (s(i) + s(i + 1))
            end if
        end do
        breaks(size(breaks)) = s(size(s))
    end subroutine maximum_breaks

    subroutine box_breaks(x, iqr_mult, qtype, legacy, breaks)
        real(dp), intent(in) :: x(:) !! Finite observations used to compute quartiles, fences, and box-map intervals.
        real(dp), intent(in) :: iqr_mult !! Nonnegative multiplier applied to the interquartile range; classInt defaults to 1.5.
        integer, intent(in) :: qtype !! R quantile algorithm number used for the five-number summary.
        logical, intent(in) :: legacy !! If true, reproduce pre-0.4-9 finite fence/end behavior instead of +/-Inf outer breaks.
        real(dp), allocatable, intent(out) :: breaks(:) !! Seven box-map break values defining six classes.
        real(dp) :: qv(5)
        real(dp) :: iqr
        real(dp) :: lowfence
        real(dp) :: upfence
        integer :: i

        if (iqr_mult < 0.0_dp) error stop "box_breaks: iqr_mult must be nonnegative"
        do i = 1, 5
            qv(i) = quantile_r(x, real(i - 1, dp) / 4.0_dp, qtype)
        end do
        iqr = iqr_mult * (qv(4) - qv(2))
        lowfence = qv(2) - iqr
        upfence = qv(4) + iqr
        allocate(breaks(7))
        if (lowfence < qv(1)) then
            if (legacy) then
                breaks(1) = lowfence
                breaks(2) = floor(qv(1))
            else
                breaks(1) = ieee_value(0.0_dp, ieee_negative_inf)
                breaks(2) = lowfence
            end if
        else
            breaks(1) = qv(1)
            breaks(2) = lowfence
        end if
        breaks(3:5) = qv(2:4)
        if (upfence > qv(5)) then
            if (legacy) then
                breaks(6) = ceiling(qv(5))
                breaks(7) = upfence
            else
                breaks(6) = upfence
                breaks(7) = ieee_value(0.0_dp, ieee_positive_inf)
            end if
        else
            breaks(6) = upfence
            breaks(7) = qv(5)
        end if
    end subroutine box_breaks

    subroutine get_hclust_class_intervals(fit, k, result)
        type(class_intervals), intent(in) :: fit !! Existing classInt fit produced with style="hclust" and retaining its hierarchy.
        integer, intent(in) :: k !! New hierarchical cut count from one through the number of finite observations.
        type(class_intervals), intent(out) :: result !! Interval copy with breaks recomputed at the requested hclust cut.
        real(dp), allocatable :: x(:)

        if (.not. fit%has_hclust) error stop "get_hclust_class_intervals: fit has no hclust model"
        result = fit
        x = finite_values(fit%values)
        call hclust_cut_breaks(x, fit%hclust, k, result%breaks)
    end subroutine get_hclust_class_intervals

    ! This remains impure until the translated e1071 bagged-clustering helpers expose pure interfaces.
    subroutine get_bclust_class_intervals(fit, k, result)
        !! Recomputes a retained bagged-clustering fit for a requested number of classes.
        type(class_intervals), intent(in) :: fit !! Existing classInt fit produced with style="bclust" and retaining its committee.
        integer, intent(in) :: k !! New bagged-hierarchy cut count supported by the retained e1071 bclust model.
        type(class_intervals), intent(out) :: result !! Interval copy with breaks recomputed at the requested bclust cut.
        real(dp), allocatable :: x(:)

        if (.not. fit%has_bclust) error stop "get_bclust_class_intervals: fit has no bclust model"
        result = fit
        x = finite_values(fit%values)
        call bclust_cut_breaks(x, fit%bclust, k, result%breaks)
    end subroutine get_bclust_class_intervals
end module classint_core
