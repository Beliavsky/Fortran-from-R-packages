! SPDX-License-Identifier: GPL-2.0-or-later
! The adaptive rejection Metropolis sampling core follows the upstream C port
! in dlm/src/arms-R.c, which credits the original ARMS code to Wally Gilks.
module dlm_arms
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use dlm_types, only : dp, dlm_success, dlm_invalid_shape, dlm_invalid_argument
    use dlm_types, only : dlm_sampling_failure
    use dlm_bounds, only : dlm_indicator, convex_bounds
    use dlm_random, only : seed_dlm_rng
    implicit none
    private

    real(dp), parameter :: xeps = 1.0e-5_dp
    real(dp), parameter :: yeps = 0.1_dp
    real(dp), parameter :: eyeps = 0.001_dp
    real(dp), parameter :: yceil = 50.0_dp
    integer, parameter :: ninit_default = 4
    integer, parameter :: npoint_default = 100

    type, abstract, public :: dlm_log_density
    contains
        procedure(dlm_log_density_evaluate), deferred, pass :: evaluate
    end type dlm_log_density

    abstract interface
        pure real(dp) function dlm_log_density_evaluate(self, x) result(log_density)
            import :: dlm_log_density, dp
            class(dlm_log_density), intent(in) :: self !! Density object holding any fixed target-distribution parameters.
            real(dp), intent(in) :: x(:) !! Point at which the natural logarithm of the unnormalized density is evaluated.
        end function dlm_log_density_evaluate
    end interface

    type :: arms_point
        real(dp) :: x = 0.0_dp
        real(dp) :: y = 0.0_dp
        real(dp) :: ey = 0.0_dp
        real(dp) :: cum = 0.0_dp
        logical :: evaluated = .false.
        integer :: left = 0
        integer :: right = 0
    end type arms_point

    type :: arms_envelope
        integer :: current_points = 0
        integer :: maximum_points = 0
        integer :: evaluations = 0
        real(dp) :: ymax = 0.0_dp
        real(dp) :: convexity = 1.0_dp
        type(arms_point), allocatable :: point(:)
    end type arms_envelope

    type :: arms_metropolis
        logical :: enabled = .true.
        real(dp) :: previous_x = 0.0_dp
        real(dp) :: previous_y = 0.0_dp
    end type arms_metropolis

    public :: arms

contains

    subroutine arms(y_start, density, indicator, n_sample, draws, info, seed, max_points)
        real(dp), intent(in) :: y_start(:) !! Starting point inside the bounded convex support.
        class(dlm_log_density), intent(in) :: density !! Log-density object evaluated by the ARMS rejection/Metropolis steps.
        class(dlm_indicator), intent(in) :: indicator !! Bounded convex support indicator used to find line-search limits.
        integer, intent(in) :: n_sample !! Number of retained samples; zero returns an empty draw matrix.
        real(dp), allocatable, intent(out) :: draws(:, :) !! Samples by row, with one column per component of `y_start`.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status when support/envelope construction fails.
        integer, intent(in), optional :: seed !! Optional deterministic scalar seed for the Fortran intrinsic RNG stream.
        integer, intent(in), optional :: max_points !! Maximum univariate envelope points; default 100 and minimum nine.

        real(dp), allocatable :: direction(:), line_draw(:), current(:), base(:)
        real(dp) :: bounds_offset(2), bounds_actual(2)
        integer :: dimension, k, point_limit, status

        info = dlm_success
        dimension = size(y_start)
        point_limit = npoint_default
        if (present(max_points)) point_limit = max_points
        if (dimension < 1 .or. n_sample < 0 .or. point_limit < 2 * ninit_default + 1) then
            allocate(draws(0, max(0, dimension)))
            info = dlm_invalid_argument
            return
        end if
        if (.not. indicator%evaluate(y_start)) then
            allocate(draws(0, dimension))
            info = dlm_invalid_argument
            return
        end if
        if (present(seed)) call seed_dlm_rng(seed)
        allocate(draws(n_sample, dimension))
        if (n_sample == 0) return

        allocate(direction(dimension), current(dimension), base(dimension))
        current = y_start
        if (dimension == 1) then
            direction = 1.0_dp
            call convex_bounds(current, direction, indicator, bounds_offset, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            bounds_actual = current(1) + bounds_offset
            if (bounds_actual(2) - bounds_actual(1) < 1.0e-7_dp) then
                draws(:, 1) = current(1)
                return
            end if
            base = 0.0_dp
            call arms_line(bounds_actual, density, base, direction, current(1), n_sample, &
                           line_draw, status, point_limit)
            if (status /= dlm_success) then
                info = status
                return
            end if
            draws(:, 1) = line_draw
            return
        end if

        do k = 1, n_sample
            call normal_direction(direction)
            call convex_bounds(current, direction, indicator, bounds_offset, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            if (bounds_offset(2) - bounds_offset(1) < 1.0e-7_dp) then
                draws(k, :) = current
                cycle
            end if
            base = current
            call arms_line(bounds_offset, density, base, direction, 0.0_dp, 1, &
                           line_draw, status, point_limit)
            if (status /= dlm_success) then
                info = status
                return
            end if
            current = current + line_draw(1) * direction
            draws(k, :) = current
        end do
    end subroutine arms

    subroutine arms_line(bounds, density, base, direction, previous, n_sample, sample_values, info, max_points)
        real(dp), intent(in) :: bounds(2) !! Lower and upper scalar limits for the current one-dimensional line target.
        class(dlm_log_density), intent(in) :: density !! Full-space log density evaluated along `base + x*direction`.
        real(dp), intent(in) :: base(:) !! Fixed full-space point defining the origin of the sampled line.
        real(dp), intent(in) :: direction(:) !! Full-space direction multiplying the scalar ARMS coordinate.
        real(dp), intent(in) :: previous !! Previous scalar chain position used by the Metropolis correction.
        integer, intent(in) :: n_sample !! Number of scalar ARMS draws requested from this line target.
        real(dp), allocatable, intent(out) :: sample_values(:) !! Retained scalar line coordinates.
        integer, intent(out) :: info !! Zero on success or `dlm_sampling_failure`/input status on failure.
        integer, intent(in) :: max_points !! Maximum number of points retained in the adaptive rejection envelope.

        type(arms_envelope) :: envelope
        type(arms_metropolis) :: metropolis
        type(arms_point) :: work
        real(dp) :: x_initial(ninit_default)
        integer :: accepted, i, status

        info = dlm_success
        allocate(sample_values(n_sample))
        if (size(base) /= size(direction) .or. size(base) < 1 .or. n_sample < 1) then
            info = dlm_invalid_shape
            return
        end if
        if (.not. bounds(1) < bounds(2) .or. previous < bounds(1) .or. previous > bounds(2)) then
            info = dlm_invalid_argument
            return
        end if
        do i = 1, ninit_default
            x_initial(i) = bounds(1) + real(i, dp) * (bounds(2) - bounds(1)) / real(ninit_default + 1, dp)
        end do

        envelope%convexity = 1.0_dp
        metropolis%enabled = .true.
        call initial_envelope(x_initial, bounds(1), bounds(2), max_points, density, &
                              base, direction, envelope, metropolis, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        metropolis%previous_x = previous
        call evaluate_line(density, base, direction, previous, envelope, metropolis%previous_y, status)
        if (status /= dlm_success) then
            info = status
            return
        end if

        accepted = 0
        do while (accepted < n_sample)
            call sample_envelope(envelope, work, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            call test_point(envelope, work, density, base, direction, metropolis, status)
            select case (status)
            case (1)
                accepted = accepted + 1
                sample_values(accepted) = work%x
            case (0)
                cycle
            case default
                info = dlm_sampling_failure
                return
            end select
        end do
    end subroutine arms_line

    subroutine initial_envelope(x_initial, lower, upper, max_points, density, base, direction, envelope, metropolis, info)
        real(dp), intent(in) :: x_initial(:) !! Ordered interior scalar points used to initialize the rejection envelope.
        real(dp), intent(in) :: lower !! Strict lower support bound for the line target.
        real(dp), intent(in) :: upper !! Strict upper support bound for the line target.
        integer, intent(in) :: max_points !! Maximum adaptive envelope storage capacity.
        class(dlm_log_density), intent(in) :: density !! Full-space log density evaluated at initial line points.
        real(dp), intent(in) :: base(:) !! Full-space line origin used for density evaluations.
        real(dp), intent(in) :: direction(:) !! Full-space line direction used for density evaluations.
        type(arms_envelope), intent(out) :: envelope !! Initialized piecewise-linear log-envelope and cumulative integral.
        type(arms_metropolis), intent(in) :: metropolis !! Metropolis configuration controlling convexity correction behavior.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status when envelope initialization fails.

        integer :: i, index, mpoint, ninit, status

        info = dlm_success
        envelope%convexity = 1.0_dp
        ninit = size(x_initial)
        mpoint = 2 * ninit + 1
        if (ninit < 3 .or. max_points < mpoint .or. lower >= upper) then
            info = dlm_invalid_argument
            return
        end if
        if (x_initial(1) <= lower .or. x_initial(ninit) >= upper) then
            info = dlm_invalid_argument
            return
        end if
        do i = 2, ninit
            if (x_initial(i) <= x_initial(i - 1)) then
                info = dlm_invalid_argument
                return
            end if
        end do
        if (envelope%convexity < 0.0_dp) then
            info = dlm_invalid_argument
            return
        end if

        envelope%maximum_points = max_points
        envelope%evaluations = 0
        allocate(envelope%point(max_points))
        envelope%point(1)%x = lower
        envelope%point(1)%evaluated = .false.
        envelope%point(1)%left = 0
        envelope%point(1)%right = 2

        do index = 2, mpoint - 1
            envelope%point(index)%left = index - 1
            envelope%point(index)%right = index + 1
            if (mod(index, 2) == 0) then
                i = index / 2
                envelope%point(index)%x = x_initial(i)
                call evaluate_line(density, base, direction, envelope%point(index)%x, envelope, &
                                   envelope%point(index)%y, status)
                if (status /= dlm_success) then
                    info = status
                    return
                end if
                envelope%point(index)%evaluated = .true.
            else
                envelope%point(index)%evaluated = .false.
            end if
        end do

        envelope%point(mpoint)%x = upper
        envelope%point(mpoint)%evaluated = .false.
        envelope%point(mpoint)%left = mpoint - 1
        envelope%point(mpoint)%right = 0
        do index = 1, mpoint, 2
            call meet_point(index, envelope, metropolis, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
        end do
        call cumulate_envelope(envelope, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        envelope%current_points = mpoint
    end subroutine initial_envelope

    subroutine sample_envelope(envelope, work, info)
        type(arms_envelope), intent(in) :: envelope !! Current adaptive rejection envelope from which a proposal is sampled.
        type(arms_point), intent(out) :: work !! Proposed scalar coordinate and envelope height with neighboring point indices.
        integer, intent(out) :: info !! Zero on success or `dlm_sampling_failure` for an invalid cumulative envelope.

        real(dp) :: probability

        call random_number(probability)
        call invert_envelope(probability, envelope, work, info)
    end subroutine sample_envelope

    pure subroutine invert_envelope(probability, envelope, work, info)
        real(dp), intent(in) :: probability !! Uniform probability in `[0,1)` used to invert the envelope CDF.
        type(arms_envelope), intent(in) :: envelope !! Piecewise exponential rejection envelope and cumulative areas.
        type(arms_point), intent(out) :: work !! Scalar proposal obtained by inverting the envelope integral.
        integer, intent(out) :: info !! Zero on success or `dlm_sampling_failure` when envelope geometry is invalid.

        real(dp) :: e_left, e_right, proportion, u, x_left, x_right, y_left, y_right
        integer :: left_index, q

        info = dlm_success
        q = 1
        do while (envelope%point(q)%right /= 0)
            q = envelope%point(q)%right
        end do
        if (envelope%point(q)%cum <= 0.0_dp) then
            info = dlm_sampling_failure
            return
        end if
        u = min(max(probability, 0.0_dp), 1.0_dp) * envelope%point(q)%cum
        do
            left_index = envelope%point(q)%left
            if (left_index == 0) then
                info = dlm_sampling_failure
                return
            end if
            if (envelope%point(left_index)%cum <= u) exit
            q = left_index
        end do

        left_index = envelope%point(q)%left
        work%left = left_index
        work%right = q
        work%evaluated = .false.
        work%cum = u
        if (envelope%point(q)%cum <= envelope%point(left_index)%cum) then
            info = dlm_sampling_failure
            return
        end if
        proportion = (u - envelope%point(left_index)%cum) / &
                     (envelope%point(q)%cum - envelope%point(left_index)%cum)
        proportion = min(1.0_dp, max(0.0_dp, proportion))

        x_left = envelope%point(left_index)%x
        x_right = envelope%point(q)%x
        y_left = envelope%point(left_index)%y
        y_right = envelope%point(q)%y
        e_left = envelope%point(left_index)%ey
        e_right = envelope%point(q)%ey
        if (same_coordinate(x_left, x_right)) then
            work%x = x_right
            work%y = y_right
            work%ey = e_right
        else if (abs(y_right - y_left) < yeps) then
            if (abs(e_right - e_left) > eyeps * abs(e_right + e_left)) then
                work%x = x_left + (x_right - x_left) / (e_right - e_left) * &
                         (-e_left + sqrt(max(0.0_dp, (1.0_dp - proportion) * e_left * e_left + &
                                                   proportion * e_right * e_right)))
            else
                work%x = x_left + (x_right - x_left) * proportion
            end if
            work%ey = (work%x - x_left) / (x_right - x_left) * (e_right - e_left) + e_left
            work%y = log_shift(work%ey, envelope%ymax)
        else
            work%x = x_left + (x_right - x_left) / (y_right - y_left) * &
                     (-y_left + log_shift((1.0_dp - proportion) * e_left + proportion * e_right, &
                                          envelope%ymax))
            work%y = (work%x - x_left) / (x_right - x_left) * (y_right - y_left) + y_left
            work%ey = exp_shift(work%y, envelope%ymax)
        end if
        if (work%x < min(x_left, x_right) - 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x_left)) .or. &
            work%x > max(x_left, x_right) + 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x_right))) then
            info = dlm_sampling_failure
        end if
    end subroutine invert_envelope

    subroutine test_point(envelope, work, density, base, direction, metropolis, decision)
        type(arms_envelope), intent(inout) :: envelope !! Adaptive rejection envelope updated after evaluated proposals.
        type(arms_point), intent(inout) :: work !! Proposed point, replaced by the previous state after Metropolis rejection.
        class(dlm_log_density), intent(in) :: density !! Full-space log density used for proposal evaluation.
        real(dp), intent(in) :: base(:) !! Full-space origin defining the currently sampled line.
        real(dp), intent(in) :: direction(:) !! Full-space direction defining the currently sampled line.
        type(arms_metropolis), intent(inout) :: metropolis !! Previous chain state and Metropolis enable flag.
        integer, intent(out) :: decision !! One for an accepted sample, zero for rejection, or negative on envelope failure.

        real(dp) :: acceptance, envelope_old, envelope_new, log_reject, log_squeeze
        real(dp) :: old_density, random_value, weight, y_new
        integer :: left_eval, q_left, q_right, right_eval, status

        decision = 0
        call random_number(random_value)
        log_reject = log_shift(random_value * work%ey, envelope%ymax)

        if (.not. metropolis%enabled .and. envelope%point(work%left)%left /= 0 .and. &
            envelope%point(work%right)%right /= 0) then
            if (envelope%point(work%left)%evaluated) then
                left_eval = work%left
            else
                left_eval = envelope%point(work%left)%left
            end if
            if (envelope%point(work%right)%evaluated) then
                right_eval = work%right
            else
                right_eval = envelope%point(work%right)%right
            end if
            if (right_eval /= 0 .and. left_eval /= 0) then
                if (.not. same_coordinate(envelope%point(right_eval)%x, envelope%point(left_eval)%x)) then
                    log_squeeze = (envelope%point(right_eval)%y * (work%x - envelope%point(left_eval)%x) + &
                                   envelope%point(left_eval)%y * (envelope%point(right_eval)%x - work%x)) / &
                                  (envelope%point(right_eval)%x - envelope%point(left_eval)%x)
                    if (log_reject <= log_squeeze) then
                        decision = 1
                        return
                    end if
                end if
            end if
        end if

        call evaluate_line(density, base, direction, work%x, envelope, y_new, status)
        if (status /= dlm_success) then
            decision = -1
            return
        end if
        if (.not. metropolis%enabled .or. log_reject >= y_new) then
            work%y = y_new
            work%ey = exp_shift(work%y, envelope%ymax)
            work%evaluated = .true.
            call update_envelope(envelope, work, density, base, direction, metropolis, status)
            if (status /= dlm_success) then
                decision = -1
                return
            end if
            if (log_reject >= y_new) then
                decision = 0
            else
                decision = 1
            end if
            return
        end if

        old_density = metropolis%previous_y
        q_left = 1
        do while (envelope%point(q_left)%left /= 0)
            q_left = envelope%point(q_left)%left
        end do
        do while (envelope%point(q_left)%right /= 0)
            q_right = envelope%point(q_left)%right
            if (envelope%point(q_right)%x >= metropolis%previous_x) exit
            q_left = q_right
        end do
        q_right = envelope%point(q_left)%right
        if (q_right == 0) then
            decision = -1
            return
        end if
        if (same_coordinate(envelope%point(q_right)%x, envelope%point(q_left)%x)) then
            decision = -1
            return
        end if
        weight = (metropolis%previous_x - envelope%point(q_left)%x) / &
                 (envelope%point(q_right)%x - envelope%point(q_left)%x)
        envelope_old = envelope%point(q_left)%y + weight * &
                       (envelope%point(q_right)%y - envelope%point(q_left)%y)
        envelope_new = work%y
        envelope_old = min(envelope_old, old_density)
        envelope_new = min(envelope_new, y_new)
        acceptance = y_new - envelope_new - old_density + envelope_old
        acceptance = min(0.0_dp, acceptance)
        if (acceptance > -yceil) then
            acceptance = exp(acceptance)
        else
            acceptance = 0.0_dp
        end if
        call random_number(random_value)
        if (random_value > acceptance) then
            work%x = metropolis%previous_x
            work%y = metropolis%previous_y
            work%ey = exp_shift(work%y, envelope%ymax)
            work%evaluated = .true.
            work%left = q_left
            work%right = q_right
        else
            metropolis%previous_x = work%x
            metropolis%previous_y = y_new
        end if
        decision = 1
    end subroutine test_point

    subroutine update_envelope(envelope, work, density, base, direction, metropolis, info)
        type(arms_envelope), intent(inout) :: envelope !! Adaptive envelope receiving the newly evaluated density point.
        type(arms_point), intent(in) :: work !! Evaluated proposal and indices of the envelope piece containing it.
        class(dlm_log_density), intent(in) :: density !! Full-space log density used if an inserted point must be moved.
        real(dp), intent(in) :: base(:) !! Full-space origin defining the sampled line.
        real(dp), intent(in) :: direction(:) !! Full-space direction defining the sampled line.
        type(arms_metropolis), intent(in) :: metropolis !! Metropolis configuration controlling convexity correction.
        integer, intent(out) :: info !! Zero on success or `dlm_sampling_failure` for invalid envelope geometry.

        integer :: far_left, far_right, left_neighbor, middle, q, right_neighbor, status
        real(dp) :: left_x, right_x

        info = dlm_success
        if (.not. work%evaluated .or. envelope%current_points > envelope%maximum_points - 2) return
        q = envelope%current_points + 1
        envelope%current_points = q
        envelope%point(q) = work
        envelope%point(q)%evaluated = .true.
        middle = envelope%current_points + 1
        envelope%current_points = middle
        envelope%point(middle) = arms_point()
        envelope%point(middle)%evaluated = .false.

        if (envelope%point(work%left)%evaluated .and. .not. envelope%point(work%right)%evaluated) then
            envelope%point(middle)%left = work%left
            envelope%point(middle)%right = q
            envelope%point(q)%left = middle
            envelope%point(q)%right = work%right
            envelope%point(work%left)%right = middle
            envelope%point(work%right)%left = q
        else if (.not. envelope%point(work%left)%evaluated .and. envelope%point(work%right)%evaluated) then
            envelope%point(middle)%right = work%right
            envelope%point(middle)%left = q
            envelope%point(q)%right = middle
            envelope%point(q)%left = work%left
            envelope%point(work%right)%left = middle
            envelope%point(work%left)%right = q
        else
            info = dlm_sampling_failure
            return
        end if

        if (envelope%point(envelope%point(q)%left)%left /= 0) then
            left_neighbor = envelope%point(envelope%point(q)%left)%left
        else
            left_neighbor = envelope%point(q)%left
        end if
        if (envelope%point(envelope%point(q)%right)%right /= 0) then
            right_neighbor = envelope%point(envelope%point(q)%right)%right
        else
            right_neighbor = envelope%point(q)%right
        end if
        left_x = envelope%point(left_neighbor)%x
        right_x = envelope%point(right_neighbor)%x
        if (envelope%point(q)%x < (1.0_dp - xeps) * left_x + xeps * right_x) then
            envelope%point(q)%x = (1.0_dp - xeps) * left_x + xeps * right_x
            call evaluate_line(density, base, direction, envelope%point(q)%x, envelope, envelope%point(q)%y, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
        else if (envelope%point(q)%x > xeps * left_x + (1.0_dp - xeps) * right_x) then
            envelope%point(q)%x = xeps * left_x + (1.0_dp - xeps) * right_x
            call evaluate_line(density, base, direction, envelope%point(q)%x, envelope, envelope%point(q)%y, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
        end if

        call meet_point(envelope%point(q)%left, envelope, metropolis, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        call meet_point(envelope%point(q)%right, envelope, metropolis, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        left_neighbor = envelope%point(q)%left
        far_left = envelope%point(left_neighbor)%left
        if (far_left /= 0) then
            far_left = envelope%point(far_left)%left
            if (far_left /= 0) then
                call meet_point(far_left, envelope, metropolis, status)
                if (status /= dlm_success) then
                    info = status
                    return
                end if
            end if
        end if
        right_neighbor = envelope%point(q)%right
        far_right = envelope%point(right_neighbor)%right
        if (far_right /= 0) then
            far_right = envelope%point(far_right)%right
            if (far_right /= 0) then
                call meet_point(far_right, envelope, metropolis, status)
                if (status /= dlm_success) then
                    info = status
                    return
                end if
            end if
        end if
        call cumulate_envelope(envelope, info)
    end subroutine update_envelope

    pure subroutine cumulate_envelope(envelope, info)
        type(arms_envelope), intent(inout) :: envelope !! Envelope whose shifted exponentials and cumulative areas are recomputed.
        integer, intent(out) :: info !! Zero on success or `dlm_sampling_failure` when an area is negative or nonfinite.

        real(dp) :: segment_area
        integer :: leftmost, q

        info = dlm_success
        leftmost = 1
        do while (envelope%point(leftmost)%left /= 0)
            leftmost = envelope%point(leftmost)%left
        end do
        envelope%ymax = envelope%point(leftmost)%y
        q = envelope%point(leftmost)%right
        do while (q /= 0)
            envelope%ymax = max(envelope%ymax, envelope%point(q)%y)
            q = envelope%point(q)%right
        end do
        q = leftmost
        do while (q /= 0)
            envelope%point(q)%ey = exp_shift(envelope%point(q)%y, envelope%ymax)
            q = envelope%point(q)%right
        end do
        envelope%point(leftmost)%cum = 0.0_dp
        q = envelope%point(leftmost)%right
        do while (q /= 0)
            segment_area = envelope_area(q, envelope)
            if (segment_area < 0.0_dp .or. .not. ieee_is_finite(segment_area)) then
                info = dlm_sampling_failure
                return
            end if
            envelope%point(q)%cum = envelope%point(envelope%point(q)%left)%cum + segment_area
            q = envelope%point(q)%right
        end do
    end subroutine cumulate_envelope

    pure subroutine meet_point(q, envelope, metropolis, info)
        integer, intent(in) :: q !! Index of an unevaluated envelope intersection or boundary point to update.
        type(arms_envelope), intent(inout) :: envelope !! Envelope updated with the computed intersection.
        type(arms_metropolis), intent(in) :: metropolis !! Metropolis flag permitting convexity-gradient adjustment.
        integer, intent(out) :: info !! Zero on success or `dlm_sampling_failure` for an invalid/violating envelope.

        real(dp) :: d_left, d_right, gradient_left, gradient_right, gradient_span
        integer :: left1, left2, left3, right1, right2, right3
        logical :: have_left, have_right, have_span

        info = dlm_success
        if (q < 1 .or. q > size(envelope%point) .or. envelope%point(q)%evaluated) then
            info = dlm_sampling_failure
            return
        end if
        left1 = envelope%point(q)%left
        left2 = 0
        left3 = 0
        if (left1 /= 0) left2 = envelope%point(left1)%left
        if (left2 /= 0) left3 = envelope%point(left2)%left
        have_left = left3 /= 0
        if (have_left) then
            if (same_coordinate(envelope%point(left1)%x, envelope%point(left3)%x)) then
                info = dlm_sampling_failure
                return
            end if
            gradient_left = (envelope%point(left1)%y - envelope%point(left3)%y) / &
                            (envelope%point(left1)%x - envelope%point(left3)%x)
        else
            gradient_left = 0.0_dp
        end if

        right1 = envelope%point(q)%right
        right2 = 0
        right3 = 0
        if (right1 /= 0) right2 = envelope%point(right1)%right
        if (right2 /= 0) right3 = envelope%point(right2)%right
        have_right = right3 /= 0
        if (have_right) then
            if (same_coordinate(envelope%point(right1)%x, envelope%point(right3)%x)) then
                info = dlm_sampling_failure
                return
            end if
            gradient_right = (envelope%point(right1)%y - envelope%point(right3)%y) / &
                             (envelope%point(right1)%x - envelope%point(right3)%x)
        else
            gradient_right = 0.0_dp
        end if

        have_span = left1 /= 0 .and. right1 /= 0
        if (have_span) then
            if (same_coordinate(envelope%point(right1)%x, envelope%point(left1)%x)) then
                info = dlm_sampling_failure
                return
            end if
            gradient_span = (envelope%point(right1)%y - envelope%point(left1)%y) / &
                            (envelope%point(right1)%x - envelope%point(left1)%x)
        else
            gradient_span = 0.0_dp
        end if

        if (have_span .and. have_left .and. gradient_left < gradient_span) then
            if (.not. metropolis%enabled) then
                info = dlm_sampling_failure
                return
            end if
            gradient_left = gradient_left + (1.0_dp + envelope%convexity) * &
                            (gradient_span - gradient_left)
        end if
        if (have_span .and. have_right .and. gradient_right > gradient_span) then
            if (.not. metropolis%enabled) then
                info = dlm_sampling_failure
                return
            end if
            gradient_right = gradient_right + (1.0_dp + envelope%convexity) * &
                             (gradient_span - gradient_right)
        end if

        d_right = 0.0_dp
        d_left = 0.0_dp
        if (have_left .and. have_span) then
            d_right = (gradient_left - gradient_span) * &
                      (envelope%point(right1)%x - envelope%point(left1)%x)
            d_right = max(d_right, yeps)
        end if
        if (have_right .and. have_span) then
            d_left = (gradient_span - gradient_right) * &
                     (envelope%point(right1)%x - envelope%point(left1)%x)
            d_left = max(d_left, yeps)
        end if

        if (have_left .and. have_right .and. have_span) then
            envelope%point(q)%x = (d_left * envelope%point(right1)%x + &
                                   d_right * envelope%point(left1)%x) / (d_left + d_right)
            envelope%point(q)%y = (d_left * envelope%point(right1)%y + &
                                   d_right * envelope%point(left1)%y + d_left * d_right) / &
                                  (d_left + d_right)
        else if (have_left .and. have_span) then
            envelope%point(q)%x = envelope%point(right1)%x
            envelope%point(q)%y = envelope%point(right1)%y + d_right
        else if (have_right .and. have_span) then
            envelope%point(q)%x = envelope%point(left1)%x
            envelope%point(q)%y = envelope%point(left1)%y + d_left
        else if (have_left) then
            envelope%point(q)%y = envelope%point(left1)%y + gradient_left * &
                                  (envelope%point(q)%x - envelope%point(left1)%x)
        else if (have_right) then
            envelope%point(q)%y = envelope%point(right1)%y - gradient_right * &
                                  (envelope%point(right1)%x - envelope%point(q)%x)
        else
            info = dlm_sampling_failure
            return
        end if
        if (left1 /= 0) then
            if (envelope%point(q)%x < envelope%point(left1)%x) info = dlm_sampling_failure
        end if
        if (right1 /= 0) then
            if (envelope%point(q)%x > envelope%point(right1)%x) info = dlm_sampling_failure
        end if
    end subroutine meet_point

    pure real(dp) function envelope_area(q, envelope) result(area_value)
        integer, intent(in) :: q !! Right-hand point index of the envelope segment being integrated.
        type(arms_envelope), intent(in) :: envelope !! Envelope supplying segment coordinates and shifted exponentials.

        integer :: left

        left = envelope%point(q)%left
        if (left == 0) then
            area_value = 0.0_dp
        else if (same_coordinate(envelope%point(left)%x, envelope%point(q)%x)) then
            area_value = 0.0_dp
        else if (abs(envelope%point(q)%y - envelope%point(left)%y) < yeps) then
            area_value = 0.5_dp * (envelope%point(q)%ey + envelope%point(left)%ey) * &
                         (envelope%point(q)%x - envelope%point(left)%x)
        else
            area_value = (envelope%point(q)%ey - envelope%point(left)%ey) / &
                         (envelope%point(q)%y - envelope%point(left)%y) * &
                         (envelope%point(q)%x - envelope%point(left)%x)
        end if
    end function envelope_area

    pure elemental logical function same_coordinate(a, b) result(equal_coordinate)
        real(dp), intent(in) :: a !! First scalar coordinate being compared for envelope degeneracy.
        real(dp), intent(in) :: b !! Second scalar coordinate being compared for envelope degeneracy.

        equal_coordinate = abs(a - b) <= 10.0_dp * epsilon(1.0_dp) * &
                           max(1.0_dp, abs(a), abs(b))
    end function same_coordinate

    pure elemental real(dp) function exp_shift(y, y0) result(value)
        real(dp), intent(in) :: y !! Unshifted log-envelope height.
        real(dp), intent(in) :: y0 !! Current maximum log-envelope height used as the exponential reference.

        if (y - y0 > -2.0_dp * yceil) then
            value = exp(y - y0 + yceil)
        else
            value = 0.0_dp
        end if
    end function exp_shift

    pure elemental real(dp) function log_shift(y, y0) result(value)
        real(dp), intent(in) :: y !! Shifted exponential height produced by `exp_shift`.
        real(dp), intent(in) :: y0 !! Current maximum log-envelope height used as the exponential reference.

        if (y > 0.0_dp) then
            value = log(y) + y0 - yceil
        else
            value = -huge(1.0_dp)
        end if
    end function log_shift

    subroutine evaluate_line(density, base, direction, line_x, envelope, value, info)
        class(dlm_log_density), intent(in) :: density !! Full-space target log density.
        real(dp), intent(in) :: base(:) !! Full-space line origin.
        real(dp), intent(in) :: direction(:) !! Full-space line direction.
        real(dp), intent(in) :: line_x !! Scalar line coordinate at which the target is evaluated.
        type(arms_envelope), intent(inout) :: envelope !! Envelope whose function-evaluation counter is incremented.
        real(dp), intent(out) :: value !! Log density at `base + line_x*direction`.
        integer, intent(out) :: info !! Zero for a finite target value or `dlm_sampling_failure` otherwise.

        value = density%evaluate(base + line_x * direction)
        envelope%evaluations = envelope%evaluations + 1
        if (ieee_is_finite(value)) then
            info = dlm_success
        else
            info = dlm_sampling_failure
        end if
    end subroutine evaluate_line

    subroutine normal_direction(direction)
        real(dp), intent(out) :: direction(:) !! Independent standard-normal components used as a random hit-and-run direction.

        real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
        real(dp) :: radius, theta, u1, u2
        integer :: i

        i = 1
        do while (i <= size(direction))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            radius = sqrt(-2.0_dp * log(u1))
            theta = two_pi * u2
            direction(i) = radius * cos(theta)
            if (i + 1 <= size(direction)) direction(i + 1) = radius * sin(theta)
            i = i + 2
        end do
    end subroutine normal_direction

end module dlm_arms
