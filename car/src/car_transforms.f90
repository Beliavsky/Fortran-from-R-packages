module car_transforms
    use car_kinds, only: dp
    use car_utils, only: mean_valid
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_quiet_nan, ieee_value
    implicit none
    private

    public :: basic_power
    public :: bc_power
    public :: bcn_power
    public :: bcn_power_inverse
    public :: box_cox_variable
    public :: logit_transform
    public :: yj_power

    interface basic_power
        module procedure basic_power_vector
        module procedure basic_power_matrix
    end interface basic_power

    interface bc_power
        module procedure bc_power_vector
        module procedure bc_power_matrix
    end interface bc_power

    interface yj_power
        module procedure yj_power_vector
        module procedure yj_power_matrix
    end interface yj_power

    interface bcn_power
        module procedure bcn_power_vector
        module procedure bcn_power_matrix
    end interface bcn_power

contains

    pure subroutine basic_power_vector(u, lambda, z, status, gamma)
        real(dp), intent(in) :: u(:) !! Input data vector; NaN values are propagated and finite values must be positive after any shift.
        real(dp), intent(in) :: lambda !! Power parameter; values near zero select the logarithm.
        real(dp), intent(out) :: z(:) !! Transformed vector with the same length as u.
        integer, intent(out) :: status !! Zero on success; nonzero for shape mismatch or nonpositive transformed-domain values.
        real(dp), intent(in), optional :: gamma !! Optional scalar location shift added before applying the power.
        real(dp) :: shifted
        integer :: i

        if (size(z) /= size(u)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        do i = 1, size(u)
            if (ieee_is_nan(u(i))) then
                z(i) = u(i)
                cycle
            end if
            shifted = u(i)
            if (present(gamma)) shifted = shifted + gamma
            if (shifted <= 0.0_dp) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 2
                return
            end if
            if (abs(lambda) <= 1.0e-6_dp) then
                z(i) = log(shifted)
            else
                z(i) = shifted**lambda
            end if
        end do
        status = 0
    end subroutine basic_power_vector

    pure subroutine basic_power_matrix(u, lambda, z, status, gamma)
        real(dp), intent(in) :: u(:, :) !! Input data matrix with variables in columns; NaNs are propagated.
        real(dp), intent(in) :: lambda(:) !! Power parameter for each column of u.
        real(dp), intent(out) :: z(:, :) !! Transformed matrix with the same shape as u.
        integer, intent(out) :: status !! Zero on success; nonzero for incompatible shapes or invalid domains.
        real(dp), intent(in), optional :: gamma(:) !! Optional location shift for each column; length must equal size(u, 2).
        integer :: j
        integer :: local_status

        if (any(shape(z) /= shape(u)) .or. size(lambda) /= size(u, 2)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (present(gamma)) then
            if (size(gamma) /= size(u, 2)) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 1
                return
            end if
        end if

        do j = 1, size(u, 2)
            if (present(gamma)) then
                call basic_power_vector(u(:, j), lambda(j), z(:, j), local_status, gamma(j))
            else
                call basic_power_vector(u(:, j), lambda(j), z(:, j), local_status)
            end if
            if (local_status /= 0) then
                status = local_status
                return
            end if
        end do
        status = 0
    end subroutine basic_power_matrix

    pure subroutine bc_power_vector(u, lambda, z, status, jacobian_adjusted, gamma)
        real(dp), intent(in) :: u(:) !! Input data vector; finite values must be positive after any optional shift.
        real(dp), intent(in) :: lambda !! Box-Cox power parameter; values near zero use log(u).
        real(dp), intent(out) :: z(:) !! Box-Cox transformed vector with the same length as u.
        integer, intent(out) :: status !! Zero on success; nonzero for shape mismatch, domain failure, or empty valid data.
        logical, intent(in), optional :: jacobian_adjusted !! If true, apply the geometric-mean Jacobian normalization used by car::bcPower.
        real(dp), intent(in), optional :: gamma !! Optional scalar shift added to u before transformation.
        real(dp), allocatable :: logs(:)
        real(dp) :: geometric_mean
        real(dp) :: mean_log
        real(dp) :: shifted
        logical :: adjust
        integer :: i
        integer :: n_valid
        integer :: mean_status

        if (size(z) /= size(u)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        adjust = .false.
        if (present(jacobian_adjusted)) adjust = jacobian_adjusted
        allocate(logs(size(u)))
        logs = ieee_value(0.0_dp, ieee_quiet_nan)

        do i = 1, size(u)
            if (ieee_is_nan(u(i))) then
                z(i) = u(i)
                cycle
            end if
            shifted = u(i)
            if (present(gamma)) shifted = shifted + gamma
            if (shifted <= 0.0_dp) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 2
                return
            end if
            logs(i) = log(shifted)
            if (abs(lambda) <= 1.0e-6_dp) then
                z(i) = logs(i)
            else
                z(i) = (shifted**lambda - 1.0_dp) / lambda
            end if
        end do

        if (adjust) then
            call mean_valid(logs, mean_log, n_valid, mean_status)
            if (mean_status /= 0) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 3
                return
            end if
            geometric_mean = exp(mean_log)
            z = z * geometric_mean**(1.0_dp - lambda)
        end if
        status = 0
    end subroutine bc_power_vector

    pure subroutine bc_power_matrix(u, lambda, z, status, jacobian_adjusted, gamma)
        real(dp), intent(in) :: u(:, :) !! Input data matrix with variables in columns and NaNs allowed.
        real(dp), intent(in) :: lambda(:) !! Box-Cox power for each column.
        real(dp), intent(out) :: z(:, :) !! Transformed matrix with the same shape as u.
        integer, intent(out) :: status !! Zero on success; nonzero for incompatible shapes or invalid domains.
        logical, intent(in), optional :: jacobian_adjusted !! If true, normalize each column by its Jacobian factor.
        real(dp), intent(in), optional :: gamma(:) !! Optional location shift for each column.
        logical :: adjust
        integer :: j
        integer :: local_status

        if (any(shape(z) /= shape(u)) .or. size(lambda) /= size(u, 2)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (present(gamma)) then
            if (size(gamma) /= size(u, 2)) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 1
                return
            end if
        end if

        adjust = .false.
        if (present(jacobian_adjusted)) adjust = jacobian_adjusted
        do j = 1, size(u, 2)
            if (present(gamma)) then
                call bc_power_vector(u(:, j), lambda(j), z(:, j), local_status, adjust, gamma(j))
            else
                call bc_power_vector(u(:, j), lambda(j), z(:, j), local_status, adjust)
            end if
            if (local_status /= 0) then
                status = local_status
                return
            end if
        end do
        status = 0
    end subroutine bc_power_matrix

    pure subroutine yj_power_vector(u, lambda, z, status, jacobian_adjusted)
        real(dp), intent(in) :: u(:) !! Input data vector on the full real line; NaNs are propagated.
        real(dp), intent(in) :: lambda !! Yeo-Johnson power parameter.
        real(dp), intent(out) :: z(:) !! Yeo-Johnson transformed vector with the same length as u.
        integer, intent(out) :: status !! Zero on success; nonzero for shape mismatch or empty valid data.
        logical, intent(in), optional :: jacobian_adjusted !! If true, apply car's Jacobian normalization.
        real(dp), allocatable :: jac_logs(:)
        real(dp) :: base
        real(dp) :: mean_jac_log
        logical :: adjust
        integer :: i
        integer :: n_valid
        integer :: mean_status

        if (size(z) /= size(u)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        adjust = .false.
        if (present(jacobian_adjusted)) adjust = jacobian_adjusted
        allocate(jac_logs(size(u)))
        jac_logs = ieee_value(0.0_dp, ieee_quiet_nan)

        do i = 1, size(u)
            if (ieee_is_nan(u(i))) then
                z(i) = u(i)
                cycle
            end if

            if (u(i) >= 0.0_dp) then
                base = u(i) + 1.0_dp
                if (abs(lambda) <= 1.0e-6_dp) then
                    z(i) = log(base)
                else
                    z(i) = (base**lambda - 1.0_dp) / lambda
                end if
                jac_logs(i) = log(base)
            else
                base = -u(i) + 1.0_dp
                if (abs(2.0_dp - lambda) <= 1.0e-6_dp) then
                    z(i) = -log(base)
                else
                    z(i) = -(base**(2.0_dp - lambda) - 1.0_dp) / (2.0_dp - lambda)
                end if
                jac_logs(i) = -log(base)
            end if
        end do

        if (adjust) then
            call mean_valid(jac_logs, mean_jac_log, n_valid, mean_status)
            if (mean_status /= 0) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 2
                return
            end if
            z = z * exp(mean_jac_log)**(1.0_dp - lambda)
        end if
        status = 0
    end subroutine yj_power_vector

    pure subroutine yj_power_matrix(u, lambda, z, status, jacobian_adjusted)
        real(dp), intent(in) :: u(:, :) !! Input data matrix with variables in columns and unrestricted finite values.
        real(dp), intent(in) :: lambda(:) !! Yeo-Johnson power for each column.
        real(dp), intent(out) :: z(:, :) !! Transformed matrix with the same shape as u.
        integer, intent(out) :: status !! Zero on success; nonzero for incompatible shapes or transformation failure.
        logical, intent(in), optional :: jacobian_adjusted !! If true, normalize each column by the Jacobian factor.
        logical :: adjust
        integer :: j
        integer :: local_status

        if (any(shape(z) /= shape(u)) .or. size(lambda) /= size(u, 2)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        adjust = .false.
        if (present(jacobian_adjusted)) adjust = jacobian_adjusted
        do j = 1, size(u, 2)
            call yj_power_vector(u(:, j), lambda(j), z(:, j), local_status, adjust)
            if (local_status /= 0) then
                status = local_status
                return
            end if
        end do
        status = 0
    end subroutine yj_power_matrix

    pure subroutine bcn_power_vector(u, lambda, gamma, z, status, jacobian_adjusted)
        real(dp), intent(in) :: u(:) !! Input data vector; with gamma near zero, finite values must be strictly positive.
        real(dp), intent(in) :: lambda !! BCN power parameter; values near zero use the logarithmic limit.
        real(dp), intent(in) :: gamma !! Nonnegative BCN location parameter.
        real(dp), intent(out) :: z(:) !! BCN-transformed vector with the same length as u.
        integer, intent(out) :: status !! Zero on success; nonzero for shape, gamma, domain, or Jacobian failures.
        logical, intent(in), optional :: jacobian_adjusted !! If true, divide by the car::bcnPower Jacobian normalization.
        real(dp), allocatable :: log_one_plus(:)
        real(dp), allocatable :: log_u_plus_s(:)
        real(dp) :: jacobian
        real(dp) :: mean_one_plus
        real(dp) :: mean_u_plus_s
        real(dp) :: s
        real(dp) :: t
        logical :: adjust
        integer :: i
        integer :: n_valid
        integer :: mean_status

        if (size(z) /= size(u)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (gamma < 0.0_dp) then
            z = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 2
            return
        end if

        adjust = .false.
        if (present(jacobian_adjusted)) adjust = jacobian_adjusted
        allocate(log_one_plus(size(u)), log_u_plus_s(size(u)))
        log_one_plus = ieee_value(0.0_dp, ieee_quiet_nan)
        log_u_plus_s = ieee_value(0.0_dp, ieee_quiet_nan)

        do i = 1, size(u)
            if (ieee_is_nan(u(i))) then
                z(i) = u(i)
                cycle
            end if
            if (abs(gamma) <= 1.0e-10_dp .and. u(i) <= 0.0_dp) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 3
                return
            end if

            s = sqrt(u(i) * u(i) + gamma * gamma)
            t = 0.5_dp * (u(i) + s)
            if (t <= 0.0_dp .or. s <= 0.0_dp) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 3
                return
            end if

            if (abs(lambda) <= 1.0e-10_dp) then
                z(i) = log(t)
            else
                z(i) = (t**lambda - 1.0_dp) / lambda
            end if
            log_u_plus_s(i) = log(u(i) + s)
            log_one_plus(i) = log(1.0_dp + u(i) / s)
        end do

        if (adjust) then
            call mean_valid(log_u_plus_s, mean_u_plus_s, n_valid, mean_status)
            if (mean_status /= 0) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 4
                return
            end if
            call mean_valid(log_one_plus, mean_one_plus, n_valid, mean_status)
            if (mean_status /= 0) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 4
                return
            end if
            jacobian = 0.5_dp**lambda * exp((lambda - 1.0_dp) * mean_u_plus_s) * exp(mean_one_plus)
            if (jacobian <= 0.0_dp .or. ieee_is_nan(jacobian)) then
                z = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 4
                return
            end if
            z = z / jacobian
        end if
        status = 0
    end subroutine bcn_power_vector

    pure subroutine bcn_power_matrix(u, lambda, gamma, z, status, jacobian_adjusted)
        real(dp), intent(in) :: u(:, :) !! Input data matrix with variables in columns and NaNs allowed.
        real(dp), intent(in) :: lambda(:) !! BCN power parameter for each column.
        real(dp), intent(in) :: gamma(:) !! Nonnegative BCN location parameter for each column.
        real(dp), intent(out) :: z(:, :) !! BCN-transformed matrix with the same shape as u.
        integer, intent(out) :: status !! Zero on success; nonzero for incompatible shapes or transformation failure.
        logical, intent(in), optional :: jacobian_adjusted !! If true, apply Jacobian normalization to each column.
        logical :: adjust
        integer :: j
        integer :: local_status

        if (any(shape(z) /= shape(u)) .or. size(lambda) /= size(u, 2) .or. size(gamma) /= size(u, 2)) then
            status = 1
            if (size(z) > 0) z = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        adjust = .false.
        if (present(jacobian_adjusted)) adjust = jacobian_adjusted
        do j = 1, size(u, 2)
            call bcn_power_vector(u(:, j), lambda(j), gamma(j), z(:, j), local_status, adjust)
            if (local_status /= 0) then
                status = local_status
                return
            end if
        end do
        status = 0
    end subroutine bcn_power_matrix

    pure elemental function bcn_power_inverse(z, lambda, gamma) result(u)
        real(dp), intent(in) :: z !! BCN-transformed scalar.
        real(dp), intent(in) :: lambda !! BCN power parameter used to obtain z.
        real(dp), intent(in) :: gamma !! Nonnegative BCN location parameter used to obtain z.
        real(dp) :: u
        real(dp) :: q

        if (ieee_is_nan(z) .or. ieee_is_nan(lambda) .or. ieee_is_nan(gamma)) then
            u = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (abs(lambda) < 1.0e-7_dp) then
            q = 2.0_dp * exp(z)
        else
            q = 2.0_dp * (lambda * z + 1.0_dp)**(1.0_dp / lambda)
        end if
        u = (q * q - gamma * gamma) / (2.0_dp * q)
    end function bcn_power_inverse

    pure subroutine logit_transform(p, transformed, status, percents, adjust, interpreted_as_percents, used_adjust)
        real(dp), intent(in) :: p(:) !! Proportions in [0,1] or percentages in [0,100]; NaNs are ignored in range checks and propagated.
        real(dp), intent(out) :: transformed(:) !! Logit-transformed values with the same length as p.
        integer, intent(out) :: status !! Zero on success; nonzero for shape mismatch, empty data, invalid range, or invalid adjustment.
        logical, intent(in), optional :: percents !! Explicitly interpret p as percentages when true; if absent infer from max(p) > 1.
        real(dp), intent(in), optional :: adjust !! Endpoint remapping amount a in [0,0.5); absent uses 0.025 only when an endpoint is present.
        logical, intent(out), optional :: interpreted_as_percents !! Returns whether percentage interpretation was used.
        real(dp), intent(out), optional :: used_adjust !! Returns the endpoint adjustment actually applied.
        real(dp) :: a
        real(dp) :: alpha
        real(dp) :: max_p
        real(dp) :: min_p
        real(dp) :: prop
        logical :: as_percent
        logical :: have_valid
        integer :: i

        if (size(transformed) /= size(p)) then
            status = 1
            if (size(transformed) > 0) transformed = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        have_valid = .false.
        min_p = huge(1.0_dp)
        max_p = -huge(1.0_dp)
        do i = 1, size(p)
            if (.not. ieee_is_nan(p(i))) then
                have_valid = .true.
                min_p = min(min_p, p(i))
                max_p = max(max_p, p(i))
            end if
        end do
        if (.not. have_valid) then
            transformed = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 2
            return
        end if

        if (present(percents)) then
            as_percent = percents
        else
            as_percent = max_p > 1.0_dp
        end if
        if (present(interpreted_as_percents)) interpreted_as_percents = as_percent

        if (as_percent) then
            if (min_p < 0.0_dp .or. max_p > 100.0_dp) then
                transformed = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 3
                return
            end if
            min_p = min_p / 100.0_dp
            max_p = max_p / 100.0_dp
        else
            if (min_p < 0.0_dp .or. max_p > 1.0_dp) then
                transformed = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 3
                return
            end if
        end if

        if (present(adjust)) then
            a = adjust
        else
            if (abs(min_p) <= 32.0_dp * epsilon(1.0_dp) .or. &
                abs(max_p - 1.0_dp) <= 32.0_dp * epsilon(1.0_dp)) then
                a = 0.025_dp
            else
                a = 0.0_dp
            end if
        end if
        if (a < 0.0_dp .or. a >= 0.5_dp) then
            transformed = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 4
            return
        end if
        if (present(used_adjust)) used_adjust = a

        alpha = 1.0_dp - 2.0_dp * a
        do i = 1, size(p)
            if (ieee_is_nan(p(i))) then
                transformed(i) = p(i)
                cycle
            end if
            prop = p(i)
            if (as_percent) prop = prop / 100.0_dp
            prop = 0.5_dp + alpha * (prop - 0.5_dp)
            transformed(i) = log(prop / (1.0_dp - prop))
        end do
        status = 0
    end subroutine logit_transform

    pure subroutine box_cox_variable(y, variable, status)
        real(dp), intent(in) :: y(:) !! Positive response values; NaNs are omitted from the geometric mean and propagated.
        real(dp), intent(out) :: variable(:) !! Constructed Box-Cox variable y * (log(y/geometric_mean) - 1).
        integer, intent(out) :: status !! Zero on success; nonzero for shape mismatch, nonpositive finite values, or empty valid input.
        real(dp), allocatable :: log_y(:)
        real(dp) :: geometric_mean
        real(dp) :: mean_log
        integer :: i
        integer :: n_valid
        integer :: mean_status

        if (size(variable) /= size(y)) then
            status = 1
            if (size(variable) > 0) variable = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        allocate(log_y(size(y)))
        log_y = ieee_value(0.0_dp, ieee_quiet_nan)

        do i = 1, size(y)
            if (ieee_is_nan(y(i))) cycle
            if (y(i) <= 0.0_dp) then
                variable = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 2
                return
            end if
            log_y(i) = log(y(i))
        end do
        call mean_valid(log_y, mean_log, n_valid, mean_status)
        if (mean_status /= 0) then
            variable = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 3
            return
        end if
        geometric_mean = exp(mean_log)

        do i = 1, size(y)
            if (ieee_is_nan(y(i))) then
                variable(i) = y(i)
            else
                variable(i) = y(i) * (log(y(i) / geometric_mean) - 1.0_dp)
            end if
        end do
        status = 0
    end subroutine box_cox_variable

end module car_transforms
