module e1071_svm_solver
    use e1071_kinds, only: dp
    implicit none
    private

    public :: solve_dual, initialize_nu_svc_alpha

contains

    subroutine solve_dual(q, p, y, alpha, cp, cn, tolerance, nu_mode, requested_max_iterations, rho, r)
        real(dp), intent(in) :: q(:, :) !! Symmetric dual Hessian matrix in LIBSVM sign convention.
        real(dp), intent(in) :: p(:) !! Linear objective term matching q and alpha lengths.
        integer, intent(in) :: y(:) !! Constraint signs, each +1 or -1.
        real(dp), intent(inout) :: alpha(:) !! Feasible initial dual variables replaced by the optimized solution.
        real(dp), intent(in) :: cp !! Upper bound for variables with positive y.
        real(dp), intent(in) :: cn !! Upper bound for variables with negative y.
        real(dp), intent(in) :: tolerance !! LIBSVM KKT stopping tolerance; must be positive.
        logical, intent(in) :: nu_mode !! If true, use the additional constant-sum working-set rule for nu formulations.
        integer, intent(in) :: requested_max_iterations !! Positive explicit iteration limit, or zero for a LIBSVM-like
        !! automatic bound.
        real(dp), intent(out) :: rho !! Bias threshold calculated from the final dual gradient.
        real(dp), intent(out) :: r !! nu-formulation scaling quantity; zero for ordinary formulations.
        real(dp), allocatable :: gradient(:)
        real(dp) :: old_i
        real(dp) :: old_j
        real(dp) :: delta
        real(dp) :: diff
        real(dp) :: sum_value
        real(dp) :: quad
        real(dp) :: ci
        real(dp) :: cj
        real(dp) :: di
        real(dp) :: dj
        integer :: i
        integer :: j
        integer :: iter
        integer :: max_iterations
        logical :: optimal

        if (size(q, 1) /= size(q, 2) .or. size(q, 1) /= size(alpha) .or. size(p) /= size(alpha) .or. size(y) /= size(alpha)) then
            error stop "solve_dual: shape mismatch"
        end if
        if (cp <= 0.0_dp .or. cn <= 0.0_dp .or. tolerance <= 0.0_dp) error stop "solve_dual: invalid bounds or tolerance"
        gradient = p + matmul(q, alpha)
        if (requested_max_iterations > 0) then
            max_iterations = requested_max_iterations
        else
            max_iterations = max(10000000, 100 * size(alpha))
        end if
        do iter = 1, max_iterations
            if (nu_mode) then
                call select_working_set_nu(q, gradient, y, alpha, cp, cn, tolerance, i, j, optimal)
            else
                call select_working_set_general(q, gradient, y, alpha, cp, cn, tolerance, i, j, optimal)
            end if
            if (optimal) exit
            ci = bound_for_sign(y(i), cp, cn)
            cj = bound_for_sign(y(j), cp, cn)
            old_i = alpha(i)
            old_j = alpha(j)
            if (y(i) /= y(j)) then
                quad = q(i, i) + q(j, j) + 2.0_dp * q(i, j)
                if (quad <= 0.0_dp) quad = 1.0e-12_dp
                delta = (-gradient(i) - gradient(j)) / quad
                diff = alpha(i) - alpha(j)
                alpha(i) = alpha(i) + delta
                alpha(j) = alpha(j) + delta
                if (diff > 0.0_dp) then
                    if (alpha(j) < 0.0_dp) then
                        alpha(j) = 0.0_dp
                        alpha(i) = diff
                    end if
                else
                    if (alpha(i) < 0.0_dp) then
                        alpha(i) = 0.0_dp
                        alpha(j) = -diff
                    end if
                end if
                if (diff > ci - cj) then
                    if (alpha(i) > ci) then
                        alpha(i) = ci
                        alpha(j) = ci - diff
                    end if
                else
                    if (alpha(j) > cj) then
                        alpha(j) = cj
                        alpha(i) = cj + diff
                    end if
                end if
            else
                quad = q(i, i) + q(j, j) - 2.0_dp * q(i, j)
                if (quad <= 0.0_dp) quad = 1.0e-12_dp
                delta = (gradient(i) - gradient(j)) / quad
                sum_value = alpha(i) + alpha(j)
                alpha(i) = alpha(i) - delta
                alpha(j) = alpha(j) + delta
                if (sum_value > ci) then
                    if (alpha(i) > ci) then
                        alpha(i) = ci
                        alpha(j) = sum_value - ci
                    end if
                else
                    if (alpha(j) < 0.0_dp) then
                        alpha(j) = 0.0_dp
                        alpha(i) = sum_value
                    end if
                end if
                if (sum_value > cj) then
                    if (alpha(j) > cj) then
                        alpha(j) = cj
                        alpha(i) = sum_value - cj
                    end if
                else
                    if (alpha(i) < 0.0_dp) then
                        alpha(i) = 0.0_dp
                        alpha(j) = sum_value
                    end if
                end if
            end if
            di = alpha(i) - old_i
            dj = alpha(j) - old_j
            gradient = gradient + q(:, i) * di + q(:, j) * dj
        end do
        if (nu_mode) then
            call calculate_rho_nu(gradient, y, alpha, cp, cn, rho, r)
        else
            call calculate_rho_general(gradient, y, alpha, cp, cn, rho)
            r = 0.0_dp
        end if
    end subroutine solve_dual

    subroutine select_working_set_general(q, gradient, y, alpha, cp, cn, tolerance, out_i, out_j, optimal)
        real(dp), intent(in) :: q(:, :) !! Symmetric dual Hessian used to estimate pairwise objective decreases.
        real(dp), intent(in) :: gradient(:) !! Current dual objective gradient.
        integer, intent(in) :: y(:) !! Constraint signs, each +1 or -1.
        real(dp), intent(in) :: alpha(:) !! Current dual variables used to determine bound status.
        real(dp), intent(in) :: cp !! Positive-sign upper bound.
        real(dp), intent(in) :: cn !! Negative-sign upper bound.
        real(dp), intent(in) :: tolerance !! KKT stopping tolerance.
        integer, intent(out) :: out_i !! Selected first working-set index when optimal is false.
        integer, intent(out) :: out_j !! Selected second working-set index when optimal is false.
        logical, intent(out) :: optimal !! True when no pair violates the KKT condition by tolerance.
        real(dp) :: gmax
        real(dp) :: gmax2
        real(dp) :: obj_min
        real(dp) :: grad_diff
        real(dp) :: quad
        real(dp) :: obj_diff
        integer :: t
        integer :: i
        integer :: j

        gmax = -huge(1.0_dp)
        gmax2 = -huge(1.0_dp)
        obj_min = huge(1.0_dp)
        i = 0
        j = 0
        do t = 1, size(alpha)
            if (y(t) == 1) then
                if (.not. is_upper(alpha(t), cp)) then
                    if (-gradient(t) >= gmax) then
                        gmax = -gradient(t)
                        i = t
                    end if
                end if
            else
                if (.not. is_lower(alpha(t))) then
                    if (gradient(t) >= gmax) then
                        gmax = gradient(t)
                        i = t
                    end if
                end if
            end if
        end do
        if (i == 0) then
            optimal = .true.
            out_i = 0
            out_j = 0
            return
        end if
        do t = 1, size(alpha)
            if (y(t) == 1) then
                if (.not. is_lower(alpha(t))) then
                    grad_diff = gmax + gradient(t)
                    gmax2 = max(gmax2, gradient(t))
                    if (grad_diff > 0.0_dp) then
                        quad = q(i, i) + q(t, t) - 2.0_dp * real(y(i), dp) * q(i, t)
                        if (quad <= 0.0_dp) quad = 1.0e-12_dp
                        obj_diff = -(grad_diff * grad_diff) / quad
                        if (obj_diff <= obj_min) then
                            obj_min = obj_diff
                            j = t
                        end if
                    end if
                end if
            else
                if (.not. is_upper(alpha(t), cn)) then
                    grad_diff = gmax - gradient(t)
                    gmax2 = max(gmax2, -gradient(t))
                    if (grad_diff > 0.0_dp) then
                        quad = q(i, i) + q(t, t) + 2.0_dp * real(y(i), dp) * q(i, t)
                        if (quad <= 0.0_dp) quad = 1.0e-12_dp
                        obj_diff = -(grad_diff * grad_diff) / quad
                        if (obj_diff <= obj_min) then
                            obj_min = obj_diff
                            j = t
                        end if
                    end if
                end if
            end if
        end do
        optimal = gmax + gmax2 < tolerance .or. j == 0
        out_i = i
        out_j = j
    end subroutine select_working_set_general

    subroutine select_working_set_nu(q, gradient, y, alpha, cp, cn, tolerance, out_i, out_j, optimal)
        real(dp), intent(in) :: q(:, :) !! Symmetric dual Hessian used to estimate pairwise objective decreases.
        real(dp), intent(in) :: gradient(:) !! Current dual objective gradient.
        integer, intent(in) :: y(:) !! Constraint signs, each +1 or -1.
        real(dp), intent(in) :: alpha(:) !! Current dual variables used to determine bound status.
        real(dp), intent(in) :: cp !! Positive-sign upper bound.
        real(dp), intent(in) :: cn !! Negative-sign upper bound.
        real(dp), intent(in) :: tolerance !! KKT stopping tolerance for both sign-specific working sets.
        integer, intent(out) :: out_i !! Selected first index sharing the second index's sign when optimal is false.
        integer, intent(out) :: out_j !! Selected second working-set index when optimal is false.
        logical, intent(out) :: optimal !! True when both sign-specific KKT gaps fall below tolerance.
        real(dp) :: gmaxp
        real(dp) :: gmaxp2
        real(dp) :: gmaxn
        real(dp) :: gmaxn2
        real(dp) :: obj_min
        real(dp) :: grad_diff
        real(dp) :: quad
        real(dp) :: obj_diff
        integer :: ip
        integer :: ineg
        integer :: j
        integer :: t

        gmaxp = -huge(1.0_dp)
        gmaxp2 = -huge(1.0_dp)
        gmaxn = -huge(1.0_dp)
        gmaxn2 = -huge(1.0_dp)
        obj_min = huge(1.0_dp)
        ip = 0
        ineg = 0
        j = 0
        do t = 1, size(alpha)
            if (y(t) == 1) then
                if (.not. is_upper(alpha(t), cp)) then
                    if (-gradient(t) >= gmaxp) then
                        gmaxp = -gradient(t)
                        ip = t
                    end if
                end if
            else
                if (.not. is_lower(alpha(t))) then
                    if (gradient(t) >= gmaxn) then
                        gmaxn = gradient(t)
                        ineg = t
                    end if
                end if
            end if
        end do
        do t = 1, size(alpha)
            if (y(t) == 1) then
                if (ip > 0 .and. .not. is_lower(alpha(t))) then
                    grad_diff = gmaxp + gradient(t)
                    gmaxp2 = max(gmaxp2, gradient(t))
                    if (grad_diff > 0.0_dp) then
                        quad = q(ip, ip) + q(t, t) - 2.0_dp * q(ip, t)
                        if (quad <= 0.0_dp) quad = 1.0e-12_dp
                        obj_diff = -(grad_diff * grad_diff) / quad
                        if (obj_diff <= obj_min) then
                            obj_min = obj_diff
                            j = t
                        end if
                    end if
                end if
            else
                if (ineg > 0 .and. .not. is_upper(alpha(t), cn)) then
                    grad_diff = gmaxn - gradient(t)
                    gmaxn2 = max(gmaxn2, -gradient(t))
                    if (grad_diff > 0.0_dp) then
                        quad = q(ineg, ineg) + q(t, t) - 2.0_dp * q(ineg, t)
                        if (quad <= 0.0_dp) quad = 1.0e-12_dp
                        obj_diff = -(grad_diff * grad_diff) / quad
                        if (obj_diff <= obj_min) then
                            obj_min = obj_diff
                            j = t
                        end if
                    end if
                end if
            end if
        end do
        optimal = max(gmaxp + gmaxp2, gmaxn + gmaxn2) < tolerance .or. j == 0
        if (optimal) then
            out_i = 0
            out_j = 0
        else
            if (y(j) == 1) then
                out_i = ip
            else
                out_i = ineg
            end if
            out_j = j
        end if
    end subroutine select_working_set_nu

    subroutine calculate_rho_general(gradient, y, alpha, cp, cn, rho)
        real(dp), intent(in) :: gradient(:) !! Final dual gradient used to identify free/bound threshold values.
        integer, intent(in) :: y(:) !! Constraint signs, each +1 or -1.
        real(dp), intent(in) :: alpha(:) !! Final dual variables.
        real(dp), intent(in) :: cp !! Positive-sign upper bound.
        real(dp), intent(in) :: cn !! Negative-sign upper bound.
        real(dp), intent(out) :: rho !! LIBSVM bias threshold calculated from free variables or bound midpoint.
        real(dp) :: upper
        real(dp) :: lower
        real(dp) :: sum_free
        real(dp) :: yg
        real(dp) :: c
        integer :: nfree
        integer :: i

        upper = huge(1.0_dp)
        lower = -huge(1.0_dp)
        sum_free = 0.0_dp
        nfree = 0
        do i = 1, size(alpha)
            c = bound_for_sign(y(i), cp, cn)
            yg = real(y(i), dp) * gradient(i)
            if (is_upper(alpha(i), c)) then
                if (y(i) == -1) then
                    upper = min(upper, yg)
                else
                    lower = max(lower, yg)
                end if
            else if (is_lower(alpha(i))) then
                if (y(i) == 1) then
                    upper = min(upper, yg)
                else
                    lower = max(lower, yg)
                end if
            else
                nfree = nfree + 1
                sum_free = sum_free + yg
            end if
        end do
        if (nfree > 0) then
            rho = sum_free / real(nfree, dp)
        else
            rho = 0.5_dp * (upper + lower)
        end if
    end subroutine calculate_rho_general

    subroutine calculate_rho_nu(gradient, y, alpha, cp, cn, rho, r)
        real(dp), intent(in) :: gradient(:) !! Final dual gradient used for the two nu-specific threshold estimates.
        integer, intent(in) :: y(:) !! Constraint signs, each +1 or -1.
        real(dp), intent(in) :: alpha(:) !! Final dual variables.
        real(dp), intent(in) :: cp !! Positive-sign upper bound.
        real(dp), intent(in) :: cn !! Negative-sign upper bound.
        real(dp), intent(out) :: rho !! Half-difference of positive/negative gradient thresholds.
        real(dp), intent(out) :: r !! Half-sum of positive/negative gradient thresholds used by nu-SVC scaling.
        real(dp) :: upper1
        real(dp) :: lower1
        real(dp) :: upper2
        real(dp) :: lower2
        real(dp) :: sum1
        real(dp) :: sum2
        real(dp) :: r1
        real(dp) :: r2
        integer :: nfree1
        integer :: nfree2
        integer :: i

        upper1 = huge(1.0_dp)
        lower1 = -huge(1.0_dp)
        upper2 = huge(1.0_dp)
        lower2 = -huge(1.0_dp)
        sum1 = 0.0_dp
        sum2 = 0.0_dp
        nfree1 = 0
        nfree2 = 0
        do i = 1, size(alpha)
            if (y(i) == 1) then
                if (is_upper(alpha(i), cp)) then
                    lower1 = max(lower1, gradient(i))
                else if (is_lower(alpha(i))) then
                    upper1 = min(upper1, gradient(i))
                else
                    nfree1 = nfree1 + 1
                    sum1 = sum1 + gradient(i)
                end if
            else
                if (is_upper(alpha(i), cn)) then
                    lower2 = max(lower2, gradient(i))
                else if (is_lower(alpha(i))) then
                    upper2 = min(upper2, gradient(i))
                else
                    nfree2 = nfree2 + 1
                    sum2 = sum2 + gradient(i)
                end if
            end if
        end do
        if (nfree1 > 0) then
            r1 = sum1 / real(nfree1, dp)
        else
            r1 = 0.5_dp * (upper1 + lower1)
        end if
        if (nfree2 > 0) then
            r2 = sum2 / real(nfree2, dp)
        else
            r2 = 0.5_dp * (upper2 + lower2)
        end if
        r = 0.5_dp * (r1 + r2)
        rho = 0.5_dp * (r1 - r2)
    end subroutine calculate_rho_nu

    subroutine initialize_nu_svc_alpha(y, nu, alpha)
        integer, intent(in) :: y(:) !! Binary class signs, each +1 or -1.
        real(dp), intent(in) :: nu !! nu-SVC parameter in (0,1], determining equal positive/negative initial alpha sums.
        real(dp), allocatable, intent(out) :: alpha(:) !! Feasible nu-SVC initial dual vector with unit upper bounds.
        real(dp) :: sum_pos
        real(dp) :: sum_neg
        integer :: i

        allocate(alpha(size(y)))
        alpha = 0.0_dp
        sum_pos = nu * real(size(y), dp) / 2.0_dp
        sum_neg = sum_pos
        do i = 1, size(y)
            if (y(i) == 1) then
                alpha(i) = min(1.0_dp, sum_pos)
                sum_pos = sum_pos - alpha(i)
            else
                alpha(i) = min(1.0_dp, sum_neg)
                sum_neg = sum_neg - alpha(i)
            end if
        end do
    end subroutine initialize_nu_svc_alpha

    pure function is_upper(alpha, c) result(value)
        real(dp), intent(in) :: alpha !! Dual variable tested against its upper bound.
        real(dp), intent(in) :: c !! Positive upper bound for the dual variable.
        logical :: value

        value = alpha >= c - 1.0e-12_dp * max(1.0_dp, c)
    end function is_upper

    pure function is_lower(alpha) result(value)
        real(dp), intent(in) :: alpha !! Dual variable tested against the zero lower bound.
        logical :: value

        value = alpha <= 1.0e-12_dp
    end function is_lower

    pure function bound_for_sign(sign_value, cp, cn) result(c)
        integer, intent(in) :: sign_value !! Constraint sign selecting the positive or negative upper bound.
        real(dp), intent(in) :: cp !! Upper bound for positive-sign variables.
        real(dp), intent(in) :: cn !! Upper bound for negative-sign variables.
        real(dp) :: c

        if (sign_value > 0) then
            c = cp
        else
            c = cn
        end if
    end function bound_for_sign

end module e1071_svm_solver
