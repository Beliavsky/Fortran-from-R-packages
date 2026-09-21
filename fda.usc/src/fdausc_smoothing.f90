module fdausc_smoothing
    use r_kinds, only : dp
    use r_linalg, only : solve_system, inverse_matrix
    use fdausc_kernels, only : kernel_value
    implicit none
    private

    public :: smoothing_nw, smoothing_knn, smoothing_llr, smoothing_lpr, smoothing_lcr
    public :: cv_score, gcv_score, gccv_score, residual_variance, fitted_variance
    public :: penalty_matrix, matrix_trace, basis_smoothing_matrix

contains

    pure subroutine smoothing_nw(t, bandwidth, s, kernel_code, weights, cv)
        real(dp), intent(in) :: t(:) !! Ordered discretization points used to form pairwise kernel arguments.
        real(dp), intent(in) :: bandwidth !! Positive kernel bandwidth.
        real(dp), intent(out) :: s(:, :) !! Nadaraya-Watson smoothing matrix with shape [size(t),size(t)].
        integer, intent(in), optional :: kernel_code !! Kernel selector matching fdausc_kernels; default normal.
        real(dp), intent(in), optional :: weights(:) !! Optional observation weights applied to smoothing-matrix columns.
        logical, intent(in), optional :: cv !! If true, leave each observation out by zeroing its diagonal kernel weight.
        integer :: i
        integer :: j
        integer :: code
        real(dp) :: rowsum
        logical :: leave_one_out

        code = 1
        if (present(kernel_code)) code = kernel_code
        leave_one_out = .false.
        if (present(cv)) leave_one_out = cv
        do i = 1, size(t)
            do j = 1, size(t)
                if (leave_one_out .and. i == j) then
                    s(i, j) = 0.0_dp
                else
                    s(i, j) = kernel_value(abs(t(i) - t(j))/bandwidth, code)
                    if (present(weights)) s(i, j) = s(i, j)*weights(j)
                end if
            end do
            rowsum = sum(s(i, :))
            if (abs(rowsum) <= tiny(1.0_dp)) rowsum = 1.0e-28_dp
            s(i, :) = s(i, :)/rowsum
        end do
    end subroutine smoothing_nw

    pure subroutine smoothing_knn(t, k_neighbors, s, kernel_code, weights, cv)
        real(dp), intent(in) :: t(:) !! Ordered discretization points used to define nearest-neighbor distances.
        integer, intent(in) :: k_neighbors !! Positive nearest-neighbor count, usually between one and size(t).
        real(dp), intent(out) :: s(:, :) !! K-nearest-neighbor smoothing matrix with shape [size(t),size(t)].
        integer, intent(in), optional :: kernel_code !! Kernel selector matching fdausc_kernels; default uniform.
        real(dp), intent(in), optional :: weights(:) !! Optional observation weights applied to smoothing-matrix columns.
        logical, intent(in), optional :: cv !! If true, exclude each target point from its own neighborhood.
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: work(:)
        real(dp) :: h
        real(dp) :: rowsum
        integer :: i
        integer :: j
        integer :: kk
        integer :: code
        logical :: leave_one_out

        code = 6
        if (present(kernel_code)) code = kernel_code
        leave_one_out = .false.
        if (present(cv)) leave_one_out = cv
        kk = max(1, min(k_neighbors, size(t)))
        allocate(d(size(t)), work(size(t)))
        do i = 1, size(t)
            d = abs(t(i) - t)
            if (leave_one_out) d(i) = huge(1.0_dp)
            work = d
            call sort_increasing(work)
            h = work(min(kk, size(work))) + 1.0e-19_dp
            do j = 1, size(t)
                if (leave_one_out .and. i == j) then
                    s(i, j) = 0.0_dp
                else
                    s(i, j) = kernel_value(d(j)/h, code)
                    if (present(weights)) s(i, j) = s(i, j)*weights(j)
                end if
            end do
            rowsum = sum(s(i, :))
            if (abs(rowsum) <= tiny(1.0_dp)) rowsum = 1.0e-28_dp
            s(i, :) = s(i, :)/rowsum
        end do
    end subroutine smoothing_knn

    pure subroutine smoothing_llr(t, bandwidth, s, kernel_code, weights, cv)
        real(dp), intent(in) :: t(:) !! Ordered discretization points for local-linear smoothing.
        real(dp), intent(in) :: bandwidth !! Positive local-linear kernel bandwidth.
        real(dp), intent(out) :: s(:, :) !! Local-linear smoothing matrix with shape [size(t),size(t)].
        integer, intent(in), optional :: kernel_code !! Kernel selector matching fdausc_kernels; default normal.
        real(dp), intent(in), optional :: weights(:) !! Optional target-row weights, matching upstream S.LLR semantics.
        logical, intent(in), optional :: cv !! If true, zero diagonal contributions for leave-one-out fitting.
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: k(:)
        real(dp) :: s1
        real(dp) :: s2
        real(dp) :: rowsum
        integer :: i
        integer :: j
        integer :: code
        logical :: leave_one_out

        code = 1
        if (present(kernel_code)) code = kernel_code
        leave_one_out = .false.
        if (present(cv)) leave_one_out = cv
        allocate(d(size(t)), k(size(t)))
        do i = 1, size(t)
            d = t(i) - t
            do j = 1, size(t)
                k(j) = kernel_value(d(j)/bandwidth, code)
            end do
            if (leave_one_out) k(i) = 0.0_dp
            s1 = sum(k*d)
            s2 = sum(k*d*d)
            s(i, :) = k*(s2 - d*s1)
            if (leave_one_out) s(i, i) = 0.0_dp
            if (present(weights)) s(i, :) = s(i, :)*weights(i)
            rowsum = sum(s(i, :))
            if (abs(rowsum) <= tiny(1.0_dp)) rowsum = 1.0e-28_dp
            s(i, :) = s(i, :)/rowsum
        end do
    end subroutine smoothing_llr

    subroutine smoothing_lpr(t, bandwidth, degree, s, kernel_code, weights, cv, info)
        real(dp), intent(in) :: t(:) !! Ordered discretization points for local-polynomial smoothing.
        real(dp), intent(in) :: bandwidth !! Positive kernel bandwidth for every target point.
        integer, intent(in) :: degree !! Polynomial degree, zero or greater.
        real(dp), intent(out) :: s(:, :) !! Local-polynomial smoothing matrix with shape [size(t),size(t)].
        integer, intent(in), optional :: kernel_code !! Kernel selector matching fdausc_kernels; default normal.
        real(dp), intent(in), optional :: weights(:) !! Optional target-row weights, matching upstream S.LPR semantics.
        logical, intent(in), optional :: cv !! If true, zero each diagonal kernel contribution before fitting.
        integer, intent(out), optional :: info !! Zero when all local systems solve; positive if a local system is singular.
        real(dp), allocatable :: xmat(:, :)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: coef(:)
        real(dp), allocatable :: kval(:)
        real(dp), allocatable :: d(:)
        real(dp) :: rowsum
        integer :: i
        integer :: j
        integer :: a
        integer :: b
        integer :: code
        integer :: stat
        logical :: leave_one_out

        code = 1
        if (present(kernel_code)) code = kernel_code
        leave_one_out = .false.
        if (present(cv)) leave_one_out = cv
        if (present(info)) info = 0
        allocate(xmat(size(t), degree + 1), gram(degree + 1, degree + 1))
        allocate(rhs(degree + 1), coef(degree + 1), kval(size(t)), d(size(t)))
        do i = 1, size(t)
            d = t - t(i)
            do j = 1, size(t)
                kval(j) = kernel_value(d(j)/bandwidth, code)/bandwidth
            end do
            if (leave_one_out) kval(i) = 0.0_dp
            xmat(:, 1) = 1.0_dp
            do a = 1, degree
                xmat(:, a + 1) = d**a
            end do
            gram = 0.0_dp
            do a = 1, degree + 1
                do b = 1, degree + 1
                    gram(a, b) = sum(kval*xmat(:, a)*xmat(:, b))
                end do
            end do
            rhs = 0.0_dp
            rhs(1) = 1.0_dp
            call solve_system(gram, rhs, coef, stat)
            if (stat /= 0) then
                s(i, :) = 0.0_dp
                if (present(info)) info = stat
                cycle
            end if
            do j = 1, size(t)
                s(i, j) = kval(j)*dot_product(coef, xmat(j, :))
            end do
            if (leave_one_out) s(i, i) = 0.0_dp
            if (present(weights)) s(i, :) = s(i, :)*weights(i)
            rowsum = sum(s(i, :))
            if (abs(rowsum) <= tiny(1.0_dp)) rowsum = 1.0e-28_dp
            s(i, :) = s(i, :)/rowsum
        end do
    end subroutine smoothing_lpr

    subroutine smoothing_lcr(t, bandwidth, s, kernel_code, weights, cv, info)
        real(dp), intent(in) :: t(:) !! Ordered discretization points for local-cubic smoothing.
        real(dp), intent(in) :: bandwidth !! Positive local-cubic kernel bandwidth.
        real(dp), intent(out) :: s(:, :) !! Local-cubic smoothing matrix with shape [size(t),size(t)].
        integer, intent(in), optional :: kernel_code !! Kernel selector matching fdausc_kernels; default normal.
        real(dp), intent(in), optional :: weights(:) !! Optional target-row weights.
        logical, intent(in), optional :: cv !! If true, construct leave-one-out local fits.
        integer, intent(out), optional :: info !! Zero on success; positive if a local polynomial system is singular.
        integer :: code
        integer :: stat
        logical :: leave_one_out

        code = 1
        if (present(kernel_code)) code = kernel_code
        leave_one_out = .false.
        if (present(cv)) leave_one_out = cv
        call smoothing_lpr(t, bandwidth, 3, s, code, weights, leave_one_out, stat)
        if (present(info)) info = stat
    end subroutine smoothing_lcr

    pure real(dp) function cv_score(y, s, weights) result(score)
        real(dp), intent(in) :: y(:) !! Observed scalar response vector.
        real(dp), intent(in) :: s(:, :) !! Linear smoothing or hat matrix for y.
        real(dp), intent(in), optional :: weights(:) !! Optional diagonal loss weights, default one.
        real(dp), allocatable :: residual(:)
        real(dp) :: denom
        real(dp) :: term
        integer :: i

        residual = y - matmul(s, y)
        score = 0.0_dp
        do i = 1, size(y)
            denom = 1.0_dp - s(i, i)
            if (abs(denom) <= tiny(1.0_dp)) then
                score = huge(1.0_dp)
                return
            end if
            term = residual(i)**2/(denom*denom)
            if (present(weights)) term = term*weights(i)
            score = score + term
        end do
        score = score/real(size(y), dp)
    end function cv_score

    pure subroutine gcv_score(y, s, score, df, criteria_code, weights)
        real(dp), intent(in) :: y(:) !! Observed scalar response vector.
        real(dp), intent(in) :: s(:, :) !! Linear smoothing or hat matrix for y.
        real(dp), intent(out) :: score !! Generalized cross-validation or related criterion value.
        real(dp), intent(out) :: df !! Effective degrees of freedom trace(S).
        integer, intent(in), optional :: criteria_code !! Criterion selector: 1 GCV, 2 AIC, 3 FPE, 4 Shibata, 5 Rice; default GCV.
        real(dp), intent(in), optional :: weights(:) !! Optional diagonal quadratic-loss weights, default one.
        real(dp), allocatable :: residual(:)
        real(dp) :: res
        real(dp) :: dmean
        real(dp) :: factor
        integer :: code
        integer :: i

        code = 1
        if (present(criteria_code)) code = criteria_code
        residual = y - matmul(s, y)
        res = 0.0_dp
        do i = 1, size(y)
            if (present(weights)) then
                res = res + weights(i)*residual(i)**2
            else
                res = res + residual(i)**2
            end if
        end do
        df = matrix_trace(s)
        dmean = df/real(size(y), dp)
        select case (code)
        case (2)
            factor = exp(2.0_dp*dmean)
        case (3)
            if (abs(1.0_dp - dmean) <= tiny(1.0_dp)) then
                factor = huge(1.0_dp)
            else
                factor = (1.0_dp + dmean)/(1.0_dp - dmean)
            end if
        case (4)
            factor = 1.0_dp + 2.0_dp*dmean
        case (5)
            if (dmean > 0.5_dp) then
                factor = huge(1.0_dp)
            else
                factor = 1.0_dp/(1.0_dp - 2.0_dp*dmean)
            end if
        case default
            if (abs(1.0_dp - dmean) <= tiny(1.0_dp)) then
                factor = huge(1.0_dp)
            else
                factor = (1.0_dp - dmean)**(-2)
            end if
        end select
        score = res*factor/real(size(y), dp)
    end subroutine gcv_score

    subroutine gccv_score(y, s, score, df, criteria_code, w, info)
        real(dp), intent(in) :: y(:) !! Observed scalar response vector.
        real(dp), intent(in) :: s(:, :) !! Linear smoothing or hat matrix for y.
        real(dp), intent(out) :: score !! Generalized correlated cross-validation value.
        real(dp), intent(out) :: df !! Criterion-specific effective degrees of freedom.
        integer, intent(in), optional :: criteria_code !! Selector: 1 GCCV1, 2 GCCV2, 3 GCCV3, 4 ordinary GCV; default 1.
        real(dp), intent(in), optional :: w(:, :) !! Positive-definite precision matrix; identity if omitted.
        integer, intent(out), optional :: info !! Zero on success; positive when inversion of w fails.
        real(dp), allocatable :: residual(:)
        real(dp), allocatable :: sigma(:, :)
        real(dp), allocatable :: sc(:, :)
        real(dp) :: mse
        integer :: code
        integer :: stat
        integer :: n

        n = size(y)
        code = 1
        if (present(criteria_code)) code = criteria_code
        if (present(info)) info = 0
        residual = y - matmul(s, y)
        mse = sum(residual*residual)/real(n, dp)
        if (code == 4) then
            df = matrix_trace(s)
        else if (present(w)) then
            call inverse_matrix(w, sigma, stat)
            if (stat /= 0) then
                score = huge(1.0_dp)
                df = huge(1.0_dp)
                if (present(info)) info = stat
                return
            end if
            select case (code)
            case (2)
                df = matrix_trace(matmul(s, sigma))
            case (3)
                df = matrix_trace(matmul(matmul(s, sigma), transpose(s)))
            case default
                sc = matmul(s, sigma)
                df = matrix_trace(2.0_dp*sc - matmul(sc, transpose(s)))
            end select
        else
            select case (code)
            case (2)
                df = matrix_trace(s)
            case (3)
                df = matrix_trace(matmul(s, transpose(s)))
            case default
                df = matrix_trace(2.0_dp*s - matmul(s, transpose(s)))
            end select
        end if
        if (abs(1.0_dp - df/real(n, dp)) <= tiny(1.0_dp)) then
            score = huge(1.0_dp)
        else
            score = mse/(1.0_dp - df/real(n, dp))**2
        end if
    end subroutine gccv_score

    pure subroutine residual_variance(y, s, var_e)
        real(dp), intent(in) :: y(:) !! Observed scalar response vector.
        real(dp), intent(in) :: s(:, :) !! Linear smoothing matrix applied to y.
        real(dp), intent(out) :: var_e(:, :) !! Estimated residual covariance, scalar variance times identity.
        real(dp), allocatable :: residual(:)
        real(dp) :: se
        real(dp) :: df
        integer :: i

        df = matrix_trace(s)
        residual = y - matmul(s, y)
        se = sum(residual*residual)/max(1.0_dp, real(size(y), dp) - df)
        var_e = 0.0_dp
        do i = 1, min(size(var_e, 1), size(var_e, 2))
            var_e(i, i) = se
        end do
    end subroutine residual_variance

    pure subroutine fitted_variance(s, var_e, var_y)
        real(dp), intent(in) :: s(:, :) !! Linear smoothing matrix.
        real(dp), intent(in) :: var_e(:, :) !! Residual covariance matrix conformable with s.
        real(dp), intent(out) :: var_y(:, :) !! Fitted covariance S*Var(e)*transpose(S).
        var_y = matmul(matmul(s, var_e), transpose(s))
    end subroutine fitted_variance

    pure subroutine penalty_matrix(t, pcoef, penalty)
        real(dp), intent(in) :: t(:) !! Ordered grid used to scale the finite-difference penalty.
        real(dp), intent(in) :: pcoef(:) !! Difference-order weights, e.g. [0,0,1] for a second-difference penalty.
        real(dp), intent(out) :: penalty(:, :) !! Symmetric finite-difference roughness penalty matrix.
        real(dp), allocatable :: d(:, :)
        real(dp), allocatable :: current(:, :)
        real(dp), allocatable :: next(:, :)
        integer :: n
        integer :: order
        integer :: i
        real(dp) :: h
        real(dp) :: scale

        n = size(t)
        order = max(0, size(pcoef) - 1)
        allocate(current(n, n))
        current = 0.0_dp
        do i = 1, n
            current(i, i) = 1.0_dp
        end do
        do i = 1, order
            if (size(current, 1) <= 1) exit
            allocate(next(size(current, 1) - 1, n))
            next = current(2:, :) - current(:size(current, 1) - 1, :)
            call move_alloc(next, current)
        end do
        d = current
        if (n > 1) then
            h = (t(n) - t(1))/real(n - 1, dp)
        else
            h = 1.0_dp
        end if
        scale = 1.0_dp
        if (order > 0 .and. abs(h) > tiny(1.0_dp)) scale = h**(-real(2*order - 1, dp))
        penalty = matmul(transpose(d), d)*scale
        if (size(pcoef) > 0) penalty = penalty*max(1.0_dp, abs(pcoef(size(pcoef))))
    end subroutine penalty_matrix

    pure real(dp) function matrix_trace(a) result(value)
        real(dp), intent(in) :: a(:, :) !! Matrix whose main diagonal is summed.
        integer :: i
        value = 0.0_dp
        do i = 1, min(size(a, 1), size(a, 2))
            value = value + a(i, i)
        end do
    end function matrix_trace

    pure subroutine sort_increasing(x)
        real(dp), intent(inout) :: x(:) !! Values sorted in ascending order in place.
        integer :: i
        integer :: j
        real(dp) :: key

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
    end subroutine sort_increasing

    subroutine basis_smoothing_matrix(phi, case_weights, penalty, lambda, s, info)
        real(dp), intent(in) :: phi(:, :) !! Basis design matrix with rows at discretization points and basis functions in columns.
        real(dp), intent(in) :: case_weights(:) !! Nonnegative weights for discretization points, one per row of phi.
        real(dp), intent(in) :: penalty(:, :) !! Square basis-coefficient roughness penalty matrix.
        real(dp), intent(in) :: lambda !! Nonnegative multiplier applied to the roughness penalty.
        real(dp), intent(out) :: s(:, :) !! Basis smoother matrix mapping observed values to fitted values.
        integer, intent(out) :: info !! Zero on success; positive when the penalized normal matrix is singular.
        real(dp), allocatable :: normal(:, :)
        real(dp), allocatable :: inverse(:, :)
        real(dp), allocatable :: weighted_phi(:, :)
        integer :: i

        allocate(weighted_phi(size(phi, 1), size(phi, 2)))
        do i = 1, size(phi, 1)
            weighted_phi(i, :) = case_weights(i)*phi(i, :)
        end do
        normal = matmul(transpose(phi), weighted_phi) + lambda*penalty
        call inverse_matrix(normal, inverse, info)
        if (info /= 0) then
            s = 0.0_dp
            return
        end if
        s = matmul(matmul(phi, inverse), transpose(weighted_phi))
    end subroutine basis_smoothing_matrix

end module fdausc_smoothing
