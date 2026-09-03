! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Heterogeneous two-level normal Gibbs sampler derived from mice.impute.2l.norm.R.
module mice_2l_norm
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mice_regression, only : chol_lower
    use mice_rng, only : mice_rng_state, rng_chisq, rng_gamma, rng_normal
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed, mice_singular
    implicit none
    private

    type, public :: mice_2l_norm_state
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: population_mean(:)
        real(dp), allocatable :: random_precision(:, :)
        real(dp), allocatable :: residual_precision(:)
        real(dp) :: sigma2_0 = 1.0_dp
        real(dp) :: theta = 1.0_dp
    end type mice_2l_norm_state

    public :: impute_2l_norm

contains

    subroutine impute_2l_norm(y, observed, cluster, random_x, where, rng, imputed, info, intercept, n_iter, ridge, state)
        real(dp), intent(in) :: y(:) !! Incomplete level-1 numeric response.
        logical, intent(in) :: observed(:) !! True for response values used in the two-level Gibbs updates.
        integer, intent(in) :: cluster(:) !! Integer class identifier for every row; every class must contain an observed response.
        real(dp), intent(in) :: random_x(:, :) !! Complete predictors treated as fixed and random effects; no intercept.
        logical, intent(in) :: where(:) !! True for rows where posterior-predictive imputations are requested.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for all Gibbs and posterior-predictive draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Two-level normal imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        logical, intent(in), optional :: intercept !! Add a random intercept when true; default true as in upstream mice.
        integer, intent(in), optional :: n_iter !! Number of Gibbs iterations; default 100 as in upstream mice.
        real(dp), intent(in), optional :: ridge !! Relative symmetric ridge; default `1e-4` as in upstream `symridge`.
        type(mice_2l_norm_state), intent(out), optional :: state !! Final Gibbs state for diagnostics or warm-start inspection.

        integer, allocatable :: groups(:), group_of_row(:), n_g(:)
        real(dp), allocatable :: design(:, :), bees(:, :), ss(:), mu(:), inv_psi(:, :), inv_sigma2(:)
        real(dp), allocatable :: xss(:, :, :), xgy(:, :), vv(:, :), bees_var(:, :), rhs(:), draw(:)
        real(dp), allocatable :: centered(:, :), scatter(:, :), psi_scale(:, :), psi(:, :), covariance(:, :)
        real(dp) :: denominator, geometric, harmonic, ridge_value, sigma2_0, theta
        integer :: c, df_wishart, i, iter, iterations, j, k, nclass, nrc, status
        logical :: add_intercept

        if (size(observed) /= size(y) .or. size(cluster) /= size(y) .or. size(where) /= size(y) .or. &
            size(random_x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        if (count(observed) < 1) then
            info = mice_no_observed
            return
        end if
        add_intercept = .true.
        if (present(intercept)) add_intercept = intercept
        iterations = 100
        if (present(n_iter)) iterations = n_iter
        ridge_value = 1.0e-4_dp
        if (present(ridge)) ridge_value = ridge
        if (iterations < 1 .or. ridge_value < 0.0_dp) then
            info = mice_invalid_argument
            return
        end if

        call unique_groups(cluster, groups)
        nclass = size(groups)
        nrc = size(random_x, 2)
        if (add_intercept) nrc = nrc + 1
        if (nclass < 3 .or. nrc < 1) then
            info = mice_invalid_argument
            return
        end if
        df_wishart = nclass - nrc - 1
        if (df_wishart < nrc) then
            info = mice_invalid_argument
            return
        end if
        allocate(group_of_row(size(y)), n_g(nclass))
        n_g = 0
        do i = 1, size(y)
            group_of_row(i) = group_position(groups, cluster(i))
            if (observed(i)) n_g(group_of_row(i)) = n_g(group_of_row(i)) + 1
        end do
        if (any(n_g < 1)) then
            info = mice_invalid_argument
            return
        end if

        allocate(design(size(y), nrc))
        k = 0
        if (add_intercept) then
            k = 1
            design(:, 1) = 1.0_dp
        end if
        if (size(random_x, 2) > 0) design(:, k + 1:) = random_x
        allocate(bees(nclass, nrc), ss(nclass), mu(nrc), inv_psi(nrc, nrc), inv_sigma2(nclass))
        allocate(xss(nrc, nrc, nclass), xgy(nrc, nclass))
        bees = 0.0_dp
        ss = 0.0_dp
        mu = 0.0_dp
        inv_psi = 0.0_dp
        do j = 1, nrc
            inv_psi(j, j) = 1.0_dp
        end do
        inv_sigma2 = 1.0_dp
        sigma2_0 = 1.0_dp
        theta = 1.0_dp
        xss = 0.0_dp
        xgy = 0.0_dp
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            c = group_of_row(i)
            xss(:, :, c) = xss(:, :, c) + outer_product(design(i, :), design(i, :))
            xgy(:, c) = xgy(:, c) + design(i, :) * y(i)
        end do

        allocate(vv(nrc, nrc), bees_var(nrc, nrc), rhs(nrc), draw(nrc))
        allocate(centered(nclass, nrc), scatter(nrc, nrc), psi_scale(nrc, nrc), psi(nrc, nrc), covariance(nrc, nrc))
        do iter = 1, iterations
            do c = 1, nclass
                vv = inv_sigma2(c) * xss(:, :, c) + inv_psi
                call symridge(vv, ridge_value)
                call inverse_matrix(vv, bees_var, status)
                if (status /= 0) then
                    info = mice_singular
                    return
                end if
                rhs = inv_sigma2(c) * xgy(:, c) + matmul(inv_psi, mu)
                bees(c, :) = matmul(bees_var, rhs)
                call mvn_zero(bees_var, rng, draw, status)
                if (status /= mice_ok) then
                    info = status
                    return
                end if
                bees(c, :) = bees(c, :) + draw
                ss(c) = 0.0_dp
                do i = 1, size(y)
                    if (.not. observed(i) .or. group_of_row(i) /= c) cycle
                    ss(c) = ss(c) + (y(i) - dot_product(design(i, :), bees(c, :)))**2
                end do
            end do

            call inverse_matrix(inv_psi, psi, status)
            if (status /= 0) then
                info = mice_singular
                return
            end if
            covariance = psi / real(nclass, dp)
            mu = sum(bees, dim=1) / real(nclass, dp)
            call mvn_zero(covariance, rng, draw, status)
            if (status /= mice_ok) then
                info = status
                return
            end if
            mu = mu + draw

            do c = 1, nclass
                centered(c, :) = bees(c, :) - mu
            end do
            scatter = matmul(transpose(centered), centered)
            call symridge(scatter, ridge_value)
            call inverse_matrix(scatter, psi_scale, status)
            if (status /= 0) then
                info = mice_singular
                return
            end if
            call wishart_draw(psi_scale, df_wishart, rng, inv_psi, status)
            if (status /= mice_ok) then
                info = status
                return
            end if

            do c = 1, nclass
                denominator = ss(c) * theta + sigma2_0
                if (denominator <= 0.0_dp .or. theta <= 0.0_dp) then
                    info = mice_invalid_argument
                    return
                end if
                inv_sigma2(c) = rng_gamma(rng, 0.5_dp * real(n_g(c), dp) + 0.5_dp / theta) * &
                                    (2.0_dp * theta / denominator)
                inv_sigma2(c) = max(inv_sigma2(c), tiny(1.0_dp))
            end do
            harmonic = 1.0_dp / (sum(inv_sigma2) / real(nclass, dp))
            sigma2_0 = rng_gamma(rng, real(nclass, dp) / (2.0_dp * theta) + 1.0_dp) * &
                         (2.0_dp * theta * harmonic / real(nclass, dp))
            sigma2_0 = max(sigma2_0, tiny(1.0_dp))
            geometric = exp(sum(log(1.0_dp / inv_sigma2)) / real(nclass, dp))
            denominator = real(nclass, dp) * &
                          (sigma2_0 / harmonic - log(sigma2_0) + log(geometric) - 1.0_dp)
            if (denominator <= 0.0_dp) then
                info = mice_invalid_argument
                return
            end if
            theta = 1.0_dp / rng_gamma(rng, 0.5_dp * real(nclass, dp) - 1.0_dp)
            theta = theta * denominator / 2.0_dp
            if (theta <= 0.0_dp) then
                info = mice_invalid_argument
                return
            end if
        end do

        allocate(imputed(count(where)))
        j = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            j = j + 1
            c = group_of_row(i)
            imputed(j) = dot_product(design(i, :), bees(c, :)) + rng_normal(rng) / sqrt(inv_sigma2(c))
        end do
        if (present(state)) then
            state%random_effects = bees
            state%population_mean = mu
            state%random_precision = inv_psi
            state%residual_precision = inv_sigma2
            state%sigma2_0 = sigma2_0
            state%theta = theta
        end if
        info = mice_ok
    end subroutine impute_2l_norm

    subroutine wishart_draw(scale, df, rng, matrix, info)
        real(dp), intent(in) :: scale(:, :) !! Symmetric positive-definite Wishart scale matrix.
        integer, intent(in), value :: df !! Wishart degrees of freedom; must be at least the matrix dimension.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for Bartlett Gaussian and chi-square draws.
        real(dp), intent(out) :: matrix(:, :) !! Wishart random matrix with expectation `df*scale`.
        integer, intent(out) :: info !! `mice_ok` on success or a shape/positive-definiteness status code.

        real(dp), allocatable :: lower(:, :), bartlett(:, :), product(:, :)
        integer :: i, j, n

        n = size(scale, 1)
        if (size(scale, 2) /= n .or. size(matrix, 1) /= n .or. size(matrix, 2) /= n .or. df < n) then
            info = mice_invalid_shape
            return
        end if
        allocate(lower(n, n), bartlett(n, n), product(n, n))
        call chol_lower(scale, lower, info)
        if (info /= mice_ok) return
        bartlett = 0.0_dp
        do i = 1, n
            bartlett(i, i) = sqrt(rng_chisq(rng, real(df - i + 1, dp)))
            do j = 1, i - 1
                bartlett(i, j) = rng_normal(rng)
            end do
        end do
        product = matmul(lower, bartlett)
        matrix = matmul(product, transpose(product))
        info = mice_ok
    end subroutine wishart_draw

    subroutine mvn_zero(covariance, rng, draw, info)
        real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite covariance matrix.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for independent standard-normal components.
        real(dp), intent(out) :: draw(:) !! Zero-mean multivariate-normal draw.
        integer, intent(out) :: info !! `mice_ok` on success or a positive-definiteness status code.

        real(dp), allocatable :: lower(:, :), z(:)
        integer :: j, n

        n = size(covariance, 1)
        if (size(covariance, 2) /= n .or. size(draw) /= n) then
            info = mice_invalid_shape
            return
        end if
        allocate(lower(n, n), z(n))
        call chol_lower(covariance, lower, info)
        if (info /= mice_ok) return
        do j = 1, n
            z(j) = rng_normal(rng)
        end do
        draw = matmul(lower, z)
        info = mice_ok
    end subroutine mvn_zero

    pure subroutine symridge(matrix, ridge)
        real(dp), intent(inout) :: matrix(:, :) !! Square matrix symmetrized and diagonally ridged in place.
        real(dp), intent(in), value :: ridge !! Relative diagonal ridge, matching upstream `symridge` semantics.
        integer :: j

        matrix = 0.5_dp * (matrix + transpose(matrix))
        if (size(matrix, 1) == 1) return
        do j = 1, min(size(matrix, 1), size(matrix, 2))
            matrix(j, j) = matrix(j, j) * (1.0_dp + ridge)
        end do
    end subroutine symridge

    pure subroutine unique_groups(cluster, groups)
        integer, intent(in) :: cluster(:) !! Row-level integer cluster identifiers.
        integer, allocatable, intent(out) :: groups(:) !! Unique cluster identifiers in first-occurrence order.
        integer, allocatable :: work(:)
        integer :: i, j, n
        logical :: found

        allocate(work(max(size(cluster), 1)))
        n = 0
        do i = 1, size(cluster)
            found = .false.
            do j = 1, n
                if (work(j) == cluster(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                n = n + 1
                work(n) = cluster(i)
            end if
        end do
        allocate(groups(n))
        if (n > 0) groups = work(1:n)
    end subroutine unique_groups

    pure integer function group_position(groups, value) result(position)
        integer, intent(in) :: groups(:) !! Unique cluster identifiers.
        integer, intent(in), value :: value !! Cluster identifier to locate.
        integer :: i

        position = 0
        do i = 1, size(groups)
            if (groups(i) == value) then
                position = i
                return
            end if
        end do
    end function group_position

    pure function outer_product(a, b) result(matrix)
        real(dp), intent(in) :: a(:) !! Left vector of an outer product.
        real(dp), intent(in) :: b(:) !! Right vector of an outer product.
        real(dp) :: matrix(size(a), size(b))
        integer :: i, j

        do j = 1, size(b)
            do i = 1, size(a)
                matrix(i, j) = a(i) * b(j)
            end do
        end do
    end function outer_product

end module mice_2l_norm
