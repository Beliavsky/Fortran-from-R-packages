! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Distance-aided donor selection derived from mice.impute.midastouch.R and auxiliary.R.
module mice_midastouch
    use r_kinds, only : dp
    use r_linalg, only : solve_spd
    use mice_rng, only : mice_rng_state, rng_integer, rng_uniform
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed, mice_singular
    implicit none
    private

    public :: impute_midastouch

contains

    subroutine impute_midastouch(y, observed, x, where, rng, imputed, info, ridge, midas_kappa, outout, neff)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response; observed entries form the donor pool.
        logical, intent(in) :: observed(:) !! True for rows used as donors and in the bootstrap regression.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept column.
        logical, intent(in) :: where(:) !! True for recipient rows requiring an imputed donor value.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for bootstrap counts and donor sampling.
        real(dp), allocatable, intent(out) :: imputed(:) !! Donor values ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative diagonal ridge; default `1e-5`.
        real(dp), intent(in), optional :: midas_kappa !! Distance exponent; absent means estimate it from bootstrap R-squared.
        logical, intent(in), optional :: outout !! Use one leave-one-out model per donor when true; default true.
        real(dp), intent(out), optional :: neff !! Mean effective donor-pool size from the unweighted distance probabilities.

        integer, allocatable :: obs_index(:), mis_index(:), omega_i(:)
        real(dp), allocatable :: design(:, :), xobs(:, :), xmis(:, :), yobs(:), omega(:)
        real(dp), allocatable :: xcx(:, :), xy(:), beta(:), yhat_obs(:), yhat_rec(:, :), yhat_don(:)
        real(dp), allocatable :: delta(:, :), probabilities(:, :), a(:, :), b(:), beta_loo(:)
        real(dp) :: denom, distance, effective, kappa, lower, mean_y, ridge_value, r2, rss, tss, upper, u
        integer :: donor, i, j, k, m, nmis, nobs, pick, status
        logical :: leave_one_out

        if (size(observed) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        nmis = count(where)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        if (ridge_value < 0.0_dp) then
            info = mice_invalid_argument
            return
        end if
        leave_one_out = .true.
        if (present(outout)) leave_one_out = outout

        allocate(obs_index(nobs), mis_index(nmis), yobs(nobs), omega_i(nobs), omega(nobs))
        j = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            j = j + 1
            obs_index(j) = i
            yobs(j) = y(i)
        end do
        j = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            j = j + 1
            mis_index(j) = i
        end do
        omega_i = 0
        do i = 1, nobs
            pick = rng_integer(rng, 1, nobs)
            omega_i(pick) = omega_i(pick) + 1
        end do
        omega = real(omega_i, dp)

        m = size(x, 2) + 1
        allocate(design(size(y), m), xobs(nobs, m), xmis(nmis, m))
        design(:, 1) = 1.0_dp
        if (size(x, 2) > 0) design(:, 2:) = x
        do i = 1, nobs
            xobs(i, :) = design(obs_index(i), :)
        end do
        do i = 1, nmis
            xmis(i, :) = design(mis_index(i), :)
        end do
        allocate(xcx(m, m), xy(m), beta(m), yhat_obs(nobs))
        call weighted_system(xobs, yobs, omega, ridge_value, xcx, xy)
        call solve_spd(xcx, xy, beta, status)
        if (status /= 0) then
            info = mice_singular
            return
        end if
        yhat_obs = matmul(xobs, beta)

        if (present(midas_kappa)) then
            if (midas_kappa <= 0.0_dp) then
                info = mice_invalid_argument
                return
            end if
            kappa = midas_kappa
        else
            mean_y = dot_product(yobs, omega) / real(nobs, dp)
            rss = dot_product(omega, (yobs - yhat_obs)**2)
            tss = dot_product(omega, (yobs - mean_y)**2)
            if (tss <= tiny(1.0_dp)) then
                kappa = 3.0_dp
            else
                r2 = 1.0_dp - rss / tss
                if (r2 <= 0.0_dp) then
                    kappa = 3.0_dp
                else if (r2 >= 1.0_dp) then
                    kappa = 100.0_dp
                else
                    kappa = min((50.0_dp * r2 / (1.0_dp - r2))**(3.0_dp / 8.0_dp), 100.0_dp)
                end if
            end if
        end if

        allocate(yhat_rec(nmis, nobs), yhat_don(nobs), delta(nobs, nmis), probabilities(nobs, nmis))
        if (leave_one_out) then
            allocate(a(m, m), b(m), beta_loo(m))
            do donor = 1, nobs
                a = xcx - omega(donor) * outer_product(xobs(donor, :), xobs(donor, :))
                do k = 2, m
                    a(k, k) = a(k, k) - ridge_value * omega(donor) * xobs(donor, k)**2
                end do
                b = xy - omega(donor) * xobs(donor, :) * yobs(donor)
                call repair_diagonal(a, ridge_value)
                call solve_spd(a, b, beta_loo, status)
                if (status /= 0) then
                    info = mice_singular
                    return
                end if
                yhat_don(donor) = dot_product(xobs(donor, :), beta_loo)
                if (nmis > 0) yhat_rec(:, donor) = matmul(xmis, beta_loo)
            end do
            do j = 1, nmis
                do donor = 1, nobs
                    distance = abs(yhat_don(donor) - yhat_rec(j, donor))
                    delta(donor, j) = bounded_inverse_distance(distance, kappa)
                end do
            end do
        else
            do j = 1, nmis
                distance = dot_product(xmis(j, :), beta)
                do donor = 1, nobs
                    delta(donor, j) = bounded_inverse_distance(abs(yhat_obs(donor) - distance), kappa)
                end do
            end do
        end if

        lower = sqrt(epsilon(1.0_dp))
        upper = sqrt(huge(1.0_dp))
        allocate(imputed(nmis))
        effective = 0.0_dp
        do j = 1, nmis
            probabilities(:, j) = delta(:, j) * omega
            denom = min(upper, max(lower, sum(probabilities(:, j))))
            probabilities(:, j) = probabilities(:, j) / denom
            if (present(neff)) then
                denom = min(upper, max(lower, sum(delta(:, j) * omega)))
                effective = effective + 1.0_dp / max(sum((delta(:, j) / denom)**2), lower)
            end if
            u = rng_uniform(rng)
            denom = 0.0_dp
            pick = nobs
            do donor = 1, nobs
                denom = denom + probabilities(donor, j)
                if (u <= denom) then
                    pick = donor
                    exit
                end if
            end do
            if (u > denom) pick = maxloc(probabilities(:, j), dim=1)
            imputed(j) = yobs(pick)
        end do
        if (present(neff)) then
            if (nmis > 0) then
                neff = effective / real(nmis, dp)
            else
                neff = 0.0_dp
            end if
        end if
        info = mice_ok
    end subroutine impute_midastouch

    pure subroutine weighted_system(x, y, weights, ridge, matrix, rhs)
        real(dp), intent(in) :: x(:, :) !! Design matrix including an intercept.
        real(dp), intent(in) :: y(:) !! Observed response vector.
        real(dp), intent(in) :: weights(:) !! Bootstrap frequency weights.
        real(dp), intent(in), value :: ridge !! Multiplicative ridge for non-intercept diagonal entries.
        real(dp), intent(out) :: matrix(:, :) !! Weighted and regularized `X'WX` matrix.
        real(dp), intent(out) :: rhs(:) !! Weighted `X'Wy` vector.
        integer :: i, j, k

        matrix = 0.0_dp
        rhs = 0.0_dp
        do i = 1, size(x, 1)
            do j = 1, size(x, 2)
                rhs(j) = rhs(j) + weights(i) * x(i, j) * y(i)
                do k = 1, size(x, 2)
                    matrix(j, k) = matrix(j, k) + weights(i) * x(i, j) * x(i, k)
                end do
            end do
        end do
        if (ridge > 0.0_dp) then
            do j = 2, size(x, 2)
                matrix(j, j) = matrix(j, j) * (1.0_dp + ridge)
            end do
        end if
        call repair_diagonal(matrix, ridge)
    end subroutine weighted_system

    pure subroutine repair_diagonal(matrix, ridge)
        real(dp), intent(inout) :: matrix(:, :) !! Symmetric cross-product matrix repaired in place when diagonal entries are zero.
        real(dp), intent(in), value :: ridge !! User ridge used as one candidate replacement value.
        real(dp) :: replacement
        integer :: j

        replacement = max(epsilon(1.0_dp)**0.25_dp, ridge)
        do j = 1, min(size(matrix, 1), size(matrix, 2))
            if (abs(matrix(j, j)) <= tiny(1.0_dp)) matrix(j, j) = replacement
        end do
    end subroutine repair_diagonal

    pure real(dp) function bounded_inverse_distance(distance, kappa) result(value)
        real(dp), intent(in), value :: distance !! Absolute predicted-distance magnitude.
        real(dp), intent(in), value :: kappa !! Positive distance-power exponent.
        real(dp) :: lower, upper

        lower = sqrt(epsilon(1.0_dp))
        upper = sqrt(huge(1.0_dp))
        if (distance <= 0.0_dp) then
            value = upper
        else
            value = exp(min(log(upper), max(log(lower), -kappa * log(distance))))
        end if
        value = min(upper, max(lower, value))
    end function bounded_inverse_distance

    pure function outer_product(a, b) result(matrix)
        real(dp), intent(in) :: a(:) !! Left vector of the outer product.
        real(dp), intent(in) :: b(:) !! Right vector of the outer product.
        real(dp) :: matrix(size(a), size(b))
        integer :: i, j

        do j = 1, size(b)
            do i = 1, size(a)
                matrix(i, j) = a(i) * b(j)
            end do
        end do
    end function outer_product

end module mice_midastouch
