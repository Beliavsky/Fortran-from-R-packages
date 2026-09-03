! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0 2lonly mean/norm/PMM helpers.
module mice_twolevel
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use r_kinds, only : dp
    use mice_impute_continuous, only : impute_pmm
    use mice_regression, only : impute_norm
    use mice_rng, only : mice_rng_state
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed
    implicit none
    private

    public :: impute_2lonly_mean
    public :: impute_2lonly_norm
    public :: impute_2lonly_pmm

contains

    pure subroutine impute_2lonly_mean(y, observed, cluster, where, imputed, info)
        real(dp), intent(in) :: y(:) !! Numeric response expected to be constant within level-2 clusters apart from missing entries.
        logical, intent(in) :: observed(:) !! True for observed response values.
        integer, intent(in) :: cluster(:) !! Integer level-2 cluster identifier for each row.
        logical, intent(in) :: where(:) !! True for rows where class-mean repair values are requested.
        real(dp), allocatable, intent(out) :: imputed(:) !! Cluster means at requested rows; empty clusters return NaN.
        integer, intent(out) :: info !! `mice_ok` on success or a shape status code.

        integer, allocatable :: groups(:)
        real(dp) :: mean_value, nan_value
        integer :: g, i, j, k, ng, nobs

        if (size(observed) /= size(y) .or. size(cluster) /= size(y) .or. size(where) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        call unique_groups(cluster, groups)
        ng = size(groups)
        allocate(imputed(count(where)))
        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        k = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            k = k + 1
            g = group_position(groups, cluster(i))
            mean_value = 0.0_dp
            nobs = 0
            do j = 1, size(y)
                if (group_position(groups, cluster(j)) /= g .or. .not. observed(j)) cycle
                mean_value = mean_value + y(j)
                nobs = nobs + 1
            end do
            if (nobs > 0) then
                imputed(k) = mean_value / real(nobs, dp)
            else
                imputed(k) = nan_value
            end if
        end do
        if (ng < 1) then
            info = mice_invalid_argument
        else
            info = mice_ok
        end if
    end subroutine impute_2lonly_mean

    subroutine impute_2lonly_norm(y, observed, cluster, x, where, rng, imputed, info, ridge)
        real(dp), intent(in) :: y(:) !! Level-2 numeric response repeated within each cluster.
        logical, intent(in) :: observed(:) !! True for observed response entries; partial cluster missingness is rejected.
        integer, intent(in) :: cluster(:) !! Integer level-2 cluster identifier for each row.
        real(dp), intent(in) :: x(:, :) !! Complete row predictors averaged to cluster level before normal imputation.
        logical, intent(in) :: where(:) !! True for row-level locations whose cluster imputation should be returned.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used by the cluster-level Bayesian normal imputer.
        real(dp), allocatable, intent(out) :: imputed(:) !! Cluster-level normal imputations expanded to requested rows.
        integer, intent(out) :: info !! `mice_ok` on success or a shape/partial-cluster status code.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default `1e-5`.

        call impute_2lonly_model(y, observed, cluster, x, where, rng, .false., imputed, info, ridge=ridge)
    end subroutine impute_2lonly_norm

    subroutine impute_2lonly_pmm(y, observed, cluster, x, where, rng, imputed, info, donors, ridge)
        real(dp), intent(in) :: y(:) !! Level-2 numeric response repeated within each cluster.
        logical, intent(in) :: observed(:) !! True for observed response entries; partial cluster missingness is rejected.
        integer, intent(in) :: cluster(:) !! Integer level-2 cluster identifier for each row.
        real(dp), intent(in) :: x(:, :) !! Complete row-level predictors averaged to cluster level before PMM.
        logical, intent(in) :: where(:) !! True for row-level locations whose cluster PMM value should be returned.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used by cluster-level PMM.
        real(dp), allocatable, intent(out) :: imputed(:) !! Cluster-level PMM donor values expanded to requested rows.
        integer, intent(out) :: info !! `mice_ok` on success or a shape/partial-cluster status code.
        integer, intent(in), optional :: donors !! PMM donor-pool size; default five.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default `1e-5`.

        call impute_2lonly_model(y, observed, cluster, x, where, rng, .true., imputed, info, donors=donors, ridge=ridge)
    end subroutine impute_2lonly_pmm

    subroutine impute_2lonly_model(y, observed, cluster, x, where, rng, use_pmm, imputed, info, donors, ridge)
        real(dp), intent(in) :: y(:) !! Level-2 response repeated within clusters.
        logical, intent(in) :: observed(:) !! True for observed level-2 response entries.
        integer, intent(in) :: cluster(:) !! Integer cluster identifier for each row.
        real(dp), intent(in) :: x(:, :) !! Complete predictors aggregated to cluster means.
        logical, intent(in) :: where(:) !! True for row-level positions whose cluster imputation is returned.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used by the selected cluster-level imputer.
        logical, intent(in), value :: use_pmm !! Select PMM when true and Bayesian normal regression when false.
        real(dp), allocatable, intent(out) :: imputed(:) !! Expanded cluster-level imputations.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        integer, intent(in), optional :: donors !! PMM donor-pool size; ignored for normal imputation.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default `1e-5`.

        integer, allocatable :: groups(:)
        real(dp), allocatable :: y2(:), x2(:, :), imp2(:)
        logical, allocatable :: observed2(:), where2(:)
        real(dp) :: ridge_value
        integer :: donor_count, g, i, j, k, ng, nobs, nmis

        if (size(observed) /= size(y) .or. size(cluster) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        call unique_groups(cluster, groups)
        ng = size(groups)
        if (ng < 1) then
            info = mice_invalid_argument
            return
        end if
        allocate(y2(ng), x2(ng, size(x, 2)), observed2(ng), where2(ng))
        do g = 1, ng
            nobs = 0
            nmis = 0
            y2(g) = 0.0_dp
            x2(g, :) = 0.0_dp
            k = 0
            do i = 1, size(y)
                if (cluster(i) /= groups(g)) cycle
                k = k + 1
                if (observed(i)) then
                    nobs = nobs + 1
                    y2(g) = y2(g) + y(i)
                else
                    nmis = nmis + 1
                end if
                if (size(x, 2) > 0) x2(g, :) = x2(g, :) + x(i, :)
            end do
            if (nobs > 0 .and. nmis > 0) then
                info = mice_invalid_argument
                return
            end if
            if (nobs > 0) y2(g) = y2(g) / real(nobs, dp)
            if (k > 0 .and. size(x, 2) > 0) x2(g, :) = x2(g, :) / real(k, dp)
            observed2(g) = nobs > 0
            where2(g) = any(where .and. cluster == groups(g))
        end do
        if (count(observed2) < 1) then
            info = mice_no_observed
            return
        end if
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        donor_count = 5
        if (present(donors)) donor_count = donors
        if (use_pmm) then
            call impute_pmm(y2, observed2, x2, where2, rng, imp2, info, donors=donor_count, ridge=ridge_value)
        else
            call impute_norm(y2, observed2, x2, where2, rng, imp2, info, ridge=ridge_value)
        end if
        if (info /= mice_ok) return
        allocate(imputed(count(where)))
        k = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            k = k + 1
            g = group_position(groups, cluster(i))
            j = count(where2(1:g))
            imputed(k) = imp2(j)
        end do
        info = mice_ok
    end subroutine impute_2lonly_model

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

end module mice_twolevel
