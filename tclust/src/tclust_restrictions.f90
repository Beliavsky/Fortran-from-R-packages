module tclust_restrictions
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use tclust_kinds, only : dp
    use tclust_linalg, only : symmetric_eigen, determinant_sym, covariance_matrix
    implicit none
    private

    public :: restrict_eigenvalues
    public :: restrict_determinants
    public :: restrict_covariances
    public :: unrestricted_factor

contains

    subroutine restrict_eigenvalues(eigenvalues, sizes, factor, zero_tol, restricted)
        real(dp), intent(in) :: eigenvalues(:, :) !! Raw eigenvalues, with variables in rows and clusters in columns.
        real(dp), intent(in) :: sizes(:) !! Current effective cluster sizes; length equals number of clusters.
        real(dp), intent(in) :: factor !! Maximum allowed ratio between any active eigenvalues; must be at least one.
        real(dp), intent(in) :: zero_tol !! Numerical tolerance used to identify degenerate scatter matrices.
        real(dp), intent(out) :: restricted(:, :) !! Restricted eigenvalues with the same shape as eigenvalues.
        real(dp), allocatable :: d(:, :)
        real(dp), allocatable :: points(:)
        real(dp), allocatable :: ed(:)
        real(dp) :: c
        real(dp) :: den
        real(dp) :: e
        real(dp) :: m
        real(dp) :: n_total
        real(dp) :: num
        real(dp) :: objective
        real(dp) :: best_objective
        real(dp) :: rel_max
        real(dp) :: rel_min
        real(dp) :: mean_active
        real(dp) :: r
        real(dp) :: ssum
        real(dp) :: tsum
        integer :: i
        integer :: j
        integer :: k
        integer :: p
        integer :: mp
        integer :: np
        integer :: idx

        p = size(eigenvalues, 1)
        k = size(eigenvalues, 2)
        if (size(sizes) /= k) error stop 'restrict_eigenvalues: bad sizes length'
        c = max(factor, 1.0_dp)
        allocate(d(p, k))
        d = max(eigenvalues, 0.0_dp)

        rel_max = 0.0_dp
        rel_min = huge(1.0_dp)
        mean_active = 0.0_dp
        idx = 0
        do j = 1, k
            if (sizes(j) > zero_tol) then
                rel_max = max(rel_max, maxval(d(:, j)))
                rel_min = min(rel_min, minval(d(:, j)))
                mean_active = mean_active + sum(d(:, j))
                idx = idx + p
            end if
        end do
        if (idx == 0 .or. rel_max <= zero_tol) then
            restricted = 0.0_dp
            return
        end if
        mean_active = mean_active / real(idx, dp)

        if (rel_min > zero_tol .and. rel_max / rel_min <= c * (1.0_dp + 100.0_dp * epsilon(1.0_dp))) then
            restricted = d
            do j = 1, k
                if (sizes(j) <= zero_tol) restricted(:, j) = mean_active
            end do
            return
        end if

        np = 2 * p * k
        allocate(points(np), ed(np + 1))
        idx = 0
        do j = 1, k
            do i = 1, p
                idx = idx + 1
                points(idx) = d(i, j)
                idx = idx + 1
                points(idx) = d(i, j) / c
            end do
        end do
        call sort_real(points)
        ed(1) = 0.5_dp * points(1)
        do mp = 2, np
            ed(mp) = 0.5_dp * (points(mp - 1) + points(mp))
        end do
        ed(np + 1) = 1.5_dp * points(np)

        n_total = sum(max(sizes, 0.0_dp))
        best_objective = -huge(1.0_dp)
        m = max(zero_tol, rel_min)
        do mp = 1, np + 1
            num = 0.0_dp
            den = 0.0_dp
            do j = 1, k
                if (sizes(j) <= zero_tol) cycle
                r = 0.0_dp
                ssum = 0.0_dp
                tsum = 0.0_dp
                do i = 1, p
                    if (d(i, j) < ed(mp)) then
                        r = r + 1.0_dp
                        ssum = ssum + d(i, j)
                    else if (d(i, j) > ed(mp) * c) then
                        r = r + 1.0_dp
                        tsum = tsum + d(i, j)
                    end if
                end do
                num = num + (sizes(j) / n_total) * (ssum + tsum / c)
                den = den + (sizes(j) / n_total) * r
            end do
            if (den <= tiny(1.0_dp)) cycle
            e = num / den
            if (e <= zero_tol) cycle
            objective = 0.0_dp
            do j = 1, k
                if (sizes(j) <= zero_tol) cycle
                do i = 1, p
                    r = min(max(d(i, j), e), c * e)
                    r = max(r, zero_tol)
                    objective = objective - 0.5_dp * (sizes(j) / n_total) * (log(r) + d(i, j) / r)
                end do
            end do
            if (objective > best_objective) then
                best_objective = objective
                m = e
            end if
        end do

        do j = 1, k
            do i = 1, p
                restricted(i, j) = min(max(d(i, j), m), c * m)
            end do
        end do
    end subroutine restrict_eigenvalues

    subroutine restrict_determinants(eigenvalues, sizes, factor, cshape, zero_tol, restricted)
        real(dp), intent(in) :: eigenvalues(:, :) !! Raw covariance eigenvalues, variables by clusters.
        real(dp), intent(in) :: sizes(:) !! Effective cluster sizes used as restriction weights.
        real(dp), intent(in) :: factor !! Maximum determinant ratio between clusters; must be at least one.
        real(dp), intent(in) :: cshape !! Maximum within-cluster eigenvalue ratio used to constrain shape.
        real(dp), intent(in) :: zero_tol !! Numerical tolerance for degenerate eigenvalues.
        real(dp), intent(out) :: restricted(:, :) !! Restricted eigenvalues, variables by clusters.
        real(dp), allocatable :: shape_ev(:, :)
        real(dp), allocatable :: gm(:, :)
        real(dp), allocatable :: scale_raw(:, :)
        real(dp), allocatable :: scale_res(:, :)
        real(dp) :: prod_ev
        integer :: i
        integer :: j
        integer :: k
        integer :: p

        p = size(eigenvalues, 1)
        k = size(eigenvalues, 2)
        if (p == 1) then
            call restrict_eigenvalues(eigenvalues, sizes, factor, zero_tol, restricted)
            return
        end if
        allocate(shape_ev(p, k), gm(p, k), scale_raw(1, k), scale_res(1, k))
        do j = 1, k
            call restrict_eigenvalues(eigenvalues(:, j:j), [1.0_dp], cshape, zero_tol, shape_ev(:, j:j))
            prod_ev = product(shape_ev(:, j))
            if (prod_ev <= zero_tol) prod_ev = 1.0_dp
            gm(:, j) = shape_ev(:, j) / prod_ev ** (1.0_dp / real(p, dp))
            do i = 1, p
                if (gm(i, j) <= 0.0_dp) gm(i, j) = 1.0_dp
            end do
            scale_raw(1, j) = sum(eigenvalues(:, j) / gm(:, j)) / real(p, dp)
            if (ieee_is_nan(scale_raw(1, j))) scale_raw(1, j) = 0.0_dp
        end do
        call restrict_eigenvalues(scale_raw, sizes, max(factor, 1.0_dp) ** (1.0_dp / real(p, dp)), &
                                  zero_tol, scale_res)
        do j = 1, k
            restricted(:, j) = gm(:, j) * scale_res(1, j)
        end do
    end subroutine restrict_determinants

    subroutine restrict_covariances(cov, sizes, restriction, factor, cshape, zero_tol, code)
        real(dp), intent(inout) :: cov(:, :, :) !! Cluster covariance matrices, modified in place after restriction.
        real(dp), intent(in) :: sizes(:) !! Effective cluster sizes corresponding to covariance slices.
        character(len=*), intent(in) :: restriction !! Restriction mode: 'eigen' or 'deter'.
        real(dp), intent(in) :: factor !! Across-cluster restriction factor.
        real(dp), intent(in) :: cshape !! Within-cluster shape restriction used for determinant mode.
        real(dp), intent(in) :: zero_tol !! Numerical tolerance used to detect a fully singular solution.
        integer, intent(out) :: code !! One for a usable restricted covariance set, zero for a degenerate set.
        real(dp), allocatable :: eig(:, :)
        real(dp), allocatable :: reig(:, :)
        real(dp), allocatable :: vec(:, :, :)
        real(dp), allocatable :: tmp(:, :)
        integer :: j
        integer :: k
        integer :: p

        p = size(cov, 1)
        k = size(cov, 3)
        allocate(eig(p, k), reig(p, k), vec(p, p, k), tmp(p, p))
        do j = 1, k
            call symmetric_eigen(cov(:, :, j), eig(:, j), vec(:, :, j))
        end do
        eig = max(eig, 0.0_dp)
        select case (trim(adjustl(restriction)))
        case ('eigen')
            call restrict_eigenvalues(eig, sizes, factor, zero_tol, reig)
        case ('deter')
            call restrict_determinants(eig, sizes, factor, cshape, zero_tol, reig)
        case default
            error stop 'restrict_covariances: supported restrictions are eigen and deter'
        end select
        if (maxval(reig) <= zero_tol) then
            code = 0
            return
        end if
        do j = 1, k
            tmp = 0.0_dp
            tmp = matmul(vec(:, :, j) * spread(reig(:, j), 1, p), transpose(vec(:, :, j)))
            cov(:, :, j) = 0.5_dp * (tmp + transpose(tmp))
        end do
        code = 1
    end subroutine restrict_covariances

    real(dp) function unrestricted_factor(x, cluster, k, determinant_mode) result(factor)
        real(dp), intent(in) :: x(:, :) !! Original observations in rows and variables in columns.
        integer, intent(in) :: cluster(:) !! Cluster labels, with zero denoting trimmed observations.
        integer, intent(in) :: k !! Number of positive cluster labels to inspect.
        logical, intent(in) :: determinant_mode !! True to compare determinants; false to compare covariance eigenvalues.
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        real(dp), allocatable :: sub(:, :)
        real(dp) :: vmax
        real(dp) :: vmin
        real(dp) :: val
        integer :: count_j
        integer :: i
        integer :: j
        integer :: p
        integer :: pos

        p = size(x, 2)
        allocate(cov(p, p), eig(p), vec(p, p))
        vmax = 0.0_dp
        vmin = huge(1.0_dp)
        do j = 1, k
            count_j = count(cluster == j)
            if (count_j <= 1) cycle
            allocate(sub(count_j, p))
            pos = 0
            do i = 1, size(cluster)
                if (cluster(i) == j) then
                    pos = pos + 1
                    sub(pos, :) = x(i, :)
                end if
            end do
            call covariance_matrix(sub, cov, unbiased=.true.)
            if (determinant_mode) then
                val = determinant_sym(cov)
                if (val > 0.0_dp) then
                    vmax = max(vmax, val)
                    vmin = min(vmin, val)
                end if
            else
                call symmetric_eigen(cov, eig, vec)
                if (maxval(eig) > 0.0_dp) vmax = max(vmax, maxval(eig))
                if (minval(eig, mask=eig > 0.0_dp) > 0.0_dp) vmin = min(vmin, minval(eig, mask=eig > 0.0_dp))
            end if
            deallocate(sub)
        end do
        if (vmax <= 0.0_dp .or. vmin >= 0.5_dp * huge(1.0_dp)) then
            factor = 1.0_dp
        else
            factor = ceiling(vmax / vmin)
        end if
    end function unrestricted_factor

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:) !! Real vector sorted in ascending order in place.
        real(dp) :: t
        integer :: i
        integer :: j
        integer :: m

        do i = 1, size(x) - 1
            m = i
            do j = i + 1, size(x)
                if (x(j) < x(m)) m = j
            end do
            if (m /= i) then
                t = x(i)
                x(i) = x(m)
                x(m) = t
            end if
        end do
    end subroutine sort_real

end module tclust_restrictions
