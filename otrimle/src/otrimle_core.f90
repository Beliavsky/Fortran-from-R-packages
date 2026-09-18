module otrimle_core
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use otrimle_kinds, only : dp
    use otrimle_types, only : otrimle_fit, otrimle_grid_point
    use mclust_hierarchical, only : hc_result, hc_fit, hclass
    use otrimle_linalg, only : pi_dp, symmetric_eigen, covariance_matrix, pchisq, qchisq, &
        sort_real, count_unique_rows
    implicit none
    private

    public :: init_clust
    public :: rimle
    public :: otrimle_fit_grid
    public :: ecm_fit
    public :: cluster_to_assign

contains

    subroutine cluster_to_assign(cluster, assign)
        integer, intent(in) :: cluster(:) !! Cluster labels, where zero denotes the noise/unassigned component.
        real(dp), allocatable, intent(out) :: assign(:, :) !! One-hot assignment matrix with noise in column 1.
        integer, allocatable :: labels(:)
        integer :: g
        integer :: i
        integer :: j
        integer :: nlabels

        allocate(labels(size(cluster)))
        nlabels = 0
        do i = 1, size(cluster)
            if (cluster(i) == 0) cycle
            if (.not. any(labels(1:nlabels) == cluster(i))) then
                nlabels = nlabels + 1
                labels(nlabels) = cluster(i)
            end if
        end do
        if (nlabels > 1) call sort_int(labels(1:nlabels))
        g = nlabels
        allocate(assign(size(cluster), g + 1))
        assign = 0.0_dp
        do i = 1, size(cluster)
            if (cluster(i) == 0) then
                assign(i, 1) = 1.0_dp
            else
                do j = 1, g
                    if (cluster(i) == labels(j)) then
                        assign(i, j + 1) = 1.0_dp
                        exit
                    end if
                end do
            end if
        end do
    end subroutine cluster_to_assign

    subroutine init_clust(data, g, cluster, k, knnd_trim, model_name)
        real(dp), intent(in) :: data(:, :) !! Numeric observations in rows and variables in columns; shape n by p.
        integer, intent(in) :: g !! Requested positive number of Gaussian clusters.
        integer, allocatable, intent(out) :: cluster(:) !! Initial labels of length n; zero marks kNN-trimmed observations.
        integer, intent(in), optional :: k !! Neighbor rank used by trimming; default 3 as in InitClust().
        real(dp), intent(in), optional :: knnd_trim !! Fraction trimmed by k-nearest-neighbor distance; default 0.5.
        character(len=*), intent(in), optional :: model_name !! mclust hierarchical covariance model; default ``VVV``.
        real(dp), allocatable :: kth(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: sorted(:)
        real(dp), allocatable :: regular(:, :)
        integer, allocatable :: reg_idx(:)
        integer, allocatable :: reg_labels(:)
        type(hc_result) :: hc
        character(len=3) :: model
        integer :: i
        integer :: info
        integer :: kk
        integer :: n
        integer :: nr
        real(dp) :: trim
        real(dp) :: cutoff

        n = size(data, 1)
        if (g < 1 .or. g >= n) error stop 'init_clust: require 1 <= g < n'
        kk = 3
        if (present(k)) kk = k
        kk = max(1, min(kk, n - 1))
        trim = 0.5_dp
        if (present(knnd_trim)) trim = knnd_trim
        trim = min(max(0.0_dp, real(max(0, n - 2 * g), dp) / real(n, dp)), trim)
        model = 'VVV'
        if (present(model_name)) model = adjustl(model_name)
        allocate(kth(n), d(n - 1), sorted(n))
        do i = 1, n
            call distances_from_row(data, i, d)
            call sort_real(d)
            kth(i) = d(kk)
        end do
        sorted = kth
        call sort_real(sorted)
        cutoff = sorted(max(1, n - int(floor(real(n, dp) * trim))))
        nr = count(kth <= cutoff)
        if (nr < g) error stop 'init_clust: too few regular observations after kNN trimming'
        allocate(reg_idx(nr), regular(nr, size(data, 2)), cluster(n))
        nr = 0
        do i = 1, n
            if (kth(i) <= cutoff) then
                nr = nr + 1
                reg_idx(nr) = i
                regular(nr, :) = data(i, :)
            end if
        end do

        cluster = 0
        info = -1
        call hc_fit(regular, model, hc, minclus=g, status=info)
        if (info == 0) call hclass(hc, g, reg_labels, status=info)
        if (info == 0 .and. allocated(reg_labels)) then
            if (size(reg_labels) == nr) then
                do i = 1, nr
                    cluster(reg_idx(i)) = reg_labels(i)
                end do
                if (valid_partition(data, cluster, g)) return
            end if
        end if

        if (allocated(reg_labels)) deallocate(reg_labels)
        allocate(reg_labels(nr))
        call average_link_manhattan(regular, g, reg_labels)
        cluster = 0
        do i = 1, nr
            cluster(reg_idx(i)) = reg_labels(i)
        end do
        if (valid_partition(data, cluster, g)) return

        call deterministic_regular_partition(regular, reg_idx, data, g, cluster)
        if (.not. valid_partition(data, cluster, g)) then
            error stop 'init_clust: could not find a valid initial partition'
        end if
    end subroutine init_clust

    subroutine rimle(data, g, fit, initial, logicd, npr_max, erc, iter_max, tol)
        real(dp), intent(in) :: data(:, :) !! Numeric observations in rows; shape n by p.
        integer, intent(in) :: g !! Positive number of Gaussian components.
        type(otrimle_fit), intent(out) :: fit !! Fitted robust improper maximum-likelihood mixture.
        integer, intent(in), optional :: initial(:) !! Optional labels of length n; zero denotes initial noise.
        real(dp), intent(in), optional :: logicd !! Log improper constant density; omit for the package heuristic.
        real(dp), intent(in), optional :: npr_max !! Maximum posterior noise proportion in [0,1); default 0.5.
        real(dp), intent(in), optional :: erc !! Maximum global covariance eigenvalue ratio; default 20.
        integer, intent(in), optional :: iter_max !! Maximum ECM iterations; default 500.
        real(dp), intent(in), optional :: tol !! Absolute improper-log-likelihood convergence tolerance; default 1e-6.
        integer, allocatable :: labels(:)
        real(dp), allocatable :: assign(:, :)
        real(dp) :: ld
        real(dp) :: npmax
        real(dp) :: ratio
        real(dp) :: eps
        real(dp) :: qc
        real(dp) :: gcost
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        real(dp), allocatable :: subset(:, :)
        integer :: imax
        integer :: j
        integer :: n
        integer :: p
        integer :: m
        real(dp) :: tolerance

        n = size(data, 1)
        p = size(data, 2)
        npmax = 0.5_dp
        if (present(npr_max)) npmax = npr_max
        ratio = 20.0_dp
        if (present(erc)) ratio = erc
        imax = 500
        if (present(iter_max)) imax = iter_max
        tolerance = 1.0e-6_dp
        if (present(tol)) tolerance = tol
        if (present(initial)) then
            allocate(labels(n))
            labels = initial
        else
            call init_clust(data, g, labels, knnd_trim=npmax)
        end if
        if (.not. valid_partition(data, labels, g)) error stop 'rimle: invalid initial partition'
        call cluster_to_assign(labels, assign)
        if (size(assign, 2) /= g + 1) error stop 'rimle: initial labels do not contain g clusters'
        if (present(logicd)) then
            ld = logicd
        else
            allocate(cov(p, p), eig(p), vec(p, p))
            eps = sqrt(epsilon(1.0_dp))
            gcost = real(p, dp) * log(2.0_dp * pi_dp)
            qc = qchisq(0.95_dp ** (1.0_dp / real(n, dp)), real(p, dp))
            ld = huge(1.0_dp)
            do j = 1, g
                m = count(labels == j)
                allocate(subset(m, p))
                call pack_rows(data, labels == j, subset)
                call covariance_matrix(subset, cov, unbiased=.true.)
                call symmetric_eigen(cov, eig, vec)
                if (product(abs(eig)) < tiny(1.0_dp) .or. any(eig <= 0.0_dp) .or. &
                    any(.not. ieee_is_finite(eig))) then
                    if (maxval(eig) > 0.0_dp) then
                        where (eig < maxval(eig) / eps) eig = maxval(eig) / eps
                    else
                        eig = 1.0_dp / eps
                    end if
                end if
                ld = min(ld, -0.5_dp * (gcost + sum(log(max(eig, tiny(1.0_dp)))) + qc))
                deallocate(subset)
            end do
        end if
        if (npmax <= 0.0_dp) ld = -huge(1.0_dp)
        call ecm_fit(data, assign, ld, npmax, ratio, 0.0_dp, imax, tolerance, fit)
    end subroutine rimle

    subroutine otrimle_fit_grid(data, g, fit, optimization, initial, logicd_grid, npr_max, erc, beta, iter_max, tol)
        real(dp), intent(in) :: data(:, :) !! Numeric observations in rows; shape n by p.
        integer, intent(in) :: g !! Positive number of Gaussian components.
        type(otrimle_fit), intent(out) :: fit !! Fit at the grid point minimizing the OTRIMLE criterion.
        type(otrimle_grid_point), allocatable, intent(out) :: optimization(:) !! Per-grid criterion and likelihood results.
        integer, intent(in), optional :: initial(:) !! Optional initial labels; zero denotes noise.
        real(dp), intent(in), optional :: logicd_grid(:) !! Candidate log improper densities; omit for the upstream default grid.
        real(dp), intent(in), optional :: npr_max !! Maximum noise proportion; default 0.5.
        real(dp), intent(in), optional :: erc !! Global eigenratio constraint; default 20.
        real(dp), intent(in), optional :: beta !! Nonnegative penalty multiplying the fitted noise proportion; default zero.
        integer, intent(in), optional :: iter_max !! Maximum ECM iterations per grid point; default 500.
        real(dp), intent(in), optional :: tol !! ECM likelihood convergence tolerance; default 1e-6.
        real(dp), allocatable :: grid(:)
        real(dp), allocatable :: assign(:, :)
        integer, allocatable :: labels(:)
        type(otrimle_fit) :: trial
        real(dp) :: npmax
        real(dp) :: ratio
        real(dp) :: bpen
        real(dp) :: tolerance
        integer :: i
        integer :: best
        integer :: imax
        integer :: ngrid

        npmax = 0.5_dp
        if (present(npr_max)) npmax = npr_max
        ratio = 20.0_dp
        if (present(erc)) ratio = erc
        bpen = 0.0_dp
        if (present(beta)) bpen = beta
        imax = 500
        if (present(iter_max)) imax = iter_max
        tolerance = 1.0e-6_dp
        if (present(tol)) tolerance = tol
        if (present(initial)) then
            allocate(labels(size(initial)))
            labels = initial
        else
            call init_clust(data, g, labels, knnd_trim=npmax)
        end if
        call cluster_to_assign(labels, assign)
        if (present(logicd_grid)) then
            allocate(grid(size(logicd_grid)))
            grid = logicd_grid
            call sort_real(grid)
        else
            call default_logicd_grid(grid)
        end if
        ngrid = size(grid)
        allocate(optimization(ngrid))
        best = 0
        do i = 1, ngrid
            call ecm_fit(data, assign, grid(i), npmax, ratio, bpen, imax, tolerance, trial)
            optimization(i)%logicd = grid(i)
            optimization(i)%criterion = trial%criterion
            optimization(i)%iloglik = trial%iloglik
            if (allocated(trial%exproportion)) optimization(i)%noise_proportion = trial%exproportion(1)
            optimization(i)%code = trial%code
            optimization(i)%flags = trial%flags
            if (trial%code > 0) then
                if (best == 0) then
                    best = i
                else if (trial%criterion < optimization(best)%criterion) then
                    best = i
                end if
            end if
        end do
        if (best == 0) then
            fit%code = 0
            fit%g = g
            return
        end if
        call ecm_fit(data, assign, grid(best), npmax, ratio, bpen, imax, tolerance, fit)
    end subroutine otrimle_fit_grid

    subroutine ecm_fit(data, initial, logicd, npr_max, erc, beta, iter_max, tol, fit)
        real(dp), intent(in) :: data(:, :) !! Observations in rows; shape n by p.
        real(dp), intent(in) :: initial(:, :) !! Initial posterior/assignment matrix; shape n by g+1 with noise first.
        real(dp), intent(in) :: logicd !! Log constant improper density; very negative values approximate no-noise fitting.
        real(dp), intent(in) :: npr_max !! Upper bound on posterior noise proportion.
        real(dp), intent(in) :: erc !! Global maximum ratio across all component covariance eigenvalues.
        real(dp), intent(in) :: beta !! OTRIMLE penalty coefficient on the fitted noise proportion.
        integer, intent(in) :: iter_max !! Maximum number of ECM iterations.
        real(dp), intent(in) :: tol !! Absolute change in improper log likelihood required for convergence.
        type(otrimle_fit), intent(out) :: fit !! Fitted parameters, posterior probabilities, diagnostics, and status.
        real(dp), allocatable :: tau_old(:, :)
        real(dp), allocatable :: tau_new(:, :)
        real(dp), allocatable :: sumtau_old(:)
        real(dp), allocatable :: sumtau_new(:)
        real(dp), allocatable :: pr(:)
        real(dp), allocatable :: phi(:, :)
        real(dp), allocatable :: psi(:)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: delta(:, :)
        real(dp), allocatable :: scatter(:, :)
        real(dp), allocatable :: eig(:, :)
        real(dp), allocatable :: vec(:, :, :)
        real(dp), allocatable :: smd(:, :)
        real(dp), allocatable :: means(:, :)
        real(dp), allocatable :: covs(:, :, :)
        real(dp), allocatable :: kd(:)
        real(dp) :: icd
        real(dp) :: gausscost
        real(dp) :: sw
        real(dp) :: loglik_old
        real(dp) :: loglik_new
        real(dp) :: dloglik
        real(dp) :: root
        real(dp) :: detpart
        real(dp) :: ratio
        real(dp) :: val
        real(dp) :: bestdist
        integer :: g
        integer :: i
        integer :: iter
        integer :: j
        integer :: k_component
        integer :: n
        integer :: p
        logical :: failed
        logical :: converged
        logical :: root_ok

        n = size(data, 1)
        p = size(data, 2)
        g = size(initial, 2) - 1
        if (size(initial, 1) /= n) error stop 'ecm_fit: initial row mismatch'
        allocate(tau_old(n, g + 1), tau_new(n, g + 1), sumtau_old(g + 1), sumtau_new(g + 1))
        allocate(pr(g + 1), phi(n, g + 1), psi(n), w(n), delta(n, p), scatter(p, p))
        allocate(eig(p, g), vec(p, p, g), smd(n, g), means(p, g), covs(p, p, g), kd(g))
        tau_old = initial
        sumtau_old = sum(tau_old, dim=1)
        icd = exp(max(logicd, log(tiny(1.0_dp))))
        if (logicd <= log(tiny(1.0_dp))) icd = 0.0_dp
        gausscost = (2.0_dp * pi_dp) ** (-0.5_dp * real(p, dp))
        fit%flags = .false.
        failed = .false.
        converged = .false.
        loglik_old = -huge(1.0_dp)
        smd = 0.0_dp
        do iter = 1, iter_max
            pr = sumtau_old / real(n, dp)
            do j = 1, g
                sw = sumtau_old(j + 1)
                if (sw <= tiny(1.0_dp)) then
                    failed = .true.
                    fit%flags(1) = .true.
                    exit
                end if
                w = tau_old(:, j + 1) / sw
                means(:, j) = matmul(transpose(data), w)
                delta = data - spread(means(:, j), 1, n)
                scatter = 0.0_dp
                do i = 1, n
                    scatter = scatter + w(i) * outer(delta(i, :), delta(i, :))
                end do
                call symmetric_eigen(scatter, eig(:, j), vec(:, :, j))
            end do
            if (failed) exit
            ratio = maxval(eig) / minval(eig)
            if (.not. ieee_is_finite(ratio) .or. ratio > erc .or. any(eig <= 0.0_dp) .or. &
                any(.not. ieee_is_finite(eig))) then
                fit%flags(4) = .true.
                call restrict_eigenvalues(eig, sumtau_old(2:), erc)
            end if
            phi(:, 1) = pr(1) * icd
            do j = 1, g
                delta = data - spread(means(:, j), 1, n)
                do i = 1, n
                    val = sum((matmul(transpose(vec(:, :, j)), delta(i, :)) ** 2) / eig(:, j))
                    smd(i, j) = val
                end do
                detpart = product(eig(:, j)) ** (-0.5_dp)
                phi(:, j + 1) = pr(j + 1) * gausscost * detpart * exp(-0.5_dp * smd(:, j))
            end do
            psi = sum(phi, dim=2)
            tau_new = 0.0_dp
            do i = 1, n
                if (psi(i) > 0.0_dp .and. ieee_is_finite(psi(i))) then
                    tau_new(i, :) = phi(i, :) / psi(i)
                else if (icd > 0.0_dp) then
                    tau_new(i, 1) = 1.0_dp
                    psi(i) = max(pr(1) * icd, tiny(1.0_dp))
                else
                    bestdist = huge(1.0_dp)
                    j = 1
                    do k_component = 1, g
                        val = sum((data(i, :) - means(:, k_component)) ** 2)
                        if (val < bestdist) then
                            bestdist = val
                            j = k_component
                        end if
                    end do
                    tau_new(i, j + 1) = 1.0_dp
                    psi(i) = tiny(1.0_dp)
                end if
            end do
            sumtau_new = sum(tau_new, dim=1)
            if (any(.not. ieee_is_finite(tau_new)) .or. any(sumtau_new(2:) <= 0.0_dp)) then
                failed = .true.
                fit%flags(1) = .true.
                exit
            end if
            if (icd > 0.0_dp .and. sumtau_new(1) / real(n, dp) > npr_max) then
                call solve_noise_prior(icd, sumtau_old, pr, psi, npr_max, root, root_ok)
                if (.not. root_ok) then
                    failed = .true.
                    fit%flags(2) = .true.
                    exit
                end if
                fit%flags(3) = .true.
                do j = 1, g
                    if (pr(j + 1) > 0.0_dp) phi(:, j + 1) = phi(:, j + 1) / pr(j + 1)
                end do
                pr(1) = root
                pr(2:) = (1.0_dp - root) * sumtau_old(2:) / max(real(n, dp) - sumtau_old(1), tiny(1.0_dp))
                do j = 1, g
                    phi(:, j + 1) = pr(j + 1) * phi(:, j + 1)
                end do
                phi(:, 1) = pr(1) * icd
                psi = sum(phi, dim=2)
                do i = 1, n
                    if (psi(i) <= 0.0_dp) then
                        failed = .true.
                        fit%flags(2) = .true.
                        exit
                    end if
                    tau_new(i, :) = phi(i, :) / psi(i)
                end do
                if (failed) exit
                sumtau_new = sum(tau_new, dim=1)
                if (real(count(maxloc(tau_new, dim=2) == 1), dp) / real(n, dp) > npr_max + 0.01_dp) then
                    failed = .true.
                    fit%flags(2) = .true.
                    exit
                end if
            end if
            loglik_new = sum(log(max(psi, tiny(1.0_dp))))
            dloglik = abs(loglik_new - loglik_old)
            tau_old = tau_new
            sumtau_old = sumtau_new
            if (ieee_is_finite(dloglik) .and. dloglik < tol) then
                converged = .true.
                loglik_old = loglik_new
                exit
            end if
            loglik_old = loglik_new
        end do
        fit%g = g
        fit%logicd = logicd
        fit%iter = min(iter, iter_max)
        if (failed) then
            fit%code = 0
            return
        end if
        if (converged) then
            fit%code = 2
        else
            fit%code = 1
        end if
        pr = sumtau_old / real(n, dp)
        kd = 0.0_dp
        do j = 1, g
            kd(j) = weighted_chisq_ks(smd(:, j), tau_old(:, j + 1), p)
        end do
        fit%criterion = sum(kd * pr(2:)) / max(sum(pr(2:)), tiny(1.0_dp)) + beta * pr(1)
        fit%iloglik = loglik_old
        allocate(fit%pi(g + 1), fit%mean(p, g), fit%cov(p, p, g), fit%tau(n, g + 1), fit%smd(n, g))
        allocate(fit%exproportion(g + 1), fit%cluster(n), fit%size(g + 1))
        fit%pi = pr
        fit%mean = means
        fit%tau = tau_old
        fit%smd = smd
        fit%exproportion = sum(tau_old, dim=1) / real(n, dp)
        fit%cluster = maxloc(tau_old, dim=2) - 1
        fit%size = 0
        do i = 1, n
            fit%size(fit%cluster(i) + 1) = fit%size(fit%cluster(i) + 1) + 1
        end do
        do j = 1, g
            fit%cov(:, :, j) = matmul(vec(:, :, j), matmul(diag_matrix(eig(:, j)), transpose(vec(:, :, j))))
        end do
    end subroutine ecm_fit

    subroutine restrict_eigenvalues(values, sumtau, erc)
        real(dp), intent(inout) :: values(:, :) !! Component covariance eigenvalues; shape p by g.
        real(dp), intent(in) :: sumtau(:) !! Effective Gaussian-component membership totals; length g.
        real(dp), intent(in) :: erc !! Maximum allowed global ratio between any two eigenvalues.
        real(dp) :: a
        real(dp) :: b
        real(dp) :: f1
        real(dp) :: f2
        real(dp) :: gr
        real(dp) :: l
        real(dp) :: x1
        real(dp) :: x2
        integer :: iter

        a = log(max(minval(values), tiny(1.0_dp)))
        b = log(min(maxval(values), huge(1.0_dp)))
        if (b <= a) then
            values = max(values, tiny(1.0_dp))
            return
        end if
        gr = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp
        x1 = a + (1.0_dp - gr) * (b - a)
        x2 = a + gr * (b - a)
        f1 = restriction_objective(values, sumtau, erc, exp(x1))
        f2 = restriction_objective(values, sumtau, erc, exp(x2))
        do iter = 1, 120
            if (abs(b - a) <= sqrt(epsilon(1.0_dp))) exit
            if (f1 > f2) then
                a = x1
                x1 = x2
                f1 = f2
                x2 = a + gr * (b - a)
                f2 = restriction_objective(values, sumtau, erc, exp(x2))
            else
                b = x2
                x2 = x1
                f2 = f1
                x1 = a + (1.0_dp - gr) * (b - a)
                f1 = restriction_objective(values, sumtau, erc, exp(x1))
            end if
        end do
        l = exp(0.5_dp * (a + b))
        where (values < l) values = l
        where (values > erc * l) values = erc * l
    end subroutine restrict_eigenvalues

    real(dp) function restriction_objective(values, sumtau, erc, lower) result(obj)
        real(dp), intent(in) :: values(:, :) !! Original positive covariance eigenvalues; shape p by g.
        real(dp), intent(in) :: sumtau(:) !! Effective membership weights for each Gaussian component.
        real(dp), intent(in) :: erc !! Maximum eigenvalue ratio.
        real(dp), intent(in) :: lower !! Candidate lower eigenvalue truncation threshold.
        real(dp), allocatable :: clipped(:, :)
        integer :: j

        allocate(clipped(size(values, 1), size(values, 2)))
        clipped = values
        where (clipped < lower) clipped = lower
        where (clipped > erc * lower) clipped = erc * lower
        obj = 0.0_dp
        do j = 1, size(values, 2)
            obj = obj + sumtau(j) * sum(log(clipped(:, j)) + values(:, j) / clipped(:, j))
        end do
    end function restriction_objective

    subroutine solve_noise_prior(icd, sumtau_old, pr, psi, npr_max, root, ok)
        real(dp), intent(in) :: icd !! Positive improper density constant.
        real(dp), intent(in) :: sumtau_old(:) !! Previous posterior component totals, noise first.
        real(dp), intent(in) :: pr(:) !! Previous component prior probabilities, noise first.
        real(dp), intent(in) :: psi(:) !! Current mixture density at every observation.
        real(dp), intent(in) :: npr_max !! Requested upper bound on average posterior noise probability.
        real(dp), intent(out) :: root !! Updated noise mixing proportion satisfying the constraint.
        logical, intent(out) :: ok !! True when a bracketed root was found and bisected successfully.
        real(dp) :: a
        real(dp) :: b
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: fm
        real(dp) :: m
        integer :: iter

        a = 0.0_dp
        b = npr_max
        fa = noise_equation(a, icd, sumtau_old, pr, psi, npr_max)
        fb = noise_equation(b, icd, sumtau_old, pr, psi, npr_max)
        if (abs(fa) <= tiny(1.0_dp)) then
            root = a
            ok = .true.
            return
        end if
        if (abs(fb) <= tiny(1.0_dp)) then
            root = b
            ok = .true.
            return
        end if
        if (fa * fb > 0.0_dp) then
            root = b
            ok = .false.
            return
        end if
        do iter = 1, 200
            m = 0.5_dp * (a + b)
            fm = noise_equation(m, icd, sumtau_old, pr, psi, npr_max)
            if (abs(fm) <= 100.0_dp * epsilon(1.0_dp) * real(size(psi), dp)) exit
            if (fa * fm <= 0.0_dp) then
                b = m
                fb = fm
            else
                a = m
                fa = fm
            end if
        end do
        root = 0.5_dp * (a + b)
        ok = .true.
    end subroutine solve_noise_prior

    real(dp) function noise_equation(x, icd, sumtau_old, pr, psi, npr_max) result(f)
        real(dp), intent(in) :: x !! Candidate noise prior probability.
        real(dp), intent(in) :: icd !! Positive improper density constant.
        real(dp), intent(in) :: sumtau_old(:) !! Previous component posterior totals, noise first.
        real(dp), intent(in) :: pr(:) !! Previous component prior probabilities, noise first.
        real(dp), intent(in) :: psi(:) !! Current mixture densities.
        real(dp), intent(in) :: npr_max !! Target mean posterior noise probability.
        real(dp) :: den
        real(dp) :: gaussian_part
        integer :: i
        integer :: n

        n = size(psi)
        if (icd * x <= tiny(1.0_dp)) then
            f = -real(n, dp) * npr_max
            return
        end if
        f = -real(n, dp) * npr_max
        do i = 1, n
            gaussian_part = max(psi(i) - pr(1) * icd, 0.0_dp)
            den = icd * x + ((1.0_dp - x) / max(real(n, dp) - sumtau_old(1), tiny(1.0_dp))) * &
                gaussian_part * real(n, dp)
            f = f + icd * x / max(den, tiny(1.0_dp))
        end do
    end function noise_equation

    real(dp) function weighted_chisq_ks(x, weights, df) result(stat)
        real(dp), intent(in) :: x(:) !! Squared Mahalanobis distances for one Gaussian component.
        real(dp), intent(in) :: weights(:) !! Posterior component weights corresponding to x.
        integer, intent(in) :: df !! Chi-square degrees of freedom, normally the data dimension.
        integer, allocatable :: ord(:)
        real(dp) :: c
        real(dp) :: sw
        real(dp) :: d
        integer :: k

        sw = sum(weights)
        if (sw <= 0.0_dp) then
            stat = 0.0_dp
            return
        end if
        allocate(ord(size(x)))
        call order_real(x, ord)
        c = 0.0_dp
        stat = 0.0_dp
        do k = 1, size(x)
            c = c + weights(ord(k)) / sw
            d = abs(c - pchisq(x(ord(k)), real(df, dp)))
            stat = max(stat, d)
        end do
    end function weighted_chisq_ks

    subroutine default_logicd_grid(grid)
        real(dp), allocatable, intent(out) :: grid(:) !! Default sorted log-improper-density grid used by otrimle().
        real(dp), allocatable :: tmp(:)
        integer :: i
        integer :: k
        integer :: n

        n = 1 + 13 + 9 + 17 + 10
        allocate(tmp(n))
        k = 1
        tmp(k) = -huge(1.0_dp)
        do i = -700, -100, 50
            k = k + 1
            tmp(k) = real(i, dp)
        end do
        do i = -95, -55, 5
            k = k + 1
            tmp(k) = real(i, dp)
        end do
        do i = 0, 16
            k = k + 1
            tmp(k) = -50.0_dp + 2.5_dp * real(i, dp)
        end do
        do i = -9, 0
            k = k + 1
            tmp(k) = real(i, dp)
        end do
        allocate(grid(k))
        grid = tmp(1:k)
    end subroutine default_logicd_grid

    subroutine distances_from_row(x, row, distances)
        real(dp), intent(in) :: x(:, :) !! Observation matrix used to compute Euclidean distances.
        integer, intent(in) :: row !! One-based row index whose distances to all other rows are requested.
        real(dp), intent(out) :: distances(:) !! Euclidean distances to every row except row; length n-1.
        integer :: i
        integer :: k

        k = 0
        do i = 1, size(x, 1)
            if (i == row) cycle
            k = k + 1
            distances(k) = sqrt(sum((x(i, :) - x(row, :)) ** 2))
        end do
    end subroutine distances_from_row

    subroutine average_link_manhattan(x, g, labels)
        real(dp), intent(in) :: x(:, :) !! Regular observations to cluster; shape n by p.
        integer, intent(in) :: g !! Desired number of average-linkage clusters.
        integer, intent(out) :: labels(:) !! Cluster labels from 1 through g for the rows of x.
        integer, allocatable :: group(:)
        integer, allocatable :: uniq(:)
        integer :: a
        integer :: b
        integer :: best_a
        integer :: best_b
        integer :: i
        integer :: j
        integer :: n
        integer :: ng
        integer :: newlabel
        real(dp) :: best
        real(dp) :: d
        real(dp) :: sumd
        integer :: cnt

        n = size(x, 1)
        if (size(labels) /= n) error stop 'average_link_manhattan: label size mismatch'
        allocate(group(n), uniq(n))
        group = [(i, i=1,n)]
        ng = n
        do while (ng > g)
            call unique_positive(group, uniq, ng)
            best = huge(1.0_dp)
            best_a = uniq(1)
            best_b = uniq(2)
            do a = 1, ng - 1
                do b = a + 1, ng
                    sumd = 0.0_dp
                    cnt = 0
                    do i = 1, n
                        if (group(i) /= uniq(a)) cycle
                        do j = 1, n
                            if (group(j) /= uniq(b)) cycle
                            d = sum(abs(x(i, :) - x(j, :)))
                            sumd = sumd + d
                            cnt = cnt + 1
                        end do
                    end do
                    if (cnt > 0 .and. sumd / real(cnt, dp) < best) then
                        best = sumd / real(cnt, dp)
                        best_a = uniq(a)
                        best_b = uniq(b)
                    end if
                end do
            end do
            where (group == best_b) group = best_a
            call unique_positive(group, uniq, ng)
        end do
        call unique_positive(group, uniq, ng)
        labels = 0
        do a = 1, ng
            newlabel = a
            where (group == uniq(a)) labels = newlabel
        end do
    end subroutine average_link_manhattan

    subroutine deterministic_regular_partition(regular, reg_idx, data, g, labels)
        real(dp), intent(in) :: regular(:, :) !! kNN-retained observations eligible for non-noise cluster labels.
        integer, intent(in) :: reg_idx(:) !! Original row indices corresponding to regular.
        real(dp), intent(in) :: data(:, :) !! Full observation matrix used to validate distinct rows.
        integer, intent(in) :: g !! Desired positive number of clusters.
        integer, intent(out) :: labels(:) !! Labels for all observations; trimmed rows remain zero.
        integer, allocatable :: ord(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: shift

        n = size(regular, 1)
        allocate(ord(n))
        call order_real(regular(:, 1), ord)
        labels = 0
        do shift = 0, min(9, max(0, n - 1))
            labels = 0
            do i = 1, n
                j = 1 + modulo(i - 1 + shift, n)
                labels(reg_idx(ord(j))) = min(g, 1 + (i - 1) * g / n)
            end do
            if (valid_partition(data, labels, g)) return
        end do
    end subroutine deterministic_regular_partition

    logical function valid_partition(data, labels, g) result(ok)
        real(dp), intent(in) :: data(:, :) !! Observation matrix whose within-cluster distinct rows are checked.
        integer, intent(in) :: labels(:) !! Cluster labels; zero is ignored as noise.
        integer, intent(in) :: g !! Expected number of non-noise clusters.
        real(dp), allocatable :: subset(:, :)
        integer :: j
        integer :: m

        ok = .true.
        do j = 1, g
            m = count(labels == j)
            if (m < 2) then
                ok = .false.
                return
            end if
            allocate(subset(m, size(data, 2)))
            call pack_rows(data, labels == j, subset)
            if (count_unique_rows(subset) < 2) ok = .false.
            deallocate(subset)
            if (.not. ok) return
        end do
    end function valid_partition

    subroutine pack_rows(x, mask, out)
        real(dp), intent(in) :: x(:, :) !! Source matrix from which selected rows are copied.
        logical, intent(in) :: mask(:) !! Row-selection mask with length size(x,1).
        real(dp), intent(out) :: out(:, :) !! Packed selected rows, with count(mask) rows.
        integer :: i
        integer :: k

        k = 0
        do i = 1, size(x, 1)
            if (mask(i)) then
                k = k + 1
                out(k, :) = x(i, :)
            end if
        end do
    end subroutine pack_rows

    pure function outer(a, b) result(c)
        real(dp), intent(in) :: a(:) !! Left vector in the outer product.
        real(dp), intent(in) :: b(:) !! Right vector in the outer product.
        real(dp) :: c(size(a), size(b))

        c = spread(a, 2, size(b)) * spread(b, 1, size(a))
    end function outer

    pure function diag_matrix(x) result(a)
        real(dp), intent(in) :: x(:) !! Values placed on the diagonal of the returned square matrix.
        real(dp) :: a(size(x), size(x))
        integer :: i

        a = 0.0_dp
        do i = 1, size(x)
            a(i, i) = x(i)
        end do
    end function diag_matrix

    subroutine order_real(x, order)
        real(dp), intent(in) :: x(:) !! Real vector whose ascending-order permutation is requested.
        integer, intent(out) :: order(:) !! Permutation of 1:size(x) into ascending value order.
        integer :: i
        integer :: j
        integer :: m
        integer :: t

        do i = 1, size(x)
            order(i) = i
        end do
        do i = 1, size(x) - 1
            m = i
            do j = i + 1, size(x)
                if (x(order(j)) < x(order(m))) m = j
            end do
            if (m /= i) then
                t = order(i)
                order(i) = order(m)
                order(m) = t
            end if
        end do
    end subroutine order_real

    subroutine unique_positive(x, values, nvalues)
        integer, intent(in) :: x(:) !! Positive integer labels whose unique sorted values are requested.
        integer, intent(out) :: values(:) !! Workspace receiving unique labels.
        integer, intent(out) :: nvalues !! Number of unique values written to values.
        integer :: i

        nvalues = 0
        do i = 1, size(x)
            if (.not. any(values(1:nvalues) == x(i))) then
                nvalues = nvalues + 1
                values(nvalues) = x(i)
            end if
        end do
        call sort_int(values(1:nvalues))
    end subroutine unique_positive

    subroutine sort_int(x)
        integer, intent(inout) :: x(:) !! Integer vector sorted into ascending order in place.
        integer :: i
        integer :: j
        integer :: m
        integer :: t

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
    end subroutine sort_int

end module otrimle_core
