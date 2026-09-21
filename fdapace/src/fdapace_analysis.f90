module fdapace_analysis
    use fdapace_fpca, only : fpca_dense
    use fdapace_kinds, only : dp
    use fdapace_math, only : jacobi_eigen_sym, linear_interp, sort_eigen_descending
    use fdapace_numerics, only : column_correlation, empirical_quantile, rank_average
    use fdapace_smoothing, only : lwls1d, trapz_rcpp
    use fdapace_types, only : cluster_result, fpca_der_result, fpca_result, fvpa_result, stringing_result, wfda_result
    implicit none
    private

    public :: fclust
    public :: fpca_der
    public :: fvpa
    public :: k_cfc
    public :: stringing
    public :: wfda

contains

    subroutine fpca_der(fpca_obj, result, derivative_order, bandwidth, kernel_type)
        type(fpca_result), intent(in) :: fpca_obj !! Previously fitted FPCA object whose mean and eigenfunctions are differentiated.
        type(fpca_der_result), intent(out) :: result !! Derivatives of the fitted mean and eigenfunctions on the FPCA work grid.
        integer, intent(in), optional :: derivative_order !! Positive derivative order; defaults to one.
        real(dp), intent(in), optional :: bandwidth !! Positive local-polynomial bandwidth; defaults to 10 percent of the work-grid range.
        character(len=*), intent(in), optional :: kernel_type !! Smoothing kernel for local-polynomial differentiation; defaults to gauss.
        character(len=16) :: kernel
        real(dp), allocatable :: weights(:)
        real(dp) :: bw
        integer :: j
        integer :: p

        if (.not. allocated(fpca_obj%work_grid) .or. .not. allocated(fpca_obj%mu) .or. .not. allocated(fpca_obj%phi)) then
            error stop "fpca_der: incomplete FPCA object"
        end if
        p = 1
        if (present(derivative_order)) p = derivative_order
        if (p < 1) error stop "fpca_der: derivative_order must be positive"
        bw = 0.1_dp * (maxval(fpca_obj%work_grid) - minval(fpca_obj%work_grid))
        if (present(bandwidth)) bw = bandwidth
        if (bw <= 0.0_dp) error stop "fpca_der: bandwidth must be positive"
        kernel = "gauss"
        if (present(kernel_type)) kernel = trim(kernel_type)
        allocate(weights(size(fpca_obj%work_grid)))
        weights = 1.0_dp
        result%mu_der = lwls1d(bw, kernel, fpca_obj%work_grid, fpca_obj%mu, fpca_obj%work_grid, &
            weights, p, p)
        allocate(result%phi_der(size(fpca_obj%phi, 1), size(fpca_obj%phi, 2)))
        do j = 1, size(fpca_obj%phi, 2)
            result%phi_der(:, j) = lwls1d(bw, kernel, fpca_obj%work_grid, fpca_obj%phi(:, j), &
                fpca_obj%work_grid, weights, p, p)
        end do
        result%derivative_order = p
    end subroutine fpca_der

    subroutine fvpa(y, t, result, q, fve_threshold)
        real(dp), intent(in) :: y(:,:) !! Dense regular functional observations used for variance-process analysis.
        real(dp), intent(in) :: t(:) !! Common regular observation grid corresponding to columns of y.
        type(fvpa_result), intent(out) :: result !! FPCA fits for the signal and log variance process, plus residual adjustment and variance.
        real(dp), intent(in), optional :: q !! Residual-square quantile used as the positive log adjustment; defaults to 0.1.
        real(dp), intent(in), optional :: fve_threshold !! FVE threshold used in both dense FPCA fits; defaults to 0.9.
        real(dp), allocatable :: residual2(:,:)
        real(dp), allocatable :: flat(:)
        real(dp), allocatable :: log_res(:,:)
        real(dp) :: prob
        real(dp) :: threshold

        if (size(y, 2) /= size(t)) error stop "fvpa: y columns must match t"
        prob = 0.1_dp
        if (present(q)) prob = q
        if (prob < 0.0_dp .or. prob > 1.0_dp) prob = 0.1_dp
        threshold = 0.9_dp
        if (present(fve_threshold)) threshold = fve_threshold
        call fpca_dense(y, t, result%fpca_y, fve_threshold=threshold, assume_error=.true.)
        allocate(residual2(size(y, 1), size(y, 2)))
        residual2 = (y - result%fpca_y%fitted_y)**2
        flat = reshape(residual2, [size(residual2)])
        result%delta = empirical_quantile(flat, prob)
        allocate(log_res(size(y, 1), size(y, 2)))
        log_res = log(max(result%delta + residual2, tiny(1.0_dp)))
        call fpca_dense(log_res, t, result%fpca_r, fve_threshold=threshold, assume_error=.true.)
        result%sigma2 = result%fpca_r%sigma2
    end subroutine fvpa

    subroutine stringing(x, result, y, standardize, distance, dis_mat)
        real(dp), intent(in) :: x(:,:) !! Input data matrix with observations in rows and variables to be ordered in columns.
        type(stringing_result), intent(out) :: result !! Variable ordering, ordered matrix, optional standardized matrix, and dissimilarities.
        real(dp), intent(in), optional :: y(:) !! Optional scalar response required by the xycor distance.
        logical, intent(in), optional :: standardize !! If true, standardize each variable before computing distances; defaults to false.
        character(len=*), intent(in), optional :: distance !! Dissimilarity name; supports euclidean, maximum, manhattan, canberra, correlation, spearman, hamming, xycor, and user.
        real(dp), intent(in), optional :: dis_mat(:,:) !! User-supplied symmetric dissimilarity matrix for distance='user'.
        character(len=16) :: metric
        logical :: do_standardize
        real(dp), allocatable :: b(:,:)
        real(dp), allocatable :: col_mean(:)
        real(dp), allocatable :: d(:,:)
        real(dp), allocatable :: eval(:)
        real(dp), allocatable :: evec(:,:)
        real(dp), allocatable :: row_mean(:)
        real(dp), allocatable :: rx(:)
        real(dp), allocatable :: ry(:)
        real(dp), allocatable :: work(:,:)
        real(dp) :: grand_mean
        real(dp) :: scale
        integer :: i
        integer :: info
        integer :: j
        integer :: n
        integer :: p

        n = size(x, 1)
        p = size(x, 2)
        if (n < 2 .or. p < 1) error stop "stringing: X must contain at least two observations"
        do_standardize = .false.
        if (present(standardize)) do_standardize = standardize
        metric = "euclidean"
        if (present(distance)) metric = trim(distance)
        allocate(work(n, p))
        work = x
        if (do_standardize) then
            do j = 1, p
                scale = sqrt(sum((work(:, j) - sum(work(:, j)) / real(n, dp))**2) / real(n - 1, dp))
                if (scale <= 0.0_dp) error stop "stringing: cannot standardize a constant column"
                work(:, j) = (work(:, j) - sum(work(:, j)) / real(n, dp)) / scale
            end do
            result%standardized_x = work
            result%standardized = .true.
        end if
        allocate(d(p, p))
        if (trim(metric) == "user") then
            if (.not. present(dis_mat)) error stop "stringing: dis_mat is required for user distance"
            if (size(dis_mat, 1) /= p .or. size(dis_mat, 2) /= p) error stop "stringing: dis_mat shape mismatch"
            d = dis_mat
        else
            d = 0.0_dp
            if (trim(metric) == "spearman") allocate(rx(n), ry(n))
            do j = 1, p
                do i = j + 1, p
                    select case (trim(metric))
                    case ("euclidean")
                        d(i, j) = sqrt(sum((work(:, i) - work(:, j))**2))
                    case ("maximum")
                        d(i, j) = maxval(abs(work(:, i) - work(:, j)))
                    case ("manhattan")
                        d(i, j) = sum(abs(work(:, i) - work(:, j)))
                    case ("canberra")
                        d(i, j) = canberra_distance(work(:, i), work(:, j))
                    case ("correlation")
                        d(i, j) = sqrt(max(0.0_dp, 2.0_dp * (1.0_dp - column_correlation(work(:, i), work(:, j)))))
                    case ("spearman")
                        call rank_average(work(:, i), rx)
                        call rank_average(work(:, j), ry)
                        d(i, j) = 1.0_dp - column_correlation(rx, ry)
                    case ("hamming", "binary")
                        d(i, j) = real(count(work(:, i) /= work(:, j)), dp) / real(n, dp)
                    case ("minkowski")
                        d(i, j) = sum(abs(work(:, i) - work(:, j))**2)**0.5_dp
                    case ("xycor")
                        if (.not. present(y)) error stop "stringing: y is required for xycor distance"
                        if (size(y) /= n) error stop "stringing: y length mismatch"
                        d(i, j) = abs(column_correlation(work(:, i), y) - column_correlation(work(:, j), y))
                    case default
                        error stop "stringing: unsupported distance"
                    end select
                    d(j, i) = d(i, j)
                end do
            end do
        end if
        result%distance = d
        allocate(row_mean(p), col_mean(p), b(p, p), eval(p), evec(p, p))
        row_mean = sum(d**2, dim=2) / real(p, dp)
        col_mean = sum(d**2, dim=1) / real(p, dp)
        grand_mean = sum(d**2) / real(p * p, dp)
        do j = 1, p
            do i = 1, p
                b(i, j) = -0.5_dp * (d(i, j)**2 - row_mean(i) - col_mean(j) + grand_mean)
            end do
        end do
        call jacobi_eigen_sym(b, eval, evec, info)
        if (info /= 0) error stop "stringing: classical MDS eigensolver failed"
        call sort_eigen_descending(eval, evec)
        allocate(result%order(p))
        call order_real(evec(:, 1), result%order)
        allocate(result%stringed_x(n, p))
        result%stringed_x = work(:, result%order)
    end subroutine stringing

    pure real(dp) function canberra_distance(x, y) result(value)
        real(dp), intent(in) :: x(:) !! First vector in a Canberra-distance calculation.
        real(dp), intent(in) :: y(:) !! Second vector of the same length in a Canberra-distance calculation.
        real(dp) :: den
        integer :: i

        if (size(x) /= size(y)) error stop "canberra_distance: shape mismatch"
        value = 0.0_dp
        do i = 1, size(x)
            den = abs(x(i)) + abs(y(i))
            if (den > 0.0_dp) value = value + abs(x(i) - y(i)) / den
        end do
    end function canberra_distance

    pure subroutine order_real(x, order)
        real(dp), intent(in) :: x(:) !! Values whose one-based ascending permutation is requested.
        integer, intent(out) :: order(:) !! Ascending one-based permutation of x.
        integer :: i
        integer :: j
        integer :: key

        if (size(order) /= size(x)) error stop "order_real: shape mismatch"
        order = [(i, i=1, size(x))]
        do i = 2, size(x)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (x(order(j)) <= x(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine order_real

    subroutine fclust(y, t, k, result, max_iter)
        real(dp), intent(in) :: y(:,:) !! Dense regular functional sample with subjects in rows.
        real(dp), intent(in) :: t(:) !! Common regular observation grid corresponding to columns of y.
        integer, intent(in) :: k !! Number of requested clusters.
        type(cluster_result), intent(out) :: result !! K-means clustering of dense FPCA scores.
        integer, intent(in), optional :: max_iter !! Maximum k-means iterations; defaults to 125.
        type(fpca_result) :: fit
        integer :: itmax

        if (k < 1 .or. k > size(y, 1)) error stop "fclust: invalid k"
        itmax = 125
        if (present(max_iter)) itmax = max_iter
        call fpca_dense(y, t, fit, assume_error=.false.)
        call kmeans_scores(fit%xi_est, k, itmax, result)
    end subroutine fclust

    subroutine k_cfc(y, t, k, result, max_iter)
        real(dp), intent(in) :: y(:,:) !! Dense regular functional sample with subjects in rows.
        real(dp), intent(in) :: t(:) !! Common regular observation grid corresponding to columns of y.
        integer, intent(in) :: k !! Number of functional clusters.
        type(cluster_result), intent(out) :: result !! Iterative cluster assignment based on cluster-specific FPCA reconstruction error.
        integer, intent(in), optional :: max_iter !! Maximum reassignment iterations; defaults to 125.
        type(cluster_result) :: initial
        type(fpca_result), allocatable :: fits(:)
        integer, allocatable :: assign(:)
        integer, allocatable :: old_assign(:)
        integer, allocatable :: idx(:)
        real(dp), allocatable :: costs(:)
        real(dp), allocatable :: subset(:,:)
        integer :: c
        integer :: i
        integer :: iter
        integer :: itmax
        integer :: n

        n = size(y, 1)
        if (k < 1 .or. k > n) error stop "k_cfc: invalid k"
        itmax = 125
        if (present(max_iter)) itmax = max_iter
        call fclust(y, t, k, initial, max_iter=min(25, itmax))
        assign = initial%cluster
        allocate(old_assign(n), fits(k), costs(k))
        do iter = 1, itmax
            old_assign = assign
            do c = 1, k
                idx = pack([(i, i=1, n)], assign == c)
                if (size(idx) < 2) then
                    result%cluster = assign
                    result%iterations = iter - 1
                    result%converged = .false.
                    return
                end if
                allocate(subset(size(idx), size(y, 2)))
                subset = y(idx, :)
                call fpca_dense(subset, t, fits(c), fve_threshold=0.7_dp, assume_error=.false.)
                deallocate(subset)
            end do
            do i = 1, n
                do c = 1, k
                    costs(c) = fpca_reconstruction_cost(y(i, :), t, fits(c))
                end do
                assign(i) = minloc(costs, dim=1)
            end do
            if (all(assign == old_assign)) exit
        end do
        result%cluster = assign
        result%iterations = iter
        result%converged = all(assign == old_assign)
        allocate(result%centers(k, size(initial%centers, 2)))
        result%centers = 0.0_dp
    end subroutine k_cfc

    subroutine kmeans_scores(scores, k, max_iter, result)
        real(dp), intent(in) :: scores(:,:) !! Feature matrix with observations in rows for deterministic MacQueen-style clustering.
        integer, intent(in) :: k !! Number of clusters.
        integer, intent(in) :: max_iter !! Positive maximum number of Lloyd reassignment iterations.
        type(cluster_result), intent(out) :: result !! Cluster labels, centers, iteration count, and convergence flag.
        integer, allocatable :: assign(:)
        integer, allocatable :: old_assign(:)
        integer :: c
        integer :: i
        integer :: iter
        integer :: n_in
        real(dp) :: best
        real(dp) :: dist

        if (k < 1 .or. k > size(scores, 1) .or. max_iter < 1) error stop "kmeans_scores: invalid controls"
        allocate(result%centers(k, size(scores, 2)), assign(size(scores, 1)), old_assign(size(scores, 1)))
        do c = 1, k
            i = 1 + (c - 1) * max(1, (size(scores, 1) - 1) / max(1, k - 1))
            i = min(i, size(scores, 1))
            result%centers(c, :) = scores(i, :)
        end do
        assign = 1
        old_assign = 0
        do iter = 1, max_iter
            old_assign = assign
            do i = 1, size(scores, 1)
                assign(i) = 1
                best = sum((scores(i, :) - result%centers(1, :))**2)
                do c = 2, k
                    dist = sum((scores(i, :) - result%centers(c, :))**2)
                    if (dist < best) then
                        best = dist
                        assign(i) = c
                    end if
                end do
            end do
            do c = 1, k
                n_in = count(assign == c)
                if (n_in > 0) then
                    result%centers(c, :) = 0.0_dp
                    do i = 1, size(scores, 1)
                        if (assign(i) == c) result%centers(c, :) = result%centers(c, :) + scores(i, :)
                    end do
                    result%centers(c, :) = result%centers(c, :) / real(n_in, dp)
                end if
            end do
            if (all(assign == old_assign)) exit
        end do
        result%cluster = assign
        result%iterations = iter
        result%converged = all(assign == old_assign)
    end subroutine kmeans_scores

    real(dp) function fpca_reconstruction_cost(y, t, fit) result(cost)
        real(dp), intent(in) :: y(:) !! New dense curve on the same grid as fit used for reconstruction-error scoring.
        real(dp), intent(in) :: t(:) !! Common support grid used for numerical integration.
        type(fpca_result), intent(in) :: fit !! Cluster-specific FPCA fit defining mean and eigenfunctions.
        real(dp), allocatable :: recon(:)
        real(dp) :: score
        integer :: j

        if (size(y) /= size(t) .or. size(y) /= size(fit%mu)) error stop "fpca_reconstruction_cost: shape mismatch"
        allocate(recon(size(y)))
        recon = fit%mu
        do j = 1, size(fit%phi, 2)
            score = trapz_rcpp(t, (y - fit%mu) * fit%phi(:, j))
            recon = recon + score * fit%phi(:, j)
        end do
        cost = trapz_rcpp(t, (y - recon)**2)
    end function fpca_reconstruction_cost

    subroutine wfda(y, t, result, lambda, subset_prop)
        real(dp), intent(in) :: y(:,:) !! Dense common-grid curves to align by monotone one-parameter warp averaging.
        real(dp), intent(in) :: t(:) !! Strictly increasing common support, rescaled internally to [0,1].
        type(wfda_result), intent(out) :: result !! Estimated warps, inverse warps, aligned curves, costs, and regularization weight.
        real(dp), intent(in), optional :: lambda !! Nonnegative warp regularization weight; defaults to 1e-4 times sample scale.
        real(dp), intent(in), optional :: subset_prop !! Retained compatibility argument in (0,1]; this dense approximation uses all curves.
        real(dp), allocatable :: grid(:)
        real(dp), allocatable :: mean_curve(:)
        real(dp), allocatable :: norm_y(:,:)
        real(dp), allocatable :: warp(:)
        real(dp) :: a
        real(dp) :: best_a
        real(dp) :: best_cost
        real(dp) :: cost
        real(dp) :: lam
        real(dp) :: max_abs
        real(dp) :: prop
        integer :: ia
        integer :: i
        integer :: j
        integer, parameter :: na = 61

        if (size(y, 2) /= size(t) .or. size(t) < 2) error stop "wfda: invalid dense dimensions"
        prop = 0.5_dp
        if (present(subset_prop)) prop = subset_prop
        if (prop <= 0.0_dp .or. prop > 1.0_dp) error stop "wfda: subset_prop must lie in (0,1]"
        allocate(grid(size(t)), norm_y(size(y, 1), size(y, 2)), mean_curve(size(t)))
        grid = (t - t(1)) / (t(size(t)) - t(1))
        do i = 1, size(y, 1)
            max_abs = maxval(abs(y(i, :)))
            if (max_abs <= 0.0_dp) then
                norm_y(i, :) = y(i, :)
            else
                norm_y(i, :) = y(i, :) / max_abs
            end if
        end do
        mean_curve = sum(norm_y, dim=1) / real(size(y, 1), dp)
        lam = 1.0e-4_dp * sqrt(sum((norm_y - spread(mean_curve, 1, size(y, 1)))**2) / &
            max(1.0_dp, real(size(y, 1) - 1, dp)))
        if (present(lambda)) lam = lambda
        if (lam < 0.0_dp) error stop "wfda: lambda must be nonnegative"
        result%lambda = lam
        allocate(result%h(size(y, 1), size(t)), result%h_inv(size(y, 1), size(t)))
        allocate(result%aligned(size(y, 1), size(t)), result%costs(size(y, 1)), warp(size(t)))
        do i = 1, size(y, 1)
            best_cost = huge(1.0_dp)
            best_a = 1.0_dp
            do ia = 1, na
                a = 0.5_dp + 1.5_dp * real(ia - 1, dp) / real(na - 1, dp)
                warp = grid**a
                cost = 0.0_dp
                do j = 1, size(t)
                    cost = cost + (linear_interp(grid, norm_y(i, :), warp(j), .true.) - mean_curve(j))**2
                end do
                cost = cost + lam * sum((warp - grid)**2)
                if (cost < best_cost) then
                    best_cost = cost
                    best_a = a
                end if
            end do
            result%h(i, :) = grid**best_a
            result%h_inv(i, :) = grid**(1.0_dp / best_a)
            do j = 1, size(t)
                result%aligned(i, j) = linear_interp(grid, y(i, :), result%h(i, j), .true.)
            end do
            result%costs(i) = best_cost / real(size(t), dp)
        end do
    end subroutine wfda

end module fdapace_analysis
