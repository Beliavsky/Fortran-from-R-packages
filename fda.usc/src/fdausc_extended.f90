module fdausc_extended
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use fdausc_integration, only : integrate_curve
    use fdausc_kernels, only : kernel_value, ker_norm
    use fdausc_metrics, only : lp_distance_matrix, combine_distance_matrices, distance_correlation
    use fdausc_depth, only : mahalanobis_depth
    use fdausc_smoothing, only : gcv_score
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: basis_project_grid, basis_reconstruct_grid, select_basis_gcv
    public :: default_bandwidth_grid, horizontal_shift_distance, horizontal_shift_distance_matrix
    public :: fourier_semimetric, basis_semimetric
    public :: likelihood_depth, bivariate_simplicial_depth
    public :: rpd_depth_scores, multicomponent_mahalanobis_depth, multicomponent_modal_depth
    public :: depth_classify, select_classification_cv
    public :: bootstrap_mean_replicates, depth_outlier_flags, lrt_outlier_statistic
    public :: local_distance_correlation_select, hetero_anova_onefactor

contains

    subroutine basis_project_grid(curves, basis, coefficients, mean_curve, center, info)
        real(dp), intent(in) :: curves(:, :) !! Functional observations by rows on a common grid.
        real(dp), intent(in) :: basis(:, :) !! Basis functions by rows evaluated on the same grid as curves.
        real(dp), intent(out) :: coefficients(:, :) !! Least-squares basis coefficients, one observation per row.
        real(dp), intent(out) :: mean_curve(:) !! Mean curve removed before projection when centering is enabled.
        logical, intent(in), optional :: center !! If true, subtract the sample mean before basis projection; default true.
        integer, intent(out) :: info !! Zero on success; positive when the basis Gram matrix is singular.
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: invgram(:, :)
        real(dp), allocatable :: centered(:, :)
        real(dp), allocatable :: rhs(:)
        logical :: do_center
        integer :: i

        do_center = .true.
        if (present(center)) do_center = center
        mean_curve = 0.0_dp
        if (do_center) mean_curve = sum(curves, dim=1)/real(size(curves, 1), dp)
        allocate(centered(size(curves, 1), size(curves, 2)))
        centered = curves - spread(mean_curve, 1, size(curves, 1))
        if (.not. do_center) centered = curves
        gram = matmul(basis, transpose(basis))
        call inverse_matrix(gram, invgram, info)
        if (info /= 0) then
            coefficients = 0.0_dp
            return
        end if
        allocate(rhs(size(basis, 1)))
        do i = 1, size(curves, 1)
            rhs = matmul(basis, centered(i, :))
            coefficients(i, :) = matmul(invgram, rhs)
        end do
    end subroutine basis_project_grid

    pure subroutine basis_reconstruct_grid(coefficients, basis, mean_curve, curves)
        real(dp), intent(in) :: coefficients(:, :) !! Basis coefficients, one functional observation per row.
        real(dp), intent(in) :: basis(:, :) !! Basis functions by rows evaluated on the target grid.
        real(dp), intent(in) :: mean_curve(:) !! Mean curve added to each reconstruction; use zeros for no centering.
        real(dp), intent(out) :: curves(:, :) !! Reconstructed functional observations on the basis grid.

        curves = matmul(coefficients, basis) + spread(mean_curve, 1, size(coefficients, 1))
    end subroutine basis_reconstruct_grid

    pure subroutine select_basis_gcv(curves, smoothers, scores, best_index)
        real(dp), intent(in) :: curves(:, :) !! Functional observations by rows used to evaluate each candidate smoother.
        real(dp), intent(in) :: smoothers(:, :, :) !! Candidate square smoothing matrices in the third dimension.
        real(dp), intent(out) :: scores(:) !! GCV score for each candidate smoothing matrix.
        integer, intent(out) :: best_index !! One-based index of the candidate with minimum GCV score.
        real(dp) :: df
        real(dp) :: one_score
        integer :: i
        integer :: j

        do i = 1, size(smoothers, 3)
            scores(i) = 0.0_dp
            do j = 1, size(curves, 1)
                call gcv_score(curves(j, :), smoothers(:, :, i), one_score, df)
                scores(i) = scores(i) + one_score
            end do
            scores(i) = scores(i)/real(size(curves, 1), dp)
        end do
        best_index = minloc(scores, dim=1)
    end subroutine select_basis_gcv

    pure subroutine default_bandwidth_grid(distances, prob_min, prob_max, bandwidths, kernel_code)
        real(dp), intent(in) :: distances(:, :) !! Pairwise distance matrix whose diagonal is ignored.
        real(dp), intent(in) :: prob_min !! Lower type-4 empirical quantile probability in [0,1].
        real(dp), intent(in) :: prob_max !! Upper type-4 empirical quantile probability in [0,1].
        real(dp), intent(out) :: bandwidths(:) !! Positive bandwidth sequence spanning the selected quantiles.
        integer, intent(in), optional :: kernel_code !! Kernel selector matching fdausc_kernels; default normal.
        real(dp), allocatable :: vals(:)
        real(dp) :: qmin
        real(dp) :: qmax
        real(dp) :: c0
        integer :: code
        integer :: i
        integer :: j
        integer :: k
        integer :: n

        n = size(distances, 1)*size(distances, 2) - min(size(distances, 1), size(distances, 2))
        allocate(vals(max(1, n)))
        k = 0
        do i = 1, size(distances, 1)
            do j = 1, size(distances, 2)
                if (i == j .and. size(distances, 1) == size(distances, 2)) cycle
                if (distances(i, j) < huge(1.0_dp) .and. distances(i, j) >= 0.0_dp) then
                    k = k + 1
                    vals(k) = distances(i, j)
                end if
            end do
        end do
        if (k == 0) then
            bandwidths = 0.0_dp
            return
        end if
        call sort_real_local(vals(:k))
        qmin = quantile_type4_sorted(vals(:k), max(0.0_dp, min(1.0_dp, prob_min)))
        qmax = quantile_type4_sorted(vals(:k), max(0.0_dp, min(1.0_dp, prob_max)))
        qmax = max(qmin, qmax)
        code = 1
        if (present(kernel_code)) code = kernel_code
        c0 = kernel_bandwidth_scale(code)
        if (size(bandwidths) == 1) then
            bandwidths(1) = qmax*c0
        else
            do i = 1, size(bandwidths)
                bandwidths(i) = (qmin + real(i - 1, dp)*(qmax - qmin)/real(size(bandwidths) - 1, dp))*c0
            end do
        end if
    end subroutine default_bandwidth_grid

    pure subroutine horizontal_shift_distance(x, y, t, distance, lag)
        real(dp), intent(in) :: x(:) !! First discretized curve.
        real(dp), intent(in) :: y(:) !! Second discretized curve on the same grid.
        real(dp), intent(in) :: t(:) !! Increasing grid used to define the maximum horizontal lag.
        real(dp), intent(out) :: distance !! Minimum root integrated squared distance after horizontal shifting.
        integer, intent(out) :: lag !! Nonnegative grid lag attaining the reported minimum distance.
        real(dp) :: rang
        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: best
        real(dp) :: denom
        integer :: lagmax
        integer :: l
        integer :: m

        rang = t(size(t)) - t(1)
        if (rang <= 0.0_dp .or. size(x) < 2) then
            distance = 0.0_dp
            lag = 0
            return
        end if
        lagmax = min(size(x) - 2, max(0, floor(0.2_dp*rang)))
        best = sum((x(1:size(x) - 1) - y(1:size(y) - 1))**2 + &
                   (x(2:size(x)) - y(2:size(y)))**2)/(2.0_dp*rang)
        lag = 0
        do l = 1, lagmax
            m = size(x) - l
            denom = 2.0_dp*(rang - 2.0_dp*real(l, dp))
            if (denom <= 0.0_dp .or. m < 2) cycle
            d1 = sum((x(l + 1:size(x) - 1) - y(1:m - 1))**2 + &
                     (x(l + 2:size(x)) - y(2:m))**2)/denom
            d2 = sum((x(1:m - 1) - y(l + 1:size(y) - 1))**2 + &
                     (x(2:m) - y(l + 2:size(y)))**2)/denom
            if (d1 < best) then
                best = d1
                lag = l
            end if
            if (d2 < best) then
                best = d2
                lag = l
            end if
        end do
        distance = sqrt(max(0.0_dp, best))
    end subroutine horizontal_shift_distance

    pure subroutine horizontal_shift_distance_matrix(curves1, curves2, t, distances, lags)
        real(dp), intent(in) :: curves1(:, :) !! First functional sample, one curve per row.
        real(dp), intent(in) :: curves2(:, :) !! Second functional sample on the same grid.
        real(dp), intent(in) :: t(:) !! Increasing grid shared by both samples.
        real(dp), intent(out) :: distances(:, :) !! Pairwise horizontal-shift distances.
        integer, intent(out), optional :: lags(:, :) !! Optional pairwise optimal nonnegative grid lags.
        integer :: i
        integer :: j
        integer :: lag0

        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                call horizontal_shift_distance(curves1(i, :), curves2(j, :), t, distances(i, j), lag0)
                if (present(lags)) lags(i, j) = lag0
            end do
        end do
    end subroutine horizontal_shift_distance_matrix

    subroutine fourier_semimetric(t, curves1, curves2, nbasis, nderiv, period, distances, info)
        real(dp), intent(in) :: t(:) !! Functional grid used for Fourier fitting and derivative evaluation.
        real(dp), intent(in) :: curves1(:, :) !! First functional sample, one curve per row.
        real(dp), intent(in) :: curves2(:, :) !! Second functional sample, one curve per row.
        integer, intent(in) :: nbasis !! Number of Fourier basis functions retained; must be positive.
        integer, intent(in) :: nderiv !! Nonnegative derivative order defining the semimetric.
        real(dp), intent(in), optional :: period !! Positive Fourier period; default is grid range.
        real(dp), intent(out) :: distances(:, :) !! Pairwise L2 distances between fitted derivative curves.
        integer, intent(out) :: info !! Zero on success; positive if a basis projection is singular.
        real(dp), allocatable :: basis(:, :)
        real(dp), allocatable :: dbasis(:, :)
        real(dp), allocatable :: coef1(:, :)
        real(dp), allocatable :: coef2(:, :)
        real(dp), allocatable :: dcurves1(:, :)
        real(dp), allocatable :: dcurves2(:, :)
        real(dp), allocatable :: mean1(:)
        real(dp), allocatable :: mean2(:)
        real(dp) :: per
        integer :: stat

        per = t(size(t)) - t(1)
        if (present(period)) per = period
        if (per <= 0.0_dp .or. nbasis < 1) then
            distances = 0.0_dp
            info = 1
            return
        end if
        allocate(basis(nbasis, size(t)), dbasis(nbasis, size(t)))
        call fourier_basis_values(t, nbasis, per, 0, basis)
        call fourier_basis_values(t, nbasis, per, nderiv, dbasis)
        allocate(coef1(size(curves1, 1), nbasis), coef2(size(curves2, 1), nbasis))
        allocate(mean1(size(t)), mean2(size(t)))
        call basis_project_grid(curves1, basis, coef1, mean1, center=.false., info=stat)
        if (stat /= 0) then
            distances = 0.0_dp
            info = stat
            return
        end if
        call basis_project_grid(curves2, basis, coef2, mean2, center=.false., info=stat)
        if (stat /= 0) then
            distances = 0.0_dp
            info = stat
            return
        end if
        allocate(dcurves1(size(curves1, 1), size(t)), dcurves2(size(curves2, 1), size(t)))
        dcurves1 = matmul(coef1, dbasis)
        dcurves2 = matmul(coef2, dbasis)
        call lp_distance_matrix(t, dcurves1, dcurves2, distances, p=2.0_dp, method=2)
        info = 0
    end subroutine fourier_semimetric

    pure subroutine basis_semimetric(t, coefficients1, coefficients2, derivative_basis, distances)
        real(dp), intent(in) :: t(:) !! Increasing grid on which derivative basis functions are supplied.
        real(dp), intent(in) :: coefficients1(:, :) !! Basis coefficients for the first functional sample.
        real(dp), intent(in) :: coefficients2(:, :) !! Basis coefficients for the second functional sample.
        real(dp), intent(in) :: derivative_basis(:, :) !! Requested derivative of each basis function by rows.
        real(dp), intent(out) :: distances(:, :) !! Pairwise L2 semimetric between reconstructed derivatives.
        real(dp), allocatable :: x1(:, :)
        real(dp), allocatable :: x2(:, :)

        x1 = matmul(coefficients1, derivative_basis)
        x2 = matmul(coefficients2, derivative_basis)
        call lp_distance_matrix(t, x1, x2, distances, p=2.0_dp, method=2)
    end subroutine basis_semimetric

    pure subroutine likelihood_depth(distances_eval_ref, bandwidth, depth, scale, distances_ref_ref)
        real(dp), intent(in) :: distances_eval_ref(:, :) !! Distances from evaluated points to the reference sample.
        real(dp), intent(in) :: bandwidth !! Positive Gaussian-kernel bandwidth in distance units.
        real(dp), intent(out) :: depth(:) !! Likelihood-depth scores proportional to Gaussian kernel sums.
        logical, intent(in), optional :: scale !! If true, normalize by the maximum reference-sample kernel sum.
        real(dp), intent(in), optional :: distances_ref_ref(:, :) !! Reference pairwise distances used for scaling.
        real(dp), allocatable :: ref_depth(:)
        real(dp) :: denom
        integer :: i

        do i = 1, size(distances_eval_ref, 1)
            depth(i) = sum(ker_norm(distances_eval_ref(i, :)/bandwidth))
        end do
        if (present(scale)) then
            if (scale) then
                if (present(distances_ref_ref)) then
                    allocate(ref_depth(size(distances_ref_ref, 1)))
                    do i = 1, size(distances_ref_ref, 1)
                        ref_depth(i) = sum(ker_norm(distances_ref_ref(i, :)/bandwidth))
                    end do
                    denom = maxval(ref_depth)
                else
                    denom = maxval(depth)
                end if
                if (denom > tiny(1.0_dp)) depth = depth/denom
            end if
        end if
    end subroutine likelihood_depth

    pure subroutine bivariate_simplicial_depth(x, reference, depth, scale, same_sample)
        real(dp), intent(in) :: x(:, :) !! Evaluated two-dimensional points, one point per row.
        real(dp), intent(in) :: reference(:, :) !! Two-dimensional reference sample defining simplices.
        real(dp), intent(out) :: depth(:) !! Simplicial depth following the upstream vertex-inclusion convention.
        logical, intent(in), optional :: scale !! If true, normalize the resulting depths by their maximum.
        logical, intent(in), optional :: same_sample !! True when x and reference are the same aligned sample.
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: nref
        integer :: count_inside
        integer :: vertex_count
        integer :: denom_count
        logical :: aligned
        logical :: do_scale

        nref = size(reference, 1)
        aligned = .false.
        if (present(same_sample)) aligned = same_sample
        do i = 1, size(x, 1)
            count_inside = 0
            do j = 1, nref - 2
                if (aligned .and. i == j) cycle
                do k = j + 1, nref - 1
                    if (aligned .and. i == k) cycle
                    do l = k + 1, nref
                        if (aligned .and. i == l) cycle
                        if (point_in_triangle(x(i, :), reference(j, :), reference(k, :), reference(l, :))) then
                            count_inside = count_inside + 1
                        end if
                    end do
                end do
            end do
            if (aligned) then
                vertex_count = (nref - 1)*(nref - 2)/2
                denom_count = nref*(nref - 1)*(nref - 2)/6
            else
                vertex_count = nref*(nref - 1)/2
                denom_count = (nref + 1)*nref*(nref - 1)/6
            end if
            if (denom_count > 0) then
                depth(i) = real(count_inside + vertex_count, dp)/real(denom_count, dp)
            else
                depth(i) = 0.0_dp
            end if
        end do
        do_scale = .false.
        if (present(scale)) do_scale = scale
        if (do_scale .and. maxval(depth) > tiny(1.0_dp)) depth = depth/maxval(depth)
    end subroutine bivariate_simplicial_depth

    subroutine rpd_depth_scores(projected_eval, projected_ref, depth, bandwidth_probability, scale)
        real(dp), intent(in) :: projected_eval(:, :, :) !! Evaluation projection scores [case,feature,projection].
        real(dp), intent(in) :: projected_ref(:, :, :) !! Reference projection scores [case,feature,projection].
        real(dp), intent(out) :: depth(:) !! Average likelihood depth across supplied projections.
        real(dp), intent(in), optional :: bandwidth_probability !! Reference distance quantile; default 0.15.
        logical, intent(in), optional :: scale !! If true, scale likelihood depth within each projection.
        real(dp), allocatable :: deval(:, :)
        real(dp), allocatable :: dref(:, :)
        real(dp), allocatable :: depj(:)
        real(dp), allocatable :: offdiag(:)
        real(dp) :: prob
        real(dp) :: h
        logical :: do_scale
        integer :: i
        integer :: j
        integer :: k
        integer :: p
        integer :: n

        prob = 0.15_dp
        if (present(bandwidth_probability)) prob = bandwidth_probability
        do_scale = .true.
        if (present(scale)) do_scale = scale
        allocate(deval(size(projected_eval, 1), size(projected_ref, 1)))
        allocate(dref(size(projected_ref, 1), size(projected_ref, 1)))
        allocate(depj(size(projected_eval, 1)))
        n = size(projected_ref, 1)*(size(projected_ref, 1) - 1)
        allocate(offdiag(max(1, n)))
        depth = 0.0_dp
        do p = 1, size(projected_eval, 3)
            do i = 1, size(projected_eval, 1)
                do j = 1, size(projected_ref, 1)
                    deval(i, j) = sqrt(sum((projected_eval(i, :, p) - projected_ref(j, :, p))**2))
                end do
            end do
            do i = 1, size(projected_ref, 1)
                do j = 1, size(projected_ref, 1)
                    dref(i, j) = sqrt(sum((projected_ref(i, :, p) - projected_ref(j, :, p))**2))
                end do
            end do
            k = 0
            do i = 1, size(dref, 1)
                do j = 1, size(dref, 2)
                    if (i == j) cycle
                    k = k + 1
                    offdiag(k) = dref(i, j)
                end do
            end do
            call sort_real_local(offdiag(:k))
            h = quantile_type4_sorted(offdiag(:k), max(0.0_dp, min(1.0_dp, prob)))
            h = max(h, sqrt(tiny(1.0_dp)))
            call likelihood_depth(deval, h, depj, scale=do_scale, distances_ref_ref=dref)
            depth = depth + depj
        end do
        depth = depth/real(size(projected_eval, 3), dp)
    end subroutine rpd_depth_scores

    subroutine multicomponent_mahalanobis_depth(values_eval, values_ref, depth, info)
        real(dp), intent(in) :: values_eval(:, :, :) !! Evaluation data [case,grid,component].
        real(dp), intent(in) :: values_ref(:, :, :) !! Reference data [case,grid,component].
        real(dp), intent(out) :: depth(:) !! Mean pointwise multivariate Mahalanobis depth.
        integer, intent(out) :: info !! Zero if all pointwise covariance inversions succeed; otherwise first failure code.
        real(dp), allocatable :: depj(:)
        real(dp), allocatable :: xe(:, :)
        real(dp), allocatable :: xr(:, :)
        integer :: j
        integer :: stat

        allocate(depj(size(values_eval, 1)))
        allocate(xe(size(values_eval, 1), size(values_eval, 3)))
        allocate(xr(size(values_ref, 1), size(values_ref, 3)))
        depth = 0.0_dp
        info = 0
        do j = 1, size(values_eval, 2)
            xe = values_eval(:, j, :)
            xr = values_ref(:, j, :)
            call mahalanobis_depth(xe, xr, depj, stat)
            if (stat /= 0 .and. info == 0) info = stat
            depth = depth + depj
        end do
        depth = depth/real(size(values_eval, 2), dp)
    end subroutine multicomponent_mahalanobis_depth

    pure subroutine multicomponent_modal_depth(distance_components, component_weights, method, bandwidth, depth, scale)
        real(dp), intent(in) :: distance_components(:, :, :) !! Component eval-to-reference distance matrices.
        real(dp), intent(in) :: component_weights(:) !! Nonnegative component distance weights.
        integer, intent(in) :: method !! Combination code used by combine_distance_matrices.
        real(dp), intent(in) :: bandwidth !! Positive bandwidth for the combined modal-depth kernel.
        real(dp), intent(out) :: depth(:) !! Modal depth based on the combined component distance.
        logical, intent(in), optional :: scale !! If true, normalize by the largest evaluated depth score.
        real(dp), allocatable :: combined(:, :)
        real(dp) :: denom
        integer :: i

        allocate(combined(size(distance_components, 1), size(distance_components, 2)))
        call combine_distance_matrices(distance_components, component_weights, method, combined)
        do i = 1, size(combined, 1)
            depth(i) = sum(ker_norm(combined(i, :)/bandwidth))
        end do
        if (present(scale)) then
            if (scale) then
                denom = maxval(depth)
                if (denom > tiny(1.0_dp)) depth = depth/denom
            end if
        end if
    end subroutine multicomponent_modal_depth

    pure subroutine depth_classify(depth_by_class, predicted, probabilities)
        real(dp), intent(in) :: depth_by_class(:, :) !! Depth score of each case relative to each class reference sample.
        integer, intent(out) :: predicted(:) !! Class with maximum depth for each case.
        real(dp), intent(out), optional :: probabilities(:, :) !! Row-normalized nonnegative depth scores by class.
        real(dp) :: s
        integer :: i

        do i = 1, size(depth_by_class, 1)
            predicted(i) = maxloc(depth_by_class(i, :), dim=1)
            if (present(probabilities)) then
                probabilities(i, :) = max(depth_by_class(i, :), 0.0_dp)
                s = sum(probabilities(i, :))
                if (s > tiny(1.0_dp)) probabilities(i, :) = probabilities(i, :)/s
            end if
        end do
    end subroutine depth_classify

    pure subroutine select_classification_cv(truth, candidate_predictions, best_index, error_rates, costs)
        integer, intent(in) :: truth(:) !! Observed positive integer class labels for cross-validation cases.
        integer, intent(in) :: candidate_predictions(:, :) !! Predictions [case,candidate] from completed CV fits.
        integer, intent(out) :: best_index !! One-based candidate index with minimum weighted error.
        real(dp), intent(out) :: error_rates(:) !! Misclassification rate for each candidate.
        real(dp), intent(in), optional :: costs(:) !! Optional per-class misclassification costs indexed by true class.
        real(dp) :: wsum
        integer :: i
        integer :: j

        do j = 1, size(candidate_predictions, 2)
            error_rates(j) = 0.0_dp
            wsum = 0.0_dp
            do i = 1, size(truth)
                if (present(costs)) then
                    if (truth(i) >= 1 .and. truth(i) <= size(costs)) then
                        wsum = wsum + costs(truth(i))
                        if (candidate_predictions(i, j) /= truth(i)) error_rates(j) = error_rates(j) + costs(truth(i))
                    end if
                else
                    wsum = wsum + 1.0_dp
                    if (candidate_predictions(i, j) /= truth(i)) error_rates(j) = error_rates(j) + 1.0_dp
                end if
            end do
            if (wsum > 0.0_dp) error_rates(j) = error_rates(j)/wsum
        end do
        best_index = minloc(error_rates, dim=1)
    end subroutine select_classification_cv

    pure subroutine bootstrap_mean_replicates(t, curves, bootstrap_indices, noise, alpha, center, replicates, distances, dband)
        real(dp), intent(in) :: t(:) !! Increasing functional grid used for L2 bootstrap distances.
        real(dp), intent(in) :: curves(:, :) !! Original functional sample, one curve per row.
        integer, intent(in) :: bootstrap_indices(:, :) !! Resampled row indices [sample_position,replicate].
        real(dp), intent(in) :: noise(:, :, :) !! Optional pre-generated smoothing noise [sample_position,grid,replicate].
        real(dp), intent(in) :: alpha !! Tail probability defining the confidence-ball radius.
        real(dp), intent(out) :: center(:) !! Original sample mean curve.
        real(dp), intent(out) :: replicates(:, :) !! Bootstrap mean curve for each replicate by rows.
        real(dp), intent(out) :: distances(:) !! L2 distance of each bootstrap mean from the original mean.
        real(dp), intent(out) :: dband !! Empirical (1-alpha) bootstrap confidence-ball radius.
        real(dp), allocatable :: sorted(:)
        real(dp), allocatable :: diff(:)
        integer :: b
        integer :: i
        integer :: idx
        integer :: k

        center = sum(curves, dim=1)/real(size(curves, 1), dp)
        allocate(diff(size(t)))
        do b = 1, size(bootstrap_indices, 2)
            replicates(b, :) = 0.0_dp
            do i = 1, size(bootstrap_indices, 1)
                idx = bootstrap_indices(i, b)
                if (idx >= 1 .and. idx <= size(curves, 1)) then
                    replicates(b, :) = replicates(b, :) + curves(idx, :) + noise(i, :, b)
                end if
            end do
            replicates(b, :) = replicates(b, :)/real(size(bootstrap_indices, 1), dp)
            diff = replicates(b, :) - center
            distances(b) = sqrt(max(0.0_dp, integrate_curve(t, diff*diff, 2)))
        end do
        sorted = distances
        call sort_real_local(sorted)
        k = max(1, min(size(sorted), floor((1.0_dp - alpha)*real(size(sorted), dp))))
        dband = sorted(k)
    end subroutine bootstrap_mean_replicates

    pure subroutine depth_outlier_flags(depth, cutoff, flags)
        real(dp), intent(in) :: depth(:) !! Depth score for each observation.
        real(dp), intent(in) :: cutoff !! Observations with depth below this cutoff are flagged as outliers.
        logical, intent(out) :: flags(:) !! True for observations below the supplied depth cutoff.

        flags = depth < cutoff
    end subroutine depth_outlier_flags

    pure subroutine lrt_outlier_statistic(t, curves, trimmed_mean, trimmed_variance, statistics, maximum)
        real(dp), intent(in) :: t(:) !! Increasing grid used for standardized L2 distances.
        real(dp), intent(in) :: curves(:, :) !! Curves to score for the likelihood-ratio-style outlier statistic.
        real(dp), intent(in) :: trimmed_mean(:) !! Robust center curve used by the statistic.
        real(dp), intent(in) :: trimmed_variance(:) !! Positive pointwise robust variance curve.
        real(dp), intent(out) :: statistics(:) !! Standardized L2 distance for each curve.
        real(dp), intent(out) :: maximum !! Largest outlier statistic across the supplied curves.
        real(dp), allocatable :: z(:)
        integer :: i

        allocate(z(size(t)))
        do i = 1, size(curves, 1)
            z = (curves(i, :) - trimmed_mean)/sqrt(max(trimmed_variance, tiny(1.0_dp)))
            statistics(i) = sqrt(max(0.0_dp, integrate_curve(t, z*z, 2)))
        end do
        maximum = maxval(statistics)
    end subroutine lrt_outlier_statistic

    pure subroutine local_distance_correlation_select(x, y_distance, tolerance, correlations, selected)
        real(dp), intent(in) :: x(:, :) !! Candidate scalar covariates in columns.
        real(dp), intent(in) :: y_distance(:, :) !! Pairwise distance matrix for the scalar or multivariate response.
        real(dp), intent(in) :: tolerance !! Minimum distance correlation required for impact-point selection.
        real(dp), intent(out) :: correlations(:) !! Distance correlation between each covariate and the response.
        logical, intent(out) :: selected(:) !! True at local maxima whose distance correlation exceeds tolerance.
        real(dp), allocatable :: dx(:, :)
        integer :: i
        integer :: j
        integer :: n

        n = size(x, 1)
        allocate(dx(n, n))
        do j = 1, size(x, 2)
            do i = 1, n
                dx(i, :) = abs(x(i, j) - x(:, j))
            end do
            correlations(j) = distance_correlation(dx, y_distance)
        end do
        selected = .false.
        if (size(correlations) == 1) then
            selected(1) = correlations(1) > tolerance
            return
        end if
        if (correlations(1) > correlations(2) .and. correlations(1) > tolerance) selected(1) = .true.
        do j = 2, size(correlations) - 1
            if (correlations(j) >= correlations(j - 1) .and. correlations(j) > correlations(j + 1) .and. &
                correlations(j) > tolerance) selected(j) = .true.
        end do
        if (correlations(size(correlations)) >= correlations(size(correlations) - 1) .and. &
            correlations(size(correlations)) > tolerance) selected(size(correlations)) = .true.
    end subroutine local_distance_correlation_select

    pure subroutine hetero_anova_onefactor(y, groups, statistic, df1, df2)
        real(dp), intent(in) :: y(:) !! Scalar response values for a one-factor heteroscedastic ANOVA.
        integer, intent(in) :: groups(:) !! Positive integer factor level for each response observation.
        real(dp), intent(out) :: statistic !! Brunner-Dette-Munk box-type statistic for the one-factor contrast.
        real(dp), intent(out) :: df1 !! Numerator degrees-of-freedom approximation.
        real(dp), intent(out) :: df2 !! Denominator degrees-of-freedom approximation.
        integer :: g
        integer :: i
        integer :: j
        integer :: k
        integer, allocatable :: counts(:)
        real(dp), allocatable :: means(:)
        real(dp), allocatable :: sig(:)
        real(dp), allocatable :: mm(:, :)
        real(dp) :: grand
        real(dp) :: tr_dsn
        real(dp) :: tr_msn2
        real(dp) :: tr_d2sn2delta
        real(dp) :: nall

        g = maxval(groups)
        allocate(counts(g), means(g), sig(g), mm(g, g))
        counts = 0
        means = 0.0_dp
        sig = 0.0_dp
        do i = 1, size(y)
            if (groups(i) >= 1 .and. groups(i) <= g) then
                counts(groups(i)) = counts(groups(i)) + 1
                means(groups(i)) = means(groups(i)) + y(i)
            end if
        end do
        do j = 1, g
            if (counts(j) > 0) means(j) = means(j)/real(counts(j), dp)
        end do
        do i = 1, size(y)
            if (groups(i) >= 1 .and. groups(i) <= g) then
                j = groups(i)
                sig(j) = sig(j) + (y(i) - means(j))**2
            end if
        end do
        do j = 1, g
            if (counts(j) > 1) then
                sig(j) = sig(j)/real(counts(j) - 1, dp)
                sig(j) = sig(j)*real(counts(j) - 1, dp)/real(counts(j)*counts(j), dp)
            else
                sig(j) = 0.0_dp
            end if
        end do
        mm = -1.0_dp/real(g, dp)
        do j = 1, g
            mm(j, j) = 1.0_dp - 1.0_dp/real(g, dp)
        end do
        nall = real(sum(counts), dp)
        grand = dot_product(means, matmul(mm, means))
        tr_dsn = 0.0_dp
        tr_msn2 = 0.0_dp
        tr_d2sn2delta = 0.0_dp
        do j = 1, g
            tr_dsn = tr_dsn + mm(j, j)*nall*sig(j)
            do k = 1, g
                tr_msn2 = tr_msn2 + mm(j, k)*nall*sig(k)*mm(k, j)*nall*sig(j)
            end do
            if (counts(j) > 1) then
                tr_d2sn2delta = tr_d2sn2delta + mm(j, j)**2*(nall*sig(j))**2/real(counts(j) - 1, dp)
            end if
        end do
        if (tr_dsn <= tiny(1.0_dp)) then
            statistic = 0.0_dp
            df1 = 0.0_dp
            df2 = 0.0_dp
            return
        end if
        statistic = nall*grand/tr_dsn
        if (tr_msn2 > tiny(1.0_dp)) then
            df1 = tr_dsn*tr_dsn/tr_msn2
        else
            df1 = 0.0_dp
        end if
        if (tr_d2sn2delta > tiny(1.0_dp)) then
            df2 = tr_dsn*tr_dsn/tr_d2sn2delta
        else
            df2 = 0.0_dp
        end if
    end subroutine hetero_anova_onefactor

    pure subroutine fourier_basis_values(t, nbasis, period, nderiv, values)
        real(dp), intent(in) :: t(:) !! Grid at which the orthonormal Fourier basis is evaluated.
        integer, intent(in) :: nbasis !! Number of Fourier basis functions requested.
        real(dp), intent(in) :: period !! Positive Fourier period.
        integer, intent(in) :: nderiv !! Nonnegative derivative order.
        real(dp), intent(out) :: values(:, :) !! Basis values [basis,grid] for the requested derivative order.
        real(dp) :: omega
        real(dp) :: phase
        real(dp) :: factor
        integer :: b
        integer :: k

        values = 0.0_dp
        if (nderiv == 0) values(1, :) = 1.0_dp/sqrt(period)
        do b = 2, nbasis
            k = (b - 1 + 1)/2
            omega = 2.0_dp*pi*real(k, dp)/period
            factor = sqrt(2.0_dp/period)*omega**nderiv
            phase = 0.5_dp*pi*real(nderiv, dp)
            if (mod(b, 2) == 0) then
                values(b, :) = factor*sin(omega*(t - t(1)) + phase)
            else
                values(b, :) = factor*cos(omega*(t - t(1)) + phase)
            end if
        end do
    end subroutine fourier_basis_values

    pure logical function point_in_triangle(p, a, b, c) result(inside)
        real(dp), intent(in) :: p(:) !! Two-coordinate point whose triangle membership is tested.
        real(dp), intent(in) :: a(:) !! First triangle vertex.
        real(dp), intent(in) :: b(:) !! Second triangle vertex.
        real(dp), intent(in) :: c(:) !! Third triangle vertex.
        real(dp) :: d1
        real(dp) :: d2
        real(dp) :: d3
        real(dp) :: tol
        logical :: has_neg
        logical :: has_pos

        tol = 100.0_dp*epsilon(1.0_dp)
        d1 = cross2(p - a, b - a)
        d2 = cross2(p - b, c - b)
        d3 = cross2(p - c, a - c)
        has_neg = d1 < -tol .or. d2 < -tol .or. d3 < -tol
        has_pos = d1 > tol .or. d2 > tol .or. d3 > tol
        inside = .not. (has_neg .and. has_pos)
    end function point_in_triangle

    pure real(dp) function cross2(a, b) result(c)
        real(dp), intent(in) :: a(:) !! First two-dimensional vector.
        real(dp), intent(in) :: b(:) !! Second two-dimensional vector.

        c = a(1)*b(2) - a(2)*b(1)
    end function cross2

    pure real(dp) function quantile_type4_sorted(sorted, probability) result(q)
        real(dp), intent(in) :: sorted(:) !! Increasing sample values.
        real(dp), intent(in) :: probability !! Quantile probability in [0,1].
        real(dp) :: h
        real(dp) :: g
        integer :: j
        integer :: n

        n = size(sorted)
        if (probability <= 0.0_dp) then
            q = sorted(1)
            return
        end if
        if (probability >= 1.0_dp) then
            q = sorted(n)
            return
        end if
        h = real(n, dp)*probability
        j = floor(h)
        g = h - real(j, dp)
        if (j <= 0) then
            q = sorted(1)
        else if (j >= n) then
            q = sorted(n)
        else
            q = (1.0_dp - g)*sorted(j) + g*sorted(j + 1)
        end if
    end function quantile_type4_sorted

    pure real(dp) function kernel_bandwidth_scale(kernel_code) result(c0)
        integer, intent(in) :: kernel_code !! Kernel selector matching fdausc_kernels.
        integer, parameter :: nstep = 2000
        real(dp) :: ck
        real(dp) :: dk
        real(dp) :: dx
        real(dp) :: u
        real(dp) :: w
        real(dp) :: kval
        integer :: i

        ck = 0.0_dp
        dk = 0.0_dp
        dx = 8.0_dp/real(nstep, dp)
        do i = 0, nstep
            u = -4.0_dp + real(i, dp)*dx
            if (i == 0 .or. i == nstep) then
                w = 0.5_dp
            else
                w = 1.0_dp
            end if
            kval = kernel_value(u, kernel_code)
            ck = ck + w*kval*kval
            dk = dk + w*u*u*kval
        end do
        ck = ck*dx
        dk = dk*dx
        if (dk > tiny(1.0_dp)) then
            c0 = (ck/(dk*dk))**0.2_dp/1.3510_dp
        else
            c0 = 1.0_dp
        end if
    end function kernel_bandwidth_scale

    pure subroutine sort_real_local(x)
        real(dp), intent(inout) :: x(:) !! Real array sorted in ascending order in place.
        real(dp) :: key
        integer :: i
        integer :: j

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
    end subroutine sort_real_local

end module fdausc_extended
