module fdausc_models
    use r_kinds, only : dp
    use r_linalg, only : thin_svd, solve_system
    use fdausc_kernels, only : aker_epa, iker_epa, kernel_value
    use fdausc_integration, only : integrate_curve
    implicit none
    private

    public :: functional_pca, functional_pls1
    public :: conditional_cdf, conditional_quantile, conditional_mode
    public :: functional_kmeans, knn_classify, kernel_classify
    public :: pca_score_distance, derivative_distance

contains

    subroutine functional_pca(t, curves, ncomp, components, scores, singular_values, mean_curve, info)
        real(dp), intent(in) :: t(:) !! Increasing functional grid shared by all input curves.
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row.
        integer, intent(in) :: ncomp !! Number of leading principal components requested.
        real(dp), intent(out) :: components(:, :) !! Leading component functions in rows, shape [ncomp,size(t)].
        real(dp), intent(out) :: scores(:, :) !! Principal-component scores, shape [size(curves,1),ncomp].
        real(dp), intent(out) :: singular_values(:) !! Leading singular values after functional-grid normalization.
        real(dp), intent(out) :: mean_curve(:) !! Pointwise sample mean removed before the decomposition.
        integer, intent(out), optional :: info !! Zero on successful SVD; nonzero value propagated from rfortran-linalg.
        real(dp), allocatable :: centered(:, :)
        real(dp), allocatable :: u(:, :)
        real(dp), allocatable :: sval(:)
        real(dp), allocatable :: vt(:, :)
        real(dp) :: delta
        real(dp) :: normv
        integer :: i
        integer :: j
        integer :: nc
        integer :: stat

        mean_curve = sum(curves, dim=1)/real(size(curves, 1), dp)
        allocate(centered(size(curves, 1), size(curves, 2)))
        do i = 1, size(curves, 1)
            centered(i, :) = curves(i, :) - mean_curve
        end do
        call thin_svd(centered, u, sval, vt, stat)
        if (present(info)) info = stat
        components = 0.0_dp
        scores = 0.0_dp
        singular_values = 0.0_dp
        if (stat /= 0) return
        nc = min(ncomp, min(size(vt, 1), size(scores, 2)))
        if (size(t) > 1) then
            delta = sqrt((t(size(t)) - t(1))/real(size(t) - 1, dp))
        else
            delta = 1.0_dp
        end if
        do j = 1, nc
            components(j, :) = vt(j, :)
            normv = sqrt(max(0.0_dp, integrate_curve(t, components(j, :)**2, 2)))
            if (normv > sqrt(tiny(1.0_dp))) components(j, :) = components(j, :)/normv
            singular_values(j) = sval(j)*delta
            do i = 1, size(curves, 1)
                scores(i, j) = integrate_curve(t, centered(i, :)*components(j, :), 2)
            end do
        end do
    end subroutine functional_pca

    subroutine functional_pls1(curves, y, ncomp, weights, scores, loadings, coefficients, xmean, ymean)
        real(dp), intent(in) :: curves(:, :) !! Predictor matrix or discretized functional sample, observations in rows.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with curve rows.
        integer, intent(in) :: ncomp !! Number of orthogonal PLS1 components requested.
        real(dp), intent(out) :: weights(:, :) !! PLS predictor weight vectors in columns.
        real(dp), intent(out) :: scores(:, :) !! PLS component scores in columns.
        real(dp), intent(out) :: loadings(:, :) !! Predictor loadings corresponding to the PLS scores.
        real(dp), intent(out) :: coefficients(:, :) !! Regression coefficient vector after each component, columns 1:ncomp.
        real(dp), intent(out) :: xmean(:) !! Predictor column means removed before PLS fitting.
        real(dp), intent(out) :: ymean !! Response mean removed before PLS fitting.
        real(dp), allocatable :: xres(:, :)
        real(dp), allocatable :: yres(:)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: tscore(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: wstar(:, :)
        real(dp), allocatable :: ptw(:, :)
        real(dp), allocatable :: invcoef(:)
        real(dp) :: wn
        real(dp) :: denom
        integer :: h
        integer :: i
        integer :: nc
        integer :: stat

        xmean = sum(curves, dim=1)/real(size(curves, 1), dp)
        ymean = sum(y)/real(size(y), dp)
        allocate(xres(size(curves, 1), size(curves, 2)), yres(size(y)))
        do i = 1, size(curves, 1)
            xres(i, :) = curves(i, :) - xmean
        end do
        yres = y - ymean
        weights = 0.0_dp
        scores = 0.0_dp
        loadings = 0.0_dp
        coefficients = 0.0_dp
        nc = min(ncomp, min(size(weights, 2), size(scores, 2)))
        allocate(w(size(curves, 2)), tscore(size(curves, 1)), p(size(curves, 2)), q(nc))
        do h = 1, nc
            w = matmul(transpose(xres), yres)
            wn = sqrt(dot_product(w, w))
            if (wn <= sqrt(tiny(1.0_dp))) exit
            w = w/wn
            tscore = matmul(xres, w)
            denom = dot_product(tscore, tscore)
            if (denom <= tiny(1.0_dp)) exit
            p = matmul(transpose(xres), tscore)/denom
            q(h) = dot_product(yres, tscore)/denom
            weights(:, h) = w
            scores(:, h) = tscore
            loadings(:, h) = p
            xres = xres - spread(tscore, 2, size(xres, 2))*spread(p, 1, size(xres, 1))
            yres = yres - q(h)*tscore
            allocate(ptw(h, h), invcoef(h), wstar(size(curves, 2), h))
            ptw = matmul(transpose(loadings(:, :h)), weights(:, :h))
            call solve_system(transpose(ptw), q(:h), invcoef, stat)
            if (stat == 0) then
                wstar = weights(:, :h)
                coefficients(:, h) = matmul(wstar, invcoef)
            end if
            deallocate(ptw, invcoef, wstar)
        end do
    end subroutine functional_pls1

    pure subroutine conditional_cdf(distances_ref_eval, y, y0, h, g, fc)
        real(dp), intent(in) :: distances_ref_eval(:, :) !! Reference-to-evaluation distances; reference rows align with y.
        real(dp), intent(in) :: y(:) !! Scalar reference responses paired with rows of distances_ref_eval.
        real(dp), intent(in) :: y0(:) !! Response grid values at which conditional CDFs are evaluated.
        real(dp), intent(in) :: h !! Positive predictor-distance bandwidth.
        real(dp), intent(in) :: g !! Positive response bandwidth in response units.
        real(dp), intent(out) :: fc(:, :) !! Conditional CDF values with shape [size(y0),size(distances_ref_eval,2)].
        real(dp), allocatable :: w(:)
        real(dp) :: denom
        integer :: i
        integer :: j
        integer :: k

        allocate(w(size(y)))
        do j = 1, size(distances_ref_eval, 2)
            do k = 1, size(y)
                w(k) = aker_epa(distances_ref_eval(k, j)/h)
            end do
            denom = sum(w)
            if (denom <= tiny(1.0_dp)) then
                fc(:, j) = 0.0_dp
                cycle
            end if
            do i = 1, size(y0)
                fc(i, j) = 0.0_dp
                do k = 1, size(y)
                    fc(i, j) = fc(i, j) + w(k)*(1.0_dp - iker_epa((y(k) - y0(i))/g))
                end do
                fc(i, j) = fc(i, j)/denom
            end do
        end do
    end subroutine conditional_cdf

    pure subroutine conditional_quantile(y_grid, fc, probability, quantiles)
        real(dp), intent(in) :: y_grid(:) !! Strictly increasing response grid associated with rows of fc.
        real(dp), intent(in) :: fc(:, :) !! Conditional CDF values, one evaluation case per column.
        real(dp), intent(in) :: probability !! Target conditional probability in [0,1].
        real(dp), intent(out) :: quantiles(:) !! Linearly interpolated conditional quantile for each fc column.
        integer :: i
        integer :: j
        real(dp) :: f

        do j = 1, size(fc, 2)
            if (probability <= fc(1, j)) then
                quantiles(j) = y_grid(1)
            else if (probability >= fc(size(y_grid), j)) then
                quantiles(j) = y_grid(size(y_grid))
            else
                quantiles(j) = y_grid(size(y_grid))
                do i = 2, size(y_grid)
                    if (fc(i, j) >= probability) then
                        if (abs(fc(i, j) - fc(i - 1, j)) <= tiny(1.0_dp)) then
                            quantiles(j) = y_grid(i)
                        else
                            f = (probability - fc(i - 1, j))/(fc(i, j) - fc(i - 1, j))
                            quantiles(j) = (1.0_dp - f)*y_grid(i - 1) + f*y_grid(i)
                        end if
                        exit
                    end if
                end do
            end if
        end do
    end subroutine conditional_quantile

    pure subroutine conditional_mode(y_grid, fc, modes, density)
        real(dp), intent(in) :: y_grid(:) !! Strictly increasing response grid associated with rows of fc.
        real(dp), intent(in) :: fc(:, :) !! Conditional CDF values, one evaluation case per column.
        real(dp), intent(out) :: modes(:) !! Grid locations of the maximum finite-difference conditional density.
        real(dp), intent(out) :: density(:, :) !! First-difference density approximation between successive CDF grid points.
        integer :: i
        integer :: j
        integer :: idx

        density = 0.0_dp
        do j = 1, size(fc, 2)
            do i = 2, size(y_grid)
                density(i, j) = (fc(i, j) - fc(i - 1, j))/(y_grid(i) - y_grid(i - 1))
            end do
            density(1, j) = density(2, j)
            idx = maxloc(density(:, j), dim=1)
            modes(j) = y_grid(idx)
        end do
    end subroutine conditional_mode

    pure subroutine functional_kmeans(t, curves, k, max_iter, centers, labels, inertia)
        real(dp), intent(in) :: t(:) !! Increasing functional grid used in L2 assignment distances.
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row.
        integer, intent(in) :: k !! Requested number of clusters, between one and number of curves.
        integer, intent(in) :: max_iter !! Maximum Lloyd iterations.
        real(dp), intent(out) :: centers(:, :) !! Cluster-center curves, shape [k,size(t)].
        integer, intent(out) :: labels(:) !! One-based cluster assignment for each input curve.
        real(dp), intent(out) :: inertia !! Sum of squared functional L2 distances to assigned centers.
        real(dp), allocatable :: new_centers(:, :)
        integer, allocatable :: counts(:)
        integer, allocatable :: old_labels(:)
        real(dp) :: dist
        real(dp) :: best
        integer :: i
        integer :: j
        integer :: iter
        integer :: bestj

        do j = 1, k
            centers(j, :) = curves(1 + (j - 1)*size(curves, 1)/k, :)
        end do
        labels = 0
        allocate(new_centers(k, size(t)), counts(k), old_labels(size(labels)))
        do iter = 1, max_iter
            old_labels = labels
            do i = 1, size(curves, 1)
                best = huge(1.0_dp)
                bestj = 1
                do j = 1, k
                    dist = integrate_curve(t, (curves(i, :) - centers(j, :))**2, 2)
                    if (dist < best) then
                        best = dist
                        bestj = j
                    end if
                end do
                labels(i) = bestj
            end do
            if (iter > 1 .and. all(labels == old_labels)) exit
            new_centers = 0.0_dp
            counts = 0
            do i = 1, size(curves, 1)
                new_centers(labels(i), :) = new_centers(labels(i), :) + curves(i, :)
                counts(labels(i)) = counts(labels(i)) + 1
            end do
            do j = 1, k
                if (counts(j) > 0) centers(j, :) = new_centers(j, :)/real(counts(j), dp)
            end do
        end do
        inertia = 0.0_dp
        do i = 1, size(curves, 1)
            inertia = inertia + integrate_curve(t, (curves(i, :) - centers(labels(i), :))**2, 2)
        end do
    end subroutine functional_kmeans

    pure subroutine knn_classify(dist_train_test, train_labels, k, predicted, probabilities)
        real(dp), intent(in) :: dist_train_test(:, :) !! Distances from training observations (rows) to test observations (columns).
        integer, intent(in) :: train_labels(:) !! Positive integer class labels aligned with distance rows.
        integer, intent(in) :: k !! Number of nearest training observations voting for each test case.
        integer, intent(out) :: predicted(:) !! Predicted one-based class label for each test observation.
        real(dp), intent(out) :: probabilities(:, :) !! Class vote fractions with shape [nclass,ntest].
        real(dp), allocatable :: d(:)
        integer, allocatable :: idx(:)
        integer :: i
        integer :: j
        integer :: c
        integer :: kk
        integer :: nclass

        nclass = size(probabilities, 1)
        kk = max(1, min(k, size(train_labels)))
        allocate(d(size(train_labels)), idx(size(train_labels)))
        do j = 1, size(dist_train_test, 2)
            d = dist_train_test(:, j)
            do i = 1, size(idx)
                idx(i) = i
            end do
            call sort_index_by_value(d, idx)
            probabilities(:, j) = 0.0_dp
            do i = 1, kk
                c = train_labels(idx(i))
                if (c >= 1 .and. c <= nclass) probabilities(c, j) = probabilities(c, j) + 1.0_dp
            end do
            probabilities(:, j) = probabilities(:, j)/real(kk, dp)
            predicted(j) = maxloc(probabilities(:, j), dim=1)
        end do
    end subroutine knn_classify

    pure subroutine kernel_classify(dist_train_test, train_labels, bandwidth, kernel_code, predicted, probabilities, weights)
        real(dp), intent(in) :: dist_train_test(:, :) !! Distances from training observations (rows) to test observations (columns).
        integer, intent(in) :: train_labels(:) !! Positive integer class labels aligned with distance rows.
        real(dp), intent(in) :: bandwidth !! Positive kernel bandwidth for distance weighting.
        integer, intent(in), optional :: kernel_code !! Symmetric kernel selector; default normal.
        integer, intent(out) :: predicted(:) !! Predicted one-based class label for each test observation.
        real(dp), intent(out) :: probabilities(:, :) !! Normalized class kernel weights with shape [nclass,ntest].
        real(dp), intent(in), optional :: weights(:) !! Optional training-observation weights.
        real(dp) :: w
        real(dp) :: total
        integer :: code
        integer :: i
        integer :: j
        integer :: c

        code = 1
        if (present(kernel_code)) code = kernel_code
        do j = 1, size(dist_train_test, 2)
            probabilities(:, j) = 0.0_dp
            do i = 1, size(train_labels)
                w = kernel_value(dist_train_test(i, j)/bandwidth, code)
                if (present(weights)) w = w*weights(i)
                c = train_labels(i)
                if (c >= 1 .and. c <= size(probabilities, 1)) probabilities(c, j) = probabilities(c, j) + w
            end do
            total = sum(probabilities(:, j))
            if (total > tiny(1.0_dp)) probabilities(:, j) = probabilities(:, j)/total
            predicted(j) = maxloc(probabilities(:, j), dim=1)
        end do
    end subroutine kernel_classify

    pure subroutine pca_score_distance(scores1, scores2, distances)
        real(dp), intent(in) :: scores1(:, :) !! First matrix of retained functional principal-component scores.
        real(dp), intent(in) :: scores2(:, :) !! Second matrix of retained functional principal-component scores.
        real(dp), intent(out) :: distances(:, :) !! Euclidean distances between score rows.
        integer :: i
        integer :: j

        do i = 1, size(scores1, 1)
            do j = 1, size(scores2, 1)
                distances(i, j) = sqrt(sum((scores1(i, :) - scores2(j, :))**2))
            end do
        end do
    end subroutine pca_score_distance

    pure subroutine derivative_distance(t, derivatives1, derivatives2, distances, p)
        real(dp), intent(in) :: t(:) !! Increasing grid on which derivative curves are supplied.
        real(dp), intent(in) :: derivatives1(:, :) !! First derivative-curve sample, one curve per row.
        real(dp), intent(in) :: derivatives2(:, :) !! Second derivative-curve sample, one curve per row.
        real(dp), intent(out) :: distances(:, :) !! Pairwise Lp distances between derivative curves.
        real(dp), intent(in), optional :: p !! Lp exponent, default two.
        real(dp) :: pp
        integer :: i
        integer :: j

        pp = 2.0_dp
        if (present(p)) pp = p
        do i = 1, size(derivatives1, 1)
            do j = 1, size(derivatives2, 1)
                distances(i, j) = integrate_curve(t, abs(derivatives1(i, :) - derivatives2(j, :))**pp, 2)**(1.0_dp/pp)
            end do
        end do
    end subroutine derivative_distance

    pure subroutine sort_index_by_value(values, index)
        real(dp), intent(in) :: values(:) !! Values used to order corresponding integer indices.
        integer, intent(inout) :: index(:) !! Indices permuted in ascending order of values(index).
        integer :: i
        integer :: j
        integer :: key

        do i = 2, size(index)
            key = index(i)
            j = i - 1
            do while (j >= 1)
                if (values(index(j)) <= values(key)) exit
                index(j + 1) = index(j)
                j = j - 1
            end do
            index(j + 1) = key
        end do
    end subroutine sort_index_by_value

end module fdausc_models
