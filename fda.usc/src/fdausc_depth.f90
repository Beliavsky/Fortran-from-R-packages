module fdausc_depth
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use fdausc_integration, only : integrate_curve
    use fdausc_kernels, only : ker_norm
    implicit none
    private

    public :: depth_mode_scores, functional_spatial_depth, kernel_functional_spatial_depth
    public :: fraiman_muniz_depth, tukey_coordinate_depth, mahalanobis_depth
    public :: multivariate_spatial_depth, multivariate_kernel_spatial_depth
    public :: projection_tukey_depth, depth_summary

contains

    pure subroutine depth_mode_scores(distances_eval_ref, bandwidth, depth, scale)
        real(dp), intent(in) :: distances_eval_ref(:, :) !! Evaluation-to-reference distances, evaluated observations in rows.
        real(dp), intent(in) :: bandwidth !! Positive kernel bandwidth used in the modal depth.
        real(dp), intent(out) :: depth(:) !! Modal depth scores, proportional to sums of normal-kernel similarities.
        logical, intent(in), optional :: scale !! If true, divide all scores by their maximum so the deepest point has depth one.
        integer :: i
        logical :: do_scale

        do i = 1, size(distances_eval_ref, 1)
            depth(i) = sum(ker_norm(distances_eval_ref(i, :)/bandwidth))
        end do
        do_scale = .false.
        if (present(scale)) do_scale = scale
        if (do_scale .and. maxval(depth) > 0.0_dp) depth = depth/maxval(depth)
    end subroutine depth_mode_scores

    pure subroutine functional_spatial_depth(xgrid, curves, reference, depth, method)
        real(dp), intent(in) :: xgrid(:) !! Increasing functional grid common to evaluated and reference curves.
        real(dp), intent(in) :: curves(:, :) !! Curves whose depth is evaluated, one curve per row.
        real(dp), intent(in) :: reference(:, :) !! Reference functional sample, one curve per row.
        real(dp), intent(out) :: depth(:) !! Functional spatial-depth values in [0,1] up to numerical error.
        integer, intent(in), optional :: method !! Integration selector passed to L2 norms, default composite Simpson.
        real(dp), allocatable :: mean_sign(:)
        real(dp), allocatable :: diff(:)
        real(dp), allocatable :: sq(:)
        real(dp) :: normdiff
        real(dp) :: normmean
        integer :: i
        integer :: j
        integer :: meth

        meth = 2
        if (present(method)) meth = method
        allocate(mean_sign(size(xgrid)), diff(size(xgrid)), sq(size(xgrid)))
        do i = 1, size(curves, 1)
            mean_sign = 0.0_dp
            do j = 1, size(reference, 1)
                diff = reference(j, :) - curves(i, :)
                sq = diff*diff
                normdiff = sqrt(max(0.0_dp, integrate_curve(xgrid, sq, meth)))
                if (normdiff > sqrt(tiny(1.0_dp))) mean_sign = mean_sign + diff/normdiff
            end do
            mean_sign = mean_sign/real(size(reference, 1), dp)
            sq = mean_sign*mean_sign
            normmean = sqrt(max(0.0_dp, integrate_curve(xgrid, sq, meth)))
            depth(i) = max(0.0_dp, 1.0_dp - normmean)
        end do
    end subroutine functional_spatial_depth

    pure subroutine kernel_functional_spatial_depth(xgrid, curves, reference, bandwidth, depth, method)
        real(dp), intent(in) :: xgrid(:) !! Increasing functional grid shared by evaluated and reference curves.
        real(dp), intent(in) :: curves(:, :) !! Curves whose kernelized spatial depth is evaluated.
        real(dp), intent(in) :: reference(:, :) !! Reference functional sample.
        real(dp), intent(in) :: bandwidth !! Positive Gaussian-kernel bandwidth in L2 distance units.
        real(dp), intent(out) :: depth(:) !! Kernelized functional spatial-depth values.
        integer, intent(in), optional :: method !! Integration selector used for L2 distances, default composite Simpson.
        real(dp), allocatable :: krr(:, :)
        real(dp), allocatable :: ker(:)
        real(dp), allocatable :: diff(:)
        real(dp) :: dist2
        real(dp) :: denom1
        real(dp) :: denom2
        real(dp) :: s
        integer :: i
        integer :: j
        integer :: k
        integer :: meth
        integer :: nr

        meth = 2
        if (present(method)) meth = method
        nr = size(reference, 1)
        allocate(krr(nr, nr), ker(nr), diff(size(xgrid)))
        do j = 1, nr
            do k = j, nr
                diff = reference(j, :) - reference(k, :)
                dist2 = integrate_curve(xgrid, diff*diff, meth)
                krr(j, k) = exp(-dist2/(bandwidth*bandwidth))
                krr(k, j) = krr(j, k)
            end do
        end do
        do i = 1, size(curves, 1)
            do j = 1, nr
                diff = curves(i, :) - reference(j, :)
                dist2 = integrate_curve(xgrid, diff*diff, meth)
                ker(j) = exp(-dist2/(bandwidth*bandwidth))
            end do
            s = 0.0_dp
            do j = 1, nr
                denom1 = sqrt(max(0.0_dp, 2.0_dp - 2.0_dp*ker(j)))
                if (denom1 <= sqrt(tiny(1.0_dp))) cycle
                do k = 1, nr
                    denom2 = sqrt(max(0.0_dp, 2.0_dp - 2.0_dp*ker(k)))
                    if (denom2 <= sqrt(tiny(1.0_dp))) cycle
                    s = s + (1.0_dp + krr(j, k) - ker(j) - ker(k))/(denom1*denom2)
                end do
            end do
            depth(i) = max(0.0_dp, 1.0_dp - sqrt(max(0.0_dp, s))/real(nr, dp))
        end do
    end subroutine kernel_functional_spatial_depth

    pure subroutine fraiman_muniz_depth(curves, reference, depth, scale)
        real(dp), intent(in) :: curves(:, :) !! Functional observations whose pointwise empirical depth is evaluated.
        real(dp), intent(in) :: reference(:, :) !! Reference functional sample defining each pointwise empirical CDF.
        real(dp), intent(out) :: depth(:) !! Mean pointwise Fraiman-Muniz depth for each evaluated curve.
        logical, intent(in), optional :: scale !! If true, divide by the maximum score among evaluated curves.
        real(dp) :: f
        real(dp) :: d
        integer :: i
        integer :: j
        integer :: k
        logical :: do_scale

        depth = 0.0_dp
        do i = 1, size(curves, 1)
            do j = 1, size(curves, 2)
                k = count(reference(:, j) <= curves(i, j))
                f = real(k, dp)/real(size(reference, 1), dp)
                d = 1.0_dp - abs(0.5_dp - f)
                depth(i) = depth(i) + d
            end do
            depth(i) = depth(i)/real(size(curves, 2), dp)
        end do
        do_scale = .false.
        if (present(scale)) do_scale = scale
        if (do_scale .and. maxval(depth) > 0.0_dp) depth = depth/maxval(depth)
    end subroutine fraiman_muniz_depth

    pure subroutine tukey_coordinate_depth(x, reference, depth, scale)
        real(dp), intent(in) :: x(:, :) !! Evaluated multivariate observations, one row per observation.
        real(dp), intent(in) :: reference(:, :) !! Reference multivariate sample defining marginal empirical CDFs.
        real(dp), intent(out) :: depth(:) !! Average coordinatewise Tukey depth of each row of x.
        logical, intent(in), optional :: scale !! If true, multiply coordinate depths by two before averaging.
        real(dp) :: f_le
        real(dp) :: f_lt
        real(dp) :: d
        integer :: i
        integer :: j
        logical :: do_scale

        do_scale = .false.
        if (present(scale)) do_scale = scale
        do i = 1, size(x, 1)
            depth(i) = 0.0_dp
            do j = 1, size(x, 2)
                f_le = real(count(reference(:, j) <= x(i, j)), dp)/real(size(reference, 1), dp)
                f_lt = real(count(reference(:, j) < x(i, j)), dp)/real(size(reference, 1), dp)
                d = min(f_le, 1.0_dp - f_lt)
                if (do_scale) d = 2.0_dp*d
                depth(i) = depth(i) + d
            end do
            depth(i) = depth(i)/real(size(x, 2), dp)
        end do
    end subroutine tukey_coordinate_depth

    subroutine mahalanobis_depth(x, reference, depth, info)
        real(dp), intent(in) :: x(:, :) !! Evaluated multivariate observations, one row per observation.
        real(dp), intent(in) :: reference(:, :) !! Reference sample used to estimate mean and covariance.
        real(dp), intent(out) :: depth(:) !! Mahalanobis depths 1/(1+squared distance).
        integer, intent(out), optional :: info !! Zero on success; positive if the sample covariance cannot be inverted.
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: invcov(:, :)
        real(dp), allocatable :: d(:)
        integer :: i
        integer :: j
        integer :: k
        integer :: stat

        mu = sum(reference, dim=1)/real(size(reference, 1), dp)
        allocate(cov(size(reference, 2), size(reference, 2)), d(size(reference, 2)))
        cov = 0.0_dp
        do i = 1, size(reference, 1)
            d = reference(i, :) - mu
            do j = 1, size(d)
                do k = 1, size(d)
                    cov(j, k) = cov(j, k) + d(j)*d(k)
                end do
            end do
        end do
        cov = cov/real(max(1, size(reference, 1) - 1), dp)
        call inverse_matrix(cov, invcov, stat)
        if (present(info)) info = stat
        if (stat /= 0) then
            depth = 0.0_dp
            return
        end if
        do i = 1, size(x, 1)
            d = x(i, :) - mu
            depth(i) = 1.0_dp/(1.0_dp + max(0.0_dp, dot_product(d, matmul(invcov, d))))
        end do
    end subroutine mahalanobis_depth

    pure subroutine multivariate_spatial_depth(x, reference, depth)
        real(dp), intent(in) :: x(:, :) !! Evaluated multivariate observations, one row per point.
        real(dp), intent(in) :: reference(:, :) !! Reference multivariate sample.
        real(dp), intent(out) :: depth(:) !! Spatial depth 1 minus norm of the mean spatial sign.
        real(dp), allocatable :: sign_mean(:)
        real(dp), allocatable :: d(:)
        real(dp) :: dn
        integer :: i
        integer :: j

        allocate(sign_mean(size(x, 2)), d(size(x, 2)))
        do i = 1, size(x, 1)
            sign_mean = 0.0_dp
            do j = 1, size(reference, 1)
                d = reference(j, :) - x(i, :)
                dn = sqrt(dot_product(d, d))
                if (dn > sqrt(tiny(1.0_dp))) sign_mean = sign_mean + d/dn
            end do
            sign_mean = sign_mean/real(size(reference, 1), dp)
            depth(i) = max(0.0_dp, 1.0_dp - sqrt(dot_product(sign_mean, sign_mean)))
        end do
    end subroutine multivariate_spatial_depth

    pure subroutine multivariate_kernel_spatial_depth(x, reference, bandwidth, depth)
        real(dp), intent(in) :: x(:, :) !! Evaluated multivariate observations.
        real(dp), intent(in) :: reference(:, :) !! Reference multivariate sample.
        real(dp), intent(in) :: bandwidth !! Positive Gaussian-kernel bandwidth in Euclidean-distance units.
        real(dp), intent(out) :: depth(:) !! Kernelized spatial-depth values.
        real(dp), allocatable :: krr(:, :)
        real(dp), allocatable :: ker(:)
        real(dp), allocatable :: d(:)
        real(dp) :: denom1
        real(dp) :: denom2
        real(dp) :: s
        integer :: i
        integer :: j
        integer :: k
        integer :: nr

        nr = size(reference, 1)
        allocate(krr(nr, nr), ker(nr), d(size(x, 2)))
        do j = 1, nr
            do k = j, nr
                d = reference(j, :) - reference(k, :)
                krr(j, k) = exp(-dot_product(d, d)/(bandwidth*bandwidth))
                krr(k, j) = krr(j, k)
            end do
        end do
        do i = 1, size(x, 1)
            do j = 1, nr
                d = x(i, :) - reference(j, :)
                ker(j) = exp(-dot_product(d, d)/(bandwidth*bandwidth))
            end do
            s = 0.0_dp
            do j = 1, nr
                denom1 = sqrt(max(0.0_dp, 2.0_dp - 2.0_dp*ker(j)))
                if (denom1 <= sqrt(tiny(1.0_dp))) cycle
                do k = 1, nr
                    denom2 = sqrt(max(0.0_dp, 2.0_dp - 2.0_dp*ker(k)))
                    if (denom2 <= sqrt(tiny(1.0_dp))) cycle
                    s = s + (1.0_dp + krr(j, k) - ker(j) - ker(k))/(denom1*denom2)
                end do
            end do
            depth(i) = max(0.0_dp, 1.0_dp - sqrt(max(0.0_dp, s))/real(nr, dp))
        end do
    end subroutine multivariate_kernel_spatial_depth

    pure subroutine projection_tukey_depth(x, reference, projections, depth, average)
        real(dp), intent(in) :: x(:, :) !! Evaluated multivariate or discretized functional observations.
        real(dp), intent(in) :: reference(:, :) !! Reference observations with the same number of columns as x.
        real(dp), intent(in) :: projections(:, :) !! Projection directions in rows; each row has size(x,2) entries.
        real(dp), intent(out) :: depth(:) !! Mean or minimum Tukey depth across the supplied projection directions.
        logical, intent(in), optional :: average !! If true average over directions; otherwise take the minimum; default true.
        real(dp), allocatable :: zr(:)
        real(dp) :: z
        real(dp) :: f_le
        real(dp) :: f_lt
        real(dp) :: d
        integer :: i
        integer :: j
        logical :: mean_depth

        mean_depth = .true.
        if (present(average)) mean_depth = average
        allocate(zr(size(reference, 1)))
        do i = 1, size(x, 1)
            if (mean_depth) then
                depth(i) = 0.0_dp
            else
                depth(i) = huge(1.0_dp)
            end if
            do j = 1, size(projections, 1)
                zr = matmul(reference, projections(j, :))
                z = dot_product(x(i, :), projections(j, :))
                f_le = real(count(zr <= z), dp)/real(size(zr), dp)
                f_lt = real(count(zr < z), dp)/real(size(zr), dp)
                d = min(f_le, 1.0_dp - f_lt)
                if (mean_depth) then
                    depth(i) = depth(i) + d/real(size(projections, 1), dp)
                else
                    depth(i) = min(depth(i), d)
                end if
            end do
        end do
    end subroutine projection_tukey_depth

    pure subroutine depth_summary(curves, depth, trim, median_index, trimmed_mean, trimmed_variance, selected)
        real(dp), intent(in) :: curves(:, :) !! Functional or multivariate observations in rows.
        real(dp), intent(in) :: depth(:) !! Depth score for each row of curves.
        real(dp), intent(in) :: trim !! Lower depth quantile defining observations retained in the trimmed summary.
        integer, intent(out) :: median_index !! One-based index of the deepest observation.
        real(dp), intent(out) :: trimmed_mean(:) !! Coordinatewise mean over observations at or above the trim-depth quantile.
        real(dp), intent(out) :: trimmed_variance(:) !! Coordinatewise population variance over the retained observations.
        logical, intent(out) :: selected(:) !! Logical mask indicating which observations are retained.
        real(dp), allocatable :: sorted(:)
        real(dp) :: cutoff
        integer :: i
        integer :: nsel
        integer :: idx

        median_index = maxloc(depth, dim=1)
        sorted = depth
        call sort_real(sorted)
        idx = max(1, min(size(sorted), 1 + int(trim*real(size(sorted) - 1, dp))))
        cutoff = sorted(idx)
        selected = depth >= cutoff
        nsel = count(selected)
        trimmed_mean = 0.0_dp
        do i = 1, size(curves, 1)
            if (selected(i)) trimmed_mean = trimmed_mean + curves(i, :)
        end do
        if (nsel > 0) trimmed_mean = trimmed_mean/real(nsel, dp)
        trimmed_variance = 0.0_dp
        do i = 1, size(curves, 1)
            if (selected(i)) trimmed_variance = trimmed_variance + (curves(i, :) - trimmed_mean)**2
        end do
        if (nsel > 0) trimmed_variance = trimmed_variance/real(nsel, dp)
    end subroutine depth_summary

    pure subroutine sort_real(x)
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
    end subroutine sort_real

end module fdausc_depth
