module fdapace_fpca
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use fdapace_kinds, only : dp
    use fdapace_math, only : jacobi_eigen_sym, linear_interp, sort_eigen_descending, thin_svd
    use fdapace_smoothing, only : trapz_rcpp
    use fdapace_types, only : fpca_result, fsvd_result, select_k_result
    implicit none
    private

    public :: fpca_dense
    public :: fsvd_dense
    public :: get_normalised_sample
    public :: get_normalized_sample
    public :: select_k_fixed
    public :: select_k_fve

contains

    subroutine fpca_dense(y, t, result, fve_threshold, max_k, assume_error, shrink)
        real(dp), intent(in) :: y(:,:) !! Dense regular functional data with subjects in rows and no missing values in this translated path.
        real(dp), intent(in) :: t(:) !! Common increasing regular observation grid corresponding to columns of y.
        type(fpca_result), intent(out) :: result !! Dense cross-sectional FPCA estimates, fitted covariance/correlation, scores, and fitted curves.
        real(dp), intent(in), optional :: fve_threshold !! Fraction-of-variance-explained target in (0,1]; defaults to 0.99.
        integer, intent(in), optional :: max_k !! Maximum number of nonnegative eigencomponents retained before FVE selection.
        logical, intent(in), optional :: assume_error !! If true, estimate measurement-error variance by second differences and remove it from the covariance diagonal.
        logical, intent(in), optional :: shrink !! If true with error estimation, apply the upstream dense IN-score shrinkage factor.
        real(dp), allocatable :: centered(:,:)
        real(dp), allocatable :: cov(:,:)
        real(dp), allocatable :: eval(:)
        real(dp), allocatable :: evec(:,:)
        real(dp), allocatable :: phi_all(:,:)
        real(dp), allocatable :: lambda_all(:)
        real(dp), allocatable :: scaled_phi(:,:)
        real(dp), allocatable :: second(:,:)
        real(dp) :: area
        real(dp) :: fve_target
        real(dp) :: grid_size
        real(dp) :: sign_sum
        real(dp) :: total_eval
        real(dp) :: range_t
        integer :: i
        integer :: info
        integer :: j
        integer :: kk
        integer :: max_components
        integer :: n
        integer :: n_positive
        integer :: p
        logical :: error_model
        logical :: do_shrink

        n = size(y, 1)
        p = size(y, 2)
        if (p /= size(t)) error stop "fpca_dense: y columns do not match t"
        if (n < 2 .or. p < 3) error stop "fpca_dense: dense FPCA requires at least two curves and three grid points"
        if (any(ieee_is_nan(y))) error stop "fpca_dense: translated dense path does not accept missing y values"
        do j = 2, p
            if (t(j) <= t(j - 1)) error stop "fpca_dense: t must be strictly increasing"
        end do
        grid_size = t(2) - t(1)
        if (maxval(abs((t(2:p) - t(1:p - 1)) - grid_size)) > 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(grid_size))) then
            error stop "fpca_dense: translated dense path requires a regular grid"
        end if

        fve_target = 0.99_dp
        if (present(fve_threshold)) fve_target = fve_threshold
        if (fve_target <= 0.0_dp .or. fve_target > 1.0_dp) error stop "fpca_dense: fve_threshold must lie in (0,1]"
        max_components = max(1, min(20, min(max(1, n - 2), max(1, p - 2))))
        if (present(max_k)) max_components = max_k
        if (max_components < 1) error stop "fpca_dense: max_k must be positive"
        error_model = .true.
        if (present(assume_error)) error_model = assume_error
        do_shrink = .false.
        if (present(shrink)) do_shrink = shrink

        allocate(result%mu(p), centered(n, p), cov(p, p))
        result%mu = sum(y, dim=1) / real(n, dp)
        do i = 1, n
            centered(i, :) = y(i, :) - result%mu
        end do
        cov = matmul(transpose(centered), centered) / real(n - 1, dp)
        cov = 0.5_dp * (cov + transpose(cov))

        if (error_model) then
            allocate(second(n, p - 2))
            do i = 1, n
                second(i, :) = y(i, 3:p) - 2.0_dp * y(i, 2:p - 1) + y(i, 1:p - 2)
            end do
            result%sigma2 = sum(second**2) / real(size(second), dp) / 6.0_dp
            result%has_sigma2 = .true.
            do j = 1, p
                cov(j, j) = cov(j, j) - result%sigma2
            end do
        else
            result%sigma2 = 0.0_dp
            result%has_sigma2 = .false.
        end if
        result%smoothed_cov = cov
        result%obs_grid = t
        result%work_grid = t

        allocate(eval(p), evec(p, p))
        call jacobi_eigen_sym(cov, eval, evec, info)
        if (info /= 0) error stop "fpca_dense: symmetric eigensolver did not converge"
        call sort_eigen_descending(eval, evec)
        n_positive = count(eval >= 0.0_dp)
        if (n_positive == 0) error stop "fpca_dense: covariance estimate has no nonnegative eigenvalue"
        n_positive = min(n_positive, max_components)
        allocate(lambda_all(n_positive), phi_all(p, n_positive), result%cum_fve(n_positive))
        total_eval = sum(eval(1:n_positive))
        if (total_eval <= 0.0_dp) error stop "fpca_dense: retained covariance variation is zero"
        result%cum_fve = [(sum(eval(1:j)) / total_eval, j=1, n_positive)]
        kk = n_positive
        do j = 1, n_positive
            if (result%cum_fve(j) >= fve_target) then
                kk = j
                exit
            end if
        end do
        do j = 1, n_positive
            area = trapz_rcpp(t, evec(:, j)**2)
            if (area <= 0.0_dp) error stop "fpca_dense: eigenfunction has zero numerical norm"
            phi_all(:, j) = evec(:, j) / sqrt(area)
            sign_sum = sum(phi_all(:, j) * result%mu)
            if (sign_sum < 0.0_dp) phi_all(:, j) = -phi_all(:, j)
            lambda_all(j) = grid_size * eval(j)
        end do

        allocate(result%lambda(kk), result%phi(p, kk))
        result%lambda = lambda_all(1:kk)
        result%phi = phi_all(:, 1:kk)
        result%select_k = kk
        result%fve = result%cum_fve(kk)

        allocate(scaled_phi(p, n_positive))
        scaled_phi = phi_all
        do j = 1, n_positive
            scaled_phi(:, j) = scaled_phi(:, j) * lambda_all(j)
        end do
        result%fitted_cov = matmul(scaled_phi, transpose(phi_all))
        result%fitted_cov = 0.5_dp * (result%fitted_cov + transpose(result%fitted_cov))
        allocate(result%fitted_corr(p, p))
        do j = 1, p
            do i = 1, p
                if (result%fitted_cov(i, i) > 0.0_dp .and. result%fitted_cov(j, j) > 0.0_dp) then
                    result%fitted_corr(i, j) = result%fitted_cov(i, j) &
                        / sqrt(result%fitted_cov(i, i) * result%fitted_cov(j, j))
                else
                    result%fitted_corr(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
                end if
            end do
        end do
        do j = 1, p
            if (.not. ieee_is_nan(result%fitted_corr(j, j))) result%fitted_corr(j, j) = 1.0_dp
        end do

        allocate(result%xi_est(n, kk), result%fitted_y(n, p))
        range_t = t(p) - t(1)
        do i = 1, n
            do j = 1, kk
                result%xi_est(i, j) = trapz_rcpp(t, centered(i, :) * result%phi(:, j))
                if (do_shrink .and. result%has_sigma2) then
                    result%xi_est(i, j) = result%xi_est(i, j) * result%lambda(j) &
                        / (result%lambda(j) + range_t * result%sigma2 / real(p, dp))
                end if
            end do
            result%fitted_y(i, :) = result%mu + matmul(result%phi, result%xi_est(i, :))
        end do
    end subroutine fpca_dense

    function select_k_fve(fpca_obj, fve_threshold) result(selected)
        type(fpca_result), intent(in) :: fpca_obj !! FPCA result containing cumulative FVE values.
        real(dp), intent(in), optional :: fve_threshold !! Fraction-of-variance-explained target in (0,1]; defaults to 0.95.
        type(select_k_result) :: selected
        real(dp) :: threshold
        real(dp) :: buff
        integer :: i

        threshold = 0.95_dp
        if (present(fve_threshold)) threshold = fve_threshold
        if (threshold <= 0.0_dp .or. threshold > 1.0_dp) error stop "select_k_fve: threshold must lie in (0,1]"
        if (.not. allocated(fpca_obj%cum_fve)) error stop "select_k_fve: cum_fve is unavailable"
        buff = epsilon(1.0_dp)
        selected%k = size(fpca_obj%cum_fve)
        do i = 1, size(fpca_obj%cum_fve)
            if (fpca_obj%cum_fve(i) > threshold - buff) then
                selected%k = i
                exit
            end if
        end do
        selected%criterion = fpca_obj%cum_fve(selected%k)
        selected%has_criterion = .true.
    end function select_k_fve

    function select_k_fixed(fpca_obj, k) result(selected)
        type(fpca_result), intent(in) :: fpca_obj !! FPCA result whose available component count bounds the fixed choice.
        integer, intent(in) :: k !! Positive fixed number of components requested by the caller.
        type(select_k_result) :: selected

        if (.not. allocated(fpca_obj%lambda)) error stop "select_k_fixed: lambda is unavailable"
        if (k < 1 .or. k > size(fpca_obj%lambda)) error stop "select_k_fixed: requested k is unavailable"
        selected%k = k
        selected%criterion = 0.0_dp
        selected%has_criterion = .false.
    end function select_k_fixed

    function get_normalised_sample(fpca_obj, y, t, error_sigma) result(ynorm)
        type(fpca_result), intent(in) :: fpca_obj !! FPCA result supplying the fitted mean, covariance diagonal, work grid, and optional sigma2.
        real(dp), intent(in) :: y(:,:) !! Functional observations with subjects in rows and coordinates given by t.
        real(dp), intent(in) :: t(:) !! Observation grid corresponding to y columns and lying within the FPCA work grid.
        logical, intent(in), optional :: error_sigma !! If true, include fpca_obj%sigma2 in the normalizing variance.
        real(dp), allocatable :: ynorm(:,:)
        real(dp), allocatable :: scale_grid(:)
        real(dp) :: mu_t
        real(dp) :: scale_t
        real(dp) :: sigmae
        integer :: i
        integer :: j

        if (size(y, 2) /= size(t)) error stop "get_normalised_sample: y columns do not match t"
        if (.not. allocated(fpca_obj%fitted_cov) .or. .not. allocated(fpca_obj%mu) &
            .or. .not. allocated(fpca_obj%work_grid)) error stop "get_normalised_sample: incomplete FPCA object"
        if (any([(fpca_obj%fitted_cov(j, j) <= 0.0_dp, j=1, size(fpca_obj%fitted_cov, 1))])) then
            error stop "get_normalised_sample: fitted covariance has nonpositive diagonal values"
        end if
        sigmae = 0.0_dp
        if (present(error_sigma)) then
            if (error_sigma) then
                if (.not. fpca_obj%has_sigma2) error stop "get_normalised_sample: sigma2 was not estimated"
                sigmae = fpca_obj%sigma2
            end if
        end if
        allocate(scale_grid(size(fpca_obj%work_grid)), ynorm(size(y, 1), size(y, 2)))
        do j = 1, size(scale_grid)
            scale_grid(j) = sqrt(sigmae + fpca_obj%fitted_cov(j, j))
        end do
        do j = 1, size(t)
            mu_t = linear_interp(fpca_obj%work_grid, fpca_obj%mu, t(j))
            scale_t = linear_interp(fpca_obj%work_grid, scale_grid, t(j))
            do i = 1, size(y, 1)
                ynorm(i, j) = (y(i, j) - mu_t) / scale_t
            end do
        end do
    end function get_normalised_sample

    function get_normalized_sample(fpca_obj, y, t, error_sigma) result(ynorm)
        type(fpca_result), intent(in) :: fpca_obj !! FPCA result passed unchanged to get_normalised_sample.
        real(dp), intent(in) :: y(:,:) !! Functional observations with subjects in rows.
        real(dp), intent(in) :: t(:) !! Observation grid corresponding to y columns.
        logical, intent(in), optional :: error_sigma !! If true, include estimated measurement-error variance in the normalizer.
        real(dp), allocatable :: ynorm(:,:)

        if (present(error_sigma)) then
            ynorm = get_normalised_sample(fpca_obj, y, t, error_sigma)
        else
            ynorm = get_normalised_sample(fpca_obj, y, t)
        end if
    end function get_normalized_sample

    subroutine fsvd_dense(y1, t1, y2, t2, result, fve_threshold, max_k, fixed_k, flip, assume_error, shrink)
        real(dp), intent(in) :: y1(:,:) !! First dense regular functional sample with subjects in rows.
        real(dp), intent(in) :: t1(:) !! Common regular grid for y1.
        real(dp), intent(in) :: y2(:,:) !! Second paired dense regular functional sample with the same number of subjects as y1.
        real(dp), intent(in) :: t2(:) !! Common regular grid for y2.
        type(fsvd_result), intent(out) :: result !! Cross-covariance SVD, singular functions/values, canonical correlations, and dense integration scores.
        real(dp), intent(in), optional :: fve_threshold !! Squared-singular-value FVE target in (0,1]; defaults to 0.99.
        integer, intent(in), optional :: max_k !! Maximum number of singular components considered; defaults to twenty.
        integer, intent(in), optional :: fixed_k !! Optional fixed number of singular components overriding FVE selection.
        logical, intent(in), optional :: flip !! If true, multiply both singular-function families by -1 as in the R option.
        logical, intent(in), optional :: assume_error !! Measurement-error option passed to both dense FPCA fits; defaults to true.
        logical, intent(in), optional :: shrink !! Dense FPCA score-shrinkage option; defaults to false.
        type(fpca_result) :: fp1
        type(fpca_result) :: fp2
        real(dp), allocatable :: c1(:,:)
        real(dp), allocatable :: c2(:,:)
        real(dp), allocatable :: u(:,:)
        real(dp), allocatable :: sval(:)
        real(dp), allocatable :: v(:,:)
        real(dp) :: area
        real(dp) :: denom1
        real(dp) :: denom2
        real(dp) :: dt1
        real(dp) :: dt2
        real(dp) :: fve_target
        real(dp) :: total_s2
        integer :: i
        integer :: info
        integer :: j
        integer :: kkeep
        integer :: kmax
        integer :: n
        integer :: npos
        logical :: do_flip
        logical :: err
        logical :: shr

        n = size(y1, 1)
        if (size(y2, 1) /= n) error stop "fsvd_dense: paired samples have different subject counts"
        if (size(y1, 2) /= size(t1) .or. size(y2, 2) /= size(t2)) error stop "fsvd_dense: sample columns do not match grids"
        if (n <= 2) error stop "fsvd_dense: dense sample cross covariance requires more than two subjects"
        fve_target = 0.99_dp
        if (present(fve_threshold)) fve_target = fve_threshold
        if (fve_target <= 0.0_dp .or. fve_target > 1.0_dp) error stop "fsvd_dense: fve_threshold must lie in (0,1]"
        kmax = 20
        if (present(max_k)) kmax = max_k
        if (kmax < 1) error stop "fsvd_dense: max_k must be positive"
        do_flip = .false.
        if (present(flip)) do_flip = flip
        err = .true.
        if (present(assume_error)) err = assume_error
        shr = .false.
        if (present(shrink)) shr = shrink

        call fpca_dense(y1, t1, fp1, fve_threshold=fve_target, max_k=min(kmax, size(t1)), assume_error=err, shrink=shr)
        call fpca_dense(y2, t2, fp2, fve_threshold=fve_target, max_k=min(kmax, size(t2)), assume_error=err, shrink=shr)
        allocate(c1(n, size(t1)), c2(n, size(t2)))
        do i = 1, n
            c1(i, :) = y1(i, :) - fp1%mu
            c2(i, :) = y2(i, :) - fp2%mu
        end do
        result%cr_cov = matmul(transpose(c1), c2) / real(n - 2, dp)
        call thin_svd(result%cr_cov, u, sval, v, info)
        if (info /= 0) error stop "fsvd_dense: singular-value decomposition failed"
        npos = count(sval > 0.0_dp)
        npos = min(npos, kmax)
        if (npos < 1) error stop "fsvd_dense: cross covariance has no positive singular values"
        total_s2 = sum(sval(1:npos)**2)
        if (present(fixed_k)) then
            kkeep = min(max(1, fixed_k), npos)
        else
            kkeep = npos
            do j = 1, npos
                if (sum(sval(1:j)**2) / total_s2 >= fve_target) then
                    kkeep = j
                    exit
                end if
            end do
        end if
        result%nsvd = kkeep
        result%fve = sum(sval(1:kkeep)**2) / total_s2
        result%grid1 = t1
        result%grid2 = t2
        allocate(result%s_fun1(size(t1), kkeep), result%s_fun2(size(t2), kkeep), result%s_values(kkeep))
        do j = 1, kkeep
            area = trapz_rcpp(t1, u(:, j)**2)
            if (area <= 0.0_dp) error stop "fsvd_dense: first singular function has zero numerical norm"
            result%s_fun1(:, j) = u(:, j) / sqrt(area)
            area = trapz_rcpp(t2, v(:, j)**2)
            if (area <= 0.0_dp) error stop "fsvd_dense: second singular function has zero numerical norm"
            result%s_fun2(:, j) = v(:, j) / sqrt(area)
            if (do_flip) then
                result%s_fun1(:, j) = -result%s_fun1(:, j)
                result%s_fun2(:, j) = -result%s_fun2(:, j)
            end if
        end do
        dt1 = t1(2) - t1(1)
        dt2 = t2(2) - t2(1)
        result%s_values = sqrt(dt1 * dt2) * sval(1:kkeep)
        allocate(result%can_corr(kkeep))
        do j = 1, kkeep
            denom1 = dt1**2 * dot_product(result%s_fun1(:, j), matmul(fp1%fitted_cov, result%s_fun1(:, j)))
            denom2 = dt2**2 * dot_product(result%s_fun2(:, j), matmul(fp2%fitted_cov, result%s_fun2(:, j)))
            if (denom1 <= 0.0_dp .or. denom2 <= 0.0_dp) then
                result%can_corr(j) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
                result%can_corr(j) = result%s_values(j) / sqrt(denom1 * denom2)
                result%can_corr(j) = sign(min(1.0_dp, abs(result%can_corr(j))), result%can_corr(j))
            end if
        end do

        allocate(result%scores1(n, kkeep), result%scores2(n, kkeep))
        do i = 1, n
            do j = 1, kkeep
                result%scores1(i, j) = trapz_rcpp(t1, c1(i, :) * result%s_fun1(:, j))
                result%scores2(i, j) = trapz_rcpp(t2, c2(i, :) * result%s_fun2(:, j))
            end do
        end do
    end subroutine fsvd_dense

end module fdapace_fpca
