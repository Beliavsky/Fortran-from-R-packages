module fdapace_covariance
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use fdapace_kinds, only : dp
    use fdapace_numerics, only : empirical_quantile, sample_indices_with_replacement
    use fdapace_smoothing, only : lwls1d
    use fdapace_types, only : cov_surface_result, mean_ci_result, mean_curve_result
    implicit none
    private

    public :: get_cov_surface
    public :: get_cr_cor_yx
    public :: get_cr_cor_yz
    public :: get_cr_cov_yx
    public :: get_cr_cov_yz
    public :: get_mean_ci
    public :: get_mean_curve

    interface get_cr_cor_yx
        module procedure get_cr_cor_yx_diag
        module procedure get_cr_cor_yx_matrix
    end interface get_cr_cor_yx

    interface get_cr_cor_yz
        module procedure get_cr_cor_yz_diag
        module procedure get_cr_cor_yz_matrix
    end interface get_cr_cor_yz

contains

    subroutine get_mean_curve(y, t, result, bw, kernel_type)
        real(dp), intent(in) :: y(:,:) !! Dense common-grid functional observations with subjects in rows.
        real(dp), intent(in) :: t(:) !! Strictly increasing common observation grid corresponding to columns of y.
        type(mean_curve_result), intent(out) :: result !! Estimated mean curve, work grid, and smoothing bandwidth metadata.
        real(dp), intent(in), optional :: bw !! Optional positive local-linear bandwidth; omission uses the upstream dense cross-sectional mean.
        character(len=*), intent(in), optional :: kernel_type !! Optional smoothing kernel name passed to lwls1d; defaults to gauss.
        character(len=16) :: kernel
        real(dp), allocatable :: raw_mu(:)
        real(dp), allocatable :: weights(:)

        if (size(y, 2) /= size(t) .or. size(y, 1) < 1) error stop "get_mean_curve: invalid dense dimensions"
        if (size(t) < 1) error stop "get_mean_curve: empty grid"
        allocate(raw_mu(size(t)))
        raw_mu = sum(y, dim=1) / real(size(y, 1), dp)
        result%work_grid = t
        result%bw_mu = 0.0_dp
        if (present(bw)) then
            if (bw <= 0.0_dp) error stop "get_mean_curve: bw must be positive"
            kernel = "gauss"
            if (present(kernel_type)) kernel = trim(kernel_type)
            allocate(weights(size(t)))
            weights = 1.0_dp
            result%mu = lwls1d(bw, kernel, t, raw_mu, t, weights, 1, 0)
            result%bw_mu = bw
        else
            result%mu = raw_mu
        end if
    end subroutine get_mean_curve

    subroutine get_cov_surface(y, t, result, assume_error)
        real(dp), intent(in) :: y(:,:) !! Dense common-grid functional observations with subjects in rows.
        real(dp), intent(in) :: t(:) !! Strictly increasing common observation grid corresponding to columns of y.
        type(cov_surface_result), intent(out) :: result !! Cross-sectional covariance estimate, grid, and optional measurement-error variance.
        logical, intent(in), optional :: assume_error !! If true, estimate measurement-error variance from second differences and subtract it from the diagonal.
        real(dp), allocatable :: centered(:,:)
        real(dp), allocatable :: second(:,:)
        real(dp), allocatable :: mu(:)
        logical :: error_model
        integer :: i
        integer :: j
        integer :: n
        integer :: p

        n = size(y, 1)
        p = size(y, 2)
        if (p /= size(t) .or. n < 2 .or. p < 1) error stop "get_cov_surface: invalid dense dimensions"
        allocate(mu(p), centered(n, p), result%cov(p, p))
        mu = sum(y, dim=1) / real(n, dp)
        do i = 1, n
            centered(i, :) = y(i, :) - mu
        end do
        result%cov = matmul(transpose(centered), centered) / real(n - 1, dp)
        result%cov = 0.5_dp * (result%cov + transpose(result%cov))
        error_model = .true.
        if (present(assume_error)) error_model = assume_error
        result%sigma2 = 0.0_dp
        if (error_model .and. p >= 3) then
            allocate(second(n, p - 2))
            do i = 1, n
                second(i, :) = y(i, 3:p) - 2.0_dp * y(i, 2:p - 1) + y(i, 1:p - 2)
            end do
            result%sigma2 = sum(second**2) / real(size(second), dp) / 6.0_dp
            do j = 1, p
                result%cov(j, j) = result%cov(j, j) - result%sigma2
            end do
        end if
        result%work_grid = t
        result%bw_cov = 0.0_dp
    end subroutine get_cov_surface

    subroutine get_mean_ci(y, t, result, level, n_boot)
        real(dp), intent(in) :: y(:,:) !! Dense common-grid sample whose mean curve is bootstrapped by subject.
        real(dp), intent(in) :: t(:) !! Common observation grid corresponding to columns of y.
        type(mean_ci_result), intent(out) :: result !! Pointwise percentile-bootstrap confidence limits and their grid.
        real(dp), intent(in), optional :: level !! Confidence level in [0,1]; defaults to 0.95.
        integer, intent(in), optional :: n_boot !! Number of bootstrap replicates; defaults to 999 and must be positive.
        integer, allocatable :: idx(:)
        real(dp), allocatable :: boot(:,:)
        real(dp) :: lev
        integer :: b
        integer :: j
        integer :: nb
        integer :: n

        n = size(y, 1)
        if (n < 1 .or. size(y, 2) /= size(t)) error stop "get_mean_ci: invalid dense dimensions"
        lev = 0.95_dp
        if (present(level)) lev = level
        if (lev < 0.0_dp .or. lev > 1.0_dp) error stop "get_mean_ci: level must lie in [0,1]"
        nb = 999
        if (present(n_boot)) nb = n_boot
        if (nb < 1) error stop "get_mean_ci: n_boot must be positive"
        allocate(idx(n), boot(nb, size(t)), result%lower(size(t)), result%upper(size(t)))
        do b = 1, nb
            call sample_indices_with_replacement(n, idx)
            boot(b, :) = sum(y(idx, :), dim=1) / real(n, dp)
        end do
        do j = 1, size(t)
            result%lower(j) = empirical_quantile(boot(:, j), 0.5_dp * (1.0_dp - lev))
            result%upper(j) = empirical_quantile(boot(:, j), 1.0_dp - 0.5_dp * (1.0_dp - lev))
        end do
        result%grid = t
        result%level = lev
    end subroutine get_mean_ci

    function get_cr_cov_yx(y1, y2) result(cross_cov)
        real(dp), intent(in) :: y1(:,:) !! First dense functional sample with paired subjects in rows.
        real(dp), intent(in) :: y2(:,:) !! Second dense functional sample with the same paired-subject count as y1.
        real(dp), allocatable :: cross_cov(:,:)
        real(dp), allocatable :: c1(:,:)
        real(dp), allocatable :: c2(:,:)
        real(dp), allocatable :: m1(:)
        real(dp), allocatable :: m2(:)
        integer :: i
        integer :: n

        n = size(y1, 1)
        if (n /= size(y2, 1) .or. n < 2) error stop "get_cr_cov_yx: paired samples require at least two subjects"
        allocate(m1(size(y1, 2)), m2(size(y2, 2)), c1(n, size(y1, 2)), c2(n, size(y2, 2)))
        m1 = sum(y1, dim=1) / real(n, dp)
        m2 = sum(y2, dim=1) / real(n, dp)
        do i = 1, n
            c1(i, :) = y1(i, :) - m1
            c2(i, :) = y2(i, :) - m2
        end do
        cross_cov = matmul(transpose(c1), c2) / real(n - 1, dp)
    end function get_cr_cov_yx

    function get_cr_cov_yz(y, z) result(cross_cov)
        real(dp), intent(in) :: y(:,:) !! Dense functional sample with subjects in rows.
        real(dp), intent(in) :: z(:) !! Paired scalar response, one value per subject.
        real(dp), allocatable :: cross_cov(:)
        real(dp), allocatable :: centered_z(:)
        real(dp), allocatable :: mu_y(:)
        real(dp) :: mu_z
        integer :: i
        integer :: n

        n = size(y, 1)
        if (size(z) /= n .or. n < 2) error stop "get_cr_cov_yz: z length must match at least two subjects"
        allocate(mu_y(size(y, 2)), centered_z(n), cross_cov(size(y, 2)))
        mu_y = sum(y, dim=1) / real(n, dp)
        mu_z = sum(z) / real(n, dp)
        centered_z = z - mu_z
        cross_cov = 0.0_dp
        do i = 1, n
            cross_cov = cross_cov + (y(i, :) - mu_y) * centered_z(i)
        end do
        cross_cov = cross_cov / real(n - 1, dp)
    end function get_cr_cov_yz

    function get_cr_cor_yx_diag(cc_xy, diag_xx, diag_yy) result(cross_cor)
        real(dp), intent(in) :: cc_xy(:,:) !! Cross-covariance matrix from the first support to the second support.
        real(dp), intent(in) :: diag_xx(:) !! Positive marginal variances for rows of cc_xy.
        real(dp), intent(in) :: diag_yy(:) !! Positive marginal variances for columns of cc_xy.
        real(dp), allocatable :: cross_cor(:,:)
        integer :: i
        integer :: j

        if (size(diag_xx) /= size(cc_xy, 1) .or. size(diag_yy) /= size(cc_xy, 2)) then
            error stop "get_cr_cor_yx: covariance dimensions are incompatible"
        end if
        if (any(diag_xx <= 1.0e-12_dp) .or. any(diag_yy <= 1.0e-12_dp)) then
            error stop "get_cr_cor_yx: marginal variances must exceed 1e-12"
        end if
        allocate(cross_cor(size(cc_xy, 1), size(cc_xy, 2)))
        do j = 1, size(cc_xy, 2)
            do i = 1, size(cc_xy, 1)
                cross_cor(i, j) = cc_xy(i, j) / sqrt(diag_xx(i) * diag_yy(j))
            end do
        end do
    end function get_cr_cor_yx_diag

    function get_cr_cor_yx_matrix(cc_xy, cc_xx, cc_yy) result(cross_cor)
        real(dp), intent(in) :: cc_xy(:,:) !! Cross-covariance matrix from the first support to the second support.
        real(dp), intent(in) :: cc_xx(:,:) !! Square auto-covariance matrix for rows of cc_xy.
        real(dp), intent(in) :: cc_yy(:,:) !! Square auto-covariance matrix for columns of cc_xy.
        real(dp), allocatable :: cross_cor(:,:)
        real(dp), allocatable :: diag_x(:)
        real(dp), allocatable :: diag_y(:)
        integer :: i

        if (size(cc_xx, 1) /= size(cc_xx, 2) .or. size(cc_yy, 1) /= size(cc_yy, 2)) then
            error stop "get_cr_cor_yx: auto-covariances must be square"
        end if
        allocate(diag_x(size(cc_xx, 1)), diag_y(size(cc_yy, 1)))
        do i = 1, size(diag_x)
            diag_x(i) = cc_xx(i, i)
        end do
        do i = 1, size(diag_y)
            diag_y(i) = cc_yy(i, i)
        end do
        cross_cor = get_cr_cor_yx_diag(cc_xy, diag_x, diag_y)
    end function get_cr_cor_yx_matrix

    function get_cr_cor_yz_diag(cc_yz, diag_yy, var_z) result(cross_cor)
        real(dp), intent(in) :: cc_yz(:) !! Functional-scalar cross-covariance vector over the functional support.
        real(dp), intent(in) :: diag_yy(:) !! Positive marginal variances of the functional process.
        real(dp), intent(in) :: var_z !! Positive scalar variance of the scalar response.
        real(dp), allocatable :: cross_cor(:)

        if (size(cc_yz) /= size(diag_yy)) error stop "get_cr_cor_yz: covariance dimensions are incompatible"
        if (any(diag_yy <= 1.0e-12_dp) .or. var_z <= 1.0e-12_dp) then
            error stop "get_cr_cor_yz: marginal variances must exceed 1e-12"
        end if
        cross_cor = cc_yz / sqrt(diag_yy * var_z)
    end function get_cr_cor_yz_diag

    function get_cr_cor_yz_matrix(cc_yz, cc_yy, var_z) result(cross_cor)
        real(dp), intent(in) :: cc_yz(:) !! Functional-scalar cross-covariance vector over the functional support.
        real(dp), intent(in) :: cc_yy(:,:) !! Square functional auto-covariance matrix.
        real(dp), intent(in) :: var_z !! Positive scalar variance of the scalar response.
        real(dp), allocatable :: cross_cor(:)
        real(dp), allocatable :: diag_y(:)
        integer :: i

        if (size(cc_yy, 1) /= size(cc_yy, 2)) error stop "get_cr_cor_yz: auto-covariance must be square"
        allocate(diag_y(size(cc_yy, 1)))
        do i = 1, size(diag_y)
            diag_y(i) = cc_yy(i, i)
        end do
        cross_cor = get_cr_cor_yz_diag(cc_yz, diag_y, var_z)
    end function get_cr_cor_yz_matrix

end module fdapace_covariance
