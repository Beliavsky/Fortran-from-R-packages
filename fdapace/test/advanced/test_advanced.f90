program test_advanced
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use fdapace, only : cluster_result, cov_surface_result, dp, fam, fam_result, fc_reg, fclust, fcreg_result, &
        flm_ci_result, flm_ci_scalar_dense, flm_result, flm_scalar_dense, f_opt_des, fopt_result, fpc_quantile, &
        fpc_quantile_result, fpca_dense, fpca_der, fpca_der_result, fpca_result, fvpa, fvpa_result, get_cov_surface, &
        get_cr_cor_yx, get_cr_cor_yz, get_cr_cov_yx, get_cr_cov_yz, get_mean_ci, get_mean_curve, k_cfc, &
        make_sparse_gp, mean_ci_result, mean_curve_result, multi_fam, multifam_result, sb_fitting, sbf_result, &
        seed_rng, sparse_gp_result, stringing, stringing_result, tvam, tvam_result, vcam, vcam_result, wfda, wfda_result
    implicit none

    integer :: failures

    failures = 0
    call test_covariance_family(failures)
    call test_sparse_and_derivatives(failures)
    call test_ordering_and_clustering(failures)
    call test_regression_family(failures)
    call test_longitudinal_family(failures)

    if (failures /= 0) then
        write (*, '(a,i0)') 'advanced fdapace tests failed: ', failures
        error stop 1
    end if
    write (*, '(a)') 'All advanced fdapace tests passed.'

contains

    subroutine check(condition, message, failures)
        logical, intent(in) :: condition !! Test condition that must be true.
        character(len=*), intent(in) :: message !! Diagnostic text printed when the assertion fails.
        integer, intent(inout) :: failures !! Running number of failed assertions.

        if (.not. condition) then
            failures = failures + 1
            write (*, '(a)') 'FAIL: ' // trim(message)
        end if
    end subroutine check

    subroutine check_close(actual, expected, tol, message, failures)
        real(dp), intent(in) :: actual !! Computed scalar value.
        real(dp), intent(in) :: expected !! Reference scalar value.
        real(dp), intent(in) :: tol !! Maximum permitted absolute error.
        character(len=*), intent(in) :: message !! Diagnostic text printed when the tolerance is exceeded.
        integer, intent(inout) :: failures !! Running number of failed assertions.

        call check(abs(actual - expected) <= tol, message, failures)
    end subroutine check_close

    subroutine make_rank_data(t, scores, y1, y2)
        real(dp), intent(out) :: t(:) !! Equally spaced common grid generated on [0,1].
        real(dp), intent(in) :: scores(:) !! Subject scores controlling the rank-one functional variation.
        real(dp), intent(out) :: y1(:,:) !! First rank-one functional sample with rows matching scores.
        real(dp), intent(out) :: y2(:,:) !! Second paired rank-one functional sample with rows matching scores.
        integer :: i

        t = [(real(i - 1, dp) / real(size(t) - 1, dp), i=1, size(t))]
        do i = 1, size(scores)
            y1(i, :) = 1.0_dp + 0.3_dp * t + scores(i) * (1.0_dp + 0.4_dp * cos(acos(-1.0_dp) * t))
            y2(i, :) = -0.5_dp + 0.2_dp * t + 1.5_dp * scores(i) * (0.2_dp + sin(0.5_dp * acos(-1.0_dp) * t))
        end do
    end subroutine make_rank_data

    subroutine test_covariance_family(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: t(9)
        real(dp) :: scores(10)
        real(dp) :: y1(10, 9)
        real(dp) :: y2(10, 9)
        real(dp) :: z(10)
        real(dp), allocatable :: cc(:,:)
        real(dp), allocatable :: cc_z(:)
        real(dp), allocatable :: cor(:,:)
        real(dp), allocatable :: cor_z(:)
        real(dp), allocatable :: var1(:)
        real(dp), allocatable :: var2(:)
        type(mean_curve_result) :: mu
        type(cov_surface_result) :: cov
        type(mean_ci_result) :: ci
        integer :: i

        scores = [(real(i, dp) - 5.5_dp, i=1, 10)]
        call make_rank_data(t, scores, y1, y2)
        z = 2.0_dp + 0.7_dp * scores
        call get_mean_curve(y1, t, mu)
        call check(all(shape(mu%mu) == [9]), 'get_mean_curve shape', failures)
        call check_close(mu%mu(1), 1.0_dp, 1.0e-12_dp, 'get_mean_curve centered-score mean', failures)
        call get_cov_surface(y1, t, cov, assume_error=.false.)
        call check(all(shape(cov%cov) == [9, 9]), 'get_cov_surface shape', failures)
        cc = get_cr_cov_yx(y1, y2)
        call check(all(shape(cc) == [9, 9]), 'get_cr_cov_yx shape', failures)
        allocate(var1(9), var2(9))
        do i = 1, 9
            var1(i) = cov%cov(i, i)
            var2(i) = sum((y2(:, i) - sum(y2(:, i)) / 10.0_dp)**2) / 9.0_dp
        end do
        cor = get_cr_cor_yx(cc, var1, var2)
        call check(maxval(abs(cor - 1.0_dp)) < 1.0e-10_dp, 'get_cr_cor_yx paired rank-one correlation', failures)
        cc_z = get_cr_cov_yz(y1, z)
        cor_z = get_cr_cor_yz(cc_z, var1, sum((z - sum(z) / 10.0_dp)**2) / 9.0_dp)
        call check(maxval(abs(cor_z - 1.0_dp)) < 1.0e-10_dp, 'get_cr_cor_yz rank-one correlation', failures)
        call seed_rng(123)
        call get_mean_ci(y1, t, ci, level=0.8_dp, n_boot=20)
        call check(all(ci%lower <= ci%upper), 'get_mean_ci ordered limits', failures)
    end subroutine test_covariance_family

    subroutine test_sparse_and_derivatives(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: t(11)
        real(dp) :: scores(12)
        real(dp) :: y1(12, 11)
        real(dp) :: y2(12, 11)
        real(dp) :: yv(12, 11)
        type(sparse_gp_result) :: sp
        type(fpca_result) :: fit
        type(fpca_der_result) :: der
        type(fvpa_result) :: vp
        integer :: i
        integer :: j

        scores = [(real(i, dp) - 6.5_dp, i=1, 12)]
        call make_rank_data(t, scores, y1, y2)
        call seed_rng(456)
        sp = make_sparse_gp(8, [3, 4], 2, [1.0_dp, 0.5_dp], sigma=0.1_dp, basis_type='cos')
        call check(size(sp%sample%ly) == 8 .and. sp%has_true, 'make_sparse_gp output and truth', failures)
        call fpca_dense(y1, t, fit, assume_error=.false.)
        call fpca_der(fit, der, derivative_order=1, bandwidth=0.5_dp)
        call check(all(shape(der%phi_der) == shape(fit%phi)), 'fpca_der eigenfunction shape', failures)
        do i = 1, 12
            do j = 1, 11
                yv(i, j) = y1(i, j) + 0.08_dp * sin(real(i * j, dp)) + 0.03_dp * cos(real(2 * i + j, dp))
            end do
        end do
        call fvpa(yv, t, vp, q=0.2_dp, fve_threshold=0.8_dp)
        call check(ieee_is_finite(vp%sigma2) .and. vp%sigma2 >= 0.0_dp, 'fvpa finite residual variance', failures)
    end subroutine test_sparse_and_derivatives

    subroutine test_ordering_and_clustering(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: x(20, 4)
        real(dp) :: t(9)
        real(dp) :: y(12, 9)
        real(dp) :: scores(12)
        real(dp) :: dummy(12, 9)
        type(stringing_result) :: st
        type(cluster_result) :: cl
        type(cluster_result) :: kc
        type(wfda_result) :: wr
        integer :: i

        do i = 1, 20
            x(i, 1) = real(i, dp)
            x(i, 2) = real(i, dp) + 0.1_dp * sin(real(i, dp))
            x(i, 3) = -real(i, dp)
            x(i, 4) = sin(real(i, dp))
        end do
        call stringing(x, st, standardize=.true., distance='correlation')
        call check(size(st%order) == 4 .and. all(sort_copy(st%order) == [1, 2, 3, 4]), 'stringing permutation', failures)
        scores = [-3.0_dp, -2.8_dp, -2.5_dp, -2.2_dp, -1.9_dp, -1.6_dp, 1.6_dp, 1.9_dp, 2.2_dp, 2.5_dp, 2.8_dp, 3.0_dp]
        call make_rank_data(t, scores, y, dummy)
        call fclust(y, t, 2, cl, max_iter=50)
        call check(count(cl%cluster == cl%cluster(1)) == 6, 'fclust separates ordered two-group scores', failures)
        call k_cfc(y, t, 2, kc, max_iter=20)
        call check(size(kc%cluster) == 12, 'k_cfc returns all cluster labels', failures)
        call wfda(y, t, wr, lambda=1.0e-4_dp)
        call check(all(shape(wr%aligned) == shape(y)) .and. all(wr%costs >= 0.0_dp), 'wfda aligned output', failures)
    end subroutine test_ordering_and_clustering

    function sort_copy(x) result(y)
        integer, intent(in) :: x(:) !! Integer values copied and sorted in ascending order for permutation checks.
        integer, allocatable :: y(:)
        integer :: i
        integer :: j
        integer :: key

        y = x
        do i = 2, size(y)
            key = y(i)
            j = i - 1
            do while (j >= 1)
                if (y(j) <= key) exit
                y(j + 1) = y(j)
                j = j - 1
            end do
            y(j + 1) = key
        end do
    end function sort_copy

    subroutine test_regression_family(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        integer, parameter :: n = 24
        integer, parameter :: m = 11
        real(dp) :: t(m)
        real(dp) :: s1(n)
        real(dp) :: s2(n)
        real(dp) :: xfun(n, m, 2)
        real(dp) :: response(n)
        real(dp) :: q(3)
        real(dp) :: xeval(21, 2)
        real(dp) :: xobs(n, 2)
        real(dp) :: yfunc(n, m)
        type(sbf_result) :: sb
        type(fam_result) :: fr
        type(multifam_result) :: mf
        type(flm_result) :: lr
        type(flm_ci_result) :: lci
        type(fcreg_result) :: cr
        type(fopt_result) :: od
        type(fpc_quantile_result) :: qr
        integer :: i
        integer :: j

        t = [(real(j - 1, dp) / real(m - 1, dp), j=1, m)]
        do i = 1, n
            s1(i) = (real(i, dp) - 12.5_dp) / 5.0_dp
            s2(i) = sin(0.7_dp * real(i, dp))
            xfun(i, :, 1) = s1(i) * (1.0_dp + 0.5_dp * cos(acos(-1.0_dp) * t))
            xfun(i, :, 2) = s2(i) * (0.3_dp + sin(0.5_dp * acos(-1.0_dp) * t))
            response(i) = 2.0_dp + 1.5_dp * s1(i) - 0.8_dp * s2(i)
            xobs(i, :) = [s1(i), s2(i)]
            yfunc(i, :) = 1.0_dp + 2.0_dp * xfun(i, :, 1) - 1.2_dp * xfun(i, :, 2)
        end do
        do i = 1, 21
            xeval(i, 1) = -2.5_dp + 5.0_dp * real(i - 1, dp) / 20.0_dp
            xeval(i, 2) = -1.1_dp + 2.2_dp * real(i - 1, dp) / 20.0_dp
        end do
        call sb_fitting(response, xeval, xobs, sb, h=[1.0_dp, 0.8_dp])
        call check(all(shape(sb%fit) == [21, 2]), 'sb_fitting component shape', failures)
        call fam(response, xfun(:, :, 1), t, fr, n_eval=21, fve_threshold=0.9_dp)
        call check(size(fr%fam, 1) == 21, 'fam evaluation size', failures)
        call multi_fam(response, xfun, t, mf, n_eval=21, fve_threshold=0.9_dp)
        call check(size(mf%sbfit, 1) == 21, 'multi_fam evaluation size', failures)
        call flm_scalar_dense(response, xfun, t, lr, fve_threshold=0.9_dp)
        call check(lr%r2 > 0.99_dp, 'flm_scalar_dense explains synthetic linear response', failures)
        call seed_rng(789)
        call flm_ci_scalar_dense(response, xfun, t, lci, level=0.8_dp, n_boot=8, fve_threshold=0.9_dp)
        call check(all(lci%beta_lower <= lci%beta_upper), 'flm_ci_scalar_dense ordered intervals', failures)
        call fc_reg(xfun, yfunc, t, cr)
        call check(maxval(abs(cr%beta(1, :) - 2.0_dp)) < 1.0e-6_dp, 'fc_reg first concurrent coefficient', failures)
        call check(maxval(abs(cr%beta(2, :) + 1.2_dp)) < 1.0e-6_dp, 'fc_reg second concurrent coefficient', failures)
        call f_opt_des(xfun(:, :, 1), t, 2, 0.05_dp, od, response=response)
        call check(size(od%indices) == 2 .and. od%r2 >= 0.0_dp .and. od%r2 <= 1.0_dp, 'f_opt_des selection', failures)
        q = [0.25_dp, 0.5_dp, 0.75_dp]
        call fpc_quantile(xfun(:, :, 1), t, response, q, qr, fve_threshold=0.9_dp)
        call check(all(shape(qr%pred_quantile) == [n, 3]), 'fpc_quantile prediction shape', failures)
        call check(all(ieee_is_finite(qr%pred_quantile)), 'fpc_quantile finite predictions', failures)
    end subroutine test_regression_family

    subroutine test_longitudinal_family(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        integer, parameter :: n = 20
        integer, parameter :: m = 9
        integer, parameter :: d = 2
        real(dp) :: t(m)
        real(dp) :: y(n, m)
        real(dp) :: x(n, d)
        real(dp) :: grid_x(15, d)
        real(dp) :: pooled_t(n * m)
        real(dp) :: pooled_y(n * m)
        real(dp) :: pooled_x(n * m, d)
        real(dp) :: grid_t(5)
        type(vcam_result) :: vr
        type(tvam_result) :: tr
        integer :: i
        integer :: j
        integer :: pos

        t = [(real(j - 1, dp) / real(m - 1, dp), j=1, m)]
        do i = 1, n
            x(i, 1) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
            x(i, 2) = sin(0.5_dp * real(i, dp))
            y(i, :) = sin(acos(-1.0_dp) * t) + x(i, 1) * (0.5_dp + t) + x(i, 2) * (0.2_dp + t**2)
        end do
        do i = 1, 15
            grid_x(i, 1) = -1.0_dp + 2.0_dp * real(i - 1, dp) / 14.0_dp
            grid_x(i, 2) = -1.0_dp + 2.0_dp * real(i - 1, dp) / 14.0_dp
        end do
        call vcam(y, t, x, grid_x, vr, bandwidth=[0.6_dp, 0.6_dp])
        call check(all(shape(vr%fitted_y) == [n, m]), 'vcam fitted shape', failures)
        pos = 0
        do i = 1, n
            do j = 1, m
                pos = pos + 1
                pooled_t(pos) = t(j)
                pooled_y(pos) = y(i, j)
                pooled_x(pos, :) = x(i, :)
            end do
        end do
        grid_t = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
        call tvam(pooled_t, pooled_y, pooled_x, grid_t, grid_x, 0.26_dp, [0.7_dp, 0.7_dp], tr)
        call check(all(shape(tr%components) == [5, 15, 2]), 'tvam component surface shape', failures)
    end subroutine test_longitudinal_family

end program test_advanced
