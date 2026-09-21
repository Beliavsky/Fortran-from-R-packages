program test_fdapace
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use fdapace, only : bwnn_result, bw_nn, convert_support, create_basis, cumtrapz_rcpp, dp, &
                        dyn_corr, dyn_test, dyn_test_result, fc_cor, fccor_result, fpca_dense, &
                        fpca_inputs, fpca_result, fsvd_dense, fsvd_result, get_normalised_sample, &
                        gp_functional_data, lwls1d, lwls2d, lwls2d_deriv, make_bw_to_zscore_02y, &
                        make_fpca_inputs_dense, make_gp_functional_data, make_hc_to_zscore_02y, &
                        make_ln_to_zscore_02y, norm_curv_to_area, real_vector, seed_rng, select_k_fve, &
                        select_k_result, sparse_sample, sparsify, trapz_rcpp, wiener
    implicit none

    integer :: failures

    failures = 0
    call test_integration(failures)
    call test_basis_and_support(failures)
    call test_smoothers(failures)
    call test_dynamic(failures)
    call test_bwnn(failures)
    call test_fccor(failures)
    call test_growth(failures)
    call test_simulation(failures)
    call test_fpca_and_fsvd(failures)

    if (failures /= 0) then
        write (*, '(a,i0)') 'fdapace tests failed: ', failures
        error stop 1
    end if
    write (*, '(a)') 'All fdapace tests passed.'

contains

    subroutine check(condition, message, failures)
        logical, intent(in) :: condition !! Test condition that must be true.
        character(len=*), intent(in) :: message !! Diagnostic text printed when the condition is false.
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

    subroutine test_integration(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: x(3)
        real(dp) :: y(3)
        real(dp), allocatable :: c(:)
        real(dp), allocatable :: scaled(:)

        x = [0.0_dp, 0.5_dp, 1.0_dp]
        y = x
        call check_close(trapz_rcpp(x, y), 0.5_dp, 1.0e-13_dp, 'trapz_rcpp integrates a line', failures)
        c = cumtrapz_rcpp(x, y)
        call check_close(c(3), 0.5_dp, 1.0e-13_dp, 'cumtrapz_rcpp terminal value', failures)
        scaled = norm_curv_to_area([1.0_dp, 1.0_dp, 1.0_dp], x, 2.0_dp)
        call check_close(trapz_rcpp(x, scaled), 2.0_dp, 1.0e-13_dp, 'norm_curv_to_area target', failures)
    end subroutine test_integration

    subroutine test_basis_and_support(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: from_grid(3)
        real(dp) :: to_grid(2)
        real(dp) :: mu(3)
        real(dp), allocatable :: basis(:,:)
        real(dp), allocatable :: vals(:)

        from_grid = [0.0_dp, 0.5_dp, 1.0_dp]
        to_grid = [0.25_dp, 0.75_dp]
        mu = [0.0_dp, 1.0_dp, 2.0_dp]
        vals = convert_support(from_grid, to_grid, mu)
        call check_close(vals(1), 0.5_dp, 1.0e-13_dp, 'convert_support first interpolation', failures)
        call check_close(vals(2), 1.5_dp, 1.0e-13_dp, 'convert_support second interpolation', failures)
        basis = create_basis(3, from_grid, 'cos')
        call check(all(abs(basis(:, 1) - 1.0_dp) < 1.0e-13_dp), 'create_basis constant cosine column', failures)
        call check(size(basis, 1) == 3 .and. size(basis, 2) == 3, 'create_basis shape', failures)
    end subroutine test_basis_and_support

    subroutine test_smoothers(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: xin1(5)
        real(dp) :: yin1(5)
        real(dp) :: xout1(2)
        real(dp) :: xin2(9, 2)
        real(dp) :: yin2(9)
        real(dp) :: gx(2)
        real(dp) :: gy(2)
        real(dp) :: bw2(2)
        real(dp), allocatable :: fit1(:)
        real(dp), allocatable :: fit2(:,:)
        real(dp), allocatable :: dx(:,:)
        integer :: i
        integer :: j
        integer :: q

        xin1 = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        yin1 = 2.0_dp + 3.0_dp * xin1
        xout1 = [0.2_dp, 0.8_dp]
        fit1 = lwls1d(2.0_dp, 'rect', xin1, yin1, xout1, npoly=1, nder=0)
        call check_close(fit1(1), 2.6_dp, 1.0e-11_dp, 'lwls1d recovers a line at first point', failures)
        call check_close(fit1(2), 4.4_dp, 1.0e-11_dp, 'lwls1d recovers a line at second point', failures)

        q = 0
        do j = 0, 2
            do i = 0, 2
                q = q + 1
                xin2(q, 1) = 0.5_dp * real(i, dp)
                xin2(q, 2) = 0.5_dp * real(j, dp)
                yin2(q) = 1.0_dp + 2.0_dp * xin2(q, 1) + 3.0_dp * xin2(q, 2)
            end do
        end do
        gx = [0.25_dp, 0.75_dp]
        gy = [0.25_dp, 0.75_dp]
        bw2 = [2.0_dp, 2.0_dp]
        fit2 = lwls2d(bw2, 'rect', xin2, yin2, gx, gy, crosscov=.true.)
        call check_close(fit2(1, 1), 2.25_dp, 1.0e-10_dp, 'lwls2d recovers a plane', failures)
        dx = lwls2d_deriv(bw2, 'rect', xin2, yin2, gx, gy, npoly=1, nder1=1, nder2=0)
        call check_close(dx(2, 1), 2.0_dp, 1.0e-10_dp, 'lwls2d_deriv recovers first derivative', failures)
    end subroutine test_smoothers

    subroutine test_dynamic(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: t(5)
        real(dp) :: x(5, 5)
        real(dp) :: y(5, 5)
        real(dp), allocatable :: corr(:)
        type(dyn_test_result) :: bt
        integer :: i

        t = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        do i = 1, 5
            x(i, :) = real(i, dp) * [0.0_dp, 1.0_dp, -0.5_dp, 0.7_dp, -0.2_dp] &
                      + 0.1_dp * real(i * i, dp) * [1.0_dp, 0.0_dp, 0.4_dp, -0.3_dp, 0.2_dp]
            y(i, :) = (0.5_dp + 0.15_dp * real(i, dp)) * x(i, :) &
                      + 0.2_dp * real(i - 3, dp) * [0.0_dp, -0.2_dp, 0.5_dp, 0.1_dp, -0.4_dp]
        end do
        corr = dyn_corr(x, y, t)
        call check(all(ieee_is_finite(corr)), 'dyn_corr produces finite individual correlations', failures)
        call seed_rng(12345)
        call dyn_test(x, y, t, bt, b=30)
        call check(ieee_is_finite(bt%stats), 'dyn_test finite statistic', failures)
        call check(bt%pval >= 0.0_dp .and. bt%pval <= 1.0_dp, 'dyn_test p-value range', failures)
    end subroutine test_dynamic

    subroutine test_bwnn(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        type(real_vector) :: lt(5)
        type(bwnn_result) :: bw

        lt(1)%v = [1.0_dp, 7.0_dp]
        lt(2)%v = [2.0_dp, 3.0_dp]
        lt(3)%v = [6.0_dp]
        lt(4)%v = [2.0_dp, 4.0_dp]
        lt(5)%v = [4.0_dp, 5.0_dp]
        bw = bw_nn(lt, k=2)
        call check_close(bw%cov_bw, 3.0_dp, 1.0e-13_dp, 'bw_nn covariance example', failures)
        call check_close(bw%mu_bw, 2.0_dp, 1.0e-13_dp, 'bw_nn mean example', failures)
    end subroutine test_bwnn

    subroutine test_fccor(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        type(real_vector) :: x(4)
        type(real_vector) :: y(4)
        type(real_vector) :: t(4)
        type(fccor_result) :: res
        integer :: i

        do i = 1, 4
            t(i)%v = [0.0_dp, 0.5_dp, 1.0_dp]
            x(i)%v = real(i, dp) + [0.0_dp, 0.5_dp, 1.0_dp] + 0.1_dp * real(i * i, dp) * [0.2_dp, -0.1_dp, 0.3_dp]
            y(i)%v = 0.7_dp * x(i)%v + 0.2_dp * real(i, dp) * [0.3_dp, -0.2_dp, 0.1_dp]
        end do
        res = fc_cor(x, y, t, [1.0_dp], 'gauss')
        call check(size(res%corr) == 3, 'fc_cor output grid size', failures)
        call check(count(.not. ieee_is_nan(res%corr)) >= 1, 'fc_cor has at least one finite correlation', failures)
    end subroutine test_fccor

    subroutine test_growth(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp), allocatable :: z(:)

        z = make_hc_to_zscore_02y('F', [0.0_dp], [33.8787_dp])
        call check_close(z(1), 0.0_dp, 1.0e-12_dp, 'head circumference reference mean gives z=0', failures)
        z = make_ln_to_zscore_02y('F', [0.0_dp], [49.1477_dp])
        call check_close(z(1), 0.0_dp, 1.0e-12_dp, 'length reference mean gives z=0', failures)
        z = make_bw_to_zscore_02y('F', [0.0_dp], [3.2322_dp])
        call check_close(z(1), 0.0_dp, 1.0e-12_dp, 'weight reference median gives z=0', failures)
    end subroutine test_growth

    subroutine test_simulation(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        type(gp_functional_data) :: gp
        type(fpca_inputs) :: inputs
        type(sparse_sample) :: sparse
        real(dp), allocatable :: w(:,:)
        real(dp) :: mu(20)
        integer :: i

        mu = 0.0_dp
        call seed_rng(2026)
        gp = make_gp_functional_data(6, 20, mu, 2, [2.0_dp, 0.5_dp], 0.01_dp, 'cos')
        call check(all(shape(gp%y) == [6, 20]), 'make_gp_functional_data Y shape', failures)
        call check(gp%has_noisy .and. all(shape(gp%yn) == [6, 20]), 'make_gp_functional_data noisy output', failures)
        call check_close(sum(gp%xi(:, 1)) / 6.0_dp, 0.0_dp, 1.0e-12_dp, 'GP scores are centered', failures)
        inputs = make_fpca_inputs_dense(gp%pts, gp%y)
        call check(size(inputs%ly) == 6 .and. size(inputs%lt(1)%v) == 20, 'make_fpca_inputs_dense shape', failures)
        call seed_rng(7)
        sparse = sparsify(gp%y, gp%pts, [3])
        do i = 1, 6
            call check(size(sparse%lt(i)%v) == 3, 'sparsify fixed retained count', failures)
        end do
        call seed_rng(11)
        w = wiener(4, gp%pts, 8)
        call check(all(shape(w) == [4, 20]), 'wiener output shape', failures)
    end subroutine test_simulation

    subroutine test_fpca_and_fsvd(failures)
        integer, intent(inout) :: failures !! Running number of failed assertions.
        real(dp) :: t(9)
        real(dp) :: y1(8, 9)
        real(dp) :: y2(8, 9)
        real(dp) :: scores(8)
        real(dp), allocatable :: ynorm(:,:)
        type(fpca_result) :: fp
        type(fsvd_result) :: sv
        type(select_k_result) :: selected
        integer :: i

        do i = 1, 9
            t(i) = real(i - 1, dp) / 8.0_dp
        end do
        scores = [-2.0_dp, -1.4_dp, -0.8_dp, -0.2_dp, 0.3_dp, 0.9_dp, 1.5_dp, 2.1_dp]
        do i = 1, 8
            y1(i, :) = 1.0_dp + 0.5_dp * t + scores(i) * (1.0_dp + 0.5_dp * cos(acos(-1.0_dp) * t))
            y2(i, :) = -0.2_dp + 0.3_dp * t + 1.7_dp * scores(i) * sin(0.5_dp * acos(-1.0_dp) * t)
        end do
        call fpca_dense(y1, t, fp, fve_threshold=0.90_dp, max_k=5, assume_error=.false.)
        call check(fp%select_k == 1, 'dense fpca selects rank-one component', failures)
        call check(fp%lambda(1) > 0.0_dp, 'dense fpca positive leading eigenvalue', failures)
        call check(maxval(abs(fp%fitted_y - y1)) < 1.0e-8_dp, 'dense fpca reconstructs rank-one sample', failures)
        selected = select_k_fve(fp, 0.90_dp)
        call check(selected%k == 1, 'select_k_fve selects first component', failures)
        ynorm = get_normalised_sample(fp, y1, t)
        call check(all(shape(ynorm) == shape(y1)), 'get_normalised_sample shape', failures)
        call fsvd_dense(y1, t, y2, t, sv, fve_threshold=0.90_dp, max_k=5, fixed_k=1, assume_error=.false.)
        call check(sv%nsvd == 1 .and. sv%s_values(1) > 0.0_dp, 'dense fsvd leading singular component', failures)
        call check(abs(sv%can_corr(1)) <= 1.0_dp + 1.0e-12_dp, 'fsvd canonical correlation bounded', failures)
    end subroutine test_fpca_and_fsvd

end program test_fdapace
