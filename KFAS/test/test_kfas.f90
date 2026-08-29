! SPDX-License-Identifier: GPL-2.0-or-later

program test_kfas
    use kfas, only: dp, kfas_approx_result, kfas_approximate_nongaussian, &
        kfas_ar_transform, kfas_binomial, kfas_filter_result, kfas_gamma, &
        kfas_gaussian_filter, kfas_gaussian_loglik, kfas_gaussian_smooth, &
        kfas_init_theta, kfas_ldl_factor, kfas_model, kfas_negative_binomial, &
        kfas_nongaussian_loglik, kfas_poisson, kfas_smooth_result, &
        kfas_weighted_mean_cov
    implicit none

    call test_ar_transform()
    call test_ldl()
    call test_weighted_moments()
    call test_gaussian_filter()
    call test_gaussian_smoother()
    call test_correlated_observation_covariance()
    call test_diffuse_filter()
    call test_init_theta()
    call test_nongaussian_fixed_poisson()
    call test_density_bridge()
    print '(a)', "All KFAS tests passed."

contains

    subroutine assert_close(actual, expected, tol, label)
        real(dp), intent(in) :: actual, expected, tol
        character(len = *), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print '(a)', "FAIL: " // trim(label)
            print '(a,es24.16)', "  actual   = ", actual
            print '(a,es24.16)', "  expected = ", expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len = *), intent(in) :: label

        if (.not. condition) then
            print '(a)', "FAIL: " // trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine test_ar_transform()
        real(dp) :: phi(2)

        phi = [0.2_dp, -0.3_dp]
        call kfas_ar_transform(phi)
        call assert_close(phi(1), 0.26_dp, 1.0e-14_dp, "AR transform phi(1)")
        call assert_close(phi(2), -0.3_dp, 1.0e-14_dp, "AR transform phi(2)")
    end subroutine test_ar_transform

    subroutine test_ldl()
        real(dp) :: a(2, 2)
        integer :: info

        a = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
        call kfas_ldl_factor(a, info = info)
        call assert_true(info == 0, "LDL status")
        call assert_close(a(1, 1), 4.0_dp, 1.0e-14_dp, "LDL d1")
        call assert_close(a(2, 1), 0.5_dp, 1.0e-14_dp, "LDL l21")
        call assert_close(a(2, 2), 2.0_dp, 1.0e-14_dp, "LDL d2")
    end subroutine test_ldl

    subroutine test_weighted_moments()
        real(dp) :: x(1, 1, 3), w(3), meanx(1, 1), covx(1, 1, 1)

        x(1, 1, :) = [1.0_dp, 2.0_dp, 4.0_dp]
        w = [0.2_dp, 0.3_dp, 0.5_dp]
        call kfas_weighted_mean_cov(x, w, meanx, covx)
        call assert_close(meanx(1, 1), 2.8_dp, 1.0e-14_dp, "weighted mean")
        call assert_close(covx(1, 1, 1), 1.56_dp, 1.0e-14_dp, "weighted variance")
    end subroutine test_weighted_moments

    subroutine test_gaussian_filter()
        type(kfas_model) :: model
        type(kfas_filter_result) :: result
        real(dp), parameter :: expected_a(6) = [ &
            0.0_dp, 0.8_dp, 1.4545454545454546_dp, 1.4766355140186915_dp, &
            1.9573835480673933_dp, 2.442383123078554_dp]
        real(dp), parameter :: expected_f(5) = [ &
            1.25_dp, 0.55_dp, 0.48636363636363633_dp, &
            0.47149532710280373_dp, 0.4674430128840436_dp]
        integer :: info, i
        real(dp) :: loglik

        allocate(model%y(5, 1), model%missing(5, 1))
        allocate(model%z(1, 1, 1), model%h(1, 1, 1))
        allocate(model%tmat(1, 1, 1), model%rmat(1, 1, 1), model%q(1, 1, 1))
        allocate(model%a1(1), model%p1(1, 1), model%p1inf(1, 1))
        model%y(:, 1) = [1.0_dp, 2.0_dp, 1.5_dp, 2.5_dp, 3.0_dp]
        model%missing = 0
        model%z = 1.0_dp
        model%h = 0.25_dp
        model%tmat = 1.0_dp
        model%rmat = 1.0_dp
        model%q = 0.1_dp
        model%a1 = 0.0_dp
        model%p1 = 1.0_dp
        model%p1inf = 0.0_dp
        model%time_varying = 0
        model%diffuse_rank = 0

        call kfas_gaussian_filter(model, result, filter_signal = .true., info = info)
        call assert_true(info == 0, "Gaussian filter status")
        do i = 1, 6
            call assert_close(result%a(1, i), expected_a(i), 2.0e-13_dp, "filtered prediction a")
        end do
        do i = 1, 5
            call assert_close(result%f(1, i), expected_f(i), 2.0e-13_dp, "innovation variance f")
        end do
        call assert_close(result%loglik, -7.275349449509625_dp, 2.0e-13_dp, "Gaussian log-likelihood")
        call assert_close(result%signal(1, 1), 0.0_dp, 2.0e-13_dp, "filtered signal first")
        call assert_close(result%signal(5, 1), expected_a(5), 2.0e-13_dp, "filtered signal fifth")

        loglik = kfas_gaussian_loglik(model, info = info)
        call assert_true(info == 0, "Gaussian loglik wrapper status")
        call assert_close(loglik, result%loglik, 1.0e-14_dp, "Gaussian loglik wrapper")
    end subroutine test_gaussian_filter

    subroutine test_gaussian_smoother()
        type(kfas_model) :: model
        type(kfas_smooth_result) :: result
        real(dp), parameter :: expected_state(5) = [ &
            1.401547757871303_dp, 1.7023216368069543_dp, 1.8840241704653873_dp, &
            2.2193363723099755_dp, 2.442383123078554_dp]
        real(dp), parameter :: expected_var(5) = [ &
            0.1043358422559101_dp, 0.08475564507579773_dp, 0.08131029364995229_dp, &
            0.08793596946888582_dp, 0.11629386197392134_dp]
        integer :: info, i

        allocate(model%y(5, 1), model%missing(5, 1))
        allocate(model%z(1, 1, 1), model%h(1, 1, 1))
        allocate(model%tmat(1, 1, 1), model%rmat(1, 1, 1), model%q(1, 1, 1))
        allocate(model%a1(1), model%p1(1, 1), model%p1inf(1, 1))
        model%y(:, 1) = [1.0_dp, 2.0_dp, 1.5_dp, 2.5_dp, 3.0_dp]
        model%missing = 0
        model%z = 1.0_dp
        model%h = 0.25_dp
        model%tmat = 1.0_dp
        model%rmat = 1.0_dp
        model%q = 0.1_dp
        model%a1 = 0.0_dp
        model%p1 = 1.0_dp
        model%p1inf = 0.0_dp
        model%time_varying = 0
        model%diffuse_rank = 0

        call kfas_gaussian_smooth(model, result, info = info)
        call assert_true(info == 0, "Gaussian smoother status")
        do i = 1, 5
            call assert_close(result%state(1, i), expected_state(i), 3.0e-13_dp, &
                "smoothed state")
            call assert_close(result%state_var(1, 1, i), expected_var(i), 3.0e-13_dp, &
                "smoothed state variance")
            call assert_close(result%signal(1, i), expected_state(i), 3.0e-13_dp, &
                "smoothed signal")
        end do
        call assert_close(result%loglik, -7.275349449509625_dp, 2.0e-13_dp, &
            "smoother log-likelihood")
    end subroutine test_gaussian_smoother

    subroutine test_correlated_observation_covariance()
        type(kfas_model) :: model
        type(kfas_filter_result) :: result
        integer :: info

        allocate(model%y(3, 2), model%missing(3, 2))
        allocate(model%z(2, 1, 1), model%h(2, 2, 1))
        allocate(model%tmat(1, 1, 1), model%rmat(1, 1, 1), model%q(1, 1, 1))
        allocate(model%a1(1), model%p1(1, 1), model%p1inf(1, 1))
        model%y = reshape([1.0_dp, 2.0_dp, 1.4_dp, 1.2_dp, 1.7_dp, 1.5_dp], [3, 2])
        model%missing = 0
        model%z = 1.0_dp
        model%h(:, :, 1) = reshape([0.4_dp, 0.1_dp, 0.1_dp, 0.3_dp], [2, 2])
        model%tmat = 1.0_dp
        model%rmat = 1.0_dp
        model%q = 0.2_dp
        model%a1 = 0.0_dp
        model%p1 = 1.0_dp
        model%p1inf = 0.0_dp
        model%time_varying = 0
        model%diffuse_rank = 0

        call kfas_gaussian_filter(model, result, filter_signal = .true., info = info)
        call assert_true(info == 0, "correlated-H filter status")
        call assert_true(result%transformed_h, "correlated-H LDL transformation")
        call assert_close(result%loglik, -5.360188606476303_dp, 3.0e-13_dp, &
            "correlated-H multivariate log-likelihood")
        call assert_close(result%att(1, 1), 0.9180327868852461_dp, 3.0e-13_dp, &
            "correlated-H filtered state")
        call assert_close(result%signal(2, 1), 0.9180327868852461_dp, 3.0e-13_dp, &
            "correlated-H original-scale signal")

        model%missing(2, 2) = 1
        call kfas_gaussian_filter(model, result, info = info)
        call assert_true(info == 0, "correlated-H partial-missing status")
        call assert_close(result%loglik, -4.9482252090428895_dp, 4.0e-13_dp, &
            "correlated-H partial-missing log-likelihood")
    end subroutine test_correlated_observation_covariance

    subroutine test_diffuse_filter()
        type(kfas_model) :: model
        type(kfas_filter_result) :: result
        integer :: info

        allocate(model%y(5, 1), model%missing(5, 1))
        allocate(model%z(1, 1, 1), model%h(1, 1, 1))
        allocate(model%tmat(1, 1, 1), model%rmat(1, 1, 1), model%q(1, 1, 1))
        allocate(model%a1(1), model%p1(1, 1), model%p1inf(1, 1))
        model%y(:, 1) = [1.0_dp, 2.0_dp, 1.5_dp, 2.5_dp, 3.0_dp]
        model%missing = 0
        model%z = 1.0_dp
        model%h = 0.25_dp
        model%tmat = 1.0_dp
        model%rmat = 1.0_dp
        model%q = 0.1_dp
        model%a1 = 0.0_dp
        model%p1 = 0.0_dp
        model%p1inf = 1.0_dp
        model%time_varying = 0
        model%diffuse_rank = 1

        call kfas_gaussian_filter(model, result, info = info)
        call assert_true(info == 0, "diffuse filter status")
        call assert_true(result%diffuse_end == 1, "diffuse phase length")
        call assert_true(result%remaining_diffuse_rank == 0, "diffuse rank exhausted")
        call assert_close(result%finf(1, 1), 1.0_dp, 1.0e-14_dp, "diffuse innovation variance")
        call assert_close(result%loglik, -5.2047353029994365_dp, 2.0e-13_dp, &
            "exact diffuse log-likelihood")
    end subroutine test_diffuse_filter

    subroutine test_init_theta()
        real(dp) :: y(2, 4), u(2, 4), theta(2, 4)
        integer :: missing(2, 4), distribution(4)

        y = 0.0_dp
        u = 1.0_dp
        missing = 0
        distribution = [kfas_poisson, kfas_binomial, kfas_gamma, kfas_negative_binomial]
        y(:, 1) = [0.0_dp, 2.0_dp]
        u(:, 1) = [1.0_dp, 4.0_dp]
        y(:, 2) = [2.0_dp, 5.0_dp]
        u(:, 2) = [4.0_dp, 10.0_dp]
        y(:, 3) = [0.5_dp, 3.0_dp]
        y(:, 4) = [0.01_dp, 2.0_dp]
        missing(2, 2) = 1

        call kfas_init_theta(y, u, distribution, missing, theta)
        call assert_close(theta(1, 1), log(0.1_dp), 1.0e-14_dp, "Poisson init theta floor")
        call assert_close(theta(2, 1), log(0.5_dp), 1.0e-14_dp, "Poisson init theta")
        call assert_close(theta(1, 2), 0.0_dp, 1.0e-14_dp, "binomial init theta")
        call assert_close(theta(2, 2), log((1.0_dp / 11.0_dp) / (10.0_dp / 11.0_dp)), &
            1.0e-14_dp, "binomial missing init theta")
        call assert_close(theta(1, 3), 0.0_dp, 1.0e-14_dp, "gamma init theta floor")
        call assert_close(theta(2, 3), log(3.0_dp), 1.0e-14_dp, "gamma init theta")
        call assert_close(theta(1, 4), log(1.0_dp / 6.0_dp), 1.0e-14_dp, &
            "negative-binomial init theta floor")
        call assert_close(theta(2, 4), log(2.0_dp), 1.0e-14_dp, &
            "negative-binomial init theta")
    end subroutine test_init_theta

    subroutine test_nongaussian_fixed_poisson()
        type(kfas_model) :: model
        type(kfas_approx_result) :: approximation
        real(dp) :: u(3, 1), theta(3, 1), loglik, expected
        integer :: distribution(1), info

        allocate(model%y(3, 1), model%missing(3, 1))
        allocate(model%z(1, 1, 1), model%h(1, 1, 1))
        allocate(model%tmat(1, 1, 1), model%rmat(1, 1, 1), model%q(1, 1, 1))
        allocate(model%a1(1), model%p1(1, 1), model%p1inf(1, 1))
        model%y(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
        model%missing = 0
        model%z = 1.0_dp
        model%h = 0.0_dp
        model%tmat = 1.0_dp
        model%rmat = 1.0_dp
        model%q = 0.0_dp
        model%a1 = log(2.0_dp)
        model%p1 = 0.0_dp
        model%p1inf = 0.0_dp
        model%time_varying = 0
        model%diffuse_rank = 0
        u = 1.0_dp
        distribution = kfas_poisson
        theta = log(2.0_dp)

        call kfas_nongaussian_loglik(model, u, distribution, loglik, theta = theta, info = info)
        expected = sum(model%y(:, 1) * log(2.0_dp) - 2.0_dp) &
            - log_gamma(2.0_dp) - log_gamma(3.0_dp) - log_gamma(4.0_dp)
        call assert_true(info == 0, "non-Gaussian loglik status")
        call assert_close(loglik, expected, 2.0e-13_dp, "fixed-state Poisson log-likelihood")

        call kfas_approximate_nongaussian(model, u, distribution, approximation, &
            theta = theta, info = info)
        call assert_true(info == 0, "non-Gaussian approximation status")
        call assert_true(approximation%iterations == 2, "non-Gaussian approximation iterations")
        call assert_close(approximation%difference, 0.0_dp, 1.0e-14_dp, &
            "non-Gaussian approximation difference")
        call assert_close(approximation%theta(1, 1), log(2.0_dp), 1.0e-14_dp, &
            "non-Gaussian approximation theta")
    end subroutine test_nongaussian_fixed_poisson

    subroutine test_density_bridge()
        interface
            subroutine dnormf(x, mu, sigma, res)
                import :: dp
                real(dp), intent(in) :: x, mu, sigma
                real(dp), intent(inout) :: res
            end subroutine dnormf
            subroutine dpoisf(x, lambda, res)
                import :: dp
                real(dp), intent(in) :: x, lambda
                real(dp), intent(inout) :: res
            end subroutine dpoisf
            subroutine dbinomf(x, n, prob, res)
                import :: dp
                real(dp), intent(in) :: x, n, prob
                real(dp), intent(inout) :: res
            end subroutine dbinomf
            subroutine dgammaf(x, shape, scale, res)
                import :: dp
                real(dp), intent(in) :: x, shape, scale
                real(dp), intent(inout) :: res
            end subroutine dgammaf
            subroutine dnbinomf(x, size, mu, res)
                import :: dp
                real(dp), intent(in) :: x, size, mu
                real(dp), intent(inout) :: res
            end subroutine dnbinomf
        end interface
        real(dp) :: res
        real(dp), parameter :: pi = acos(-1.0_dp)

        res = 0.0_dp
        call dnormf(0.0_dp, 0.0_dp, 1.0_dp, res)
        call assert_close(res, 0.5_dp * log(2.0_dp * pi), 1.0e-14_dp, "normal bridge")

        res = 0.0_dp
        call dpoisf(2.0_dp, 3.0_dp, res)
        call assert_close(res, 2.0_dp * log(3.0_dp) - 3.0_dp - log(2.0_dp), &
            1.0e-14_dp, "Poisson bridge")

        res = 0.0_dp
        call dbinomf(2.0_dp, 5.0_dp, 0.4_dp, res)
        call assert_close(res, log(10.0_dp) + 2.0_dp * log(0.4_dp) + &
            3.0_dp * log(0.6_dp), 1.0e-14_dp, "binomial bridge")

        res = 0.0_dp
        call dgammaf(2.0_dp, 3.0_dp, 2.0_dp, res)
        call assert_close(res, 2.0_dp * log(2.0_dp) - 1.0_dp - log_gamma(3.0_dp) &
            - 3.0_dp * log(2.0_dp), 1.0e-14_dp, "gamma bridge")

        res = 0.0_dp
        call dnbinomf(2.0_dp, 4.0_dp, 3.0_dp, res)
        call assert_close(res, log_gamma(6.0_dp) - log_gamma(4.0_dp) - log_gamma(3.0_dp) &
            + 4.0_dp * log(4.0_dp / 7.0_dp) + 2.0_dp * log(3.0_dp / 7.0_dp), &
            1.0e-14_dp, "negative-binomial bridge")
    end subroutine test_density_bridge

end program test_kfas
