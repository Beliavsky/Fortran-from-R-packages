program test_generalized
    use fdausc
    implicit none

    call test_basis_regression()
    call test_glm()
    call test_gls_and_partial()
    call test_additive_and_deviance()
    call test_influence_and_projection()
    call test_penalized_and_additive()
    call test_selection_and_bootstrap()
    print '(a)', 'All fda.usc generalized tests passed'

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Actual scalar value produced by the routine under test.
        real(dp), intent(in) :: expected !! Expected deterministic scalar reference value.
        real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference.
        character(len=*), intent(in) :: label !! Test label printed if the comparison fails.

        if (abs(actual - expected) > tolerance) then
            print '(a,2(1x,es16.8))', trim(label), actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Condition that must evaluate true for the test to pass.
        character(len=*), intent(in) :: label !! Test label printed if the condition fails.

        if (.not. condition) then
            print '(a)', trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine test_basis_regression()
        real(dp) :: curves(5, 2)
        real(dp) :: basis(2, 2)
        real(dp) :: y(5)
        real(dp) :: coefficients(3)
        real(dp) :: beta(2)
        real(dp) :: fitted(5)
        real(dp) :: residuals(5)
        real(dp) :: hat(5, 5)
        integer :: info

        curves = reshape([0.0_dp, 0.0_dp, &
                          1.0_dp, 0.0_dp, &
                          0.0_dp, 1.0_dp, &
                          1.0_dp, 1.0_dp, &
                          2.0_dp, 1.0_dp], shape(curves), order=[2, 1])
        basis = 0.0_dp
        basis(1, 1) = 1.0_dp
        basis(2, 2) = 1.0_dp
        y = 1.0_dp + 2.0_dp*curves(:, 1) + 3.0_dp*curves(:, 2)
        call basis_scalar_regression(curves, basis, y, coefficients, beta, fitted, residuals, hat, center=.false., info=info)
        call assert_true(info == 0, 'basis_scalar_regression info')
        call assert_close(coefficients(1), 1.0_dp, 1.0e-10_dp, 'basis regression intercept')
        call assert_close(beta(1), 2.0_dp, 1.0e-10_dp, 'basis regression beta1')
        call assert_close(beta(2), 3.0_dp, 1.0e-10_dp, 'basis regression beta2')
        call assert_close(maxval(abs(residuals)), 0.0_dp, 1.0e-10_dp, 'basis regression residuals')
    end subroutine test_basis_regression

    subroutine test_glm()
        real(dp) :: x(6, 2)
        real(dp) :: y(6)
        real(dp) :: coefficients(2)
        real(dp) :: fitted(6)
        real(dp) :: residuals(6)
        real(dp) :: hat(6, 6)
        real(dp) :: deviance
        integer :: iterations
        integer :: info
        integer :: i
        logical :: converged

        do i = 1, 6
            x(i, 1) = 1.0_dp
            x(i, 2) = real(i - 3, dp)
        end do
        y = [0.10_dp, 0.20_dp, 0.35_dp, 0.60_dp, 0.78_dp, 0.90_dp]
        call glm_irls(x, y, 2, coefficients, fitted, residuals, hat, deviance, iterations, converged, info=info)
        call assert_true(info == 0, 'glm_irls info')
        call assert_true(converged, 'glm_irls convergence')
        call assert_true(coefficients(2) > 0.0_dp, 'glm_irls positive slope')
        call assert_true(all(fitted > 0.0_dp .and. fitted < 1.0_dp), 'glm_irls fitted range')
    end subroutine test_glm

    subroutine test_gls_and_partial()
        real(dp) :: x(5, 2)
        real(dp) :: y(5)
        real(dp) :: covariance(5, 5)
        real(dp) :: coefficients(2)
        real(dp) :: fitted(5)
        real(dp) :: residuals(5)
        real(dp) :: hat(5, 5)
        real(dp) :: smoother(5, 5)
        real(dp) :: nonparametric(5)
        integer :: i
        integer :: info

        covariance = 0.0_dp
        smoother = 0.0_dp
        do i = 1, 5
            x(i, 1) = 1.0_dp
            x(i, 2) = real(i - 1, dp)
            covariance(i, i) = 1.0_dp
        end do
        y = 2.0_dp + 0.5_dp*x(:, 2)
        call gls_fit(x, y, covariance, coefficients, fitted, residuals, hat, info=info)
        call assert_true(info == 0, 'gls_fit info')
        call assert_close(coefficients(1), 2.0_dp, 1.0e-10_dp, 'gls intercept')
        call assert_close(coefficients(2), 0.5_dp, 1.0e-10_dp, 'gls slope')
        call partially_linear_fit(smoother, x, y, coefficients, nonparametric, fitted, residuals, hat, info)
        call assert_true(info == 0, 'partially_linear_fit info')
        call assert_close(maxval(abs(residuals)), 0.0_dp, 1.0e-10_dp, 'partial residuals')
    end subroutine test_gls_and_partial

    subroutine test_additive_and_deviance()
        real(dp) :: smoothers(4, 4, 1)
        real(dp) :: y(4)
        real(dp) :: components(4, 1)
        real(dp) :: fitted(4)
        real(dp) :: residuals(4)
        real(dp) :: offset(4)
        real(dp) :: score
        real(dp) :: intercept
        integer :: i
        integer :: j
        integer :: iterations
        logical :: converged

        smoothers = -0.25_dp
        do i = 1, 4
            smoothers(i, i, 1) = 0.75_dp
        end do
        y = [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp]
        call additive_backfit(smoothers, y, intercept, components, fitted, residuals, iterations, converged)
        call assert_true(converged, 'additive_backfit convergence')
        call assert_close(maxval(abs(residuals)), 0.0_dp, 1.0e-10_dp, 'additive residuals')
        smoothers = 0.0_dp
        offset = 0.0_dp
        call deviance_smoother_score(y, smoothers(:, :, 1), y, 1, offset, 0.0_dp, 1, score)
        call assert_close(score, sum(y*y)/4.0_dp, 1.0e-10_dp, 'deviance smoother score')
        j = iterations
        call assert_true(j >= 1, 'additive iterations')
    end subroutine test_additive_and_deviance

    subroutine test_influence_and_projection()
        real(dp) :: dcp(2)
        real(dp) :: dce(2)
        real(dp) :: dp_stat(2)
        real(dp) :: boot(4)
        real(dp) :: q1(2)
        real(dp) :: q2(2)
        real(dp) :: q3(2)
        real(dp) :: x(3, 2)
        real(dp) :: y(3, 2)
        real(dp) :: projections(2, 2)
        real(dp) :: statistics(2)
        real(dp) :: pvalues(2)
        real(dp) :: fdr_pvalue
        real(dp) :: anova_stats(2)
        integer :: groups(6)
        real(dp) :: curves(6, 2)

        dcp = [1.0_dp, 3.0_dp]
        dce = [1.0_dp, 3.0_dp]
        dp_stat = [1.0_dp, 3.0_dp]
        boot = [0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        call influence_quantile_summary(dcp, dce, dp_stat, boot, boot, boot, q1, q2, q3)
        call assert_close(q1(1), 0.5_dp, 1.0e-12_dp, 'influence quantile 1')
        call assert_close(q1(2), 0.75_dp, 1.0e-12_dp, 'influence quantile 2')
        x = reshape([0.0_dp, 0.0_dp, 0.2_dp, 0.1_dp, 0.4_dp, 0.2_dp], shape(x), order=[2, 1])
        y = reshape([1.0_dp, 1.0_dp, 1.2_dp, 1.1_dp, 1.4_dp, 1.2_dp], shape(y), order=[2, 1])
        projections = 0.0_dp
        projections(1, 1) = 1.0_dp
        projections(2, 2) = 1.0_dp
        call projected_two_sample_ks(x, y, projections, statistics, pvalues, fdr_pvalue)
        call assert_true(all(statistics > 0.0_dp), 'projected KS statistics')
        call assert_true(fdr_pvalue >= 0.0_dp .and. fdr_pvalue <= 1.0_dp, 'projected KS FDR p-value')
        curves(1:3, :) = x
        curves(4:6, :) = y
        groups = [1, 1, 1, 2, 2, 2]
        call random_projection_anova(curves, projections, groups, anova_stats, heteroscedastic=.false.)
        call assert_true(all(anova_stats > 0.0_dp), 'random projection ANOVA')
    end subroutine test_influence_and_projection


    subroutine test_penalized_and_additive()
        real(dp) :: x(6, 2)
        real(dp) :: y(6)
        real(dp) :: penalty(2, 2)
        real(dp) :: coefficients(2)
        real(dp) :: fitted(6)
        real(dp) :: residuals(6)
        real(dp) :: deviance
        real(dp) :: smoothers(6, 6, 1)
        real(dp) :: components(6, 1)
        real(dp) :: intercept
        integer :: iterations
        integer :: info
        integer :: i
        logical :: converged

        do i = 1, 6
            x(i, 1) = 1.0_dp
            x(i, 2) = real(i - 3, dp)
        end do
        y = [0.10_dp, 0.20_dp, 0.35_dp, 0.60_dp, 0.78_dp, 0.90_dp]
        penalty = 0.0_dp
        penalty(2, 2) = 0.1_dp
        call penalized_glm_irls(x, y, 2, penalty, coefficients, fitted, residuals, deviance, iterations, &
                                converged, info=info)
        call assert_true(info == 0, 'penalized_glm_irls info')
        call assert_true(converged, 'penalized_glm_irls convergence')
        call assert_true(coefficients(2) > 0.0_dp, 'penalized_glm_irls slope')
        smoothers = -1.0_dp/6.0_dp
        do i = 1, 6
            smoothers(i, i, 1) = 5.0_dp/6.0_dp
        end do
        call generalized_additive_backfit(smoothers, y, 2, intercept, components, fitted, residuals, iterations, converged)
        call assert_true(all(fitted > 0.0_dp .and. fitted < 1.0_dp), 'generalized additive fitted range')
        call assert_true(iterations >= 1, 'generalized additive iterations')
    end subroutine test_penalized_and_additive

    subroutine test_selection_and_bootstrap()
        real(dp) :: x(6, 3)
        real(dp) :: y(6)
        integer :: groups(3)
        logical :: selected(2)
        integer :: order(2)
        integer :: nselected
        real(dp) :: deviance
        integer :: info
        integer :: i
        integer :: obs_index(2, 6)
        integer :: res_index(2, 6)
        real(dp) :: multipliers(2, 6)
        real(dp) :: fitted(6)
        real(dp) :: residuals(6)
        real(dp) :: coefficient_replicates(2, 2)
        real(dp) :: xlin(6, 2)

        do i = 1, 6
            x(i, 1) = 1.0_dp
            x(i, 2) = real(i - 1, dp)
            x(i, 3) = real(mod(i, 2), dp)
        end do
        y = 1.0_dp + 2.0_dp*x(:, 2)
        groups = [0, 1, 2]
        call forward_glm_select(x, y, groups, 1, selected, order, nselected, deviance, info=info)
        call assert_true(info == 0, 'forward_glm_select info')
        call assert_true(nselected == 1, 'forward_glm_select count')
        call assert_true(selected(1), 'forward_glm_select informative group')
        call assert_true(.not. selected(2), 'forward_glm_select noise group')
        xlin = x(:, 1:2)
        fitted = y
        residuals = [0.5_dp, -0.5_dp, 0.25_dp, -0.25_dp, 0.1_dp, -0.1_dp]
        do i = 1, 6
            obs_index(1, i) = i
            obs_index(2, i) = i
            res_index(1, i) = i
            res_index(2, i) = 7 - i
        end do
        multipliers = 0.0_dp
        call bootstrap_linear_refits(xlin, fitted, residuals, obs_index, res_index, multipliers, &
                                     coefficient_replicates, info=info)
        call assert_true(info == 0, 'bootstrap_linear_refits info')
        call assert_close(coefficient_replicates(1, 1), 1.0_dp, 1.0e-10_dp, 'bootstrap intercept')
        call assert_close(coefficient_replicates(1, 2), 2.0_dp, 1.0e-10_dp, 'bootstrap slope')
    end subroutine test_selection_and_bootstrap

end program test_generalized
