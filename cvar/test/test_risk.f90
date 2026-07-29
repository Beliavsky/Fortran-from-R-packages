! SPDX-License-Identifier: GPL-2.0-or-later
module test_risk_distributions
    use cvar, only : dp, student_t_cdf, student_t_quantile
    implicit none
contains
    pure function logistic_quantile(probability) result(value)
        real(dp), intent(in) :: probability
        real(dp) :: value
        value = log(probability / (1.0_dp - probability))
    end function logistic_quantile

    pure function logistic_cdf(value_in) result(value)
        real(dp), intent(in) :: value_in
        real(dp) :: value
        if (value_in >= 0.0_dp) then
            value = 1.0_dp / (1.0_dp + exp(-value_in))
        else
            value = exp(value_in) / (1.0_dp + exp(value_in))
        end if
    end function logistic_cdf

    pure function logistic_pdf(value_in) result(value)
        real(dp), intent(in) :: value_in
        real(dp) :: value, cdf_value
        cdf_value = logistic_cdf(value_in)
        value = cdf_value * (1.0_dp - cdf_value)
    end function logistic_pdf

    pure function t4_quantile(probability) result(value)
        real(dp), intent(in) :: probability
        real(dp) :: value
        value = student_t_quantile(probability, 4.0_dp)
    end function t4_quantile

    pure function t4_cdf(value_in) result(value)
        real(dp), intent(in) :: value_in
        real(dp) :: value
        value = student_t_cdf(value_in, 4.0_dp)
    end function t4_cdf
end module test_risk_distributions

program test_risk
    use test_risk_distributions, only : logistic_quantile, logistic_cdf, logistic_pdf, &
                                         t4_quantile, t4_cdf
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use cvar, only : dp, cvar_ok, cvar_invalid_probability, &
                     normal_pdf, normal_cdf, normal_quantile, &
                     student_t_cdf, student_t_quantile, &
                     std_student_t_cdf, std_student_t_quantile, &
                     ged_pdf, ged_cdf, ged_quantile, &
                     var_qf, var_cdf, var_sample, &
                     es_qf, es_cdf, es_pdf, es_sample
    implicit none

    real(dp), parameter :: p = 0.05_dp
    real(dp), parameter :: z_ref = -1.6448536269514729_dp
    real(dp), parameter :: es_ref = 2.0627128075074253_dp
    real(dp), parameter :: mu = 0.006408553_dp
    real(dp), parameter :: variance = 0.0004018977_dp
    real(dp), parameter :: sd = sqrt(variance)
    real(dp), parameter :: transformed_var_ref = 0.02621667904226832_dp
    real(dp), parameter :: transformed_es_ref = 0.034313278401780334_dp
    real(dp), parameter :: sample(10) = [ -4.0_dp, -3.0_dp, -2.0_dp, -1.0_dp, &
                                          0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
                                          4.0_dp, 5.0_dp ]
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: matrix_values(:, :)
    real(dp) :: x, y, z
    real(dp) :: sample_matrix(10, 2)
    integer :: status

    call assert_close(normal_quantile(p), z_ref, 2.0e-14_dp, "normal quantile")
    call assert_close(normal_cdf(z_ref), p, 2.0e-15_dp, "normal cdf")
    call assert_close(normal_pdf(z_ref) / p, es_ref, 2.0e-14_dp, "normal ES identity")

    x = var_qf(normal_quantile, p, status=status)
    call assert_status(status, cvar_ok, "VaR qf status")
    call assert_close(x, -z_ref, 3.0e-14_dp, "VaR from quantile")

    y = var_cdf(normal_cdf, p, tol=1.0e-13_dp, status=status)
    call assert_status(status, cvar_ok, "VaR cdf status")
    call assert_close(y, x, 2.0e-12_dp, "VaR cdf versus qf")

    x = es_qf(normal_quantile, p, tol=2.0e-11_dp, status=status)
    call assert_status(status, cvar_ok, "ES qf status")
    call assert_close(x, es_ref, 2.0e-9_dp, "ES from quantile")

    y = es_cdf(normal_cdf, p, tol=2.0e-8_dp, status=status)
    call assert_status(status, cvar_ok, "ES cdf status")
    call assert_close(y, es_ref, 2.0e-7_dp, "ES from cdf")

    z = es_pdf(normal_pdf, normal_quantile, p, tol=2.0e-11_dp, status=status)
    call assert_status(status, cvar_ok, "ES pdf status")
    call assert_close(z, es_ref, 3.0e-9_dp, "ES from pdf")

    x = var_qf(normal_quantile, p, intercept=mu, slope=sd)
    call assert_close(x, 0.02656646317059502_dp, 2.0e-14_dp, "location-scale VaR")
    y = var_qf(normal_quantile, p, intercept=mu, slope=sd, transf=.true.)
    call assert_close(y, transformed_var_ref, 2.0e-14_dp, "transformed VaR")
    z = es_qf(normal_quantile, p, intercept=mu, slope=sd, transf=.true., tol=2.0e-11_dp)
    call assert_close(z, transformed_es_ref, 2.0e-9_dp, "transformed ES")

    z = es_pdf(normal_pdf, normal_quantile, p, intercept=mu, slope=sd, &
               transf=.true., tol=2.0e-11_dp)
    call assert_close(z, transformed_es_ref, 2.0e-9_dp, "transformed PDF ES")

    x = var_sample(sample, 0.2_dp)
    y = es_sample(sample, 0.2_dp)
    call assert_close(x, 2.2_dp, 1.0e-14_dp, "empirical type-7 VaR")
    call assert_close(y, 3.5_dp, 1.0e-14_dp, "empirical ES")

    sample_matrix(:, 1) = sample
    sample_matrix(:, 2) = 2.0_dp * sample
    values = var_sample(sample_matrix, 0.2_dp)
    call assert_close(values(1), 2.2_dp, 1.0e-14_dp, "matrix VaR first column")
    call assert_close(values(2), 4.4_dp, 1.0e-14_dp, "matrix VaR second column")
    matrix_values = es_sample(sample_matrix, [0.2_dp, 0.3_dp])
    call assert_close(matrix_values(1, 1), 3.5_dp, 1.0e-14_dp, "matrix ES first column")
    call assert_close(matrix_values(1, 2), 7.0_dp, 1.0e-14_dp, "matrix ES second column")

    values = var_qf(normal_quantile, [0.10_dp, 0.05_dp, 0.01_dp])
    call assert_close(values(1), 1.2815515655446004_dp, 3.0e-14_dp, "vector VaR 10 percent")
    call assert_close(values(2), -z_ref, 3.0e-14_dp, "vector VaR 5 percent")
    call assert_close(values(3), 2.3263478740408408_dp, 3.0e-14_dp, "vector VaR 1 percent")

    x = es_qf(logistic_quantile, 0.05_dp, tol=2.0e-11_dp)
    y = es_cdf(logistic_cdf, 0.05_dp, tol=2.0e-8_dp)
    z = es_pdf(logistic_pdf, logistic_quantile, 0.05_dp, tol=2.0e-11_dp)
    call assert_close(y, x, 3.0e-7_dp, "callback cdf ES")
    call assert_close(z, x, 3.0e-9_dp, "callback pdf ES")

    x = es_qf(t4_quantile, 0.05_dp, tol=2.0e-10_dp)
    y = es_cdf(t4_cdf, 0.05_dp, tol=5.0e-8_dp)
    call assert_close(y, x, 2.0e-6_dp, "Student-t qf/cdf ES")

    call assert_close(student_t_quantile(0.025_dp, 4.0_dp), &
                      -2.7764451051977943_dp, 3.0e-11_dp, "Student-t quantile")
    call assert_close(student_t_cdf(-2.7764451051977943_dp, 4.0_dp), &
                      0.025_dp, 2.0e-13_dp, "Student-t cdf")
    call assert_close(std_student_t_quantile(0.025_dp, 5.0_dp), &
                      -1.9911641278965484_dp, 3.0e-11_dp, "standardized-t quantile")
    call assert_close(std_student_t_cdf(-1.9911641278965484_dp, 5.0_dp), &
                      0.025_dp, 2.0e-13_dp, "standardized-t cdf")

    call assert_close(ged_quantile(0.025_dp, 1.5_dp), &
                      -2.033146704578768_dp, 3.0e-11_dp, "GED quantile")
    call assert_close(ged_cdf(-2.033146704578768_dp, 1.5_dp), &
                      0.025_dp, 3.0e-13_dp, "GED cdf")
    call assert_close(ged_pdf(0.0_dp, 1.5_dp), &
                      0.47596665240714864_dp, 3.0e-14_dp, "GED density")

    x = var_qf(normal_quantile, 0.0_dp, status=status)
    call assert_status(status, cvar_invalid_probability, "invalid probability status")
    if (.not. ieee_is_nan(x)) error stop "invalid probability should return NaN"

    print '(a)', "test_risk: all tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        if (abs(actual - expected) > tolerance) then
            print '(a)', trim(label)
            print '(a,es24.16)', "actual:   ", actual
            print '(a,es24.16)', "expected: ", expected
            print '(a,es24.16)', "tolerance:", tolerance
            error stop "assert_close failed"
        end if
    end subroutine assert_close

    subroutine assert_status(actual, expected, label)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        if (actual /= expected) then
            print '(a,2(1x,i0))', trim(label), actual, expected
            error stop "assert_status failed"
        end if
    end subroutine assert_status

end program test_risk
