module test_probplot_distribution
    use e1071_kinds, only: dp
    implicit none
contains
    pure function logistic_quantile(p) result(q)
        real(dp), intent(in) :: p !! Probability in (0,1) mapped to the standard logistic quantile.
        real(dp) :: q

        q = log(p / (1.0_dp - p))
    end function logistic_quantile
end module test_probplot_distribution

program test_probplot_custom
    use e1071, only: dp, probplot_distribution, probplot_result, probplot_custom
    use test_probplot_distribution, only: logistic_quantile
    implicit none

    type(probplot_distribution) :: distribution
    type(probplot_result) :: result
    real(dp) :: x(8)
    real(dp) :: expected_quartiles(2)

    x = [-3.0_dp, -1.0_dp, -0.5_dp, 0.0_dp, 0.4_dp, 0.9_dp, 1.7_dp, 3.0_dp]
    distribution%quantile => logistic_quantile
    call probplot_custom(x, distribution, result, [0.25_dp, 0.5_dp, 0.75_dp])
    expected_quartiles = [log(0.25_dp / 0.75_dp), log(0.75_dp / 0.25_dp)]
    if (maxval(abs(result%probability_quantiles([1, 3]) - expected_quartiles)) > 1.0e-12_dp) then
        error stop 'custom probability-plot quantiles failed'
    end if
    if (abs(result%probability_quantiles(2)) > 1.0e-14_dp) error stop 'custom probability-plot median failed'
    if (any(result%theoretical(2:) <= result%theoretical(:size(result%theoretical) - 1))) then
        error stop 'custom probability-plot theoretical values are not increasing'
    end if
    print '(a)', 'test_probplot_custom: PASS'
end program test_probplot_custom
