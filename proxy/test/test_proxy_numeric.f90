program test_proxy_numeric
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use proxy, only: dp, euclidean_distance, manhattan_distance, supremum_distance, &
                     minkowski_distance, canberra_distance, wave_hedges_distance, &
                     soergel_distance, cosine_similarity, &
                     extended_jaccard_similarity, extended_dice_similarity, &
                     mahalanobis_distance, levenshtein_distance
    implicit none

    real(dp) :: x(3)
    real(dp) :: y(3)
    real(dp) :: covariance(3, 3)
    real(dp) :: nan

    x = [1.0_dp, 2.0_dp, 3.0_dp]
    y = [2.0_dp, 4.0_dp, 1.0_dp]
    call assert_close(euclidean_distance(x, y), 3.0_dp, 1.0e-12_dp, 'euclidean')
    call assert_close(manhattan_distance(x, y), 5.0_dp, 1.0e-12_dp, 'manhattan')
    call assert_close(supremum_distance(x, y), 2.0_dp, 1.0e-12_dp, 'supremum')
    call assert_close(minkowski_distance(x, y, 1.0_dp), 5.0_dp, 1.0e-12_dp, 'minkowski')
    call assert_close(minkowski_distance(x, y, ieee_value(0.0_dp, ieee_positive_inf)), &
                      1.0_dp, 1.0e-12_dp, 'minkowski large exponent')
    call assert_close(canberra_distance(x, y), 7.0_dp / 6.0_dp, 1.0e-12_dp, 'canberra')
    call assert_close(wave_hedges_distance(x, y), 0.75_dp, 1.0e-12_dp, 'wave hedges scalar min max')
    call assert_close(soergel_distance(x, y), 1.25_dp, 1.0e-12_dp, 'soergel scalar max')
    call assert_close(cosine_similarity(x, y), 13.0_dp / sqrt(14.0_dp * 21.0_dp), 1.0e-12_dp, 'cosine')
    call assert_close(extended_jaccard_similarity(x, y), 13.0_dp / 22.0_dp, 1.0e-12_dp, 'extended jaccard')
    call assert_close(extended_dice_similarity(x, y), 26.0_dp / 35.0_dp, 1.0e-12_dp, 'extended dice')

    covariance = 0.0_dp
    covariance(1, 1) = 1.0_dp
    covariance(2, 2) = 1.0_dp
    covariance(3, 3) = 1.0_dp
    call assert_close(mahalanobis_distance(x, y, covariance), 3.0_dp, 1.0e-12_dp, 'mahalanobis')

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    x = [1.0_dp, nan, 3.0_dp]
    y = [2.0_dp, 5.0_dp, 1.0_dp]
    call assert_close(euclidean_distance(x, y), sqrt(7.5_dp), 1.0e-12_dp, 'euclidean missing compensation')
    call assert_close(manhattan_distance(x, y), 4.5_dp, 1.0e-12_dp, 'manhattan missing compensation')

    call assert_close(levenshtein_distance('kitten', 'sitting'), 3.0_dp, 1.0e-12_dp, 'levenshtein')

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Calculated scalar result to compare with the deterministic reference value.
        real(dp), intent(in) :: expected !! Reference value expected from the corresponding upstream proxy formula.
        real(dp), intent(in) :: tolerance !! Maximum allowed absolute error for this deterministic numeric assertion.
        character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.

        if (abs(actual - expected) > tolerance) then
            write (*, '(a,2es24.14)') 'FAIL '//trim(label)//': ', actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_proxy_numeric
