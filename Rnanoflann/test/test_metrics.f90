program test_metrics
    use rnanoflann, only: dp, metric_distance
    implicit none
    real(dp) :: x(3), y(3)

    x = [0.1_dp, 0.3_dp, 0.6_dp]
    y = [0.2_dp, 0.5_dp, 0.3_dp]

    call check(metric_distance(x, y, "euclidean"), 0.37416573867739417_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "euclidean", square=.true.), 0.14_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "hellinger"), 0.21683228641053510_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "manhattan"), 0.6_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "canberra"), 0.9166666666666666_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "kullback_leibler"), 0.3794239969771763_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "jensen_shannon"), 0.09322676863972956_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "itakura_saito"), 0.6108256237659906_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "bhattacharyya"), 0.04815741684776861_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "jeffries_matusita"), 0.30664716020214594_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "minimum"), 0.1_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "maximum"), 0.3_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "total_variation"), 0.3_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "sorensen"), 0.9166666666666666_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "cosine"), 0.8371385281611448_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "gower"), 0.2_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "minkowski", p=3.0_dp), 0.33019272488946266_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "soergel"), 0.4615384615384617_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "kulczynski"), 0.8571428571428573_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "wave_hedges"), 1.4_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "motyka"), 0.65_dp, 2.0e-14_dp)
    call check(metric_distance(x, y, "harmonic_mean"), 0.35_dp, 2.0e-14_dp)

    print '(a)', 'test_metrics: PASS'
contains
    subroutine check(value, expected, tol)
        real(dp), intent(in) :: value, expected, tol
        if (abs(value - expected) > tol) then
            print *, 'value=', value, ' expected=', expected
            error stop 1
        end if
    end subroutine check
end program test_metrics
