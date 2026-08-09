program test_reference_values
    use global_opt_tests, only : dp, n_benchmarks, benchmark_names, get_default_bounds, go_test
    implicit none
    real(dp), parameter :: expected(n_benchmarks) = [ &
        1.78915039087759844e+01_dp, 5.51135200399999619e+01_dp, 2.59920000000000240e+00_dp, &
        7.83244176458245761e+02_dp, 7.82782002272460431e+02_dp, 1.83352439960671845e+01_dp, &
        1.07601371125180250e+02_dp, 5.58182087090360710e+02_dp, 1.55291310699213647e+00_dp, &
        3.10582621398427294e+00_dp, 6.95665585895600915e+06_dp, 2.65629466581454744e-65_dp, &
        -1.01408826582952716e+00_dp, -4.42516722543033337e-33_dp, 3.96994642343062242e+03_dp, &
        6.62105068726344967e+01_dp, 3.17455361571281749e+01_dp, -4.62214460704315100e-01_dp, &
        -1.16863553356085914e+00_dp, -9.03488852583760993e-01_dp, 9.98459111668805605e-03_dp, &
        1.05034148278620858e+01_dp, 5.12880590541980865e+01_dp, 2.46465111901105551e+01_dp, &
        -2.65831324788571965e-01_dp, 1.20053860613229233e+01_dp, 1.68734913641134865e+00_dp, &
        3.63061545974124082e-13_dp, 5.07786127247360309e+03_dp, -6.75961157799564588e-02_dp, &
        1.32547604742947569e+04_dp, 1.97570250000000306e+03_dp, 1.33347818247153214e+01_dp, &
        1.51669524124699429e+00_dp, 5.09369140625000000e+03_dp, 1.18338079618551233e+06_dp, &
        1.99821973455273866e+05_dp, 4.07996533800000101e+07_dp, 1.00123026988186563e+01_dp, &
        5.22475040893504605e-01_dp, 6.29325936360478444e+00_dp, 9.48249433286350141e-01_dp, &
        -1.19422787861972643e+03_dp, -2.66488448280410894e+00_dp, -2.31101742928505960e+00_dp, &
        -2.54159435869462902e+00_dp, -7.66943887043436634e-01_dp, 1.88123563238400122e+05_dp, &
        -1.41919123561646803e-02_dp, -8.05641505282735174e-05_dp &
    ]
    real(dp), allocatable :: lower(:), upper(:), x(:)
    real(dp) :: got, tol
    integer :: i, status

    do i = 1, n_benchmarks
        call get_default_bounds(trim(benchmark_names(i)), lower, upper, status)
        if (status /= 0) error stop 'bounds lookup failed'
        allocate(x(size(lower)))
        x = lower + 0.37_dp*(upper-lower)
        got = go_test(x, trim(benchmark_names(i)), status=status)
        if (status /= 0) error stop 'go_test failed'
        tol = 2.0e-10_dp*max(1.0_dp, abs(expected(i)))
        if (abs(got-expected(i)) > tol) then
            write(*,'(a,1x,a,2(1x,es24.16))') 'mismatch', trim(benchmark_names(i)), got, expected(i)
            error stop 'reference mismatch'
        end if
        deallocate(lower, upper, x)
    end do
    print '(a)', 'test_reference_values: PASS'
end program test_reference_values
