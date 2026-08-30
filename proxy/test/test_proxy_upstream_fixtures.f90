program test_proxy_upstream_fixtures
    use proxy, only: dp, proxy_simil_auto, proxy_simil_cross, proxy_row_sums_dist, proxy_row_dist, proxy_col_dist
    implicit none

    real(dp) :: x(5, 4)
    real(dp) :: y(5, 4)
    real(dp), allocatable :: auto(:, :)
    real(dp), allocatable :: cross(:, :)
    real(dp), allocatable :: sums(:)
    integer, allocatable :: rows(:)
    integer, allocatable :: cols(:)
    real(dp) :: packed(10)

    x(1, :) = [0.4595829_dp, 0.246283120_dp, 0.6662151_dp, 0.2658000_dp]
    x(2, :) = [0.8390691_dp, 0.002929016_dp, 0.3449861_dp, 0.4751930_dp]
    x(3, :) = [0.0377036_dp, 0.300615279_dp, 0.1399626_dp, 0.5143060_dp]
    x(4, :) = [0.7621676_dp, 0.963274350_dp, 0.4358104_dp, 0.7934745_dp]
    x(5, :) = [0.1211588_dp, 0.452265898_dp, 0.8005272_dp, 0.3062947_dp]

    y(1, :) = [0.3435311_dp, 0.45003283_dp, 0.1931588_dp, 0.22999466_dp]
    y(2, :) = [0.4417228_dp, 0.78191550_dp, 0.9028450_dp, 0.01728483_dp]
    y(3, :) = [0.3872093_dp, 0.45154308_dp, 0.7324979_dp, 0.52502253_dp]
    y(4, :) = [0.9119128_dp, 0.21720445_dp, 0.3986998_dp, 0.72623146_dp]
    y(5, :) = [0.9374529_dp, 0.04343485_dp, 0.8170873_dp, 0.69974714_dp]

    call proxy_simil_auto(x, 'cosine', auto)
    call assert_close(auto(2, 1), 0.8175510_dp, 2.0e-7_dp, 'upstream cosine auto 2,1')
    call assert_close(auto(5, 4), 0.7514756_dp, 2.0e-7_dp, 'upstream cosine auto 5,4')

    call proxy_simil_cross(x, y, 'cosine', cross)
    call assert_close(cross(1, 1), 0.8068091_dp, 2.0e-7_dp, 'upstream cosine cross 1,1')
    call assert_close(cross(4, 2), 0.7703657_dp, 2.0e-7_dp, 'upstream cosine cross 4,2')
    call assert_close(cross(5, 5), 0.7184058_dp, 2.0e-7_dp, 'upstream cosine cross 5,5')

    packed = [0.8390691_dp, 0.0377036_dp, 0.7621676_dp, 0.1211588_dp, 0.3006153_dp, &
              0.9632744_dp, 0.4522659_dp, 0.4358104_dp, 0.8005272_dp, 0.3062947_dp]
    call proxy_row_sums_dist(packed, 5, sums)
    call assert_close(sums(1), 1.760099_dp, 1.0e-6_dp, 'upstream packed row sum 1')
    call assert_close(sums(5), 1.680247_dp, 1.0e-6_dp, 'upstream packed row sum 5')

    call proxy_row_dist(5, rows)
    call proxy_col_dist(5, cols)
    if (any(rows /= [2, 3, 4, 5, 3, 4, 5, 4, 5, 5])) error stop 'upstream row.dist fixture failed'
    if (any(cols /= [1, 1, 1, 1, 2, 2, 2, 3, 3, 4])) error stop 'upstream col.dist fixture failed'

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Calculated scalar value being checked against an upstream saved-output fixture.
        real(dp), intent(in) :: expected !! Reference scalar copied from the corresponding upstream proxy `.Rout.save` output.
        real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference, allowing for printed fixture precision.
        character(len=*), intent(in) :: label !! Human-readable assertion label used to identify a parity failure.

        if (abs(actual - expected) > tolerance) then
            write (*, '(a,2es24.14)') 'FAIL '//trim(label)//': ', actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_proxy_upstream_fixtures
