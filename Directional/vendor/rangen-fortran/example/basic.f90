program basic
    use rangen
    implicit none
    real(dp), allocatable :: x(:), m(:,:)
    integer, allocatable :: ix(:)

    call seed_all(42_i8)
    x = rnorm(10000, 2.0_dp, 3.0_dp)
    print '(a,f10.5)', 'normal mean = ', sum(x) / real(size(x), dp)

    m = col_rgamma(5, 2, [2.0_dp, 5.0_dp], [1.0_dp, 2.0_dp])
    print '(a,2f10.5)', 'gamma first row = ', m(1,:)

    ix = sample_int(10, 5, .false.)
    print '(a,5i4)', 'sample = ', ix
end program basic
