program test_seed
    use rangen
    implicit none
    real(dp), allocatable :: a(:), b(:)

    call seed_all(123_i8)
    a = runif(20)
    call seed_all(123_i8)
    b = runif(20)
    if (maxval(abs(a - b)) > epsilon(1.0_dp)) error stop "seed reproducibility"

    call seed_all(123_i8)
    a = rnorm(20)
    call seed_all(123_i8)
    b = rnorm(20)
    if (maxval(abs(a - b)) > epsilon(1.0_dp)) error stop "normal seed reproducibility"

    print *, "test_seed: PASS"
end program test_seed
