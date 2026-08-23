program test_matrix
    use rangen
    implicit none
    real(dp), allocatable :: x(:,:), y(:,:)
    real(dp), parameter :: mins(2) = [0.0_dp, 10.0_dp]
    real(dp), parameter :: maxs(2) = [1.0_dp, 12.0_dp]
    real(dp), parameter :: means(2) = [0.0_dp, 5.0_dp]
    real(dp), parameter :: sds(2) = [1.0_dp, 2.0_dp]
    integer, parameter :: n = 30000

    call seed_all(9876_i8)
    x = runif_mat(7, 3, -2.0_dp, 4.0_dp)
    if (size(x,1) /= 7 .or. size(x,2) /= 3) error stop "matrix dimensions"
    if (minval(x) < -2.0_dp .or. maxval(x) > 4.0_dp) error stop "matrix range"

    x = col_runif(n, 2, mins, maxs)
    if (abs(sum(x(:,1)) / n - 0.5_dp) > 0.01_dp) error stop "col uniform 1"
    if (abs(sum(x(:,2)) / n - 11.0_dp) > 0.02_dp) error stop "col uniform 2"

    y = col_rnorm(n, 2, means, sds)
    if (abs(sum(y(:,1)) / n) > 0.02_dp) error stop "col normal 1"
    if (abs(sum(y(:,2)) / n - 5.0_dp) > 0.03_dp) error stop "col normal 2"

    print *, "test_matrix: PASS"
end program test_matrix
