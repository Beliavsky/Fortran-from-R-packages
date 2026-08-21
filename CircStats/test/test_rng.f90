program test_rng
    use circstats
    implicit none
    integer, parameter :: n=20000
    real(dp), allocatable :: x(:)
    real(dp) :: rho
    integer :: nseed
    integer, allocatable :: seed(:)
    allocate(x(n))
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed=314159
    call random_seed(put=seed)
    call rvm(n,1.0_dp,2.5_dp,x)
    if (minval(x)<0.0_dp .or. maxval(x)>=twopi) error stop 1
    rho=est_rho(x)
    if (abs(rho-a1(2.5_dp))>0.025_dp) error stop 1
    call rwrpcauchy(n,0.7_dp,0.6_dp,x)
    if (abs(est_rho(x)-0.6_dp)>0.025_dp) error stop 1
    call rwrpnorm(n,0.4_dp,0.7_dp,x)
    if (abs(est_rho(x)-0.7_dp)>0.025_dp) error stop 1
    call rcard(n,1.2_dp,-0.25_dp,x)
    if (minval(x)<0.0_dp .or. maxval(x)>=twopi) error stop 1
    call rtri(n,0.2_dp,x)
    if (minval(x)<-pi .or. maxval(x)>pi) error stop 1
    print *, "test_rng: PASS"
end program
