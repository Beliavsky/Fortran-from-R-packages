program test_stable_bootstrap
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use circstats
    implicit none
    integer, parameter :: n=30000
    real(dp), allocatable :: x(:)
    real(dp) :: m,v
    real(dp) :: angles(10)
    type(vm_bootstrap_result) :: ci
    integer :: nseed
    integer, allocatable :: seed(:)
    allocate(x(n))
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed=271828
    call random_seed(put=seed)
    call rstable(n,1.3_dp,2.0_dp,0.8_dp,x)
    m=sum(x)/real(n,dp)
    v=sum((x-m)**2)/real(n-1,dp)
    if (abs(m)>0.08_dp) error stop 1
    if (abs(v-2.0_dp*1.3_dp**2)>0.15_dp) error stop 1
    call rstable(100,1.2_dp,1.0_dp,0.5_dp,x(1:100))
    if (any(ieee_is_nan(x(1:100)))) error stop 1
    angles=[0.1_dp,0.2_dp,0.4_dp,0.5_dp,0.7_dp,0.9_dp,1.0_dp,1.1_dp,0.3_dp,0.6_dp]
    ci=vm_bootstrap_ci(angles,reps=200)
    if (size(ci%mu_reps)/=200 .or. size(ci%kappa_reps)/=200) error stop 1
    if (ci%kappa_ci(1)>ci%kappa_ci(2)) error stop 1
    print *, "test_stable_bootstrap: PASS"
end program
