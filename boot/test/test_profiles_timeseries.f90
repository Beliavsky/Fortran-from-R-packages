program test_profiles_timeseries
    use boot_kinds, only : dp
    use boot_profiles
    use boot_timeseries
    implicit none
    real(dp)::u(3),ll,lam,eef,theta(7),like(7),lo,hi
    integer::info,idx(3,10),r,j
    u=[-1.0_dp,0.0_dp,1.0_dp]
    call empirical_loglikelihood(u,ll,lam,info)
    if(info/=0 .or. abs(lam)>1.0e-10_dp .or. abs(ll)>1.0e-10_dp)error stop 1
    call eef_loglikelihood(u,ll,eef,lam,info)
    if(abs(lam)>1.0e-10_dp .or. abs(eef)>1.0e-10_dp)error stop 2
    theta=[-2.0_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,2.0_dp]
    like=-theta*theta
    call likelihood_ci(theta,like,-1.0_dp,lo,hi)
    if(abs(lo+1.0_dp)>1.0e-10_dp .or. abs(hi-1.0_dp)>1.0e-10_dp)error stop 3
    call seed_fixed()
    call fixed_block_indices(7,10,3,3,.true.,idx)
    if(any(idx<1).or.any(idx>7))error stop 4
    do r=1,3
        do j=2,3
            if(idx(r,j)/=1+mod(idx(r,j-1),7))error stop 5
        end do
    end do
    print '(a)', 'test_profiles_timeseries: PASS'
contains
    subroutine seed_fixed()
        integer::m,i
        integer,allocatable::s(:)
        call random_seed(size=m)
        allocate(s(m))
        s=[(1234+53*i,i=1,m)]
        call random_seed(put=s)
    end subroutine seed_fixed
end program test_profiles_timeseries
