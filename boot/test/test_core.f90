program test_core
    use boot_kinds, only : dp
    use boot_core
    use boot_statistics, only : mean_dp
    implicit none
    real(dp)::data(5,1)
    type(bootstrap_result)::res
    data(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    call seed_fixed()
    call bootstrap_weighted(data,stat,200,'balanced',res)
    if(abs(res%t0-3.0_dp)>1.0e-12_dp)error stop 1
    if(abs(mean_dp(res%t)-3.0_dp)>1.0e-12_dp)error stop 2
    if(any(sum(res%frequencies,dim=2)/=5))error stop 3
    print '(a)', 'test_core: PASS'
contains
    function stat(x,w) result(v)
        real(dp),intent(in)::x(:,:),w(:)
        real(dp)::v
        v=sum(x(:,1)*w)/sum(w)
    end function stat
    subroutine seed_fixed()
        integer::m,i
        integer,allocatable::s(:)
        call random_seed(size=m)
        allocate(s(m))
        s=[(9001+17*i,i=1,m)]
        call random_seed(put=s)
    end subroutine seed_fixed
end program test_core
