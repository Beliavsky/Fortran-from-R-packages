program test_resampling
    use boot_resampling
    implicit none
    integer,parameter::n=6,r=12
    integer::strata(n),idx(r,n),freq(r,n),j
    strata=[1,1,1,2,2,2]
    call seed_fixed()
    call balanced_array(n,r,strata,idx)
    call frequency_array(idx,freq)
    if(any(sum(freq,dim=2)/=n))error stop 1
    do j=1,n
        if(sum(freq(:,j))/=r)error stop 2
    end do
    call permutation_array(n,r,strata,idx)
    do j=1,r
        if(sum(idx(j,1:3))/=6 .or. product(idx(j,1:3))/=6)error stop 3
        if(sum(idx(j,4:6))/=15 .or. product(idx(j,4:6))/=120)error stop 4
    end do
    print '(a)', 'test_resampling: PASS'
contains
    subroutine seed_fixed()
        integer::m,i
        integer,allocatable::s(:)
        call random_seed(size=m)
        allocate(s(m))
        s=[(7919+37*i,i=1,m)]
        call random_seed(put=s)
    end subroutine seed_fixed
end program test_resampling
