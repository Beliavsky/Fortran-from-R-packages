program test_nested
    use boot_kinds, only : dp
    use boot_nested
    implicit none
    real(dp)::d(8,2),w(8),z,p
    integer::i
    do i=1,8
        d(i,1)=real(i,dp)
        d(i,2)=0.7_dp*real(i,dp)+0.2_dp*real(mod(i,3),dp)
    end do
    w=1.0_dp/8.0_dp
    call seed_fixed()
    call nested_correlation(d,w,0.0_dp,49,z,p)
    if(.not.(p>=0.0_dp .and. p<=1.0_dp))error stop 1
    if(z<=0.0_dp)error stop 2
    print '(a)', 'test_nested: PASS'
contains
    subroutine seed_fixed()
        integer::m,j
        integer,allocatable::s(:)
        call random_seed(size=m)
        allocate(s(m))
        s=[(501+j*11,j=1,m)]
        call random_seed(put=s)
    end subroutine seed_fixed
end program test_nested
