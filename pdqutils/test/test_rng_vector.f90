program test_rng_vector
    use pdqutils, only : dp, rapx_cf, qapx_cf_vec, dapx_edgeworth_vec, gca_basis_from_name, gca_gamma
    implicit none
    real(dp), allocatable :: x(:),q(:),d(:)
    real(dp) :: meanv,varv
    integer :: n,fails
    fails=0
    call set_seed(24681357)
    n=50000
    allocate(x(n))
    call rapx_cf([2.0_dp,9.0_dp,0.0_dp,0.0_dp,0.0_dp],x)
    meanv=sum(x)/real(n,dp)
    varv=sum((x-meanv)**2)/real(n-1,dp)
    call check(abs(meanv-2.0_dp)<0.05_dp,'rng mean')
    call check(abs(sqrt(varv)-3.0_dp)<0.05_dp,'rng sd')
    q=qapx_cf_vec([0.1_dp,0.5_dp,0.9_dp],[0.0_dp,1.0_dp,0.0_dp,0.0_dp])
    call check(maxval(abs(q-[-1.2815515655446004_dp,0.0_dp,1.2815515655446004_dp]))<2e-14_dp,'vector q')
    d=dapx_edgeworth_vec([-1.0_dp,0.0_dp,1.0_dp],[0.0_dp,1.0_dp,0.0_dp,0.0_dp])
    call check(maxval(abs(d-[0.24197072451914337_dp,0.3989422804014327_dp,0.24197072451914337_dp]))<2e-14_dp,'vector d')
    call check(gca_basis_from_name('Gamma')==gca_gamma,'basis name')
    if(fails==0) then
        print '(a)','test_rng_vector: PASS'
    else
        print '(a,i0)','test_rng_vector: FAIL ',fails
        error stop 1
    end if
contains
    subroutine set_seed(seed)
        integer,intent(in)::seed
        integer :: nseed,i
        integer, allocatable :: put(:)
        call random_seed(size=nseed)
        allocate(put(nseed))
        do i=1,nseed
            put(i)=mod(seed+104729*i,huge(1)-1)
            if(put(i)<=0) put(i)=i
        end do
        call random_seed(put=put)
    end subroutine set_seed
    subroutine check(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok) then
            print '(a,a)','FAIL: ',trim(name); fails=fails+1
        end if
    end subroutine check
end program test_rng_vector
