program test_distribution
    use circstats
    implicit none
    real(dp), parameter :: tol=3.0e-7_dp
    call check(i0e(2.3_dp),0.28369298571049134_dp,2.0e-7_dp,"i0e")
    call check(i1e(2.3_dp),0.21032300512056257_dp,2.0e-7_dp,"i1e")
    call check(a1(2.3_dp),0.7413754153766677_dp,3.0e-7_dp,"a1")
    call check(ip(2,2.3_dp),1.0054316636559144_dp,3.0e-12_dp,"ip")
    call check(dvm(1.2_dp,0.4_dp,2.3_dp),0.2792664393490445_dp,tol,"dvm")
    call check(pvm(1.2_dp,0.4_dp,2.3_dp),0.5728606953990438_dp,2.0e-7_dp,"pvm")
    call check(dcard(1.1_dp,0.3_dp,0.2_dp),0.20350866976305126_dp,1.0e-13_dp,"dcard")
    call check(dtri(1.1_dp,0.25_dp),0.18857971351657637_dp,1.0e-13_dp,"dtri")
    call check(dwrpcauchy(1.1_dp,0.3_dp,0.6_dp),0.19440554389639456_dp,1.0e-13_dp,"dwrpcauchy")
    call check(dwrpnorm(1.1_dp,0.3_dp,0.7_dp,tol=1.0e-13_dp),0.30160573549651737_dp,2.0e-12_dp,"dwrpnorm")
    call check(pvm(twopi,0.7_dp,4.0_dp),1.0_dp,1.0e-12_dp,"pvm full circle")
    print *, "test_distribution: PASS"
contains
    subroutine check(got,want,eps,label)
        real(dp),intent(in)::got,want,eps
        character(*),intent(in)::label
        if (abs(got-want)>eps) then
            print *, trim(label),got,want
            error stop 1
        end if
    end subroutine
end program
