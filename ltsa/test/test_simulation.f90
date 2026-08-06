program test_simulation
    use iso_fortran_env, only : int64
    use ltsa, only : dp, ltsa_error, set_ltsa_seed, sim_glp, dh_condition, dh_simulate
    implicit none
    type(ltsa_error) :: error
    real(dp), allocatable :: z(:), r(:), z_source(:)
    real(dp) :: psi(3), innovations(6), expected(4), mean_value, variance, lag1
    logical :: valid
    integer :: i, n

    psi = [1.0_dp,0.5_dp,-0.25_dp]
    innovations = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
    expected = [3.75_dp,5.0_dp,6.25_dp,7.5_dp]
    call sim_glp(psi,innovations,z,error)
    call assert_true(error%ok(),'SimGLP failed')
    call assert_vector_close(z,expected,1.0e-12_dp,'SimGLP values')

    n = 4096
    allocate(r(n))
    r = [(0.7_dp**(i-1)/0.51_dp,i=1,n)]
    valid = dh_condition(n,r)
    call assert_true(valid,'AR(1) should satisfy Davies-Harte condition')
    call set_ltsa_seed(1234567_int64)
    call dh_simulate(n,r,z,error)
    call assert_true(error%ok(),'DHSimulate failed')
    mean_value = sum(z)/real(n,dp)
    variance = sum((z-mean_value)**2)/real(n,dp)
    lag1 = dot_product(z(1:n-1)-mean_value,z(2:n)-mean_value)/real(n-1,dp)
    call assert_close(mean_value,0.0_dp,0.12_dp,'DH mean')
    call assert_close(variance,r(1),0.22_dp,'DH variance')
    call assert_close(lag1/r(1),0.7_dp,0.09_dp,'DH lag one correlation')

    call set_ltsa_seed(1234567_int64)
    call dh_simulate(64,r,z_source,error,source_compatible=.true.)
    call assert_true(error%ok(),'source-compatible DH failed')
    call assert_true(all(z_source==z_source),'source-compatible DH finite')

    deallocate(r)
    allocate(r(3))
    r = [1.0_dp,2.0_dp,2.0_dp]
    valid = dh_condition(8,r)
    call assert_true(.not.valid,'invalid covariance should fail DH condition')

    print '(a)', 'test_simulation: PASS'
contains
    subroutine assert_true(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) error stop message
    end subroutine assert_true
    subroutine assert_close(actual,expected_value,tolerance,message)
        real(dp),intent(in)::actual,expected_value,tolerance
        character(len=*),intent(in)::message
        if(abs(actual-expected_value)>tolerance) then
            print *,message,actual,expected_value
            error stop
        end if
    end subroutine assert_close
    subroutine assert_vector_close(actual,expected_value,tolerance,message)
        real(dp),intent(in)::actual(:),expected_value(:),tolerance
        character(len=*),intent(in)::message
        if(maxval(abs(actual-expected_value))>tolerance) then
            print *,message,maxval(abs(actual-expected_value))
            error stop
        end if
    end subroutine assert_vector_close
end program test_simulation
