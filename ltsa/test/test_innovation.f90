program test_innovation
    use iso_fortran_env, only : int64
    use ltsa, only : dp, ltsa_error, innovation_variance_result, set_ltsa_seed, tacvf_arma, dl_simulate, &
                     innovation_variance
    implicit none
    type(ltsa_error) :: error
    type(innovation_variance_result) :: ar_result, k_result, bad_result
    real(dp), allocatable :: r(:), z(:)
    real(dp) :: phi(1)

    phi = 0.65_dp
    call tacvf_arma(phi,[real(dp)::],1199,1.0_dp,r,error)
    call assert_true(error%ok(),'tacvf setup failed')
    call set_ltsa_seed(991827_int64)
    call dl_simulate(1200,r,z,error)
    call assert_true(error%ok(),'simulation setup failed')
    ar_result = innovation_variance(z,'AR',max_order=15)
    call assert_true(ar_result%error%ok(),'AR innovation variance failed')
    call assert_close(ar_result%variance,1.0_dp,0.18_dp,'AR innovation variance')
    call assert_true(ar_result%selected_order>=1 .and. ar_result%selected_order<=15,'selected AR order')

    k_result = innovation_variance(z,'Kolmogoroff',smooth_span=5)
    call assert_true(k_result%error%ok(),'Kolmogoroff innovation variance failed')
    call assert_true(k_result%variance>0.0_dp .and. k_result%variance<5.0_dp,'Kolmogoroff variance range')

    bad_result = innovation_variance(z(1:5),'AR')
    call assert_true(.not.bad_result%error%ok(),'short series should fail')

    print '(a)', 'test_innovation: PASS'
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
end program test_innovation
