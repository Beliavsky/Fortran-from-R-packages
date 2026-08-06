program test_toeplitz
    use ltsa, only : dp, ltsa_error, exact_likelihood_result, toeplitz_matrix, trench_inverse, &
                     toeplitz_inverse_update, trench_mean, trench_loglikelihood, dl_loglikelihood, exact_loglikelihood
    implicit none
    type(ltsa_error) :: error
    type(exact_likelihood_result) :: exact
    real(dp), allocatable :: r(:), g(:,:), gi(:,:), expected(:,:), updated(:,:), direct(:,:)
    real(dp) :: phi, d, mean_est, ll1, ll2
    real(dp) :: z(4)
    integer :: i

    phi = 0.6_dp
    r = [(phi**(i-1), i=1,4)]
    g = toeplitz_matrix(r)
    call trench_inverse(g,gi,error)
    call assert_true(error%ok(),'TrenchInverse failed')
    allocate(expected(4,4),source=0.0_dp)
    d = 1.0_dp/(1.0_dp-phi*phi)
    expected(1,1)=d; expected(4,4)=d
    expected(2,2)=(1.0_dp+phi*phi)*d; expected(3,3)=expected(2,2)
    expected(1,2)=-phi*d; expected(2,1)=expected(1,2)
    expected(2,3)=-phi*d; expected(3,2)=expected(2,3)
    expected(3,4)=-phi*d; expected(4,3)=expected(3,4)
    call assert_matrix_close(gi,expected,2.0e-11_dp,'Trench inverse')
    call assert_matrix_close(matmul(g,gi),identity(4),2.0e-11_dp,'inverse product')

    call toeplitz_inverse_update(gi,r,phi**4,updated,error)
    call assert_true(error%ok(),'ToeplitzInverseUpdate failed')
    call trench_inverse(toeplitz_matrix([(phi**(i-1),i=1,5)]),direct,error)
    call assert_true(error%ok(),'direct updated inverse failed')
    call assert_matrix_close(updated,direct,2.0e-10_dp,'updated inverse')

    z = 3.5_dp
    mean_est = trench_mean(r,z,error)
    call assert_true(error%ok(),'TrenchMean failed')
    call assert_close(mean_est,3.5_dp,1.0e-12_dp,'TrenchMean constant')

    z = [0.2_dp,-0.3_dp,0.8_dp,0.1_dp]
    ll1 = trench_loglikelihood(r,z,error)
    call assert_true(error%ok(),'TrenchLoglikelihood failed')
    ll2 = dl_loglikelihood(r,z,error)
    call assert_true(error%ok(),'DLLoglikelihood comparison failed')
    call assert_close(ll1,ll2,2.0e-10_dp,'Trench and DL loglikelihood')
    exact = exact_loglikelihood(r,z)
    call assert_true(exact%error%ok(),'exactLoglikelihood failed')
    call assert_true(exact%sigma_sq > 0.0_dp,'exact variance must be positive')

    print '(a)', 'test_toeplitz: PASS'
contains
    function identity(n) result(a)
        integer,intent(in)::n
        real(dp)::a(n,n)
        integer::j
        a=0.0_dp
        do j=1,n
            a(j,j)=1.0_dp
        end do
    end function identity
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
    subroutine assert_matrix_close(actual,expected_value,tolerance,message)
        real(dp),intent(in)::actual(:,:),expected_value(:,:),tolerance
        character(len=*),intent(in)::message
        if(maxval(abs(actual-expected_value))>tolerance) then
            print *,message,maxval(abs(actual-expected_value))
            error stop
        end if
    end subroutine assert_matrix_close
end program test_toeplitz
