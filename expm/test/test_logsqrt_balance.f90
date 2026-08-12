program test_logsqrt_balance
    use expm_module
    implicit none
    real(dp) :: a(3,3),u(2,2)
    complex(dp), allocatable :: s(:,:),l(:,:),ea(:,:)
    type(balance_real_result) :: br
    a=reshape([2.0_dp,0.3_dp,-0.2_dp, 0.1_dp,1.5_dp,0.4_dp, -0.1_dp,0.2_dp,1.2_dp],[3,3])
    s=sqrtm(a)
    if(maxval(abs(matmul(s,s)-cmplx(a,0.0_dp,dp)))>2.0e-8_dp) error stop "sqrtm failed"
    ea=cmplx(expm(a,balancing=.false.),0.0_dp,dp); l=logm(real(ea,dp))
    if(maxval(abs(l-cmplx(a,0.0_dp,dp)))>2.0e-8_dp) error stop "logm(expm(A)) failed"
    u=reshape([1.0_dp,1.0_dp,0.0_dp,1.0_dp],[2,2]); l=logm(u)
    if(abs(l(2,1)-1.0_dp)>1.0e-12_dp .or. maxval(abs([l(1,1),l(1,2),l(2,2)]))>1.0e-12_dp) &
        error stop "logm Jordan block failed"
    br=balance_real(a,'B'); if(br%info/=0) error stop "balance failed"
    print *, "test_logsqrt_balance: PASS"
end program test_logsqrt_balance
