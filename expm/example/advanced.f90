program advanced
    use expm_module
    implicit none
    real(dp) :: a(3,3),e(3,3),v(3),c1,cf
    real(dp), allocatable :: ea(:,:),l(:,:)
    complex(dp), allocatable :: sa(:,:),la(:,:)
    type(exp_action_result) :: av
    a=reshape([0.2_dp,-0.1_dp,0.3_dp, 0.4_dp,0.1_dp,-0.2_dp, &
               0.0_dp,0.5_dp,-0.3_dp],[3,3])
    e=0.0_dp; e(1,2)=1.0_dp; v=[1.0_dp,2.0_dp,-1.0_dp]
    call expm_frechet_sps(a,e,ea,l)
    av=exp_at_v(a,v,t=0.5_dp)
    call expm_cond_exact(a,c1,cf)
    sa=sqrtm(matmul(transpose(a),a)+reshape([1.0_dp,0.0_dp,0.0_dp, &
        0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp],[3,3]))
    la=logm(ea)
    print '(a,es12.4)', 'cond_1(exp,A) = ', c1
    print '(a,es12.4)', 'cond_F(exp,A) = ', cf
    print '(a,3f12.6)', 'exp(0.5*A)v = ', av%value
    print '(a,es12.4)', 'sqrt residual = ', maxval(abs(matmul(sa,sa)-cmplx( &
        matmul(transpose(a),a)+reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp, &
        0.0_dp,0.0_dp,0.0_dp,1.0_dp],[3,3]),0.0_dp,dp)))
    print '(a,es12.4)', 'log(exp(A))-A = ', maxval(abs(la-cmplx(a,0.0_dp,dp)))
end program advanced
