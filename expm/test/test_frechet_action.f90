program test_frechet_action
    use expm_module
    implicit none
    real(dp) :: a(3,3),e(3,3),v(3),cond1,condf
    real(dp), allocatable :: x1(:,:),l1(:,:),x2(:,:),l2(:,:),direct(:)
    type(exp_action_result) :: ar
    a=reshape([0.2_dp,-0.1_dp,0.3_dp, 0.4_dp,0.1_dp,-0.2_dp, 0.0_dp,0.5_dp,-0.3_dp],[3,3])
    e=reshape([0.1_dp,0.2_dp,-0.1_dp, 0.0_dp,-0.3_dp,0.4_dp, 0.2_dp,0.0_dp,0.1_dp],[3,3])
    call expm_frechet_sps(a,e,x1,l1); call expm_frechet_block(a,e,x2,l2)
    if(maxval(abs(x1-x2))>2.0e-13_dp) error stop "Frechet expm mismatch"
    if(maxval(abs(l1-l2))>5.0e-13_dp) error stop "Frechet derivative mismatch"
    v=[1.0_dp,-2.0_dp,0.5_dp]; ar=exp_at_v(a,v,t=0.7_dp,tol=1.0e-10_dp)
    direct=matmul(expm(0.7_dp*a,balancing=.false.),v)
    if(maxval(abs(ar%value-direct))>2.0e-9_dp) error stop "exp_at_v mismatch"
    call expm_cond_exact(a,cond1,condf)
    if(.not.(cond1>0.0_dp .and. condf>0.0_dp)) error stop "condition number failed"
    print *, "test_frechet_action: PASS"
end program test_frechet_action
