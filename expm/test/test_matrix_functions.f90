program test_matrix_functions
    use expm_module
    implicit none
    real(dp) :: a(3,3),truth(3,3),err
    real(dp), allocatable :: x(:,:),xp(:,:),xt(:,:),xa(:,:),xr(:,:),xw(:,:),p4(:,:)
    a = reshape([4.0_dp,1.0_dp,1.0_dp, 2.0_dp,4.0_dp,1.0_dp, 0.0_dp,1.0_dp,4.0_dp],[3,3])
    truth=reshape([147.866622446369_dp,127.781085523181_dp,127.781085523182_dp, &
        183.765138646367_dp,183.765138646366_dp,163.679601723179_dp, &
        71.797032399996_dp,91.8825693231832_dp,111.968106246371_dp],[3,3])
    x=expm(a); xp=expm_pade(a,8); xt=expm_taylor(a,18); xa=expm_almohy09(a,13); xr=expm_rbs(a,6); xw=expm_ward77(a,8)
    err=maxval(abs(x-truth))/maxval(abs(truth)); if(err>2.0e-12_dp) error stop "Higham expm regression failed"
    if(maxval(abs(xp-truth))/maxval(abs(truth))>2.0e-12_dp) error stop "Pade regression failed"
    if(maxval(abs(xt-truth))/maxval(abs(truth))>2.0e-12_dp) error stop "Taylor regression failed"
    if(maxval(abs(xa-truth))/maxval(abs(truth))>2.0e-12_dp) error stop "AlMohy09 regression failed"
    if(maxval(abs(xr-truth))/maxval(abs(truth))>2.0e-12_dp) error stop "RBS regression failed"
    if(maxval(abs(xw-truth))/maxval(abs(truth))>2.0e-12_dp) error stop "Ward regression failed"
    p4=matrix_power(a,4); if(maxval(abs(p4-matmul(matmul(a,a),matmul(a,a))))>1.0e-12_dp) error stop "matrix_power failed"
    print *, "test_matrix_functions: PASS"
end program test_matrix_functions
