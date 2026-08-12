program basic
    use expm_module
    implicit none
    real(dp) :: a(2,2)
    real(dp), allocatable :: ea(:,:)
    a=reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
    ea=expm(a)
    print '(a)', 'exp(A) ='
    print '(2f16.10)', ea(1,:)
    print '(2f16.10)', ea(2,:)
end program basic
