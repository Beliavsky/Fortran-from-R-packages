program test_complex_and_roots
    use expm_module
    implicit none
    complex(dp) :: a(2,2),truth(2,2)
    real(dp) :: r(2,2),j3(3,3)
    complex(dp), allocatable :: x(:,:),s(:,:),l(:,:)
    a=cmplx(0.0_dp,0.0_dp,dp)
    a(1,1)=cmplx(0.2_dp,0.7_dp,dp); a(2,2)=cmplx(-0.3_dp,-0.2_dp,dp)
    truth=cmplx(0.0_dp,0.0_dp,dp); truth(1,1)=exp(a(1,1)); truth(2,2)=exp(a(2,2))
    x=expm(a)
    if(maxval(abs(x-truth))>2.0e-13_dp) error stop "complex expm failed"
    r=0.0_dp; r(1,1)=-4.0_dp; r(2,2)=9.0_dp; s=sqrtm(r)
    if(maxval(abs(matmul(s,s)-cmplx(r,0.0_dp,dp)))>2.0e-12_dp) error stop "complex sqrt branch failed"
    j3=reshape([1.0_dp,1.0_dp,0.0_dp, 0.0_dp,1.0_dp,1.0_dp, 0.0_dp,0.0_dp,1.0_dp],[3,3])
    l=logm(j3)
    if(maxval(abs(expm(real(l,dp),balancing=.false.)-j3))>1.0e-10_dp) &
        error stop "Jordan log roundtrip failed"
    print *, "test_complex_and_roots: PASS"
end program test_complex_and_roots
