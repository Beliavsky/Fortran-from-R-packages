program test_random_regression
    use expm_module
    implicit none
    integer :: n,i,j
    real(dp), allocatable :: a(:,:),x(:,:),xm(:,:),id(:,:),e(:,:),xp(:,:),lfr(:,:),xb(:,:),lbr(:,:),spd(:,:)
    complex(dp), allocatable :: s(:,:),lg(:,:)
    real(dp) :: err
    do n=2,8
        allocate(a(n,n),id(n,n),e(n,n)); id=0.0_dp
        do i=1,n
            id(i,i)=1.0_dp
            do j=1,n
                a(i,j)=0.12_dp*sin(real(17*i+11*j+n,dp))
                e(i,j)=0.08_dp*cos(real(13*i-7*j+n,dp))
            end do
        end do
        x=expm(a,balancing=.false.); xm=expm(-a,balancing=.false.)
        err=maxval(abs(matmul(x,xm)-id)); if(err>3.0e-12_dp) error stop "random expm inverse identity failed"
        call expm_frechet_sps(a,e,xp,lfr); call expm_frechet_block(a,e,xb,lbr)
        if(maxval(abs(lfr-lbr))>3.0e-12_dp) error stop "random Frechet mismatch"
        lg=logm(x); if(maxval(abs(lg-cmplx(a,0.0_dp,dp)))>2.0e-10_dp) error stop "random log(expm) mismatch"
        spd=matmul(transpose(a),a)+id; s=sqrtm(spd)
        if(maxval(abs(matmul(s,s)-cmplx(spd,0.0_dp,dp)))>2.0e-10_dp) error stop "random sqrtm mismatch"
        deallocate(a,id,e)
    end do
    print *, "test_random_regression: PASS"
end program test_random_regression
