program test_eigen_methods
    use expm_module
    implicit none
    real(dp) :: a(2,2)
    complex(dp), allocatable :: xe(:,:),le(:,:)
    real(dp), allocatable :: xh(:,:)
    integer :: info
    a=reshape([1.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2])
    xe=expm_eigen(a,info); if(info/=0) error stop "expm_eigen failed"
    if(abs(xe(1,1)-exp(1.0_dp))>1.0e-13_dp .or. abs(xe(2,2)-exp(2.0_dp))>1.0e-13_dp) error stop "expm_eigen wrong"
    le=logm_eigen(a,info); if(info/=0) error stop "logm_eigen failed"
    if(abs(le(1,1))>1.0e-13_dp .or. abs(le(2,2)-log(2.0_dp))>1.0e-13_dp) error stop "logm_eigen wrong"
    xh=expm_hybrid_eigen_ward(a); if(maxval(abs(xh-real(xe,dp)))>1.0e-13_dp) error stop "hybrid failed"
    print *, "test_eigen_methods: PASS"
end program test_eigen_methods
