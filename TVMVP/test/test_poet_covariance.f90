program test_poet_covariance
  use tvmvp, only : dp, local_pca_result, poet_result, tvmvp_error, local_pca_all, &
                    estimate_residual_cov_poet_local, cholesky_factor, silverman
  implicit none
  integer, parameter :: n=48,p=6
  real(dp) :: x(n,p),f,grid(8)
  real(dp), allocatable :: l(:,:)
  integer :: i,j
  type(local_pca_result) :: local
  type(poet_result) :: poet
  type(tvmvp_error) :: err
  do i=1,n
    f=0.03_dp*sin(0.15_dp*real(i,dp))+0.02_dp*cos(0.04_dp*real(i,dp))
    do j=1,p
      x(i,j)=f*(0.5_dp+0.2_dp*real(j,dp))+0.005_dp*sin(real(i*j,dp))+ &
             0.002_dp*cos(real(i+3*j,dp))
    end do
  end do
  do i=1,size(grid)
    grid(i)=0.01_dp+0.12_dp*real(i-1,dp)
  end do
  call local_pca_all(x,silverman(x),1,local)
  call check(.not.local%error%failed(),'local PCA status')
  call estimate_residual_cov_poet_local(local,x,poet,m0=4,rho_grid=grid)
  call check(.not.poet%error%failed(),'POET status')
  call check(maxval(abs(poet%total_cov-transpose(poet%total_cov)))<1.0e-11_dp,'covariance symmetry')
  call cholesky_factor(poet%total_cov,l,err)
  call check(.not.err%failed(),'positive definite covariance')
  call check(any(abs(poet%best_rho-grid)<1.0e-14_dp),'rho selected from grid')
  call check(maxval(abs([(poet%residual_cov(i,i)-poet%naive_residual_cov(i,i),i=1,p)]))<1.0e-14_dp, &
             'residual diagonal retained')
  print *, 'test_poet_covariance: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',msg
      error stop 1
    end if
  end subroutine check
end program test_poet_covariance
