program test_sparse
use fields, only: dp,csr_matrix,wendland_covariance,wendland_covariance_sparse,krig_fit,sparse_krig_fit, &
                  krig_fit_covariance,sparse_krig_fit_covariance
implicit none
integer,parameter::n=10
real(dp)::x(n,1),y(n),t(n,1)
real(dp),allocatable::kd(:,:)
type(csr_matrix)::ks
type(krig_fit)::fd
type(sparse_krig_fit)::fs
integer::i
do i=1,n;x(i,1)=real(i-1,dp)/real(n-1,dp);y(i)=cos(3*x(i,1));end do
t=1.0_dp;kd=wendland_covariance(x,x,0.8_dp,2);ks=wendland_covariance_sparse(x,x,0.8_dp,2)
fd=krig_fit_covariance(y,kd,0.1_dp,t);fs=sparse_krig_fit_covariance(y,ks,0.1_dp,t,pivot='none')
call check(fd%info==0 .and. fs%info==0,'fits')
call check(maxval(abs(fd%fitted-fs%fitted))<2e-8_dp,'sparse/dense fitted')
call check(maxval(abs(fd%c-fs%c))<2e-8_dp,'sparse/dense coefficients')
print *,'test_sparse: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(len=*),intent(in)::msg
if(.not.ok)then;print *,'FAIL: ',msg;error stop 1;end if
end subroutine
end program
