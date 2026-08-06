program test_pca_selection
  use tvmvp, only : dp, local_pca_result, factor_selection_result, local_pca_all, &
                    determine_factors, silverman
  implicit none
  integer, parameter :: n=32,p=5
  real(dp) :: x(n,p),f(n),load(p),h
  integer :: i,j
  type(local_pca_result) :: local
  type(factor_selection_result) :: sel
  do i=1,n
    f(i)=sin(0.23_dp*real(i,dp))+0.4_dp*cos(0.07_dp*real(i,dp))
  end do
  load=[1.0_dp,0.8_dp,-0.5_dp,1.3_dp,-0.9_dp]
  do j=1,p
    do i=1,n
      x(i,j)=f(i)*load(j)+0.015_dp*sin(real(i*(j+2),dp))
    end do
  end do
  h=silverman(x)
  call local_pca_all(x,h,1,local)
  call check(.not.local%error%failed(),'local PCA status')
  call check(all(shape(local%f_hat)==[n,1]),'f_hat shape')
  call check(all(shape(local%loadings)==[p,1,n]),'loading shape')
  call check(minval(local%weights)>=-1.0e-14_dp,'nonnegative Epanechnikov weights')
  call determine_factors(x,2,sel,h)
  call check(.not.sel%error%failed(),'factor selection status')
  call check(sel%optimal_m>=1 .and. sel%optimal_m<=2,'selected factor range')
  call check(all(sel%ic_values<huge(1.0_dp)),'finite information criterion')
  call check(sel%residual_variances(2)<=sel%residual_variances(1)+1.0e-10_dp,'residual variance monotonic')
  print *, 'test_pca_selection: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',msg
      error stop 1
    end if
  end subroutine check
end program test_pca_selection
