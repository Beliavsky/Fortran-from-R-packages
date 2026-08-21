program test_factor
  use survey, only : dp, factor_result_t, factor_ml_cov, varimax_rotate
  implicit none
  real(dp) :: cov(4,4), l(4), psi(4), load2(6,2), before(6,6), after(6,6)
  type(factor_result_t) :: fit
  integer :: i,j
  l=[0.8_dp,0.7_dp,0.6_dp,0.5_dp];psi=1.0_dp-l*l
  do j=1,4
    do i=1,4
      cov(i,j)=l(i)*l(j)
    end do
    cov(j,j)=cov(j,j)+psi(j)
  end do
  call factor_ml_cov(cov,1,fit,500.0_dp,.false.,maxfun=12000)
  call assert_true(fit%converged,'factor convergence')
  call assert_close(maxval(abs(abs(fit%loadings(:,1))-l)),0.0_dp,3e-3_dp,'factor loading recovery')
  call assert_close(maxval(abs(fit%uniqueness-psi)),0.0_dp,4e-3_dp,'factor uniqueness recovery')
  call assert_true(fit%criterion<1e-6_dp,'factor exact-model criterion')
  call assert_true(fit%p_value>0.99_dp,'factor exact-model p value')

  load2=reshape([0.8_dp,0.7_dp,0.6_dp,0.1_dp,0.1_dp,0.2_dp, &
                 0.1_dp,0.2_dp,0.1_dp,0.8_dp,0.7_dp,0.6_dp],[6,2],order=[2,1])
  before=matmul(load2,transpose(load2));call varimax_rotate(load2);after=matmul(load2,transpose(load2))
  call assert_close(maxval(abs(before-after)),0.0_dp,2e-12_dp,'varimax preserves reproduced covariance')
  print '(a)','test_factor: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
  end subroutine assert_true
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(*),intent(in)::msg
    if(abs(a-b)>tol)then;write(*,'(a,2es24.14)')'FAIL: '//trim(msg)//' ',a,b;error stop 1;end if
  end subroutine assert_close
end program test_factor
