! SPDX-License-Identifier: GPL-3.0-or-later
program test_skewt
  use nvmix
  implicit none
  type(integration_control) :: ctrl
  type(sample_result) :: draw
  real(dp) :: scale(1,1),v1,v2,mean
  ctrl%samples=1024; scale(1,1)=1.4_dp
  v1=dskewt([0.3_dp],[0.0_dp],9.0_dp,[0.0_dp],scale,ctrl)
  v2=dstudent_mv([0.3_dp],9.0_dp,[0.0_dp],scale)
  call assert_close(v1,v2,4.0e-4_dp,'zero-skew density')
  call assert_close(pskewt1d(0.7_dp,0.0_dp,9.0_dp,0.0_dp,1.4_dp,ctrl),&
    student_cdf(0.7_dp/sqrt(1.4_dp),9.0_dp),2.0e-13_dp,'zero-skew cdf')
  draw=rskewt(100000,[0.25_dp],10.0_dp,[0.0_dp],reshape([1.0_dp],[1,1]),777_i8)
  mean=sum(draw%x(:,1))/real(size(draw%x,1),dp)
  call assert_close(mean,0.25_dp*10.0_dp/8.0_dp,0.025_dp,'skew-t mean')
  call assert_close(pskewt1d(qskewt1d(0.8_dp,0.2_dp,8.0_dp,0.0_dp,1.0_dp,ctrl),&
    0.2_dp,8.0_dp,0.0_dp,1.0_dp,ctrl),0.8_dp,3.0e-4_dp,'skew-t inversion')
  print '(a)','test_skewt: PASS'
contains
  subroutine assert_close(a,b,tol,label)
    real(dp), intent(in) :: a,b,tol
    character(*), intent(in) :: label
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      write(*,'(a,3es24.15)')trim(label)//' mismatch: ',a,b,abs(a-b); error stop 1
    end if
  end subroutine
end program
