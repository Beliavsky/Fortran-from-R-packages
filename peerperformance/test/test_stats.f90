! SPDX-License-Identifier: GPL-2.0-or-later
program test_stats
  use peerperformance, only: dp, sharpe, modified_sharpe, alpha_coefficients
  implicit none
  integer, parameter :: n=80, p=4
  real(dp) :: x(n,p), f(n,2), sr(p), msr(p)
  real(dp), allocatable :: beta(:,:)
  integer :: nobs(p), status
  integer, allocatable :: nobs_alpha(:)

  call make_data(x,f)
  call sharpe(x,sr,nobs,status)
  call assert_true(status == 0, 'sharpe status')
  call assert_true(all(nobs == n), 'sharpe counts')
  call assert_vector(sr,[0.460915768788_dp,0.180186993478_dp,-0.066074786910_dp, &
                         0.147767442128_dp],2.0e-8_dp,'sharpe')

  call modified_sharpe(x,0.95_dp,msr,nobs,.false.,status)
  call assert_true(status == 0, 'modified Sharpe status')
  call assert_vector(msr,[0.364756072327_dp,0.124764111401_dp,-0.038488219676_dp, &
                          0.097084722682_dp],2.0e-8_dp,'modified Sharpe')

  call alpha_coefficients(x,beta,nobs_alpha,f,status)
  call assert_true(status == 0, 'alpha coefficients status')
  call assert_true(all(shape(beta) == [3,p]), 'alpha coefficient shape')
  call assert_vector(beta(1,:),[0.006638796002_dp,0.002154761192_dp,-0.001247466633_dp, &
                                0.001813507150_dp],2.0e-8_dp,'alphas')
  call assert_vector(beta(2,:),[0.717769356_dp,-0.037094966_dp,-1.04596528_dp,0.062311161_dp], &
                     2.0e-8_dp,'factor one coefficients')
  print '(a)', 'test_stats: PASS'
contains
  subroutine make_data(a,fac)
    real(dp), intent(out) :: a(:,:), fac(:,:)
    integer :: i
    real(dp) :: t
    do i=1,size(a,1)
      t=real(i,dp)
      fac(i,1)=0.01_dp*sin(0.17_dp*t)
      fac(i,2)=0.008_dp*cos(0.11_dp*t)
      a(i,1)=0.006_dp+0.8_dp*fac(i,1)+0.2_dp*fac(i,2)+0.020_dp*sin(0.37_dp*t)
      a(i,2)=0.003_dp+0.4_dp*fac(i,1)-0.1_dp*fac(i,2)+0.018_dp*cos(0.29_dp*t)
      a(i,3)=-0.001_dp-0.2_dp*fac(i,1)+0.3_dp*fac(i,2)+0.022_dp*sin(0.23_dp*t+0.4_dp)
      a(i,4)=0.002_dp+0.1_dp*fac(i,1)+0.5_dp*fac(i,2)+0.019_dp*cos(0.31_dp*t+0.2_dp)
    end do
  end subroutine make_data
  subroutine assert_vector(actual,expected,tol,label)
    real(dp), intent(in) :: actual(:),expected(:),tol
    character(len=*), intent(in) :: label
    if (size(actual)/=size(expected) .or. any(abs(actual-expected)>tol*(1.0_dp+abs(expected)))) then
      print '(a)', 'mismatch: '//trim(label)
      print '(*(es24.16,1x))', actual
      print '(*(es24.16,1x))', expected
      error stop 1
    end if
  end subroutine assert_vector
  subroutine assert_true(value,label)
    logical, intent(in) :: value
    character(len=*), intent(in) :: label
    if (.not.value) then
      print '(a)', 'failed: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true
end program test_stats
