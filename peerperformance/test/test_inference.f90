! SPDX-License-Identifier: GPL-2.0-or-later
program test_inference
  use peerperformance, only: dp, test_result, alpha_testing, &
       sharpe_testing_asymptotic, modified_sharpe_testing_asymptotic, &
       sharpe_testing_bootstrap, modified_sharpe_testing_bootstrap, &
       sharpe_block_size, modified_sharpe_block_size
  implicit none
  integer, parameter :: n=80, p=4
  real(dp) :: x(n,p), f(n,2)
  type(test_result) :: result, repeat
  integer :: block

  call make_data(x,f)
  call alpha_testing(x(:,1),x(:,2),result,factors=f,hac=.false.,screen_beta=.true.,min_obs=20)
  call assert_true(result%status==0,'alpha test')
  call assert_true(size(result%difference)==3,'alpha test coefficient count')
  call assert_close(result%difference(1),result%estimate(1,1)-result%estimate(1,2),1.0e-12_dp,'alpha difference')
  call assert_probability(result%pvalue,'alpha p-values')

  call alpha_testing(x(:,1),x(:,2),result,factors=f,hac=.true.,screen_beta=.false.,min_obs=20)
  call assert_true(result%status==0,'HAC alpha test')
  call assert_probability(result%pvalue,'HAC alpha p-value')

  call sharpe_testing_asymptotic(x(:,1),x(:,2),result,hac=.false.,ttype=1,min_obs=20)
  call assert_true(result%status==0,'Sharpe asymptotic')
  call assert_close(result%difference(1),result%estimate(1,1)-result%estimate(1,2),1.0e-12_dp,'Sharpe difference')
  call assert_probability(result%pvalue,'Sharpe p-value')

  call modified_sharpe_testing_asymptotic(x(:,1),x(:,2),0.95_dp,result,.false.,.true.,2,20)
  call assert_true(result%status==0,'modified Sharpe asymptotic')
  call assert_probability(result%pvalue,'modified Sharpe p-value')

  call sharpe_testing_bootstrap(x(:,1),x(:,2),result,n_boot=59,block_length=4, &
       ttype=2,p_boot=2,seed=73,min_obs=20)
  call sharpe_testing_bootstrap(x(:,1),x(:,2),repeat,n_boot=59,block_length=4, &
       ttype=2,p_boot=2,seed=73,min_obs=20)
  call assert_true(result%status==0 .and. repeat%status==0,'Sharpe bootstrap')
  call assert_close(result%pvalue(1),repeat%pvalue(1),0.0_dp,'bootstrap reproducibility')
  call assert_probability(result%pvalue,'bootstrap p-value')

  call modified_sharpe_testing_bootstrap(x(:,1),x(:,2),0.95_dp,result,.false., &
       n_boot=39,block_length=3,ttype=1,p_boot=1,seed=91,min_obs=20)
  call assert_true(result%status==0,'modified Sharpe bootstrap')
  call assert_probability(result%pvalue,'modified bootstrap p-value')

  block=sharpe_block_size(x(:,1),x(:,2),candidates=[1,3],alpha=0.10_dp, &
       m_boot=9,k_sim=3,average_block=3,start_length=8,ttype=2,seed=17)
  call assert_true(block==1 .or. block==3,'Sharpe block selector')
  block=modified_sharpe_block_size(x(:,1),x(:,2),0.95_dp,candidates=[1,3], &
       alpha=0.10_dp,m_boot=9,k_sim=3,average_block=3,start_length=8, &
       ttype=2,seed=23,na_negative=.false.)
  call assert_true(block==1 .or. block==3,'modified Sharpe block selector')
  print '(a)', 'test_inference: PASS'
contains
  subroutine make_data(a,fac)
    real(dp), intent(out) :: a(:,:), fac(:,:)
    integer :: i
    real(dp) :: t
    do i=1,size(a,1)
      t=real(i,dp); fac(i,1)=0.01_dp*sin(0.17_dp*t); fac(i,2)=0.008_dp*cos(0.11_dp*t)
      a(i,1)=0.006_dp+0.8_dp*fac(i,1)+0.2_dp*fac(i,2)+0.020_dp*sin(0.37_dp*t)
      a(i,2)=0.003_dp+0.4_dp*fac(i,1)-0.1_dp*fac(i,2)+0.018_dp*cos(0.29_dp*t)
      a(i,3)=-0.001_dp-0.2_dp*fac(i,1)+0.3_dp*fac(i,2)+0.022_dp*sin(0.23_dp*t+0.4_dp)
      a(i,4)=0.002_dp+0.1_dp*fac(i,1)+0.5_dp*fac(i,2)+0.019_dp*cos(0.31_dp*t+0.2_dp)
    end do
  end subroutine make_data
  subroutine assert_probability(values,label)
    real(dp), intent(in) :: values(:)
    character(len=*), intent(in) :: label
    call assert_true(all(values>=0.0_dp .and. values<=1.0_dp),label)
  end subroutine assert_probability
  subroutine assert_close(actual,expected,tol,label)
    real(dp), intent(in) :: actual,expected,tol
    character(len=*), intent(in) :: label
    if (abs(actual-expected)>tol*(1.0_dp+abs(expected))) then
      print '(a,3(1x,es24.16))', trim(label),actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(value,label)
    logical, intent(in) :: value
    character(len=*), intent(in) :: label
    if (.not.value) then
      print '(a)', 'failed: '//trim(label); error stop 1
    end if
  end subroutine assert_true
end program test_inference
