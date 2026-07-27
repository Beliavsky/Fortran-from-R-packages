! SPDX-License-Identifier: GPL-2.0-or-later
program test_screening
  use peerperformance, only: dp, peer_control, screening_result, rolling_result, &
       alpha_screening, sharpe_screening, modified_sharpe_screening, &
       target_peer_performance, roll_screening, exposure_heterogeneity
  implicit none
  integer, parameter :: n=80, p=4
  real(dp) :: x(n,p), f(n,2)
  type(peer_control) :: control, boot_control
  type(screening_result) :: alpha, sh, ms, target, direct
  type(rolling_result) :: rolling
  real(dp), allocatable :: equal(:), hetero(:)
  integer :: status

  call make_data(x,f)
  control%has_lambda=.true.; control%lambda=0.5_dp; control%n_boot=39
  control%min_obs=20; control%min_obs_pi=1; control%screen_beta=.true.
  call alpha_screening(x,control,alpha,factors=f)
  call assert_true(alpha%status==0,'alpha screening')
  call assert_true(alpha%ncoef==3 .and. alpha%n_focal==p,'alpha shape')
  call assert_symmetric(alpha%pvalue,'alpha p-value symmetry')
  call assert_ratios(alpha,'alpha ratios')
  call exposure_heterogeneity(alpha,equal,hetero,status)
  call assert_true(status==0 .and. all(abs(equal+hetero-1.0_dp)<1.0e-13_dp), &
       'exposure heterogeneity')

  control%screen_beta=.false.
  call sharpe_screening(x,control,sh)
  call assert_true(sh%status==0,'Sharpe screening')
  call assert_symmetric(sh%pvalue,'Sharpe symmetry')
  call assert_ratios(sh,'Sharpe ratios')
  call modified_sharpe_screening(x,0.95_dp,control,ms,.false.)
  call assert_true(ms%status==0,'modified Sharpe screening')
  call assert_ratios(ms,'modified Sharpe ratios')

  call target_peer_performance(x,[1,3],'sharpe',control,target)
  call sharpe_screening(x(:,[1,3]),control,direct,peers=x)
  call assert_true(target%status==0 .and. direct%status==0,'target screening')
  call assert_true(all(abs(target%pizero-direct%pizero)<1.0e-14_dp),'target equality')

  call roll_screening(x,'alpha',40,20,control,rolling,factors=f)
  call assert_true(rolling%status==0 .and. rolling%nwindow==3,'rolling screening')
  call assert_true(all(abs(rolling%heterogeneity+rolling%pizero-1.0_dp)<1.0e-13_dp), &
       'rolling heterogeneity')

  boot_control=control; boot_control%test_type=2; boot_control%n_boot=19
  boot_control%block_length=4; boot_control%p_boot=2
  call sharpe_screening(x(:,1:3),boot_control,sh)
  call assert_true(sh%status==0,'bootstrap screening')
  call assert_ratios(sh,'bootstrap ratios')
  print '(a)', 'test_screening: PASS'
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
  subroutine assert_symmetric(a,label)
    real(dp), intent(in) :: a(:,:,:)
    character(len=*), intent(in) :: label
    integer :: i,j
    do i=1,size(a,2)
      do j=1,size(a,3)
        if (i==j) cycle
        if (abs(a(1,i,j)-a(1,j,i))>1.0e-13_dp) then
          print '(a)', 'failed: '//trim(label); error stop 1
        end if
      end do
    end do
  end subroutine assert_symmetric
  subroutine assert_ratios(r,label)
    type(screening_result), intent(in) :: r
    character(len=*), intent(in) :: label
    if (any(abs(r%pizero+r%pipos+r%pineg-1.0_dp)>2.0e-12_dp)) then
      print '(a)', 'failed: '//trim(label); error stop 1
    end if
  end subroutine assert_ratios
  subroutine assert_true(value,label)
    logical, intent(in) :: value
    character(len=*), intent(in) :: label
    if (.not.value) then
      print '(a)', 'failed: '//trim(label); error stop 1
    end if
  end subroutine assert_true
end program test_screening
