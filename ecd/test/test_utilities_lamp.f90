! SPDX-License-Identifier: Artistic-2.0
program test_utilities_lamp
  use ecd_api
  implicit none
  real(dp) :: x(6),q
  real(dp), allocatable :: d(:),cur(:),lagged(:),acf(:),acov(:),trimmed(:),lt(:),ut(:)
  integer :: labels(6),st
  type(sample_statistics) :: s
  type(rng_state) :: rng
  type(lamp_model) :: lm
  type(lamp_result) :: lr

  x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
  s=sample_stats(x)
  call close(s%mean,3.5_dp,1.0e-14_dp,'sample mean')
  call close(s%variance,3.5_dp,1.0e-14_dp,'sample variance')
  call difference_series(x,2,d)
  call check(all(abs(d-2.0_dp)<1.0e-14_dp),'difference')
  call lag_series(x,2,cur,lagged)
  call check(all(cur==[3.0_dp,4.0_dp,5.0_dp,6.0_dp]),'current lag series')
  call check(all(lagged==[1.0_dp,2.0_dp,3.0_dp,4.0_dp]),'lagged series')
  call lag_stats(x,2,acf,acov,st)
  call check(st==ecd_ok,'lag stats status')
  call close(acf(0),1.0_dp,1.0e-14_dp,'acf zero')
  q=empirical_quantile(x,0.25_dp)
  call close(q,2.25_dp,1.0e-14_dp,'type 7 quantile')
  call quantilize(x,3,labels,status=st)
  call check(st==ecd_ok .and. all(labels==[1,1,2,2,3,3]),'quantilize')
  call manage_hist_tails(x,0.2_dp,0.8_dp,trimmed,lt,ut)
  call check(size(trimmed)==4 .and. size(lt)==1 .and. size(ut)==1,'tail split')

  call rng_seed(rng,123456_i8)
  lm=lamp_new(lambda=4.0_dp,beta=0.0_dp,random_walk=22,t_infinity=50.0_dp, &
    random_count=2000,n_lower=0.0_dp,n_upper=1000.0_dp,status=st)
  call check(st==ecd_ok,'lamp constructor')
  call lamp_simulate_once(lm,rng,lr,drop=5,keep_tau=0)
  call check(lr%status==ecd_ok,'lamp simulation status')
  call check(size(lr%z)>10,'lamp generated observations')
  call check(all(lr%n>=0.0_dp),'lamp nonnegative counts')
  call check(all(abs(lr%z-lr%b*lr%n)<1.0e-12_dp),'lamp product identity')
  print '(a)', 'test_utilities_lamp: PASS'
contains
  subroutine close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol
    character(len=*),intent(in)::msg
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      write(*,*)trim(msg),a,b; error stop 1
    end if
  end subroutine close
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,*)trim(msg);error stop 1;end if
  end subroutine check
end program test_utilities_lamp
