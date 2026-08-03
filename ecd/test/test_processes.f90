! SPDX-License-Identifier: Artistic-2.0
program test_processes
  use ecd_api
  implicit none
  type(cumulants4) :: k
  type(sld_model) :: d
  type(rng_state) :: rng
  real(dp) :: x(50000),m,v
  integer :: st

  call close(stdlap_pdf(0.25_dp,0.7_dp,2.0_dp,0.3_dp,0.1_dp), &
    0.5906335493931782_dp,3.0e-11_dp,'stdlap pdf')
  call close(stdlap_cdf(0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp),0.5_dp,2.0e-9_dp,'stdlap symmetry')
  call close(stdlap_quantile(0.5_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,st),0.0_dp,2.0e-8_dp,'stdlap median')
  call check(st==ecd_ok,'stdlap quantile status')
  k=stdlap_cumulants(1.0_dp,1.0_dp,0.0_dp,0.0_dp)
  call close(k%variance,1.0_dp,1.0e-14_dp,'stdlap variance')
  call close(k%kurtosis,6.0_dp,1.0e-14_dp,'stdlap kurtosis')

  call close(stable_count_pdf(7.0_dp,0.5_dp,1.0_dp,2.0_dp), &
    0.05769987105204573_dp,2.0e-12_dp,'stable count pdf')
  call close(stable_count_cdf(7.0_dp,0.5_dp,1.0_dp,2.0_dp), &
    0.31772966966378746_dp,2.0e-12_dp,'stable count cdf')
  call close(stable_count_quantile(0.8_dp,0.5_dp,1.0_dp,2.0_dp,status=st), &
    19.5665107043498_dp,2.0e-9_dp,'stable count quantile')
  call check(st==ecd_ok,'stable count quantile status')
  k=stable_count_cumulants(theta=1.0_dp,status=st)
  call close(k%mean,6.0_dp,1.0e-14_dp,'stable count mean')
  call close(k%variance,24.0_dp,1.0e-14_dp,'stable count variance')
  call close(k%skewness,sqrt(8.0_dp/3.0_dp),1.0e-14_dp,'stable count skewness')
  call close(k%kurtosis,7.0_dp,1.0e-14_dp,'stable count kurtosis')

  d=sld_new(t=1.0_dp,nu0=0.0_dp,theta=1.0_dp,convo=1.0_dp,beta_a=0.0_dp,mu=0.0_dp)
  k=sld_cumulants(d,st)
  call check(st==ecd_ok,'sld cumulant status')
  call close(k%mean,0.0_dp,1.0e-14_dp,'sld mean')
  call close(k%variance,60.0_dp,1.0e-12_dp,'sld variance')
  call close(real(sld_cf(d,0.0_dp),dp),1.0_dp,2.0e-6_dp,'sld cf zero')
  call close(qsl_variance_analytic(),60.0_dp,1.0e-13_dp,'qsl variance')

  call rng_seed(rng,20260729_i8)
  call stdlap_random(rng,x,1.0_dp,1.0_dp,0.0_dp,0.0_dp)
  m=sum(x)/real(size(x),dp); v=sum((x-m)**2)/real(size(x)-1,dp)
  call check(abs(m)<0.025_dp,'stdlap random mean')
  call check(abs(v-1.0_dp)<0.05_dp,'stdlap random variance')

  call close(levy_dlambda(2.0_dp,4.0_dp),exp(-sqrt(2.0_dp))/4.0_dp,1.0e-14_dp,'levy dlambda')
  print '(a)', 'test_processes: PASS'
contains
  subroutine close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol*max(1.0_dp,abs(y)))then
      write(*,*)trim(msg),x,y; error stop 1
    end if
  end subroutine close
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,*)trim(msg);error stop 1;end if
  end subroutine check
end program test_processes
