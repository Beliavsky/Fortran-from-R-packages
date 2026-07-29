! SPDX-License-Identifier: GPL-3.0-or-later
program test_data_diagnostics
  use acdm
  implicit none
  integer :: failures
  failures=0
  call test_duration_construction
  call test_diurnal
  call test_diagnostics
  call test_specification
  if(failures>0) error stop 'test_data_diagnostics failed'
  print '(a)','test_data_diagnostics: PASS'
contains
  subroutine test_duration_construction
    integer,parameter::n=10
    integer::yy(n),mm(n),dd(n)
    real(dp)::tt(n),pr(n),vol(n)
    type(duration_result)::r
    yy=2025;mm=7;dd=[1,1,1,1,1,2,2,2,2,2]
    tt=[36010._dp,36010._dp,36020._dp,36040._dp,70000._dp, &
        36005._dp,36015._dp,36025._dp,36025._dp,36050._dp]
    pr=[100._dp,100.02_dp,100.15_dp,100.18_dp,100.4_dp, &
        101._dp,101.05_dp,101.25_dp,101.27_dp,101.5_dp]
    vol=[40._dp,60._dp,70._dp,40._dp,100._dp,30._dp,80._dp,50._dp,70._dp,100._dp]
    call compute_durations(yy,mm,dd,tt,r,DURATION_TRADE,36000._dp,66600._dp,.true.,pr,vol)
    call assert_true(r%status==ACDM_SUCCESS,'trade duration status')
    call assert_true(size(r%durations)==7,'trade duration count')
    call assert_close(r%durations(1),10._dp,1e-13_dp,'trade first duration')
    call assert_true(r%transaction_count(1)==2,'trade aggregation')
    call compute_durations(yy,mm,dd,tt,r,DURATION_PRICE,36000._dp,66600._dp,.true.,pr,vol,0.1_dp)
    call assert_true(r%status==ACDM_SUCCESS.and.size(r%durations)>=3,'price durations')
    call compute_durations(yy,mm,dd,tt,r,DURATION_VOLUME,36000._dp,66600._dp,.true.,pr,vol,cumulative_volume=100._dp)
    call assert_true(r%status==ACDM_SUCCESS.and.size(r%durations)>=3,'volume durations')
    call assert_true(all(r%durations>0._dp),'positive aggregated durations')
  end subroutine
  subroutine test_diurnal
    integer,parameter::n=240
    real(dp)::time(n),dur(n),nodes(9)
    integer::i,meth,grp(n)
    type(diurnal_result)::r
    nodes=[36000._dp,36600._dp,37200._dp,37800._dp,38400._dp,39000._dp,39600._dp,40200._dp,40800._dp]
    do i=1,n
      time(i)=36000._dp+real(mod(i-1,120),dp)*(4800._dp/119._dp)
      dur(i)=2._dp+0.7_dp*cos(2._dp*acos(-1._dp)*(time(i)-36000._dp)/4800._dp)+0.05_dp*sin(real(i,dp))
      grp(i)=1+mod((i-1)/120,2)
    end do
    do meth=DIURNAL_CUBIC_SPLINE,DIURNAL_FFF
      if(meth<=DIURNAL_SMOOTH_SPLINE)then
        call diurnal_adjust(time,dur,r,meth,nodes=nodes)
      else
        call diurnal_adjust(time,dur,r,meth)
      end if
      call assert_true(r%status==ACDM_SUCCESS,'diurnal method')
      if(r%status==ACDM_SUCCESS)then
        call assert_true(all(r%fitted>0._dp).and.all(r%adjusted>0._dp),'diurnal positive')
        call assert_true(abs(sum(r%adjusted)/real(n,dp)-1._dp)<0.15_dp,'diurnal normalization')
      end if
    end do
    call diurnal_adjust(time,dur,r,DIURNAL_FFF,group=grp,fourier_order=3)
    call assert_true(r%status==ACDM_SUCCESS.and.size(r%grid_fit,2)==2,'grouped FFF')
    do i=1,n
      time(i)=real(i-1,dp)
      dur(i)=2.0_dp+0.2_dp*time(i)/real(n-1,dp)+0.3_dp*cos(2.0_dp*acos(-1.0_dp)*time(i)/real(n-1,dp))-0.1_dp*sin(2.0_dp*acos(-1.0_dp)*time(i)/real(n-1,dp))
    end do
    call diurnal_adjust(time,dur,r,DIURNAL_FFF,fourier_order=1)
    call assert_true(r%status==ACDM_SUCCESS.and.maxval(abs(r%fitted-dur))<1e-9_dp,'exact FFF recovery')
  end subroutine
  subroutine test_diagnostics
    real(dp)::res(5),pit(5),cs(5)
    type(acf_result)::ar
    type(density_result)::dr
    type(qq_result)::qr
    type(summary_result)::sr
    integer::st
    res=[0.2_dp,0.5_dp,1._dp,1.5_dp,2._dp]
    call standardize_residuals(res,DIST_EXPONENTIAL,[real(dp)::],.true.,pit,st)
    call assert_true(st==ACDM_SUCCESS,'PIT status')
    call assert_close(pit(3),1._dp-exp(-1._dp),1e-13_dp,'PIT value')
    call standardize_residuals(res,DIST_EXPONENTIAL,[real(dp)::],.true.,cs,st,cox_snell=.true.)
    call assert_true(maxval(abs(cs-res))<1e-13_dp,'Cox-Snell exponential')
    call acf_acd(res,3,ar)
    call assert_true(ar%status==ACDM_SUCCESS.and.size(ar%acf)==3,'ACF')
    call residual_density_acd(res,DIST_EXPONENTIAL,[real(dp)::],.true.,dr,ngrid=51)
    call assert_true(dr%status==ACDM_SUCCESS.and.size(dr%x)==51,'density')
    call qqplot_acd(res,DIST_EXPONENTIAL,[real(dp)::],.true.,qr,npoints=5)
    call assert_true(qr%status==ACDM_SUCCESS.and.all(qr%theoretical>0._dp),'QQ')
    call summarize_durations(res,sr,window=3)
    call assert_true(sr%status==ACDM_SUCCESS.and.size(sr%rolling_mean)==3,'summary')
  end subroutine
  subroutine test_specification
    integer,parameter::n=350
    real(dp)::x(n),mu(n),rr(n),par(3),time(n)
    type(acd_order)::o
    type(rng_state)::rng
    type(lm_test_result)::lr
    integer::st,i
    o=acd_order(1,0,1);par=[0.2_dp,0.15_dp,0.7_dp]
    call seed_rng(rng,4281)
    call simulate_acd(n,MODEL_ACD,o,par,DIST_EXPONENTIAL,[real(dp)::],x,st,rng,burn=200)
    call filter_acd(x,MODEL_ACD,o,par,mu,rr,st)
    time=[(real(i,dp),i=1,n)]
    call test_rm_acd(x,mu,par,1,1,2,.true.,lr)
    call assert_test(lr,'remaining ACD robust')
    call test_rm_acd(x,mu,par,1,1,2,.false.,lr)
    call assert_test(lr,'remaining ACD nonrobust')
    call test_st_acd(x,mu,par,1,1,2,.true.,lr)
    call assert_test(lr,'STACD')
    call test_tv_acd(x,mu,par,1,1,time,2,.true.,lr)
    call assert_test(lr,'TVACD')
  end subroutine
  subroutine assert_test(r,label)
    type(lm_test_result),intent(in)::r;character(*),intent(in)::label
    call assert_true(r%status==ACDM_SUCCESS,trim(label)//' status')
    call assert_true(r%statistic>=0._dp.and.r%p_value>=0._dp.and.r%p_value<=1._dp,trim(label)//' values')
  end subroutine
  subroutine assert_close(a,b,tol,label)
    real(dp),intent(in)::a,b,tol;character(*),intent(in)::label
    if(abs(a-b)>tol*max(1._dp,abs(b)))then;failures=failures+1;print *,'FAIL ',trim(label),a,b;end if
  end subroutine
  subroutine assert_true(ok,label)
    logical,intent(in)::ok;character(*),intent(in)::label
    if(.not.ok)then;failures=failures+1;print *,'FAIL ',trim(label);end if
  end subroutine
end program
