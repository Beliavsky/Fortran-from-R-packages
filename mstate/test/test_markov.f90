program test_markov
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use mstate
    implicit none
    type(transition_map)::tr
    type(msdata_type)::ms
    type(markov_test_result)::mt
    real(dp)::times(12,3),nanv,grid(3)
    integer::status(12,3),info,i

    call trans_illdeath(tr)
    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    times(:,1)=nanv
    status(:,1)=0
    do i=1,12
        times(i,2)=1.0_dp+0.6_dp*real(i,dp)
        times(i,3)=times(i,2)+1.0_dp+0.15_dp*real(mod(i,4),dp)
        status(i,2)=merge(1,0,mod(i,3)/=0)
        status(i,3)=1
    end do
    call msprep(times,status,tr,ms,info=info)
    call check(info==0,'msprep')
    grid=[1.0_dp,3.0_dp,5.0_dp]
    call markov_test(ms,tr,3,grid,30,mt,dist='normal',seed=1234,info=info)
    call check(info==0,'markov info')
    call check(mt%transition==3.and.mt%b==30,'markov metadata')
    call check(size(mt%qualset)>=2,'qualifying states')
    call check(all(shape(mt%zbar)==[3,size(mt%qualset)]),'trace shape')
    call check(all(ieee_is_finite(mt%obs_chisq_trace)),'finite chi trace')
    call check(all(mt%p_stat_wb>=0.0_dp.and.mt%p_stat_wb<=1.0_dp),'state p values')
    call check(mt%p_ch_stat_wb>=0.0_dp.and.mt%p_ch_stat_wb<=1.0_dp,'overall p value')
    print '(a)','test_markov: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program test_markov
