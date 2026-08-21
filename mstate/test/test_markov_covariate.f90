program test_markov_covariate
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use mstate
    implicit none
    type(transition_map)::tr
    type(msdata_type)::ms
    type(markov_test_result)::mt
    real(dp)::times(18,3),nanv,grid(3)
    real(dp),allocatable::cov(:,:)
    integer::status(18,3),info,i

    call trans_illdeath(tr)
    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    times(:,1)=nanv;status(:,1)=0
    do i=1,18
        times(i,2)=0.5_dp+0.35_dp*real(i,dp)
        times(i,3)=times(i,2)+0.8_dp+0.1_dp*real(mod(i,5),dp)
        status(i,2)=merge(1,0,mod(i,4)/=0)
        status(i,3)=merge(1,0,mod(i,6)/=0)
    end do
    call msprep(times,status,tr,ms,info=info)
    call check(info==0,'msprep')
    allocate(cov(ms%n,1))
    do i=1,ms%n
        cov(i,1)=sin(0.7_dp*real(ms%id(i),dp))
    end do
    grid=[1.0_dp,3.0_dp,5.0_dp]
    call markov_test(ms,tr,3,grid,20,mt,covariates=cov,dist='poisson',seed=4321,info=info)
    call check(info==0,'markov covariate info')
    call check(size(mt%beta)==1.and.size(mt%beta_vcov,1)==1,'cox fit stored')
    call check(ieee_is_finite(mt%beta(1)),'finite beta')
    call check(all(ieee_is_finite(mt%zbar)),'finite normalized traces')
    print '(a)','test_markov_covariate: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program test_markov_covariate
