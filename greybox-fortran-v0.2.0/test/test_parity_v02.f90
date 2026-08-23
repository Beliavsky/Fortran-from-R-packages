program test_parity_v02
    use greybox, only: dp, alm_model, alm_fit, alm_predict, alm_occurrence_model, &
        alm_dynamic_model, alm_fit_occurrence, alm_fit_arima_errors, rmcb_result, &
        rmcb_test, dsrboot_result, dsr_bootstrap, aid_result, aidcat_result, &
        aid_fit, aid_cat, calm_model, calm_fit, calm_predict, calm_dynamic_model, lm_dynamic_fit
    implicit none
    integer, parameter :: n=60
    real(dp) :: x(n,2), y(n), pred(n), z, err, data(20,4), yi(24), catdata(24,2)
    integer :: i, j, seed_size
    integer, allocatable :: seed(:)
    type(alm_model) :: m, mloss
    type(alm_occurrence_model) :: hm
    type(alm_dynamic_model) :: dm
    type(rmcb_result) :: rt
    type(dsrboot_result) :: boot, iboot
    type(aid_result) :: ar
    type(aidcat_result) :: ac
    type(calm_model) :: cm
    type(calm_dynamic_model) :: cdm
    real(dp) :: pb(n)

    call random_seed(size=seed_size)
    allocate(seed(seed_size)); seed=24681357; call random_seed(put=seed)

    do i=1,n
        z=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
        x(i,1)=1.0_dp; x(i,2)=z
        y(i)=1.0_dp/(1.0_dp+exp(-(0.35_dp+0.9_dp*z))) + 0.025_dp*sin(0.7_dp*real(i,dp))
        y(i)=min(0.98_dp,max(0.02_dp,y(i)))
    end do
    call alm_fit(x,y,'dbeta',m,max_iter=900,tol=2.0e-6_dp)
    call check(m%converged,'beta alm converged')
    call check(allocated(m%scale_beta),'beta scale coefficients allocated')
    call alm_predict(m,x,pred)
    call check(all(pred>0.0_dp .and. pred<1.0_dp),'beta predictions in support')
    call check(sum(abs(pred-y))/real(n,dp)<0.12_dp,'beta fit accuracy')

    call alm_fit(x,y,'dnorm',mloss,loss='LASSO',lambda=0.15_dp,max_iter=500,tol=1.0e-6_dp)
    call check(mloss%converged,'lasso alm')
    call alm_fit(x,y,'dnorm',mloss,loss='RIDGE',lambda=0.15_dp,max_iter=500,tol=1.0e-6_dp)
    call check(mloss%converged,'ridge alm')
    call alm_fit(x,y,'dnorm',mloss,loss='ROLE',trim_fraction=0.1_dp,max_iter=500,tol=1.0e-6_dp)
    call check(mloss%converged,'role alm')
    call alm_fit(x,y,'dnorm',mloss,loss='QUALE',lambda=0.5_dp,max_iter=500,tol=1.0e-6_dp)
    call check(mloss%converged,'quale alm')

    do i=1,n
        z=x(i,2)
        y(i)=2.0_dp+0.7_dp*z+0.08_dp*sin(real(i,dp))
        if(mod(i,4)==0)y(i)=0.0_dp
    end do
    call alm_fit_occurrence(x,y,'dnorm','plogis',hm,max_iter=700,tol=1.0e-6_dp)
    call check(hm%converged,'occurrence model converged')
    call check(all(hm%fitted>=0.0_dp),'occurrence predictions nonnegative')
    call check(hm%loglik>-huge(1.0_dp)/10.0_dp,'occurrence finite likelihood')

    y=0.0_dp
    y(1)=1.0_dp
    do i=2,n
        y(i)=1.0_dp+0.025_dp*real(i,dp)+0.68_dp*(y(i-1)-1.0_dp-0.025_dp*real(i-1,dp)) + &
            0.04_dp*sin(0.4_dp*real(i,dp))
    end do
    call alm_fit_arima_errors(x,y,'dnorm',[1,0,0],dm,max_iter=500,tol=1.0e-6_dp,n_outer=6)
    call check(dm%converged,'dynamic alm converged')
    call check(size(dm%ar)==1,'dynamic ar size')
    call check(abs(dm%ar(1))<1.2_dp,'dynamic ar stable-ish')

    do i=1,20
        do j=1,4
            data(i,j)=real(j-1,dp)+0.03_dp*sin(real(i*j,dp))
        end do
    end do
    call rmcb_test(data,0.95_dp,rt,'tukey')
    call check(rt%selected==1,'rmcb identifies best')
    call check(rt%p_value<0.05_dp,'rmcb significant')
    call check(all([(rt%groups(i,i),i=1,4)]),'rmcb diagonal groups')
    call rmcb_test(data,0.95_dp,rt,'dlaplace')
    call check(rt%selected==1,'rmcb ALM branch identifies best')

    yi=[1.0_dp,2.0_dp,1.5_dp,2.2_dp,3.0_dp,2.8_dp,2.0_dp,1.8_dp,2.3_dp,2.7_dp, &
        3.1_dp,2.9_dp,2.4_dp,2.1_dp,1.7_dp,1.9_dp,2.6_dp,3.2_dp,3.0_dp,2.5_dp, &
        2.2_dp,2.0_dp,2.4_dp,2.8_dp]
    call dsr_bootstrap(yi,25,boot,multiplicative=.false.,parametric=.false.)
    call check(all(shape(boot%boot)==[24,25]),'dsrboot dimensions')
    call check(all(abs(boot%boot)<huge(1.0_dp)/100.0_dp),'dsrboot finite')

    yi=[0.0_dp,2.0_dp,0.0_dp,0.0_dp,3.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,2.0_dp, &
        0.0_dp,0.0_dp,1.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp,0.0_dp,3.0_dp,0.0_dp, &
        1.0_dp,0.0_dp,2.0_dp,0.0_dp]
    call dsr_bootstrap(yi,12,iboot,intermittent=.true.,multiplicative=.false.,parametric=.true.)
    call check(all(shape(iboot%boot)==[24,12]),'intermittent bootstrap dimensions')
    call check(all(iboot%boot>=0.0_dp),'intermittent bootstrap support')

    call aid_fit(yi,ar)
    call check(len_trim(ar%name)>0,'aid category')
    call check(trim(ar%type2)=='intermittent','aid identifies intermittent')
    catdata(:,1)=yi
    catdata(:,2)=3.0_dp+0.1_dp*[(real(i,dp),i=1,24)]
    call aid_cat(catdata,ac)
    call check(sum(ac%types)==2,'aidCat counts series')

    do i=1,n
        z=x(i,2)
        y(i)=1.0_dp/(1.0_dp+exp(-(0.2_dp+0.7_dp*z))) + 0.02_dp*cos(real(i,dp))
        y(i)=min(0.98_dp,max(0.02_dp,y(i)))
    end do
    call calm_fit(x,y,'dbeta',cm,max_predictors=4,max_iter=500,tol=2.0e-6_dp)
    call check(cm%converged,'beta calm converged')
    call check(allocated(cm%scale_beta),'beta calm scale coefficients')
    call calm_predict(cm,x,pb)
    err=sum(abs(pb-y))/real(n,dp)
    call check(err<0.16_dp,'beta calm prediction')


    ! Point-IC dynamic model averaging with LOWESS-style smoothing.
    do i=1,n
        z=x(i,2)
        y(i)=0.8_dp+1.4_dp*z+0.05_dp*sin(0.3_dp*real(i,dp))
    end do
    call lm_dynamic_fit(x,y,'dnorm',cdm,max_predictors=4,max_iter=400,tol=1.0e-6_dp)
    call check(cdm%converged,'lmDynamic core converged')
    call check(all(shape(cdm%weights)==[n,2]),'lmDynamic weight shape')
    call check(maxval(abs(sum(cdm%weights,dim=2)-1.0_dp))<1.0e-10_dp,'lmDynamic weights sum to one')
    call check(sum(abs(cdm%fitted-y))/real(n,dp)<0.12_dp,'lmDynamic fitted accuracy')

    print '(a)', 'test_parity_v02: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok)then
            print '(a)', 'FAIL: '//trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_parity_v02
