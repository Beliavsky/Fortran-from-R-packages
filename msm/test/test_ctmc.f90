program test_ctmc
    use msm, only : dp, make_generator, transition_matrix, transition_derivatives, observation_kernel, &
        death_transition_density, ctmc_minus2loglik, obs_exact, obs_death, obs_panel, &
        sojourn_times, next_state_probabilities, expected_total_time, expected_first_passage, &
        eventual_passage_probability, passage_probability_matrix, expected_visits, state_prevalence, &
        ctmc_aggregate_minus2loglik, ctmc_censored_minus2loglik, ctmc_gradient
    implicit none
    real(dp) :: off(2,2), q(2,2), a,b,t,e,expected,nll,eps
    real(dp), allocatable :: p(:,:),dpq(:,:,:),dq(:,:,:),kp(:,:),tau(:),r(:,:),et(:,:),ef(:),hit(:),pp(:,:),ev(:,:),prev(:,:)
    integer :: states(3),ot(3),fr(2),to(2),cnt(2)
    logical :: allowed(2,3)
    real(dp), allocatable :: q3(:,:,:),grad(:),qrep(:,:),drep(:,:,:),fdrep(:,:,:),pplus(:,:),pminus(:,:)
    real(dp) :: times(3)
    logical :: target(2),target3(3)
    a=0.4_dp; b=0.2_dp; t=1.7_dp; off=0.0_dp; off(1,2)=a; off(2,1)=b; q=make_generator(off)
    p=transition_matrix(q,t); e=exp(-(a+b)*t)
    call check(abs(p(1,1)-(b/(a+b)+a/(a+b)*e))<2e-13_dp,"two-state p11")
    call check(abs(p(1,2)-(a/(a+b)*(1.0_dp-e)))<2e-13_dp,"two-state p12")
    allocate(dq(2,2,1)); dq=0.0_dp; dq(1,1,1)=-1.0_dp; dq(1,2,1)=1.0_dp
    call transition_derivatives(q,dq,t,dpq,obs_panel)
    eps=1e-6_dp; off=0.0_dp; off(1,2)=a+eps; off(2,1)=b
    kp=transition_matrix(make_generator(off),t); off(1,2)=a-eps
    kp=(kp-transition_matrix(make_generator(off),t))/(2.0_dp*eps)
    call check(maxval(abs(kp-dpq(:,:,1)))<2e-8_dp,"Frechet dP/da")
    off=0.0_dp; off(1,2)=a; q=make_generator(off); kp=observation_kernel(q,t,obs_exact)
    call check(abs(kp(1,2)-a*exp(-a*t))<1e-14_dp,"exact transition density")
    p=transition_matrix(q,t)
    call check(abs(death_transition_density(1,2,p,q)-a*exp(-a*t))<2e-14_dp,"death density")
    states=[1,1,2]; times=[0.0_dp,1.0_dp,2.0_dp]; ot=[obs_panel,obs_panel,obs_panel]
    nll=ctmc_minus2loglik(states,times,q,ot)
    expected=-2.0_dp*(log(exp(-a))+log(1.0_dp-exp(-a)))
    call check(abs(nll-expected)<2e-13_dp,"panel likelihood")
    ot=[obs_panel,obs_panel,obs_death]
    nll=ctmc_minus2loglik(states,times,q,ot)
    expected=-2.0_dp*(log(exp(-a))+log(a*exp(-a)))
    call check(abs(nll-expected)<2e-13_dp,"death likelihood")
    tau=sojourn_times(q); call check(abs(tau(1)-1.0_dp/a)<1e-14_dp,"sojourn")
    r=next_state_probabilities(q); call check(abs(r(1,2)-1.0_dp)<1e-14_dp,"next state")
    et=expected_total_time(q,2.0_dp)
    call check(abs(et(1,1)-(1.0_dp-exp(-2.0_dp*a))/a)<2e-13_dp,"finite total time transient")
    call check(abs(et(1,2)-(2.0_dp-(1.0_dp-exp(-2.0_dp*a))/a))<2e-13_dp,"finite total time absorbing")
    target=[.false.,.true.]; ef=expected_first_passage(q,target); hit=eventual_passage_probability(q,target)
    call check(abs(ef(1)-1.0_dp/a)<2e-13_dp,"expected first passage")
    call check(abs(hit(1)-1.0_dp)<2e-13_dp,"eventual passage")
    pp=passage_probability_matrix(q,2.0_dp); call check(abs(pp(1,2)-(1.0_dp-exp(-2.0_dp*a)))<2e-13_dp,"finite passage")
    ev=expected_visits(q,2.0_dp); call check(abs(ev(1,2)-(1.0_dp-exp(-2.0_dp*a)))<2e-13_dp,"expected visits")
    prev=state_prevalence([1.0_dp,0.0_dp],q,[0.0_dp,2.0_dp]); call check(abs(prev(2,2)-pp(1,2))<2e-13_dp,"prevalence")

    allocate(q3(2,2,1)); q3(:,:,1)=q
    fr=[1,1]; to=[1,2]; cnt=[2,3]
    nll=ctmc_aggregate_minus2loglik(fr,to,[1.0_dp,1.0_dp],cnt,q3,[obs_panel,obs_panel])
    expected=-2.0_dp*(2.0_dp*log(exp(-a))+3.0_dp*log(1.0_dp-exp(-a)))
    call check(abs(nll-expected)<3e-13_dp,"aggregate likelihood")
    allowed=.true.; allowed(:,1)=[.true.,.false.]; allowed(:,3)=[.false.,.true.]
    nll=ctmc_censored_minus2loglik(times,q3,allowed)
    expected=-2.0_dp*log(1.0_dp-exp(-2.0_dp*a))
    call check(abs(nll-expected)<3e-13_dp,"censored likelihood")
    states=[1,1,2]; ot=[obs_panel,obs_panel,obs_panel]; grad=ctmc_gradient(states,times,q,dq,ot)
    eps=1e-6_dp; off=0.0_dp; off(1,2)=a+eps; nll=ctmc_minus2loglik(states,times,make_generator(off),ot)
    off(1,2)=a-eps; expected=ctmc_minus2loglik(states,times,make_generator(off),ot)
    call check(abs(grad(1)-(nll-expected)/(2.0_dp*eps))<3e-7_dp,"CTMC gradient")

    ! Repeated nonzero eigenvalue: progressive 1->2->3 chain with equal rates.
    allocate(qrep(3,3),drep(3,3,1)); qrep=0.0_dp
    qrep(1,2)=0.3_dp; qrep(2,3)=0.3_dp
    qrep(1,1)=-0.3_dp; qrep(2,2)=-0.3_dp
    drep=0.0_dp; drep(1,1,1)=-1.0_dp; drep(1,2,1)=1.0_dp
    call transition_derivatives(qrep,drep,1.3_dp,fdrep,obs_panel)
    pplus=qrep; pplus(1,1)=-(0.3_dp+eps); pplus(1,2)=0.3_dp+eps
    pminus=qrep; pminus(1,1)=-(0.3_dp-eps); pminus(1,2)=0.3_dp-eps
    pplus=(transition_matrix(pplus,1.3_dp)-transition_matrix(pminus,1.3_dp))/(2.0_dp*eps)
    call check(maxval(abs(pplus-fdrep(:,:,1)))<3e-8_dp,"Frechet repeated-eigen derivative")

    ! Hitting a target need not be certain when another absorbing class exists.
    qrep=0.0_dp; qrep(1,2)=0.25_dp; qrep(1,3)=0.25_dp; qrep(1,1)=-0.5_dp
    target3=[.false.,.true.,.false.]
    hit=eventual_passage_probability(qrep,target3)
    ef=expected_first_passage(qrep,target3)
    call check(abs(hit(1)-0.5_dp)<3e-13_dp,"partial eventual passage")
    call check(ef(1)>0.5_dp*huge(1.0_dp),"infinite first passage when hit probability below one")
    call check(abs(hit(3))<3e-13_dp,"unreachable target probability")
    print '(a)', 'test_ctmc: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; write(*,'(a)') 'FAIL: '//msg; error stop 1; end if
    end subroutine check
end program test_ctmc
