program test_v02
    use actuar
    implicit none
    real(dp),parameter::tol=2.0e-8_dp
    real(dp)::x,p,meanv,varv,orders(2),gm(2)
    real(dp)::sev(2),freq(5),ratios(4,4),weights(4,4),design(5,2)
    real(dp)::hratios(5,5),hweights(5,5),ha(5),hb(5),noise(5)
    real(dp),allocatable::epc(:),etc(:,:),epw(:),etw(:,:)
    real(dp)::prob_c(1),tc(1,1),prob_w(1),tw(1,1)
    real(dp)::bounds(4),counts(3)
    integer::class1(4,1),i
    type(aggregate_dist_t)::agg
    type(ruin_result_t)::rr
    type(hierarchical_credibility_result_t)::hc
    type(hachemeister_result_t)::hr

    p=pfpareto(6.0_dp,1.0_dp,2.5_dp,1.4_dp,1.7_dp,3.0_dp)
    x=qfpareto(p,1.0_dp,2.5_dp,1.4_dp,1.7_dp,3.0_dp)
    call check_close(x,6.0_dp,2.0e-9_dp,'Feller-Pareto roundtrip')
    call check_close(mfpareto(1.0_dp,0.0_dp,3.0_dp,2.0_dp,1.5_dp,2.0_dp), &
                     mtrbeta(1.0_dp,3.0_dp,2.0_dp,1.5_dp,2.0_dp),tol,'Feller min=0 moment')
    call check_close(pinvpareto(qinvpareto(0.61_dp,2.2_dp,4.0_dp),2.2_dp,4.0_dp), &
                     0.61_dp,2.0e-11_dp,'inverse Pareto roundtrip')
    call check_close(pinvexp(3.0_dp,2.0_dp),exp(-2.0_dp/3.0_dp),1.0e-14_dp,'inverse exponential CDF')
    call check_close(dinvexp(3.0_dp,2.0_dp),2.0_dp*exp(-2.0_dp/3.0_dp)/9.0_dp,1.0e-14_dp, &
                     'inverse exponential density')
    call check_close(qinvexp(pinvexp(3.0_dp,2.0_dp),2.0_dp),3.0_dp,1.0e-13_dp, &
                     'inverse exponential quantile')
    call check_close(minvexp(0.5_dp,2.0_dp),sqrt(2.0_dp)*gamma(0.5_dp),1.0e-13_dp, &
                     'inverse exponential moment')
    call check_close(pinvtrgamma(qinvtrgamma(0.72_dp,3.0_dp,1.4_dp,2.0_dp),3.0_dp,1.4_dp,2.0_dp), &
                     0.72_dp,2.0e-9_dp,'inverse transformed gamma roundtrip')
    call check_close(levinvweibull(1.0e6_dp,3.5_dp,2.0_dp,1.0_dp), &
                     minvweibull(1.0_dp,3.5_dp,2.0_dp),2.0e-6_dp,'inverse Weibull limited moment')
    call check_close(levtrgamma(1.0e4_dp,2.5_dp,1.2_dp,2.0_dp,1.0_dp), &
                     mtrgamma(1.0_dp,2.5_dp,1.2_dp,2.0_dp),2.0e-7_dp,'trgamma limited moment')
    call check_close(betaint_raw(0.4_dp,2.2_dp,3.1_dp), &
                     gamma(2.2_dp)*gamma(3.1_dp)*reg_beta(0.4_dp,2.2_dp,3.1_dp),1.0e-11_dp,'beta integral')

    sev=[0.0_dp,1.0_dp]
    freq=[0.10_dp,0.20_dp,0.30_dp,0.25_dp,0.15_dp]
    agg=aggregate_exact(sev,freq,100.0_dp)
    call check_true(size(agg%pmf)==size(freq),'exact aggregate size')
    call check_true(maxval(abs(agg%pmf-freq))<1.0e-13_dp,'exact aggregate degenerate severity')
    call check_close(aggregate_normal_cdf(5.0_dp,5.0_dp,4.0_dp),0.5_dp,1.0e-14_dp,'normal aggregate')
    p=aggregate_npower_cdf(8.0_dp,5.0_dp,4.0_dp,0.8_dp)
    call check_close(aggregate_npower_cdf(aggregate_npower_quantile(p,5.0_dp,4.0_dp,0.8_dp), &
                     5.0_dp,4.0_dp,0.8_dp),p,2.0e-10_dp,'normal-power inversion')

    prob_c=[1.0_dp]
    tc(1,1)=-2.0_dp
    prob_w=[1.0_dp]
    tw(1,1)=-0.5_dp
    rr=ruin_phase_type(prob_c,tc,prob_w,tw,1.0_dp)
    call check_true(rr%converged,'ruin Cramer-Lundberg convergence')
    call check_close(rr%probability(0.0_dp),0.25_dp,2.0e-11_dp,'ruin at zero')
    call check_close(rr%probability(1.0_dp),0.25_dp*exp(-1.5_dp),3.0e-10_dp,'ruin exponential exact')
    call make_erlang_phase([2],[2.0_dp],[1.0_dp],epc,etc)
    call make_erlang_phase([2],[1.5_dp],[1.0_dp],epw,etw)
    rr=ruin_phase_type(epc,etc,epw,etw,1.5_dp,tol=1.0e-11_dp,maxit=1000)
    call check_true(rr%converged,'ruin Erlang fixed point')
    call check_true(rr%probability(2.0_dp)<rr%probability(0.0_dp),'ruin decreases with surplus')

    bounds=[0.0_dp,2.0_dp,5.0_dp,9.0_dp]
    counts=[2.0_dp,3.0_dp,5.0_dp]
    meanv=grouped_mean(bounds,counts)
    call check_close(meanv,(2.0_dp*1.0_dp+3.0_dp*3.5_dp+5.0_dp*7.0_dp)/10.0_dp,1.0e-13_dp,'grouped mean')
    call check_close(ogive_eval(5.0_dp,bounds,counts),0.5_dp,1.0e-13_dp,'ogive boundary')
    call check_close(grouped_quantile(0.5_dp,bounds,counts),5.0_dp,1.0e-13_dp,'grouped median')
    orders=[1.0_dp,2.0_dp]
    gm=grouped_moments(bounds,counts,orders)
    call check_true(gm(2)>gm(1)**2,'grouped second moment')
    call check_true(elev_grouped(4.0_dp,bounds,counts)<=4.0_dp,'grouped LEV bound')

    ratios=reshape([1.0_dp,1.1_dp,0.9_dp,1.0_dp, 1.3_dp,1.2_dp,1.4_dp,1.3_dp, &
                    0.8_dp,0.9_dp,0.7_dp,0.8_dp, 1.5_dp,1.4_dp,1.6_dp,1.5_dp],shape(ratios))
    weights=1.0_dp
    class1(:,1)=[1,2,3,4]
    hc=hierarchical_credibility(ratios,weights,class1)
    call check_true(hc%nlevels==1,'hierarchical one level')
    call check_true(all(hc%level(1)%credibility>=0.0_dp .and. hc%level(1)%credibility<=1.0_dp), &
                    'hierarchical credibility range')
    call check_true(size(hc%level(1)%premium)==4,'hierarchical premiums')

    design(:,1)=1.0_dp
    design(:,2)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    ha=[0.90_dp,1.15_dp,0.78_dp,1.28_dp,1.02_dp]
    hb=[0.38_dp,0.53_dp,0.34_dp,0.47_dp,0.59_dp]
    noise=[0.02_dp,-0.015_dp,0.01_dp,-0.005_dp,-0.01_dp]
    hweights=1.0_dp
    do i=1,5
        hratios(i,:)=ha(i)+hb(i)*design(:,2)+cshift(noise,i-1)
    end do
    hr=hachemeister_fit(hratios,hweights,design,maxit=80)
    call check_true(hr%converged,'Hachemeister convergence')
    call check_close(hr%collective(1),sum(ha)/5.0_dp,5.0e-2_dp,'Hachemeister collective intercept')
    call check_close(hr%collective(2),sum(hb)/5.0_dp,5.0e-2_dp,'Hachemeister collective slope')
    call check_true(size(hr%adjusted,2)==5,'Hachemeister adjusted size')
    call check_true(hr%process_variance>=0.0_dp,'Hachemeister process variance')

    call check_close(mbeta_act(1.0_dp,2.0_dp,3.0_dp),0.4_dp,1.0e-13_dp,'beta raw moment')
    call check_close(levbeta_act(2.0_dp,2.0_dp,3.0_dp,1.0_dp),0.4_dp,1.0e-13_dp,'beta limited moment')
    call check_close(mbeta(1.0_dp,2.0_dp,3.0_dp),0.4_dp,1.0e-13_dp,'beta compatibility alias')
    call check_close(levbeta(2.0_dp,2.0_dp,3.0_dp,1.0_dp),0.4_dp,1.0e-13_dp,'beta limited compatibility alias')
    call check_close(munif(1.0_dp,2.0_dp,6.0_dp),4.0_dp,1.0e-13_dp,'uniform compatibility alias')
    call check_close(levunif(6.0_dp,2.0_dp,6.0_dp,1.0_dp),4.0_dp,1.0e-13_dp,'uniform limited compatibility alias')
    call check_close(mchisq_act(1.0_dp,4.0_dp,0.0_dp),4.0_dp,1.0e-13_dp,'chi-square mean')
    call check_close(levchisq_act(1.0e4_dp,4.0_dp,0.0_dp,1.0_dp),4.0_dp,1.0e-10_dp, &
                     'chi-square limited mean')
    call check_close(mgfchisq_act(0.1_dp,4.0_dp,0.0_dp),1.5625_dp,1.0e-13_dp,'chi-square mgf')
    call check_close(levinvexp_act(1.0e8_dp,2.0_dp,0.5_dp),minvexp(0.5_dp,2.0_dp),2.0e-4_dp, &
                     'inverse exponential limited moment')
    call check_close(levinvgauss_act(1.0e6_dp,2.0_dp,0.3_dp),2.0_dp,2.0e-5_dp, &
                     'inverse Gaussian limited mean')

    print '(a)','test_v02: PASS'
contains
    pure function cshift(base,k) result(v)
        real(dp),intent(in)::base(:)
        integer,intent(in)::k
        real(dp)::v(size(base))
        integer::j,n
        n=size(base)
        do j=1,n
            v(j)=base(mod(j+k-1,n)+1)
        end do
    end function cshift
    subroutine check_close(got,expected,eps,name)
        real(dp),intent(in)::got,expected,eps
        character(*),intent(in)::name
        if(abs(got-expected)>eps*max(1.0_dp,abs(expected))) then
            print *,'FAIL ',trim(name),got,expected
            error stop 1
        end if
    end subroutine check_close
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(*),intent(in)::name
        if(.not.ok) then
            print *,'FAIL ',trim(name)
            error stop 1
        end if
    end subroutine check_true
end program test_v02
