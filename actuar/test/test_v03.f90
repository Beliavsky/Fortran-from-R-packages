module test_v03_callbacks
    use actuar_kinds, only: dp
    implicit none
contains
    function exp_cdf_cb(x,par) result(v)
        real(dp),intent(in)::x,par(:)
        real(dp)::v
        if(x<=0.0_dp)then;v=0.0_dp;else;v=1.0_dp-exp(-par(1)*x);end if
    end function exp_cdf_cb
    function exp_pdf_cb(x,par) result(v)
        real(dp),intent(in)::x,par(:)
        real(dp)::v
        if(x<0.0_dp)then;v=0.0_dp;else;v=par(1)*exp(-par(1)*x);end if
    end function exp_pdf_cb
    function exp_lev_cb(x,par) result(v)
        real(dp),intent(in)::x,par(:)
        real(dp)::v
        if(x<=0.0_dp)then;v=0.0_dp;else;v=(1.0_dp-exp(-par(1)*x))/par(1);end if
    end function exp_lev_cb
    subroutine freq_mix_draw(anc,values)
        real(dp),intent(in)::anc(:,:)
        real(dp),intent(out)::values(:)
        integer::i
        do i=1,size(values);values(i)=real(i,dp);end do
    end subroutine freq_mix_draw
    subroutine freq_final_draw(anc,values)
        real(dp),intent(in)::anc(:,:)
        integer,intent(out)::values(:)
        integer::i
        do i=1,size(values);values(i)=max(0,nint(anc(i,1)));end do
    end subroutine freq_final_draw
    subroutine sev_mix_draw(anc,values)
        real(dp),intent(in)::anc(:,:)
        real(dp),intent(out)::values(:)
        integer::i
        do i=1,size(values);values(i)=10.0_dp*real(i,dp);end do
    end subroutine sev_mix_draw
    subroutine sev_final_draw(anc,values)
        real(dp),intent(in)::anc(:,:)
        real(dp),intent(out)::values(:)
        values=anc(:,1)
    end subroutine sev_final_draw
end module test_v03_callbacks

program test_v03
    use actuar
    use test_v03_callbacks
    implicit none
    type(mde_result_t)::mc,mg,ml
    type(coverage_spec_t)::cs
    type(hache_barycenter_result_t)::hb
    type(hierarc_exact_result_t)::he
    type(bstraub_result_t)::bs
    type(hierarchy_node_counts_t)::nodes(2)
    type(hierarchy_real_model_t)::fmix(1),smix(1),sfinal
    type(hierarchy_count_model_t)::ffinal
    type(hierarchical_portfolio_t)::port
    real(dp)::bounds(5),counts(4),par(1),y,d,r,u,expected,raw(20)
    real(dp)::ratios(3,5),weights(3,5),design(5,2),ints(3),slopes(3)
    real(dp)::crat(4,4),cwt(4,4)
    integer::class1(4,1),i,j

    bounds=[0.0_dp,-log(0.75_dp)/2.0_dp,-log(0.5_dp)/2.0_dp,-log(0.25_dp)/2.0_dp,10.0_dp]
    counts=25.0_dp;par=[1.2_dp]
    do i=1,20;raw(i)=-log(1.0_dp-real(i,dp)/21.0_dp)/2.0_dp;end do
    mc=mde_cvm(raw,exp_cdf_cb,[1.0_dp],lower=[0.1_dp],upper=[5.0_dp])
    call check_true(mc%converged .and. mc%estimate(1)>1.5_dp .and. mc%estimate(1)<2.5_dp, &
                    'MDE individual CvM')
    mg=mde_grouped_chisq(bounds,counts,exp_cdf_cb,par,lower=[0.1_dp],upper=[5.0_dp])
    call check_true(mg%converged,'MDE chi-square convergence')
    call check_close(mg%estimate(1),2.0_dp,2.0e-4_dp,'MDE chi-square exponential')
    mc=mde_grouped_cvm(bounds,counts,exp_cdf_cb,par,lower=[0.1_dp],upper=[5.0_dp])
    call check_true(mc%converged,'MDE grouped CvM convergence')
    call check_close(mc%estimate(1),2.0_dp,2.0e-4_dp,'MDE grouped CvM exponential')
    ml=mde_grouped_las(bounds,counts,exp_lev_cb,[1.5_dp],lower=[0.1_dp],upper=[5.0_dp])
    call check_true(ml%converged .and. ml%estimate(1)>0.0_dp,'MDE LAS convergence')

    cs%deductible=1.0_dp;cs%limit=3.0_dp;cs%coinsurance=0.8_dp;cs%inflation=0.1_dp;cs%per_loss=.true.
    par=[2.0_dp];r=1.1_dp;d=1.0_dp/r;u=3.0_dp/r
    call check_close(coverage_cdf(0.0_dp,par,exp_cdf_cb,cs),exp_cdf_cb(d,par),1.0e-13_dp,'coverage zero mass')
    call check_close(coverage_pdf(0.0_dp,par,exp_pdf_cb,exp_cdf_cb,cs),exp_cdf_cb(d,par),1.0e-13_dp, &
                     'coverage pdf zero mass')
    y=0.6_dp;expected=exp_cdf_cb((y/0.8_dp+1.0_dp)/r,par)
    call check_close(coverage_cdf(y,par,exp_cdf_cb,cs),expected,1.0e-13_dp,'coverage interior CDF')
    expected=exp_pdf_cb((y/0.8_dp+1.0_dp)/r,par)/(0.8_dp*r)
    call check_close(coverage_pdf(y,par,exp_pdf_cb,exp_cdf_cb,cs),expected,1.0e-13_dp,'coverage interior density')
    y=cs%maximum_payment();call check_close(coverage_pdf(y,par,exp_pdf_cb,exp_cdf_cb,cs), &
                     1.0_dp-exp_cdf_cb(u,par),1.0e-13_dp,'coverage limit mass')
    cs%franchise=.true.;cs%per_loss=.true.
    call check_close(coverage_cdf(0.4_dp,par,exp_cdf_cb,cs),exp_cdf_cb(d,par),1.0e-13_dp, &
                     'franchise per-loss zero branch')
    cs%per_loss=.false.
    call check_close(coverage_cdf(0.4_dp,par,exp_cdf_cb,cs),0.0_dp,1.0e-14_dp, &
                     'franchise per-payment gap')
    cs%franchise=.false.;y=0.6_dp
    expected=(exp_cdf_cb((y/0.8_dp+1.0_dp)/r,par)-exp_cdf_cb(d,par))/(1.0_dp-exp_cdf_cb(d,par))
    call check_close(coverage_cdf(y,par,exp_cdf_cb,cs),expected,1.0e-13_dp,'coverage per-payment CDF')

    do j=1,5;design(j,1)=1.0_dp;design(j,2)=real(j-1,dp);end do
    ints=[1.0_dp,2.0_dp,3.0_dp];slopes=[0.5_dp,-0.2_dp,1.0_dp];weights=1.0_dp
    do i=1,3;do j=1,5;ratios(i,j)=ints(i)+slopes(i)*design(j,2);end do;end do
    hb=hachemeister_barycenter_fit(ratios,weights,design)
    call check_true(hb%converged,'Hachemeister barycenter convergence')
    do i=1,3
        call check_close(hb%individual(1,i),ints(i),2.0e-10_dp,'Hachemeister bary intercept')
        call check_close(hb%individual(2,i),slopes(i),2.0e-10_dp,'Hachemeister bary slope')
        call check_true(maxval(abs(hb%adjusted(:,i)-hb%individual(:,i)))<1.0e-8_dp,'Hachemeister zero-noise credibility')
    end do

    crat(1,:)=[0.95_dp,1.00_dp,1.05_dp,1.00_dp]
    crat(2,:)=[1.95_dp,2.00_dp,2.05_dp,2.00_dp]
    crat(3,:)=[2.95_dp,3.00_dp,3.05_dp,3.00_dp]
    crat(4,:)=[3.95_dp,4.00_dp,4.05_dp,4.00_dp]
    cwt=1.0_dp;class1(:,1)=[1,2,3,4]
    bs=bstraub_fit(crat,cwt,iterative=.true.)
    he=hierarc_exact_fit(crat,cwt,class1,HIERARC_ITERATIVE)
    call check_true(he%converged,'exact hierarc convergence')
    call check_close(he%variance(2),bs%process_variance,2.0e-12_dp,'hierarc process variance')
    call check_close(he%variance(1),bs%between_variance,2.0e-9_dp,'hierarc one-level between variance')
    call check_close(he%level(0)%mean(1),bs%collective_mean,2.0e-9_dp,'hierarc collective mean')
    call check_true(maxval(abs(he%level(1)%credibility-bs%credibility))<2.0e-9_dp,'hierarc credibility oracle')

    allocate(nodes(1)%count(1),nodes(2)%count(2));nodes(1)%count=[2];nodes(2)%count=[2,1]
    fmix(1)%draw=>freq_mix_draw;ffinal%draw=>freq_final_draw
    smix(1)%draw=>sev_mix_draw;sfinal%draw=>sev_final_draw
    port=rcomphierarc_simulate(nodes,fmix,ffinal,smix,sfinal)
    call check_true(all(port%frequency==[1,1,2]),'hierarchical simulation frequencies')
    call check_true(size(port%claims)==4,'hierarchical simulation claim count')
    call check_true(maxval(abs(port%aggregate-[10.0_dp,10.0_dp,40.0_dp]))<1.0e-13_dp,'hierarchical simulation aggregate')

    print '(a)', 'test_v03: PASS'
contains
    subroutine check_close(got,expected,eps,name)
        real(dp),intent(in)::got,expected,eps
        character(*),intent(in)::name
        if(abs(got-expected)>eps*max(1.0_dp,abs(expected)))then
            print *,'FAIL ',trim(name),got,expected;error stop 1
        end if
    end subroutine check_close
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(*),intent(in)::name
        if(.not.ok)then;print *,'FAIL ',trim(name);error stop 1;end if
    end subroutine check_true
end program test_v03
