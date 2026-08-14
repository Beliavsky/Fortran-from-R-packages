! Statistical network models translated from R/sna models.R and src/likelihood.c.
! Upstream copyright (C) 2004-2024 Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_models
    use sna_kinds, only : dp, sna_nan, is_missing, sna_eps
    use sna_types, only : regression_result, brokerage_result
    use sna_linalg, only : least_squares, logistic_irls, inverse_matrix, dominant_eigenvector, jacobi_eigen_symmetric
    use sna_prep, only : gvectorize, log_mean, log_sub
    use sna_graph, only : dyad_census, triad_census, gden, grecip, gtrans, mutuality, &
        connectedness, hierarchy, lubness, efficiency, neighborhood, degree
    use sna_multivariate, only : centralgraph, graph_correlation
    use sna_centrality, only : betweenness, closeness, centralization_from_scores
    use sna_permutation, only : rmperm
    use sna_random, only : rgraph
    use sna_bn_triad, only : bn_triad_stats, bn_nll_triad
    implicit none
    private

    type, public :: bbnam_result
        real(dp), allocatable :: net(:,:,:)
        real(dp), allocatable :: em(:,:)
        real(dp), allocatable :: ep(:,:)
        integer :: nactors=0, nobservers=0, draws=0
        character(len=12) :: model=''
    end type bbnam_result


    type, public :: bayes_factor_result
        real(dp) :: integrated_loglik(3)=0.0_dp
        real(dp) :: log_bayes_factor(3,3)=0.0_dp
        real(dp) :: log_mc_se(3)=0.0_dp
        integer :: reps=0
    end type bayes_factor_result

    type, public :: bn_result
        real(dp) :: pi=0.0_dp, sigma=0.0_dp, rho=0.0_dp, d=0.0_dp
        real(dp) :: g_square=0.0_dp
        real(dp), allocatable :: triads(:), triads_pred(:)
        real(dp) :: dyads(3)=0.0_dp, dyads_pred(3)=0.0_dp
        logical :: converged=.false.
        integer :: iterations=0
        character(len=16) :: method=''
    end type bn_result

    type, public :: netcancor_result
        real(dp), allocatable :: cor(:), xcoef(:,:), ycoef(:,:)
        real(dp), allocatable :: cor_distribution(:,:)
        real(dp), allocatable :: xcoef_distribution(:,:,:), ycoef_distribution(:,:,:)
        real(dp), allocatable :: cor_p_lower(:), cor_p_upper(:)
        real(dp), allocatable :: xcoef_p_lower(:,:), xcoef_p_upper(:,:)
        real(dp), allocatable :: ycoef_p_lower(:,:), ycoef_p_upper(:,:)
        character(len=16) :: nullhyp='qap'
    end type netcancor_result

    type, public :: network_regression_result
        type(regression_result) :: fit
        real(dp), allocatable :: p_lower(:), p_upper(:), p_two_sided(:)
        real(dp), allocatable :: permutation_distribution(:,:)
        character(len=16) :: nullhyp='classical'
    end type network_regression_result

    type, public :: pstar_result
        type(regression_result) :: fit
        real(dp), allocatable :: tie_data(:,:)
    end type pstar_result

    type, public :: lnam_result
        real(dp), allocatable :: beta(:), rho1(:), rho2(:)
        real(dp), allocatable :: fitted(:), residual(:), disturbances(:)
        real(dp), allocatable :: acvm(:,:), beta_se(:), rho1_se(:), rho2_se(:)
        real(dp) :: sigmasq=0.0_dp, sigma=0.0_dp, loglik=0.0_dp
        real(dp) :: sigmasq_se=0.0_dp, sigma_se=0.0_dp
        logical :: converged=.false.
        integer :: iterations=0, hessian_info=0
    end type lnam_result

    public :: bbnam_probtie, bbnam_fixed_posterior, bbnam_fixed_draws
    public :: bbnam_pooled, bbnam_actor, bbnam_joint_loglik, bbnam_bayes_factor, potscalered_mcmc
    public :: bn_ptriad, bn_dyad_stats, bn_nll_dyad, bn_nll_edge, bn_fit
    public :: brokerage, consensus, nacf, netcancor, netlm, netlogit, lnam
    public :: pstar, pstar_basic, eval_edgeperturbation, npostpred_scalar

abstract interface
        function scalar_graph_function(g) result(v)
            import dp
            real(dp), intent(in) :: g(:,:)
            real(dp) :: v
        end function scalar_graph_function
    end interface

contains

    real(dp) function bbnam_probtie(obs,nprior,em,ep) result(p)
        real(dp),intent(in)::obs(:),nprior
        real(dp),intent(in)::em(:),ep(:)
        real(dp)::num,den
        integer::k
        num=nprior
        den=1.0_dp-nprior
        do k=1,size(obs)
            if(is_missing(obs(k)))cycle
            num=num*(obs(k)*(1.0_dp-em(k))+(1.0_dp-obs(k))*em(k))
            den=den*(obs(k)*ep(k)+(1.0_dp-obs(k))*(1.0_dp-ep(k)))
        end do
        if(num+den<=0.0_dp)then
        p=0.5_dp
        else
        p=num/(num+den)
        end if
    end function bbnam_probtie

    function bbnam_fixed_posterior(dat,nprior,em,ep,diag,mode) result(post)
        real(dp),intent(in)::dat(:,:,:),nprior(:,:),em(:),ep(:)
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        real(dp),allocatable::post(:,:)
        real(dp),allocatable::obs(:)
        integer::n,m,i,j
        logical::dg,undir
        character(len=12)::md
        m=size(dat,1)
        n=size(dat,2)
        allocate(post(n,n),obs(m))
        if(size(em)/=m.or.size(ep)/=m)error stop 'bbnam_fixed_posterior: error-vector length mismatch'
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        do i=1,n
        do j=1,n
            if(.not.dg.and.i==j)then
            post(i,j)=0.0_dp
            cycle
            end if
            obs=dat(:,i,j)
            post(i,j)=bbnam_probtie(obs,nprior(i,j),em,ep)
        end do
        end do
        if(undir)then
            do i=1,n
            do j=i+1,n
            post(j,i)=post(i,j)
            end do
            end do
        end if
    end function bbnam_fixed_posterior

    function bbnam_fixed_draws(dat,nprior,em,ep,draws,diag,mode) result(out)
        real(dp),intent(in)::dat(:,:,:),nprior(:,:),em(:),ep(:)
        integer,intent(in)::draws
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        type(bbnam_result)::out
        real(dp),allocatable::post(:,:)
        integer::s,i,j,n,m
        real(dp)::u
        logical::dg,undir
        character(len=12)::md
        post=bbnam_fixed_posterior(dat,nprior,em,ep,diag,mode)
        n=size(post,1)
        m=size(dat,1)
        allocate(out%net(draws,n,n),out%em(draws,m),out%ep(draws,m))
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        out%em=spread(em,1,draws)
        out%ep=spread(ep,1,draws)
        do s=1,draws
            out%net(s,:,:)=0.0_dp
            do i=1,n
            do j=1,n
                if(.not.dg.and.i==j)cycle
                if(undir.and.j>i)cycle
                call random_number(u)
                out%net(s,i,j)=merge(1.0_dp,0.0_dp,u<post(i,j))
                if(undir)out%net(s,j,i)=out%net(s,i,j)
            end do
            end do
        end do
        out%nactors=n
        out%nobservers=m
        out%draws=draws
        out%model='fixed'
    end function bbnam_fixed_draws

    function bbnam_pooled(dat,nprior,emprior,epprior,draws,burntime,diag,mode) result(out)
        real(dp),intent(in)::dat(:,:,:),nprior(:,:),emprior(2),epprior(2)
        integer,intent(in)::draws,burntime
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        type(bbnam_result)::out
        real(dp),allocatable::a(:,:),post(:,:),emv(:),epv(:),obs(:)
        real(dp)::em,ep,u
        integer::m,n,it,s,i,j,ne0,ne1,np0,np1
        logical::dg,undir
        character(len=12)::md
        m=size(dat,1)
        n=size(dat,2)
        allocate(a(n,n),obs(m),emv(m),epv(m))
        a=0.0_dp
        allocate(out%net(draws,n,n),out%em(draws,1),out%ep(draws,1))
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        call random_number(u)
        em=0.5_dp*u
        call random_number(u)
        ep=0.5_dp*u
        s=0
        do it=1,burntime+draws
            emv=em
            epv=ep
            post=bbnam_fixed_posterior(dat,nprior,emv,epv,dg,md)
            do i=1,n
            do j=1,n
                if(.not.dg.and.i==j)then
                a(i,j)=0
                cycle
                end if
                if(undir.and.j>i)cycle
                call random_number(u)
                a(i,j)=merge(1.0_dp,0.0_dp,u<post(i,j))
                if(undir)a(j,i)=a(i,j)
            end do
            end do
            ne0=0
            ne1=0
            np0=0
            np1=0
            do i=1,n
            do j=1,n
                if(.not.dg.and.i==j)cycle
                if(undir.and.j>i)cycle
                do m=1,size(dat,1)
                    if(is_missing(dat(m,i,j)))cycle
                    if(a(i,j)>0.5_dp)then
                        if(dat(m,i,j)<0.5_dp)then
                        ne0=ne0+1
                        else
                        ne1=ne1+1
                        end if
                    else
                        if(dat(m,i,j)>0.5_dp)then
                        np0=np0+1
                        else
                        np1=np1+1
                        end if
                    end if
                end do
            end do
            end do
            em=beta_rng(emprior(1)+real(ne0,dp),emprior(2)+real(ne1,dp))
            ep=beta_rng(epprior(1)+real(np0,dp),epprior(2)+real(np1,dp))
            if(it>burntime)then
            s=s+1
            out%net(s,:,:)=a
            out%em(s,1)=em
            out%ep(s,1)=ep
            end if
        end do
        out%nactors=n
        out%nobservers=size(dat,1)
        out%draws=draws
        out%model='pooled'
    end function bbnam_pooled

    function bbnam_actor(dat,nprior,emprior,epprior,draws,burntime,diag,mode) result(out)
        real(dp),intent(in)::dat(:,:,:),nprior(:,:),emprior(:,:),epprior(:,:)
        integer,intent(in)::draws,burntime
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        type(bbnam_result)::out
        real(dp),allocatable::a(:,:),post(:,:),em(:),ep(:),obs(:)
        real(dp)::u
        integer::m,n,it,s,i,j,k,ne0,ne1,np0,np1
        logical::dg,undir
        character(len=12)::md
        m=size(dat,1)
        n=size(dat,2)
        if(size(emprior,1)/=m.or.size(epprior,1)/=m)error stop 'bbnam_actor: prior rows must match observers'
        allocate(a(n,n),em(m),ep(m),obs(m))
        a=0.0_dp
        allocate(out%net(draws,n,n),out%em(draws,m),out%ep(draws,m))
        do k=1,m
        call random_number(u)
        em(k)=0.5_dp*u
        call random_number(u)
        ep(k)=0.5_dp*u
        end do
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        s=0
        do it=1,burntime+draws
            allocate(post(n,n))
            do i=1,n
            do j=1,n
                if(.not.dg.and.i==j)then
                post(i,j)=0
                cycle
                end if
                obs=dat(:,i,j)
                post(i,j)=bbnam_probtie(obs,nprior(i,j),em,ep)
            end do
            end do
            if(undir)then
            do i=1,n
            do j=i+1,n
            post(j,i)=post(i,j)
            end do
            end do
            end if
            do i=1,n
            do j=1,n
                if(.not.dg.and.i==j)then
                a(i,j)=0
                cycle
                end if
                if(undir.and.j>i)cycle
                call random_number(u)
                a(i,j)=merge(1.0_dp,0.0_dp,u<post(i,j))
                if(undir)a(j,i)=a(i,j)
            end do
            end do
            deallocate(post)
            do k=1,m
                ne0=0
                ne1=0
                np0=0
                np1=0
                do i=1,n
                do j=1,n
                    if(.not.dg.and.i==j)cycle
                    if(undir.and.j>i)cycle
                    if(is_missing(dat(k,i,j)))cycle
                    if(a(i,j)>0.5_dp)then
                        if(dat(k,i,j)<0.5_dp)then
                        ne0=ne0+1
                        else
                        ne1=ne1+1
                        end if
                    else
                        if(dat(k,i,j)>0.5_dp)then
                        np0=np0+1
                        else
                        np1=np1+1
                        end if
                    end if
                end do
                end do
                em(k)=beta_rng(emprior(k,1)+real(ne0,dp),emprior(k,2)+real(ne1,dp))
                ep(k)=beta_rng(epprior(k,1)+real(np0,dp),epprior(k,2)+real(np1,dp))
            end do
            if(it>burntime)then
            s=s+1
            out%net(s,:,:)=a
            out%em(s,:)=em
            out%ep(s,:)=ep
            end if
        end do
        out%nactors=n
        out%nobservers=m
        out%draws=draws
        out%model='actor'
    end function bbnam_actor

    real(dp) function bbnam_joint_loglik(dat,a,em,ep) result(ll)
        real(dp),intent(in)::dat(:,:,:),a(:,:),em(:),ep(:)
        integer::s,i,j
        real(dp)::p
        ll=0.0_dp
        do s=1,size(dat,1)
        do i=1,size(dat,2)
        do j=1,size(dat,3)
            if(is_missing(dat(s,i,j)))cycle
            p=(1.0_dp-a(i,j))*(dat(s,i,j)*ep(s)+(1.0_dp-dat(s,i,j))*(1.0_dp-ep(s))) + &
              a(i,j)*(dat(s,i,j)*(1.0_dp-em(s))+(1.0_dp-dat(s,i,j))*em(s))
            ll=ll+log(max(p,tiny(1.0_dp)))
        end do
        end do
        end do
    end function bbnam_joint_loglik

    function bbnam_bayes_factor(dat,nprior,em_fp,ep_fp,emprior_pooled,epprior_pooled,emprior_actor,epprior_actor,reps,diag, &
        & mode) result(out)
        real(dp),intent(in)::dat(:,:,:),nprior(:,:),em_fp,ep_fp,emprior_pooled(2),epprior_pooled(2)
        real(dp),intent(in)::emprior_actor(:,:),epprior_actor(:,:)
        integer,intent(in),optional::reps
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        type(bayes_factor_result)::out
        real(dp),allocatable::dd(:,:,:),a(:,:),lf(:),lp(:),la(:),emv(:),epv(:),tmp(:)
        real(dp)::emp,epp,lm2(3),twolm(3)
        integer::m,n,nr,r,k,i,j
        logical::dg,undir
        character(len=12)::md
        m=size(dat,1)
        n=size(dat,2)
        if(size(nprior,1)/=n.or.size(nprior,2)/=n)error stop 'bbnam_bayes_factor: nprior size mismatch'
        if(size(emprior_actor,1)/=m.or.size(epprior_actor,1)/=m.or.size(emprior_actor,2)/=2.or.size(epprior_actor,2)/= &
            & 2)error stop 'bbnam_bayes_factor: actor-prior size mismatch'
        nr=1000
        if(present(reps))nr=reps
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        allocate(dd(m,n,n))
        dd=dat
        if(.not.dg)then
        do k=1,m
        do i=1,n
        dd(k,i,i)=sna_nan()
        end do
        end do
        end if
        if(undir)then
        do k=1,m
        do i=1,n
        do j=i+1,n
        dd(k,i,j)=sna_nan()
        end do
        end do
        end do
        end if
        allocate(lf(nr),lp(nr),la(nr),emv(m),epv(m))
        do r=1,nr
            a=rgraph(n,nprior,md,dg)
            emv=em_fp
            epv=ep_fp
            lf(r)=bbnam_joint_loglik(dd,a,emv,epv)
            emp=beta_rng(emprior_pooled(1),emprior_pooled(2))
            epp=beta_rng(epprior_pooled(1),epprior_pooled(2))
            emv=emp
            epv=epp
            lp(r)=bbnam_joint_loglik(dd,a,emv,epv)
            do k=1,m
            emv(k)=beta_rng(emprior_actor(k,1),emprior_actor(k,2))
            epv(k)=beta_rng(epprior_actor(k,1),epprior_actor(k,2))
            end do
            la(r)=bbnam_joint_loglik(dd,a,emv,epv)
        end do
        out%integrated_loglik=[log_mean(lf),log_mean(lp),log_mean(la)]
        do i=1,3
        do j=1,3
        if(i==j)then
        out%log_bayes_factor(i,j)=out%integrated_loglik(i)
        else
        out%log_bayes_factor(i,j)=out%integrated_loglik(i)-out%integrated_loglik(j)
        end if
        end do
        end do
        lm2=[log_mean(2.0_dp*lf),log_mean(2.0_dp*lp),log_mean(2.0_dp*la)]
        twolm=2.0_dp*out%integrated_loglik
        tmp=log_sub(lm2,twolm)
        out%log_mc_se=(tmp-log(real(nr,dp)))/2.0_dp
        out%reps=nr
    end function bbnam_bayes_factor

    real(dp) function potscalered_mcmc(psi) result(rhat)
        real(dp),intent(in)::psi(:,:)
        real(dp),allocatable::means(:)
        real(dp)::b,w,varp,mtot
        integer::j,nc,nd
        nd=size(psi,1)
        nc=size(psi,2)
        if(nc<2.or.nd<2)then
        rhat=sna_nan()
        return
        end if
        allocate(means(nc))
        do j=1,nc
        means(j)=sum(psi(:,j))/real(nd,dp)
        end do
        mtot=sum(means)/real(nc,dp)
        b=real(nd,dp)/real(nc-1,dp)*sum((means-mtot)**2)
        w=0.0_dp
        do j=1,nc
        w=w+sum((psi(:,j)-means(j))**2)/real(nd-1,dp)
        end do
        w=w/real(nc,dp)
        if(w<=0)then
        rhat=sna_nan()
        return
        end if
        varp=(real(nd-1,dp)/real(nd,dp))*w+b/real(nd,dp)
        rhat=sqrt(varp/w)
    end function potscalered_mcmc

    function bn_ptriad(pi,sigma,rho,d) result(pt)
        real(dp),intent(in)::pi,sigma,rho,d
        real(dp)::pt(16),m0,a0,n0,m1,a1,n1,mp1,ap1,np1,sr
        m0=d*(pi+(1.0_dp-pi)*d)
        a0=d*(1.0_dp-d)*(1.0_dp-pi)
        n0=(1.0_dp-d)*(1.0_dp-d*(1.0_dp-pi))
        m1=(sigma+(1.0_dp-sigma)*d)*(1.0_dp-(1.0_dp-pi)*(1.0_dp-sigma)*(1.0_dp-rho)*(1.0_dp-d))
        a1=(sigma+(1.0_dp-sigma)*d)*(1.0_dp-pi)*(1.0_dp-sigma)*(1.0_dp-rho)*(1.0_dp-d)
        n1=1.0_dp-(sigma+(1.0_dp-sigma)*d)*(1.0_dp+(1.0_dp-pi)*(1.0_dp-sigma)*(1.0_dp-rho)*(1.0_dp-d))
        mp1=sigma*(1.0_dp-(1.0_dp-sigma)*(1.0_dp-rho))
        ap1=sigma*(1.0_dp-sigma)*(1.0_dp-rho)
        np1=1.0_dp-sigma*(1.0_dp-(1.0_dp-sigma)*(1.0_dp-rho)+2.0_dp*(1.0_dp-sigma)*(1.0_dp-rho))
        sr=1.0_dp-(1.0_dp-sigma)*(1.0_dp-rho)
        pt(1)=n0**3
        pt(2)=6*a0*n0*n0
        pt(3)=3*m0*n0*n0
        pt(4)=a0*a0*(n1+2*n0*np1)
        pt(5)=3*a0*a0*n0
        pt(6)=6*a0*a0*n0
        pt(7)=6*m0*a0*n0
        pt(8)=2*m0*a0*(n1+2*n0*np1)
        pt(9)=2*a0*a0*(a1+2*a0*(1-sr)+2*n0*ap1)
        pt(10)=2*a0**3
        pt(11)=m0*m0*(n1+2*n0*np1)
        pt(12)=a0*a0*(m1+2*m0+2*n0*mp1+4*a0*sr)
        pt(13)=m0*a0*(1-sr)*(2*a1+a0*(1-sr)+4*n0*ap1)
        pt(14)=2*m0*a0*(a1+2*a0*(1-sr)+2*n0*ap1)
        pt(15)=m0*(2*m0*a1+2*a0*m1*(1-sr)+2*m0*a0*(1-sr)+4*a0*n0*ap1*sr+4*a0*n0*mp1*(1-sr)+4*m0*n0*ap1+2*a0*a1*sr+6*a0*a0*sr*(1-sr))
        pt(16)=m0*(m0*m1+4*a0*n0*mp1*sr+2*m0*n0*mp1+5*a0*a0*sr*sr+2*a0*m1*sr+2*m0*a0*sr)
    end function bn_ptriad

    function bn_dyad_stats(g) result(stats)
        real(dp),intent(in)::g(:,:)
        real(dp),allocatable::stats(:,:)
        integer::n,i,j,k,parents
        n=size(g,1)
        allocate(stats(n-1,4))
        stats=0.0_dp
        do i=1,n-1
        stats(i,1)=real(i-1,dp)
        end do
        do i=1,n-1
        do j=i+1,n
            parents=0
            do k=1,n
            if(g(k,i)>0.0_dp.and.g(k,j)>0.0_dp)parents=parents+1
            end do
            if(g(i,j)>0.0_dp.and.g(j,i)>0.0_dp)then
            stats(parents+1,2)=stats(parents+1,2)+1
            else if(g(i,j)>0.0_dp.or.g(j,i)>0.0_dp)then
            stats(parents+1,3)=stats(parents+1,3)+1
            else
            stats(parents+1,4)=stats(parents+1,4)+1
            end if
        end do
        end do
    end function bn_dyad_stats

    real(dp) function bn_nll_dyad(params,stats) result(nll)
        real(dp),intent(in)::params(4),stats(:,:)
        real(dp)::pi,sigma,rho,d,lm,la,ln,k,p1,p2
        integer::i
        pi=params(1)
        sigma=params(2)
        rho=params(3)
        d=params(4)
        nll=0.0_dp
        do i=1,size(stats,1)
            k=stats(i,1)
            lm=log(max(1.0_dp-(1.0_dp-pi)*(1.0_dp-rho)**k*(1.0_dp-sigma)**k*(1.0_dp-d),tiny(1.0_dp))) + &
               log(max(1.0_dp-(1.0_dp-sigma)**k*(1.0_dp-d),tiny(1.0_dp)))
            la=log(max(1.0_dp-(1.0_dp-sigma)**k*(1.0_dp-d),tiny(1.0_dp))) + log(max((1.0_dp-pi)*(1.0_dp-rho)**k*(1.0_dp- &
                & sigma)**k*(1.0_dp-d),tiny(1.0_dp)))
            p1=1.0_dp-(1.0_dp-sigma)**k*(1.0_dp-d)
            p2=1.0_dp+(1.0_dp-pi)*(1.0_dp-sigma)**k*(1.0_dp-rho)**k*(1.0_dp-d)
            ln=log(max(1.0_dp-p1*p2,tiny(1.0_dp)))
            nll=nll-(stats(i,2)*lm+stats(i,3)*la+stats(i,4)*ln)
        end do
    end function bn_nll_dyad

    real(dp) function bn_nll_edge(params,stats) result(nll)
        real(dp),intent(in)::params(4),stats(:,:)
        real(dp)::pi,sigma,rho,d,p(4),k
        integer::i
        pi=params(1)
        sigma=params(2)
        rho=params(3)
        d=params(4)
        nll=0
        do i=1,size(stats,1)
            k=stats(i,1)
            p(1)=1-(1-pi)*(1-rho)**k*(1-sigma)**k*(1-d)
            p(2)=1-(1-sigma)**k*(1-d)
            p=max(min(p,1.0_dp-1e-15_dp),1e-15_dp)
            p(3)=1-p(2)
            p(4)=1-p(1)
            nll=nll-(2*stats(i,2)*log(p(1))+stats(i,3)*log(p(2))+stats(i,3)*log(p(3))+2*stats(i,4)*log(p(4)))
        end do
    end function bn_nll_edge

    function bn_fit(g,method,epsilon,maxiter,tol) result(out)
        real(dp),intent(in)::g(:,:)
        character(len=*),intent(in),optional::method
        real(dp),intent(in),optional::epsilon,tol
        integer,intent(in),optional::maxiter
        type(bn_result)::out
        real(dp),allocatable::stats(:,:),tc(:)
        integer,allocatable::tstats(:,:)
        real(dp)::p(4),cand(4),best,v,step,eps,atol,dc(3),choose3,choose2
        integer::i,it,mit
        character(len=16)::meth
        logical::changed
        meth='mple.triad'
        if(present(method))meth=trim(method)
        eps=1e-5_dp
        if(present(epsilon))eps=epsilon
        mit=500
        if(present(maxiter))mit=maxiter
        atol=1e-7_dp
        if(present(tol))atol=tol
        ! sna names are pi,sigma,rho,d; use density as d and reciprocity as pi.
        p=[max(eps,min(1-eps,grecip(g,'edgewise'))),max(eps,min(1-eps,gtrans(g,'weak'))), &
           max(eps,min(1-eps,gtrans(g,'weak'))),max(eps,min(1-eps,gden(g,.false.,'digraph',.true.)))]
        stats=bn_dyad_stats(g)
        tc=triad_census(g,.true.)
        tstats=bn_triad_stats(g)
        best=bn_objective(p,meth,stats,tc,g,tstats)
        step=0.15_dp
        do it=1,mit
            changed=.false.
            do i=1,4
                cand=p
                cand(i)=min(1-eps,p(i)+step)
                v=bn_objective(cand,meth,stats,tc,g,tstats)
                if(v<best)then
                p=cand
                best=v
                changed=.true.
                cycle
                end if
                cand=p
                cand(i)=max(eps,p(i)-step)
                v=bn_objective(cand,meth,stats,tc,g,tstats)
                if(v<best)then
                p=cand
                best=v
                changed=.true.
                end if
            end do
            if(.not.changed)step=step*0.5_dp
            if(step<atol)exit
        end do
        out%pi=p(1)
        out%sigma=p(2)
        out%rho=p(3)
        out%d=p(4)
        out%g_square=2*best
        out%triads=tc
        out%triads_pred=bn_ptriad(p(1),p(2),p(3),p(4))
        dc=dyad_census(g)
        out%dyads=dc
        choose3=real(size(g,1)*(size(g,1)-1)*(size(g,1)-2),dp)/6.0_dp
        choose2=real(size(g,1)*(size(g,1)-1),dp)/2.0_dp
        if(size(g,1)>2.and.choose2>0)then
            out%dyads_pred(1)=sum(out%triads_pred*[0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp,1._dp,0._dp,0._dp,2._dp,1._dp, &
                & 1._dp,1._dp,2._dp,3._dp])*choose3/choose2/real(size(g,1)-2,dp)
            out%dyads_pred(2)=sum(out%triads_pred*[0._dp,1._dp,0._dp,2._dp,2._dp,2._dp,1._dp,1._dp,3._dp,3._dp,0._dp,2._dp, &
                & 2._dp,2._dp,1._dp,0._dp])*choose3/choose2/real(size(g,1)-2,dp)
            out%dyads_pred(3)=sum(out%triads_pred*[3._dp,2._dp,2._dp,1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,0._dp,1._dp,0._dp, &
                & 0._dp,0._dp,0._dp,0._dp])*choose3/choose2/real(size(g,1)-2,dp)
        end if
        out%converged=step<atol
        out%iterations=it
        out%method=meth
    end function bn_fit

    real(dp) function bn_objective(p,method,stats,tc,g,tstats) result(v)
        real(dp),intent(in)::p(4),stats(:,:),tc(:),g(:,:)
        integer,intent(in)::tstats(:,:)
        character(len=*),intent(in)::method
        real(dp)::pt(16)
        if(trim(method)=='mtle')then
            pt=max(bn_ptriad(p(1),p(2),p(3),p(4)),1e-300_dp)
            v=-sum(tc*log(pt))
        else if(trim(method)=='mple.edge')then
            v=bn_nll_edge(p,stats)
        else if(trim(method)=='mple.triad')then
            v=bn_nll_triad(p,g,tstats)
        else
            v=bn_nll_dyad(p,stats)
        end if
    end function bn_objective

    function brokerage(g,cl) result(out)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::cl(:)
        type(brokerage_result)::out
        integer::n,i,j,k,t,ng,gi,gj,gk,m,ii,jj,kk
        integer,allocatable::gid(:)
        real(dp)::d,nn,term
        real(dp),allocatable::ebr(:,:),vbr(:,:),egbr(:),vgbr(:)
        n=size(g,1)
        if(size(g,2)/=n.or.size(cl)/=n)error stop 'brokerage: size mismatch'
        allocate(out%raw(n,6),out%aggregate(1,6))
        out%raw=0.0_dp
        do i=1,n
        do j=1,n
            if(i==j.or.g(i,j)<=0.0_dp)cycle
            do k=1,n
                if(k==i.or.k==j.or.g(j,k)<=0.0_dp.or.g(i,k)>0.0_dp)cycle
                if(cl(j)==cl(i))then
                    if(cl(j)==cl(k))then
                    t=1
                    else
                    t=3
                    end if
                else if(cl(j)==cl(k))then
                t=4
                else if(cl(i)==cl(k))then
                t=2
                else
                t=5
                end if
                out%raw(j,t)=out%raw(j,t)+1.0_dp
            end do
        end do
        end do
        out%raw(:,6)=sum(out%raw(:,1:5),dim=2)
        out%aggregate(1,:)=sum(out%raw,dim=1)

        ! Group identities in first-occurrence order, as in unique(cl).
        allocate(gid(n),out%class_ids(n),out%class_sizes(n))
        gid=0
        out%class_ids=0
        out%class_sizes=0
        ng=0
        do i=1,n
            gi=0
            do j=1,ng
            if(out%class_ids(j)==cl(i))then
            gi=j
            exit
            end if
            end do
            if(gi==0)then
            ng=ng+1
            out%class_ids(ng)=cl(i)
            gi=ng
            end if
            gid(i)=gi
            out%class_sizes(gi)=out%class_sizes(gi)+1
        end do
        out%class_ids=out%class_ids(:ng)
        out%class_sizes=out%class_sizes(:ng)
        allocate(ebr(ng,6),vbr(ng,6),egbr(6),vgbr(6))
        ebr=0.0_dp
        vbr=0.0_dp
        egbr=0.0_dp
        vgbr=0.0_dp
        m=0
        do i=1,n
        do j=1,n
        if(i/=j.and.g(i,j)>0.0_dp)m=m+1
        end do
        end do
        if(n>1)then
        d=real(m,dp)/real(n*(n-1),dp)
        else
        d=0.0_dp
        end if
        nn=real(n,dp)

        do gi=1,ng
            block
                real(dp)::ni,s2,s3,lia
                ni=real(out%class_sizes(gi),dp)
                ebr(gi,1)=d**2*(1-d)*(ni-1)*(ni-2)
                vbr(gi,1)=ebr(gi,1)*(1-d**2*(1-d))+2*(ni-1)*(ni-2)*(ni-3)*d**3*(1-d)**3
                s2=0.0_dp
                s3=0.0_dp
                do gj=1,ng
                if(gj==gi)cycle
                term=real(out%class_sizes(gj),dp)
                s2=s2+term*(term-1)
                s3=s3+term*(term-1)*(term-2)
                end do
                ebr(gi,2)=d**2*(1-d)*s2
                vbr(gi,2)=ebr(gi,2)*(1-d**2*(1-d))+2*s3*d**3*(1-d)**3
                ebr(gi,3)=d**2*(1-d)*(nn-ni)*(ni-1)
                vbr(gi,3)=ebr(gi,3)*(1-d**2*(1-d))+2*((ni-1)*(nn-ni)*(nn-ni-1)/2+(nn-ni)*(ni-1)*(ni-2)/2)*d**3*(1-d)**3
                ebr(gi,4)=ebr(gi,3)
                vbr(gi,4)=vbr(gi,3)
                lia=0.0_dp
                do gj=1,ng
                if(gj==gi)cycle
                do gk=1,ng
                if(gk==gi.or.gk==gj)cycle
                lia=lia+real(out%class_sizes(gj)*out%class_sizes(gk),dp)
                end do
                end do
                ebr(gi,5)=d**2*(1-d)*lia
                term=0.0_dp
                do gj=1,ng
                if(gj==gi)cycle
                term=term+real(out%class_sizes(gj),dp)*real((n-out%class_sizes(gj)-out%class_sizes(gi))*(n- &
                    & out%class_sizes(gj)-out%class_sizes(gi)-1),dp)/2.0_dp
                end do
                vbr(gi,5)=ebr(gi,5)*(1-d**2*(1-d))+4*term*d**3*(1-d)**3
                ebr(gi,6)=d**2*(1-d)*(nn-1)*(nn-2)
                vbr(gi,6)=ebr(gi,6)*(1-d**2*(1-d))+2*(nn-1)*(nn-2)*(nn-3)*d**3*(1-d)**3
            end block
        end do

        allocate(out%expected(n,6),out%sd(n,6),out%z(n,6))
        do i=1,n
        out%expected(i,:)=ebr(gid(i),:)
        out%sd(i,:)=sqrt(max(vbr(gid(i),:),0.0_dp))
        end do
        out%z=0.0_dp
        where(out%sd>0.0_dp)out%z=(out%raw-out%expected)/out%sd

        ! Global expectations and variances from Gould-Fernandez null moments in sna::brokerage.
        do gi=1,ng
            term=real(out%class_sizes(gi),dp)
            egbr(1)=egbr(1)+term*(term-1)*(term-2)
            egbr(2)=egbr(2)+term*(nn-term)*(term-1)
            egbr(3)=egbr(3)+term*(nn-term)*(term-1)
        end do
        egbr(1)=d**2*(1-d)*egbr(1)
        egbr(2)=d**2*(1-d)*egbr(2)
        egbr(3)=d**2*(1-d)*egbr(3)
        egbr(4)=egbr(3)
        term=0.0_dp
        do gi=1,ng
        do gj=1,ng
        term=term+real(out%class_sizes(gi)*out%class_sizes(gj)*(n-out%class_sizes(gi)-out%class_sizes(gj)),dp)
        end do
        end do
        do gi=1,ng
        term=term-real(out%class_sizes(gi)*out%class_sizes(gi)*(n-2*out%class_sizes(gi)),dp)
        end do
        egbr(5)=d**2*(1-d)*term
        egbr(6)=d**2*(1-d)*nn*(nn-1)*(nn-2)

        term=0.0_dp
        do gi=1,ng
            block
                real(dp)::ni
                ni=real(out%class_sizes(gi),dp)
                term=term+ni*(ni-1)*(ni-2)*((4*ni-10)*d**3*(1-d)**3-4*(ni-3)*d**4*(1-d)**2+(ni-3)*d**5*(1-d))
            end block
        end do
        vgbr(1)=egbr(1)*(1-d**2*(1-d))+term
        term=0.0_dp
        do gi=1,ng
        do gj=1,ng
            term=term+real(out%class_sizes(gi)*out%class_sizes(gj)*(out%class_sizes(gi)-1),dp)*((2.0_dp*out%class_sizes(gi)+ &
                & 2.0_dp*out%class_sizes(gj)-6.0_dp)*d**3*(1-d)**3+(nn-out%class_sizes(gi)-1.0_dp)*d**5*(1-d))
        end do
        end do
        do gi=1,ng
            term=term-real(out%class_sizes(gi)**2*(out%class_sizes(gi)-1),dp)*((4.0_dp*out%class_sizes(gi)-6.0_dp)*d**3*(1- &
                & d)**3+(nn-out%class_sizes(gi)-1.0_dp)*d**5*(1-d))
        end do
        vgbr(2)=egbr(2)*(1-d**2*(1-d))+term
        term=0.0_dp
        do gi=1,ng
            block
                real(dp)::ni
                ni=real(out%class_sizes(gi),dp)
                term=term+ni*(nn-ni)*(ni-1)*((nn-3)*d**3*(1-d)**3+(ni-2)*d**5*(1-d))
            end block
        end do
        vgbr(3)=egbr(3)*(1-d**2*(1-d))+term
        vgbr(4)=vgbr(3)
        vgbr(5)=egbr(5)*(1-d**2*(1-d))
        do ii=1,ng
        do jj=1,ng
        do kk=1,ng
            if(ii==jj.or.jj==kk.or.ii==kk)cycle
            block
                real(dp)::ni,nj,nk
                ni=real(out%class_sizes(ii),dp)
                nj=real(out%class_sizes(jj),dp)
                nk=real(out%class_sizes(kk),dp)
                vgbr(5)=vgbr(5)+ni*nj*nk*((4*(nn-nj)-2*(ni+nk+1))*d**3*(1-d)**3-(4*(nn-nk)-2*(ni+nj+1))*d**4*(1-d)**2+(nn- &
                    & (ni+nk+1))*d**5*(1-d))
            end block
        end do
        end do
        end do
        vgbr(6)=egbr(6)*(1-d**2*(1-d))+nn*(nn-1)*(nn-2)*((4*nn-10)*d**3*(1-d)**3-4*(nn-3)*d**4*(1-d)**2+(nn-3)*d**5*(1-d))

        allocate(out%global_expected(6),out%global_sd(6),out%global_z(6),out%group_expected(ng,6),out%group_sd(ng,6))
        out%global_expected=egbr
        out%global_sd=sqrt(max(vgbr,0.0_dp))
        out%global_z=0.0_dp
        where(out%global_sd>0.0_dp)out%global_z=(out%aggregate(1,:)-egbr)/out%global_sd
        out%group_expected=ebr
        out%group_sd=sqrt(max(vbr,0.0_dp))
    end function brokerage

    function consensus(dat,method,mode,diag,tol,maxiter,no_bias) result(cong)
        real(dp),intent(in)::dat(:,:,:)
        character(len=*),intent(in),optional::method,mode
        logical,intent(in),optional::diag,no_bias
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter
        real(dp),allocatable::cong(:,:),gc(:,:),comp(:),bias(:),ocomp(:),weights(:),vec(:)
        real(dp)::atol,cdiff,ll1,ll0,s1,s0,correct,drate,eigval
        integer::m,n,i,j,k,it,mit
        character(len=24)::meth,md
        logical::dg,nb,ans
        meth='central.graph'
        if(present(method))meth=trim(method)
        md='digraph'
        if(present(mode))md=trim(mode)
        dg=.false.
        if(present(diag))dg=diag
        nb=.false.
        if(present(no_bias))nb=no_bias
        atol=1e-6_dp
        if(present(tol))atol=tol
        mit=1000
        if(present(maxiter))mit=maxiter
        m=size(dat,1)
        n=size(dat,2)
        if(trim(meth)=='central.graph')then
            cong=centralgraph(dat,.false.)
        else if(trim(meth)=='single.reweight'.or.trim(meth)=='PCA.reweight')then
            allocate(gc(m,m))
            do i=1,m
            do j=1,m
            gc(i,j)=graph_correlation(dat(i,:,:),dat(j,:,:),dg,md)
            if(is_missing(gc(i,j)))gc(i,j)=0
            end do
            gc(i,i)=1
            end do
            allocate(weights(m))
            if(trim(meth)=='single.reweight')then
                weights=sum(gc,dim=2)
            else
                allocate(vec(m))
                call dominant_eigenvector(gc,vec,eigval)
                weights=abs(vec)
            end if
            if(abs(sum(weights))<=sna_eps)weights=1
            weights=weights/sum(weights)
            allocate(cong(n,n))
            cong=0
            do k=1,m
            cong=cong+weights(k)*dat(k,:,:)
            end do
        else if(trim(meth)=='iterative.reweight'.or.trim(meth)=='romney.batchelder')then
            cong=centralgraph(dat,.false.)
            allocate(comp(m),bias(m),ocomp(m))
            comp=0.5_dp
            bias=0.5_dp
            do k=1,m
                correct=0
                drate=0
                it=0
                do i=1,n
                do j=1,n
                if(.not.dg.and.i==j)cycle
                if(is_missing(dat(k,i,j)))cycle
                it=it+1
                if((dat(k,i,j)>0.5_dp).eqv.(cong(i,j)>0.5_dp))correct=correct+1
                drate=drate+dat(k,i,j)
                end do
                end do
                if(it>0)then
                correct=correct/it
                drate=drate/it
                comp(k)=max(0.5_dp,correct)
                if(trim(meth)=='romney.batchelder')comp(k)=max(0.0_dp,2*correct-1)
                if(.not.nb)bias(k)=drate
                end if
            end do
            do it=1,mit
                ocomp=comp
                do i=1,n
                do j=1,n
                    if(.not.dg.and.i==j)then
                    cong(i,j)=0
                    cycle
                    end if
                    ll1=0
                    ll0=0
                    do k=1,m
                        if(is_missing(dat(k,i,j)))cycle
                        ll1=ll1+dat(k,i,j)*log(max(comp(k)+(1-comp(k))*bias(k),1e-15_dp))+(1-dat(k,i,j))*log(max((1-comp(k))* &
                            & (1-bias(k)),1e-15_dp))
                        ll0=ll0+(1-dat(k,i,j))*log(max(comp(k)+(1-comp(k))*(1-bias(k)),1e-15_dp))+dat(k,i,j)*log(max((1- &
                            & comp(k))*bias(k),1e-15_dp))
                    end do
                    cong(i,j)=merge(1.0_dp,0.0_dp,ll1>ll0)
                end do
                end do
                do k=1,m
                    correct=0
                    drate=0
                    j=0
                    do i=1,n
                    do j=1,n
                        if(.not.dg.and.i==j)cycle
                        if(is_missing(dat(k,i,j)))cycle
                        correct=correct+merge(1.0_dp,0.0_dp,(dat(k,i,j)>0.5_dp).eqv.(cong(i,j)>0.5_dp))
                        drate=drate+dat(k,i,j)
                    end do
                    end do
                    ! Use all admissible dyads; this matches the iterative-reweight update closely.
                    correct=correct/real(max(1,n*n-merge(0,n,dg)),dp)
                    comp(k)=max(0.5_dp,min(1.0_dp,correct))
                end do
                cdiff=sum(abs(ocomp-comp))
                if(cdiff<=atol)exit
            end do
        else if(trim(meth)=='LAS.intersection'.or.trim(meth)=='LAS.union')then
            allocate(cong(n,n))
            cong=0
            do i=1,n
            do j=1,n
                if(i<=m.and.j<=m)then
                    if(trim(meth)=='LAS.intersection')then
                    cong(i,j)=merge(1.0_dp,0.0_dp,dat(i,i,j)>0.0_dp.and.dat(j,i,j)>0.0_dp)
                    else
                    cong(i,j)=merge(1.0_dp,0.0_dp,dat(i,i,j)>0.0_dp.or.dat(j,i,j)>0.0_dp)
                    end if
                end if
            end do
            end do
        else if(trim(meth)=='OR.row')then
            allocate(cong(n,n))
            cong=0.0_dp
            do i=1,min(n,m)
            cong(i,:)=dat(i,i,:)
            end do
        else if(trim(meth)=='OR.col')then
            allocate(cong(n,n))
            cong=0.0_dp
            do i=1,min(n,m)
            cong(:,i)=dat(i,:,i)
            end do
        else
            cong=centralgraph(dat,.false.)
        end if
        if(trim(md)=='graph')then
        do i=1,n
        do j=i+1,n
        cong(i,j)=cong(j,i)
        end do
        end do
        end if
        if(.not.dg)then
        do i=1,n
        cong(i,i)=0
        end do
        end if
    end function consensus

    function nacf(net,y,lag_max,type,neighborhood_type,partial,mode,diag,thresh,demean) result(v)
        real(dp),intent(in)::net(:,:),y(:)
        integer,intent(in),optional::lag_max
        character(len=*),intent(in),optional::type,neighborhood_type,mode
        logical,intent(in),optional::partial,diag,demean
        real(dp),intent(in),optional::thresh
        real(dp),allocatable::v(:),yy(:),nh(:,:)
        integer::lag,n,lm,i,j,ec
        real(dp)::vary,num,den,thr
        character(len=16)::tp,nt,md
        logical::pt,dg,dm
        n=size(y)
        lm=n-1
        if(present(lag_max))lm=min(lag_max,n-1)
        tp='correlation'
        if(present(type))tp=trim(type)
        nt='in'
        if(present(neighborhood_type))nt=trim(neighborhood_type)
        md='digraph'
        if(present(mode))md=trim(mode)
        pt=.true.
        if(present(partial))pt=partial
        dg=.false.
        if(present(diag))dg=diag
        thr=0
        if(present(thresh))thr=thresh
        dm=.true.
        if(present(demean))dm=demean
        allocate(yy(n))
        yy=y
        if(dm.or.trim(tp)=='moran')yy=yy-sum(yy)/real(n,dp)
        vary=sum((yy-sum(yy)/real(n,dp))**2)/real(max(1,n-1),dp)
        allocate(v(lm+1))
        v=0
        if(trim(tp)=='covariance')v(1)=dot_product(yy,yy)/real(n,dp)
        if(trim(tp)=='correlation'.or.trim(tp)=='moran')v(1)=1
        do lag=1,lm
            nh=neighborhood(net,lag,nt,md,dg,thr,pt)
            ec=count(nh>0.0_dp)
            if(ec==0)cycle
            select case(trim(tp))
            case('covariance')
            v(lag+1)=dot_product(yy,matmul(nh,yy))/real(ec,dp)
            case('correlation')
            v(lag+1)=dot_product(yy,matmul(nh,yy))/real(ec,dp)/max(vary,sna_eps)
            case('moran')
            v(lag+1)=real(n,dp)/real(ec,dp)*sum(spread(yy,2,n)*spread(yy,1,n)*nh)/max(sum(yy**2),sna_eps)
            case('geary')
                num=0
                do i=1,n
                do j=1,n
                if(nh(i,j)>0)num=num+nh(i,j)*(yy(i)-yy(j))**2
                end do
                end do
                den=sum((yy-sum(yy)/real(n,dp))**2)
                v(lag+1)=real(n-1,dp)/(2*real(ec,dp))*num/max(den,sna_eps)
            end select
        end do
    end function nacf

    function netcancor(y,x,mode,diag,nullhyp,reps) result(out)
        ! Canonical correlations among network variables, including the Monte Carlo
        ! nulls used by sna::netcancor.  Graph stacks are (variables,n,n).
        real(dp),intent(in)::y(:,:,:),x(:,:,:)
        character(len=*),intent(in),optional::mode,nullhyp
        logical,intent(in),optional::diag
        integer,intent(in),optional::reps
        type(netcancor_result)::out
        real(dp),allocatable::ymat(:,:),xmat(:,:),yr(:,:,:),xr(:,:,:),p(:,:),gr(:,:)
        real(dp),allocatable::cor(:),xc(:,:),yc(:,:)
        integer::my,mx,n,q,nr,r,k,i,j
        logical::dg
        character(len=16)::md,nh
        my=size(y,1)
        mx=size(x,1)
        n=size(y,2)
        if(size(y,3)/=n.or.size(x,2)/=n.or.size(x,3)/=n)error stop 'netcancor: graph orders differ'
        q=min(mx,my)
        nr=1000
        if(present(reps))nr=reps
        md='digraph'
        if(present(mode))md=trim(mode)
        dg=.false.
        if(present(diag))dg=diag
        nh='cugtie'
        if(present(nullhyp))nh=trim(nullhyp)
        out%nullhyp=nh
        call stack_design(x,xmat)
        call stack_design(y,ymat)
        call canonical_fit(xmat,ymat,cor,xc,yc)
        allocate(out%cor(q),out%xcoef(mx,q),out%ycoef(my,q))
        out%cor=cor
        out%xcoef=xc
        out%ycoef=yc
        allocate(out%cor_distribution(nr,q),out%xcoef_distribution(nr,mx,q),out%ycoef_distribution(nr,my,q))
        allocate(xr(mx,n,n),yr(my,n,n),p(n,n))
        do r=1,nr
            do k=1,mx
                select case(trim(nh))
                case('qap')
                xr(k,:,:)=rmperm(x(k,:,:))
                case('cug')
                    p=0.5_dp
                    xr(k,:,:)=rgraph(n,p,md,dg)
                case('cugden')
                    p=gden(x(k,:,:),dg,md,.true.)
                    xr(k,:,:)=rgraph(n,p,md,dg)
                case default
                    xr(k,:,:)=resample_graph_values(x(k,:,:),md,dg)
                end select
            end do
            do k=1,my
                select case(trim(nh))
                case('qap')
                yr(k,:,:)=rmperm(y(k,:,:))
                case('cug')
                    p=0.5_dp
                    yr(k,:,:)=rgraph(n,p,md,dg)
                case('cugden')
                    p=gden(y(k,:,:),dg,md,.true.)
                    yr(k,:,:)=rgraph(n,p,md,dg)
                case default
                    yr(k,:,:)=resample_graph_values(y(k,:,:),md,dg)
                end select
            end do
            call stack_design(xr,xmat)
            call stack_design(yr,ymat)
            call canonical_fit(xmat,ymat,cor,xc,yc)
            out%cor_distribution(r,:)=cor
            out%xcoef_distribution(r,:,:)=xc
            out%ycoef_distribution(r,:,:)=yc
        end do
        allocate(out%cor_p_lower(q),out%cor_p_upper(q),out%xcoef_p_lower(mx,q),out%xcoef_p_upper(mx,q),out%ycoef_p_lower(my, &
            & q),out%ycoef_p_upper(my,q))
        do k=1,q
            out%cor_p_lower(k)=real(count(out%cor_distribution(:,k)<=out%cor(k)),dp)/real(nr,dp)
            out%cor_p_upper(k)=real(count(out%cor_distribution(:,k)>=out%cor(k)),dp)/real(nr,dp)
        end do
        do i=1,mx
        do k=1,q
            out%xcoef_p_lower(i,k)=real(count(out%xcoef_distribution(:,i,k)<=out%xcoef(i,k)),dp)/real(nr,dp)
            out%xcoef_p_upper(i,k)=real(count(out%xcoef_distribution(:,i,k)>=out%xcoef(i,k)),dp)/real(nr,dp)
        end do
        end do
        do i=1,my
        do k=1,q
            out%ycoef_p_lower(i,k)=real(count(out%ycoef_distribution(:,i,k)<=out%ycoef(i,k)),dp)/real(nr,dp)
            out%ycoef_p_upper(i,k)=real(count(out%ycoef_distribution(:,i,k)>=out%ycoef(i,k)),dp)/real(nr,dp)
        end do
        end do
    end function netcancor

    subroutine stack_design(g,z)
        real(dp),intent(in)::g(:,:,:)
        real(dp),allocatable,intent(out)::z(:,:)
        integer::m,n,k,i,j,c
        m=size(g,1)
        n=size(g,2)
        allocate(z(n*n,m))
        do k=1,m
            c=0
            do j=1,n
            do i=1,n
            c=c+1
            z(c,k)=g(k,i,j)
            end do
            end do
        end do
    end subroutine stack_design

    subroutine canonical_fit(x,y,cor,xcoef,ycoef)
        real(dp),intent(in)::x(:,:),y(:,:)
        real(dp),allocatable,intent(out)::cor(:),xcoef(:,:),ycoef(:,:)
        real(dp),allocatable::xc(:,:),yc(:,:),sxx(:,:),syy(:,:),sxy(:,:),ixh(:,:),iyh(:,:),mm(:,:),aa(:,:),eval(:),u(:,:),vv(:)
        integer::n,px,py,q,i,info
        n=size(x,1)
        px=size(x,2)
        py=size(y,2)
        q=min(px,py)
        allocate(xc(n,px),yc(n,py))
        xc=x-spread(sum(x,dim=1)/real(n,dp),1,n)
        yc=y-spread(sum(y,dim=1)/real(n,dp),1,n)
        sxx=matmul(transpose(xc),xc)
        syy=matmul(transpose(yc),yc)
        sxy=matmul(transpose(xc),yc)
        call symmetric_invsqrt(sxx,ixh)
        call symmetric_invsqrt(syy,iyh)
        mm=matmul(matmul(ixh,sxy),iyh)
        aa=matmul(mm,transpose(mm))
        allocate(eval(px),u(px,px))
        call jacobi_eigen_symmetric(aa,eval,u,info)
        allocate(cor(q),xcoef(px,q),ycoef(py,q),vv(py))
        do i=1,q
            cor(i)=sqrt(max(eval(i),0.0_dp))
            xcoef(:,i)=matmul(ixh,u(:,i))
            if(cor(i)>1.0e-12_dp)then
            vv=matmul(transpose(mm),u(:,i))/cor(i)
            else
            vv=0.0_dp
            end if
            ycoef(:,i)=matmul(iyh,vv)
            ! R's cancor coefficient signs are arbitrary; standardize to a reproducible orientation.
            if(sum(xcoef(:,i))<0.0_dp)then
            xcoef(:,i)=-xcoef(:,i)
            ycoef(:,i)=-ycoef(:,i)
            end if
        end do
    end subroutine canonical_fit

    subroutine symmetric_invsqrt(a,aih)
        real(dp),intent(in)::a(:,:)
        real(dp),allocatable,intent(out)::aih(:,:)
        real(dp),allocatable::ev(:),vec(:,:),d(:,:)
        real(dp)::mx,tol
        integer::n,i,info
        n=size(a,1)
        allocate(ev(n),vec(n,n),d(n,n),aih(n,n))
        call jacobi_eigen_symmetric(a,ev,vec,info)
        d=0.0_dp
        mx=max(1.0_dp,maxval(abs(ev)))
        tol=1.0e-12_dp*mx
        do i=1,n
        if(ev(i)>tol)d(i,i)=1.0_dp/sqrt(ev(i))
        end do
        aih=matmul(matmul(vec,d),transpose(vec))
    end subroutine symmetric_invsqrt

    function resample_graph_values(g,mode,diag) result(out)
        real(dp),intent(in)::g(:,:)
        character(len=*),intent(in)::mode
        logical,intent(in)::diag
        real(dp),allocatable::out(:,:),vals(:)
        integer::n,i,j,k,ix
        real(dp)::u
        n=size(g,1)
        allocate(out(n,n),vals(n*n))
        k=0
        do j=1,n
        do i=1,n
        k=k+1
        vals(k)=g(i,j)
        end do
        end do
        do j=1,n
        do i=1,n
        call random_number(u)
        ix=1+min(n*n-1,int(u*real(n*n,dp)))
        out(i,j)=vals(ix)
        end do
        end do
        if(.not.diag)then
        do i=1,n
        out(i,i)=0.0_dp
        end do
        end if
        if(trim(mode)/='digraph')then
        do i=1,n
        do j=i+1,n
        out(i,j)=out(j,i)
        end do
        end do
        end if
    end function resample_graph_values

    function netlm(y,x,intercept,mode,diag,nullhyp,reps) result(out)
        real(dp),intent(in)::y(:,:),x(:,:,:)
        logical,intent(in),optional::intercept,diag
        character(len=*),intent(in),optional::mode,nullhyp
        integer,intent(in),optional::reps
        type(network_regression_result)::out
        real(dp),allocatable::yv(:),xv(:,:),b(:),cov(:,:),yr(:,:),xr(:,:,:),v(:),cv(:,:)
        real(dp)::sig
        integer::nobs,p,k,info,nr,r
        logical::ic,dg
        character(len=16)::md,nh
        ic=.true.
        if(present(intercept))ic=intercept
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        nh='classical'
        if(present(nullhyp))nh=trim(nullhyp)
        nr=1000
        if(present(reps))nr=reps
        call graph_design(y,x,ic,md,dg,yv,xv)
        nobs=size(yv)
        p=size(xv,2)
        allocate(b(p),cov(p,p))
        call least_squares(xv,yv,b,cov,sig,info)
        call fill_regression(out%fit,xv,yv,b,cov,sig,info,.false.)
        out%nullhyp=nh
        allocate(out%p_lower(p),out%p_upper(p),out%p_two_sided(p))
        out%p_lower=0
        out%p_upper=0
        out%p_two_sided=0
        if(trim(nh)/='classical')then
            allocate(out%permutation_distribution(nr,p),xr(size(x,1),size(x,2),size(x,3)))
            do r=1,nr
                xr=x
                if(trim(nh)=='qapy'.or.trim(nh)=='qap')then
                    if(allocated(yr))deallocate(yr)
                    yr=rmperm(y)
                else
                    if(allocated(yr))deallocate(yr)
                    allocate(yr(size(y,1),size(y,2)))
                    yr=y
                    do k=1,size(x,1)
                    xr(k,:,:)=rmperm(x(k,:,:))
                    end do
                end if
                call graph_design(yr,xr,ic,md,dg,yv,xv)
                allocate(v(size(xv,2)),cv(size(xv,2),size(xv,2)))
                call least_squares(xv,yv,v,cv,sig,info)
                out%permutation_distribution(r,:)=v
                deallocate(v,cv,yr)
            end do
            do k=1,p
                out%p_lower(k)=real(count(out%permutation_distribution(:,k)<=out%fit%coef(k)),dp)/nr
                out%p_upper(k)=real(count(out%permutation_distribution(:,k)>=out%fit%coef(k)),dp)/nr
                out%p_two_sided(k)=real(count(abs(out%permutation_distribution(:,k))>=abs(out%fit%coef(k))),dp)/nr
            end do
        end if
    end function netlm

    function netlogit(y,x,intercept,mode,diag,nullhyp,reps) result(out)
        real(dp),intent(in)::y(:,:),x(:,:,:)
        logical,intent(in),optional::intercept,diag
        character(len=*),intent(in),optional::mode,nullhyp
        integer,intent(in),optional::reps
        type(network_regression_result)::out
        real(dp),allocatable::yv(:),xv(:,:),b(:),cov(:,:),yr(:,:),xr(:,:,:),v(:),cv(:,:)
        real(dp)::ll
        integer::p,k,info,nr,r
        logical::ic,dg
        character(len=16)::md,nh
        ic=.true.
        if(present(intercept))ic=intercept
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        nh='classical'
        if(present(nullhyp))nh=trim(nullhyp)
        nr=1000
        if(present(reps))nr=reps
        call graph_design(y,x,ic,md,dg,yv,xv)
        p=size(xv,2)
        allocate(b(p),cov(p,p))
        call logistic_irls(xv,yv,b,cov,ll,info)
        call fill_regression(out%fit,xv,yv,b,cov,0.0_dp,info,.true.)
        out%fit%loglik=ll
        out%nullhyp=nh
        allocate(out%p_lower(p),out%p_upper(p),out%p_two_sided(p))
        out%p_lower=0
        out%p_upper=0
        out%p_two_sided=0
        if(trim(nh)/='classical')then
            allocate(out%permutation_distribution(nr,p),xr(size(x,1),size(x,2),size(x,3)),yr(size(y,1),size(y,2)))
            do r=1,nr
                xr=x
                yr=y
                if(trim(nh)=='qapy'.or.trim(nh)=='qap')then
                yr=rmperm(y)
                else
                do k=1,size(x,1)
                xr(k,:,:)=rmperm(x(k,:,:))
                end do
                end if
                call graph_design(yr,xr,ic,md,dg,yv,xv)
                allocate(v(size(xv,2)),cv(size(xv,2),size(xv,2)))
                call logistic_irls(xv,yv,v,cv,ll,info)
                out%permutation_distribution(r,:)=v
                deallocate(v,cv)
            end do
            do k=1,p
                out%p_lower(k)=real(count(out%permutation_distribution(:,k)<=out%fit%coef(k)),dp)/nr
                out%p_upper(k)=real(count(out%permutation_distribution(:,k)>=out%fit%coef(k)),dp)/nr
                out%p_two_sided(k)=real(count(abs(out%permutation_distribution(:,k))>=abs(out%fit%coef(k))),dp)/nr
            end do
        end if
    end function netlogit

    function lnam(y,x,w1,w2,tol,maxiter) result(out)
        real(dp),intent(in)::y(:)
        real(dp),intent(in),optional::x(:,:),w1(:,:,:),w2(:,:,:)
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter
        type(lnam_result)::out
        real(dp),allocatable::rho1(:),rho2(:),beta(:),a1(:,:),a2(:,:),z(:),xx(:,:),cov(:,:),res(:),cand(:),se(:)
        real(dp)::sig,olddev,dev,step,atol,v
        integer::n,nx,nw1,nw2,info,it,mit,k
        logical::changed
        n=size(y)
        nx=0
        if(present(x))nx=size(x,2)
        nw1=0
        if(present(w1))nw1=size(w1,1)
        nw2=0
        if(present(w2))nw2=size(w2,1)
        allocate(rho1(nw1),rho2(nw2),beta(nx))
        rho1=0
        rho2=0
        beta=0
        atol=1e-8_dp
        if(present(tol))atol=tol
        mit=200
        if(present(maxiter))mit=maxiter
        olddev=huge(1.0_dp)
        step=0.1_dp
        do it=1,mit
            call lnam_mats(n,w1,w2,rho1,rho2,a1,a2)
            if(nx>0)then
                if(nw2>0)then
                xx=matmul(a2,x)
                z=matmul(a2,matmul(a1,y))
                else
                xx=x
                z=matmul(a1,y)
                end if
                allocate(cov(nx,nx))
                call least_squares(xx,z,beta,cov,sig,info)
                deallocate(cov)
            end if
            z=matmul(a1,y)
            if(nx>0)z=z-matmul(x,beta)
            if(nw2>0)z=matmul(a2,z)
            sig=sum(z*z)/real(n,dp)
            dev=lnam_deviance(y,x,w1,w2,beta,rho1,rho2,max(sig,1e-15_dp))
            changed=.false.
            do k=1,nw1
                v=rho1(k)
                rho1(k)=v+step
                if(lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)<dev)then
                dev=lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)
                changed=.true.
                cycle
                end if
                rho1(k)=v-step
                if(lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)<dev)then
                dev=lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)
                changed=.true.
                else
                rho1(k)=v
                end if
            end do
            do k=1,nw2
                v=rho2(k)
                rho2(k)=v+step
                if(lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)<dev)then
                dev=lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)
                changed=.true.
                cycle
                end if
                rho2(k)=v-step
                if(lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)<dev)then
                dev=lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)
                changed=.true.
                else
                rho2(k)=v
                end if
            end do
            if(.not.changed)step=step*0.5_dp
            if(abs(olddev-dev)<atol.and.step<sqrt(atol))exit
            olddev=dev
        end do
        call lnam_mats(n,w1,w2,rho1,rho2,a1,a2)
        z=matmul(a1,y)
        if(nx>0)z=z-matmul(x,beta)
        if(nw2>0)z=matmul(a2,z)
        sig=sum(z*z)/real(n,dp)
        out%beta=beta
        out%rho1=rho1
        out%rho2=rho2
        out%sigmasq=sig
        out%sigma=sqrt(sig)
        out%loglik=-0.5_dp*lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sig)
        call lnam_information(y,x,w1,w2,beta,rho1,rho2,sig,out%acvm,se,out%hessian_info)
        allocate(out%beta_se(nx),out%rho1_se(nw1),out%rho2_se(nw2))
        if(nx>0)out%beta_se=se(1:nx)
        if(nw1>0)out%rho1_se=se(nx+1:nx+nw1)
        if(nw2>0)out%rho2_se=se(nx+nw1+1:nx+nw1+nw2)
        out%sigmasq_se=se(nx+nw1+nw2+1)
        if(sig>0.0_dp)out%sigma_se=out%sigmasq_se**2/(4.0_dp*sig)
        allocate(out%fitted(n),out%residual(n),out%disturbances(n))
        if(nx>0)then
        out%fitted=matmul(x,beta)
        if(nw1>0)then
        allocate(res(n))
        call solve_fit(a1,out%fitted,res)
        out%fitted=res
        deallocate(res)
        end if
        else
        out%fitted=0
        end if
        out%residual=y-out%fitted
        out%disturbances=z
        out%converged=step<sqrt(atol)
        out%iterations=it
    end function lnam

    function pstar(dat,effects,attr,memb,diag,mode) result(out)
        real(dp),intent(in)::dat(:,:)
        character(len=*),intent(in)::effects(:)
        real(dp),intent(in),optional::attr(:,:)
        integer,intent(in),optional::memb(:,:)
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        type(pstar_result)::out
        real(dp),allocatable::x(:,:),y(:),g0(:,:),g1(:,:),b(:),cov(:,:),chg(:)
        integer::n,i,j,k,row,p,nobs,info,col,w,na,nm
        real(dp)::ll
        logical::dg,undir
        character(len=12)::md
        n=size(dat,1)
        if(size(dat,2)/=n)error stop 'pstar: square graph required'
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        if(present(attr))then
        if(size(attr,1)/=n)error stop 'pstar: attr row count mismatch'
        na=size(attr,2)
        else
        na=0
        end if
        if(present(memb))then
        if(size(memb,1)/=n)error stop 'pstar: memb row count mismatch'
        nm=size(memb,2)
        else
        nm=0
        end if
        p=na+nm
        do k=1,size(effects)
        p=p+pstar_effect_width(trim(effects(k)),n)
        end do
        nobs=0
        do i=1,n
        do j=1,n
        if(.not.dg.and.i==j)cycle
        if(undir.and.j>i)cycle
        if(.not.is_missing(dat(i,j)))nobs=nobs+1
        end do
        end do
        allocate(x(nobs,p),y(nobs),out%tie_data(nobs,p+1))
        row=0
        do i=1,n
        do j=1,n
            if(.not.dg.and.i==j)cycle
            if(undir.and.j>i)cycle
            if(is_missing(dat(i,j)))cycle
            row=row+1
            y(row)=merge(1.0_dp,0.0_dp,dat(i,j)>0.0_dp)
            allocate(g0(n,n),g1(n,n))
            g0=dat
            g1=dat
            g0(i,j)=0
            g1(i,j)=1
            if(undir)then
            g0(j,i)=0
            g1(j,i)=1
            end if
            col=0
            do k=1,size(effects)
                w=pstar_effect_width(trim(effects(k)),n)
                allocate(chg(w))
                call pstar_effect_change(g0,g1,trim(effects(k)),dg,md,chg)
                x(row,col+1:col+w)=chg
                col=col+w
                deallocate(chg)
            end do
            if(present(attr))then
            x(row,col+1:col+na)=abs(attr(i,:)-attr(j,:))
            col=col+na
            end if
            if(present(memb))then
                do k=1,nm
                x(row,col+k)=merge(1.0_dp,0.0_dp,memb(i,k)==memb(j,k))
                end do
                col=col+nm
            end if
            deallocate(g0,g1)
        end do
        end do
        out%tie_data(:,1)=y
        out%tie_data(:,2:)=x
        allocate(b(p),cov(p,p))
        call logistic_irls(x,y,b,cov,ll,info)
        call fill_regression(out%fit,x,y,b,cov,0.0_dp,info,.true.)
        out%fit%loglik=ll
    end function pstar

    integer function pstar_effect_width(effect,n) result(w)
        character(len=*),intent(in)::effect
        integer,intent(in)::n
        select case(trim(effect))
        case('outdegree','indegree','betweenness','closeness')
        w=n
        case default
        w=1
        end select
    end function pstar_effect_width

    subroutine pstar_effect_change(g0,g1,effect,diag,mode,chg)
        real(dp),intent(in)::g0(:,:),g1(:,:)
        character(len=*),intent(in)::effect,mode
        logical,intent(in)::diag
        real(dp),intent(out)::chg(:)
        real(dp),allocatable::v0(:),v1(:)
        integer::n
        n=size(g0,1)
        select case(trim(effect))
        case('outdegree')
            if(trim(mode)=='graph')then
            v0=degree(g0,'indegree',diag,.false.)
            v1=degree(g1,'indegree',diag,.false.)
            else
            v0=degree(g0,'outdegree',diag,.false.)
            v1=degree(g1,'outdegree',diag,.false.)
            end if
            chg=v1-v0
        case('indegree')
            v0=degree(g0,'indegree',diag,.false.)
            v1=degree(g1,'indegree',diag,.false.)
            chg=v1-v0
        case('betweenness')
            if(trim(mode)=='graph')then
            v0=betweenness(g0,'undirected')
            v1=betweenness(g1,'undirected')
            else
            v0=betweenness(g0,'directed')
            v1=betweenness(g1,'directed')
            end if
            chg=v1-v0
        case('closeness')
            if(trim(mode)=='graph')then
            v0=closeness(g0,'undirected')
            v1=closeness(g1,'undirected')
            else
            v0=closeness(g0,'directed')
            v1=closeness(g1,'directed')
            end if
            chg=v1-v0
        case('degcentralization')
            chg(1)=degree_centralization(g1,diag,mode)-degree_centralization(g0,diag,mode)
        case('betcentralization')
            chg(1)=bet_centralization(g1,mode)-bet_centralization(g0,mode)
        case('clocentralization')
            chg(1)=clo_centralization(g1,mode)-clo_centralization(g0,mode)
        case default
            chg(1)=effect_value(g1,effect,diag,mode)-effect_value(g0,effect,diag,mode)
        end select
    end subroutine pstar_effect_change

    real(dp) function degree_centralization(g,diag,mode) result(v)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::diag
        character(len=*),intent(in)::mode
        real(dp),allocatable::sc(:)
        real(dp)::tm
        integer::n
        n=size(g,1)
        if(trim(mode)=='graph')then
        sc=degree(g,'indegree',diag,.false.)
        tm=real((n-1)*(n-2+merge(1,0,diag)),dp)
        else
        sc=degree(g,'freeman',diag,.false.)
        tm=real((n-1)*(2*(n-1)-2+merge(1,0,diag)),dp)
        end if
        v=centralization_from_scores(sc,tm)
    end function degree_centralization

    real(dp) function bet_centralization(g,mode) result(v)
        real(dp),intent(in)::g(:,:)
        character(len=*),intent(in)::mode
        real(dp),allocatable::sc(:)
        real(dp)::tm
        integer::n
        n=size(g,1)
        if(trim(mode)=='graph')then
        sc=betweenness(g,'undirected')
        tm=real((n-1)*(n-1)*(n-2),dp)/2.0_dp
        else
        sc=betweenness(g,'directed')
        tm=real((n-1)*(n-1)*(n-2),dp)
        end if
        v=centralization_from_scores(sc,tm)
    end function bet_centralization

    real(dp) function clo_centralization(g,mode) result(v)
        real(dp),intent(in)::g(:,:)
        character(len=*),intent(in)::mode
        real(dp),allocatable::sc(:)
        real(dp)::tm
        integer::n
        n=size(g,1)
        if(trim(mode)=='graph')then
        sc=closeness(g,'undirected')
        tm=real((n-2)*(n-1),dp)/real(max(1,2*n-3),dp)
        else
        sc=closeness(g,'directed')
        tm=real(n-1,dp)*(1.0_dp-1.0_dp/real(max(1,n),dp))
        end if
        v=centralization_from_scores(sc,tm)
    end function clo_centralization

    function pstar_basic(dat,effects,diag,mode) result(fit)
        real(dp),intent(in)::dat(:,:)
        character(len=*),intent(in)::effects(:)
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        type(regression_result)::fit
        real(dp),allocatable::x(:,:),y(:),g0(:,:),g1(:,:),b(:),cov(:,:)
        integer::n,i,j,k,row,p,nobs,info
        real(dp)::ll
        logical::dg,undir
        character(len=12)::md
        n=size(dat,1)
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        nobs=0
        do i=1,n
        do j=1,n
        if(.not.dg.and.i==j)cycle
        if(undir.and.j>i)cycle
        if(.not.is_missing(dat(i,j)))nobs=nobs+1
        end do
        end do
        p=size(effects)
        allocate(x(nobs,p),y(nobs))
        row=0
        do i=1,n
        do j=1,n
            if(.not.dg.and.i==j)cycle
            if(undir.and.j>i)cycle
            if(is_missing(dat(i,j)))cycle
            row=row+1
            y(row)=merge(1.0_dp,0.0_dp,dat(i,j)>0.0_dp)
            allocate(g0(n,n),g1(n,n))
            g0=dat
            g1=dat
            g0(i,j)=0
            g1(i,j)=1
            if(undir)then
            g0(j,i)=0
            g1(j,i)=1
            end if
            do k=1,p
            x(row,k)=effect_value(g1,trim(effects(k)),dg,md)-effect_value(g0,trim(effects(k)),dg,md)
            end do
            deallocate(g0,g1)
        end do
        end do
        allocate(b(p),cov(p,p))
        call logistic_irls(x,y,b,cov,ll,info)
        call fill_regression(fit,x,y,b,cov,0.0_dp,info,.true.)
        fit%loglik=ll
    end function pstar_basic

    real(dp) function effect_value(g,effect,diag,mode) result(v)
        real(dp),intent(in)::g(:,:)
        character(len=*),intent(in)::effect,mode
        logical,intent(in)::diag
        select case(trim(effect))
        case('choice')
        v=sum(g)
        case('mutuality')
        v=mutuality(g)
        case('density')
        v=gden(g,diag,mode,.true.)
        case('reciprocity')
        v=grecip(g,'edgewise')
        case('stransitivity')
        v=gtrans(g,'strong')
        case('wtransitivity')
        v=gtrans(g,'weak')
        case('stranstri')
        v=gtrans(g,'strongcensus')
        case('wtranstri')
        v=gtrans(g,'weakcensus')
        case('connectedness')
        v=connectedness(g)
        case('hierarchy')
        v=hierarchy(g,'reciprocity')
        case('lubness')
        v=lubness(g)
        case('efficiency')
        v=efficiency(g,diag)
        case default
        v=0.0_dp
        end select
    end function effect_value

    subroutine graph_design(y,x,intercept,mode,diag,yv,xv)
        real(dp),intent(in)::y(:,:),x(:,:,:)
        logical,intent(in)::intercept,diag
        character(len=*),intent(in)::mode
        real(dp),allocatable,intent(out)::yv(:),xv(:,:)
        real(dp),allocatable::yt(:),xt(:,:),v(:)
        logical,allocatable::ok(:)
        integer::nv,p,k,i,c
        yt=gvectorize(y,mode,diag,.true.)
        nv=size(yt)
        p=size(x,1)+merge(1,0,intercept)
        allocate(xt(nv,p))
        c=0
        if(intercept)then
        c=1
        xt(:,1)=1.0_dp
        end if
        do k=1,size(x,1)
        v=gvectorize(x(k,:,:),mode,diag,.true.)
        xt(:,c+k)=v
        end do
        allocate(ok(nv))
        ok=.true.
        do i=1,nv
        if(is_missing(yt(i)).or.any(is_missing(xt(i,:))))ok(i)=.false.
        end do
        allocate(yv(count(ok)),xv(count(ok),p))
        c=0
        do i=1,nv
        if(ok(i))then
        c=c+1
        yv(c)=yt(i)
        xv(c,:)=xt(i,:)
        end if
        end do
    end subroutine graph_design

    subroutine fill_regression(fit,x,y,b,cov,sigma2,info,is_logit)
        type(regression_result),intent(out)::fit
        real(dp),intent(in)::x(:,:),y(:),b(:),cov(:,:),sigma2
        integer,intent(in)::info
        logical,intent(in)::is_logit
        integer::i,p,n
        real(dp),allocatable::eta(:)
        n=size(y)
        p=size(b)
        allocate(fit%coef(p),fit%se(p),fit%statistic(p),fit%fitted(n),fit%residual(n))
        fit%coef=b
        do i=1,p
        fit%se(i)=sqrt(max(cov(i,i),0.0_dp))
        fit%statistic(i)=merge(b(i)/fit%se(i),0.0_dp,fit%se(i)>0)
        end do
        eta=matmul(x,b)
        if(is_logit)then
        fit%fitted=1.0_dp/(1.0_dp+exp(-max(min(eta,700.0_dp),-700.0_dp)))
        else
        fit%fitted=eta
        end if
        fit%residual=y-fit%fitted
        fit%sigma2=sigma2
        fit%nobs=n
        fit%rank=p
        fit%converged=info==0
    end subroutine fill_regression

    real(dp) function eval_edgeperturbation(dat,i,j,fun) result(delta)
        real(dp), intent(in) :: dat(:,:)
        integer, intent(in) :: i,j
        procedure(scalar_graph_function) :: fun
        real(dp), allocatable :: present_graph(:,:), absent_graph(:,:)
        present_graph=dat
        absent_graph=dat
        present_graph(i,j)=1.0_dp
        absent_graph(i,j)=0.0_dp
        delta=fun(present_graph)-fun(absent_graph)
    end function eval_edgeperturbation

    function npostpred_scalar(b,fun) result(draws)
        type(bbnam_result), intent(in) :: b
        procedure(scalar_graph_function) :: fun
        real(dp), allocatable :: draws(:)
        integer :: s
        allocate(draws(size(b%net,1)))
        do s=1,size(b%net,1)
            draws(s)=fun(b%net(s,:,:))
        end do
    end function npostpred_scalar

    subroutine lnam_information(y,x,w1,w2,beta,rho1,rho2,sigmasq,acvm,se,info)
        real(dp), intent(in) :: y(:), beta(:), rho1(:), rho2(:), sigmasq
        real(dp), intent(in), optional :: x(:,:), w1(:,:,:), w2(:,:,:)
        real(dp), allocatable, intent(out) :: acvm(:,:), se(:)
        integer, intent(out) :: info
        real(dp), allocatable :: par(:), h(:), hess(:,:), invh(:,:), pp(:), pm(:), mp(:), mm(:)
        real(dp) :: f0, fp, fm, fpp, fpm, fmp, fmm
        integer :: m, nx, nw1, nw2, i, j

        nx=size(beta)
        nw1=size(rho1)
        nw2=size(rho2)
        m=nx+nw1+nw2+1
        allocate(par(m),h(m),hess(m,m),invh(m,m),pp(m),pm(m),mp(m),mm(m),se(m))
        if(nx>0)par(1:nx)=beta
        if(nw1>0)par(nx+1:nx+nw1)=rho1
        if(nw2>0)par(nx+nw1+1:nx+nw1+nw2)=rho2
        par(m)=sigmasq
        h=1.0e-4_dp*(1.0_dp+abs(par))
        h(m)=min(h(m),max(1.0e-8_dp,0.25_dp*max(sigmasq,1.0e-8_dp)))
        f0=lnam_full_nll(par,y,x,w1,w2,nx,nw1,nw2)
        hess=0.0_dp
        do i=1,m
            pp=par
            pm=par
            pp(i)=pp(i)+h(i)
            pm(i)=pm(i)-h(i)
            fp=lnam_full_nll(pp,y,x,w1,w2,nx,nw1,nw2)
            fm=lnam_full_nll(pm,y,x,w1,w2,nx,nw1,nw2)
            hess(i,i)=(fp-2.0_dp*f0+fm)/(h(i)*h(i))
            do j=i+1,m
                pp=par
                pm=par
                mp=par
                mm=par
                pp(i)=pp(i)+h(i)
                pp(j)=pp(j)+h(j)
                pm(i)=pm(i)+h(i)
                pm(j)=pm(j)-h(j)
                mp(i)=mp(i)-h(i)
                mp(j)=mp(j)+h(j)
                mm(i)=mm(i)-h(i)
                mm(j)=mm(j)-h(j)
                fpp=lnam_full_nll(pp,y,x,w1,w2,nx,nw1,nw2)
                fpm=lnam_full_nll(pm,y,x,w1,w2,nx,nw1,nw2)
                fmp=lnam_full_nll(mp,y,x,w1,w2,nx,nw1,nw2)
                fmm=lnam_full_nll(mm,y,x,w1,w2,nx,nw1,nw2)
                hess(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*h(i)*h(j))
                hess(j,i)=hess(i,j)
            end do
        end do
        call inverse_matrix(hess,invh,info)
        allocate(acvm(m,m))
        if(info==0)then
            acvm=invh
            do i=1,m
                se(i)=sqrt(max(acvm(i,i),0.0_dp))
            end do
        else
            acvm=0.0_dp
            se=sna_nan()
        end if
    end subroutine lnam_information

    real(dp) function lnam_full_nll(par,y,x,w1,w2,nx,nw1,nw2) result(v)
        real(dp), intent(in) :: par(:), y(:)
        real(dp), intent(in), optional :: x(:,:), w1(:,:,:), w2(:,:,:)
        integer, intent(in) :: nx, nw1, nw2
        real(dp), allocatable :: beta(:), rho1(:), rho2(:)
        integer :: m
        m=nx+nw1+nw2+1
        allocate(beta(nx),rho1(nw1),rho2(nw2))
        if(nx>0)beta=par(1:nx)
        if(nw1>0)rho1=par(nx+1:nx+nw1)
        if(nw2>0)rho2=par(nx+nw1+1:nx+nw1+nw2)
        v=0.5_dp*lnam_deviance(y,x,w1,w2,beta,rho1,rho2,max(par(m),tiny(1.0_dp)))
    end function lnam_full_nll

    subroutine lnam_mats(n,w1,w2,rho1,rho2,a1,a2)
        integer,intent(in)::n
        real(dp),intent(in),optional::w1(:,:,:),w2(:,:,:)
        real(dp),intent(in)::rho1(:),rho2(:)
        real(dp),allocatable,intent(out)::a1(:,:),a2(:,:)
        integer::i,k
        allocate(a1(n,n),a2(n,n))
        a1=0
        a2=0
        do i=1,n
        a1(i,i)=1
        a2(i,i)=1
        end do
        if(present(w1)) then
            do k=1,size(rho1)
                a1=a1-rho1(k)*w1(k,:,:)
            end do
        end if
        if(present(w2)) then
            do k=1,size(rho2)
                a2=a2-rho2(k)*w2(k,:,:)
            end do
        end if
    end subroutine lnam_mats

    real(dp) function lnam_deviance(y,x,w1,w2,beta,rho1,rho2,sigmasq) result(dev)
        real(dp),intent(in)::y(:),beta(:),rho1(:),rho2(:),sigmasq
        real(dp),intent(in),optional::x(:,:),w1(:,:,:),w2(:,:,:)
        real(dp),allocatable::a1(:,:),a2(:,:),z(:)
        real(dp)::ld1,ld2
        integer::n
        n=size(y)
        call lnam_mats(n,w1,w2,rho1,rho2,a1,a2)
        z=matmul(a1,y)
        if(present(x).and.size(beta)>0)z=z-matmul(x,beta)
        if(present(w2))z=matmul(a2,z)
        ld1=log(max(abs(determinant(a1)),tiny(1.0_dp)))
        ld2=log(max(abs(determinant(a2)),tiny(1.0_dp)))
        dev=real(n,dp)*(log(2.0_dp*acos(-1.0_dp))+log(max(sigmasq,tiny(1.0_dp))))+sum(z*z)/max(sigmasq,tiny(1.0_dp))-2*(ld1+ld2)
    end function lnam_deviance

    real(dp) function determinant(a) result(det)
        real(dp),intent(in)::a(:,:)
        real(dp),allocatable::b(:,:)
        real(dp)::mx,fac,t
        integer::n,i,j,k,p,sgn
        n=size(a,1)
        allocate(b(n,n))
        b=a
        det=1
        sgn=1
        do k=1,n
            p=k
            mx=abs(b(k,k))
            do i=k+1,n
            if(abs(b(i,k))>mx)then
            mx=abs(b(i,k))
            p=i
            end if
            end do
            if(mx<=tiny(1.0_dp))then
            det=0
            return
            end if
            if(p/=k)then
            do j=1,n
            t=b(k,j)
            b(k,j)=b(p,j)
            b(p,j)=t
            end do
            sgn=-sgn
            end if
            det=det*b(k,k)
            do i=k+1,n
            fac=b(i,k)/b(k,k)
            b(i,k:n)=b(i,k:n)-fac*b(k,k:n)
            end do
        end do
        det=det*real(sgn,dp)
    end function determinant

    subroutine solve_fit(a,b,x)
        real(dp),intent(in)::a(:,:),b(:)
        real(dp),intent(out)::x(:)
        real(dp),allocatable::inv(:,:)
        integer::info
        allocate(inv(size(a,1),size(a,2)))
        call inverse_matrix(a,inv,info)
        if(info==0)then
        x=matmul(inv,b)
        else
        x=b
        end if
    end subroutine solve_fit

    real(dp) function beta_rng(a,b) result(x)
        real(dp),intent(in)::a,b
        real(dp)::ga,gb
        ga=gamma_rng(max(a,1e-12_dp))
        gb=gamma_rng(max(b,1e-12_dp))
        x=ga/(ga+gb)
    end function beta_rng

    recursive real(dp) function gamma_rng(shape) result(x)
        real(dp),intent(in)::shape
        real(dp)::d,c,u,v,z
        if(shape<1.0_dp)then
            call random_number(u)
            x=gamma_rng(shape+1.0_dp)*u**(1.0_dp/shape)
            return
        end if
        d=shape-1.0_dp/3.0_dp
        c=1.0_dp/sqrt(9.0_dp*d)
        do
            z=normal_rng()
            v=(1.0_dp+c*z)**3
            if(v<=0)cycle
            call random_number(u)
            if(u<1.0_dp-0.0331_dp*z**4)exit
            if(log(u)<0.5_dp*z*z+d*(1.0_dp-v+log(v)))exit
        end do
        x=d*v
    end function gamma_rng

    real(dp) function normal_rng() result(z)
        real(dp)::u1,u2
        call random_number(u1)
        call random_number(u2)
        u1=max(u1,tiny(1.0_dp))
        z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
    end function normal_rng

end module sna_models
