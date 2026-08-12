! SPDX-License-Identifier: GPL-2.0-or-later
! HMM emission distributions translated from src/hmm.c and src/hmmderiv.c.
module msm_emissions
    use msm_kinds, only : dp
    use msm_stats, only : normal_pdf, lognormal_pdf, exponential_pdf, gamma_pdf, weibull_pdf, &
        poisson_pmf, binomial_pmf, beta_binomial_pmf, negbinomial_pmf, beta_pdf, student_t_pdf, &
        truncated_normal_pdf, me_truncated_normal_pdf, me_uniform_pdf, digamma_approx, &
        rand_uniform, rand_normal, rand_exponential, rand_gamma
    implicit none
    private
    integer, parameter, public :: hmm_cat=1, hmm_ident=2, hmm_unif=3, hmm_norm=4, hmm_lnorm=5
    integer, parameter, public :: hmm_exp=6, hmm_gamma=7, hmm_weibull=8, hmm_pois=9, hmm_binom=10
    integer, parameter, public :: hmm_beta_binom=11, hmm_tnorm=12, hmm_metnorm=13, hmm_meunif=14
    integer, parameter, public :: hmm_nbinom=15, hmm_beta=16, hmm_t=17
    type, public :: emission_model
        integer :: kind = hmm_ident
        real(dp), allocatable :: pars(:)
    end type emission_model
    public :: emission_probability, emission_derivative, emission_parameter_count, simulate_emission
contains
    pure function emission_parameter_count(model) result(n)
        type(emission_model), intent(in) :: model
        integer :: n
        if (allocated(model%pars)) then
            n=size(model%pars)
        else
            n=0
        end if
    end function emission_parameter_count

    pure function emission_probability(model,x) result(f)
        type(emission_model), intent(in) :: model
        real(dp), intent(in) :: x
        real(dp) :: f
        integer :: k
        if (.not. allocated(model%pars)) then
            f=0.0_dp; return
        end if
        select case(model%kind)
        case(hmm_cat)
            k=nint(x)
            if(k>=1 .and. k<=size(model%pars)) then; f=model%pars(k); else; f=0.0_dp; end if
        case(hmm_ident)
            f=merge(1.0_dp,0.0_dp,abs(x-model%pars(1)) <= epsilon(1.0_dp)*max(1.0_dp,abs(x)))
        case(hmm_unif)
            if(x>=model%pars(1) .and. x<=model%pars(2) .and. model%pars(2)>model%pars(1)) then
                f=1.0_dp/(model%pars(2)-model%pars(1))
            else; f=0.0_dp; end if
        case(hmm_norm)
            f=normal_pdf(x,model%pars(1),model%pars(2))
        case(hmm_lnorm)
            f=lognormal_pdf(x,model%pars(1),model%pars(2))
        case(hmm_exp)
            f=exponential_pdf(x,model%pars(1))
        case(hmm_gamma)
            f=gamma_pdf(x,model%pars(1),model%pars(2))
        case(hmm_weibull)
            f=weibull_pdf(x,model%pars(1),model%pars(2))
        case(hmm_pois)
            f=poisson_pmf(x,model%pars(1))
        case(hmm_binom)
            f=binomial_pmf(x,model%pars(1),model%pars(2))
        case(hmm_beta_binom)
            f=beta_binomial_pmf(x,model%pars(1),model%pars(2),model%pars(3))
        case(hmm_tnorm)
            f=truncated_normal_pdf(x,model%pars(1),model%pars(2),model%pars(3),model%pars(4))
        case(hmm_metnorm)
            f=me_truncated_normal_pdf(x,model%pars(1),model%pars(2),model%pars(3),model%pars(4), &
                model%pars(5),model%pars(6))
        case(hmm_meunif)
            f=me_uniform_pdf(x,model%pars(1),model%pars(2),model%pars(3),model%pars(4))
        case(hmm_nbinom)
            f=negbinomial_pmf(x,model%pars(1),model%pars(2))
        case(hmm_beta)
            f=beta_pdf(x,model%pars(1),model%pars(2))
        case(hmm_t)
            f=student_t_pdf(x,model%pars(1),model%pars(2),model%pars(3))
        case default
            f=0.0_dp
        end select
    end function emission_probability

    function emission_derivative(model,x) result(d)
        ! Derivatives with respect to model%pars in the same parameterization as msm.
        ! As in msm, derivatives for truncated-normal measurement-error families are unsupported.
        type(emission_model), intent(in) :: model
        real(dp), intent(in) :: x
        real(dp), allocatable :: d(:)
        real(dp) :: f, mean, sd, shape, rate, scale, rp, lambda, sizep, prob
        real(dp) :: meanp,sdp,a,b,pda,pdb,j11,j12,j21,j22,xmsq,df
        integer :: n,k
        n=emission_parameter_count(model); allocate(d(n)); d=0.0_dp
        if(n==0) return
        f=emission_probability(model,x)
        select case(model%kind)
        case(hmm_cat)
            k=nint(x); if(k>=1 .and. k<=n) d(k)=1.0_dp
        case(hmm_ident,hmm_unif,hmm_tnorm,hmm_metnorm,hmm_meunif)
            d=0.0_dp
        case(hmm_norm)
            mean=model%pars(1); sd=model%pars(2)
            d(1)=f*(x-mean)/(sd*sd)
            d(2)=f*(((x-mean)/sd)**2-1.0_dp)/sd
        case(hmm_lnorm)
            if(x<=0.0_dp) return
            mean=model%pars(1); sd=model%pars(2)
            d(1)=f*(log(x)-mean)/(sd*sd)
            d(2)=f*(((log(x)-mean)/sd)**2-1.0_dp)/sd
        case(hmm_exp)
            rate=model%pars(1); d(1)=(1.0_dp-rate*x)*exp(-rate*x)
        case(hmm_gamma)
            shape=model%pars(1); rate=model%pars(2)
            if(x>0.0_dp) then
                d(1)=f*(log(rate)+log(x)-digamma_approx(shape))
                d(2)=f*(shape/rate-x)
            end if
        case(hmm_weibull)
            shape=model%pars(1); scale=model%pars(2)
            if(x>0.0_dp) then
                rp=(x/scale)**shape
                d(1)=f*(1.0_dp/shape+log(x/scale)*(1.0_dp-rp))
                d(2)=f*(shape/scale*(rp-1.0_dp))
            end if
        case(hmm_pois)
            lambda=model%pars(1); if(lambda>0.0_dp) d(1)=(x/lambda-1.0_dp)*f
        case(hmm_binom)
            sizep=model%pars(1); prob=model%pars(2); d(1)=0.0_dp
            if(prob>0.0_dp .and. prob<1.0_dp) d(2)=f*(x/prob-(sizep-x)/(1.0_dp-prob))
        case(hmm_beta_binom)
            sizep=model%pars(1); meanp=model%pars(2); sdp=model%pars(3)
            if(x>=0.0_dp .and. x<=sizep .and. sdp>0.0_dp) then
                a=meanp/sdp; b=(1.0_dp-meanp)/sdp
                j11=1.0_dp/sdp; j12=-meanp/(sdp*sdp)
                j21=-1.0_dp/sdp; j22=-(1.0_dp-meanp)/(sdp*sdp)
                pda=f*(digamma_approx(x+a)-digamma_approx(sizep+a+b)-digamma_approx(a)+digamma_approx(a+b))
                pdb=f*(digamma_approx(sizep-x+b)-digamma_approx(sizep+a+b)-digamma_approx(b)+digamma_approx(a+b))
                d(1)=0.0_dp; d(2)=pda*j11+pdb*j21; d(3)=pda*j12+pdb*j22
            end if
        case(hmm_nbinom)
            sizep=model%pars(1); prob=model%pars(2)
            if(prob>0.0_dp .and. prob<1.0_dp) then
                d(1)=f*(digamma_approx(x+sizep)-digamma_approx(sizep)+log(prob))
                d(2)=f*(sizep/prob-x/(1.0_dp-prob))
            end if
        case(hmm_beta)
            a=model%pars(1); b=model%pars(2)
            if(x>0.0_dp .and. x<1.0_dp) then
                d(1)=f*(digamma_approx(a+b)-digamma_approx(a)+log(x))
                d(2)=f*(digamma_approx(a+b)-digamma_approx(b)+log(1.0_dp-x))
            end if
        case(hmm_t)
            mean=model%pars(1); scale=model%pars(2); df=model%pars(3); xmsq=(x-mean)**2
            d(1)=f*(x-mean)*(df+1.0_dp)/(df*scale*scale+xmsq)
            d(2)=f*(-1.0_dp/scale+(df+1.0_dp)*xmsq/(df*scale**3+scale*xmsq))
            d(3)=0.5_dp*f*(digamma_approx((df+1.0_dp)/2.0_dp)-digamma_approx(df/2.0_dp)-1.0_dp/df - &
                log(1.0_dp+xmsq/(df*scale*scale))+(df+1.0_dp)*xmsq/((df*scale)**2+df*xmsq))
        end select
    end function emission_derivative

    function simulate_emission(model) result(x)
        type(emission_model), intent(in) :: model
        real(dp) :: x, u, cum, p, a, b, z
        integer :: k, n, count
        select case(model%kind)
        case(hmm_cat)
            u=rand_uniform(); cum=0.0_dp; x=real(size(model%pars),dp)
            do k=1,size(model%pars)
                cum=cum+model%pars(k); if(u<=cum) then; x=real(k,dp); exit; end if
            end do
        case(hmm_ident)
            x=model%pars(1)
        case(hmm_unif)
            x=model%pars(1)+(model%pars(2)-model%pars(1))*rand_uniform()
        case(hmm_norm)
            x=model%pars(1)+model%pars(2)*rand_normal()
        case(hmm_lnorm)
            x=exp(model%pars(1)+model%pars(2)*rand_normal())
        case(hmm_exp)
            x=rand_exponential(model%pars(1))
        case(hmm_gamma)
            x=rand_gamma(model%pars(1),model%pars(2))
        case(hmm_weibull)
            x=model%pars(2)*(-log(rand_uniform()))**(1.0_dp/model%pars(1))
        case(hmm_pois)
            p=exp(-model%pars(1)); cum=p; u=rand_uniform(); k=0
            do while(u>cum)
                k=k+1; p=p*model%pars(1)/real(k,dp); cum=cum+p
            end do
            x=real(k,dp)
        case(hmm_binom)
            n=nint(model%pars(1)); count=0
            do k=1,n; if(rand_uniform()<model%pars(2)) count=count+1; end do
            x=real(count,dp)
        case(hmm_beta_binom)
            a=model%pars(2)/model%pars(3); b=(1.0_dp-model%pars(2))/model%pars(3)
            p=rand_gamma(a,1.0_dp); z=rand_gamma(b,1.0_dp); p=p/(p+z); n=nint(model%pars(1)); count=0
            do k=1,n; if(rand_uniform()<p) count=count+1; end do; x=real(count,dp)
        case(hmm_nbinom)
            ! Poisson-gamma mixture with R's size/prob parameterization.
            p=rand_gamma(model%pars(1),model%pars(2)/(1.0_dp-model%pars(2)))
            u=exp(-p); cum=u; z=rand_uniform(); k=0
            do while(z>cum); k=k+1; u=u*p/real(k,dp); cum=cum+u; end do; x=real(k,dp)
        case(hmm_beta)
            a=rand_gamma(model%pars(1),1.0_dp); b=rand_gamma(model%pars(2),1.0_dp); x=a/(a+b)
        case(hmm_t)
            z=rand_normal(); p=rand_gamma(0.5_dp*model%pars(3),0.5_dp)
            x=model%pars(1)+model%pars(2)*z/sqrt(p/model%pars(3))
        case(hmm_tnorm)
            do; x=model%pars(1)+model%pars(2)*rand_normal(); if(x>=model%pars(3).and.x<=model%pars(4)) exit; end do
        case default
            x=0.0_dp
        end select
    end function simulate_emission
end module msm_emissions
