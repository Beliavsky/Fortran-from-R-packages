! SPDX-License-Identifier: LGPL-3.0-or-later
! Translation of the computational distributions in the R package sadists
! by Steven E. Pav (LGPL-3.0-or-later).
module sadists_distributions
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use sadists_kinds, only : dp
    use sadists_special, only : nan_dp, pos_inf_dp, neg_inf_dp, normal_moments, &
        chisq_log_moment, moments_to_cumulants, cumulants_to_moments, &
        lognoncentral_chisq_cumulants, random_normal, random_chisq
    use sadists_approximations, only : edgeworth_pdf, edgeworth_cdf, cornish_fisher_quantile
    implicit none
    private

    public :: ddnf, pdnf, qdnf, rdnf, ddnt, pdnt, qdnt, rdnt
    public :: ddnbeta, pdnbeta, qdnbeta, rdnbeta
    public :: ddneta, pdneta, qdneta, rdneta
    public :: dlambdap, plambdap, qlambdap, rlambdap
    public :: dkprime, pkprime, qkprime, rkprime
    public :: dupsilon, pupsilon, qupsilon, rupsilon
    public :: dsumchisqpow, psumchisqpow, qsumchisqpow, rsumchisqpow
    public :: dsumlogchisq, psumlogchisq, qsumlogchisq, rsumlogchisq
    public :: dprodchisqpow, pprodchisqpow, qprodchisqpow, rprodchisqpow
    public :: dproddnf, pproddnf, qproddnf, rproddnf
    public :: dprodnormal, pprodnormal, qprodnormal, rprodnormal
    public :: ddnf_vec, pdnf_vec, qdnf_vec, ddnt_vec, pdnt_vec, qdnt_vec
    public :: ddnbeta_vec, pdnbeta_vec, qdnbeta_vec, ddneta_vec, pdneta_vec, qdneta_vec
    public :: norm_cumulants, chipow_cumulants, sumchisqpow_cumulants
    public :: dnf_moments, dnt_moments, lambdap_cumulants, kprime_cumulants
    public :: upsilon_cumulants, sumlogchisq_cumulants, prodnormal_cumulants

contains

    pure real(dp) function recycled(a,i) result(v)
        real(dp), intent(in) :: a(:)
        integer, intent(in) :: i
        v=a(mod(i-1,size(a))+1)
    end function recycled

    pure real(dp) function recycled_optional(a,i,default) result(v)
        real(dp), intent(in), optional :: a(:)
        integer, intent(in) :: i
        real(dp), intent(in) :: default
        if(present(a)) then
            v=recycled(a,i)
        else
            v=default
        end if
    end function recycled_optional



    pure integer function vector_terms(wts,df,ncp,pow) result(n)
        real(dp), intent(in) :: wts(:),df(:)
        real(dp), intent(in), optional :: ncp(:),pow(:)
        n=max(size(wts),size(df))
        if(present(ncp)) n=max(n,size(ncp))
        if(present(pow)) n=max(n,size(pow))
    end function vector_terms

    pure integer function get_order(order_max,default_order) result(n)
        integer, intent(in), optional :: order_max
        integer, intent(in) :: default_order
        n=default_order
        if(present(order_max)) n=order_max
        n=max(2,n)
    end function get_order

    pure logical function get_log(flag) result(v)
        logical, intent(in), optional :: flag
        v=.false.; if(present(flag)) v=flag
    end function get_log

    pure logical function get_lower(flag) result(v)
        logical, intent(in), optional :: flag
        v=.true.; if(present(flag)) v=flag
    end function get_lower

    pure subroutine norm_cumulants(mean,sd,kappa)
        real(dp),intent(in)::mean,sd
        real(dp),intent(out)::kappa(:)
        kappa=0.0_dp
        if(size(kappa)>=1) kappa(1)=mean
        if(size(kappa)>=2) kappa(2)=sd*sd
    end subroutine norm_cumulants

    pure subroutine chipow_moments(df,ncp,pow,moms)
        real(dp),intent(in)::df,ncp,pow
        real(dp),intent(out)::moms(:)
        integer::r
        do r=1,size(moms)
            moms(r)=exp(chisq_log_moment(df,ncp,pow*real(r,dp)))
        end do
    end subroutine chipow_moments

    pure subroutine chipow_cumulants(df,ncp,pow,kappa)
        real(dp),intent(in)::df,ncp,pow
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::moms(:)
        integer::r
        if(pow==1.0_dp) then
            do r=1,size(kappa)
                kappa(r)=2.0_dp**(r-1)*gamma(real(r,dp))*(df+real(r,dp)*ncp)
            end do
        else
            allocate(moms(size(kappa)))
            call chipow_moments(df,ncp,pow,moms)
            call moments_to_cumulants(moms,kappa)
        end if
    end subroutine chipow_cumulants

    pure subroutine sumchisqpow_cumulants(wts,df,ncp,pow,kappa)
        real(dp),intent(in)::wts(:),df(:)
        real(dp),intent(in),optional::ncp(:),pow(:)
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::sub(:)
        real(dp)::w,dd,nc,pw
        integer::i,r,n
        n=vector_terms(wts,df,ncp,pow)
        allocate(sub(size(kappa))); kappa=0.0_dp
        do i=1,n
            w=recycled(wts,i); dd=recycled(df,i)
            nc=recycled_optional(ncp,i,0.0_dp); pw=recycled_optional(pow,i,1.0_dp)
            call chipow_cumulants(dd,nc,pw,sub)
            do r=1,size(kappa)
                kappa(r)=kappa(r)+w**r*sub(r)
            end do
        end do
    end subroutine sumchisqpow_cumulants

    pure subroutine sumchisqpow_support(wts,lo,hi)
        real(dp),intent(in)::wts(:)
        real(dp),intent(out)::lo,hi
        if(minval(wts)<0.0_dp) then; lo=neg_inf_dp(); else; lo=0.0_dp; end if
        if(maxval(wts)>0.0_dp) then; hi=pos_inf_dp(); else; hi=0.0_dp; end if
    end subroutine sumchisqpow_support

    pure subroutine dnf_moments(df1,df2,ncp1,ncp2,moms)
        real(dp),intent(in)::df1,df2,ncp1,ncp2
        real(dp),intent(out)::moms(:)
        integer::r
        real(dp)::rr,lm
        do r=1,size(moms)
            rr=real(r,dp)
            lm=chisq_log_moment(df1,ncp1,rr)-rr*log(df1) &
                +chisq_log_moment(df2,ncp2,-rr)+rr*log(df2)
            moms(r)=exp(lm)
        end do
    end subroutine dnf_moments

    pure subroutine dnf_cumulants(df1,df2,ncp1,ncp2,kappa)
        real(dp),intent(in)::df1,df2,ncp1,ncp2
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::moms(:)
        allocate(moms(size(kappa)))
        call dnf_moments(df1,df2,ncp1,ncp2,moms)
        call moments_to_cumulants(moms,kappa)
    end subroutine dnf_cumulants

    pure subroutine dnt_moments(df,ncp1,ncp2,moms)
        real(dp),intent(in)::df,ncp1,ncp2
        real(dp),intent(out)::moms(:)
        real(dp),allocatable::nm(:)
        real(dp)::rr
        integer::r
        allocate(nm(size(moms))); call normal_moments(ncp1,1.0_dp,nm)
        do r=1,size(moms)
            rr=real(r,dp)
            moms(r)=nm(r)*exp(chisq_log_moment(df,ncp2,-0.5_dp*rr)+0.5_dp*rr*log(df))
        end do
    end subroutine dnt_moments

    pure subroutine dnt_cumulants(df,ncp1,ncp2,kappa)
        real(dp),intent(in)::df,ncp1,ncp2
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::moms(:)
        allocate(moms(size(kappa))); call dnt_moments(df,ncp1,ncp2,moms)
        call moments_to_cumulants(moms,kappa)
    end subroutine dnt_cumulants

    pure subroutine upsilon_cumulants(df,t,kappa)
        real(dp),intent(in)::df(:),t(:)
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::sub(:)
        real(dp)::dd,tt
        integer::i,r,n
        kappa=0.0_dp
        if(size(kappa)>=2) kappa(2)=1.0_dp
        allocate(sub(size(kappa)))
        n=max(size(df),size(t))
        do i=1,n
            dd=recycled(df,i); tt=recycled(t,i)
            if(.not.ieee_is_finite(dd)) then
                kappa(1)=kappa(1)+tt
            else
                call chipow_cumulants(dd,0.0_dp,0.5_dp,sub)
                do r=1,size(kappa)
                    kappa(r)=kappa(r)+(tt/sqrt(dd))**r*sub(r)
                end do
            end if
        end do
    end subroutine upsilon_cumulants

    pure subroutine lambdap_cumulants(df,t,kappa)
        real(dp),intent(in)::df,t
        real(dp),intent(out)::kappa(:)
        real(dp)::dfa(1),ta(1)
        dfa=[df]; ta=[t]
        call upsilon_cumulants(dfa,ta,kappa)
    end subroutine lambdap_cumulants

    pure subroutine lambdap_moments(df,t,moms)
        real(dp),intent(in)::df,t
        real(dp),intent(out)::moms(:)
        real(dp),allocatable::kap(:)
        allocate(kap(size(moms))); call lambdap_cumulants(df,t,kap)
        call cumulants_to_moments(kap,moms)
    end subroutine lambdap_moments

    pure subroutine kprime_cumulants(v1,v2,a,b,kappa)
        real(dp),intent(in)::v1,v2,a,b
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::moms(:),lm(:),nm(:)
        real(dp)::rr
        integer::r
        allocate(moms(size(kappa)))
        if(.not.ieee_is_finite(v1)) then
            allocate(nm(size(kappa))); call normal_moments(a,abs(b),nm)
            do r=1,size(kappa)
                rr=real(r,dp)
                if(ieee_is_finite(v2)) then
                    moms(r)=nm(r)*exp(chisq_log_moment(v2,0.0_dp,-0.5_dp*rr)+0.5_dp*rr*log(v2))
                else
                    moms(r)=nm(r)
                end if
            end do
        else if(b/=0.0_dp) then
            allocate(lm(size(kappa))); call lambdap_moments(v1,a/b,lm)
            do r=1,size(kappa)
                rr=real(r,dp)
                if(ieee_is_finite(v2)) then
                    moms(r)=(b*sqrt(v2))**r*lm(r)*exp(chisq_log_moment(v2,0.0_dp,-0.5_dp*rr))
                else
                    moms(r)=b**r*lm(r)
                end if
            end do
        else
            do r=1,size(kappa)
                rr=real(r,dp)
                moms(r)=(a/sqrt(v1))**r*exp(chisq_log_moment(v1,0.0_dp,0.5_dp*rr))
                if(ieee_is_finite(v2)) then
                    moms(r)=moms(r)*v2**(0.5_dp*rr)*exp(chisq_log_moment(v2,0.0_dp,-0.5_dp*rr))
                end if
            end do
        end if
        call moments_to_cumulants(moms,kappa)
    end subroutine kprime_cumulants

    pure subroutine sumlogchisq_cumulants(wts,df,ncp,kappa)
        real(dp),intent(in)::wts(:),df(:)
        real(dp),intent(in),optional::ncp(:)
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::sub(:)
        real(dp)::w,dd,nc
        integer::i,r,n
        n=max(size(wts),size(df)); if(present(ncp)) n=max(n,size(ncp))
        allocate(sub(size(kappa))); kappa=0.0_dp
        do i=1,n
            w=recycled(wts,i); dd=recycled(df,i); nc=recycled_optional(ncp,i,0.0_dp)
            call lognoncentral_chisq_cumulants(dd,nc,sub)
            do r=1,size(kappa)
                kappa(r)=kappa(r)+w**r*sub(r)
            end do
        end do
    end subroutine sumlogchisq_cumulants

    pure subroutine prodnormal_cumulants(mu,sigma,kappa)
        real(dp),intent(in)::mu(:),sigma(:)
        real(dp),intent(out)::kappa(:)
        real(dp),allocatable::sub(:),moms(:)
        integer::i,n
        n=max(size(mu),size(sigma)); allocate(sub(size(kappa)),moms(size(kappa))); moms=1.0_dp
        do i=1,n
            call normal_moments(recycled(mu,i),recycled(sigma,i),sub)
            moms=moms*sub
        end do
        call moments_to_cumulants(moms,kappa)
    end subroutine prodnormal_cumulants

    pure real(dp) function ddnf(x,df1,df2,ncp1,ncp2,log_density,order_max) result(v)
        real(dp),intent(in)::x,df1,df2,ncp1,ncp2
        logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call dnf_cumulants(df1,df2,ncp1,ncp2,k)
        v=edgeworth_pdf(x,k,0.0_dp,pos_inf_dp(),get_log(log_density))
    end function ddnf

    pure real(dp) function pdnf(q,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df1,df2,ncp1,ncp2
        logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call dnf_cumulants(df1,df2,ncp1,ncp2,k)
        v=edgeworth_cdf(q,k,0.0_dp,pos_inf_dp(),get_lower(lower_tail),get_log(log_p))
    end function pdnf

    pure real(dp) function qdnf(p,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df1,df2,ncp1,ncp2
        logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call dnf_cumulants(df1,df2,ncp1,ncp2,k)
        v=cornish_fisher_quantile(p,k,0.0_dp,pos_inf_dp(),get_lower(lower_tail),get_log(log_p))
    end function qdnf

    subroutine rdnf(x,df1,df2,ncp1,ncp2)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df1,df2,ncp1,ncp2
        integer::i
        do i=1,size(x); x(i)=(random_chisq(df1,ncp1)/df1)/(random_chisq(df2,ncp2)/df2); end do
    end subroutine rdnf

    pure real(dp) function ddnt(x,df,ncp1,ncp2,log_density,order_max) result(v)
        real(dp),intent(in)::x,df,ncp1,ncp2; logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call dnt_cumulants(df,ncp1,ncp2,k)
        v=edgeworth_pdf(x,k,log_density=get_log(log_density))
    end function ddnt

    pure real(dp) function pdnt(q,df,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df,ncp1,ncp2; logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call dnt_cumulants(df,ncp1,ncp2,k)
        v=edgeworth_cdf(q,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function pdnt

    pure real(dp) function qdnt(p,df,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df,ncp1,ncp2; logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call dnt_cumulants(df,ncp1,ncp2,k)
        v=cornish_fisher_quantile(p,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function qdnt

    subroutine rdnt(x,df,ncp1,ncp2)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df,ncp1,ncp2; integer::i
        do i=1,size(x); x(i)=random_normal(ncp1,1.0_dp)/sqrt(random_chisq(df,ncp2)/df); end do
    end subroutine rdnt

    pure real(dp) function ddnbeta(x,df1,df2,ncp1,ncp2,log_density,order_max) result(v)
        real(dp),intent(in)::x,df1,df2,ncp1,ncp2; logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max
        real(dp)::xf; logical::ld
        ld=get_log(log_density)
        if(x<0.0_dp .or. x>=1.0_dp) then; v=nan_dp(); return; end if
        xf=(df2/df1)*x/(1.0_dp-x)
        v=ddnf(xf,df1,df2,ncp1,ncp2,ld,order_max)
        if(ld) then; v=log(df2/df1)+v-2.0_dp*log(1.0_dp-x)
        else; v=(df2/df1)*v/(1.0_dp-x)**2; end if
    end function ddnbeta

    pure real(dp) function pdnbeta(q,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df1,df2,ncp1,ncp2; logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        logical::lower,lp; real(dp)::qf
        lower=get_lower(lower_tail); lp=get_log(log_p)
        if(q<0.0_dp) then; if(lower) then; v=0.0_dp; else; v=1.0_dp; end if
        else if(q>=1.0_dp) then; if(lower) then; v=1.0_dp; else; v=0.0_dp; end if
        else; qf=(df2/df1)*q/(1.0_dp-q); v=pdnf(qf,df1,df2,ncp1,ncp2,lower,lp,order_max); return; end if
        if(lp) then; if(v==0.0_dp) then; v=neg_inf_dp(); else; v=log(v); end if; end if
    end function pdnbeta

    pure real(dp) function qdnbeta(p,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df1,df2,ncp1,ncp2; logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp)::qf
        if((.not.get_log(log_p) .and. (p<0.0_dp .or. p>1.0_dp)) .or. &
           (get_log(log_p) .and. p>0.0_dp)) then; v=nan_dp(); return; end if
        qf=(df1/df2)*qdnf(p,df1,df2,ncp1,ncp2,get_lower(lower_tail),get_log(log_p),order_max)
        if (.not. ieee_is_finite(qf)) then
            if (qf > 0.0_dp) then; v=1.0_dp; else; v=nan_dp(); end if
        else
            v=qf/(1.0_dp+qf)
        end if
    end function qdnbeta

    subroutine rdnbeta(x,df1,df2,ncp1,ncp2)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df1,df2,ncp1,ncp2
        real(dp)::a,b; integer::i
        do i=1,size(x); a=random_chisq(df1,ncp1); b=random_chisq(df2,ncp2); x(i)=a/(a+b); end do
    end subroutine rdnbeta

    pure real(dp) function ddneta(x,df,ncp1,ncp2,log_density,order_max) result(v)
        real(dp),intent(in)::x,df,ncp1,ncp2; logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max
        real(dp)::xf; logical::ld
        ld=get_log(log_density)
        if(x< -1.0_dp .or. x>=1.0_dp) then; v=nan_dp(); return; end if
        xf=sqrt(df)*x/sqrt(1.0_dp-x*x)
        v=ddnt(xf,df,ncp1,ncp2,ld,order_max)
        if(ld) then; v=0.5_dp*log(df)+v-1.5_dp*log(1.0_dp-x*x)
        else; v=sqrt(df)*v/(1.0_dp-x*x)**1.5_dp; end if
    end function ddneta

    pure real(dp) function pdneta(q,df,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df,ncp1,ncp2; logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        logical::lower,lp; real(dp)::qf
        lower=get_lower(lower_tail); lp=get_log(log_p)
        if(q< -1.0_dp) then; if(lower) then; v=0.0_dp; else; v=1.0_dp; end if
        else if(q>=1.0_dp) then; if(lower) then; v=1.0_dp; else; v=0.0_dp; end if
        else; qf=sqrt(df)*q/sqrt(1.0_dp-q*q); v=pdnt(qf,df,ncp1,ncp2,lower,lp,order_max); return; end if
        if(lp) then; if(v==0.0_dp) then; v=neg_inf_dp(); else; v=log(v); end if; end if
    end function pdneta

    pure real(dp) function qdneta(p,df,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df,ncp1,ncp2; logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp)::qf
        qf=qdnt(p,df,ncp1,ncp2,get_lower(lower_tail),get_log(log_p),order_max)/sqrt(df)
        if (.not. ieee_is_finite(qf)) then
            if (qf > 0.0_dp) then; v=1.0_dp
            else if (qf < 0.0_dp) then; v=-1.0_dp
            else; v=nan_dp(); end if
        else
            v=qf/sqrt(1.0_dp+qf*qf)
        end if
    end function qdneta

    subroutine rdneta(x,df,ncp1,ncp2)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df,ncp1,ncp2
        real(dp)::z,y; integer::i
        do i=1,size(x); z=random_normal(ncp1,1.0_dp); y=random_chisq(df,ncp2); x(i)=z/sqrt(z*z+y); end do
    end subroutine rdneta

    pure real(dp) function dlambdap(x,df,t,log_density,order_max) result(v)
        real(dp),intent(in)::x,df,t; logical,intent(in),optional::log_density; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call lambdap_cumulants(df,t,k)
        v=edgeworth_pdf(x,k,log_density=get_log(log_density))
    end function dlambdap

    pure real(dp) function plambdap(q,df,t,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df,t; logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call lambdap_cumulants(df,t,k)
        v=edgeworth_cdf(q,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function plambdap

    pure real(dp) function qlambdap(p,df,t,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df,t; logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call lambdap_cumulants(df,t,k)
        v=cornish_fisher_quantile(p,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function qlambdap

    subroutine rlambdap(x,df,t)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df,t; real(dp)::s; integer::i
        do i=1,size(x)
            if(ieee_is_finite(df)) then; s=sqrt(random_chisq(df,0.0_dp)/df); else; s=1.0_dp; end if
            x(i)=random_normal(t*s,1.0_dp)
        end do
    end subroutine rlambdap

    pure real(dp) function dkprime(x,v1,v2,a,b,order_max,log_density) result(v)
        real(dp),intent(in)::x,v1,v2,a; real(dp),intent(in),optional::b
        integer,intent(in),optional::order_max; logical,intent(in),optional::log_density
        real(dp)::bb; real(dp),allocatable::k(:); integer::om
        bb=1.0_dp; if(present(b)) bb=b
        om=get_order(order_max,6); allocate(k(om)); call kprime_cumulants(v1,v2,a,bb,k)
        v=edgeworth_pdf(x,k,log_density=get_log(log_density))
    end function dkprime

    pure real(dp) function pkprime(q,v1,v2,a,b,order_max,lower_tail,log_p) result(v)
        real(dp),intent(in)::q,v1,v2,a; real(dp),intent(in),optional::b
        integer,intent(in),optional::order_max; logical,intent(in),optional::lower_tail,log_p
        real(dp)::bb; real(dp),allocatable::k(:); integer::om
        bb=1.0_dp; if(present(b)) bb=b
        om=get_order(order_max,6); allocate(k(om)); call kprime_cumulants(v1,v2,a,bb,k)
        v=edgeworth_cdf(q,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function pkprime

    pure real(dp) function qkprime(p,v1,v2,a,b,order_max,lower_tail,log_p) result(v)
        real(dp),intent(in)::p,v1,v2,a; real(dp),intent(in),optional::b
        integer,intent(in),optional::order_max; logical,intent(in),optional::lower_tail,log_p
        real(dp)::bb; real(dp),allocatable::k(:); integer::om
        bb=1.0_dp; if(present(b)) bb=b
        om=get_order(order_max,6); allocate(k(om)); call kprime_cumulants(v1,v2,a,bb,k)
        v=cornish_fisher_quantile(p,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function qkprime

    subroutine rkprime(x,v1,v2,a,b)
        real(dp),intent(out)::x(:); real(dp),intent(in)::v1,v2,a; real(dp),intent(in),optional::b
        real(dp)::bb,s1,s2; integer::i
        bb=1.0_dp; if(present(b)) bb=b
        do i=1,size(x)
            if(ieee_is_finite(v1)) then; s1=sqrt(random_chisq(v1,0.0_dp)/v1); else; s1=1.0_dp; end if
            if(ieee_is_finite(v2)) then; s2=sqrt(random_chisq(v2,0.0_dp)/v2); else; s2=1.0_dp; end if
            x(i)=(random_normal(0.0_dp,abs(bb))+a*s1)/s2
        end do
    end subroutine rkprime

    pure real(dp) function dupsilon(x,df,t,log_density,order_max) result(v)
        real(dp),intent(in)::x,df(:),t(:); logical,intent(in),optional::log_density; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call upsilon_cumulants(df,t,k)
        v=edgeworth_pdf(x,k,log_density=get_log(log_density))
    end function dupsilon

    pure real(dp) function pupsilon(q,df,t,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df(:),t(:); logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call upsilon_cumulants(df,t,k)
        v=edgeworth_cdf(q,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function pupsilon

    pure real(dp) function qupsilon(p,df,t,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df(:),t(:); logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call upsilon_cumulants(df,t,k)
        v=cornish_fisher_quantile(p,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function qupsilon

    subroutine rupsilon(x,df,t)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df(:),t(:)
        real(dp)::s,dd,tt; integer::i,j,n
        n=max(size(df),size(t))
        do i=1,size(x)
            x(i)=random_normal()
            do j=1,n
                dd=recycled(df,j); tt=recycled(t,j)
                if(ieee_is_finite(dd)) then; s=sqrt(random_chisq(dd,0.0_dp)/dd); else; s=1.0_dp; end if
                x(i)=x(i)+tt*s
            end do
        end do
    end subroutine rupsilon

    pure real(dp) function dsumchisqpow(x,wts,df,ncp,pow,log_density,order_max) result(v)
        real(dp),intent(in)::x,wts(:),df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        logical,intent(in),optional::log_density; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); real(dp)::lo,hi; integer::om
        om=get_order(order_max,6); allocate(k(om)); call sumchisqpow_cumulants(wts,df,ncp,pow,k)
        call sumchisqpow_support(wts,lo,hi); v=edgeworth_pdf(x,k,lo,hi,get_log(log_density))
    end function dsumchisqpow

    pure real(dp) function psumchisqpow(q,wts,df,ncp,pow,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,wts(:),df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); real(dp)::lo,hi; integer::om
        om=get_order(order_max,6); allocate(k(om)); call sumchisqpow_cumulants(wts,df,ncp,pow,k)
        call sumchisqpow_support(wts,lo,hi); v=edgeworth_cdf(q,k,lo,hi,get_lower(lower_tail),get_log(log_p))
    end function psumchisqpow

    pure real(dp) function qsumchisqpow(p,wts,df,ncp,pow,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,wts(:),df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); real(dp)::lo,hi; integer::om
        om=get_order(order_max,6); allocate(k(om)); call sumchisqpow_cumulants(wts,df,ncp,pow,k)
        call sumchisqpow_support(wts,lo,hi); v=cornish_fisher_quantile(p,k,lo,hi,get_lower(lower_tail),get_log(log_p))
    end function qsumchisqpow

    subroutine rsumchisqpow(x,wts,df,ncp,pow)
        real(dp),intent(out)::x(:); real(dp),intent(in)::wts(:),df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        real(dp)::w,dd,nc,pw; integer::i,j,n
        n=vector_terms(wts,df,ncp,pow); x=0.0_dp
        do j=1,n
            w=recycled(wts,j); dd=recycled(df,j); nc=recycled_optional(ncp,j,0.0_dp); pw=recycled_optional(pow,j,1.0_dp)
            do i=1,size(x); x(i)=x(i)+w*random_chisq(dd,nc)**pw; end do
        end do
    end subroutine rsumchisqpow

    pure real(dp) function dsumlogchisq(x,wts,df,ncp,log_density,order_max) result(v)
        real(dp),intent(in)::x,wts(:),df(:); real(dp),intent(in),optional::ncp(:)
        logical,intent(in),optional::log_density; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call sumlogchisq_cumulants(wts,df,ncp,k)
        v=edgeworth_pdf(x,k,log_density=get_log(log_density))
    end function dsumlogchisq

    pure real(dp) function psumlogchisq(q,wts,df,ncp,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,wts(:),df(:); real(dp),intent(in),optional::ncp(:)
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call sumlogchisq_cumulants(wts,df,ncp,k)
        v=edgeworth_cdf(q,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function psumlogchisq

    pure real(dp) function qsumlogchisq(p,wts,df,ncp,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,wts(:),df(:); real(dp),intent(in),optional::ncp(:)
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,6); allocate(k(om)); call sumlogchisq_cumulants(wts,df,ncp,k)
        v=cornish_fisher_quantile(p,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function qsumlogchisq

    subroutine rsumlogchisq(x,wts,df,ncp)
        real(dp),intent(out)::x(:); real(dp),intent(in)::wts(:),df(:); real(dp),intent(in),optional::ncp(:)
        real(dp)::w,dd,nc; integer::i,j,n
        n=max(size(wts),size(df)); if(present(ncp)) n=max(n,size(ncp)); x=0.0_dp
        do j=1,n
            w=recycled(wts,j); dd=recycled(df,j); nc=recycled_optional(ncp,j,0.0_dp)
            do i=1,size(x); x(i)=x(i)+w*log(random_chisq(dd,nc)); end do
        end do
    end subroutine rsumlogchisq

    pure real(dp) function dprodchisqpow(x,df,ncp,pow,log_density,order_max) result(v)
        real(dp),intent(in)::x,df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        logical,intent(in),optional::log_density; integer,intent(in),optional::order_max
        real(dp),allocatable::w(:); integer::n,i,om; logical::ld
        ld=get_log(log_density); om=get_order(order_max,5)
        if(x<=0.0_dp) then; if(ld) then; v=neg_inf_dp(); else; v=0.0_dp; end if; return; end if
        n=size(df); if(present(ncp)) n=max(n,size(ncp)); if(present(pow)) n=max(n,size(pow)); allocate(w(n))
        do i=1,n; w(i)=recycled_optional(pow,i,1.0_dp); end do
        v=dsumlogchisq(log(x),w,df,ncp,ld,om)
        if(ld) then; v=v-log(x); else; v=v/x; end if
    end function dprodchisqpow

    pure real(dp) function pprodchisqpow(q,df,ncp,pow,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::w(:); integer::n,i,om; logical::lower,lp
        lower=get_lower(lower_tail); lp=get_log(log_p); om=get_order(order_max,5)
        if(q<=0.0_dp) then; if(lower) then; v=0.0_dp; else; v=1.0_dp; end if
            if(lp) then; if(v==0.0_dp) then; v=neg_inf_dp(); else; v=0.0_dp; end if; end if; return; end if
        n=size(df); if(present(ncp)) n=max(n,size(ncp)); if(present(pow)) n=max(n,size(pow)); allocate(w(n))
        do i=1,n; w(i)=recycled_optional(pow,i,1.0_dp); end do
        v=psumlogchisq(log(q),w,df,ncp,lower,lp,om)
    end function pprodchisqpow

    pure real(dp) function qprodchisqpow(p,df,ncp,pow,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::w(:); integer::n,i,om
        om=get_order(order_max,5); n=size(df)
        if(present(ncp)) n=max(n,size(ncp))
        if(present(pow)) n=max(n,size(pow))
        allocate(w(n))
        do i=1,n; w(i)=recycled_optional(pow,i,1.0_dp); end do
        v=exp(qsumlogchisq(p,w,df,ncp,get_lower(lower_tail),get_log(log_p),om))
    end function qprodchisqpow

    subroutine rprodchisqpow(x,df,ncp,pow)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df(:); real(dp),intent(in),optional::ncp(:),pow(:)
        real(dp)::dd,nc,pw; integer::i,j,n
        n=size(df); if(present(ncp)) n=max(n,size(ncp)); if(present(pow)) n=max(n,size(pow)); x=1.0_dp
        do j=1,n
            dd=recycled(df,j); nc=recycled_optional(ncp,j,0.0_dp); pw=recycled_optional(pow,j,1.0_dp)
            do i=1,size(x); x(i)=x(i)*random_chisq(dd,nc)**pw; end do
        end do
    end subroutine rprodchisqpow

    pure subroutine proddnf_log_cumulants(df1,df2,ncp1,ncp2,kappa,shift)
        real(dp),intent(in)::df1(:),df2(:),ncp1(:),ncp2(:)
        real(dp),intent(out)::kappa(:),shift
        real(dp),allocatable::sub(:)
        real(dp)::d1,d2,n1,n2
        integer::i,r,n
        n=max(size(df1),size(df2),size(ncp1),size(ncp2)); allocate(sub(size(kappa))); kappa=0.0_dp; shift=0.0_dp
        do i=1,n
            d1=recycled(df1,i); d2=recycled(df2,i); n1=recycled(ncp1,i); n2=recycled(ncp2,i)
            call lognoncentral_chisq_cumulants(d1,n1,sub)
            do r=1,size(kappa); kappa(r)=kappa(r)+sub(r); end do
            call lognoncentral_chisq_cumulants(d2,n2,sub)
            do r=1,size(kappa); kappa(r)=kappa(r)+(-1.0_dp)**r*sub(r); end do
            shift=shift+log(d1)-log(d2)
        end do
    end subroutine proddnf_log_cumulants

    pure real(dp) function dproddnf(x,df1,df2,ncp1,ncp2,log_density,order_max) result(v)
        real(dp),intent(in)::x,df1(:),df2(:),ncp1(:),ncp2(:); logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); real(dp)::shift; integer::om; logical::ld
        ld=get_log(log_density); if(x<=0.0_dp) then; if(ld) then; v=neg_inf_dp(); else; v=0.0_dp; end if; return; end if
        om=get_order(order_max,4); allocate(k(om)); call proddnf_log_cumulants(df1,df2,ncp1,ncp2,k,shift)
        v=edgeworth_pdf(log(x)+shift,k,log_density=ld); if(ld) then; v=v-log(x); else; v=v/x; end if
    end function dproddnf

    pure real(dp) function pproddnf(q,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,df1(:),df2(:),ncp1(:),ncp2(:); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); real(dp)::shift; integer::om; logical::lower,lp
        lower=get_lower(lower_tail);lp=get_log(log_p)
        if(q<=0.0_dp) then; if(lower) then; v=0.0_dp; else; v=1.0_dp; end if
            if(lp) then; if(v==0.0_dp) then; v=neg_inf_dp(); else; v=0.0_dp; end if; end if; return; end if
        om=get_order(order_max,4); allocate(k(om)); call proddnf_log_cumulants(df1,df2,ncp1,ncp2,k,shift)
        v=edgeworth_cdf(log(q)+shift,k,lower_tail=lower,log_p=lp)
    end function pproddnf

    pure real(dp) function qproddnf(p,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,df1(:),df2(:),ncp1(:),ncp2(:); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); real(dp)::shift; integer::om
        om=get_order(order_max,4); allocate(k(om)); call proddnf_log_cumulants(df1,df2,ncp1,ncp2,k,shift)
        v=exp(cornish_fisher_quantile(p,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))-shift)
    end function qproddnf

    subroutine rproddnf(x,df1,df2,ncp1,ncp2)
        real(dp),intent(out)::x(:); real(dp),intent(in)::df1(:),df2(:),ncp1(:),ncp2(:)
        real(dp)::d1,d2,n1,n2; integer::i,j,n
        n=max(size(df1),size(df2),size(ncp1),size(ncp2)); x=1.0_dp
        do j=1,n
            d1=recycled(df1,j);d2=recycled(df2,j);n1=recycled(ncp1,j);n2=recycled(ncp2,j)
            do i=1,size(x); x(i)=x(i)*(random_chisq(d1,n1)/d1)/(random_chisq(d2,n2)/d2); end do
        end do
    end subroutine rproddnf

    pure real(dp) function dprodnormal(x,mu,sigma,log_density,order_max) result(v)
        real(dp),intent(in)::x,mu(:),sigma(:); logical,intent(in),optional::log_density; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,5); allocate(k(om)); call prodnormal_cumulants(mu,sigma,k)
        v=edgeworth_pdf(x,k,log_density=get_log(log_density))
    end function dprodnormal

    pure real(dp) function pprodnormal(q,mu,sigma,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::q,mu(:),sigma(:); logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,5); allocate(k(om)); call prodnormal_cumulants(mu,sigma,k)
        v=edgeworth_cdf(q,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function pprodnormal

    pure real(dp) function qprodnormal(p,mu,sigma,lower_tail,log_p,order_max) result(v)
        real(dp),intent(in)::p,mu(:),sigma(:); logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max
        real(dp),allocatable::k(:); integer::om
        om=get_order(order_max,5); allocate(k(om)); call prodnormal_cumulants(mu,sigma,k)
        v=cornish_fisher_quantile(p,k,lower_tail=get_lower(lower_tail),log_p=get_log(log_p))
    end function qprodnormal

    subroutine rprodnormal(x,mu,sigma)
        real(dp),intent(out)::x(:); real(dp),intent(in)::mu(:),sigma(:); integer::i,j,n
        n=max(size(mu),size(sigma)); x=1.0_dp
        do j=1,n; do i=1,size(x); x(i)=x(i)*random_normal(recycled(mu,j),recycled(sigma,j)); end do; end do
    end subroutine rprodnormal

    ! Vector evaluation helpers for the four scalar-parameter doubly-noncentral families.
    pure subroutine ddnf_vec(x,out,df1,df2,ncp1,ncp2,log_density,order_max)
        real(dp),intent(in)::x(:),df1,df2,ncp1,ncp2; real(dp),intent(out)::out(size(x))
        logical,intent(in),optional::log_density; integer,intent(in),optional::order_max; integer::i
        do i=1,size(x); out(i)=ddnf(x(i),df1,df2,ncp1,ncp2,log_density,order_max); end do
    end subroutine ddnf_vec
    pure subroutine pdnf_vec(q,out,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::q(:),df1,df2,ncp1,ncp2; real(dp),intent(out)::out(size(q))
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max; integer::i
        do i=1,size(q); out(i)=pdnf(q(i),df1,df2,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine pdnf_vec
    pure subroutine qdnf_vec(p,out,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::p(:),df1,df2,ncp1,ncp2; real(dp),intent(out)::out(size(p))
        logical,intent(in),optional::lower_tail,log_p; integer,intent(in),optional::order_max; integer::i
        do i=1,size(p); out(i)=qdnf(p(i),df1,df2,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine qdnf_vec
    pure subroutine ddnt_vec(x,out,df,ncp1,ncp2,log_density,order_max)
        real(dp),intent(in)::x(:),df,ncp1,ncp2; real(dp),intent(out)::out(size(x)); logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(x); out(i)=ddnt(x(i),df,ncp1,ncp2,log_density,order_max); end do
    end subroutine ddnt_vec
    pure subroutine pdnt_vec(q,out,df,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::q(:),df,ncp1,ncp2; real(dp),intent(out)::out(size(q)); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(q); out(i)=pdnt(q(i),df,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine pdnt_vec
    pure subroutine qdnt_vec(p,out,df,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::p(:),df,ncp1,ncp2; real(dp),intent(out)::out(size(p)); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(p); out(i)=qdnt(p(i),df,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine qdnt_vec
    pure subroutine ddnbeta_vec(x,out,df1,df2,ncp1,ncp2,log_density,order_max)
        real(dp),intent(in)::x(:),df1,df2,ncp1,ncp2; real(dp),intent(out)::out(size(x)); logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(x); out(i)=ddnbeta(x(i),df1,df2,ncp1,ncp2,log_density,order_max); end do
    end subroutine ddnbeta_vec
    pure subroutine pdnbeta_vec(q,out,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::q(:),df1,df2,ncp1,ncp2
        real(dp),intent(out)::out(size(q)); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(q); out(i)=pdnbeta(q(i),df1,df2,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine pdnbeta_vec
    pure subroutine qdnbeta_vec(p,out,df1,df2,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::p(:),df1,df2,ncp1,ncp2
        real(dp),intent(out)::out(size(p)); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(p); out(i)=qdnbeta(p(i),df1,df2,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine qdnbeta_vec
    pure subroutine ddneta_vec(x,out,df,ncp1,ncp2,log_density,order_max)
        real(dp),intent(in)::x(:),df,ncp1,ncp2; real(dp),intent(out)::out(size(x)); logical,intent(in),optional::log_density
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(x); out(i)=ddneta(x(i),df,ncp1,ncp2,log_density,order_max); end do
    end subroutine ddneta_vec
    pure subroutine pdneta_vec(q,out,df,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::q(:),df,ncp1,ncp2; real(dp),intent(out)::out(size(q)); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(q); out(i)=pdneta(q(i),df,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine pdneta_vec
    pure subroutine qdneta_vec(p,out,df,ncp1,ncp2,lower_tail,log_p,order_max)
        real(dp),intent(in)::p(:),df,ncp1,ncp2; real(dp),intent(out)::out(size(p)); logical,intent(in),optional::lower_tail,log_p
        integer,intent(in),optional::order_max; integer::i
        do i=1,size(p); out(i)=qdneta(p(i),df,ncp1,ncp2,lower_tail,log_p,order_max); end do
    end subroutine qdneta_vec

end module sadists_distributions
