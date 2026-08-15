module actuar_credibility
    use actuar_kinds, only: dp
    implicit none
    private
    public :: bstraub_result_t, bstraub_fit, bayes_linear_result_t
    public :: bayes_poisson_gamma, bayes_bernoulli_beta, bayes_normal_normal

    type :: bstraub_result_t
        real(dp) :: collective_mean = 0.0_dp
        real(dp) :: process_variance = 0.0_dp
        real(dp) :: between_variance = 0.0_dp
        real(dp), allocatable :: individual_mean(:)
        real(dp), allocatable :: credibility(:)
        real(dp), allocatable :: premium(:)
    end type bstraub_result_t

    type :: bayes_linear_result_t
        real(dp) :: collective_mean = 0.0_dp
        real(dp) :: k = 0.0_dp
        real(dp) :: credibility = 0.0_dp
        real(dp) :: individual_mean = 0.0_dp
        real(dp) :: premium = 0.0_dp
    end type bayes_linear_result_t

contains

    function bstraub_fit(ratios,weights,iterative,tol,maxit) result(res)
        real(dp),intent(in)::ratios(:,:),weights(:,:)
        logical,intent(in),optional::iterative
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(bstraub_result_t)::res
        integer::n,m,i,j,ntotal,ncontracts,it,mi
        real(dp),allocatable::ws(:),wm(:),cred(:)
        real(dp)::s2,a,old,wsum,xbar,den,t
        logical::iter
        n=size(ratios,1);m=size(ratios,2);allocate(ws(n),wm(n),cred(n))
        ws=0.0_dp;wm=0.0_dp;ntotal=0;ncontracts=0
        do i=1,n
            do j=1,m
                if(weights(i,j)>0.0_dp) then
                    ws(i)=ws(i)+weights(i,j);wm(i)=wm(i)+weights(i,j)*ratios(i,j);ntotal=ntotal+1
                end if
            end do
            if(ws(i)>0.0_dp) then;wm(i)=wm(i)/ws(i);ncontracts=ncontracts+1;end if
        end do
        s2=0.0_dp
        do i=1,n;do j=1,m;if(weights(i,j)>0.0_dp) s2=s2+weights(i,j)*(ratios(i,j)-wm(i))**2;end do;end do
        if(ntotal>ncontracts) s2=s2/real(ntotal-ncontracts,dp)
        wsum=sum(ws);xbar=dot_product(ws,wm)/wsum
        den=wsum**2-sum(ws**2)
        if(den>0.0_dp) then
            a=wsum*(dot_product(ws,(wm-xbar)**2)-real(ncontracts-1,dp)*s2)/den
        else
            a=0.0_dp
        end if
        a=max(0.0_dp,a)
        iter=.false.;if(present(iterative)) iter=iterative
        t=sqrt(epsilon(1.0_dp));if(present(tol)) t=tol
        mi=100;if(present(maxit)) mi=maxit
        if(iter .and. a>0.0_dp) then
            do it=1,mi
                old=a
                do i=1,n;cred(i)=1.0_dp/(1.0_dp+s2/(max(ws(i),tiny(1.0_dp))*a));end do
                xbar=dot_product(cred,wm)/sum(cred)
                a=dot_product(cred,(wm-xbar)**2)/real(max(1,ncontracts-1),dp)
                if(abs(a-old)<=t*max(old,tiny(1.0_dp))) exit
            end do
        end if
        if(a>0.0_dp) then
            do i=1,n;cred(i)=1.0_dp/(1.0_dp+s2/(max(ws(i),tiny(1.0_dp))*a));end do
            xbar=dot_product(cred,wm)/sum(cred)
        else
            cred=0.0_dp;xbar=dot_product(ws,wm)/wsum
        end if
        res%collective_mean=xbar;res%process_variance=s2;res%between_variance=a
        allocate(res%individual_mean(n),res%credibility(n),res%premium(n))
        res%individual_mean=wm;res%credibility=cred;res%premium=xbar+cred*(wm-xbar)
    end function bstraub_fit

    pure function bayes_poisson_gamma(x,shape,scale) result(res)
        real(dp),intent(in)::x(:),shape,scale
        type(bayes_linear_result_t)::res
        integer::n
        n=size(x);res%collective_mean=shape*scale;res%k=1.0_dp/scale
        if(n>0) res%individual_mean=sum(x)/real(n,dp)
        res%credibility=real(n,dp)/(real(n,dp)+res%k)
        res%premium=res%collective_mean+res%credibility*(res%individual_mean-res%collective_mean)
    end function bayes_poisson_gamma

    pure function bayes_bernoulli_beta(x,shape1,shape2) result(res)
        real(dp),intent(in)::x(:),shape1,shape2
        type(bayes_linear_result_t)::res
        integer::n
        n=size(x);res%k=shape1+shape2;res%collective_mean=shape1/res%k
        if(n>0) res%individual_mean=sum(x)/real(n,dp)
        res%credibility=real(n,dp)/(real(n,dp)+res%k)
        res%premium=res%collective_mean+res%credibility*(res%individual_mean-res%collective_mean)
    end function bayes_bernoulli_beta

    pure function bayes_normal_normal(x,prior_mean,prior_sd,lik_sd) result(res)
        real(dp),intent(in)::x(:),prior_mean,prior_sd,lik_sd
        type(bayes_linear_result_t)::res
        integer::n
        n=size(x);res%collective_mean=prior_mean;res%k=lik_sd**2/prior_sd**2
        if(n>0) res%individual_mean=sum(x)/real(n,dp)
        res%credibility=real(n,dp)/(real(n,dp)+res%k)
        res%premium=res%collective_mean+res%credibility*(res%individual_mean-res%collective_mean)
    end function bayes_normal_normal

end module actuar_credibility
