module actuar_risk
    use actuar_kinds, only: dp
    use actuar_discrete, only: dztpois, dztnbinom, dztbinom, dlogarithmic
    implicit none
    private
    public :: aggregate_dist_t, discretize_upper, discretize_lower
    public :: discretize_rounding, discretize_unbiased
    public :: panjer_ab, panjer_poisson, panjer_nbinom, panjer_binomial, panjer_logarithmic
    public :: convolve_pmf, aggregate_var, aggregate_cte
    public :: adjustment_coefficient_discrete, ruin_bound

    type :: aggregate_dist_t
        real(dp), allocatable :: pmf(:)
        real(dp), allocatable :: cdf(:)
        real(dp) :: x_scale = 1.0_dp
    contains
        procedure :: mean => aggregate_mean
        procedure :: variance => aggregate_variance
        procedure :: quantile => aggregate_quantile
    end type aggregate_dist_t

contains

    pure function discretize_upper(fvals) result(fx)
        real(dp),intent(in)::fvals(:)
        real(dp),allocatable::fx(:)
        integer::n
        n=size(fvals);allocate(fx(max(0,n-1)))
        if(n>1) fx=fvals(2:n)-fvals(1:n-1)
    end function discretize_upper

    pure function discretize_lower(fvals) result(fx)
        real(dp),intent(in)::fvals(:)
        real(dp),allocatable::fx(:)
        integer::n
        n=size(fvals);allocate(fx(n));fx=0.0_dp
        if(n>1) fx(2:n)=fvals(2:n)-fvals(1:n-1)
    end function discretize_lower

    pure function discretize_rounding(fvals) result(fx)
        real(dp),intent(in)::fvals(:)
        real(dp),allocatable::fx(:)
        integer::n
        n=size(fvals);allocate(fx(max(0,n-1)))
        if(n>1) fx=fvals(2:n)-fvals(1:n-1)
    end function discretize_rounding

    pure function discretize_unbiased(levvals, f_from, f_to, step) result(fx)
        real(dp),intent(in)::levvals(:),f_from,f_to,step
        real(dp),allocatable::fx(:)
        integer::n,i
        n=size(levvals);allocate(fx(n));fx=0.0_dp
        if(n<2) return
        fx(1)=-(levvals(2)-levvals(1))/step+1.0_dp-f_from
        do i=2,n-1
            fx(i)=(2.0_dp*levvals(i)-levvals(i-1)-levvals(i+1))/step
        end do
        fx(n)=(levvals(n)-levvals(n-1))/step-1.0_dp+f_to
        where(abs(fx)<1.0e-15_dp) fx=0.0_dp
    end function discretize_unbiased

    function panjer_ab(fx,a,b,fs0,tol,maxit,p0,p1) result(dist)
        real(dp),intent(in)::fx(:),a,b,fs0
        real(dp),intent(in),optional::tol,p0,p1
        integer,intent(in),optional::maxit
        type(aggregate_dist_t)::dist
        real(dp),allocatable::fs(:),tmp(:)
        real(dp)::target,cumul,norm,term,p0v,p1v
        integer::limit,x,k,m,nmax,newn
        logical::ab1
        nmax=500;if(present(maxit)) nmax=maxit
        target=1.0_dp-sqrt(epsilon(1.0_dp));if(present(tol)) target=1.0_dp-tol
        limit=size(fx)-1
        allocate(fs(0:nmax));fs=0.0_dp;fs(0)=fs0;cumul=fs0
        ab1=present(p0);p0v=0.0_dp;p1v=0.0_dp
        if(ab1) p0v=p0
        if(present(p1)) p1v=p1
        norm=1.0_dp-a*fx(1)
        if(abs(norm)<tiny(1.0_dp)) then
            allocate(dist%pmf(1),dist%cdf(1));dist%pmf=fs0;dist%cdf=fs0;return
        end if
        term=p1v-(a+b)*p0v
        x=1
        do while(cumul<target .and. x<=nmax)
            m=min(x,limit)
            do k=1,m
                fs(x)=fs(x)+(a+b*real(k,dp)/real(x,dp))*fx(k+1)*fs(x-k)
            end do
            if(ab1 .and. x<=limit) fs(x)=fs(x)+fx(x+1)*term
            fs(x)=fs(x)/norm
            if(fs(x)<0.0_dp .and. abs(fs(x))<1.0e-14_dp) fs(x)=0.0_dp
            cumul=cumul+fs(x);x=x+1
        end do
        newn=x
        allocate(dist%pmf(newn),dist%cdf(newn));dist%pmf=fs(0:newn-1)
        dist%cdf(1)=dist%pmf(1)
        do k=2,newn;dist%cdf(k)=min(1.0_dp,dist%cdf(k-1)+dist%pmf(k));end do
    end function panjer_ab

    function panjer_poisson(fx,lambda,tol,maxit) result(dist)
        real(dp),intent(in)::fx(:),lambda
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(aggregate_dist_t)::dist
        real(dp)::fs0,t
        integer::mi
        t=sqrt(epsilon(1.0_dp));if(present(tol)) t=tol
        mi=500;if(present(maxit)) mi=maxit
        fs0=exp(lambda*(fx(1)-1.0_dp))
        dist=panjer_ab(fx,0.0_dp,lambda,fs0,tol=t,maxit=mi)
    end function panjer_poisson

    function panjer_nbinom(fx,size,prob,tol,maxit) result(dist)
        real(dp),intent(in)::fx(:),size,prob
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(aggregate_dist_t)::dist
        real(dp)::a,b,fs0,t
        integer::mi
        t=sqrt(epsilon(1.0_dp));if(present(tol)) t=tol
        mi=500;if(present(maxit)) mi=maxit
        a=1.0_dp-prob;b=(size-1.0_dp)*a
        fs0=(prob/(1.0_dp-a*fx(1)))**size
        dist=panjer_ab(fx,a,b,fs0,tol=t,maxit=mi)
    end function panjer_nbinom

    function panjer_binomial(fx,n,prob,tol,maxit) result(dist)
        real(dp),intent(in)::fx(:),prob
        integer,intent(in)::n
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(aggregate_dist_t)::dist
        real(dp)::a,b,fs0,t
        integer::mi
        t=sqrt(epsilon(1.0_dp));if(present(tol)) t=tol
        mi=500;if(present(maxit)) mi=maxit
        a=-prob/(1.0_dp-prob);b=-real(n+1,dp)*a
        fs0=(1.0_dp+prob*(fx(1)-1.0_dp))**n
        dist=panjer_ab(fx,a,b,fs0,tol=t,maxit=mi)
    end function panjer_binomial

    function panjer_logarithmic(fx,prob,tol,maxit) result(dist)
        real(dp),intent(in)::fx(:),prob
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(aggregate_dist_t)::dist
        real(dp)::fs0,t
        integer::mi
        t=sqrt(epsilon(1.0_dp));if(present(tol)) t=tol
        mi=500;if(present(maxit)) mi=maxit
        fs0=log(1.0_dp-prob*fx(1))/log(1.0_dp-prob)
        dist=panjer_ab(fx,prob,-prob,fs0,tol=t,maxit=mi)
    end function panjer_logarithmic

    pure function convolve_pmf(a,b) result(c)
        real(dp),intent(in)::a(:),b(:)
        real(dp),allocatable::c(:)
        integer::i,j
        allocate(c(size(a)+size(b)-1));c=0.0_dp
        do i=1,size(a);do j=1,size(b);c(i+j-1)=c(i+j-1)+a(i)*b(j);end do;end do
    end function convolve_pmf

    pure real(dp) function aggregate_mean(self) result(m)
        class(aggregate_dist_t),intent(in)::self
        integer::i
        m=0.0_dp
        do i=1,size(self%pmf);m=m+real(i-1,dp)*self%x_scale*self%pmf(i);end do
    end function aggregate_mean

    pure real(dp) function aggregate_variance(self) result(v)
        class(aggregate_dist_t),intent(in)::self
        real(dp)::m,x
        integer::i
        m=self%mean();v=0.0_dp
        do i=1,size(self%pmf);x=real(i-1,dp)*self%x_scale;v=v+(x-m)**2*self%pmf(i);end do
    end function aggregate_variance

    pure real(dp) function aggregate_quantile(self,p) result(q)
        class(aggregate_dist_t),intent(in)::self
        real(dp),intent(in)::p
        integer::i
        q=0.0_dp
        do i=1,size(self%cdf)
            if(self%cdf(i)>=p) then;q=real(i-1,dp)*self%x_scale;return;end if
        end do
        q=real(size(self%cdf)-1,dp)*self%x_scale
    end function aggregate_quantile

    pure real(dp) function aggregate_var(dist,conf) result(v)
        type(aggregate_dist_t),intent(in)::dist
        real(dp),intent(in)::conf
        v=dist%quantile(conf)
    end function aggregate_var

    pure real(dp) function aggregate_cte(dist,conf) result(cte)
        type(aggregate_dist_t),intent(in)::dist
        real(dp),intent(in)::conf
        real(dp)::v,x,tailp,tailx
        integer::i
        v=dist%quantile(conf);tailp=0.0_dp;tailx=0.0_dp
        do i=1,size(dist%pmf)
            x=real(i-1,dp)*dist%x_scale
            if(x>v) then;tailp=tailp+dist%pmf(i);tailx=tailx+x*dist%pmf(i);end if
        end do
        if(tailp<=0.0_dp) then;cte=v;else;cte=(tailx+v*max(0.0_dp,1.0_dp-conf-tailp))/(1.0_dp-conf);end if
    end function aggregate_cte

    pure real(dp) function adjustment_coefficient_discrete(pmf,step,premium_rate) result(r)
        real(dp),intent(in)::pmf(:),step,premium_rate
        real(dp)::lo,hi,mid,g,meanx
        integer::i,it
        meanx=0.0_dp
        do i=1,size(pmf);meanx=meanx+real(i-1,dp)*step*pmf(i);end do
        if(premium_rate<=meanx) then;r=0.0_dp;return;end if
        lo=0.0_dp;hi=1.0_dp/max(step,1.0e-12_dp)
        do
            g=-premium_rate*hi
            do i=1,size(pmf);g=g+pmf(i)*(exp(hi*real(i-1,dp)*step)-1.0_dp);end do
            if(g>0.0_dp .or. hi>100.0_dp/step) exit
            hi=2.0_dp*hi
        end do
        do it=1,100
            mid=0.5_dp*(lo+hi);g=-premium_rate*mid
            do i=1,size(pmf);g=g+pmf(i)*(exp(mid*real(i-1,dp)*step)-1.0_dp);end do
            if(g>0.0_dp) then;hi=mid;else;lo=mid;end if
        end do
        r=0.5_dp*(lo+hi)
    end function adjustment_coefficient_discrete

    pure real(dp) function ruin_bound(initial_surplus,adjustment_coefficient) result(bound)
        real(dp),intent(in)::initial_surplus,adjustment_coefficient
        bound=exp(-adjustment_coefficient*initial_surplus)
    end function ruin_bound

end module actuar_risk
