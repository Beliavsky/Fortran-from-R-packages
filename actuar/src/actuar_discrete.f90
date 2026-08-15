module actuar_discrete
    use actuar_kinds, only: dp
    use actuar_special, only: poisson_pmf, poisson_cdf, poisson_quantile, &
                              nbinom_pmf, nbinom_cdf, nbinom_quantile, &
                              binom_pmf, binom_cdf, binom_quantile, &
                              random_poisson, random_negative_binomial, random_binomial
    implicit none
    private
    public :: dlogarithmic, plogarithmic, qlogarithmic, rlogarithmic
    public :: dztpois, pztpois, qztpois, rztpois
    public :: dztgeom, pztgeom, qztgeom, rztgeom
    public :: dztnbinom, pztnbinom, qztnbinom, rztnbinom
    public :: dztbinom, pztbinom, qztbinom, rztbinom
    public :: dzmpois, pzmpois, qzmpois, rzmpois
    public :: dzmgeom, pzmgeom, qzmgeom, rzmgeom
    public :: dzmnbinom, pzmnbinom, qzmnbinom, rzmnbinom
    public :: dzmbinom, pzmbinom, qzmbinom, rzmbinom
    public :: dzmlogarithmic, pzmlogarithmic, qzmlogarithmic, rzmlogarithmic
    public :: dpoisinvgauss, ppoisinvgauss, qpoisinvgauss, rpoisinvgauss

contains

    pure real(dp) function dlogarithmic(k,prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob
        real(dp)::a
        if(k<1 .or. prob<0.0_dp .or. prob>=1.0_dp) then;p=0.0_dp;return;end if
        if(prob==0.0_dp) then;p=merge(1.0_dp,0.0_dp,k==1);return;end if
        a=-1.0_dp/log(1.0_dp-prob)
        p=a*prob**k/real(k,dp)
    end function dlogarithmic

    pure real(dp) function plogarithmic(k,prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob
        integer::j
        if(k<1) then;p=0.0_dp;return;end if
        p=0.0_dp
        do j=1,k;p=p+dlogarithmic(j,prob);end do
        p=min(1.0_dp,p)
    end function plogarithmic

    pure integer function qlogarithmic(q,prob) result(k)
        real(dp),intent(in)::q,prob
        if(q<=0.0_dp) then;k=1;return;end if
        k=1
        do while(plogarithmic(k,prob)<q .and. k<100000);k=k+1;end do
    end function qlogarithmic

    integer function rlogarithmic(prob) result(k)
        real(dp),intent(in)::prob
        real(dp)::u
        call random_number(u);k=qlogarithmic(u,prob)
    end function rlogarithmic

    pure real(dp) function dztpois(k,lambda) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::lambda
        real(dp)::z
        if(k<1) then;p=0.0_dp;return;end if
        if(lambda==0.0_dp) then;p=merge(1.0_dp,0.0_dp,k==1);return;end if
        z=1.0_dp-exp(-lambda);p=poisson_pmf(k,lambda)/z
    end function dztpois
    pure real(dp) function pztpois(k,lambda) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::lambda
        real(dp)::p0
        if(k<1) then;p=0.0_dp;return;end if
        if(lambda==0.0_dp) then;p=1.0_dp;return;end if
        p0=exp(-lambda);p=(poisson_cdf(k,lambda)-p0)/(1.0_dp-p0);p=max(0.0_dp,min(1.0_dp,p))
    end function pztpois
    pure integer function qztpois(q,lambda) result(k)
        real(dp),intent(in)::q,lambda
        real(dp)::p0
        if(lambda==0.0_dp) then;k=1;return;end if
        p0=exp(-lambda);k=max(1,poisson_quantile(p0+(1.0_dp-p0)*q,lambda))
    end function qztpois
    integer function rztpois(lambda) result(k)
        real(dp),intent(in)::lambda;real(dp)::u;call random_number(u);k=qztpois(u,lambda)
    end function rztpois

    pure real(dp) function dztgeom(k,prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob
        if(k<1 .or. prob<=0.0_dp .or. prob>1.0_dp) then;p=0.0_dp;else;p=prob*(1.0_dp-prob)**(k-1);end if
    end function dztgeom
    pure real(dp) function pztgeom(k,prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob
        if(k<1) then;p=0.0_dp;else;p=1.0_dp-(1.0_dp-prob)**k;end if
    end function pztgeom
    pure integer function qztgeom(q,prob) result(k)
        real(dp),intent(in)::q,prob
        if(q<=0.0_dp) then;k=1;else if(q>=1.0_dp) then;k=huge(k);else;k=max(1,ceiling(log(1.0_dp-q)/log(1.0_dp-prob)));end if
    end function qztgeom
    integer function rztgeom(prob) result(k)
        real(dp),intent(in)::prob;real(dp)::u;call random_number(u);k=qztgeom(u,prob)
    end function rztgeom

    pure real(dp) function dztnbinom(k,size,prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::size,prob
        real(dp)::p0
        if(k<1) then;p=0.0_dp;return;end if
        if(size==0.0_dp) then;p=dlogarithmic(k,1.0_dp-prob);return;end if
        p0=prob**size;p=nbinom_pmf(k,size,prob)/(1.0_dp-p0)
    end function dztnbinom
    pure real(dp) function pztnbinom(k,size,prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::size,prob
        real(dp)::p0
        if(k<1) then;p=0.0_dp;return;end if
        if(size==0.0_dp) then;p=plogarithmic(k,1.0_dp-prob);return;end if
        p0=prob**size;p=(nbinom_cdf(k,size,prob)-p0)/(1.0_dp-p0);p=max(0.0_dp,min(1.0_dp,p))
    end function pztnbinom
    pure integer function qztnbinom(q,size,prob) result(k)
        real(dp),intent(in)::q,size,prob
        real(dp)::p0
        if(size==0.0_dp) then;k=qlogarithmic(q,1.0_dp-prob);return;end if
        p0=prob**size;k=max(1,nbinom_quantile(p0+(1.0_dp-p0)*q,size,prob))
    end function qztnbinom
    integer function rztnbinom(size,prob) result(k)
        real(dp),intent(in)::size,prob;real(dp)::u;call random_number(u);k=qztnbinom(u,size,prob)
    end function rztnbinom

    pure real(dp) function dztbinom(k,n,prob) result(p)
        integer,intent(in)::k,n
        real(dp),intent(in)::prob
        real(dp)::p0
        if(k<1) then;p=0.0_dp;return;end if
        p0=(1.0_dp-prob)**n
        if(p0>=1.0_dp) then;p=merge(1.0_dp,0.0_dp,k==1);else;p=binom_pmf(k,n,prob)/(1.0_dp-p0);end if
    end function dztbinom
    pure real(dp) function pztbinom(k,n,prob) result(p)
        integer,intent(in)::k,n
        real(dp),intent(in)::prob
        real(dp)::p0
        if(k<1) then;p=0.0_dp;return;end if
        p0=(1.0_dp-prob)**n
        if(p0>=1.0_dp) then;p=1.0_dp;else;p=(binom_cdf(k,n,prob)-p0)/(1.0_dp-p0);end if
    end function pztbinom
    pure integer function qztbinom(q,n,prob) result(k)
        real(dp),intent(in)::q,prob
        integer,intent(in)::n
        real(dp)::p0
        p0=(1.0_dp-prob)**n
        if(p0>=1.0_dp) then;k=1;else;k=max(1,binom_quantile(p0+(1.0_dp-p0)*q,n,prob));end if
    end function qztbinom
    integer function rztbinom(n,prob) result(k)
        integer,intent(in)::n;real(dp),intent(in)::prob;real(dp)::u;call random_number(u);k=qztbinom(u,n,prob)
    end function rztbinom

    pure real(dp) function dzmpois(k,lambda,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::lambda,p0m
        real(dp)::p0
        if(k<0) then;p=0.0_dp;return;end if
        if(k==0) then;p=p0m;return;end if
        p0=exp(-lambda)
        if(p0>=1.0_dp) then;p=merge(1.0_dp-p0m,0.0_dp,k==1);else;p=(1.0_dp-p0m)*poisson_pmf(k,lambda)/(1.0_dp-p0);end if
    end function dzmpois
    pure real(dp) function pzmpois(k,lambda,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::lambda,p0m
        real(dp)::p0
        if(k<0) then;p=0.0_dp;return;end if
        if(k==0) then;p=p0m;return;end if
        p0=exp(-lambda)
        if(p0>=1.0_dp) then;p=1.0_dp;else;p=p0m+(1.0_dp-p0m)*(poisson_cdf(k,lambda)-p0)/(1.0_dp-p0);end if
        p=min(1.0_dp,p)
    end function pzmpois
    pure integer function qzmpois(q,lambda,p0m) result(k)
        real(dp),intent(in)::q,lambda,p0m
        if(q<=p0m) then;k=0;else;k=qztpois((q-p0m)/(1.0_dp-p0m),lambda);end if
    end function qzmpois
    integer function rzmpois(lambda,p0m) result(k)
        real(dp),intent(in)::lambda,p0m;real(dp)::u;call random_number(u);k=qzmpois(u,lambda,p0m)
    end function rzmpois

    pure real(dp) function dzmgeom(k,prob,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=(1.0_dp-p0m)*dztgeom(k,prob);end if
    end function dzmgeom
    pure real(dp) function pzmgeom(k,prob,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=p0m+(1.0_dp-p0m)*pztgeom(k,prob);end if
    end function pzmgeom
    pure integer function qzmgeom(q,prob,p0m) result(k)
        real(dp),intent(in)::q,prob,p0m
        if(q<=p0m) then;k=0;else;k=qztgeom((q-p0m)/(1.0_dp-p0m),prob);end if
    end function qzmgeom
    integer function rzmgeom(prob,p0m) result(k)
        real(dp),intent(in)::prob,p0m;real(dp)::u;call random_number(u);k=qzmgeom(u,prob,p0m)
    end function rzmgeom

    pure real(dp) function dzmnbinom(k,size,prob,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::size,prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=(1.0_dp-p0m)*dztnbinom(k,size,prob);end if
    end function dzmnbinom
    pure real(dp) function pzmnbinom(k,size,prob,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::size,prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=p0m+(1.0_dp-p0m)*pztnbinom(k,size,prob);end if
    end function pzmnbinom
    pure integer function qzmnbinom(q,size,prob,p0m) result(k)
        real(dp),intent(in)::q,size,prob,p0m
        if(q<=p0m) then;k=0;else;k=qztnbinom((q-p0m)/(1.0_dp-p0m),size,prob);end if
    end function qzmnbinom
    integer function rzmnbinom(size,prob,p0m) result(k)
        real(dp),intent(in)::size,prob,p0m;real(dp)::u;call random_number(u);k=qzmnbinom(u,size,prob,p0m)
    end function rzmnbinom

    pure real(dp) function dzmbinom(k,n,prob,p0m) result(p)
        integer,intent(in)::k,n
        real(dp),intent(in)::prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=(1.0_dp-p0m)*dztbinom(k,n,prob);end if
    end function dzmbinom
    pure real(dp) function pzmbinom(k,n,prob,p0m) result(p)
        integer,intent(in)::k,n
        real(dp),intent(in)::prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=p0m+(1.0_dp-p0m)*pztbinom(k,n,prob);end if
    end function pzmbinom
    pure integer function qzmbinom(q,n,prob,p0m) result(k)
        real(dp),intent(in)::q,prob,p0m
        integer,intent(in)::n
        if(q<=p0m) then;k=0;else;k=qztbinom((q-p0m)/(1.0_dp-p0m),n,prob);end if
    end function qzmbinom
    integer function rzmbinom(n,prob,p0m) result(k)
        integer,intent(in)::n;real(dp),intent(in)::prob,p0m;real(dp)::u;call random_number(u);k=qzmbinom(u,n,prob,p0m)
    end function rzmbinom

    pure real(dp) function dzmlogarithmic(k,prob,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=(1.0_dp-p0m)*dlogarithmic(k,prob);end if
    end function dzmlogarithmic
    pure real(dp) function pzmlogarithmic(k,prob,p0m) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::prob,p0m
        if(k<0) then;p=0.0_dp;else if(k==0) then;p=p0m;else;p=p0m+(1.0_dp-p0m)*plogarithmic(k,prob);end if
    end function pzmlogarithmic
    pure integer function qzmlogarithmic(q,prob,p0m) result(k)
        real(dp),intent(in)::q,prob,p0m
        if(q<=p0m) then;k=0;else;k=qlogarithmic((q-p0m)/(1.0_dp-p0m),prob);end if
    end function qzmlogarithmic
    integer function rzmlogarithmic(prob,p0m) result(k)
        real(dp),intent(in)::prob,p0m;real(dp)::u;call random_number(u);k=qzmlogarithmic(u,prob,p0m)
    end function rzmlogarithmic

    pure real(dp) function dpoisinvgauss(k,mu,phi) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::mu,phi
        integer::i
        real(dp)::p0,p1,a,b,pm2,pm1,cur,two
        if(k<0 .or. mu<=0.0_dp .or. phi<=0.0_dp) then;p=0.0_dp;return;end if
        two=2.0_dp*phi
        p0=exp((1.0_dp-sqrt(1.0_dp+two*mu*mu))/(phi*mu))
        if(k==0) then;p=p0;return;end if
        p1=mu*p0/sqrt(1.0_dp+two*mu*mu)
        if(k==1) then;p=p1;return;end if
        a=1.0_dp/(1.0_dp+1.0_dp/(two*mu*mu))
        b=mu*mu/(1.0_dp+two*mu*mu)
        pm2=p0;pm1=p1
        do i=2,k
            cur=a*(1.0_dp-1.5_dp/real(i,dp))*pm1+b*pm2/(real(i,dp)*real(i-1,dp))
            pm2=pm1;pm1=cur
        end do
        p=pm1
    end function dpoisinvgauss
    pure real(dp) function ppoisinvgauss(k,mu,phi) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::mu,phi
        integer::j
        if(k<0) then;p=0.0_dp;return;end if
        p=0.0_dp;do j=0,k;p=p+dpoisinvgauss(j,mu,phi);end do;p=min(1.0_dp,p)
    end function ppoisinvgauss
    pure integer function qpoisinvgauss(q,mu,phi) result(k)
        real(dp),intent(in)::q,mu,phi
        k=0;do while(ppoisinvgauss(k,mu,phi)<q .and. k<100000);k=k+1;end do
    end function qpoisinvgauss
    integer function rpoisinvgauss(mu,phi) result(k)
        real(dp),intent(in)::mu,phi
        real(dp)::v,y,x,u
        v=mu*mu*phi
        y=random_normal_local()**2
        x=mu+0.5_dp*mu*v*y-0.5_dp*mu*sqrt(4.0_dp*v*y+v*v*y*y)
        call random_number(u)
        if(u>mu/(mu+x)) x=mu*mu/x
        k=random_poisson(x)
    end function rpoisinvgauss
    real(dp) function random_normal_local() result(z)
        real(dp)::u1,u2
        call random_number(u1);call random_number(u2)
        z=sqrt(-2.0_dp*log(max(u1,tiny(1.0_dp))))*cos(6.2831853071795864769_dp*u2)
    end function random_normal_local

end module actuar_discrete
