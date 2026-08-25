module greybox_distributions
    use greybox_kinds, only: dp, pi
    use greybox_special, only: normal_pdf, normal_cdf, normal_quantile, normal_rng, &
        gamma_p, gamma_quantile, nan_dp, inf_dp
    implicit none
    private
    public :: dlaplace, plaplace, qlaplace, rlaplace
    public :: dalaplace, palaplace, qalaplace, ralaplace
    public :: dgnorm, pgnorm, qgnorm, rgnorm
    public :: ds, ps, qs, rs
    public :: dfnorm, pfnorm, qfnorm, rfnorm
    public :: dbcnorm, pbcnorm, qbcnorm, rbcnorm
    public :: dlogitnorm, plogitnorm, qlogitnorm, rlogitnorm
    public :: drectnorm, prectnorm, qrectnorm, rrectnorm
    public :: dtplnorm, ptplnorm, qtplnorm, rtplnorm

contains

    pure elemental real(dp) function dlaplace(q,mu,scale,log_density) result(v)
        real(dp),intent(in)::q,mu,scale
        logical,intent(in),optional::log_density
        logical::lg
        lg=.false.; if(present(log_density))lg=log_density
        if(scale<=0.0_dp)then;v=nan_dp();return;end if
        if(lg)then
            v=-log(2.0_dp*scale)-abs(q-mu)/scale
        else
            v=exp(-abs(q-mu)/scale)/(2.0_dp*scale)
        end if
    end function dlaplace

    pure elemental real(dp) function plaplace(q,mu,scale) result(v)
        real(dp),intent(in)::q,mu,scale
        real(dp)::z
        if(scale<=0.0_dp)then;v=nan_dp();return;end if
        z=(q-mu)/scale
        if(z<0.0_dp)then;v=0.5_dp*exp(z);else;v=1.0_dp-0.5_dp*exp(-z);end if
    end function plaplace

    pure elemental real(dp) function qlaplace(p,mu,scale) result(v)
        real(dp),intent(in)::p,mu,scale
        if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;v=nan_dp();return;end if
        if(p<=0.0_dp)then;v=-inf_dp();else if(p>=1.0_dp)then;v=inf_dp(); &
        else if(p<0.5_dp)then;v=mu+scale*log(2.0_dp*p); &
        else;v=mu-scale*log(2.0_dp*(1.0_dp-p));end if
    end function qlaplace

    subroutine rlaplace(x,mu,scale)
        real(dp),intent(out)::x(:)
        real(dp),intent(in)::mu,scale
        real(dp)::u
        integer::i
        do i=1,size(x);call random_number(u);x(i)=qlaplace(u,mu,scale);end do
    end subroutine rlaplace

    pure elemental real(dp) function dalaplace(q,mu,scale,alpha,log_density) result(v)
        real(dp),intent(in)::q,mu,scale,alpha
        logical,intent(in),optional::log_density
        logical::lg
        real(dp)::l
        lg=.false.;if(present(log_density))lg=log_density
        if(scale<=0.0_dp.or.alpha<=0.0_dp.or.alpha>=1.0_dp)then;v=nan_dp();return;end if
        l=log(alpha)+log(1.0_dp-alpha)-log(scale)-(q-mu)/scale*(alpha-merge(1.0_dp,0.0_dp,q<=mu))
        if(lg)then;v=l;else;v=exp(l);end if
    end function dalaplace

    pure elemental real(dp) function palaplace(q,mu,scale,alpha) result(v)
        real(dp),intent(in)::q,mu,scale,alpha
        if(scale<=0.0_dp.or.alpha<=0.0_dp.or.alpha>=1.0_dp)then;v=nan_dp();return;end if
        if(q<=mu)then
            v=alpha*exp((1.0_dp-alpha)*(q-mu)/scale)
        else
            v=1.0_dp-(1.0_dp-alpha)*exp(-alpha*(q-mu)/scale)
        end if
    end function palaplace

    pure elemental real(dp) function qalaplace(p,mu,scale,alpha) result(v)
        real(dp),intent(in)::p,mu,scale,alpha
        if(scale<=0.0_dp.or.alpha<=0.0_dp.or.alpha>=1.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
            v=nan_dp();return
        end if
        if(p<=0.0_dp)then;v=-inf_dp();else if(p>=1.0_dp)then;v=inf_dp(); &
        else if(p<=alpha)then;v=mu+scale/(1.0_dp-alpha)*log(p/alpha); &
        else;v=mu-scale/alpha*log((1.0_dp-p)/(1.0_dp-alpha));end if
    end function qalaplace

    subroutine ralaplace(x,mu,scale,alpha)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,scale,alpha
        real(dp)::u;integer::i
        do i=1,size(x);call random_number(u);x(i)=qalaplace(u,mu,scale,alpha);end do
    end subroutine ralaplace

    pure elemental real(dp) function dgnorm(q,mu,scale,shape,log_density) result(v)
        real(dp),intent(in)::q,mu,scale,shape
        logical,intent(in),optional::log_density
        logical::lg;real(dp)::l,z
        lg=.false.;if(present(log_density))lg=log_density
        if(scale<=0.0_dp.or.shape<=0.0_dp)then;v=nan_dp();return;end if
        z=abs(q-mu)/scale
        l=log(shape)-log(2.0_dp*scale)-log_gamma(1.0_dp/shape)-z**shape
        if(lg)then;v=l;else;v=exp(l);end if
    end function dgnorm

    pure elemental real(dp) function pgnorm(q,mu,scale,shape,lower_tail,log_p) result(v)
        real(dp),intent(in)::q,mu,scale,shape
        logical,intent(in),optional::lower_tail,log_p
        logical::lower,lp;real(dp)::p,z
        lower=.true.;lp=.false.;if(present(lower_tail))lower=lower_tail;if(present(log_p))lp=log_p
        if(scale<=0.0_dp.or.shape<=0.0_dp)then;v=nan_dp();return;end if
        if(shape>100.0_dp)then
            p=min(1.0_dp,max(0.0_dp,(q-(mu-scale))/(2.0_dp*scale)))
        else
            z=(abs(q-mu)/scale)**shape
            p=0.5_dp+0.5_dp*sign(1.0_dp,q-mu)*gamma_p(1.0_dp/shape,z)
            if(abs(q-mu)<=epsilon(1.0_dp)*max(1.0_dp,abs(mu)))p=0.5_dp
        end if
        if(.not.lower)p=1.0_dp-p
        if(lp)then;v=log(p);else;v=p;end if
    end function pgnorm

    pure elemental real(dp) function qgnorm(p,mu,scale,shape,lower_tail,log_p) result(v)
        real(dp),intent(in)::p,mu,scale,shape
        logical,intent(in),optional::lower_tail,log_p
        logical::lower,lp;real(dp)::pp,g
        lower=.true.;lp=.false.;if(present(lower_tail))lower=lower_tail;if(present(log_p))lp=log_p
        if(scale<=0.0_dp.or.shape<=0.0_dp)then;v=nan_dp();return;end if
        pp=p;if(lp)pp=exp(pp);if(.not.lower)pp=1.0_dp-pp
        if(pp<0.0_dp.or.pp>1.0_dp)then;v=nan_dp();return;end if
        if(pp<=0.0_dp)then;v=-inf_dp();return;else if(pp>=1.0_dp)then;v=inf_dp();return;end if
        if(shape>100.0_dp)then
            v=mu-scale+2.0_dp*scale*pp;return
        end if
        if(abs(pp-0.5_dp)<=epsilon(1.0_dp))then;v=mu;return;end if
        g=gamma_quantile(2.0_dp*abs(pp-0.5_dp),1.0_dp/shape,1.0_dp)
        v=mu+sign(1.0_dp,pp-0.5_dp)*scale*g**(1.0_dp/shape)
    end function qgnorm

    subroutine rgnorm(x,mu,scale,shape)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,scale,shape
        real(dp)::u;integer::i
        do i=1,size(x);call random_number(u);x(i)=qgnorm(u,mu,scale,shape);end do
    end subroutine rgnorm

    pure elemental real(dp) function ds(q,mu,scale,log_density) result(v)
        real(dp),intent(in)::q,mu,scale;logical,intent(in),optional::log_density
        v=dgnorm(q,mu,scale*scale,0.5_dp,log_density)
    end function ds
    pure elemental real(dp) function ps(q,mu,scale) result(v)
        real(dp),intent(in)::q,mu,scale
        v=pgnorm(q,mu,scale*scale,0.5_dp)
    end function ps
    pure elemental real(dp) function qs(p,mu,scale) result(v)
        real(dp),intent(in)::p,mu,scale
        v=qgnorm(p,mu,scale*scale,0.5_dp)
    end function qs
    subroutine rs(x,mu,scale)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,scale
        call rgnorm(x,mu,scale*scale,0.5_dp)
    end subroutine rs

    pure elemental real(dp) function dfnorm(q,mu,sigma,log_density) result(v)
        real(dp),intent(in)::q,mu,sigma;logical,intent(in),optional::log_density
        logical::lg;real(dp)::f
        lg=.false.;if(present(log_density))lg=log_density
        if(q<0.0_dp)then;f=0.0_dp;else;f=normal_pdf(q,mu,sigma)+normal_pdf(q,-mu,sigma);end if
        if(lg)then;if(f>0.0_dp)then;v=log(f);else;v=-inf_dp();end if;else;v=f;end if
    end function dfnorm
    pure elemental real(dp) function pfnorm(q,mu,sigma) result(v)
        real(dp),intent(in)::q,mu,sigma
        if(q<0.0_dp)then;v=0.0_dp;else;v=normal_cdf(q,mu,sigma)+normal_cdf(q,-mu,sigma)-1.0_dp;end if
    end function pfnorm
    pure elemental real(dp) function qfnorm(p,mu,sigma) result(v)
        real(dp),intent(in)::p,mu,sigma
        real(dp)::lo,hi,mid;integer::i
        if(sigma<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;v=nan_dp();return;end if
        if(p<=0.0_dp)then;v=0.0_dp;return;else if(p>=1.0_dp)then;v=inf_dp();return;end if
        lo=0.0_dp;hi=max(abs(mu)+8.0_dp*sigma,sigma)
        do while(pfnorm(hi,mu,sigma)<p);hi=2.0_dp*hi;end do
        do i=1,100;mid=0.5_dp*(lo+hi);if(pfnorm(mid,mu,sigma)<p)then;lo=mid;else;hi=mid;end if;end do
        v=0.5_dp*(lo+hi)
    end function qfnorm
    subroutine rfnorm(x,mu,sigma)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,sigma;integer::i
        do i=1,size(x);x(i)=abs(normal_rng(mu,sigma));end do
    end subroutine rfnorm

    pure elemental real(dp) function dbcnorm(q,mu,sigma,lambda,log_density) result(v)
        real(dp),intent(in)::q,mu,sigma,lambda;logical,intent(in),optional::log_density
        logical::lg;real(dp)::f,t
        lg=.false.;if(present(log_density))lg=log_density
        if(sigma<=0.0_dp)then;v=nan_dp();return;end if
        if(abs(lambda)<=epsilon(1.0_dp))then
            if(q<=0.0_dp)then;f=0.0_dp;else;f=normal_pdf(log(q),mu,sigma)/q;end if
        else if(abs(lambda-1.0_dp)<=epsilon(1.0_dp))then
            f=normal_pdf(q,mu+1.0_dp,sigma);if(q<=0.0_dp)f=0.0_dp
        else if(q<=0.0_dp)then
            f=0.0_dp
        else
            t=(q**lambda-1.0_dp)/lambda
            f=q**(lambda-1.0_dp)*normal_pdf(t,mu,sigma)
        end if
        if(lg)then;if(f>0.0_dp)then;v=log(f);else;v=-inf_dp();end if;else;v=f;end if
    end function dbcnorm
    pure elemental real(dp) function pbcnorm(q,mu,sigma,lambda) result(v)
        real(dp),intent(in)::q,mu,sigma,lambda
        if(q<=0.0_dp)then;v=0.0_dp;else if(abs(lambda)<=epsilon(1.0_dp))then;v=normal_cdf(log(q),mu,sigma); &
        else;v=normal_cdf((q**lambda-1.0_dp)/lambda,mu,sigma);end if
    end function pbcnorm
    pure elemental real(dp) function qbcnorm(p,mu,sigma,lambda) result(v)
        real(dp),intent(in)::p,mu,sigma,lambda;real(dp)::z,b
        if(sigma<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;v=nan_dp();return;end if
        z=normal_quantile(p,mu,sigma)
        if (abs(lambda) <= epsilon(1.0_dp)) then
            v = exp(z)
        else
            b = lambda*z + 1.0_dp
            if (b <= 0.0_dp) then
                v = 0.0_dp
            else
                v = b**(1.0_dp/lambda)
            end if
        end if
    end function qbcnorm
    subroutine rbcnorm(x,mu,sigma,lambda)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,sigma,lambda;real(dp)::u;integer::i
        do i=1,size(x);call random_number(u);x(i)=qbcnorm(u,mu,sigma,lambda);end do
    end subroutine rbcnorm

    pure elemental real(dp) function dlogitnorm(q,mu,sigma,log_density) result(v)
        real(dp),intent(in)::q,mu,sigma;logical,intent(in),optional::log_density
        logical::lg;real(dp)::l,z
        lg=.false.;if(present(log_density))lg=log_density
        if(q<=0.0_dp.or.q>=1.0_dp.or.sigma<=0.0_dp)then
            if(lg)then;v=-inf_dp();else;v=0.0_dp;end if;return
        end if
        z=log(q/(1.0_dp-q))
        l=log(normal_pdf(z,mu,sigma))-log(q)-log(1.0_dp-q)
        if(lg)then;v=l;else;v=exp(l);end if
    end function dlogitnorm
    pure elemental real(dp) function plogitnorm(q,mu,sigma) result(v)
        real(dp),intent(in)::q,mu,sigma
        if(q<=0.0_dp)then;v=0.0_dp;else if(q>=1.0_dp)then;v=1.0_dp;else;v=normal_cdf(log(q/(1.0_dp-q)),mu,sigma);end if
    end function plogitnorm
    pure elemental real(dp) function qlogitnorm(p,mu,sigma) result(v)
        real(dp),intent(in)::p,mu,sigma;real(dp)::z
        z=normal_quantile(p,mu,sigma)
        if(z>=0.0_dp)then;v=1.0_dp/(1.0_dp+exp(-z));else;v=exp(z)/(1.0_dp+exp(z));end if
    end function qlogitnorm
    subroutine rlogitnorm(x,mu,sigma)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,sigma;real(dp)::u;integer::i
        do i=1,size(x);call random_number(u);x(i)=qlogitnorm(u,mu,sigma);end do
    end subroutine rlogitnorm

    pure elemental real(dp) function drectnorm(q,mu,sigma,log_density) result(v)
        real(dp),intent(in)::q,mu,sigma;logical,intent(in),optional::log_density
        logical::lg;real(dp)::f
        lg=.false.;if(present(log_density))lg=log_density
        if(q>0.0_dp)then;f=normal_pdf(q,mu,sigma);else;f=normal_cdf(0.0_dp,mu,sigma);end if
        if(lg)then;v=log(f);else;v=f;end if
    end function drectnorm
    pure elemental real(dp) function prectnorm(q,mu,sigma) result(v)
        real(dp),intent(in)::q,mu,sigma
        if(q>0.0_dp)then;v=normal_cdf(q,mu,sigma);else;v=normal_cdf(0.0_dp,mu,sigma);end if
    end function prectnorm
    pure elemental real(dp) function qrectnorm(p,mu,sigma) result(v)
        real(dp),intent(in)::p,mu,sigma
        v=max(normal_quantile(p,mu,sigma),0.0_dp)
    end function qrectnorm
    subroutine rrectnorm(x,mu,sigma)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,sigma;integer::i
        do i=1,size(x);x(i)=max(normal_rng(mu,sigma),0.0_dp);end do
    end subroutine rrectnorm

    pure elemental real(dp) function dtplnorm(q,mu,sigma,shift,log_density) result(v)
        real(dp),intent(in)::q,mu,sigma,shift;logical,intent(in),optional::log_density
        logical::lg;real(dp)::z,f
        lg=.false.;if(present(log_density))lg=log_density
        if(q<=shift.or.sigma<=0.0_dp)then
            if(lg)then;v=-inf_dp();else;v=0.0_dp;end if;return
        end if
        z=q-shift;f=normal_pdf(log(z),mu,sigma)/z
        if(lg)then;v=log(f);else;v=f;end if
    end function dtplnorm
    pure elemental real(dp) function ptplnorm(q,mu,sigma,shift) result(v)
        real(dp),intent(in)::q,mu,sigma,shift
        if(q<=shift)then;v=0.0_dp;else;v=normal_cdf(log(q-shift),mu,sigma);end if
    end function ptplnorm
    pure elemental real(dp) function qtplnorm(p,mu,sigma,shift) result(v)
        real(dp),intent(in)::p,mu,sigma,shift
        v=exp(normal_quantile(p,mu,sigma))+shift
    end function qtplnorm
    subroutine rtplnorm(x,mu,sigma,shift)
        real(dp),intent(out)::x(:);real(dp),intent(in)::mu,sigma,shift;integer::i
        do i=1,size(x);x(i)=exp(normal_rng(mu,sigma))+shift;end do
    end subroutine rtplnorm

end module greybox_distributions
