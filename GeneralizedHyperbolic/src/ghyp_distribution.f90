module ghyp_distribution
  use gh_math, only: dp, pi, bessel_k, normal_rand
  use gig_distribution, only: rgig, gig_moment
  implicit none
  private
  public :: dghyp,pghyp,qghyp,rghyp,ghyp_mean,ghyp_var,ghyp_mode,ghyp_change_pars,ghyp_scale
  public :: dhyperb, phyperb, qhyperb, rhyperb, dnig, pnig, qnig, rnig
contains
  function dghyp(x,mu,delta,alpha,beta,lambda) result(d)
    real(dp),intent(in)::x,mu,delta,alpha,beta,lambda
    real(dp)::d,gamma,y,knum,kden,logd
    if(delta<=0 .or. alpha<=abs(beta))then;d=0;return;endif
    gamma=sqrt(alpha*alpha-beta*beta); y=alpha*sqrt(delta*delta+(x-mu)**2)
    knum=bessel_k(lambda-0.5_dp,y); kden=bessel_k(lambda,delta*gamma)
    if(knum<=0.or.kden<=0)then;d=0;return;endif
    logd=(lambda-0.5_dp)*log(sqrt(delta*delta+(x-mu)**2)) + lambda*log(gamma/delta) &
      + beta*(x-mu)+log(knum)-log(kden)-0.5_dp*log(2*pi)- (lambda-0.5_dp)*log(alpha)
    d=exp(logd)
  end function

  function pghyp(q,mu,delta,alpha,beta,lambda) result(p)
    real(dp),intent(in)::q,mu,delta,alpha,beta,lambda
    real(dp)::p,lo,h,x,s,sd
    integer::i,n
    sd=sqrt(max(ghyp_var(mu,delta,alpha,beta,lambda),1.0e-12_dp))
    lo=min(mu-20*sd,q-20*sd); n=400; h=(q-lo)/n
    if(h<=0)then;p=0;return;endif
    s=0
    do i=0,n
      x=lo+h*i
      if(i==0.or.i==n)then;s=s+dghyp(x,mu,delta,alpha,beta,lambda)
      elseif(mod(i,2)==0)then;s=s+2*dghyp(x,mu,delta,alpha,beta,lambda)
      else;s=s+4*dghyp(x,mu,delta,alpha,beta,lambda);endif
    enddo
    p=max(0.0_dp,min(1.0_dp,h*s/3))
  end function

  function qghyp(p,mu,delta,alpha,beta,lambda) result(q)
    real(dp),intent(in)::p,mu,delta,alpha,beta,lambda
    real(dp)::q,lo,hi,mid,sd
    integer::i
    if(p<=0)then;q=-huge(1.0_dp);return;endif
    if(p>=1)then;q=huge(1.0_dp);return;endif
    sd=sqrt(max(ghyp_var(mu,delta,alpha,beta,lambda),1.0e-12_dp));lo=mu-12*sd;hi=mu+12*sd
    do i=1,45
      mid=(lo+hi)/2
      if(pghyp(mid,mu,delta,alpha,beta,lambda)<p)then;lo=mid;else;hi=mid;endif
    enddo
    q=(lo+hi)/2
  end function

  subroutine rghyp(x,mu,delta,alpha,beta,lambda)
    real(dp),intent(out)::x(:);real(dp),intent(in)::mu,delta,alpha,beta,lambda
    real(dp),allocatable::w(:);integer::i
    allocate(w(size(x)));call rgig(w,delta**2,alpha**2-beta**2,lambda)
    do i=1,size(x);x(i)=mu+beta*w(i)+sqrt(w(i))*normal_rand();enddo
  end subroutine

  function ghyp_mean(mu,delta,alpha,beta,lambda) result(m)
    real(dp),intent(in)::mu,delta,alpha,beta,lambda;real(dp)::m
    m=mu+beta*gig_moment(1.0_dp,delta**2,alpha**2-beta**2,lambda)
  end function
  function ghyp_var(mu,delta,alpha,beta,lambda) result(v)
    real(dp),intent(in)::mu,delta,alpha,beta,lambda;real(dp)::v,mw,vw
    mw=gig_moment(1.0_dp,delta**2,alpha**2-beta**2,lambda)
    vw=gig_moment(2.0_dp,delta**2,alpha**2-beta**2,lambda)-mw*mw
    v=mw+beta*beta*vw
  end function
  function ghyp_mode(mu,delta,alpha,beta,lambda) result(mode)
    real(dp),intent(in)::mu,delta,alpha,beta,lambda;real(dp)::mode,a,b,x1,x2,f1,f2
    integer::i
    a=mu-5*sqrt(ghyp_var(mu,delta,alpha,beta,lambda));b=mu+5*sqrt(ghyp_var(mu,delta,alpha,beta,lambda))
    do i=1,100
      x1=a+(b-a)/3;x2=b-(b-a)/3;f1=dghyp(x1,mu,delta,alpha,beta,lambda);f2=dghyp(x2,mu,delta,alpha,beta,lambda)
      if(f1<f2)then;a=x1;else;b=x2;endif
    enddo;mode=(a+b)/2
  end function

  subroutine ghyp_change_pars(from,to,param,out)
    integer,intent(in)::from,to;real(dp),intent(in)::param(5);real(dp),intent(out)::out(5)
    real(dp)::mu,delta,alpha,beta,lambda,rho,zeta,xi,chi,abar,bbar,gpi,gamma
    mu=param(1);delta=param(2);lambda=param(5)
    select case(from)
    case(1);alpha=param(3);beta=param(4)
    case(2);rho=param(3);zeta=param(4);alpha=zeta/(delta*sqrt(1-rho*rho));beta=rho*alpha
    case(3);xi=param(3);chi=param(4);alpha=(1-xi*xi)/(delta*xi*sqrt(xi*xi-chi*chi));beta=alpha*chi/xi
    case(4);alpha=param(3)/delta;beta=param(4)/delta
    case(5);gpi=param(3);zeta=param(4);alpha=zeta*sqrt(1+gpi*gpi)/delta;beta=gpi*zeta/delta
    end select
    gamma=sqrt(alpha*alpha-beta*beta)
    select case(to)
    case(1);out=[mu,delta,alpha,beta,lambda]
    case(2);out=[mu,delta,beta/alpha,delta*gamma,lambda]
    case(3);out=[mu,delta,1/sqrt(1+delta*gamma),beta/(alpha*sqrt(1+delta*gamma)),lambda]
    case(4);out=[mu,delta,delta*alpha,delta*beta,lambda]
    case(5);out=[mu,delta,beta/gamma,delta*gamma,lambda]
    end select
  end subroutine
  subroutine ghyp_scale(newmean,newsd,param,out)
    real(dp),intent(in)::newmean,newsd,param(5);real(dp),intent(out)::out(5)
    real(dp)::m,s,a,b
    m=ghyp_mean(param(1),param(2),param(3),param(4),param(5));s=sqrt(ghyp_var(param(1),param(2),param(3),param(4),param(5)))
    a=newsd/s;b=newmean-a*m
    out=[b+a*param(1),a*param(2),param(3)/a,param(4)/a,param(5)]
  end subroutine

  function dhyperb(x,mu,delta,alpha,beta) result(d)
    real(dp),intent(in)::x,mu,delta,alpha,beta
    real(dp)::d
    d=dghyp(x,mu,delta,alpha,beta,1.0_dp)
  end function dhyperb

  function phyperb(x,mu,delta,alpha,beta) result(p)
    real(dp),intent(in)::x,mu,delta,alpha,beta
    real(dp)::p
    p=pghyp(x,mu,delta,alpha,beta,1.0_dp)
  end function phyperb

  function qhyperb(p,mu,delta,alpha,beta) result(q)
    real(dp),intent(in)::p,mu,delta,alpha,beta
    real(dp)::q
    q=qghyp(p,mu,delta,alpha,beta,1.0_dp)
  end function qhyperb

  subroutine rhyperb(x,mu,delta,alpha,beta)
    real(dp),intent(out)::x(:)
    real(dp),intent(in)::mu,delta,alpha,beta
    call rghyp(x,mu,delta,alpha,beta,1.0_dp)
  end subroutine rhyperb

  function dnig(x,mu,delta,alpha,beta) result(d)
    real(dp),intent(in)::x,mu,delta,alpha,beta
    real(dp)::d
    d=dghyp(x,mu,delta,alpha,beta,-0.5_dp)
  end function dnig

  function pnig(x,mu,delta,alpha,beta) result(p)
    real(dp),intent(in)::x,mu,delta,alpha,beta
    real(dp)::p
    p=pghyp(x,mu,delta,alpha,beta,-0.5_dp)
  end function pnig

  function qnig(p,mu,delta,alpha,beta) result(q)
    real(dp),intent(in)::p,mu,delta,alpha,beta
    real(dp)::q
    q=qghyp(p,mu,delta,alpha,beta,-0.5_dp)
  end function qnig

  subroutine rnig(x,mu,delta,alpha,beta)
    real(dp),intent(out)::x(:)
    real(dp),intent(in)::mu,delta,alpha,beta
    call rghyp(x,mu,delta,alpha,beta,-0.5_dp)
  end subroutine rnig
end module ghyp_distribution
