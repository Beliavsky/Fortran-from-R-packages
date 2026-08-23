module gig_distribution
  use gh_math, only: dp, pi, bessel_k
  implicit none
  private
  public :: dgig, pgig, qgig, rgig, gig_mean, gig_var, gig_mode, gig_moment, gig_change_pars
contains
  function dgig(x,chi,psi,lambda) result(d)
    real(dp),intent(in)::x,chi,psi,lambda
    real(dp)::d,omega,c
    if (x<=0 .or. chi<=0 .or. psi<=0) then; d=0; return; end if
    omega=sqrt(chi*psi)
    c=(psi/chi)**(lambda/2.0_dp)/(2.0_dp*bessel_k(lambda,omega))
    d=c*x**(lambda-1.0_dp)*exp(-0.5_dp*(chi/x+psi*x))
  end function

  function pgig(q,chi,psi,lambda) result(p)
    real(dp),intent(in)::q,chi,psi,lambda
    real(dp)::p, lo, hi, h, t, s
    integer::i,n
    if(q<=0)then;p=0;return;endif
    ! integrate after x=exp(t): f(exp(t))*exp(t)
    lo=min(log(q)-30.0_dp, log(sqrt(chi/psi))-30.0_dp)
    hi=log(q); n=500; h=(hi-lo)/real(n,dp); s=0
    do i=0,n
      t=lo+h*i
      if(i==0.or.i==n)then
        s=s+dgig(exp(t),chi,psi,lambda)*exp(t)
      elseif(mod(i,2)==0)then
        s=s+2*dgig(exp(t),chi,psi,lambda)*exp(t)
      else
        s=s+4*dgig(exp(t),chi,psi,lambda)*exp(t)
      endif
    enddo
    p=max(0.0_dp,min(1.0_dp,h*s/3.0_dp))
  end function

  function qgig(p,chi,psi,lambda) result(q)
    real(dp),intent(in)::p,chi,psi,lambda
    real(dp)::q,lo,hi,mid
    integer::i
    if(p<=0)then;q=0;return;endif
    if(p>=1)then;q=huge(1.0_dp);return;endif
    lo=log(sqrt(chi/psi))-20; hi=log(sqrt(chi/psi))+20
    do while(pgig(exp(lo),chi,psi,lambda)>p); lo=lo-5; enddo
    do while(pgig(exp(hi),chi,psi,lambda)<p); hi=hi+5; enddo
    do i=1,45
      mid=0.5_dp*(lo+hi)
      if(pgig(exp(mid),chi,psi,lambda)<p)then;lo=mid;else;hi=mid;endif
    enddo
    q=exp(0.5_dp*(lo+hi))
  end function

  subroutine rgig(x,chi,psi,lambda)
    real(dp),intent(out)::x(:)
    real(dp),intent(in)::chi,psi,lambda
    real(dp)::alpha,beta,m,upper,ym,yp,a,b,c,r1,r2,y
    integer::i
    alpha=sqrt(psi/chi); beta=sqrt(psi*chi)
    m=(lambda-1.0_dp+sqrt((lambda-1.0_dp)**2+beta**2))/beta
    upper=m
    do while(gfun(upper,beta,lambda,m)<=0.0_dp)
      upper=2.0_dp*upper
    enddo
    ym=root_bisect(0.0_dp,m,beta,lambda,m)
    yp=root_bisect(m,upper,beta,lambda,m)
    a=(yp-m)*(yp/m)**(0.5_dp*(lambda-1.0_dp))* &
      exp(-0.25_dp*beta*(yp+1.0_dp/yp-m-1.0_dp/m))
    b=(ym-m)*(ym/m)**(0.5_dp*(lambda-1.0_dp))* &
      exp(-0.25_dp*beta*(ym+1.0_dp/ym-m-1.0_dp/m))
    c=-0.25_dp*beta*(m+1.0_dp/m)+0.5_dp*(lambda-1.0_dp)*log(m)
    do i=1,size(x)
      do
        call random_number(r1); call random_number(r2)
        r1=max(r1,tiny(1.0_dp))
        y=m+a*r2/r1+b*(1.0_dp-r2)/r1
        if(y>0.0_dp)then
          if(-log(r1)>=-0.5_dp*(lambda-1.0_dp)*log(y)+ &
             0.25_dp*beta*(y+1.0_dp/y)+c)exit
        endif
      enddo
      x(i)=y/alpha
    enddo
  end subroutine rgig

  pure function gfun(y,beta,lambda,m) result(g)
    real(dp),intent(in)::y,beta,lambda,m
    real(dp)::g
    g=0.5_dp*beta*y**3-y**2*(0.5_dp*beta*m+lambda+1.0_dp)+ &
      y*((lambda-1.0_dp)*m-0.5_dp*beta)+0.5_dp*beta*m
  end function gfun

  function root_bisect(lo0,hi0,beta,lambda,m) result(r)
    real(dp),intent(in)::lo0,hi0,beta,lambda,m
    real(dp)::r,lo,hi,mid,flo,fmid
    integer::k
    lo=lo0;hi=hi0;flo=gfun(lo,beta,lambda,m)
    do k=1,80
      mid=0.5_dp*(lo+hi);fmid=gfun(mid,beta,lambda,m)
      if(flo*fmid<=0.0_dp)then
        hi=mid
      else
        lo=mid;flo=fmid
      endif
    enddo
    r=0.5_dp*(lo+hi)
  end function root_bisect

  function gig_moment(order,chi,psi,lambda) result(m)
    real(dp),intent(in)::order,chi,psi,lambda
    real(dp)::m,omega
    omega=sqrt(chi*psi)
    m=(chi/psi)**(order/2.0_dp)*bessel_k(lambda+order,omega)/bessel_k(lambda,omega)
  end function
  function gig_mean(chi,psi,lambda) result(m)
    real(dp),intent(in)::chi,psi,lambda; real(dp)::m
    m=gig_moment(1.0_dp,chi,psi,lambda)
  end function
  function gig_var(chi,psi,lambda) result(v)
    real(dp),intent(in)::chi,psi,lambda; real(dp)::v,m
    m=gig_mean(chi,psi,lambda); v=gig_moment(2.0_dp,chi,psi,lambda)-m*m
  end function
  pure function gig_mode(chi,psi,lambda) result(m)
    real(dp),intent(in)::chi,psi,lambda; real(dp)::m
    m=(lambda-1+sqrt((lambda-1)**2+chi*psi))/psi
  end function

  subroutine gig_change_pars(from,to,param,out)
    integer,intent(in)::from,to; real(dp),intent(in)::param(3); real(dp),intent(out)::out(3)
    real(dp)::chi,psi,delta,gamma,alpha,beta,omega,eta,l
    l=param(3)
    select case(from)
    case(1); chi=param(1);psi=param(2)
    case(2); delta=param(1);gamma=param(2);chi=delta**2;psi=gamma**2
    case(3); alpha=param(1);beta=param(2);chi=beta/alpha;psi=alpha*beta
    case(4); omega=param(1);eta=param(2);chi=omega*eta;psi=omega/eta
    end select
    select case(to)
    case(1);out=[chi,psi,l]
    case(2);out=[sqrt(chi),sqrt(psi),l]
    case(3);out=[sqrt(psi/chi),sqrt(chi*psi),l]
    case(4);out=[sqrt(chi*psi),sqrt(chi/psi),l]
    end select
  end subroutine
end module gig_distribution
