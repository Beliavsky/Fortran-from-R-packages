module bzinb_distributions
  use bzinb_kinds, only : dp
  use bzinb_rng, only : random_gamma, random_poisson, categorical4
  implicit none
  private
  public :: log_poisson_pmf, log_nb_pmf, bp_pmf, bp_logpmf, bzip_a_pmf, bzip_a_logpmf
  public :: bzip_b_pmf, bzip_b_logpmf, bnb_pmf, bnb_logpmf, bzinb_pmf, bzinb_logpmf
  public :: rbp_sample, rbzip_a_sample, rbzip_b_sample, rbnb_sample, rbzinb_sample
  public :: true_correlation, nondropout_weight
contains
  pure real(dp) function log_poisson_pmf(k, lambda) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    if (k < 0 .or. lambda < 0.0_dp) then
      v=-huge(1.0_dp)
    else if (lambda <= tiny(1.0_dp)) then
      if(k==0)then;v=0.0_dp;else;v=-huge(1.0_dp);end if
    else
      v=real(k,dp)*log(lambda)-lambda-log_gamma(real(k+1,dp))
    end if
  end function log_poisson_pmf

  pure real(dp) function log_nb_pmf(k, shape, prob) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: shape, prob
    if(k<0 .or. shape<=0.0_dp .or. prob<=0.0_dp .or. prob>1.0_dp)then
      v=-huge(1.0_dp)
    else
      v=log_gamma(real(k,dp)+shape)-log_gamma(shape)-log_gamma(real(k+1,dp)) + &
        shape*log(prob)+real(k,dp)*log(max(1.0_dp-prob,tiny(1.0_dp)))
    end if
  end function log_nb_pmf

  pure real(dp) function logadd(a,b) result(v)
    real(dp),intent(in)::a,b
    real(dp)::m
    if(a < -0.25_dp*huge(1.0_dp))then;v=b;return;end if
    if(b < -0.25_dp*huge(1.0_dp))then;v=a;return;end if
    m=max(a,b);v=m+log(exp(a-m)+exp(b-m))
  end function logadd

  pure real(dp) function bp_logpmf(x,y,m0,m1,m2) result(v)
    integer,intent(in)::x,y
    real(dp),intent(in)::m0,m1,m2
    integer::s,mm
    real(dp)::lt,ls
    if(x<0.or.y<0.or.m0<0.or.m1<0.or.m2<0)then;v=-huge(1.0_dp);return;end if
    mm=min(x,y);ls=-huge(1.0_dp)
    do s=0,mm
      lt=log_poisson_pmf(s,m0)+log_poisson_pmf(x-s,m1)+log_poisson_pmf(y-s,m2)
      ls=logadd(ls,lt)
    end do
    v=ls
  end function bp_logpmf

  pure real(dp) function bp_pmf(x,y,m0,m1,m2) result(v)
    integer,intent(in)::x,y;real(dp),intent(in)::m0,m1,m2
    real(dp)::l
    l=bp_logpmf(x,y,m0,m1,m2);if(l < log(tiny(1.0_dp)))then;v=0.0_dp;else;v=exp(l);end if
  end function bp_pmf

  pure real(dp) function bzip_a_logpmf(x,y,m0,m1,m2,p) result(v)
    integer,intent(in)::x,y;real(dp),intent(in)::m0,m1,m2,p
    real(dp)::a,b
    a=log(max(1.0_dp-p,tiny(1.0_dp)))+bp_logpmf(x,y,m0,m1,m2);b=-huge(1.0_dp)
    if(x==0.and.y==0.and.p>0.0_dp)b=log(p)
    v=logadd(a,b)
  end function bzip_a_logpmf
  pure real(dp) function bzip_a_pmf(x,y,m0,m1,m2,p) result(v)
    integer,intent(in)::x,y;real(dp),intent(in)::m0,m1,m2,p;v=exp(bzip_a_logpmf(x,y,m0,m1,m2,p))
  end function bzip_a_pmf

  pure real(dp) function bzip_b_logpmf(x,y,m0,m1,m2,p1,p2,p3,p4) result(v)
    integer,intent(in)::x,y;real(dp),intent(in)::m0,m1,m2,p1,p2,p3,p4
    v=-huge(1.0_dp)
    if(p1>0)v=logadd(v,log(p1)+bp_logpmf(x,y,m0,m1,m2))
    if(y==0.and.p2>0)v=logadd(v,log(p2)+log_poisson_pmf(x,m0+m1))
    if(x==0.and.p3>0)v=logadd(v,log(p3)+log_poisson_pmf(y,m0+m2))
    if(x==0.and.y==0.and.p4>0)v=logadd(v,log(p4))
  end function bzip_b_logpmf
  pure real(dp) function bzip_b_pmf(x,y,m0,m1,m2,p1,p2,p3,p4) result(v)
    integer,intent(in)::x,y;real(dp),intent(in)::m0,m1,m2,p1,p2,p3,p4
    v=exp(bzip_b_logpmf(x,y,m0,m1,m2,p1,p2,p3,p4))
  end function bzip_b_pmf

  pure real(dp) function bnb_logpmf(x,y,a0,a1,a2,b1,b2) result(v)
    integer,intent(in)::x,y
    real(dp),intent(in)::a0,a1,a2,b1,b2
    integer::k,m
    real(dp)::t1,t2,lt,ls,c
    if(x<0.or.y<0.or.min(a0,a1,a2,b1,b2)<=0.0_dp)then;v=-huge(1.0_dp);return;end if
    t1=(b1+b2+1.0_dp)/(b1+1.0_dp);t2=(b1+b2+1.0_dp)/(b2+1.0_dp)
    ls=-huge(1.0_dp)
    do m=0,y
      do k=0,x
        lt=log_gamma(a1+real(k,dp))-log_gamma(real(k+1,dp))-log_gamma(a1) + &
          log_gamma(real(x+y-m-k,dp)+a0)-log_gamma(real(x-k+1,dp))-log_gamma(a0+real(y-m,dp)) + &
          log_gamma(real(m,dp)+a2)-log_gamma(real(m+1,dp))-log_gamma(a2) + &
          log_gamma(real(y-m,dp)+a0)-log_gamma(real(y-m+1,dp))-log_gamma(a0) + &
          real(k,dp)*log(t1)+real(m,dp)*log(t2)
        ls=logadd(ls,lt)
      end do
    end do
    c=-(real(x+y,dp)+a0)*log(1.0_dp+b1+b2)+real(x,dp)*log(b1)+real(y,dp)*log(b2) - &
      a1*log(1.0_dp+b1)-a2*log(1.0_dp+b2)
    v=ls+c
  end function bnb_logpmf
  pure real(dp) function bnb_pmf(x,y,a0,a1,a2,b1,b2) result(v)
    integer,intent(in)::x,y;real(dp),intent(in)::a0,a1,a2,b1,b2
    real(dp)::l;l=bnb_logpmf(x,y,a0,a1,a2,b1,b2);if(l<log(tiny(1.0_dp)))then;v=0;else;v=exp(l);end if
  end function bnb_pmf

  pure real(dp) function bzinb_logpmf(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4) result(v)
    integer,intent(in)::x,y
    real(dp),intent(in)::a0,a1,a2,b1,b2,p1,p2,p3,p4
    real(dp) :: l
    v=-huge(1.0_dp)
    if(p1>0)v=logadd(v,log(p1)+bnb_logpmf(x,y,a0,a1,a2,b1,b2))
    if(y==0.and.p2>0)then;l=log_nb_pmf(x,a0+a1,1.0_dp/(1.0_dp+b1));v=logadd(v,log(p2)+l);end if
    if(x==0.and.p3>0)then;l=log_nb_pmf(y,a0+a2,1.0_dp/(1.0_dp+b2));v=logadd(v,log(p3)+l);end if
    if(x==0.and.y==0.and.p4>0)v=logadd(v,log(p4))
  end function bzinb_logpmf
  pure real(dp) function bzinb_pmf(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4) result(v)
    integer,intent(in)::x,y;real(dp),intent(in)::a0,a1,a2,b1,b2,p1,p2,p3,p4
    v=exp(bzinb_logpmf(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4))
  end function bzinb_pmf

  pure real(dp) function true_correlation(a0,a1,a2,b1,b2) result(rho)
    real(dp),intent(in)::a0,a1,a2,b1,b2
    rho=a0*sqrt(b1*b2/((b1+1.0_dp)*(b2+1.0_dp)*(a0+a1)*(a0+a2)))
  end function true_correlation

  pure real(dp) function nondropout_weight(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4) result(w)
    integer,intent(in)::x,y;real(dp),intent(in)::a0,a1,a2,b1,b2,p1,p2,p3,p4
    real(dp)::lnum,lden
    lnum=log(max(p1,tiny(1.0_dp)))+bnb_logpmf(x,y,a0,a1,a2,b1,b2)
    lden=bzinb_logpmf(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4)
    w=exp(lnum-lden)
  end function nondropout_weight

  subroutine rbp_sample(n,m0,m1,m2,x,y)
    integer,intent(in)::n;real(dp),intent(in)::m0,m1,m2;integer,intent(out)::x(n),y(n)
    integer::i,u,v,w
    do i=1,n;u=random_poisson(m0);v=random_poisson(m1);w=random_poisson(m2);x(i)=u+v;y(i)=u+w;end do
  end subroutine rbp_sample

  subroutine rbzip_a_sample(n,m0,m1,m2,p,x,y)
    integer,intent(in)::n;real(dp),intent(in)::m0,m1,m2,p;integer,intent(out)::x(n),y(n)
    integer::i;real(dp)::u
    call rbp_sample(n,m0,m1,m2,x,y)
    do i=1,n;call random_number(u);if(u<p)then;x(i)=0;y(i)=0;end if;end do
  end subroutine rbzip_a_sample

  subroutine rbzip_b_sample(n,m0,m1,m2,p,x,y)
    integer,intent(in)::n;real(dp),intent(in)::m0,m1,m2,p(4);integer,intent(out)::x(n),y(n)
    integer::i,e,u,v,w
    do i=1,n
      u=random_poisson(m0);v=random_poisson(m1);w=random_poisson(m2);e=categorical4(p)
      select case(e);case(1);x(i)=u+v;y(i)=u+w;case(2);x(i)=u+v;y(i)=0;case(3);x(i)=0;y(i)=u+w;case default;x(i)=0;y(i)=0;end select
    end do
  end subroutine rbzip_b_sample

  subroutine rbnb_sample(n,a0,a1,a2,b1,b2,x,y)
    integer,intent(in)::n;real(dp),intent(in)::a0,a1,a2,b1,b2;integer,intent(out)::x(n),y(n)
    integer::i;real(dp)::r0,r1,r2
    do i=1,n
      r0=random_gamma(a0,b1);r1=random_gamma(a1,b1);r2=random_gamma(a2,b1)
      x(i)=random_poisson(r0+r1);y(i)=random_poisson((r0+r2)*b2/b1)
    end do
  end subroutine rbnb_sample

  subroutine rbzinb_sample(n,a0,a1,a2,b1,b2,p,x,y)
    integer,intent(in)::n;real(dp),intent(in)::a0,a1,a2,b1,b2,p(4);integer,intent(out)::x(n),y(n)
    integer::i,e;integer::xx(1),yy(1)
    do i=1,n
      call rbnb_sample(1,a0,a1,a2,b1,b2,xx,yy);e=categorical4(p)
      select case(e)
      case(1)
        x(i)=xx(1); y(i)=yy(1)
      case(2)
        x(i)=xx(1); y(i)=0
      case(3)
        x(i)=0; y(i)=yy(1)
      case default
        x(i)=0; y(i)=0
      end select
    end do
  end subroutine rbzinb_sample
end module bzinb_distributions
