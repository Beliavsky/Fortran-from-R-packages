module mev_tailindex
  use mev_kinds, only: dp
  use mev_math, only: sort_descending, sort_ascending, normal_quantile
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: shape_hill, shape_pickands, shape_moment, shape_vries, shape_genjack
  public :: shape_osz, pickands_xu, shape_genquant
  public :: rho_dk, rho_fagh, rho_ghp, rho_gbw
  public :: qweissman, qweissman_ci
contains

  pure real(dp) function shape_hill(xdat,k) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),allocatable::pos(:),y(:)
    integer::i,n
    n=count(xdat>0.0_dp)
    if(k<1.or.k>=n)then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(pos(n),y(n));i=0
    do n=1,size(xdat)
      if(xdat(n)>0.0_dp)then;i=i+1;pos(i)=log(xdat(n));end if
    end do
    call sort_descending(pos,y)
    xi=sum(y(1:k))/real(k,dp)-y(k+1)
  end function shape_hill

  pure real(dp) function shape_vries(xdat,k) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),allocatable::pos(:),y(:)
    real(dp)::m1,m2
    integer::i,j,n
    n=count(xdat>0.0_dp)
    if(k<5.or.k>=n)then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(pos(n),y(n));i=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;i=i+1;pos(i)=log(xdat(j));end if;end do
    call sort_descending(pos,y)
    m1=sum(y(1:k))/real(k,dp)-y(k+1)
    m2=sum((y(1:k)-y(k+1))**2)/real(k,dp)
    xi=0.5_dp*m2/m1
  end function shape_vries

  pure real(dp) function shape_genjack(xdat,k) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp)::h,v
    h=shape_hill(xdat,k);v=shape_vries(xdat,k)
    xi=2.0_dp*v-h
  end function shape_genjack

  pure real(dp) function shape_moment(xdat,k) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),allocatable::pos(:),y(:)
    real(dp)::m1,m2
    integer::i,j,n
    n=count(xdat>0.0_dp)
    if(k<5.or.k>=n)then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(pos(n),y(n));i=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;i=i+1;pos(i)=log(xdat(j));end if;end do
    call sort_descending(pos,y)
    m1=sum(y(1:k))/real(k,dp)-y(k+1)
    m2=sum((y(1:k)-y(k+1))**2)/real(k,dp)
    xi=m1+1.0_dp-0.5_dp/(1.0_dp-m1*m1/m2)
  end function shape_moment

  pure real(dp) function shape_pickands(xdat,k) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),allocatable::pos(:),y(:)
    integer::i,j,n,a,b
    n=count(xdat>0.0_dp)
    if(k<5.or.k>=n)then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(pos(n),y(n));i=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;i=i+1;pos(i)=xdat(j);end if;end do
    call sort_ascending(pos,y)
    a=floor(real(k,dp)/4.0_dp);b=floor(real(k,dp)/2.0_dp)
    if(a<1.or.b<1)then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    xi=(log(y(n-a)-y(n-b))-log(y(n-b)-y(n-k)))/log(2.0_dp)
  end function shape_pickands

  pure real(dp) function shape_osz(xdat,k) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),allocatable::y(:)
    real(dp)::shape0,mcst,sval
    integer::j,n
    n=size(xdat)
    if(k<3.or.k>n)then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(y(n));call sort_descending(xdat,y)
    if(y(1)<=y(2))then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    shape0=(2.0_dp*real(n-1,dp)/real(k-2,dp)-2.0_dp)*log(y(1)-y(2))
    mcst=1.0_dp
    do j=3,n-k+3
      mcst=mcst*real(n-j-k+4,dp)/real(n-j+1,dp)
      if(any(y(1:j-1)<=y(j)))then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
      sval=sum(log(y(1:j-1)-y(j)))
      shape0=shape0+mcst*(2.0_dp*real(n-j+1,dp)/real(k-2,dp)-real(j,dp))*sval
    end do
    xi=shape0*real(k*(k-1)*(k-2),dp)/real(n*(n-1)*(n-k+1),dp)
  end function shape_osz

  pure real(dp) function pickands_xu(xdat,m) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::m
    xi=shape_osz(xdat,m)
  end function pickands_xu

  pure real(dp) function shape_genquant(xdat,k,type,p,weight) result(xi)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    character(len=*),intent(in),optional::type
    real(dp),intent(in),optional::p,weight(:)
    real(dp),allocatable::x(:),xs(:),es(:),w(:)
    real(dp)::pp,num,den,hill
    integer::n,j,kk,idx
    character(len=16)::typ
    n=count(xdat>0.0_dp)
    if(k<5.or.k>=n-1)then;xi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(x(n),xs(n));idx=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;idx=idx+1;x(idx)=xdat(j);end if;end do
    call sort_descending(x,xs);x=xs
    typ='genmean';if(present(type))typ=trim(type)
    pp=0.5_dp;if(present(p))pp=p
    allocate(es(k+1),w(k))
    select case(trim(typ))
    case('genmean')
      do kk=1,k+1
        if(kk>=n)then;es(kk)=log(max(x(kk),tiny(1.0_dp)));cycle;end if
        hill=sum(log(x(1:kk)))/real(kk,dp)-log(x(kk+1))
        es(kk)=log(max(x(kk)*hill,tiny(1.0_dp)))
      end do
    case('genmed')
      do kk=1,k+1
        idx=floor(pp*real(kk,dp))+1
        es(kk)=log(max(x(idx)-x(kk+1),tiny(1.0_dp)))
      end do
    case('trimmean')
      do kk=1,k+1
        idx=floor(pp*real(kk,dp))+1
        es(kk)=log(max(sum(x(idx:kk))/real(kk-idx+1,dp)-x(kk+1),tiny(1.0_dp)))
      end do
    case default
      xi=ieee_value(0.0_dp,ieee_quiet_nan);return
    end select
    if(present(weight))then
      if(size(weight)==1)then;w=weight(1);else if(size(weight)>=k)then;w=weight(1:k);else;w=1.0_dp;end if
    else
      w=1.0_dp
    end if
    num=0.0_dp;den=0.0_dp
    do j=1,k
      num=num+w(j)*(log(real(k+1,dp))-log(real(j,dp)))*(es(j)-es(k+1))
      den=den+w(j)*(log(real(k+1,dp))-log(real(j,dp)))**2
    end do
    xi=num/den
  end function shape_genquant

  pure real(dp) function rho_dk(xdat,k,tau) result(rho)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),intent(in),optional::tau
    real(dp),allocatable::x(:),y(:)
    real(dp)::tt,h1,h2,h3,r
    integer::n,j,i,k1,k2
    tt=0.5_dp;if(present(tau))tt=tau
    n=count(xdat>0.0_dp);if(k<5.or.k>=n.or.tt<=0.0_dp.or.tt>=1.0_dp)then;rho=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(x(n),y(n));i=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;i=i+1;x(i)=log(xdat(j));end if;end do
    call sort_descending(x,y)
    k1=max(1,floor(tt*real(k,dp)));k2=max(1,floor(tt*tt*real(k,dp)))
    h2=sum(y(1:k1))/real(k1,dp)-y(k1+1)
    h1=sum(y(1:k2))/real(k2,dp)-y(k2+1)
    h3=sum(y(1:k))/real(k,dp)-y(k+1)
    r=abs((h1-h2)/(h2-h3))
    if(r<=0.0_dp)then;rho=0.0_dp;else;rho=min(0.0_dp,-log(r)/log(tt));end if
  end function rho_dk

  pure real(dp) function rho_fagh(xdat,k,tau) result(rho)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),intent(in),optional::tau
    real(dp),allocatable::x(:),y(:),z(:)
    real(dp)::tt,m1,m2,m3,w
    integer::n,j,i
    tt=0.0_dp;if(present(tau))tt=tau
    n=count(xdat>0.0_dp);if(k<5.or.k>=n)then;rho=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(x(n),y(n),z(k));i=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;i=i+1;x(i)=log(xdat(j));end if;end do
    call sort_descending(x,y);z=y(1:k)-y(k+1)
    m1=sum(z)/real(k,dp);m2=sum(z**2)/real(k,dp);m3=sum(z**3)/real(k,dp)
    if(abs(tt)>1.0e-14_dp)then
      w=(m1**tt-(0.5_dp*m2)**(0.5_dp*tt))/((0.5_dp*m2)**(0.5_dp*tt)-(m3/6.0_dp)**(tt/3.0_dp))
    else
      w=(log(m1)-0.5_dp*log(0.5_dp*m2))/(0.5_dp*log(0.5_dp*m2)-log(m3/6.0_dp)/3.0_dp)
    end if
    rho=min(0.0_dp,3.0_dp*(w-1.0_dp)/(w-3.0_dp))
  end function rho_fagh

  pure real(dp) function rho_ghp(xdat,k,alpha) result(rho)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),intent(in),optional::alpha
    real(dp),allocatable::x(:),y(:),z(:)
    real(dp)::a,m1,m2,qa1,qa2,sa,b1,b2,lo,hi,mid,fmid,flo
    integer::n,j,i,it
    a=2.0_dp;if(present(alpha))a=alpha
    n=count(xdat>0.0_dp);if(k<5.or.k>=n.or.a<=0.0_dp.or.abs(a-0.5_dp)<1e-12_dp.or.abs(a-1.0_dp)<1e-12_dp)then
      rho=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    allocate(x(n),y(n),z(k));i=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;i=i+1;x(i)=log(xdat(j));end if;end do
    call sort_descending(x,y);z=y(1:k)-y(k+1)
    m1=sum(z)/real(k,dp);m2=sum(z**2)/real(k,dp)
    qa1=(sum(z**(a+1.0_dp))/real(k,dp)-gamma(a+2.0_dp)*m1**(a+1.0_dp))/(m2-2.0_dp*m1*m1)
    qa2=(sum(z**(2.0_dp*a))/real(k,dp)-gamma(2.0_dp*a+1.0_dp)*m1**(2.0_dp*a))/(m2-2.0_dp*m1*m1)
    sa=exp(log(a)+2.0_dp*log(a+1.0_dp)+2.0_dp*log_gamma(a)-log(4.0_dp)-log_gamma(2.0_dp*a))*qa2/(qa1*qa1)
    b1=min((2.0_dp*a-1.0_dp)/(a*a),4.0_dp*(2.0_dp*a-1.0_dp)/(a*(a+1.0_dp)**2))
    b2=max((2.0_dp*a-1.0_dp)/(a*a),4.0_dp*(2.0_dp*a-1.0_dp)/(a*(a+1.0_dp)**2))
    if(sa<=b1.or.sa>=b2)then;rho=0.0_dp;return;end if
    if(abs(a-2.0_dp)<1.0e-12_dp)then
      rho=-(2.0_dp*(3.0_dp*sa-2.0_dp)+sqrt(3.0_dp*sa-2.0_dp))/(3.0_dp-4.0_dp*sa);return
    end if
    lo=-25.0_dp;hi=-1.0e-8_dp;flo=sa_rho(lo,sa,a)
    do it=1,150
      mid=0.5_dp*(lo+hi);fmid=sa_rho(mid,sa,a)
      if(flo*fmid<=0.0_dp)then;hi=mid;else;lo=mid;flo=fmid;end if
    end do
    rho=0.5_dp*(lo+hi)
  contains
    pure real(dp) function sa_rho(r,s,aa) result(v)
      real(dp),intent(in)::r,s,aa
      v=s-r*r*(1.0_dp-(1.0_dp-r)**(2.0_dp*aa)-2.0_dp*aa*r*(1.0_dp-r)**(2.0_dp*aa-1.0_dp)) / &
          (1.0_dp-(1.0_dp-r)**(aa+1.0_dp)-(aa+1.0_dp)*r*(1.0_dp-r)**aa)**2
    end function
  end function rho_ghp

  pure real(dp) function rho_gbw(xdat,k) result(rho)
    real(dp),intent(in)::xdat(:)
    integer,intent(in)::k
    real(dp),allocatable::x(:),y(:),z(:)
    real(dp)::t1,t2,u
    integer::n,j,i
    n=count(xdat>0.0_dp);if(k<1.or.k>=n)then;rho=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    allocate(x(n),y(n),z(k));i=0
    do j=1,size(xdat);if(xdat(j)>0.0_dp)then;i=i+1;x(i)=xdat(j);end if;end do
    call sort_descending(x,y)
    do j=1,k;z(j)=real(j,dp)*(log(y(j+1))-log(y(j)));end do
    t1=0.0_dp;t2=0.0_dp
    do j=1,k
      u=real(j,dp)/real(k+1,dp)
      t1=t1+(-1.0_dp-log(u))*z(j)
      t2=t2+(u-0.5_dp)*z(j)
    end do
    t1=t1/real(k,dp);t2=t2/real(k,dp)
    rho=min(0.0_dp,(4.0_dp*t2+t1)/(2.0_dp*t2+t1))
  end function rho_gbw

  pure real(dp) function qweissman(p,k,n,thresh,shape) result(q)
    real(dp),intent(in)::p,thresh,shape
    integer,intent(in)::k,n
    q=thresh*(real(k,dp)/(real(n,dp)*p))**shape
  end function qweissman

  subroutine qweissman_ci(p,k,n,thresh,shape,level,method,q,lower,upper)
    real(dp),intent(in)::p,thresh,shape,level
    integer,intent(in)::k,n
    character(len=*),intent(in)::method
    real(dp),intent(out)::q,lower,upper
    real(dp)::r,z,ql,qu
    q=qweissman(p,k,n,thresh,shape)
    r=real(k,dp)/(real(n,dp)*p)
    select case(trim(method))
    case('bbw1')
      z=normal_quantile((1.0_dp+level)/2.0_dp)
      lower=q*exp(-shape*log(r)*z/sqrt(real(k,dp)))
      upper=q*exp(shape*log(r)*z/sqrt(real(k,dp)))
    case('bbw2')
      ql=normal_quantile((1.0_dp+level)/2.0_dp)/sqrt(real(k,dp));qu=-ql
      lower=q**(1.0_dp/(1.0_dp+ql))*thresh**(ql/(1.0_dp+ql))
      upper=q**(1.0_dp/(1.0_dp+qu))*thresh**(qu/(1.0_dp+qu))
    case default
      lower=q;upper=q
    end select
  end subroutine qweissman_ci

end module mev_tailindex
