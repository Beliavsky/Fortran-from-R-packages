! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_quantiles
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, rep_design_t, quantile_result_t, &
    QRULE_MATH,QRULE_SCHOOL,QRULE_SHAHVAISH,QRULE_HF1,QRULE_HF2,QRULE_HF3,QRULE_HF4,QRULE_HF5,QRULE_HF6,QRULE_HF7,QRULE_HF8,QRULE_HF9
  use survey_estimators, only : svy_mean
  use survey_replicates, only : svr_var
  use survey_design, only : design_degf
  implicit none
  private
  public :: weighted_quantile, svy_quantile, rep_quantile, normal_quantile, student_t_quantile
contains

  real(dp) function weighted_quantile(x,w,p,rule) result(q)
    real(dp), intent(in) :: x(:),w(:),p
    integer, intent(in), optional :: rule
    real(dp), allocatable :: xs(:),ws(:),cum(:),pk(:),xu(:),wu(:)
    real(dp) :: total,wlow,wup,gamma,wbar
    integer :: r,n,pos,posnext,m
    logical, allocatable :: keep(:)
    if(size(x)/=size(w)) error stop 'weighted_quantile: size mismatch'
    if(p<0 .or. p>1) error stop 'weighted_quantile: p outside [0,1]'
    r=QRULE_MATH; if(present(rule)) r=rule
    keep=abs(w)>tiny(1.0_dp); n=count(keep); if(n==0) error stop 'weighted_quantile: all weights zero'
    allocate(xs(n),ws(n)); xs=pack(x,keep); ws=pack(w,keep); call sort_pairs(xs,ws)
    if(n==1) then; q=xs(1); return; end if
    total=sum(ws); if(total<=0) error stop 'weighted_quantile: nonpositive total weight'
    select case(r)
    case(QRULE_MATH,QRULE_HF1,QRULE_SCHOOL,QRULE_HF2,QRULE_HF4)
      allocate(cum(n)); cum=cumsum(ws); pos=last_le(cum,p*total); posnext=min(n,pos+1)
      wlow=p-cum(pos)/total; wup=cum(posnext)/total-p
      select case(r)
      case(QRULE_MATH,QRULE_HF1)
        q=merge(xs(pos),xs(posnext),wlow<=0)
      case(QRULE_SCHOOL,QRULE_HF2)
        if(wlow<=0) then; q=(xs(pos)+xs(posnext))/2; else; q=xs(posnext); end if
      case(QRULE_HF4)
        if(abs(wup+wlow)<=tiny(1.0_dp)) then; q=xs(pos); else; gamma=wlow/(wup+wlow); q=xs(pos)*(1-gamma)+xs(posnext)*gamma; end if
      end select
    case(QRULE_HF3)
      call aggregate_ties(xs,ws,xu,wu); m=size(xu); allocate(cum(m)); cum=cumsum(wu); pos=last_le(cum,p*sum(wu)); posnext=min(m,pos+1)
      wlow=p-cum(pos)/sum(wu); if(wlow<=0 .and. mod(pos,2)==0) then; q=xu(pos); else; q=xu(posnext); end if
    case(QRULE_HF5)
      allocate(cum(n),pk(n)); cum=cumsum(ws); pk=(cum-ws/2)/cum(n); q=linear_approx(pk,xs,p)
    case(QRULE_HF6)
      allocate(cum(n),pk(n)); cum=cumsum(ws); pk=cum/(cum(n)+ws(n)); q=linear_approx(pk,xs,p)
    case(QRULE_SHAHVAISH)
      allocate(cum(n),pk(n)); wbar=sum(ws)/real(n,dp); cum=cumsum(ws/wbar); pk=(cum+0.5_dp-ws/(2*wbar))/real(n+1,dp); q=constant_approx(pk,xs,p)
    case(QRULE_HF7)
      allocate(cum(n),pk(n)); cum=cumsum(ws); pk(1)=0
      if(n>1) pk(2:n)=cum(1:n-1)/cum(n-1); q=linear_approx(pk,xs,p)
    case(QRULE_HF8)
      allocate(cum(n),pk(n)); cum=cumsum(ws); pk(1)=(2.0_dp/3.0_dp)*cum(1)/(cum(n)+ws(n)/3)
      if(n>1) pk(2:n)=((cum(1:n-1))/3+2*cum(2:n)/3)/(cum(n)+ws(n)/3); q=linear_approx(pk,xs,p)
    case(QRULE_HF9)
      allocate(cum(n),pk(n)); cum=cumsum(ws); pk(1)=(5.0_dp/8.0_dp)*cum(1)/(cum(n)+ws(n)/4)
      if(n>1) pk(2:n)=((3.0_dp/8.0_dp)*cum(1:n-1)+(5.0_dp/8.0_dp)*cum(2:n))/(cum(n)+ws(n)/4); q=linear_approx(pk,xs,p)
    case default
      error stop 'weighted_quantile: unknown rule'
    end select
  end function weighted_quantile

  function svy_quantile(x,probs,design,rule,alpha) result(ans)
    real(dp), intent(in) :: x(:),probs(:)
    type(survey_design_t), intent(in) :: design
    integer, intent(in), optional :: rule
    real(dp), intent(in), optional :: alpha
    type(quantile_result_t) :: ans
    real(dp) :: a,qhat,phat,sep,crit,pl,pu
    integer :: j,r,df
    ! local workaround avoided: directly use svy_mean result below
    if(size(x)/=design%n) error stop 'svy_quantile: size mismatch'
    r=QRULE_MATH; if(present(rule)) r=rule; a=0.05_dp; if(present(alpha)) a=alpha
    allocate(ans%quantile(size(probs)),ans%se(size(probs)),ans%lower(size(probs)),ans%upper(size(probs)))
    df=design_degf(design); crit=student_t_quantile(1-a/2,max(df,1))
    do j=1,size(probs)
      qhat=weighted_quantile(x,design%weight,probs(j),r); ans%quantile(j)=qhat
      call cdf_point_uncertainty(x,qhat,design,phat,sep)
      pl=max(0.0_dp,phat-crit*sep); pu=min(1.0_dp,phat+crit*sep)
      ans%lower(j)=weighted_quantile(x,design%weight,pl,r); ans%upper(j)=weighted_quantile(x,design%weight,pu,r)
      ans%se(j)=(ans%upper(j)-ans%lower(j))/(2*crit)
    end do
  end function svy_quantile

  function rep_quantile(x,probs,design,rule,alpha) result(ans)
    real(dp), intent(in) :: x(:),probs(:)
    type(rep_design_t), intent(in) :: design
    integer, intent(in), optional :: rule
    real(dp), intent(in), optional :: alpha
    type(quantile_result_t) :: ans
    real(dp), allocatable :: theta(:,:), coef(:), vv(:,:)
    real(dp) :: a,crit
    integer :: j,r,qr
    if(size(x)/=design%n) error stop 'rep_quantile: size mismatch'
    qr=QRULE_MATH; if(present(rule)) qr=rule; a=0.05_dp; if(present(alpha)) a=alpha; crit=normal_quantile(1-a/2)
    allocate(ans%quantile(size(probs)),ans%se(size(probs)),ans%lower(size(probs)),ans%upper(size(probs)), &
             theta(design%r,size(probs)),coef(size(probs)),vv(size(probs),size(probs)))
    do j=1,size(probs); coef(j)=weighted_quantile(x,design%weight,probs(j),qr); end do
    do r=1,design%r; do j=1,size(probs); theta(r,j)=weighted_quantile(x,design%repweights(:,r),probs(j),qr); end do; end do
    vv=svr_var(theta,design%scale,design%rscales,design%mse,coef)
    ans%quantile=coef
    do j=1,size(probs); ans%se(j)=sqrt(max(0.0_dp,vv(j,j))); ans%lower(j)=coef(j)-crit*ans%se(j); ans%upper(j)=coef(j)+crit*ans%se(j); end do
  end function rep_quantile

  subroutine cdf_point_uncertainty(x,q,design,p,se)
    use survey_types, only : svystat_t
    real(dp),intent(in)::x(:),q; type(survey_design_t),intent(in)::design; real(dp),intent(out)::p,se
    real(dp),allocatable::z(:,:); type(svystat_t)::m
    allocate(z(size(x),1)); z(:,1)=merge(1.0_dp,0.0_dp,x<=q); m=svy_mean(z,design); p=m%estimate(1); se=sqrt(max(0.0_dp,m%variance(1,1)))
  end subroutine cdf_point_uncertainty

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a1=-3.969683028665376e1_dp,a2=2.209460984245205e2_dp,a3=-2.759285104469687e2_dp, &
      a4=1.383577518672690e2_dp,a5=-3.066479806614716e1_dp,a6=2.506628277459239_dp, &
      b1=-5.447609879822406e1_dp,b2=1.615858368580409e2_dp,b3=-1.556989798598866e2_dp,b4=6.680131188771972e1_dp,b5=-1.328068155288572e1_dp, &
      c1=-7.784894002430293e-3_dp,c2=-3.223964580411365e-1_dp,c3=-2.400758277161838_dp,c4=-2.549732539343734_dp,c5=4.374664141464968_dp,c6=2.938163982698783_dp, &
      d1=7.784695709041462e-3_dp,d2=3.224671290700398e-1_dp,d3=2.445134137142996_dp,d4=3.754408661907416_dp
    real(dp)::q,r
    if(p<=0) then; x=-huge(1.0_dp); return; else if(p>=1) then; x=huge(1.0_dp); return; end if
    if(p<0.02425_dp) then
      q=sqrt(-2*log(p)); x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
    else if(p>0.97575_dp) then
      q=sqrt(-2*log(1-p)); x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
    else
      q=p-0.5_dp; r=q*q; x=((((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q)/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
    end if
  end function normal_quantile

  pure real(dp) function student_t_quantile(p,df) result(t)
    real(dp),intent(in)::p;integer,intent(in)::df
    real(dp)::z,v,z2,z3,z5,z7,z9
    if(df<=0) then; t=normal_quantile(p); return; end if
    z=normal_quantile(p); v=real(df,dp); z2=z*z;z3=z*z2;z5=z3*z2;z7=z5*z2;z9=z7*z2
    ! Cornish-Fisher expansion, highly accurate for survey df except at very tiny df.
    t=z+(z3+z)/(4*v)+(5*z5+16*z3+3*z)/(96*v*v)+(3*z7+19*z5+17*z3-15*z)/(384*v**3) &
      +(79*z9+776*z7+1482*z5-1920*z3-945*z)/(92160*v**4)
  end function student_t_quantile

  subroutine sort_pairs(x,w)
    real(dp),intent(inout)::x(:),w(:);integer::i,j;real(dp)::kx,kw
    do i=2,size(x);kx=x(i);kw=w(i);j=i-1;do while(j>=1 .and. x(j)>kx);x(j+1)=x(j);w(j+1)=w(j);j=j-1;end do;x(j+1)=kx;w(j+1)=kw;end do
  end subroutine sort_pairs
  function cumsum(x) result(y)
    real(dp),intent(in)::x(:);real(dp)::y(size(x));integer::i;y(1)=x(1);do i=2,size(x);y(i)=y(i-1)+x(i);end do
  end function cumsum
  integer function last_le(cum,target) result(pos)
    real(dp),intent(in)::cum(:),target;integer::i;pos=1;do i=1,size(cum);if(cum(i)<=target)pos=i;end do
  end function last_le
  pure real(dp) function linear_approx(px,y,p) result(q)
    real(dp),intent(in)::px(:),y(:),p;integer::i
    if(p<=px(1))then;q=y(1);return;else if(p>=px(size(px)))then;q=y(size(y));return;end if
    do i=1,size(px)-1;if(p>=px(i).and.p<=px(i+1))then;if(abs(px(i+1)-px(i))<=epsilon(1.0_dp)*max(1.0_dp,abs(px(i))))then;q=y(i+1);else;q=y(i)+(y(i+1)-y(i))*(p-px(i))/(px(i+1)-px(i));end if;return;end if;end do;q=y(size(y))
  end function linear_approx
  pure real(dp) function constant_approx(px,y,p) result(q)
    real(dp),intent(in)::px(:),y(:),p;integer::i
    if(p<=px(1))then;q=y(1);return;end if
    q=y(size(y));do i=1,size(px)-1;if(p<px(i+1))then;q=y(i);return;end if;end do
  end function constant_approx
  subroutine aggregate_ties(x,w,xu,wu)
    real(dp),intent(in)::x(:),w(:);real(dp),allocatable,intent(out)::xu(:),wu(:);real(dp),allocatable::tx(:),tw(:);integer::i,m
    allocate(tx(size(x)),tw(size(x)));m=0;do i=1,size(x);if(m==0.or.abs(x(i)-tx(m))>epsilon(1.0_dp)*max(1.0_dp,abs(x(i)),abs(tx(m))))then;m=m+1;tx(m)=x(i);tw(m)=w(i);else;tw(m)=tw(m)+w(i);end if;end do;allocate(xu(m),wu(m));xu=tx(1:m);wu=tw(1:m)
  end subroutine aggregate_ties

end module survey_quantiles
