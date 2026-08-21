! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_empirical
  use gb2_kinds, only : dp
  implicit none
  private
  public :: weighted_quantile, weighted_gini, main_emp, robust_weights
contains
  real(dp) function weighted_quantile(x,w,p) result(q)
    real(dp), intent(in) :: x(:),w(:),p
    real(dp), allocatable :: xs(:),ws(:)
    real(dp) :: target,cum
    integer :: i,j,n
    if(size(x)/=size(w) .or. p<0.0_dp .or. p>1.0_dp .or. sum(w)<=0.0_dp) error stop 'weighted_quantile: invalid input'
    n=size(x)
    allocate(xs(n),ws(n))
    xs=x
    ws=w
    do i=2,n
      q=xs(i)
      cum=ws(i)
      j=i-1
      do while(j>=1 .and. xs(j)>q)
      xs(j+1)=xs(j)
      ws(j+1)=ws(j)
      j=j-1
      end do
      xs(j+1)=q
      ws(j+1)=cum
    end do
    if(p<=0.0_dp) then
    q=xs(1)
    return
    end if
    if(p>=1.0_dp) then
    q=xs(n)
    return
    end if
    target=p*sum(ws)
    cum=0.0_dp
    do i=1,n
      cum=cum+ws(i)
      if(cum>=target) then
      q=xs(i)
      return
      end if
    end do
    q=xs(n)
  end function weighted_quantile

  real(dp) function weighted_gini(x,w) result(g)
    real(dp), intent(in) :: x(:),w(:)
    real(dp), allocatable :: xs(:),ws(:)
    real(dp) :: sw,sy,cp,cy,pp,py,pn,yn,tmpw,tmpx
    integer :: i,j,n
    n=size(x)
    if(size(w)/=n .or. any(w<0.0_dp)) error stop 'weighted_gini: invalid input'
    allocate(xs(n),ws(n))
    xs=x
    ws=w
    do i=2,n
      tmpx=xs(i)
      tmpw=ws(i)
      j=i-1
      do while(j>=1 .and. xs(j)>tmpx)
      xs(j+1)=xs(j)
      ws(j+1)=ws(j)
      j=j-1
      end do
      xs(j+1)=tmpx
      ws(j+1)=tmpw
    end do
    sw=sum(ws)
    sy=sum(ws*xs)
    if(sw<=0.0_dp .or. sy<=0.0_dp) then
    g=0.0_dp
    return
    end if
    cp=0.0_dp
    cy=0.0_dp
    pp=0.0_dp
    py=0.0_dp
    g=0.0_dp
    do i=1,n
      cp=cp+ws(i)
      cy=cy+ws(i)*xs(i)
      pn=cp/sw
      yn=cy/sy
      g=g+(yn+py)*(pn-pp)
      pp=pn
      py=yn
    end do
    g=max(0.0_dp,min(1.0_dp,1.0_dp-g))
  end function weighted_gini

  subroutine main_emp(x,w,values,prop)
    real(dp), intent(in) :: x(:),w(:)
    real(dp), intent(out) :: values(6)
    real(dp), intent(in), optional :: prop
    real(dp), allocatable :: xp(:),wp(:)
    real(dp) :: med,meanv,thr,ar,pr,poor_med,bottom,top,sw
    integer :: npoor
    pr=0.6_dp
    if(present(prop)) pr=prop
    sw=sum(w)
    med=weighted_quantile(x,w,0.5_dp)
    meanv=dot_product(w,x)/sw
    thr=pr*med
    ar=sum(pack(w,x<thr))/sw
    npoor=count(x<thr)
    if(npoor>0) then
      allocate(xp(npoor),wp(npoor))
      xp=pack(x,x<thr)
      wp=pack(w,x<thr)
      poor_med=weighted_quantile(xp,wp,0.5_dp)
    else
      poor_med=thr
    end if
    if(thr>0.0_dp) then
    pr=1.0_dp-poor_med/thr
    else
    pr=0.0_dp
    end if
    call quintile_income_totals(x,w,bottom,top)
    if(bottom>0.0_dp) then
    top=top/bottom
    else
    top=huge(1.0_dp)
    end if
    values=[med,meanv,100.0_dp*ar,100.0_dp*pr,top,weighted_gini(x,w)]
  end subroutine main_emp

  subroutine quintile_income_totals(x,w,bottom,top)
    real(dp), intent(in) :: x(:),w(:)
    real(dp), intent(out) :: bottom,top
    real(dp), allocatable :: xs(:),ws(:)
    real(dp) :: sw,cum,lo,hi,frac,tmpx,tmpw
    integer :: i,j,n
    n=size(x)
    allocate(xs(n),ws(n))
    xs=x
    ws=w
    do i=2,n
      tmpx=xs(i)
      tmpw=ws(i)
      j=i-1
      do while(j>=1 .and. xs(j)>tmpx)
      xs(j+1)=xs(j)
      ws(j+1)=ws(j)
      j=j-1
      end do
      xs(j+1)=tmpx
      ws(j+1)=tmpw
    end do
    sw=sum(ws)
    bottom=0.0_dp
    top=0.0_dp
    cum=0.0_dp
    do i=1,n
      lo=cum
      hi=cum+ws(i)
      cum=hi
      frac=max(0.0_dp,min(hi,0.2_dp*sw)-lo)/max(ws(i),tiny(1.0_dp))
      bottom=bottom+frac*ws(i)*xs(i)
      frac=max(0.0_dp,hi-max(lo,0.8_dp*sw))/max(ws(i),tiny(1.0_dp))
      top=top+frac*ws(i)*xs(i)
    end do
  end subroutine quintile_income_totals

  subroutine robust_weights(x,w,c,alpha,corr,adjusted)
    real(dp), intent(in) :: x(:),w(:)
    real(dp), intent(in), optional :: c,alpha
    real(dp), intent(out) :: corr(:),adjusted(:)
    real(dp) :: cc,aa,sw,mlz,vlz,a,b,num,d1,d2
    integer :: i
    if(any(shape(corr)/=[size(x)]) .or. any(shape(adjusted)/=[size(x)])) error stop 'robust_weights: shape mismatch'
    cc=0.01_dp
    if(present(c)) cc=c
    aa=0.001_dp
    if(present(alpha)) aa=alpha
    sw=sum(w)
    mlz=dot_product(w,log(x))/sw
    vlz=dot_product(w,(log(x)-mlz)**2)/sw
    a=acos(-1.0_dp)/sqrt(3.0_dp*vlz)
    b=exp(mlz)
    num=abs(((1.0_dp-aa)/aa)**(1.0_dp/a)-(aa/(1.0_dp-aa))**(1.0_dp/a))
    do i=1,size(x)
      d1=abs(b/x(i)-1.0_dp)
      d2=abs(x(i)/b-1.0_dp)
      corr(i)=max(cc,min(1.0_dp,min(num/max(d1,tiny(1.0_dp)),num/max(d2,tiny(1.0_dp)))))
      adjusted(i)=corr(i)*w(i)
    end do
  end subroutine robust_weights
end module gb2_empirical
