module l1pack_lad
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use l1pack_base, only: dp, sqrt2, normal_quantile_l1, sample_quantile_l1, inverse_matrix
  use l1pack_distributions, only: rlaplace, plaplace
  use fastmatrix_regression, only: ols_result, ols_fit
  implicit none
  private

  type, public :: l1fit_result
    real(dp), allocatable :: coefficients(:), fitted(:), residuals(:)
    real(dp) :: minimum=0.0_dp
    integer :: rank=0, iterations=0, info=0
  end type l1fit_result

  type, public :: lad_result
    real(dp), allocatable :: coefficients(:), fitted(:), residuals(:), weights(:)
    real(dp), allocatable :: cov_beta(:,:)
    integer, allocatable :: basic(:)
    real(dp) :: scale=0.0_dp, sad=0.0_dp, loglik=0.0_dp, lambda=0.0_dp
    integer :: n=0, p=0, rank=0, iterations=0, info=0
    logical :: converged=.false.
    character(len=2) :: method='BR'
  end type lad_result

  public :: l1fit, lad_fit, lad_fit_br, lad_fit_em
  public :: nuisance_vcov, vcov_lad, confint_lad, predict_lad, simulate_lad
  public :: lad_quantile_residuals, lad_deviance, lad_loglik

contains

  pure real(dp) function lad_loglik(scale,n) result(ans)
    real(dp), intent(in) :: scale
    integer, intent(in) :: n
    if(scale<=0.0_dp) then
      ans=-huge(1.0_dp)
    else
      ans=-real(n,dp)*(0.5_dp*log(2.0_dp)+1.0_dp+log(scale))
    end if
  end function lad_loglik

  subroutine l1fit(x,y,res,intercept,tolerance)
    real(dp), intent(in) :: x(:,:),y(:)
    type(l1fit_result), intent(out) :: res
    logical, intent(in), optional :: intercept
    real(dp), intent(in), optional :: tolerance
    logical :: add_intercept
    real(dp) :: tol
    real(dp), allocatable :: xx(:,:),coef(:),resid(:)
    integer :: n,p,rank,iter,info
    real(dp) :: minimum

    n=size(x,1)
    if(size(y)/=n) error stop 'l1fit: x and y have different numbers of observations'
    add_intercept=.true.;if(present(intercept))add_intercept=intercept
    tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
    if(add_intercept) then
      p=size(x,2)+1
      allocate(xx(n,p));xx(:,1)=1.0_dp;xx(:,2:p)=x
    else
      p=size(x,2)
      allocate(xx(n,p));xx=x
    end if
    if(n<=p) error stop 'l1fit: more variables than observations'
    allocate(coef(p),resid(n))
    call l1br_core(xx,y,tol,coef,resid,minimum,iter,rank,info)
    allocate(res%coefficients(p),res%fitted(n),res%residuals(n))
    res%coefficients=coef
    res%fitted=matmul(xx,coef)
    res%residuals=resid
    res%minimum=minimum;res%rank=rank;res%iterations=iter;res%info=info
  end subroutine l1fit

  subroutine lad_fit(x,y,res,method,tol,maxiter,level)
    real(dp), intent(in) :: x(:,:),y(:)
    type(lad_result), intent(out) :: res
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: tol,level
    integer, intent(in), optional :: maxiter
    character(len=2) :: meth
    meth='BR';if(present(method))meth=method(1:min(2,len_trim(method)))
    if(meth=='BR'.or.meth=='br') then
      call lad_fit_br(x,y,res,tol,level)
    else if(meth=='EM'.or.meth=='em') then
      call lad_fit_em(x,y,res,tol,maxiter,level)
    else
      error stop 'lad_fit: method must be BR or EM'
    end if
  end subroutine lad_fit

  subroutine lad_fit_br(x,y,res,tol,level)
    real(dp), intent(in) :: x(:,:),y(:)
    type(lad_result), intent(out) :: res
    real(dp), intent(in), optional :: tol,level
    real(dp) :: tt,lev,minimum
    real(dp), allocatable :: coef(:),resid(:),nonzero(:),xtx(:,:),inv(:,:)
    integer, allocatable :: mask(:)
    integer :: n,p,rank,iter,info,i,m,ier,nb

    n=size(x,1);p=size(x,2)
    if(size(y)/=n.or.n<=p) error stop 'lad_fit_br: invalid dimensions'
    tt=1.0e-7_dp;if(present(tol))tt=tol
    lev=0.95_dp;if(present(level))lev=level
    if(lev<=0.0_dp.or.lev>=1.0_dp) error stop 'lad_fit_br: level must be in (0,1)'
    allocate(coef(p),resid(n))
    call l1br_core(x,y,tt,coef,resid,minimum,iter,rank,info)

    res%n=n;res%p=p;res%method='BR';res%rank=rank;res%iterations=iter;res%info=info
    res%converged=(info==0.or.info==1)
    allocate(res%coefficients(p),res%fitted(n),res%residuals(n),res%weights(n))
    res%coefficients=coef;res%fitted=matmul(x,coef);res%residuals=resid
    res%sad=minimum;res%scale=sqrt2*minimum/real(n,dp);res%loglik=lad_loglik(res%scale,n)
    res%weights=1.0_dp

    nb=count(abs(resid)<=10.0_dp*epsilon(1.0_dp)*max(1.0_dp,minimum))
    allocate(res%basic(nb),mask(n));mask=0
    nb=0
    do i=1,n
      if(abs(resid(i))<=10.0_dp*epsilon(1.0_dp)*max(1.0_dp,minimum)) then
        nb=nb+1;res%basic(nb)=i;mask(i)=1
      end if
    end do
    m=n-count(mask==1)
    if(m>1) then
      allocate(nonzero(m));m=0
      do i=1,n
        if(mask(i)==0) then;m=m+1;nonzero(m)=resid(i);end if
      end do
      res%lambda=nuisance_vcov(nonzero,1.0_dp-lev,.true.)
    else
      res%lambda=0.0_dp
    end if

    allocate(xtx(p,p),inv(p,p),res%cov_beta(p,p))
    xtx=matmul(transpose(x),x)
    call inverse_matrix(xtx,inv,ier)
    if(ier==0) then
      res%cov_beta=res%lambda**2*inv
    else
      res%cov_beta=0.0_dp
    end if
  end subroutine lad_fit_br

  subroutine lad_fit_em(x,y,res,tol,maxiter,level)
    real(dp), intent(in) :: x(:,:),y(:)
    type(lad_result), intent(out) :: res
    real(dp), intent(in), optional :: tol,level
    integer, intent(in), optional :: maxiter
    type(ols_result) :: ols
    real(dp) :: tt,lev,eps,old_sad,new_sad,conv
    real(dp), allocatable :: coef(:),fitted(:),resid(:),weights(:),a(:,:),b(:),delta(:)
    real(dp), allocatable :: xtx(:,:),inv(:,:),nonzero(:)
    integer, allocatable :: mask(:)
    integer :: n,p,mi,iter,i,j,k,info,m,nb

    n=size(x,1);p=size(x,2)
    if(size(y)/=n.or.n<=p) error stop 'lad_fit_em: invalid dimensions'
    tt=1.0e-7_dp;if(present(tol))tt=tol
    lev=0.95_dp;if(present(level))lev=level
    mi=200;if(present(maxiter))mi=maxiter
    call ols_fit(x,y,ols,info)
    if(info/=0) error stop 'lad_fit_em: OLS initialization failed'
    allocate(coef(p),fitted(n),resid(n),weights(n),a(p,p),b(p),delta(p))
    coef=ols%coefficients;fitted=ols%fitted;resid=ols%residuals
    eps=sqrt(epsilon(1.0_dp));old_sad=sum(abs(resid));new_sad=old_sad

    do iter=1,mi
      do i=1,n
        if(abs(resid(i))<eps) then
          weights(i)=1.0_dp
        else
          weights(i)=1.0_dp/abs(resid(i))
        end if
      end do
      a=0.0_dp;b=0.0_dp
      do i=1,n
        do j=1,p
          b(j)=b(j)+weights(i)*x(i,j)*resid(i)
          do k=1,p
            a(j,k)=a(j,k)+weights(i)*x(i,j)*x(i,k)
          end do
        end do
      end do
      call solve_sym(a,b,delta,info)
      if(info/=0) exit
      coef=coef+delta
      fitted=matmul(x,coef);resid=y-fitted
      new_sad=sum(abs(resid))
      conv=abs((new_sad-old_sad)/(new_sad+1.0e-50_dp))
      if(conv<tt) exit
      old_sad=new_sad
    end do

    res%n=n;res%p=p;res%method='EM';res%rank=p;res%iterations=min(iter,mi);res%info=info
    res%converged=(iter<=mi.and.info==0)
    allocate(res%coefficients(p),res%fitted(n),res%residuals(n),res%weights(n))
    res%coefficients=coef;res%fitted=fitted;res%residuals=resid
    res%sad=new_sad;res%scale=sqrt2*new_sad/real(n,dp);res%loglik=lad_loglik(res%scale,n)

    eps=epsilon(1.0_dp)**0.35_dp
    allocate(mask(n));mask=0;nb=0
    do i=1,n
      if(abs(resid(i))<eps*max(res%sad,tiny(1.0_dp))) then
        res%weights(i)=1.0_dp;mask(i)=1;nb=nb+1
      else
        res%weights(i)=eps/abs(resid(i))
      end if
    end do
    allocate(res%basic(nb));k=0
    do i=1,n
      if(mask(i)==1) then;k=k+1;res%basic(k)=i;end if
    end do
    m=n-nb
    if(m>1) then
      allocate(nonzero(m));k=0
      do i=1,n
        if(mask(i)==0) then;k=k+1;nonzero(k)=resid(i);end if
      end do
      res%lambda=nuisance_vcov(nonzero,1.0_dp-lev,.true.)
    else
      res%lambda=0.0_dp
    end if
    allocate(xtx(p,p),inv(p,p),res%cov_beta(p,p))
    xtx=matmul(transpose(x),x)
    call inverse_matrix(xtx,inv,info)
    if(info==0)then;res%cov_beta=res%lambda**2*inv;else;res%cov_beta=0.0_dp;end if
  end subroutine lad_fit_em

  real(dp) function nuisance_vcov(resid,alpha,mckean) result(value)
    real(dp), intent(in) :: resid(:),alpha
    logical, intent(in), optional :: mckean
    logical :: job
    integer :: m,lo,hi
    real(dp) :: z,k
    m=size(resid);job=.true.;if(present(mckean))job=mckean
    if(m<2.or.alpha<=0.0_dp.or.alpha>=1.0_dp) then
      value=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    z=normal_quantile_l1(1.0_dp-alpha)
    if(z<=0.0_dp) then
      value=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    value=0.5_dp*sqrt(real(m,dp))/z
    if(job) then
      k=(real(m+1,dp))/2.0_dp-z*sqrt(real(m,dp)/4.0_dp)
      lo=int(k);hi=m-lo+1
    else
      k=(real(m+1,dp))/2.0_dp-sqrt(real(m,dp));lo=int(k)
      k=(real(m+1,dp))/2.0_dp+sqrt(real(m,dp));hi=int(k)
    end if
    lo=max(1,min(m,lo));hi=max(1,min(m,hi))
    value=value*(sample_quantile_l1(resid,hi)-sample_quantile_l1(resid,lo))
  end function nuisance_vcov

  function vcov_lad(res) result(v)
    type(lad_result), intent(in) :: res
    real(dp), allocatable :: v(:,:)
    allocate(v(res%p,res%p));v=res%cov_beta
  end function vcov_lad

  subroutine confint_lad(res,level,ci)
    type(lad_result), intent(in) :: res
    real(dp), intent(in), optional :: level
    real(dp), intent(out) :: ci(:,:)
    real(dp) :: lev,z,se
    integer :: j
    lev=0.95_dp;if(present(level))lev=level
    if(size(ci,1)/=res%p.or.size(ci,2)/=2) error stop 'confint_lad: bad output shape'
    z=normal_quantile_l1(0.5_dp+0.5_dp*lev)
    do j=1,res%p
      se=sqrt(max(res%cov_beta(j,j),0.0_dp))
      ci(j,1)=res%coefficients(j)-z*se
      ci(j,2)=res%coefficients(j)+z*se
    end do
  end subroutine confint_lad

  function predict_lad(res,newx) result(yhat)
    type(lad_result), intent(in) :: res
    real(dp), intent(in) :: newx(:,:)
    real(dp) :: yhat(size(newx,1))
    yhat=matmul(newx,res%coefficients)
  end function predict_lad

  subroutine simulate_lad(res,nsim,y)
    type(lad_result), intent(in) :: res
    integer, intent(in) :: nsim
    real(dp), intent(out) :: y(:,:)
    integer :: i,j
    if(size(y,1)/=res%n.or.size(y,2)/=nsim) error stop 'simulate_lad: bad output shape'
    do j=1,nsim
      do i=1,res%n
        y(i,j)=res%fitted(i)+rlaplace(0.0_dp,res%scale)
      end do
    end do
  end subroutine simulate_lad

  function lad_quantile_residuals(res,y) result(z)
    type(lad_result), intent(in) :: res
    real(dp), intent(in) :: y(:)
    real(dp) :: z(size(y)),p
    integer :: i
    do i=1,size(y)
      p=plaplace(y(i),res%fitted(i),res%scale)
      z(i)=normal_quantile_l1(p)
    end do
  end function lad_quantile_residuals

  pure real(dp) function lad_deviance(res) result(v)
    type(lad_result), intent(in) :: res
    v=sum(abs(res%residuals))
  end function lad_deviance

  subroutine solve_sym(a,b,x,info)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp) :: aa(size(a,1),size(a,2))
    aa=a
    call solve_linear_local(aa,b,x,info)
  end subroutine solve_sym

  subroutine solve_linear_local(a,b,x,info)
    real(dp), intent(inout) :: a(:,:)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp) :: rhs(size(b)),row(size(b)),fac,piv,tmp
    integer :: n,i,k,ip
    n=size(b);rhs=b;info=0
    do k=1,n-1
      ip=k;piv=abs(a(k,k))
      do i=k+1,n
        if(abs(a(i,k))>piv)then;ip=i;piv=abs(a(i,k));end if
      end do
      if(piv<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a))))then;info=k;x=0.0_dp;return;end if
      if(ip/=k)then
        row=a(k,:);a(k,:)=a(ip,:);a(ip,:)=row
        tmp=rhs(k);rhs(k)=rhs(ip);rhs(ip)=tmp
      end if
      do i=k+1,n
        fac=a(i,k)/a(k,k);a(i,k)=0.0_dp
        a(i,k+1:n)=a(i,k+1:n)-fac*a(k,k+1:n);rhs(i)=rhs(i)-fac*rhs(k)
      end do
    end do
    if(abs(a(n,n))<=epsilon(1.0_dp))then;info=n;x=0.0_dp;return;end if
    x(n)=rhs(n)/a(n,n)
    do i=n-1,1,-1
      x(i)=(rhs(i)-dot_product(a(i,i+1:n),x(i+1:n)))/a(i,i)
    end do
  end subroutine solve_linear_local

  subroutine l1br_core(xmat,y,tol,coef,resid,minimum,iter,rank,info)
    real(dp), intent(in) :: xmat(:,:),y(:),tol
    real(dp), intent(out) :: coef(:),resid(:),minimum
    integer, intent(out) :: iter,rank,info
    integer :: m,n,m1,n1,m2,n2,i,j,k,kount,kr,kl,ienter,iout,l,jsel
    integer, allocatable :: s(:)
    real(dp), allocatable :: a(:,:),b(:)
    real(dp) :: d,sumv,pivot,maxv,minv
    logical :: stage,test
    real(dp), parameter :: big=1.0e75_dp

    m=size(xmat,1);n=size(xmat,2);m2=m+2;n2=n+2;m1=m+1;n1=n+1
    allocate(a(m2,n2),b(m),s(m));a=0.0_dp;b=y;a(1:m,1:n)=xmat
    coef=0.0_dp;resid=0.0_dp
    do j=1,n;a(m2,j)=real(j,dp);end do
    do i=1,m
      a(i,n2)=real(n+i,dp);a(i,n1)=b(i)
      if(b(i)<0.0_dp)a(i,1:n2)=-a(i,1:n2)
    end do
    do j=1,n1;a(m1,j)=sum(a(1:m,j));end do
    stage=.true.;ienter=0;iout=0;kount=0;kr=1;kl=1

70  continue
    maxv=-1.0_dp
    do j=kr,n
      if(abs(a(m2,j))>real(n,dp))cycle
      d=abs(a(m1,j))
      if(d>maxv)then;maxv=d;ienter=j;end if
    end do
    if(ienter<1)then;info=2;goto 350;end if
    if(a(m1,ienter)<0.0_dp)a(:,ienter)=-a(:,ienter)
100 continue
    k=0
    do i=kl,m
      d=a(i,ienter)
      if(d<=tol)cycle
      k=k+1;b(k)=a(i,n1)/d;s(k)=i;test=.true.
    end do
120 continue
    if(k>0)then
      test=.true.
    else
      test=.false.;goto 150
    end if
    minv=big;jsel=1
    do i=1,k
      if(b(i)<minv)then;jsel=i;minv=b(i);iout=s(i);end if
    end do
    b(jsel)=b(k);s(jsel)=s(k);k=k-1
150 continue
    if(test.or.(.not.stage))goto 170
    do i=1,m2
      d=a(i,kr);a(i,kr)=a(i,ienter);a(i,ienter)=d
    end do
    kr=kr+1;goto 260
170 continue
    if(.not.test)then;a(m2,n1)=2.0_dp;goto 350;end if
    pivot=a(iout,ienter)
    if(a(m1,ienter)-2.0_dp*pivot<=tol)goto 200
    do j=kr,n1
      d=a(iout,j);a(m1,j)=a(m1,j)-2.0_dp*d;a(iout,j)=-d
    end do
    a(iout,n2)=-a(iout,n2);goto 120
200 continue
    do j=kr,n1
      if(j/=ienter)a(iout,j)=a(iout,j)/pivot
    end do
    do i=1,m1
      if(i==iout)cycle
      d=a(i,ienter)
      do j=kr,n1
        if(j/=ienter)a(i,j)=a(i,j)-d*a(iout,j)
      end do
    end do
    do i=1,m1
      if(i/=iout)a(i,ienter)=-a(i,ienter)/pivot
    end do
    a(iout,ienter)=1.0_dp/pivot
    d=a(iout,n2);a(iout,n2)=a(m2,ienter);a(m2,ienter)=d
    kount=kount+1
    if(.not.stage)goto 270
    kl=kl+1
    do j=kr,n2
      d=a(iout,j);a(iout,j)=a(kount,j);a(kount,j)=d
    end do
260 continue
    if(kount+kr/=n1)goto 70
    stage=.false.
270 continue
    maxv=-big
    do j=kr,n
      d=a(m1,j)
      if(d<0.0_dp)then
        if(d>-2.0_dp)cycle
        d=-d-2.0_dp
      end if
      if(d>maxv)then;maxv=d;ienter=j;end if
    end do
    if(maxv<=tol)goto 310
    if(a(m1,ienter)>0.0_dp)goto 100
    a(:,ienter)=-a(:,ienter);a(m1,ienter)=a(m1,ienter)-2.0_dp;goto 100
310 continue
    l=kl-1
    do i=1,l
      if(a(i,n1)<0.0_dp)a(i,kr:n2)=-a(i,kr:n2)
    end do
    a(m2,n1)=0.0_dp
    if(kr/=1)goto 350
    do j=1,n
      d=abs(a(m1,j))
      if(d<=tol.or.2.0_dp-d<=tol)goto 350
    end do
    a(m2,n1)=1.0_dp
350 continue
    do i=1,m
      k=int(a(i,n2));d=a(i,n1)
      if(k<=0)then;k=-k;d=-d;end if
      if(i<kl)then
        if(k>=1.and.k<=n)coef(k)=d
      else
        k=k-n
        if(k>=1.and.k<=m)resid(k)=d
      end if
    end do
    a(m2,n2)=real(kount,dp);a(m1,n2)=real(n1-kr,dp)
    sumv=0.0_dp
    do i=kl,m;sumv=sumv+a(i,n1);end do
    a(m1,n1)=sumv
    minimum=a(m1,n1);rank=int(a(m1,n2));info=int(a(m2,n1));iter=int(a(m2,n2))
  end subroutine l1br_core

end module l1pack_lad
