! SPDX-License-Identifier: Artistic-2.0
module ecd_fitting
  use ecd_kinds, only : dp, ecd_ok, ecd_invalid, ecd_no_convergence
  use ecd_math, only : nan_dp
  use ecd_core, only : ecd_model, ecd_new, ecd_pdf, ecd_statistics
  use ecld_models, only : ecld_model, ecld_new, ecld_from_sd, ecld_const, ecld_solve
  use ecd_processes, only : sld_model, sld_new, sld_pdf, sld_cumulants
  use ecd_timeseries, only : sample_statistics, sample_stats
  implicit none
  private

  abstract interface
    function objective_vector(x) result(v)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: v
    end function objective_vector
  end interface

  type, public :: optimization_result
    real(dp), allocatable :: parameters(:)
    real(dp) :: objective=huge(1.0_dp)
    integer :: iterations=0, evaluations=0, status=ecd_no_convergence
  end type optimization_result

  type, public :: ecd_fit_result
    type(ecd_model) :: model
    type(optimization_result) :: optimization
  end type ecd_fit_result

  type, public :: ecld_fit_result
    type(ecld_model) :: model
    type(optimization_result) :: optimization
  end type ecld_fit_result

  type, public :: sld_fit_result
    type(sld_model) :: model
    type(optimization_result) :: optimization
  end type sld_fit_result

  interface fit_ecld_mle
    module procedure fit_ecld_moments
  end interface fit_ecld_mle

  interface ecd_standardfit
    module procedure fit_ecd_mle
  end interface ecd_standardfit

  interface qsld_fit
    module procedure fit_sld_mle
  end interface qsld_fit

  public :: nelder_mead, fit_ecd_mle, fit_ecld_moments, fit_ecld_mle, fit_sld_mle
  public :: ecd_standardfit, qsld_fit, ecd_estimate_const

contains

  function ecd_estimate_const(d) result(c)
    type(ecd_model), intent(in) :: d
    real(dp) :: c,r,y0,var3
    r=sqrt(d%alpha*d%alpha+d%gamma*d%gamma)
    y0=0.0_dp
    if(r>0.0_dp) y0=-r**(1.0_dp/3.0_dp)
    var3=1.5_dp*r**(2.0_dp/3.0_dp)-4.5_dp*r**(1.0_dp/3.0_dp)+63.0_dp/8.0_dp
    c=2.0_dp*sqrt(acos(-1.0_dp)/2.0_dp)*exp(y0)*sqrt(max(var3,0.0_dp))*d%sigma
  end function ecd_estimate_const

  subroutine nelder_mead(fun,x0,result,lower,upper,step,max_iter,tolerance)
    procedure(objective_vector) :: fun
    real(dp), intent(in) :: x0(:)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:),upper(:),step(:),tolerance
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: simplex(:,:),f(:),centroid(:),xr(:),xe(:),xc(:),st(:)
    real(dp) :: tol,spread,alpha,gamma_v,rho,sigma_v
    integer :: n,mx,i,j,best,worst,second,evaluations
    n=size(x0); mx=1000; tol=1.0e-8_dp
    if(present(max_iter))mx=max_iter
    if(present(tolerance))tol=tolerance
    allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),st(n))
    st=max(0.05_dp*max(abs(x0),1.0_dp),1.0e-3_dp)
    if(present(step))st=step
    simplex(:,1)=project(x0)
    do j=1,n
      simplex(:,j+1)=simplex(:,1); simplex(j,j+1)=simplex(j,j+1)+st(j)
      simplex(:,j+1)=project(simplex(:,j+1))
    end do
    evaluations=0
    do j=1,n+1; f(j)=safe_fun(simplex(:,j)); evaluations=evaluations+1; end do
    alpha=1.0_dp; gamma_v=2.0_dp; rho=0.5_dp; sigma_v=0.5_dp
    result%status=ecd_no_convergence
    do i=1,mx
      call order_indices(f,best,worst,second)
      spread=maxval(abs(f-f(best)))/(1.0_dp+abs(f(best)))
      if(spread<tol .and. maxval(abs(simplex-spread_cols(simplex(:,best),n+1)))<sqrt(tol)*(1+maxval(abs(simplex(:,best))))) then
        result%status=ecd_ok; exit
      end if
      centroid=(sum(simplex,dim=2)-simplex(:,worst))/real(n,dp)
      xr=project(centroid+alpha*(centroid-simplex(:,worst)))
      if(safe_fun(xr)<f(best)) then
        xe=project(centroid+gamma_v*(xr-centroid))
        if(safe_fun(xe)<safe_fun(xr)) then; simplex(:,worst)=xe; f(worst)=safe_fun(xe); evaluations=evaluations+3
        else; simplex(:,worst)=xr; f(worst)=safe_fun(xr); evaluations=evaluations+3; end if
      else if(safe_fun(xr)<f(second)) then
        simplex(:,worst)=xr; f(worst)=safe_fun(xr); evaluations=evaluations+2
      else
        if(safe_fun(xr)<f(worst)) then; xc=project(centroid+rho*(xr-centroid))
        else; xc=project(centroid-rho*(centroid-simplex(:,worst))); end if
        if(safe_fun(xc)<min(f(worst),safe_fun(xr))) then
          simplex(:,worst)=xc; f(worst)=safe_fun(xc); evaluations=evaluations+3
        else
          do j=1,n+1
            if(j/=best) then
              simplex(:,j)=project(simplex(:,best)+sigma_v*(simplex(:,j)-simplex(:,best)))
              f(j)=safe_fun(simplex(:,j)); evaluations=evaluations+1
            end if
          end do
        end if
      end if
      result%iterations=i
    end do
    call order_indices(f,best,worst,second)
    allocate(result%parameters(n)); result%parameters=simplex(:,best)
    result%objective=f(best); result%evaluations=evaluations
  contains
    function project(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y=x
      if(present(lower))y=max(y,lower)
      if(present(upper))y=min(y,upper)
    end function project
    function safe_fun(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v
      v=fun(x)
      if(.not.(v<huge(1.0_dp)))v=huge(1.0_dp)/100.0_dp
    end function safe_fun
  end subroutine nelder_mead

  subroutine order_indices(f,best,worst,second)
    real(dp), intent(in) :: f(:)
    integer, intent(out) :: best,worst,second
    integer :: i
    best=minloc(f,dim=1); worst=maxloc(f,dim=1); second=best
    do i=1,size(f)
      if(i/=worst .and. (second==best .or. f(i)>f(second)))second=i
    end do
  end subroutine order_indices

  pure function spread_cols(x,ncol) result(a)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: ncol
    real(dp) :: a(size(x),ncol)
    integer :: j
    do j=1,ncol;a(:,j)=x;end do
  end function spread_cols

  function fit_ecd_mle(data,initial,max_iter) result(fit)
    real(dp), intent(in) :: data(:)
    type(ecd_model), intent(in), optional :: initial
    integer, intent(in), optional :: max_iter
    type(ecd_fit_result) :: fit
    type(ecd_model) :: d0
    type(sample_statistics) :: s
    real(dp) :: x0(5),lo(5),hi(5)
    integer :: mx
    s=sample_stats(data); mx=300; if(present(max_iter))mx=max_iter
    if(present(initial)) then; d0=initial
    else; d0=ecd_new(alpha=0.0_dp,gamma=0.0_dp,sigma=max(s%sd,1e-4_dp),beta=0.0_dp,mu=s%mean); end if
    x0=[d0%mu,log(max(d0%sigma,1e-8_dp)),d0%alpha,d0%gamma,d0%beta]
    lo=[s%mean-5*s%sd,log(max(s%sd*0.02_dp,1e-8_dp)),-100.0_dp,-100.0_dp,-9.99_dp]
    hi=[s%mean+5*s%sd,log(max(s%sd*20.0_dp,1e-6_dp)),100.0_dp,100.0_dp,9.99_dp]
    call nelder_mead(obj,x0,fit%optimization,lo,hi,max_iter=mx,tolerance=1e-7_dp)
    fit%model=ecd_new(alpha=fit%optimization%parameters(3),gamma=fit%optimization%parameters(4), &
      sigma=exp(fit%optimization%parameters(2)),beta=fit%optimization%parameters(5), &
      mu=fit%optimization%parameters(1))
  contains
    function obj(p) result(v)
      real(dp), intent(in) :: p(:)
      real(dp) :: v,dens
      type(ecd_model) :: d
      integer :: i
      d=ecd_new(alpha=p(3),gamma=p(4),sigma=exp(p(2)),beta=p(5),mu=p(1))
      v=0.0_dp
      do i=1,size(data)
        dens=ecd_pdf(d,data(i))
        if(dens<=tiny(1.0_dp)) then; v=v+700.0_dp; else; v=v-log(dens); end if
      end do
      v=v/real(size(data),dp)
    end function obj
  end function fit_ecd_mle

  function fit_ecld_moments(data,initial,max_iter) result(fit)
    real(dp), intent(in) :: data(:)
    type(ecld_model), intent(in), optional :: initial
    integer, intent(in), optional :: max_iter
    type(ecld_fit_result) :: fit
    type(ecld_model) :: d0
    type(sample_statistics) :: s
    real(dp) :: x0(4),lo(4),hi(4)
    integer :: mx
    s=sample_stats(data); mx=500; if(present(max_iter))mx=max_iter
    if(present(initial)) then; d0=initial
    else; d0=ecld_from_sd(3.0_dp,max(s%sd,1e-6_dp),0.0_dp,s%mean); end if
    x0=[d0%lambda,d0%beta,d0%mu,log(max(d0%sigma,1e-8_dp))]
    lo=[0.5_dp,-5.0_dp,s%mean-5*s%sd,log(max(s%sd*0.005_dp,1e-10_dp))]
    hi=[10.0_dp,5.0_dp,s%mean+5*s%sd,log(max(s%sd*20.0_dp,1e-6_dp))]
    call nelder_mead(obj,x0,fit%optimization,lo,hi,max_iter=mx,tolerance=1e-8_dp)
    fit%model=ecld_new(lambda=fit%optimization%parameters(1), &
      beta=fit%optimization%parameters(2),mu=fit%optimization%parameters(3), &
      sigma=exp(fit%optimization%parameters(4)))
  contains
    function obj(p) result(v)
      real(dp), intent(in) :: p(:)
      real(dp) :: v,c,y,dens
      type(ecld_model) :: d
      integer :: i
      d=ecld_new(lambda=p(1),beta=p(2),mu=p(3),sigma=exp(p(4)))
      c=ecld_const(d)
      if(.not.(c>0.0_dp) .or. c>=huge(1.0_dp)) then; v=huge(1.0_dp)/100.0_dp; return; end if
      v=0.0_dp
      do i=1,size(data)
        y=ecld_solve(d,data(i))
        if(y<log(tiny(1.0_dp))) then; dens=0.0_dp; else; dens=exp(y)/c; end if
        if(dens<=tiny(1.0_dp)) then; v=v+700.0_dp; else; v=v-log(dens); end if
      end do
      v=v/real(size(data),dp)
    end function obj
  end function fit_ecld_moments

  function fit_sld_mle(data,initial,fit_mu,fit_convo,fit_beta,max_iter) result(fit)
    real(dp), intent(in) :: data(:)
    type(sld_model), intent(in) :: initial
    logical, intent(in), optional :: fit_mu,fit_convo,fit_beta
    integer, intent(in), optional :: max_iter
    type(sld_fit_result) :: fit
    logical :: fm,fc,fb
    integer :: npar,mx,j
    real(dp), allocatable :: x0(:),lo(:),hi(:)
    fm=.true.; fc=.true.; fb=.true.; mx=300
    if(present(fit_mu))fm=fit_mu
    if(present(fit_convo))fc=fit_convo
    if(present(fit_beta))fb=fit_beta
    if(present(max_iter))mx=max_iter
    npar=2+merge(1,0,fm)+merge(1,0,fc)+merge(1,0,fb)
    allocate(x0(npar),lo(npar),hi(npar)); j=0
    j=j+1; x0(j)=log(max(initial%nu0,1e-6_dp)); lo(j)=log(1e-6_dp); hi(j)=log(100.0_dp)
    j=j+1; x0(j)=log(max(initial%theta,1e-6_dp)); lo(j)=log(1e-6_dp); hi(j)=log(100.0_dp)
    if(fc)then;j=j+1;x0(j)=log(max(initial%convo,0.01_dp));lo(j)=log(0.01_dp);hi(j)=log(20.0_dp);end if
    if(fb)then;j=j+1;x0(j)=initial%beta_a;lo(j)=-9.99_dp;hi(j)=9.99_dp;end if
    if(fm)then;j=j+1;x0(j)=initial%mu;lo(j)=minval(data)-5*sqrt(variance(data));hi(j)=maxval(data)+5*sqrt(variance(data));end if
    call nelder_mead(obj,x0,fit%optimization,lo,hi,max_iter=mx,tolerance=1e-7_dp)
    fit%model=unpack_parameters(fit%optimization%parameters)
  contains
    function unpack_parameters(p) result(d)
      real(dp), intent(in) :: p(:)
      type(sld_model) :: d
      integer :: k
      d=initial; k=0
      k=k+1; d%nu0=exp(p(k)); k=k+1; d%theta=exp(p(k))
      if(fc)then;k=k+1;d%convo=exp(p(k));end if
      if(fb)then;k=k+1;d%beta_a=p(k);end if
      if(fm)then;k=k+1;d%mu=p(k);end if
    end function unpack_parameters
    function obj(p) result(v)
      real(dp), intent(in) :: p(:)
      real(dp) :: v,pdf
      type(sld_model) :: d
      integer :: i
      d=unpack_parameters(p); v=0.0_dp
      do i=1,size(data)
        pdf=sld_pdf(d,data(i))
        if(pdf<=tiny(1.0_dp))then;v=v+700.0_dp;else;v=v-log(pdf);end if
      end do
      v=v/real(size(data),dp)
    end function obj
    pure function variance(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v,m
      m=sum(x)/real(size(x),dp); v=sum((x-m)**2)/real(max(1,size(x)-1),dp)
    end function variance
  end function fit_sld_mle

end module ecd_fitting
