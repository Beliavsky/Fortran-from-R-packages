! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_models
  use yc_kinds, only : dp
  use yc_types, only : curve_t
  use yc_utils, only : valid_positive_vector, valid_finite_vector, sort_pairs, lower_string
  use yc_linalg, only : weighted_least_squares
  use yc_splines, only : fit_cubic_spline_coefficients
  implicit none
  private
  public :: yc_curve, yc_nelson_siegel, yc_svensson, yc_cubic_spline, yc_fit
  public :: ns_rate_scalar, sv_rate_scalar, ns_forward_scalar, sv_forward_scalar
  public :: ns_loadings_matrix, sv_loadings_matrix

contains

  pure real(dp) function stable_slope_loading(x) result(v)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-7_dp) then
      v = 1.0_dp - x/2.0_dp + x*x/6.0_dp - x*x*x/24.0_dp
    else
      v = (1.0_dp - exp(-x)) / x
    end if
  end function stable_slope_loading

  pure real(dp) function ns_rate_scalar(m, beta0, beta1, beta2, tau) result(r)
    real(dp), intent(in) :: m, beta0, beta1, beta2, tau
    real(dp) :: x, slope, curvature
    x = m / tau
    slope = stable_slope_loading(x)
    curvature = slope - exp(-x)
    r = beta0 + beta1*slope + beta2*curvature
  end function ns_rate_scalar

  pure real(dp) function sv_rate_scalar(m, beta0, beta1, beta2, beta3, tau1, tau2) result(r)
    real(dp), intent(in) :: m, beta0, beta1, beta2, beta3, tau1, tau2
    real(dp) :: x1, x2, s1, c1, c2
    x1 = m/tau1
    x2 = m/tau2
    s1 = stable_slope_loading(x1)
    c1 = s1 - exp(-x1)
    c2 = stable_slope_loading(x2) - exp(-x2)
    r = beta0 + beta1*s1 + beta2*c1 + beta3*c2
  end function sv_rate_scalar

  pure real(dp) function ns_forward_scalar(m, beta0, beta1, beta2, tau) result(f)
    real(dp), intent(in) :: m, beta0, beta1, beta2, tau
    real(dp) :: x
    x = m/tau
    f = beta0 + beta1*exp(-x) + beta2*x*exp(-x)
  end function ns_forward_scalar

  pure real(dp) function sv_forward_scalar(m, beta0, beta1, beta2, beta3, tau1, tau2) result(f)
    real(dp), intent(in) :: m, beta0, beta1, beta2, beta3, tau1, tau2
    real(dp) :: x1, x2
    x1 = m/tau1
    x2 = m/tau2
    f = beta0 + beta1*exp(-x1) + beta2*x1*exp(-x1) + beta3*x2*exp(-x2)
  end function sv_forward_scalar

  function ns_loadings_matrix(m, tau) result(xmat)
    real(dp), intent(in) :: m(:), tau
    real(dp), allocatable :: xmat(:,:)
    real(dp) :: x, s
    integer :: i
    allocate(xmat(size(m),3))
    do i = 1, size(m)
      x = m(i)/tau
      s = stable_slope_loading(x)
      xmat(i,:) = [1.0_dp, s, s-exp(-x)]
    end do
  end function ns_loadings_matrix

  function sv_loadings_matrix(m, tau1, tau2) result(xmat)
    real(dp), intent(in) :: m(:), tau1, tau2
    real(dp), allocatable :: xmat(:,:)
    real(dp) :: x1, x2, s1
    integer :: i
    allocate(xmat(size(m),4))
    do i = 1, size(m)
      x1 = m(i)/tau1
      x2 = m(i)/tau2
      s1 = stable_slope_loading(x1)
      xmat(i,:) = [1.0_dp, s1, s1-exp(-x1), stable_slope_loading(x2)-exp(-x2)]
    end do
  end function sv_loadings_matrix

  function yc_curve(maturities, rates, rate_type) result(curve)
    real(dp), intent(in) :: maturities(:), rates(:)
    character(len=*), intent(in), optional :: rate_type
    type(curve_t) :: curve
    real(dp), allocatable :: m(:), r(:)
    if (size(maturities) /= size(rates) .or. .not. valid_positive_vector(maturities) .or. &
        .not. valid_finite_vector(rates)) then
      curve%ok = .false.
      curve%message = 'Invalid maturities or rates.'
      return
    end if
    allocate(m, source=maturities)
    allocate(r, source=rates)
    call sort_pairs(m, r)
    curve%maturities = m
    curve%rates = r
    curve%method = 'observed'
    if (present(rate_type)) curve%rate_type = trim(lower_string(rate_type))
  end function yc_curve

  subroutine profile_ns(m, r, w, tau, beta, sse, ok)
    real(dp), intent(in) :: m(:), r(:), w(:), tau
    real(dp), allocatable, intent(out) :: beta(:)
    real(dp), intent(out) :: sse
    logical, intent(out) :: ok
    real(dp), allocatable :: x(:,:)
    x = ns_loadings_matrix(m, tau)
    call weighted_least_squares(x, r, w, beta, sse, ok)
  end subroutine profile_ns

  subroutine profile_sv(m, r, w, tau1, tau2, beta, sse, ok)
    real(dp), intent(in) :: m(:), r(:), w(:), tau1, tau2
    real(dp), allocatable, intent(out) :: beta(:)
    real(dp), intent(out) :: sse
    logical, intent(out) :: ok
    real(dp), allocatable :: x(:,:)
    if (abs(tau1-tau2) < 1.0e-5_dp) then
      allocate(beta(4))
      beta = 0.0_dp
      sse = huge(1.0_dp)
      ok = .false.
      return
    end if
    x = sv_loadings_matrix(m, tau1, tau2)
    call weighted_least_squares(x, r, w, beta, sse, ok)
  end subroutine profile_sv

  function yc_nelson_siegel(maturities, rates, tau_init, weights, rate_type) result(curve)
    real(dp), intent(in) :: maturities(:), rates(:)
    real(dp), intent(in), optional :: tau_init
    real(dp), intent(in), optional :: weights(:)
    character(len=*), intent(in), optional :: rate_type
    type(curve_t) :: curve
    real(dp), allocatable :: m(:), r(:), w(:), beta(:), beta_try(:)
    real(dp) :: a, b, c, d, fc, fd, sse, sse_try, tol, tau_best
    real(dp), parameter :: tau_grid(7) = [0.5_dp,1.0_dp,2.0_dp,3.0_dp,5.0_dp,8.0_dp,12.0_dp]
    integer :: i, iter
    logical :: ok

    curve = yc_curve(maturities, rates, rate_type)
    if (.not. curve%ok) return
    m = curve%maturities
    r = curve%rates
    allocate(w(size(m)))
    w = 1.0_dp
    if (present(weights)) then
      if (size(weights) /= size(m) .or. any(weights < 0.0_dp)) then
        curve%ok = .false.; curve%message = 'Invalid weights.'; return
      end if
      w = weights
      call sort_pairs_copy_weights(maturities, rates, weights, m, r, w)
    end if
    tau_best = 1.0_dp
    if (present(tau_init)) tau_best = tau_init
    sse = huge(1.0_dp)
    do i = 1, size(tau_grid)
      call profile_ns(m,r,w,tau_grid(i),beta_try,sse_try,ok)
      if (ok .and. sse_try < sse) then
        sse = sse_try
        tau_best = tau_grid(i)
        beta = beta_try
      end if
    end do
    a = 0.01_dp
    b = 30.0_dp
    tol = 1.0e-10_dp
    c = b - 0.6180339887498948482_dp*(b-a)
    d = a + 0.6180339887498948482_dp*(b-a)
    fc = ns_sse_at(c,m,r,w)
    fd = ns_sse_at(d,m,r,w)
    do iter = 1, 250
      if (abs(b-a) <= tol*(1.0_dp+abs(a)+abs(b))) exit
      if (fc < fd) then
        b=d; d=c; fd=fc; c=b-0.6180339887498948482_dp*(b-a); fc=ns_sse_at(c,m,r,w)
      else
        a=c; c=d; fc=fd; d=a+0.6180339887498948482_dp*(b-a); fd=ns_sse_at(d,m,r,w)
      end if
    end do
    tau_best = 0.5_dp*(a+b)
    call profile_ns(m,r,w,tau_best,beta,sse,ok)
    if (.not. ok) then
      curve%ok=.false.; curve%message='Nelson-Siegel fit failed.'; return
    end if
    curve%method='nelson_siegel'
    curve%beta0=beta(1); curve%beta1=beta(2); curve%beta2=beta(3); curve%tau=tau_best
    curve%objective=sse; curve%iterations=iter
    allocate(curve%fitted(size(m)),curve%residuals(size(m)))
    do i=1,size(m)
      curve%fitted(i)=ns_rate_scalar(m(i),curve%beta0,curve%beta1,curve%beta2,curve%tau)
    end do
    curve%residuals=r-curve%fitted
  contains
    real(dp) function ns_sse_at(t,mv,rv,wv) result(v)
      real(dp), intent(in)::t,mv(:),rv(:),wv(:)
      real(dp), allocatable :: bt(:)
      logical :: lok
      call profile_ns(mv,rv,wv,t,bt,v,lok)
      if(.not.lok) v=huge(1.0_dp)
    end function ns_sse_at
  end function yc_nelson_siegel

  function yc_svensson(maturities, rates, tau1_init, tau2_init, weights, rate_type) result(curve)
    real(dp), intent(in) :: maturities(:), rates(:)
    real(dp), intent(in), optional :: tau1_init, tau2_init
    real(dp), intent(in), optional :: weights(:)
    character(len=*), intent(in), optional :: rate_type
    type(curve_t) :: curve
    real(dp), allocatable :: m(:),r(:),w(:),beta(:),vertices(:,:),fval(:),centroid(:),trial(:)
    real(dp) :: t1_grid(4),t2_grid(4),sse,sse_try,t1best,t2best,step1,step2
    integer :: i,j,iter,ilo,ihi,isecond
    logical :: ok

    curve=yc_curve(maturities,rates,rate_type)
    if(.not.curve%ok)return
    m=curve%maturities; r=curve%rates
    allocate(w(size(m))); w=1.0_dp
    if(present(weights))then
      if(size(weights)/=size(m).or.any(weights<0.0_dp))then
        curve%ok=.false.; curve%message='Invalid weights.'; return
      end if
      call sort_pairs_copy_weights(maturities,rates,weights,m,r,w)
    end if
    t1_grid=[1.0_dp,2.0_dp,5.0_dp,8.0_dp]
    t2_grid=[0.5_dp,1.0_dp,3.0_dp,5.0_dp]
    t1best=1.0_dp
    t2best=5.0_dp
    if(present(tau1_init))t1best=tau1_init
    if(present(tau2_init))t2best=tau2_init
    sse=huge(1.0_dp)
    do i=1,4
      do j=1,4
        if(abs(t1_grid(i)-t2_grid(j))<0.1_dp)cycle
        call profile_sv(m,r,w,t1_grid(i),t2_grid(j),beta,sse_try,ok)
        if(ok.and.sse_try<sse)then
          sse=sse_try; t1best=t1_grid(i); t2best=t2_grid(j)
        end if
      end do
    end do
    allocate(vertices(2,3),fval(3),centroid(2),trial(2))
    step1=max(0.25_dp,0.15_dp*t1best); step2=max(0.25_dp,0.15_dp*t2best)
    vertices(:,1)=[t1best,t2best]
    vertices(:,2)=[min(30.0_dp,t1best+step1),t2best]
    vertices(:,3)=[t1best,min(30.0_dp,t2best+step2)]
    do i=1,3; fval(i)=sv_sse_at(vertices(1,i),vertices(2,i),m,r,w); end do
    do iter=1,600
      call order_three(fval,ilo,isecond,ihi)
      if(maxval(abs(vertices(:,2)-vertices(:,1)))<1.0e-8_dp.and. &
         maxval(abs(vertices(:,3)-vertices(:,1)))<1.0e-8_dp)exit
      centroid=0.5_dp*(vertices(:,ilo)+vertices(:,isecond))
      trial=clamp_tau(centroid+(centroid-vertices(:,ihi)))
      sse_try=sv_sse_at(trial(1),trial(2),m,r,w)
      if(sse_try<fval(ilo))then
        trial=clamp_tau(centroid+2.0_dp*(centroid-vertices(:,ihi)))
        if(sv_sse_at(trial(1),trial(2),m,r,w)<sse_try)then
          vertices(:,ihi)=trial; fval(ihi)=sv_sse_at(trial(1),trial(2),m,r,w)
        else
          vertices(:,ihi)=clamp_tau(centroid+(centroid-vertices(:,ihi))); fval(ihi)=sse_try
        end if
      else if(sse_try<fval(isecond))then
        vertices(:,ihi)=trial; fval(ihi)=sse_try
      else
        trial=clamp_tau(centroid+0.5_dp*(vertices(:,ihi)-centroid))
        sse_try=sv_sse_at(trial(1),trial(2),m,r,w)
        if(sse_try<fval(ihi))then
          vertices(:,ihi)=trial; fval(ihi)=sse_try
        else
          vertices(:,isecond)=vertices(:,ilo)+0.5_dp*(vertices(:,isecond)-vertices(:,ilo))
          vertices(:,ihi)=vertices(:,ilo)+0.5_dp*(vertices(:,ihi)-vertices(:,ilo))
          fval(isecond)=sv_sse_at(vertices(1,isecond),vertices(2,isecond),m,r,w)
          fval(ihi)=sv_sse_at(vertices(1,ihi),vertices(2,ihi),m,r,w)
        end if
      end if
    end do
    ilo=minloc(fval,dim=1)
    t1best=vertices(1,ilo); t2best=vertices(2,ilo)
    call profile_sv(m,r,w,t1best,t2best,beta,sse,ok)
    if(.not.ok)then; curve%ok=.false.; curve%message='Svensson fit failed.'; return; end if
    curve%method='svensson'; curve%beta0=beta(1); curve%beta1=beta(2); curve%beta2=beta(3); curve%beta3=beta(4)
    curve%tau1=t1best; curve%tau2=t2best; curve%objective=sse; curve%iterations=iter
    allocate(curve%fitted(size(m)),curve%residuals(size(m)))
    do i=1,size(m)
      curve%fitted(i)=sv_rate_scalar(m(i),curve%beta0,curve%beta1,curve%beta2,curve%beta3,curve%tau1,curve%tau2)
    end do
    curve%residuals=r-curve%fitted
  contains
    real(dp) function sv_sse_at(a,b,mv,rv,wv) result(v)
      real(dp),intent(in)::a,b,mv(:),rv(:),wv(:)
      real(dp),allocatable::bt(:)
      logical::lok
      call profile_sv(mv,rv,wv,a,b,bt,v,lok)
      if(.not.lok)v=huge(1.0_dp)
    end function sv_sse_at
    pure function clamp_tau(v) result(out)
      real(dp),intent(in)::v(2); real(dp)::out(2)
      out=max(0.01_dp,min(30.0_dp,v))
      if(abs(out(1)-out(2))<1.0e-5_dp)out(2)=min(30.0_dp,out(2)+1.0e-3_dp)
    end function clamp_tau
    subroutine order_three(f,lo,mid,hi)
      real(dp),intent(in)::f(3); integer,intent(out)::lo,mid,hi
      integer::idx(3),a,b,tmp
      idx=[1,2,3]
      do a=1,2; do b=a+1,3; if(f(idx(b))<f(idx(a)))then; tmp=idx(a);idx(a)=idx(b);idx(b)=tmp;end if;end do;end do
      lo=idx(1);mid=idx(2);hi=idx(3)
    end subroutine order_three
  end function yc_svensson

  function yc_cubic_spline(maturities,rates,spline_method,rate_type) result(curve)
    real(dp),intent(in)::maturities(:),rates(:)
    character(len=*),intent(in),optional::spline_method,rate_type
    type(curve_t)::curve
    character(len=16)::meth
    logical::ok
    curve=yc_curve(maturities,rates,rate_type)
    if(.not.curve%ok)return
    if(size(curve%maturities)<3)then;curve%ok=.false.;curve%message='Cubic spline requires at least 3 points.';return;end if
    meth='natural';if(present(spline_method))meth=trim(lower_string(spline_method))
    if(meth/='natural'.and.meth/='fmm')then;curve%ok=.false.;curve%message='Unknown spline method.';return;end if
    call fit_cubic_spline_coefficients(curve%maturities,curve%rates,meth,curve%spline_b,curve%spline_c,curve%spline_d,ok)
    if(.not.ok)then;curve%ok=.false.;curve%message='Spline fit failed.';return;end if
    curve%method='cubic_spline';curve%spline_method=meth
    allocate(curve%fitted(size(curve%rates)),curve%residuals(size(curve%rates)))
    curve%fitted=curve%rates;curve%residuals=0.0_dp
  end function yc_cubic_spline

  function yc_fit(maturities,rates,method,rate_type) result(curve)
    real(dp),intent(in)::maturities(:),rates(:)
    character(len=*),intent(in)::method
    character(len=*),intent(in),optional::rate_type
    type(curve_t)::curve
    character(len=:),allocatable::meth
    meth=trim(lower_string(method))
    select case(meth)
    case('nelson_siegel'); curve=yc_nelson_siegel(maturities,rates,rate_type=rate_type)
    case('svensson'); curve=yc_svensson(maturities,rates,rate_type=rate_type)
    case('cubic_spline'); curve=yc_cubic_spline(maturities,rates,rate_type=rate_type)
    case default; curve%ok=.false.;curve%message='Unknown fitting method.'
    end select
  end function yc_fit

  subroutine sort_pairs_copy_weights(mi,ri,wi,mo,ro,wo)
    real(dp),intent(in)::mi(:),ri(:),wi(:)
    real(dp),intent(out)::mo(:),ro(:),wo(:)
    mo=mi;ro=ri;wo=wi;call sort_pairs(mo,ro,wo)
  end subroutine sort_pairs_copy_weights

end module yc_models
