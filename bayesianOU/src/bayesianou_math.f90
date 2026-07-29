! SPDX-License-Identifier: MIT
module bayesianou_math
  use bayesianou_kinds, only : dp, pi, log_two_pi, status_ok, status_bad_input, &
                               status_not_converged, status_singular
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  type, public :: rng_state
    integer(int64) :: state = 88172645463325252_int64
  end type rng_state

  abstract interface
    function scalar_objective(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function scalar_objective
  end interface

  public :: rng_seed, rng_uniform, rng_normal, rng_gamma, rng_student_t
  public :: normal_logpdf, student_t_logpdf, normal_cdf, normal_quantile
  public :: log_sum_exp, logistic, logit, sample_mean, sample_sd, covariance_matrix
  public :: percentile, median_value, sort_real, ols_fit, solve_linear, invert_spd
  public :: symmetric_eigen, first_principal_component, nelder_mead
  public :: numerical_gradient, numerical_hessian, finite_all, clamp
  public :: autocorrelation, effective_sample_size, split_rhat, t_quantile_approx

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import dp
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*), b(ldb,*)
    end subroutine dgesv
    subroutine dpotrf(uplo, n, a, lda, info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
    end subroutine dpotrf
    subroutine dpotri(uplo, n, a, lda, info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
    end subroutine dpotri
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: w(*), work(*)
    end subroutine dsyev
  end interface

contains

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: seed
    integer(int64) :: z
    z = int(seed, int64)
    if (z == 0_int64) z = 104729_int64
    rng%state = ieor(z, int(z'9E3779B97F4A7C15', int64))
    call rng_warm(rng)
  end subroutine rng_seed

  subroutine rng_warm(rng)
    type(rng_state), intent(inout) :: rng
    integer :: i
    real(dp) :: u
    do i = 1, 16
      u = rng_uniform(rng)
    end do
  end subroutine rng_warm

  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(int64) :: x
    x = rng%state
    x = ieor(x, ishft(x, 13))
    x = ieor(x, ishft(x, -7))
    x = ieor(x, ishft(x, 17))
    rng%state = x
    u = real(iand(x, int(z'001FFFFFFFFFFFFF', int64)), dp) / real(int(z'0020000000000000', int64), dp)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
  end function rng_uniform

  function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z, u1, u2
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function rng_normal

  recursive function rng_gamma(rng, shape) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape
    real(dp) :: x, d, c, z, u, v
    if (shape <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = rng_gamma(rng, shape + 1.0_dp) * rng_uniform(rng)**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      z = rng_normal(rng)
      v = 1.0_dp + c*z
      if (v <= 0.0_dp) cycle
      v = v**3
      u = rng_uniform(rng)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
    end do
    x = d*v
  end function rng_gamma

  function rng_student_t(rng, nu) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: nu
    real(dp) :: x
    x = rng_normal(rng)/sqrt(2.0_dp*rng_gamma(rng, 0.5_dp*nu)/nu)
  end function rng_student_t

  pure elemental function normal_logpdf(x, mu, sigma) result(v)
    real(dp), intent(in) :: x, mu, sigma
    real(dp) :: v, z
    if (sigma <= 0.0_dp) then
      v = -huge(1.0_dp)
    else
      z = (x-mu)/sigma
      v = -0.5_dp*log_two_pi - log(sigma) - 0.5_dp*z*z
    end if
  end function normal_logpdf

  pure elemental function student_t_logpdf(x, nu, mu, sigma) result(v)
    real(dp), intent(in) :: x, nu, mu, sigma
    real(dp) :: v, z
    if (sigma <= 0.0_dp .or. nu <= 0.0_dp) then
      v = -huge(1.0_dp)
    else
      z = (x-mu)/sigma
      v = log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) &
          -0.5_dp*log(nu*pi) - log(sigma) - 0.5_dp*(nu+1.0_dp)*log(1.0_dp+z*z/nu)
    end if
  end function student_t_logpdf

  pure elemental function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure elemental function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ 7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp ]
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  pure function log_sum_exp(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v, m
    if (size(x) == 0) then
      v = -huge(1.0_dp)
      return
    end if
    m = maxval(x)
    if (.not. ieee_is_finite_local(m)) then
      v = m
    else
      v = m + log(sum(exp(x-m)))
    end if
  end function log_sum_exp

  pure elemental function logistic(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if (x >= 0.0_dp) then
      y = 1.0_dp/(1.0_dp+exp(-x))
    else
      y = exp(x)/(1.0_dp+exp(x))
    end if
  end function logistic

  pure elemental function logit(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x
    x = log(p/(1.0_dp-p))
  end function logit

  pure elemental function clamp(x, lo, hi) result(y)
    real(dp), intent(in) :: x, lo, hi
    real(dp) :: y
    y = min(max(x,lo),hi)
  end function clamp

  pure function sample_mean(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) == 0) then
      m = 0.0_dp
    else
      m = sum(x)/real(size(x),dp)
    end if
  end function sample_mean

  pure function sample_sd(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: s, m
    if (size(x) < 2) then
      s = 0.0_dp
    else
      m = sample_mean(x)
      s = sqrt(sum((x-m)**2)/real(size(x)-1,dp))
    end if
  end function sample_sd

  subroutine covariance_matrix(x, cov)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: cov(size(x,2),size(x,2))
    real(dp) :: mu(size(x,2))
    integer :: i
    mu = sum(x,dim=1)/real(size(x,1),dp)
    cov = 0.0_dp
    do i=1,size(x,1)
      cov = cov + spread(x(i,:)-mu,2,size(x,2))*spread(x(i,:)-mu,1,size(x,2))
    end do
    cov = cov/real(max(1,size(x,1)-1),dp)
  end subroutine covariance_matrix

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if (x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  function percentile(x,p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp) :: q, h, frac
    real(dp), allocatable :: y(:)
    integer :: k, n
    n=size(x)
    if (n==0) then
      q=0.0_dp; return
    end if
    allocate(y(n)); y=x; call sort_real(y)
    h=1.0_dp+(real(n-1,dp))*clamp(p,0.0_dp,1.0_dp)
    k=int(floor(h)); frac=h-real(k,dp)
    if (k>=n) then
      q=y(n)
    else
      q=(1.0_dp-frac)*y(k)+frac*y(k+1)
    end if
  end function percentile

  function median_value(x) result(q)
    real(dp), intent(in) :: x(:)
    real(dp) :: q
    q=percentile(x,0.5_dp)
  end function median_value

  subroutine solve_linear(a,b,x,status)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(size(b))
    integer, intent(out) :: status
    real(dp), allocatable :: aa(:,:), bb(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n, info
    n=size(b); allocate(aa(n,n),bb(n,1),ipiv(n)); aa=a; bb(:,1)=b
    call dgesv(n,1,aa,n,ipiv,bb,n,info)
    if(info==0) then; x=bb(:,1); status=status_ok
    else; x=0.0_dp; status=status_singular; end if
  end subroutine solve_linear

  subroutine ols_fit(x,y,beta,resid,cov,status,ridge)
    real(dp), intent(in) :: x(:,:), y(:)
    real(dp), intent(out) :: beta(size(x,2)), resid(size(y)), cov(size(x,2),size(x,2))
    integer, intent(out) :: status
    real(dp), intent(in), optional :: ridge
    real(dp) :: xtx(size(x,2),size(x,2)), xty(size(x,2)), rr, sig2
    integer :: j
    xtx=matmul(transpose(x),x); xty=matmul(transpose(x),y)
    rr=1.0e-10_dp; if(present(ridge)) rr=ridge
    do j=1,size(x,2); xtx(j,j)=xtx(j,j)+rr; end do
    call solve_linear(xtx,xty,beta,status)
    resid=y-matmul(x,beta)
    sig2=sum(resid**2)/real(max(1,size(y)-size(beta)),dp)
    call invert_spd(xtx,cov,status)
    cov=sig2*cov
  end subroutine ols_fit

  subroutine invert_spd(a,ainv,status,logdet)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    integer, intent(out) :: status
    real(dp), intent(out), optional :: logdet
    integer :: n,info,i,j
    n=size(a,1); ainv=0.5_dp*(a+transpose(a))
    call dpotrf('L',n,ainv,n,info)
    if(info/=0) then; status=status_singular; ainv=0.0_dp; if(present(logdet)) logdet=huge(1.0_dp); return; end if
    if(present(logdet)) then
      logdet=0.0_dp; do i=1,n; logdet=logdet+2.0_dp*log(ainv(i,i)); end do
    end if
    call dpotri('L',n,ainv,n,info)
    do j=1,n; do i=1,j-1; ainv(i,j)=ainv(j,i); end do; end do
    status=merge(status_ok,status_singular,info==0)
  end subroutine invert_spd

  subroutine symmetric_eigen(a,values,vectors,status)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: values(size(a,1)), vectors(size(a,1),size(a,2))
    integer, intent(out) :: status
    real(dp), allocatable :: work(:)
    integer :: n,lwork,info
    n=size(a,1); vectors=0.5_dp*(a+transpose(a)); lwork=max(1,3*n-1); allocate(work(lwork))
    call dsyev('V','U',n,vectors,n,values,work,lwork,info)
    status=merge(status_ok,status_singular,info==0)
  end subroutine symmetric_eigen

  subroutine first_principal_component(x, score, loading, status)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: score(size(x,1)), loading(size(x,2))
    integer, intent(out) :: status
    real(dp) :: cov(size(x,2),size(x,2)), vals(size(x,2)), vecs(size(x,2),size(x,2))
    call covariance_matrix(x,cov); call symmetric_eigen(cov,vals,vecs,status)
    if(status/=status_ok) then; score=0.0_dp; loading=0.0_dp; return; end if
    loading=vecs(:,size(vals)); score=matmul(x,loading)
  end subroutine first_principal_component

  subroutine numerical_gradient(fun,x,g,step)
    procedure(scalar_objective) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(size(x))
    real(dp), intent(in), optional :: step
    real(dp) :: h, xp(size(x)), xm(size(x)), fp, fm
    integer :: i
    h=1.0e-5_dp; if(present(step)) h=step
    do i=1,size(x)
      xp=x; xm=x; xp(i)=xp(i)+h*(1.0_dp+abs(x(i))); xm(i)=xm(i)-h*(1.0_dp+abs(x(i)))
      fp=fun(xp); fm=fun(xm)
      if(finite_all([fp,fm])) then; g(i)=(fp-fm)/(xp(i)-xm(i)); else; g(i)=0.0_dp; end if
    end do
  end subroutine numerical_gradient

  subroutine numerical_hessian(fun,x,hess,step)
    procedure(scalar_objective) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hess(size(x),size(x))
    real(dp), intent(in), optional :: step
    real(dp) :: h, f0, fpp, fpm, fmp, fmm
    real(dp) :: xp(size(x)), xm(size(x)), xpp(size(x)), xpm(size(x)), xmp(size(x)), xmm(size(x))
    integer :: i,j
    h=1.0e-4_dp; if(present(step)) h=step
    f0=fun(x); hess=0.0_dp
    do i=1,size(x)
      xp=x; xm=x; xp(i)=xp(i)+h*(1+abs(x(i))); xm(i)=xm(i)-h*(1+abs(x(i)))
      fpp=fun(xp); fmm=fun(xm)
      if(finite_all([fpp,fmm,f0])) hess(i,i)=(fpp-2*f0+fmm)/(xp(i)-x(i))**2
      do j=i+1,size(x)
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+h*(1+abs(x(i))); xpp(j)=xpp(j)+h*(1+abs(x(j)))
        xpm(i)=xpm(i)+h*(1+abs(x(i))); xpm(j)=xpm(j)-h*(1+abs(x(j)))
        xmp(i)=xmp(i)-h*(1+abs(x(i))); xmp(j)=xmp(j)+h*(1+abs(x(j)))
        xmm(i)=xmm(i)-h*(1+abs(x(i))); xmm(j)=xmm(j)-h*(1+abs(x(j)))
        fpp=fun(xpp); fpm=fun(xpm); fmp=fun(xmp); fmm=fun(xmm)
        if(finite_all([fpp,fpm,fmp,fmm])) then
          hess(i,j)=(fpp-fpm-fmp+fmm)/((xpp(i)-xmp(i))*(xpp(j)-xpm(j)))
          hess(j,i)=hess(i,j)
        end if
      end do
    end do
  end subroutine numerical_hessian

  subroutine nelder_mead(fun,x,f,status,max_iter,tol,step)
    procedure(scalar_objective) :: fun
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: f
    integer, intent(out) :: status
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol, step
    integer :: n,it,j,limit
    real(dp) :: tolerance,delta,fr,fe,fc,centroid(size(x)),xr(size(x)),xe(size(x)),xc(size(x))
    real(dp), allocatable :: simp(:,:), fv(:)
    n=size(x); limit=2000; if(present(max_iter)) limit=max_iter
    tolerance=1e-8_dp; if(present(tol)) tolerance=tol
    delta=0.1_dp; if(present(step)) delta=step
    allocate(simp(n,n+1),fv(n+1)); simp(:,1)=x
    do j=1,n; simp(:,j+1)=x; simp(j,j+1)=simp(j,j+1)+delta*(1+abs(x(j))); end do
    do j=1,n+1; fv(j)=fun(simp(:,j)); end do
    status=status_not_converged
    do it=1,limit
      call order_simplex(simp,fv)
      if(maxval(abs(fv-fv(1)))<tolerance .and. maxval(abs(simp-spread(simp(:,1),2,n+1)))<sqrt(tolerance)) then
        status=status_ok; exit
      end if
      centroid=sum(simp(:,1:n),dim=2)/real(n,dp)
      xr=centroid+(centroid-simp(:,n+1)); fr=fun(xr)
      if(fr<fv(1)) then
        xe=centroid+2*(xr-centroid); fe=fun(xe)
        if(fe<fr) then; simp(:,n+1)=xe; fv(n+1)=fe; else; simp(:,n+1)=xr; fv(n+1)=fr; end if
      else if(fr<fv(n)) then
        simp(:,n+1)=xr; fv(n+1)=fr
      else
        if(fr<fv(n+1)) then; xc=centroid+0.5_dp*(xr-centroid); else; xc=centroid+0.5_dp*(simp(:,n+1)-centroid); end if
        fc=fun(xc)
        if(fc<min(fr,fv(n+1))) then
          simp(:,n+1)=xc; fv(n+1)=fc
        else
          do j=2,n+1
            simp(:,j)=simp(:,1)+0.5_dp*(simp(:,j)-simp(:,1)); fv(j)=fun(simp(:,j))
          end do
        end if
      end if
    end do
    call order_simplex(simp,fv); x=simp(:,1); f=fv(1)
  contains
    subroutine order_simplex(s,fs)
      real(dp), intent(inout) :: s(:,:),fs(:)
      integer :: a,b,k
      real(dp) :: tf,col(size(s,1))
      do a=1,size(fs)-1
        k=a
        do b=a+1,size(fs); if(fs(b)<fs(k)) k=b; end do
        if(k/=a) then; tf=fs(a);fs(a)=fs(k);fs(k)=tf;col=s(:,a);s(:,a)=s(:,k);s(:,k)=col;end if
      end do
    end subroutine order_simplex
  end subroutine nelder_mead

  pure function finite_all(x) result(ok)
    real(dp), intent(in) :: x(:)
    logical :: ok
    integer :: i
    ok=.true.
    do i=1,size(x); if(.not.ieee_is_finite_local(x(i))) then;ok=.false.;return;end if;end do
  end function finite_all

  pure function ieee_is_finite_local(x) result(ok)
    real(dp), intent(in) :: x
    logical :: ok
    ok=ieee_is_finite(x)
  end function ieee_is_finite_local

  function autocorrelation(x,lag) result(r)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    real(dp) :: r,m,den
    if(lag<0 .or. lag>=size(x)) then;r=0;return;end if
    m=sample_mean(x); den=sum((x-m)**2)
    if(den<=0) then;r=0;else;r=sum((x(1:size(x)-lag)-m)*(x(1+lag:size(x))-m))/den;end if
  end function autocorrelation

  function effective_sample_size(x) result(ess)
    real(dp), intent(in) :: x(:)
    real(dp) :: ess,s,r1,r2
    integer :: k,n
    n=size(x); s=0.0_dp
    do k=1,(n-1)/2
      r1=autocorrelation(x,2*k-1); r2=autocorrelation(x,2*k)
      if(r1+r2<0) exit
      s=s+r1+r2
    end do
    ess=real(n,dp)/max(1.0_dp,1.0_dp+2.0_dp*s)
  end function effective_sample_size

  function split_rhat(chains) result(rhat)
    real(dp), intent(in) :: chains(:,:)
    real(dp) :: rhat
    integer :: n,m,j
    real(dp), allocatable :: halves(:,:), means(:), vars(:)
    real(dp) :: W,B,varhat
    n=size(chains,1)/2; m=2*size(chains,2)
    if(n<2 .or. size(chains,2)<2) then;rhat=1.0_dp;return;end if
    allocate(halves(n,m),means(m),vars(m))
    do j=1,size(chains,2)
      halves(:,2*j-1)=chains(1:n,j); halves(:,2*j)=chains(size(chains,1)-n+1:size(chains,1),j)
    end do
    do j=1,m; means(j)=sample_mean(halves(:,j)); vars(j)=sample_sd(halves(:,j))**2; end do
    W=sample_mean(vars); B=real(n,dp)*sample_sd(means)**2
    if(W<=0) then;rhat=1.0_dp;else;varhat=(real(n-1,dp)/n)*W+B/n;rhat=sqrt(varhat/W);end if
  end function split_rhat

  pure function t_quantile_approx(p,df) result(q)
    real(dp), intent(in) :: p,df
    real(dp) :: q,z
    z=normal_quantile(p)
    if(df>1.0e5_dp) then;q=z;else;q=z+(z**3+z)/(4*df)+(5*z**5+16*z**3+3*z)/(96*df**2);end if
  end function t_quantile_approx

end module bayesianou_math
