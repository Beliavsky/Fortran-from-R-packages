module mev_math
  use mev_kinds, only: dp, pi
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
  implicit none
  private
  public :: normal_cdf, normal_quantile, reg_gamma_p, gamma_quantile, reg_beta, beta_quantile
  public :: digamma_mev, trigamma_mev, rng_normal, rng_gamma, rng_beta
  public :: sort_ascending, sort_descending, mean_real, variance_real, median_real
  public :: chol_lower, solve_linear, inverse_matrix, logdet_spd, mvnormal_sample
  public :: outer_product, covariance_matrix, empirical_quantile
  public :: pattern_minimize, finite_diff_hessian, finite_diff_gradient
  public :: normal_pdf, student_t_cdf, student_t_quantile, student_t_pdf
  public :: expm1_mev, log1p_mev

  abstract interface
    function scalar_objective(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function scalar_objective
  end interface

contains


  pure real(dp) function expm1_mev(x) result(v)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-5_dp) then
      v = x*(1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x*(1.0_dp/24.0_dp + x/120.0_dp))))
    else
      v = exp(x)-1.0_dp
    end if
  end function expm1_mev

  pure real(dp) function log1p_mev(x) result(v)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-5_dp) then
      v = x*(1.0_dp + x*(-0.5_dp + x*(1.0_dp/3.0_dp + x*(-0.25_dp + x/5.0_dp))))
    else
      v = log(1.0_dp+x)
    end if
  end function log1p_mev

  pure real(dp) function normal_pdf(x) result(v)
    real(dp), intent(in) :: x
    v = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: lo, hi, mid
    integer :: it
    if (p <= 0.0_dp) then
      x = -ieee_value(0.0_dp, ieee_positive_inf)
      return
    else if (p >= 1.0_dp) then
      x = ieee_value(0.0_dp, ieee_positive_inf)
      return
    end if
    lo = -12.0_dp
    hi = 12.0_dp
    do it = 1, 120
      mid = 0.5_dp*(lo + hi)
      if (normal_cdf(mid) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo + hi)
  end function normal_quantile

  pure real(dp) function reg_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: itmax = 800
    real(dp), parameter :: eps = 1.0e-14_dp, fpmin = 1.0e-300_dp
    real(dp) :: sumv, del, ap, b, c, d, h, an
    integer :: n
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, itmax
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x + a*log(x) - log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/max(abs(b), fpmin)
      if (b < 0.0_dp) d = -d
      h = d
      do n = 1, itmax
        an = -real(n,dp)*(real(n,dp)-a)
        b = b + 2.0_dp
        d = an*d + b
        if (abs(d) < fpmin) d = sign(fpmin,d)
        c = b + an/c
        if (abs(c) < fpmin) c = sign(fpmin,c)
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) < eps) exit
      end do
      p = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function reg_gamma_p


  pure real(dp) function gamma_quantile(prob,a,scale) result(x)
    real(dp), intent(in) :: prob,a
    real(dp), intent(in), optional :: scale
    real(dp) :: lo,hi,mid,sc
    integer :: it
    sc=1.0_dp; if(present(scale)) sc=scale
    if(prob<=0.0_dp) then
      x=0.0_dp; return
    else if(prob>=1.0_dp) then
      x=ieee_value(0.0_dp,ieee_positive_inf); return
    end if
    if(a<=0.0_dp .or. sc<=0.0_dp) then
      x=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    lo=0.0_dp; hi=max(1.0_dp,a)
    do while(reg_gamma_p(a,hi)<prob)
      hi=2.0_dp*hi
      if(hi>huge(1.0_dp)/4.0_dp) exit
    end do
    do it=1,160
      mid=0.5_dp*(lo+hi)
      if(reg_gamma_p(a,mid)<prob) then
        lo=mid
      else
        hi=mid
      end if
    end do
    x=sc*0.5_dp*(lo+hi)
  end function gamma_quantile

  pure real(dp) function beta_cf(a, b, x) result(h)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 800
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    integer :: m, m2
    real(dp) :: aa, c, d, del, qab, qam, qap
    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, maxit
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x / &
           ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
           ((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del-1.0_dp) < eps) exit
    end do
  end function beta_cf

  pure real(dp) function reg_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      p = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
             a*log(x) + b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      p = bt*beta_cf(a,b,x)/a
    else
      p = 1.0_dp - bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    p = max(0.0_dp, min(1.0_dp,p))
  end function reg_beta

  pure real(dp) function beta_quantile(prob, a, b) result(x)
    real(dp), intent(in) :: prob, a, b
    real(dp) :: lo, hi, mid
    integer :: it
    if (prob <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (prob >= 1.0_dp) then
      x = 1.0_dp
      return
    end if
    lo = 0.0_dp
    hi = 1.0_dp
    do it = 1, 150
      mid = 0.5_dp*(lo+hi)
      if (reg_beta(mid,a,b) < prob) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function beta_quantile

  pure real(dp) function digamma_mev(x0) result(y)
    real(dp), intent(in) :: x0
    real(dp) :: x, r
    x = x0
    y = 0.0_dp
    if (x <= 0.0_dp) then
      y = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    do while (x < 8.0_dp)
      y = y - 1.0_dp/x
      x = x + 1.0_dp
    end do
    r = 1.0_dp/x
    y = y + log(x) - 0.5_dp*r - r*r*(1.0_dp/12.0_dp - &
        r*r*(1.0_dp/120.0_dp-r*r/252.0_dp))
  end function digamma_mev

  pure real(dp) function trigamma_mev(x0) result(y)
    real(dp), intent(in) :: x0
    real(dp) :: x, r
    x = x0
    y = 0.0_dp
    if (x <= 0.0_dp) then
      y = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    do while (x < 8.0_dp)
      y = y + 1.0_dp/(x*x)
      x = x + 1.0_dp
    end do
    r = 1.0_dp/x
    y = y + r + 0.5_dp*r*r + r**3/6.0_dp - r**5/30.0_dp + r**7/42.0_dp
  end function trigamma_mev

  real(dp) function rng_normal() result(z)
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function rng_normal

  recursive real(dp) function rng_gamma(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d, c, z, u, v
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      x = rng_gamma(shape+1.0_dp, scale)*u**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      z = rng_normal()
      v = 1.0_dp + c*z
      if (v <= 0.0_dp) cycle
      v = v**3
      call random_number(u)
      if (u < 1.0_dp-0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
    end do
    x = scale*d*v
  end function rng_gamma

  real(dp) function rng_beta(a, b) result(x)
    real(dp), intent(in) :: a, b
    real(dp) :: u, v
    u = rng_gamma(a,1.0_dp)
    v = rng_gamma(b,1.0_dp)
    x = u/(u+v)
  end function rng_beta

  pure subroutine sort_ascending(x, y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(size(x))
    real(dp) :: key
    integer :: i, j
    y = x
    do i = 2, size(y)
      key = y(i)
      j = i-1
      do while (j >= 1)
        if (y(j) <= key) exit
        y(j+1) = y(j)
        j = j-1
      end do
      y(j+1) = key
    end do
  end subroutine sort_ascending

  pure subroutine sort_descending(x, y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(size(x))
    real(dp) :: key
    integer :: i, j
    y = x
    do i = 2, size(y)
      key = y(i)
      j = i-1
      do while (j >= 1)
        if (y(j) >= key) exit
        y(j+1) = y(j)
        j = j-1
      end do
      y(j+1) = key
    end do
  end subroutine sort_descending

  pure real(dp) function mean_real(x) result(v)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = sum(x)/real(size(x),dp)
    end if
  end function mean_real

  pure real(dp) function variance_real(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) < 2) then
      v = 0.0_dp
      return
    end if
    m = mean_real(x)
    v = sum((x-m)**2)/real(size(x)-1,dp)
  end function variance_real

  pure real(dp) function median_real(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: y(size(x))
    call sort_ascending(x,y)
    if (mod(size(y),2) == 1) then
      v = y((size(y)+1)/2)
    else
      v = 0.5_dp*(y(size(y)/2)+y(size(y)/2+1))
    end if
  end function median_real

  pure real(dp) function empirical_quantile(x,p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp) :: y(size(x)), h, frac
    integer :: j, n
    n = size(x)
    if (n == 0) then
      q = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    call sort_ascending(x,y)
    if (p <= 0.0_dp) then
      q = y(1)
    else if (p >= 1.0_dp) then
      q = y(n)
    else
      h = 1.0_dp + real(n-1,dp)*p
      j = floor(h)
      frac = h-real(j,dp)
      q = (1.0_dp-frac)*y(j)+frac*y(min(j+1,n))
    end if
  end function empirical_quantile

  pure function outer_product(x,y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i=1,size(x)
      a(i,:) = x(i)*y
    end do
  end function outer_product

  pure subroutine covariance_matrix(x, cov, unbiased)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: cov(size(x,2),size(x,2))
    logical, intent(in), optional :: unbiased
    real(dp) :: m(size(x,2)), denom
    integer :: i
    logical :: ub
    ub = .true.
    if (present(unbiased)) ub = unbiased
    m = sum(x,dim=1)/real(size(x,1),dp)
    cov = 0.0_dp
    do i=1,size(x,1)
      cov = cov + outer_product(x(i,:)-m,x(i,:)-m)
    end do
    if (ub .and. size(x,1)>1) then
      denom = real(size(x,1)-1,dp)
    else
      denom = real(max(size(x,1),1),dp)
    end if
    cov = cov/denom
  end subroutine covariance_matrix

  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(size(b))
    integer, intent(out), optional :: info
    real(dp) :: m(size(b),size(b)), rhs(size(b)), row(size(b))
    real(dp) :: fac, piv, tmp
    integer :: n, i, k, p
    n = size(b)
    m = a
    rhs = b
    if (present(info)) info=0
    do k=1,n-1
      p=k
      piv=abs(m(k,k))
      do i=k+1,n
        if (abs(m(i,k)) > piv) then
          p=i; piv=abs(m(i,k))
        end if
      end do
      if (piv <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(m)))) then
        x=0.0_dp
        if (present(info)) info=k
        return
      end if
      if (p /= k) then
        row=m(k,:); m(k,:)=m(p,:); m(p,:)=row
        tmp=rhs(k); rhs(k)=rhs(p); rhs(p)=tmp
      end if
      do i=k+1,n
        fac=m(i,k)/m(k,k)
        m(i,k)=0.0_dp
        m(i,k+1:n)=m(i,k+1:n)-fac*m(k,k+1:n)
        rhs(i)=rhs(i)-fac*rhs(k)
      end do
    end do
    if (abs(m(n,n)) <= epsilon(1.0_dp)) then
      x=0.0_dp
      if (present(info)) info=n
      return
    end if
    x(n)=rhs(n)/m(n,n)
    do i=n-1,1,-1
      x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
    end do
  end subroutine solve_linear

  subroutine inverse_matrix(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    integer, intent(out), optional :: info
    real(dp) :: e(size(a,1)), x(size(a,1))
    integer :: j, ier
    ainv=0.0_dp
    ier=0
    do j=1,size(a,1)
      e=0.0_dp; e(j)=1.0_dp
      call solve_linear(a,e,x,ier)
      if (ier /= 0) exit
      ainv(:,j)=x
    end do
    if (present(info)) info=ier
  end subroutine inverse_matrix

  subroutine chol_lower(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(size(a,1),size(a,2))
    integer, intent(out), optional :: info
    integer :: i,j,k,n
    real(dp) :: s
    n=size(a,1)
    l=0.0_dp
    if (present(info)) info=0
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1
          s=s-l(i,k)*l(j,k)
        end do
        if (i == j) then
          if (s <= 0.0_dp) then
            if (present(info)) info=i
            return
          end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine chol_lower

  subroutine logdet_spd(a, logdet, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: logdet
    integer, intent(out), optional :: info
    real(dp) :: l(size(a,1),size(a,2))
    integer :: i, ier
    call chol_lower(a,l,ier)
    if (present(info)) info=ier
    if (ier /= 0) then
      logdet=ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    logdet=0.0_dp
    do i=1,size(a,1)
      logdet=logdet+2.0_dp*log(l(i,i))
    end do
  end subroutine logdet_spd

  subroutine mvnormal_sample(n,mu,sigma,y,info)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu(:),sigma(:,:)
    real(dp), intent(out) :: y(n,size(mu))
    integer, intent(out), optional :: info
    real(dp) :: l(size(mu),size(mu)), z(size(mu))
    integer :: i,j,ier
    call chol_lower(sigma,l,ier)
    if (present(info)) info=ier
    if (ier /= 0) then
      y=ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    do i=1,n
      do j=1,size(mu)
        z(j)=rng_normal()
      end do
      y(i,:)=mu+matmul(l,z)
    end do
  end subroutine mvnormal_sample

  subroutine finite_diff_gradient(fun,x,g,step)
    procedure(scalar_objective) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(size(x))
    real(dp), intent(in), optional :: step
    real(dp) :: h, xp(size(x)), xm(size(x))
    integer :: j
    do j=1,size(x)
      h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x(j)))
      if (present(step)) h=step*max(1.0_dp,abs(x(j)))
      xp=x; xm=x
      xp(j)=xp(j)+h; xm(j)=xm(j)-h
      g(j)=(fun(xp)-fun(xm))/(2.0_dp*h)
    end do
  end subroutine finite_diff_gradient

  subroutine finite_diff_hessian(fun,x,hess,step)
    procedure(scalar_objective) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hess(size(x),size(x))
    real(dp), intent(in), optional :: step
    real(dp) :: hi,hj, fpp,fpm,fmp,fmm,f0
    real(dp) :: xx(size(x))
    integer :: i,j
    f0=fun(x)
    hess=0.0_dp
    do i=1,size(x)
      hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
      if (present(step)) hi=step*max(1.0_dp,abs(x(i)))
      xx=x;xx(i)=x(i)+hi;fpp=fun(xx)
      xx=x;xx(i)=x(i)-hi;fmm=fun(xx)
      hess(i,i)=(fpp-2.0_dp*f0+fmm)/(hi*hi)
      do j=i+1,size(x)
        hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
        if (present(step)) hj=step*max(1.0_dp,abs(x(j)))
        xx=x;xx(i)=x(i)+hi;xx(j)=x(j)+hj;fpp=fun(xx)
        xx=x;xx(i)=x(i)+hi;xx(j)=x(j)-hj;fpm=fun(xx)
        xx=x;xx(i)=x(i)-hi;xx(j)=x(j)+hj;fmp=fun(xx)
        xx=x;xx(i)=x(i)-hi;xx(j)=x(j)-hj;fmm=fun(xx)
        hess(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
        hess(j,i)=hess(i,j)
      end do
    end do
  end subroutine finite_diff_hessian

  subroutine pattern_minimize(fun,x,fval,info,maxiter,tol,initial_step)
    procedure(scalar_objective) :: fun
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: fval
    integer, intent(out), optional :: info
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol, initial_step
    real(dp) :: best, cand, step, xt(size(x)), scale(size(x)), tol_use
    integer :: it,j,mx
    logical :: improved
    mx=1000; if(present(maxiter)) mx=maxiter
    step=0.25_dp; if(present(initial_step)) step=initial_step
    tol_use=1.0e-7_dp; if(present(tol)) tol_use=tol
    scale=max(1.0_dp,abs(x))
    best=fun(x)
    do it=1,mx
      improved=.false.
      do j=1,size(x)
        xt=x;xt(j)=xt(j)+step*scale(j)
        cand=fun(xt)
        if (ieee_is_finite(cand) .and. cand<best) then
          x=xt;best=cand;improved=.true.;cycle
        end if
        xt=x;xt(j)=xt(j)-step*scale(j)
        cand=fun(xt)
        if (ieee_is_finite(cand) .and. cand<best) then
          x=xt;best=cand;improved=.true.
        end if
      end do
      if (.not.improved) step=0.5_dp*step
      scale=max(1.0_dp,abs(x))
      if (step < tol_use) exit
    end do
    fval=best
    if(present(info)) then
      if(it<=mx) then;info=0;else;info=1;end if
    end if
  end subroutine pattern_minimize

  pure real(dp) function student_t_pdf(x,df) result(v)
    real(dp), intent(in) :: x,df
    if (df <= 0.0_dp) then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      v=exp(log_gamma((df+1.0_dp)/2.0_dp)-log_gamma(df/2.0_dp) &
          -0.5_dp*log(df*pi)-0.5_dp*(df+1.0_dp)*log(1.0_dp+x*x/df))
    end if
  end function student_t_pdf

  pure real(dp) function student_t_cdf(x,df) result(p)
    real(dp), intent(in) :: x,df
    real(dp) :: z, ib
    if (df <= 0.0_dp) then
      p=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      p=0.5_dp;return
    end if
    z=df/(df+x*x)
    ib=reg_beta(z,df/2.0_dp,0.5_dp)
    if (x>0.0_dp) then
      p=1.0_dp-0.5_dp*ib
    else
      p=0.5_dp*ib
    end if
  end function student_t_cdf

  pure real(dp) function student_t_quantile(prob,df) result(x)
    real(dp), intent(in) :: prob,df
    real(dp) :: lo,hi,mid
    integer :: it
    if(prob<=0.0_dp)then;x=-ieee_value(0.0_dp,ieee_positive_inf);return;end if
    if(prob>=1.0_dp)then;x=ieee_value(0.0_dp,ieee_positive_inf);return;end if
    lo=-1.0_dp;hi=1.0_dp
    do while(student_t_cdf(lo,df)>prob);lo=2.0_dp*lo;end do
    do while(student_t_cdf(hi,df)<prob);hi=2.0_dp*hi;end do
    do it=1,160
      mid=0.5_dp*(lo+hi)
      if(student_t_cdf(mid,df)<prob)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function student_t_quantile

end module mev_math
