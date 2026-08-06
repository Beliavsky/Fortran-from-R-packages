! SPDX-License-Identifier: GPL-3.0-only
module mass_math
  use rrcov_kinds, only : dp
  use rrcov_types, only : rrcov_success
  use rrcov_linalg, only : general_inverse, solve_linear, symmetric_inverse, identity_matrix
  use rrcov_stats, only : normal_cdf, regularized_beta, regularized_gamma_p
  use mass_types, only : mass_success, mass_invalid_argument, mass_no_convergence
  implicit none
  private
  public :: pi_dp, normal_pdf, normal_quantile, logistic_cdf, logistic_pdf
  public :: student_t_pdf, student_t_cdf, chi_square_cdf_mass
  public :: digamma_mass, trigamma_mass, log1pexp, softplus
  public :: least_squares, weighted_least_squares, matrix_pseudoinverse
  public :: bfgs_minimize, numerical_hessian, covariance_from_hessian
  public :: sample_variance, sample_sd, type7_quantile, median_mass, mad_mass
  public :: sort_real_mass, rank_average, isotonic_increasing

  real(dp), parameter :: pi_dp = acos(-1.0_dp)

  abstract interface
    function objective_function(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function objective_function
  end interface

contains

  pure elemental function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = exp(-0.5_dp * x * x) / sqrt(2.0_dp * pi_dp)
  end function normal_pdf

  pure elemental function logistic_cdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    if (x >= 0.0_dp) then
      value = 1.0_dp / (1.0_dp + exp(-x))
    else
      value = exp(x) / (1.0_dp + exp(x))
    end if
  end function logistic_cdf

  pure elemental function logistic_pdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value, p
    p = logistic_cdf(x)
    value = p * (1.0_dp - p)
  end function logistic_pdf

  pure elemental function log1pexp(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    if (x > 35.0_dp) then
      value = x
    else if (x < -35.0_dp) then
      value = exp(x)
    else
      value = log(1.0_dp + exp(x))
    end if
  end function log1pexp

  pure elemental function softplus(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = log1pexp(x)
  end function softplus

  function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
      4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
      7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp ]
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  function student_t_pdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    real(dp) :: value
    if (nu <= 0.0_dp) then
      value = 0.0_dp
    else
      value = exp(log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu)) / &
        sqrt(nu*pi_dp) * (1.0_dp + x*x/nu)**(-0.5_dp*(nu+1.0_dp))
    end if
  end function student_t_pdf

  function student_t_cdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    real(dp) :: value, z, ib
    if (nu <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      value = 0.5_dp
      return
    end if
    z = nu / (nu + x*x)
    ib = regularized_beta(z, 0.5_dp*nu, 0.5_dp)
    if (x > 0.0_dp) then
      value = 1.0_dp - 0.5_dp*ib
    else
      value = 0.5_dp*ib
    end if
  end function student_t_cdf

  function chi_square_cdf_mass(x, df) result(value)
    real(dp), intent(in) :: x, df
    real(dp) :: value
    if (x <= 0.0_dp .or. df <= 0.0_dp) then
      value = 0.0_dp
    else
      value = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
    end if
  end function chi_square_cdf_mass

  pure elemental function digamma_mass(xin) result(value)
    real(dp), intent(in) :: xin
    real(dp) :: value, x, inv, inv2
    x = xin
    value = 0.0_dp
    if (x <= 0.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    do while (x < 8.0_dp)
      value = value - 1.0_dp/x
      x = x + 1.0_dp
    end do
    inv = 1.0_dp/x
    inv2 = inv*inv
    value = value + log(x) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - &
      inv2*(1.0_dp/120.0_dp - inv2*(1.0_dp/252.0_dp - inv2/240.0_dp)))
  end function digamma_mass

  pure elemental function trigamma_mass(xin) result(value)
    real(dp), intent(in) :: xin
    real(dp) :: value, x, inv, inv2
    x = xin
    value = 0.0_dp
    if (x <= 0.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    do while (x < 8.0_dp)
      value = value + 1.0_dp/(x*x)
      x = x + 1.0_dp
    end do
    inv = 1.0_dp/x
    inv2 = inv*inv
    value = value + inv + 0.5_dp*inv2 + inv*inv2/6.0_dp - &
      inv*inv2*inv2/30.0_dp + inv*inv2*inv2*inv2/42.0_dp - &
      inv*inv2*inv2*inv2*inv2/30.0_dp
  end function trigamma_mass

  subroutine least_squares(x, y, beta, residuals, rank, status)
    real(dp), intent(in) :: x(:, :), y(:)
    real(dp), allocatable, intent(out) :: beta(:), residuals(:)
    integer, intent(out) :: rank, status
    real(dp), allocatable :: xtx(:, :), xty(:), pinv(:, :)
    integer :: p, st
    if (size(x,1) /= size(y) .or. size(x,1) == 0 .or. size(x,2) == 0) then
      allocate(beta(0), residuals(0)); rank = 0; status = mass_invalid_argument; return
    end if
    p = size(x,2)
    xtx = matmul(transpose(x), x)
    xty = matmul(transpose(x), y)
    pinv = general_inverse(xtx, st)
    beta = matmul(pinv, xty)
    residuals = y - matmul(x, beta)
    rank = count(abs([(xtx(p,p),p=1,size(xtx,1))]) > sqrt(epsilon(1.0_dp)))
    if (st == rrcov_success) then
      status = mass_success
    else
      status = mass_no_convergence
    end if
  end subroutine least_squares

  subroutine weighted_least_squares(x, y, w, beta, residuals, covariance, rank, status)
    real(dp), intent(in) :: x(:, :), y(:), w(:)
    real(dp), allocatable, intent(out) :: beta(:), residuals(:), covariance(:, :)
    integer, intent(out) :: rank, status
    real(dp), allocatable :: xw(:, :), yw(:), xtx(:, :), xty(:)
    real(dp) :: sw
    integer :: i, st
    if (size(x,1) /= size(y) .or. size(w) /= size(y) .or. any(w < 0.0_dp)) then
      allocate(beta(0), residuals(0), covariance(0,0)); rank=0; status=mass_invalid_argument; return
    end if
    allocate(xw(size(x,1),size(x,2)), yw(size(y)))
    do i=1,size(y)
      sw = sqrt(w(i))
      xw(i,:) = sw*x(i,:)
      yw(i) = sw*y(i)
    end do
    xtx = matmul(transpose(xw), xw)
    xty = matmul(transpose(xw), yw)
    covariance = general_inverse(xtx, st)
    beta = matmul(covariance, xty)
    residuals = y - matmul(x,beta)
    rank = size(beta)
    status = merge(mass_success,mass_no_convergence,st==rrcov_success)
  end subroutine weighted_least_squares

  function matrix_pseudoinverse(a, tolerance, status) result(value)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: value(:, :)
    integer :: st
    value = general_inverse(a, st, tolerance)
    if (present(status)) status = merge(mass_success,mass_no_convergence,st==rrcov_success)
  end function matrix_pseudoinverse

  subroutine bfgs_minimize(fun, x, fval, status, iterations, maxit, tolerance)
    procedure(objective_function) :: fun
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: fval
    integer, intent(out) :: status, iterations
    integer, intent(in), optional :: maxit
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: h(:, :), g(:), gnew(:), p(:), xnew(:), s(:), y(:)
    real(dp) :: tol, alpha, fnew, c1, sy, rho
    integer :: n, i, it, mit
    n = size(x)
    if (n == 0) then
      fval = fun(x); status=mass_invalid_argument; iterations=0; return
    end if
    mit=300; if (present(maxit)) mit=max(1,maxit)
    tol=1.0e-7_dp; if (present(tolerance)) tol=max(tolerance,epsilon(1.0_dp))
    allocate(h(n,n),g(n),gnew(n),p(n),xnew(n),s(n),y(n))
    h=identity_matrix(n)
    call numerical_gradient(fun,x,g)
    fval=fun(x); c1=1.0e-4_dp
    do it=1,mit
      if (maxval(abs(g)) <= tol*(1.0_dp+abs(fval))) exit
      p=-matmul(h,g)
      if (dot_product(p,g) >= 0.0_dp) p=-g
      alpha=1.0_dp
      do i=1,40
        xnew=x+alpha*p
        fnew=fun(xnew)
        if (fnew <= fval+c1*alpha*dot_product(g,p)) exit
        alpha=0.5_dp*alpha
      end do
      if (alpha < 1.0e-12_dp) exit
      call numerical_gradient(fun,xnew,gnew)
      s=xnew-x; y=gnew-g; sy=dot_product(s,y)
      if (sy > 1.0e-12_dp*sqrt(max(dot_product(s,s)*dot_product(y,y),tiny(1.0_dp)))) then
        rho=1.0_dp/sy
        h = matmul(identity_matrix(n)-rho*outer(s,y), &
             matmul(h,identity_matrix(n)-rho*outer(y,s))) + rho*outer(s,s)
      else
        h=identity_matrix(n)
      end if
      x=xnew; g=gnew; fval=fnew
    end do
    iterations=min(it,mit)
    if (maxval(abs(g)) <= 10.0_dp*tol*(1.0_dp+abs(fval))) then
      status=mass_success
    else
      status=mass_no_convergence
    end if
  contains
    pure function outer(a,b) result(c)
      real(dp), intent(in) :: a(:),b(:)
      real(dp) :: c(size(a),size(b))
      integer :: k
      do k=1,size(a); c(k,:)=a(k)*b; end do
    end function outer
  end subroutine bfgs_minimize

  subroutine numerical_gradient(fun, x, g)
    procedure(objective_function) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    real(dp), allocatable :: xp(:), xm(:)
    real(dp) :: h
    integer :: i
    allocate(xp(size(x)),xm(size(x)))
    do i=1,size(x)
      h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
      xp=x; xm=x; xp(i)=xp(i)+h; xm(i)=xm(i)-h
      g(i)=(fun(xp)-fun(xm))/(2.0_dp*h)
    end do
  end subroutine numerical_gradient

  subroutine numerical_hessian(fun, x, hessian)
    procedure(objective_function) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: hessian(:, :)
    real(dp), allocatable :: xpp(:),xpm(:),xmp(:),xmm(:)
    real(dp) :: hi,hj,f0
    integer :: i,j,n
    n=size(x); allocate(hessian(n,n),xpp(n),xpm(n),xmp(n),xmm(n)); f0=fun(x)
    do i=1,n
      hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
      xpp=x; xpm=x; xpp(i)=x(i)+hi; xpm(i)=x(i)-hi
      hessian(i,i)=(fun(xpp)-2.0_dp*f0+fun(xpm))/(hi*hi)
      do j=i+1,n
        hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=x(i)+hi; xpp(j)=x(j)+hj
        xpm(i)=x(i)+hi; xpm(j)=x(j)-hj
        xmp(i)=x(i)-hi; xmp(j)=x(j)+hj
        xmm(i)=x(i)-hi; xmm(j)=x(j)-hj
        hessian(i,j)=(fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*hi*hj)
        hessian(j,i)=hessian(i,j)
      end do
    end do
  end subroutine numerical_hessian

  function covariance_from_hessian(hessian, status) result(covariance)
    real(dp), intent(in) :: hessian(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: covariance(:, :)
    integer :: st
    covariance=symmetric_inverse(hessian,st)
    status=merge(mass_success,mass_no_convergence,st==rrcov_success)
  end function covariance_from_hessian

  pure function sample_variance(x, unbiased) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: unbiased
    real(dp) :: value,mu,den
    logical :: ub
    ub=.true.; if(present(unbiased)) ub=unbiased
    if(size(x)<1) then; value=0.0_dp; return; end if
    mu=sum(x)/real(size(x),dp)
    den=real(size(x),dp); if(ub .and. size(x)>1) den=real(size(x)-1,dp)
    value=sum((x-mu)**2)/den
  end function sample_variance

  pure function sample_sd(x, unbiased) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: unbiased
    real(dp) :: value
    value=sqrt(max(0.0_dp,sample_variance(x,unbiased)))
  end function sample_sd

  function type7_quantile(x,p) result(value)
    real(dp), intent(in) :: x(:),p
    real(dp) :: value,h,g
    real(dp), allocatable :: w(:)
    integer :: j,n
    n=size(x); if(n==0) then; value=0.0_dp; return; end if
    w=x; call sort_real_mass(w)
    if(p<=0.0_dp) then; value=w(1); return; else if(p>=1.0_dp) then; value=w(n); return; end if
    h=1.0_dp+real(n-1,dp)*p; j=floor(h); g=h-real(j,dp)
    if(j>=n) then; value=w(n); else; value=(1.0_dp-g)*w(j)+g*w(j+1); end if
  end function type7_quantile

  function median_mass(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value=type7_quantile(x,0.5_dp)
  end function median_mass

  function mad_mass(x, constant) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: constant
    real(dp) :: value,c,m
    c=1.482602218505602_dp; if(present(constant)) c=constant
    m=median_mass(x); value=c*median_mass(abs(x-m))
  end function mad_mass

  subroutine sort_real_mass(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real_mass

  subroutine rank_average(x,ranks)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: ranks(size(x))
    integer, allocatable :: idx(:)
    real(dp), allocatable :: w(:)
    real(dp) :: r
    integer :: i,j,k,n,tmp
    n=size(x); allocate(idx(n),w(n)); idx=[(i,i=1,n)]; w=x
    do i=2,n
      r=w(i); tmp=idx(i); j=i-1
      do while(j>=1)
        if(w(j)<=r) exit
        w(j+1)=w(j); idx(j+1)=idx(j); j=j-1
      end do
      w(j+1)=r; idx(j+1)=tmp
    end do
    i=1
    do while(i<=n)
      j=i
      do while(j<n)
        if(abs(w(j+1)-w(i)) > epsilon(1.0_dp)*max(1.0_dp,abs(w(i)))) exit
        j=j+1
      end do
      r=0.5_dp*real(i+j,dp)
      do k=i,j; ranks(idx(k))=r; end do
      i=j+1
    end do
  end subroutine rank_average

  subroutine isotonic_increasing(y,w,fit)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in), optional :: w(:)
    real(dp), intent(out) :: fit(size(y))
    real(dp), allocatable :: level(:),weight(:)
    integer, allocatable :: start(:),finish(:)
    integer :: n,m,i,j
    n=size(y); allocate(level(n),weight(n),start(n),finish(n)); m=0
    do i=1,n
      m=m+1; level(m)=y(i); weight(m)=1.0_dp; if(present(w)) weight(m)=w(i)
      start(m)=i; finish(m)=i
      do while(m>1)
        if(level(m-1)<=level(m)) exit
        level(m-1)=(weight(m-1)*level(m-1)+weight(m)*level(m))/(weight(m-1)+weight(m))
        weight(m-1)=weight(m-1)+weight(m); finish(m-1)=finish(m); m=m-1
      end do
    end do
    do i=1,m
      do j=start(i),finish(i); fit(j)=level(i); end do
    end do
  end subroutine isotonic_increasing

end module mass_math
