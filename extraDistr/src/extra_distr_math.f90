! SPDX-License-Identifier: GPL-2.0-only
module extra_distr_math
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf, ieee_is_finite, ieee_is_nan
  use extra_distr_kinds, only : dp, pi, sqrt_two_pi
  implicit none
  private

  public :: nan_dp, pos_inf, neg_inf, is_integer_value
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: regularized_beta, beta_pdf, beta_quantile
  public :: regularized_gamma_p, gamma_pdf, gamma_quantile
  public :: student_pdf, student_cdf, student_quantile
  public :: log_choose, binom_pmf, binom_cdf, binom_quantile
  public :: poisson_pmf, poisson_cdf, poisson_quantile
  public :: nbinom_pmf, nbinom_cdf, nbinom_quantile
  public :: log_beta, log_factorial, apply_density_log, apply_tail, decode_probability
  public :: bessel_i_integer, clamp_probability, log_sum_exp

contains

  pure real(dp) function nan_dp() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  pure real(dp) function pos_inf() result(x)
    x = ieee_value(0.0_dp, ieee_positive_inf)
  end function pos_inf

  pure real(dp) function neg_inf() result(x)
    x = ieee_value(0.0_dp, ieee_negative_inf)
  end function neg_inf

  pure logical function is_integer_value(x) result(ok)
    real(dp), intent(in) :: x
    ok = ieee_is_finite(x) .and. x == aint(x)
  end function is_integer_value

  pure real(dp) function clamp_probability(p) result(value)
    real(dp), intent(in) :: p
    if (ieee_is_nan(p)) then
      value = p
    else
      value = min(1.0_dp, max(0.0_dp, p))
    end if
  end function clamp_probability


  pure real(dp) function decode_probability(prob, lower_tail, log_p) result(p)
    real(dp), intent(in) :: prob
    logical, intent(in), optional :: lower_tail, log_p
    logical :: lower, use_log
    lower = .true.
    use_log = .false.
    if (present(lower_tail)) lower = lower_tail
    if (present(log_p)) use_log = log_p
    p = prob
    if (use_log) p = exp(p)
    if (.not. lower) p = 1.0_dp-p
  end function decode_probability

  pure real(dp) function apply_density_log(value, log_p) result(out)
    real(dp), intent(in) :: value
    logical, intent(in), optional :: log_p
    logical :: use_log
    use_log = .false.
    if (present(log_p)) use_log = log_p
    if (use_log) then
      if (value > 0.0_dp) then
        out = log(value)
      else if (value == 0.0_dp) then
        out = neg_inf()
      else
        out = nan_dp()
      end if
    else
      out = value
    end if
  end function apply_density_log

  pure real(dp) function apply_tail(value, lower_tail, log_p) result(out)
    real(dp), intent(in) :: value
    logical, intent(in), optional :: lower_tail, log_p
    logical :: lower, use_log
    real(dp) :: p
    lower = .true.
    use_log = .false.
    if (present(lower_tail)) lower = lower_tail
    if (present(log_p)) use_log = log_p
    p = value
    if (.not. lower) p = 1.0_dp-p
    p = clamp_probability(p)
    if (use_log) then
      if (p > 0.0_dp) then
        out = log(p)
      else
        out = neg_inf()
      end if
    else
      out = p
    end if
  end function apply_tail

  pure elemental real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp*x*x)/sqrt_two_pi
  end function normal_pdf

  pure elemental real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, -2.759285104469687e2_dp, &
       1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, -1.556989798598866e2_dp, &
       6.680131188771972e1_dp, -1.328068155288572e1_dp]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, -2.400758277161838_dp, &
      -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp]
    real(dp), parameter :: plow = 0.02425_dp
    real(dp), parameter :: phigh = 1.0_dp-plow
    if (p <= 0.0_dp) then
      x = neg_inf()
    else if (p >= 1.0_dp) then
      x = pos_inf()
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    if (ieee_is_finite(x)) x = x-(normal_cdf(x)-p)/max(normal_pdf(x),tiny(1.0_dp))
  end function normal_quantile

  pure real(dp) function log_beta(a,b) result(value)
    real(dp), intent(in) :: a,b
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      value = nan_dp()
    else
      value = log_gamma(a)+log_gamma(b)-log_gamma(a+b)
    end if
  end function log_beta

  pure real(dp) function betacf(a,b,x) result(value)
    real(dp), intent(in) :: a,b,x
    integer :: m,m2
    real(dp) :: aa,c,d,del,h,qab,qam,qap
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp; d=1.0_dp-qab*x/qap
    if(abs(d)<fpmin) d=fpmin
    d=1.0_dp/d; h=d
    do m=1,400
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<=eps) exit
    end do
    value=h
  end function betacf

  pure real(dp) function regularized_beta(x,a,b) result(value)
    real(dp), intent(in) :: x,a,b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      value=nan_dp(); return
    else if(x<=0.0_dp) then
      value=0.0_dp; return
    else if(x>=1.0_dp) then
      value=1.0_dp; return
    end if
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if(x<(a+1.0_dp)/(a+b+2.0_dp)) then
      value=bt*betacf(a,b,x)/a
    else
      value=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
    end if
    value=clamp_probability(value)
  end function regularized_beta

  pure real(dp) function beta_pdf(x,a,b) result(value)
    real(dp), intent(in) :: x,a,b
    if(a<=0.0_dp .or. b<=0.0_dp) then
      value=nan_dp()
    else if(x<0.0_dp .or. x>1.0_dp) then
      value=0.0_dp
    else if(x==0.0_dp) then
      if(a<1.0_dp) then; value=pos_inf(); else if(a==1.0_dp) then; value=b; else; value=0.0_dp; end if
    else if(x==1.0_dp) then
      if(b<1.0_dp) then; value=pos_inf(); else if(b==1.0_dp) then; value=a; else; value=0.0_dp; end if
    else
      value=exp((a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x)-log_beta(a,b))
    end if
  end function beta_pdf

  real(dp) function beta_quantile(p,a,b) result(value)
    real(dp), intent(in) :: p,a,b
    real(dp) :: lo,hi,mid
    integer :: iter
    if(a<=0.0_dp .or. b<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp) then
      value=nan_dp(); return
    else if(p==0.0_dp) then; value=0.0_dp; return
    else if(p==1.0_dp) then; value=1.0_dp; return
    end if
    lo=0.0_dp; hi=1.0_dp
    do iter=1,160
      mid=0.5_dp*(lo+hi)
      if(regularized_beta(mid,a,b)<p) then; lo=mid; else; hi=mid; end if
      if(hi-lo <= 2.0e-13_dp) exit
    end do
    value=0.5_dp*(lo+hi)
  end function beta_quantile

  pure real(dp) function regularized_gamma_p(a,x) result(value)
    real(dp), intent(in) :: a,x
    real(dp), parameter :: eps=3.0e-14_dp,fpmin=1.0e-300_dp
    real(dp)::ap,del,sumv,b,c,d,h,an
    integer::n
    if(a<=0.0_dp .or. x<0.0_dp) then; value=nan_dp(); return
    else if(x==0.0_dp) then; value=0.0_dp; return
    end if
    if(x<a+1.0_dp) then
      ap=a; sumv=1.0_dp/a; del=sumv
      do n=1,500
        ap=ap+1.0_dp; del=del*x/ap; sumv=sumv+del
        if(abs(del)<=abs(sumv)*eps) exit
      end do
      value=sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/max(abs(b),fpmin)
      if(b<0.0_dp)d=-d
      h=d
      do n=1,500
        an=-real(n,dp)*(real(n,dp)-a); b=b+2.0_dp
        d=an*d+b; if(abs(d)<fpmin)d=fpmin
        c=b+an/c; if(abs(c)<fpmin)c=fpmin
        d=1.0_dp/d; del=d*c; h=h*del
        if(abs(del-1.0_dp)<=eps) exit
      end do
      value=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    value=clamp_probability(value)
  end function regularized_gamma_p

  pure real(dp) function gamma_pdf(x,shape,scale) result(value)
    real(dp), intent(in)::x,shape,scale
    if(shape<=0.0_dp .or. scale<=0.0_dp) then
      value=nan_dp()
    else if(x<0.0_dp) then
      value=0.0_dp
    else if(x==0.0_dp) then
      if(shape<1.0_dp) then; value=pos_inf(); else if(shape==1.0_dp) then; value=1.0_dp/scale; else; value=0.0_dp; end if
    else
      value=exp((shape-1.0_dp)*log(x)-x/scale-log_gamma(shape)-shape*log(scale))
    end if
  end function gamma_pdf

  real(dp) function gamma_quantile(p,shape,scale) result(value)
    real(dp), intent(in)::p,shape,scale
    real(dp)::lo,hi,mid
    integer::iter
    if(shape<=0.0_dp .or. scale<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp) then; value=nan_dp(); return
    else if(p==0.0_dp) then; value=0.0_dp; return
    else if(p==1.0_dp) then; value=pos_inf(); return
    end if
    lo=0.0_dp; hi=max(1.0_dp,shape+8.0_dp*sqrt(shape)+8.0_dp)
    do while(regularized_gamma_p(shape,hi)<p .and. hi<1.0e300_dp); hi=2.0_dp*hi; end do
    do iter=1,180
      mid=0.5_dp*(lo+hi)
      if(regularized_gamma_p(shape,mid)<p) then; lo=mid; else; hi=mid; end if
      if(abs(hi-lo)<=2.0e-12_dp*max(1.0_dp,abs(mid))) exit
    end do
    value=scale*0.5_dp*(lo+hi)
  end function gamma_quantile

  pure real(dp) function student_pdf(x,nu) result(value)
    real(dp), intent(in)::x,nu
    if(nu<=0.0_dp) then; value=nan_dp(); else
      value=exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)-0.5_dp*log(pi*nu)- &
        0.5_dp*(nu+1.0_dp)*log(1.0_dp+x*x/nu))
    end if
  end function student_pdf

  pure real(dp) function student_cdf(x,nu) result(value)
    real(dp), intent(in)::x,nu
    real(dp)::z
    if(nu<=0.0_dp) then; value=nan_dp(); return
    else if(x==0.0_dp) then; value=0.5_dp; return
    end if
    z=nu/(nu+x*x)
    if(x>0.0_dp) then; value=1.0_dp-0.5_dp*regularized_beta(z,0.5_dp*nu,0.5_dp)
    else; value=0.5_dp*regularized_beta(z,0.5_dp*nu,0.5_dp); end if
  end function student_cdf

  real(dp) function student_quantile(p,nu) result(value)
    real(dp), intent(in)::p,nu
    real(dp)::lo,hi,mid
    integer::iter
    if(nu<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp) then; value=nan_dp(); return
    else if(p==0.0_dp) then; value=neg_inf(); return
    else if(p==1.0_dp) then; value=pos_inf(); return
    end if
    lo=-1.0_dp; hi=1.0_dp
    do while(student_cdf(lo,nu)>p); lo=2.0_dp*lo; end do
    do while(student_cdf(hi,nu)<p); hi=2.0_dp*hi; end do
    do iter=1,180
      mid=0.5_dp*(lo+hi)
      if(student_cdf(mid,nu)<p) then; lo=mid; else; hi=mid; end if
      if(abs(hi-lo)<=2.0e-12_dp*max(1.0_dp,abs(mid))) exit
    end do
    value=0.5_dp*(lo+hi)
  end function student_quantile

  pure real(dp) function log_factorial(n) result(value)
    integer, intent(in)::n
    if(n<0) then; value=nan_dp(); else; value=log_gamma(real(n+1,dp)); end if
  end function log_factorial

  pure real(dp) function log_choose(n,k) result(value)
    integer, intent(in)::n,k
    if(n<0 .or. k<0 .or. k>n) then; value=neg_inf(); else
      value=log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))
    end if
  end function log_choose

  pure real(dp) function binom_pmf(k,n,p) result(value)
    integer, intent(in)::k,n
    real(dp), intent(in)::p
    if(n<0 .or. p<0.0_dp .or. p>1.0_dp) then; value=nan_dp()
    else if(k<0 .or. k>n) then; value=0.0_dp
    else if(p==0.0_dp) then; value=merge(1.0_dp,0.0_dp,k==0)
    else if(p==1.0_dp) then; value=merge(1.0_dp,0.0_dp,k==n)
    else; value=exp(log_choose(n,k)+real(k,dp)*log(p)+real(n-k,dp)*log(1.0_dp-p)); end if
  end function binom_pmf

  pure real(dp) function binom_cdf(k,n,p) result(value)
    integer, intent(in)::k,n
    real(dp), intent(in)::p
    if(n<0 .or. p<0.0_dp .or. p>1.0_dp) then; value=nan_dp()
    else if(k<0) then; value=0.0_dp
    else if(k>=n) then; value=1.0_dp
    else; value=regularized_beta(1.0_dp-p,real(n-k,dp),real(k+1,dp)); end if
  end function binom_cdf

  pure integer function binom_quantile(prob,n,p) result(k)
    real(dp), intent(in)::prob,p
    integer, intent(in)::n
    integer::lo,hi,mid
    if(prob<=0.0_dp) then; k=0; return
    else if(prob>=1.0_dp) then; k=n; return
    end if
    lo=0; hi=n
    do while(lo<hi)
      mid=(lo+hi)/2
      if(binom_cdf(mid,n,p)>=prob) then; hi=mid; else; lo=mid+1; end if
    end do
    k=lo
  end function binom_quantile

  pure real(dp) function poisson_pmf(k,lambda) result(value)
    integer, intent(in)::k
    real(dp), intent(in)::lambda
    if(lambda<0.0_dp) then; value=nan_dp()
    else if(k<0) then; value=0.0_dp
    else if(lambda==0.0_dp) then; value=merge(1.0_dp,0.0_dp,k==0)
    else; value=exp(real(k,dp)*log(lambda)-lambda-log_factorial(k)); end if
  end function poisson_pmf

  pure real(dp) function poisson_cdf(k,lambda) result(value)
    integer, intent(in)::k
    real(dp), intent(in)::lambda
    if(lambda<0.0_dp) then; value=nan_dp()
    else if(k<0) then; value=0.0_dp
    else if(lambda==0.0_dp) then; value=1.0_dp
    else; value=1.0_dp-regularized_gamma_p(real(k+1,dp),lambda); end if
  end function poisson_cdf

  pure integer function poisson_quantile(prob,lambda) result(k)
    real(dp), intent(in)::prob,lambda
    integer::lo,hi,mid
    if(prob<=0.0_dp) then; k=0; return
    else if(prob>=1.0_dp) then; k=huge(1); return
    end if
    lo=0; hi=max(1,int(lambda+10.0_dp*sqrt(max(lambda,1.0_dp))+10.0_dp))
    do while(poisson_cdf(hi,lambda)<prob .and. hi<1073741823); hi=2*hi; end do
    do while(lo<hi)
      mid=lo+(hi-lo)/2
      if(poisson_cdf(mid,lambda)>=prob) then; hi=mid; else; lo=mid+1; end if
    end do
    k=lo
  end function poisson_quantile

  pure real(dp) function nbinom_pmf(k,size,p) result(value)
    integer, intent(in)::k
    real(dp), intent(in)::size,p
    if(size<=0.0_dp .or. p<=0.0_dp .or. p>1.0_dp) then; value=nan_dp()
    else if(k<0) then; value=0.0_dp
    else if(p==1.0_dp) then; value=merge(1.0_dp,0.0_dp,k==0)
    else; value=exp(log_gamma(real(k,dp)+size)-log_gamma(size)-log_gamma(real(k+1,dp))+ &
      size*log(p)+real(k,dp)*log(1.0_dp-p)); end if
  end function nbinom_pmf

  pure real(dp) function nbinom_cdf(k,size,p) result(value)
    integer, intent(in)::k
    real(dp), intent(in)::size,p
    if(size<=0.0_dp .or. p<=0.0_dp .or. p>1.0_dp) then; value=nan_dp()
    else if(k<0) then; value=0.0_dp
    else; value=regularized_beta(p,size,real(k+1,dp)); end if
  end function nbinom_cdf

  pure integer function nbinom_quantile(prob,size,p) result(k)
    real(dp), intent(in)::prob,size,p
    integer::lo,hi,mid
    if(prob<=0.0_dp) then; k=0; return
    else if(prob>=1.0_dp) then; k=huge(1); return
    end if
    lo=0; hi=max(1,int(size*(1.0_dp-p)/p+10.0_dp*sqrt(size*(1.0_dp-p)/(p*p))+10.0_dp))
    do while(nbinom_cdf(hi,size,p)<prob .and. hi<1073741823); hi=2*hi; end do
    do while(lo<hi)
      mid=lo+(hi-lo)/2
      if(nbinom_cdf(mid,size,p)>=prob) then; hi=mid; else; lo=mid+1; end if
    end do
    k=lo
  end function nbinom_quantile

  pure real(dp) function bessel_i_integer(n,x) result(value)
    integer, intent(in)::n
    real(dp), intent(in)::x
    integer::k,nn
    real(dp)::term,sumv,half
    nn=abs(n)
    if(x==0.0_dp) then; value=merge(1.0_dp,0.0_dp,nn==0); return; end if
    half=0.5_dp*x
    term=exp(real(nn,dp)*log(abs(half))-log_gamma(real(nn+1,dp)))
    sumv=term
    do k=1,10000
      term=term*(half*half)/(real(k,dp)*real(k+nn,dp))
      sumv=sumv+term
      if(abs(term)<=epsilon(sumv)*abs(sumv)) exit
    end do
    value=sumv
    if(x<0.0_dp .and. mod(nn,2)==1) value=-value
  end function bessel_i_integer

  pure real(dp) function log_sum_exp(x) result(value)
    real(dp), intent(in)::x(:)
    real(dp)::m
    if(size(x)==0) then; value=neg_inf(); return; end if
    m=maxval(x)
    if(.not.ieee_is_finite(m)) then; value=m; else; value=m+log(sum(exp(x-m))); end if
  end function log_sum_exp

end module extra_distr_math
