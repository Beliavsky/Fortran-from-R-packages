! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_special
  use survey_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: regularized_beta, regularized_gamma_q, f_survival, chisq_survival
  public :: weighted_chisq_survival, weighted_f_survival, saddle_survival_df
contains
  real(dp) function regularized_beta(x,a,b) result(v)
    real(dp), intent(in) :: x,a,b
    real(dp) :: bt
    if(a<=0.0_dp .or. b<=0.0_dp) error stop 'regularized_beta: shape parameters must be positive'
    if(x<=0.0_dp) then
      v=0.0_dp; return
    else if(x>=1.0_dp) then
      v=1.0_dp; return
    end if
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if(x < (a+1.0_dp)/(a+b+2.0_dp)) then
      v=bt*betacf(a,b,x)/a
    else
      v=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
    end if
    v=max(0.0_dp,min(1.0_dp,v))
  end function regularized_beta

  real(dp) function betacf(a,b,x) result(h)
    real(dp), intent(in) :: a,b,x
    integer, parameter :: maxit=300
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    integer :: m,m2
    real(dp) :: aa,c,d,del,qab,qam,qap
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp; d=1.0_dp-qab*x/qap
    if(abs(d)<fpmin) d=fpmin
    d=1.0_dp/d; h=d
    do m=1,maxit
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin; c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin; d=1.0_dp/d; h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin; c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin; d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<=eps) return
    end do
  end function betacf

  real(dp) function regularized_gamma_q(a,x) result(q)
    real(dp), intent(in) :: a,x
    integer, parameter :: maxit=500
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    integer :: i
    real(dp) :: ap,del,sum,b,c,d,h,an
    if(a<=0.0_dp) error stop 'regularized_gamma_q: a must be positive'
    if(x<=0.0_dp) then; q=1.0_dp; return; end if
    if(x<a+1.0_dp) then
      ap=a; sum=1.0_dp/a; del=sum
      do i=1,maxit
        ap=ap+1.0_dp; del=del*x/ap; sum=sum+del
        if(abs(del)<=abs(sum)*eps) exit
      end do
      q=1.0_dp-sum*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/b; h=d
      do i=1,maxit
        an=-real(i,dp)*(real(i,dp)-a); b=b+2.0_dp
        d=an*d+b; if(abs(d)<fpmin)d=fpmin
        c=b+an/c; if(abs(c)<fpmin)c=fpmin
        d=1.0_dp/d; del=d*c; h=h*del
        if(abs(del-1.0_dp)<=eps) exit
      end do
      q=exp(-x+a*log(x)-log_gamma(a))*h
    end if
    q=max(0.0_dp,min(1.0_dp,q))
  end function regularized_gamma_q

  real(dp) function f_survival(f,df1,df2) result(p)
    real(dp), intent(in) :: f,df1,df2
    real(dp) :: x
    if(f<=0.0_dp) then; p=1.0_dp; return; end if
    if(df1<=0.0_dp .or. df2<=0.0_dp) error stop 'f_survival: degrees of freedom must be positive'
    x=df2/(df2+df1*f)
    p=regularized_beta(x,0.5_dp*df2,0.5_dp*df1)
  end function f_survival

  real(dp) function chisq_survival(x,df) result(p)
    real(dp), intent(in) :: x,df
    if(df<=0.0_dp) error stop 'chisq_survival: df must be positive'
    if(x<=0.0_dp) then; p=1.0_dp; else; p=regularized_gamma_q(0.5_dp*df,0.5_dp*x); end if
  end function chisq_survival

  real(dp) function weighted_chisq_survival(x, lambda, saddlepoint) result(p)
    real(dp), intent(in) :: x, lambda(:)
    logical, intent(in), optional :: saddlepoint
    real(dp), allocatable :: dfs(:)
    real(dp) :: m1, m2, scale, ndf, psad
    logical :: dosad
    if (size(lambda) < 1) error stop 'weighted_chisq_survival: empty lambda'
    if (sum(lambda) <= 0.0_dp) error stop 'weighted_chisq_survival: positive mean required'
    allocate(dfs(size(lambda))); dfs = 1.0_dp
    m1 = sum(lambda)
    m2 = sum(lambda*lambda)
    scale = m2/m1
    ndf = m1*m1/m2
    p = chisq_survival(x/scale,ndf)
    dosad = .true.; if (present(saddlepoint)) dosad = saddlepoint
    if (dosad) then
      psad = saddle_survival_df(x,lambda,dfs)
      if (.not.ieee_is_nan(psad) .and. psad >= 0.0_dp .and. psad <= 1.0_dp) p = psad
    end if
  end function weighted_chisq_survival

  real(dp) function weighted_f_survival(x, lambda, ddf, saddlepoint) result(p)
    real(dp), intent(in) :: x, lambda(:), ddf
    logical, intent(in), optional :: saddlepoint
    real(dp), allocatable :: a(:), dfs(:)
    real(dp) :: m1, m2, scale, ndf, psad
    logical :: dosad
    integer :: m
    if (size(lambda) < 1 .or. ddf <= 0.0_dp) error stop 'weighted_f_survival: invalid arguments'
    m1 = sum(lambda)
    m2 = sum(lambda*lambda)
    if (m1 <= 0.0_dp .or. m2 <= 0.0_dp) error stop 'weighted_f_survival: positive weights required'
    scale = m2/m1
    ndf = m1*m1/m2
    p = f_survival(x/(ndf*scale),ndf,ddf)
    dosad = .true.; if (present(saddlepoint)) dosad = saddlepoint
    if (dosad) then
      m = size(lambda)
      allocate(a(m+1),dfs(m+1))
      a(1:m) = lambda; dfs(1:m) = 1.0_dp
      a(m+1) = -x/ddf; dfs(m+1) = ddf
      psad = saddle_survival_df(0.0_dp,a,dfs)
      if (.not.ieee_is_nan(psad) .and. psad >= 0.0_dp .and. psad <= 1.0_dp) p = psad
    end if
  end function weighted_f_survival

  real(dp) function saddle_survival_df(x, lambda, dfs) result(p)
    real(dp), intent(in) :: x, lambda(:), dfs(:)
    real(dp) :: d, xx, lo, hi, mid, flo, fhi, fm, zeta, k0, kpp, w, v, arg
    real(dp), allocatable :: lam(:)
    integer :: it
    if (size(lambda) /= size(dfs) .or. size(lambda) < 1 .or. any(dfs <= 0.0_dp)) then
      p = ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    d = maxval(abs(lambda))
    if (d <= tiny(1.0_dp)) then
      p = merge(1.0_dp,0.0_dp,x < 0.0_dp); return
    end if
    allocate(lam(size(lambda))); lam = lambda/d; xx = x/d
    if (any(lam < 0.0_dp)) then
      lo = maxval(pack(1.0_dp/(2.0_dp*lam),lam < 0.0_dp))*0.99999_dp
    else if (xx > sum(dfs*lam)) then
      lo = -0.01_dp
    else
      if (xx <= tiny(1.0_dp)) then
        lo = -1.0e4_dp
      else
        lo = -sum(dfs)/(2.0_dp*xx)
      end if
    end if
    if (.not.any(lam > 0.0_dp)) then
      hi = 1.0e4_dp
    else
      hi = minval(pack(1.0_dp/(2.0_dp*lam),lam > 0.0_dp))*0.99999_dp
    end if
    flo = saddle_kprime(lo,lam,dfs)-xx
    fhi = saddle_kprime(hi,lam,dfs)-xx
    if (ieee_is_nan(flo) .or. ieee_is_nan(fhi) .or. flo*fhi > 0.0_dp) then
      p = ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    do it = 1, 200
      mid = 0.5_dp*(lo+hi)
      fm = saddle_kprime(mid,lam,dfs)-xx
      if (abs(fm) < 1.0e-12_dp*(1.0_dp+abs(xx)) .or. abs(hi-lo) < 1.0e-11_dp*(1.0_dp+abs(mid))) exit
      if (flo*fm <= 0.0_dp) then
        hi = mid; fhi = fm
      else
        lo = mid; flo = fm
      end if
    end do
    zeta = 0.5_dp*(lo+hi)
    if (abs(zeta) < 1.0e-5_dp) then
      p = ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    k0 = -0.5_dp*sum(dfs*log(1.0_dp-2.0_dp*zeta*lam))
    kpp = 2.0_dp*sum(dfs*lam*lam/(1.0_dp-2.0_dp*zeta*lam)**2)
    if (2.0_dp*(zeta*xx-k0) <= 0.0_dp .or. kpp <= 0.0_dp) then
      p = ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    w = sign(1.0_dp,zeta)*sqrt(2.0_dp*(zeta*xx-k0))
    v = zeta*sqrt(kpp)
    if (abs(w) <= tiny(1.0_dp) .or. v/w <= 0.0_dp) then
      p = ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    arg = w+log(v/w)/w
    p = 0.5_dp*erfc(arg/sqrt(2.0_dp))
    p = max(0.0_dp,min(1.0_dp,p))
  end function saddle_survival_df

  pure real(dp) function saddle_kprime(zeta,lambda,dfs) result(v)
    real(dp), intent(in) :: zeta,lambda(:),dfs(:)
    v = sum(dfs*lambda/(1.0_dp-2.0_dp*zeta*lambda))
  end function saddle_kprime

end module survey_special
