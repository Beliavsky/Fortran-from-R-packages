! SPDX-License-Identifier: GPL-2.0-only
module compound_cox_math
  use compound_cox_kinds, only : dp
  implicit none
  private
  public :: normal_cdf, normal_quantile, chisq_cdf, chisq_quantile, median_real
  public :: invert_matrix, pseudoinverse_sym, symmetric_eigenvalues, covariance_matrix
  public :: bfgs_minimize, scalar_second_derivative, randn, shuffle_int

  abstract interface
    function objective_fn(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_fn
  end interface
contains
  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: q, r
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
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 1.0_dp-0.02425_dp) then
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

  real(dp) function gammap(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer :: n
    real(dp) :: ap, del, sumv, b, c, d, h, an
    real(dp), parameter :: eps=1.0e-14_dp, fpmin=1.0e-300_dp
    if (x <= 0.0_dp) then
    p=0.0_dp
    return
    end if
    if (x < a+1.0_dp) then
      ap=a
      sumv=1.0_dp/a
      del=sumv
      do n=1,10000
        ap=ap+1.0_dp
        del=del*x/ap
        sumv=sumv+del
        if (abs(del) < abs(sumv)*eps) exit
      end do
      p=sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a
      c=1.0_dp/fpmin
      d=1.0_dp/b
      h=d
      do n=1,10000
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp
        d=an*d+b
        if(abs(d)<fpmin)d=fpmin
        c=b+an/c
        if(abs(c)<fpmin)c=fpmin
        d=1.0_dp/d
        del=d*c
        h=h*del
        if(abs(del-1.0_dp)<eps) exit
      end do
      p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p=max(0.0_dp,min(1.0_dp,p))
  end function gammap

  real(dp) function chisq_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    if (x <= 0.0_dp) then
    p=0.0_dp
    else
    p=gammap(0.5_dp*df,0.5_dp*x)
    end if
  end function chisq_cdf

  real(dp) function chisq_quantile(p, df) result(x)
    real(dp), intent(in) :: p, df
    real(dp) :: lo, hi, mid
    integer :: k
    if(p<=0.0_dp)then
    x=0.0_dp
    return
    end if
    if(p>=1.0_dp)then
    x=huge(1.0_dp)
    return
    end if
    lo=0.0_dp
    hi=max(df+10.0_dp*sqrt(2.0_dp*df),1.0_dp)
    do while(chisq_cdf(hi,df)<p)
    hi=2.0_dp*hi
    end do
    do k=1,120
      mid=0.5_dp*(lo+hi)
      if(chisq_cdf(mid,df)<p)then
      lo=mid
      else
      hi=mid
      end if
    end do
    x=0.5_dp*(lo+hi)
  end function chisq_quantile

  real(dp) function median_real(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    integer :: n
    y=x
    call sort_inplace(y)
    n=size(y)
    if(mod(n,2)==1)then
    m=y((n+1)/2)
    else
    m=0.5_dp*(y(n/2)+y(n/2+1))
    end if
  end function median_real

  subroutine sort_inplace(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp)::v
    do i=2,size(x)
    v=x(i)
    j=i-1
    do while(j>=1)
    if(x(j)<=v)exit
    x(j+1)=x(j)
    j=j-1
    end do
    x(j+1)=v
    end do
  end subroutine sort_inplace

  subroutine invert_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:)
    real(dp) :: piv, fac
    integer :: n,i,k,imax
    n=size(a,1)
    allocate(aug(n,2*n),ainv(n,n))
    aug=0.0_dp
    aug(:,1:n)=a
    do i=1,n
    aug(i,n+i)=1.0_dp
    end do
    ok=.true.
    do k=1,n
      imax=k
      do i=k+1,n
      if(abs(aug(i,k))>abs(aug(imax,k)))imax=i
      end do
      if(abs(aug(imax,k))<1.0e-12_dp)then
      ok=.false.
      ainv=0.0_dp
      return
      end if
      if(imax/=k)call swap_rows(aug,k,imax)
      piv=aug(k,k)
      aug(k,:)=aug(k,:)/piv
      do i=1,n
        if(i==k)cycle
        fac=aug(i,k)
        aug(i,:)=aug(i,:)-fac*aug(k,:)
      end do
    end do
    ainv=aug(:,n+1:2*n)
  end subroutine invert_matrix

  subroutine swap_rows(a,i,j)
    real(dp),intent(inout)::a(:,:)
    integer,intent(in)::i,j
    real(dp)::tmp(size(a,2))
    tmp=a(i,:)
    a(i,:)=a(j,:)
    a(j,:)=tmp
  end subroutine swap_rows

  subroutine jacobi_sym(a, eval, evec)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::eval(:),evec(:,:)
    real(dp),allocatable::b(:,:)
    real(dp)::app,aqq,apq,phi,c,s,bip,biq,maxoff
    integer::n,i,j,p,q,it
    n=size(a,1)
    allocate(b(n,n),eval(n),evec(n,n))
    b=0.5_dp*(a+transpose(a))
    evec=0.0_dp
    do i=1,n
    evec(i,i)=1.0_dp
    end do
    do it=1,100*n*n
      maxoff=0.0_dp
      p=1
      q=min(2,n)
      do i=1,n-1
      do j=i+1,n
      if(abs(b(i,j))>maxoff)then
      maxoff=abs(b(i,j))
      p=i
      q=j
      end if
      end do
      end do
      if(maxoff<1.0e-12_dp)exit
      app=b(p,p)
      aqq=b(q,q)
      apq=b(p,q)
      phi=0.5_dp*atan2(2.0_dp*apq,aqq-app)
      c=cos(phi)
      s=sin(phi)
      do i=1,n
        if(i/=p .and. i/=q)then
        bip=b(i,p)
        biq=b(i,q)
        b(i,p)=c*bip-s*biq
        b(p,i)=b(i,p)
        b(i,q)=s*bip+c*biq
        b(q,i)=b(i,q)
        end if
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
      b(p,q)=0
      b(q,p)=0
      do i=1,n
      bip=evec(i,p)
      biq=evec(i,q)
      evec(i,p)=c*bip-s*biq
      evec(i,q)=s*bip+c*biq
      end do
    end do
    do i=1,n
    eval(i)=b(i,i)
    end do
  end subroutine jacobi_sym

  subroutine pseudoinverse_sym(a, pinv)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::pinv(:,:)
    real(dp),allocatable::v(:,:),e(:)
    real(dp)::tol
    integer::i
    call jacobi_sym(a,e,v)
    tol=max(1.0_dp,maxval(abs(e)))*1.0e-10_dp
    allocate(pinv(size(a,1),size(a,2)))
    pinv=0
    do i=1,size(e)
    if(abs(e(i))>tol)pinv=pinv+matmul(reshape(v(:,i),[size(a,1),1]),reshape(v(:,i),[1,size(a,1)]))/e(i)
    end do
  end subroutine pseudoinverse_sym

  subroutine symmetric_eigenvalues(a,e)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::e(:)
    real(dp),allocatable::v(:,:)
    call jacobi_sym(a,e,v)
  end subroutine symmetric_eigenvalues

  subroutine covariance_matrix(x,cov)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::cov(:,:)
    real(dp),allocatable::mu(:)
    integer::n,i
    n=size(x,1)
    allocate(mu(size(x,2)),cov(size(x,2),size(x,2)))
    mu=sum(x,dim=1)/real(n,dp)
    cov=0
    do i=1,n
    cov=cov+outer(x(i,:)-mu,x(i,:)-mu)
    end do
    if(n>1)cov=cov/real(n-1,dp)
  end subroutine covariance_matrix

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i
    do i=1,size(a)
    c(i,:)=a(i)*b
    end do
  end function outer

  subroutine numerical_gradient(fun,x,g)
    procedure(objective_fn)::fun
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    real(dp)::xp(size(x)),xm(size(x)),h
    integer::i
    do i=1,size(x)
    h=1.0e-5_dp*max(1.0_dp,abs(x(i)))
    xp=x
    xm=x
    xp(i)=xp(i)+h
    xm(i)=xm(i)-h
    g(i)=(fun(xp)-fun(xm))/(2*h)
    end do
  end subroutine numerical_gradient

  subroutine bfgs_minimize(fun,x,fval,converged,maxiter,tol)
    procedure(objective_fn)::fun
    real(dp),intent(inout)::x(:)
    real(dp),intent(out)::fval
    logical,intent(out)::converged
    integer,intent(in),optional::maxiter
    real(dp),intent(in),optional::tol
    real(dp),allocatable::hessinv(:,:),g(:),gnew(:),p(:),s(:),y(:),xn(:)
    real(dp)::f,fn,alpha,ys,rho,tolerance
    integer::n,it,mit,i
    n=size(x)
    mit=300
    if(present(maxiter))mit=maxiter
    tolerance=1e-7_dp
    if(present(tol))tolerance=tol
    allocate(hessinv(n,n),g(n),gnew(n),p(n),s(n),y(n),xn(n))
    hessinv=0
    do i=1,n
    hessinv(i,i)=1
    end do
    f=fun(x)
    call numerical_gradient(fun,x,g)
    converged=.false.
    do it=1,mit
      if(maxval(abs(g))<tolerance)then
      converged=.true.
      exit
      end if
      p=-matmul(hessinv,g)
      if(dot_product(p,g)>=0)p=-g
      alpha=1
      do
        xn=x+alpha*p
        fn=fun(xn)
        if(fn<=f+1.0e-4_dp*alpha*dot_product(g,p) .or. alpha<1e-8_dp)exit
        alpha=0.5_dp*alpha
      end do
      s=xn-x
      x=xn
      call numerical_gradient(fun,x,gnew)
      y=gnew-g
      ys=dot_product(y,s)
      if(ys>1e-12_dp)then
        rho=1/ys
        hessinv=matmul(mat_ident(n)-rho*outer(s,y),matmul(hessinv,mat_ident(n)-rho*outer(y,s)))+rho*outer(s,s)
      else
        hessinv=0
        do i=1,n
        hessinv(i,i)=1
        end do
      end if
      g=gnew
      f=fn
      if(maxval(abs(s))<tolerance*(1+maxval(abs(x))))then
      converged=.true.
      exit
      end if
    end do
    fval=f
  end subroutine bfgs_minimize

  pure function mat_ident(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0
    do i=1,n
    a(i,i)=1
    end do
  end function mat_ident

  real(dp) function scalar_second_derivative(fun,x,h) result(v)
    interface
    function fun(x) result(f)
    import dp
    real(dp),intent(in)::x
    real(dp)::f
    end function
    end interface
    real(dp),intent(in)::x
    real(dp),intent(in),optional::h
    real(dp)::hh
    hh=1.0e-3_dp*max(1.0_dp,abs(x))
    if(present(h))hh=h
    v=(fun(x+hh)-2*fun(x)+fun(x-hh))/(hh*hh)
  end function scalar_second_derivative

  real(dp) function randn() result(z)
    real(dp)::u1,u2
    call random_number(u1)
    call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    z=sqrt(-2*log(u1))*cos(2*acos(-1.0_dp)*u2)
  end function randn

  subroutine shuffle_int(x)
    integer,intent(inout)::x(:)
    integer::i,j,tmp
    real(dp)::u
    do i=size(x),2,-1
    call random_number(u)
    j=1+int(u*real(i,dp))
    if(j>i)j=i
    tmp=x(i)
    x(i)=x(j)
    x(j)=tmp
    end do
  end subroutine shuffle_int
end module compound_cox_math
