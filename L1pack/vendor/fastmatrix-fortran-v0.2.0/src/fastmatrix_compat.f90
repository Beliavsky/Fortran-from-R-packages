module fastmatrix_compat
  use fastmatrix_base, only: dp, gamma_p, normal_cdf, inverse_matrix, solve_linear
  implicit none
  private

  abstract interface
    function scalar_fun(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun
  end interface

  type, public :: mardia_result
    real(dp) :: skewness = 0.0_dp
    real(dp) :: kurtosis = 0.0_dp
    real(dp) :: skew_stat = 0.0_dp
    real(dp) :: kurt_stat = 0.0_dp
    real(dp) :: skew_pvalue = 1.0_dp
    real(dp) :: kurt_pvalue = 1.0_dp
    integer :: skew_df = 0
  end type mardia_result

  type, public :: harris_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: pvalue = 1.0_dp
    integer :: df = 0
    real(dp), allocatable :: covariance(:,:)
  end type harris_result

  type, public :: ols_compat_result
    real(dp), allocatable :: coefficients(:), fitted(:), residuals(:), cov_unscaled(:,:)
    real(dp) :: rss = 0.0_dp
    integer :: rank = 0
  end type ols_compat_result

  public :: schur_decomp, svd_decomp, matrix_fun, ols_fit_qr, ols_fit_svd
  public :: mardia_coefficients, mardia_test, harris_test

  interface
    subroutine dgees(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info)
      import dp
      character(len=1), intent(in) :: jobvs, sort
      external :: select
      logical :: select
      integer, intent(in) :: n, lda, ldvs, lwork
      integer, intent(out) :: sdim, info
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: wr(*), wi(*), vs(ldvs,*), work(*)
      logical, intent(out) :: bwork(*)
    end subroutine dgees
    subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobu, jobvt
      integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: s(*), u(ldu,*), vt(ldvt,*), work(*)
    end subroutine dgesvd
    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      import dp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*), b(ldb,*), work(*)
    end subroutine dgels
  end interface
contains

  logical function schur_select(wr, wi)
    real(dp), intent(in) :: wr, wi
    schur_select = (wr + wi) < -huge(1.0_dp)
  end function schur_select

  subroutine schur_decomp(a, t, wr, wi, q, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: t(:,:), wr(:), wi(:), q(:,:)
    integer, intent(out), optional :: info
    integer :: n, sdim, lwork, ier
    real(dp) :: workq(1)
    real(dp), allocatable :: work(:)
    logical, allocatable :: bwork(:)
    n = size(a,1)
    if (size(a,2) /= n .or. size(t,1) /= n .or. size(t,2) /= n) then
      if (present(info)) info = -1
      return
    end if
    t = a
    allocate(bwork(max(1,n)))
    lwork = -1
    call dgees('V','N',schur_select,n,t,n,sdim,wr,wi,q,n,workq,lwork,bwork,ier)
    if (ier /= 0) then
      if (present(info)) info = ier
      return
    end if
    lwork = max(1,int(workq(1)))
    allocate(work(lwork))
    t = a
    call dgees('V','N',schur_select,n,t,n,sdim,wr,wi,q,n,work,lwork,bwork,ier)
    if (present(info)) info = ier
  end subroutine schur_decomp

  subroutine svd_decomp(a, u, s, vt, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: u(:,:), s(:), vt(:,:)
    integer, intent(out), optional :: info
    integer :: m, n, k, lwork, ier
    real(dp) :: aq(size(a,1),size(a,2)), workq(1)
    real(dp), allocatable :: work(:)
    m=size(a,1); n=size(a,2); k=min(m,n)
    aq=a
    lwork=-1
    call dgesvd('S','S',m,n,aq,m,s,u,m,vt,k,workq,lwork,ier)
    if(ier/=0) then
      if(present(info)) info=ier
      return
    end if
    lwork=max(1,int(workq(1)))
    allocate(work(lwork))
    aq=a
    call dgesvd('S','S',m,n,aq,m,s,u,m,vt,k,work,lwork,ier)
    if(present(info)) info=ier
  end subroutine svd_decomp

  subroutine matrix_fun(a, fun, f, info)
    real(dp), intent(in) :: a(:,:)
    procedure(scalar_fun) :: fun
    real(dp), intent(out) :: f(:,:)
    integer, intent(out), optional :: info
    integer :: n, i, j, k, p
    real(dp) :: accum, den
    n=size(a,1)
    if(size(a,2)/=n) then
      if(present(info))info=-1
      return
    end if
    if(maxval(abs([( (a(i,j), i=j+1,n), j=1,n-1)])) > 100*epsilon(1.0_dp)) then
      if(present(info))info=-2
      return
    end if
    f=0.0_dp
    do i=1,n
      f(i,i)=fun(a(i,i))
    end do
    do p=1,n-1
      do i=1,n-p
        j=i+p
        accum=a(i,j)*(f(j,j)-f(i,i))
        do k=i+1,j-1
          accum=accum+a(i,k)*f(k,j)-f(i,k)*a(k,j)
        end do
        den=a(j,j)-a(i,i)
        if(abs(den)<=100*epsilon(1.0_dp)*max(1.0_dp,abs(a(i,i)),abs(a(j,j)))) then
          if(present(info))info=j
          return
        end if
        f(i,j)=accum/den
      end do
    end do
    if(present(info))info=0
  end subroutine matrix_fun

  subroutine ols_fit_qr(x,y,res,info)
    real(dp), intent(in) :: x(:,:), y(:)
    type(ols_compat_result), intent(out) :: res
    integer, intent(out), optional :: info
    integer :: m,n,ldb,lwork,ier,ier2
    real(dp), allocatable :: a(:,:),b(:,:),work(:),xtx(:,:)
    real(dp) :: workq(1)
    m=size(x,1); n=size(x,2); ldb=max(m,n)
    allocate(a(m,n),b(ldb,1),xtx(n,n))
    a=x; b=0.0_dp; b(1:m,1)=y
    lwork=-1
    call dgels('N',m,n,1,a,m,b,ldb,workq,lwork,ier)
    if(ier/=0) then
      if(present(info))info=ier
      return
    end if
    lwork=max(1,int(workq(1))); allocate(work(lwork)); a=x; b=0.0_dp; b(1:m,1)=y
    call dgels('N',m,n,1,a,m,b,ldb,work,lwork,ier)
    allocate(res%coefficients(n),res%fitted(m),res%residuals(m),res%cov_unscaled(n,n))
    res%coefficients=b(1:n,1); res%fitted=matmul(x,res%coefficients); res%residuals=y-res%fitted
    res%rss=sum(res%residuals**2); res%rank=n
    xtx=matmul(transpose(x),x); call inverse_matrix(xtx,res%cov_unscaled,ier2)
    if(present(info)) info=merge(ier,ier2,ier/=0)
  end subroutine ols_fit_qr

  subroutine ols_fit_svd(x,y,res,info,tol)
    real(dp), intent(in) :: x(:,:), y(:)
    type(ols_compat_result), intent(out) :: res
    integer, intent(out), optional :: info
    real(dp), intent(in), optional :: tol
    integer :: m,n,k,i,ier
    real(dp) :: eps
    real(dp), allocatable :: u(:,:),s(:),vt(:,:),tmp(:),v(:,:)
    m=size(x,1); n=size(x,2); k=min(m,n)
    allocate(u(m,k),s(k),vt(k,n),tmp(k),v(n,k))
    call svd_decomp(x,u,s,vt,ier)
    eps=max(m,n)*epsilon(1.0_dp)*maxval(s)
    if(present(tol))eps=tol
    tmp=matmul(transpose(u),y)
    do i=1,k
      if(s(i)>eps) then
        tmp(i)=tmp(i)/s(i)
      else
        tmp(i)=0.0_dp
      end if
    end do
    v=transpose(vt)
    allocate(res%coefficients(n),res%fitted(m),res%residuals(m),res%cov_unscaled(n,n))
    res%coefficients=matmul(v,tmp); res%fitted=matmul(x,res%coefficients); res%residuals=y-res%fitted
    res%rss=sum(res%residuals**2); res%rank=count(s>eps); res%cov_unscaled=0.0_dp
    do i=1,k
      if(s(i)>eps) res%cov_unscaled=res%cov_unscaled+spread(v(:,i),2,n)*spread(v(:,i),1,n)/(s(i)*s(i))
    end do
    if(present(info))info=ier
  end subroutine ols_fit_svd

  subroutine mardia_coefficients(x, skew, kurt, info)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: skew, kurt
    integer, intent(out), optional :: info
    integer :: n,p,i,j,ier
    real(dp) :: center(size(x,2)), cov(size(x,2),size(x,2)), ci(size(x,2),size(x,2))
    real(dp) :: zi(size(x,2)), zj(size(x,2)), d
    n=size(x,1); p=size(x,2)
    center=sum(x,dim=1)/real(n,dp); cov=0.0_dp
    do i=1,n
      zi=x(i,:)-center
      cov=cov+spread(zi,2,p)*spread(zi,1,p)
    end do
    cov=cov/real(n,dp)
    call inverse_matrix(cov,ci,ier)
    if(ier/=0) then
      skew=0.0_dp; kurt=0.0_dp; if(present(info))info=ier; return
    end if
    skew=0.0_dp; kurt=0.0_dp
    do i=1,n
      zi=x(i,:)-center; d=dot_product(zi,matmul(ci,zi)); kurt=kurt+d*d
      do j=1,n
        zj=x(j,:)-center; d=dot_product(zi,matmul(ci,zj)); skew=skew+d**3
      end do
    end do
    skew=skew/real(n*n,dp); kurt=kurt/real(n,dp)
    if(present(info))info=0
  end subroutine mardia_coefficients

  subroutine mardia_test(x,res,info)
    real(dp), intent(in) :: x(:,:)
    type(mardia_result), intent(out) :: res
    integer, intent(out), optional :: info
    integer :: n,p,ier
    real(dp) :: ex
    n=size(x,1); p=size(x,2)
    call mardia_coefficients(x,res%skewness,res%kurtosis,ier)
    ex=real(p*(p+2),dp); res%skew_stat=real(n,dp)*res%skewness/6.0_dp
    res%kurt_stat=(res%kurtosis-ex)/sqrt(8.0_dp*ex/real(n,dp))
    res%skew_df=p*(p+1)*(p+2)/6
    res%skew_pvalue=1.0_dp-gamma_p(0.5_dp*real(res%skew_df,dp),0.5_dp*res%skew_stat)
    res%kurt_pvalue=2.0_dp*(1.0_dp-normal_cdf(abs(res%kurt_stat)))
    if(present(info))info=ier
  end subroutine mardia_test

  subroutine harris_test(x,method,res,info)
    real(dp), intent(in) :: x(:,:)
    character(len=*), intent(in) :: method
    type(harris_result), intent(out) :: res
    integer, intent(out), optional :: info
    integer :: n,p,i,j,ier
    real(dp) :: center(size(x,2)), s(size(x,2)), h(size(x,2)-1), g(size(x,2)-1)
    real(dp) :: cov(size(x,2),size(x,2)), psi(size(x,2),size(x,2)), gg(size(x,2)-1,size(x,2)-1)
    real(dp) :: hh(size(x,2)-1,size(x,2)), z(size(x,2)), prodv
    n=size(x,1); p=size(x,2); center=sum(x,dim=1)/real(n,dp); cov=0.0_dp
    do i=1,n
      z=x(i,:)-center; cov=cov+spread(z,2,p)*spread(z,1,p)
    end do
    cov=cov/real(n-1,dp); allocate(res%covariance(p,p)); res%covariance=cov
    hh=0.0_dp
    do i=1,p-1
      hh(i,1)=-1.0_dp; hh(i,i+1)=1.0_dp
    end do
    s=[(cov(i,i),i=1,p)]
    select case(trim(method))
    case('Wald','wald')
      psi=cov*cov; h=matmul(hh,s); gg=matmul(matmul(hh,psi),transpose(hh))
      call solve_linear(gg,h,g,ier); res%statistic=0.5_dp*real(n,dp)*dot_product(g,h)
    case('log')
      psi=0.0_dp
      do i=1,p
        do j=1,p
          psi(i,j)=(cov(i,j)/sqrt(s(i)*s(j)))**2
        end do
      end do
      h=matmul(hh,log(s)); gg=matmul(matmul(hh,psi),transpose(hh))
      call solve_linear(gg,h,g,ier); res%statistic=0.5_dp*real(n,dp)*dot_product(g,h)
    case('robust')
      psi=0.0_dp
      do i=1,p
        do j=1,p
          psi(i,j)=sum((x(:,i)-center(i))**2*(x(:,j)-center(j))**2)/real(n,dp)
        end do
      end do
      h=matmul(hh,s); gg=matmul(matmul(hh,psi),transpose(hh)); call solve_linear(gg,h,g,ier)
      prodv=dot_product(g,h); res%statistic=real(n,dp)*prodv/(1.0_dp-prodv)
    case('log-robust')
      psi=0.0_dp
      do i=1,p
        do j=1,p
          psi(i,j)=sum((x(:,i)-center(i))**2*(x(:,j)-center(j))**2)/real(n,dp)/(s(i)*s(j))
        end do
      end do
      h=matmul(hh,log(s)); gg=matmul(matmul(hh,psi),transpose(hh)); call solve_linear(gg,h,g,ier)
      res%statistic=real(n,dp)*dot_product(g,h)
    case default
      if(present(info))info=-2; return
    end select
    res%df=p-1; res%pvalue=1.0_dp-gamma_p(0.5_dp*real(res%df,dp),0.5_dp*res%statistic)
    if(present(info))info=ier
  end subroutine harris_test
end module fastmatrix_compat
