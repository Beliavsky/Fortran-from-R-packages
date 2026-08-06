! SPDX-License-Identifier: GPL-3.0-only
module spantest_classical
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use spantest_kinds, only : dp
  use spantest_types, only : span_result, span_ok, span_invalid_input, &
    span_singular, span_insufficient_df
  use spantest_linalg, only : inverse_matrix, solve_vector, sample_covariance, &
    column_means, ols_fit_vector, ols_fit_matrix, quadratic_form
  use spantest_probability, only : f_upper_tail, normal_cdf, normal_quantile
  implicit none
  private
  public :: span_bj, span_f1, span_f2, span_grs, span_hk, span_km, span_py

contains

  pure real(dp) function qnan() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function qnan

  subroutine set_failure(res, h0, status, message)
    type(span_result), intent(out) :: res
    character(len=*), intent(in) :: h0, message
    integer, intent(in) :: status
    res%pval = qnan()
    res%stat = qnan()
    res%h0 = h0
    res%status = status
    res%message = message
  end subroutine set_failure

  logical function valid_pair(r1,r2)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    valid_pair = size(r1,1) == size(r2,1) .and. size(r1,1) > 0 .and. &
                 size(r1,2) > 0 .and. size(r2,2) > 0
  end function valid_pair

  function span_bj(r1,r2) result(res)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    type(span_result) :: res
    real(dp), allocatable :: r(:,:), one(:), b(:), e(:), b1(:), e1(:)
    real(dp) :: rssu, rssr
    integer :: n,k,m,p,info
    character(len=*), parameter :: h0 = 'alpha = 0'
    if (.not. valid_pair(r1,r2)) then
      call set_failure(res,h0,span_invalid_input,'incompatible or empty return matrices')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2); p=k+m
    if (n-k-m < 1) then
      call set_failure(res,h0,span_insufficient_df,'T-K-N must be positive')
      return
    end if
    allocate(r(n,p),one(n),b(p),e(n),b1(k),e1(n))
    r(:,1:k)=r1; r(:,k+1:p)=r2; one=1.0_dp
    call ols_fit_vector(r,one,b,e,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular full return cross-product')
      return
    end if
    call ols_fit_vector(r1,one,b1,e1,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular benchmark cross-product')
      return
    end if
    rssu=dot_product(e,e); rssr=dot_product(e1,e1)
    res%stat=((rssr-rssu)/real(m,dp))/(rssu/real(n-k-m,dp))
    res%pval=f_upper_tail(max(0.0_dp,res%stat),real(m,dp),real(n-k-m,dp))
    res%h0=h0; res%status=span_ok
  end function span_bj

  function span_grs(r1,r2) result(res)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    type(span_result) :: res
    real(dp), allocatable :: xx(:,:), beta(:,:), e(:,:), sig(:,:), sig_inv(:,:)
    real(dp), allocatable :: ssqm(:,:), ssqm_inv(:,:), mu(:), alpha(:)
    real(dp) :: denom, num
    integer :: n,k,m,info
    character(len=*), parameter :: h0 = 'alpha = 0'
    if (.not. valid_pair(r1,r2)) then
      call set_failure(res,h0,span_invalid_input,'incompatible or empty return matrices')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2)
    if (n-m-k < 1) then
      call set_failure(res,h0,span_insufficient_df,'T-N-K must be positive')
      return
    end if
    allocate(xx(n,k+1),beta(k+1,m),e(n,m),sig(m,m),sig_inv(m,m), &
             ssqm(k,k),ssqm_inv(k,k),mu(k),alpha(m))
    xx(:,1)=1.0_dp; xx(:,2:k+1)=r1
    call ols_fit_matrix(xx,r2,beta,e,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular benchmark design')
      return
    end if
    sig=matmul(transpose(e),e)/real(n,dp)
    call inverse_matrix(sig,sig_inv,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular residual covariance')
      return
    end if
    mu=column_means(r1)
    ssqm=sample_covariance(r1,.true.)
    call inverse_matrix(ssqm,ssqm_inv,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular benchmark covariance')
      return
    end if
    alpha=beta(1,:)
    denom=1.0_dp+quadratic_form(mu,ssqm_inv)
    num=real(n-m-k,dp)/real(m,dp)*quadratic_form(alpha,sig_inv)
    res%stat=num/denom
    res%pval=f_upper_tail(max(0.0_dp,res%stat),real(m,dp),real(n-m-k,dp))
    res%h0=h0; res%status=span_ok
  end function span_grs

  function span_f1(r1,r2) result(res)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    type(span_result) :: res
    real(dp), allocatable :: r(:,:), mu(:), mu1(:), s(:,:), s1(:,:), si(:,:), s1i(:,:)
    real(dp) :: a,a1
    integer :: n,k,m,p,df2,info
    character(len=*), parameter :: h0 = 'alpha = 0'
    if (.not. valid_pair(r1,r2)) then
      call set_failure(res,h0,span_invalid_input,'incompatible or empty return matrices')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2); p=k+m; df2=n-k-m
    if (df2 < 1) then
      call set_failure(res,h0,span_insufficient_df,'T-K-N must be positive')
      return
    end if
    allocate(r(n,p),mu(p),mu1(k),s(p,p),s1(k,k),si(p,p),s1i(k,k))
    r(:,1:k)=r1; r(:,k+1:p)=r2
    mu=column_means(r); mu1=column_means(r1)
    s=sample_covariance(r); s1=sample_covariance(r1)
    call inverse_matrix(s,si,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular full covariance')
      return
    end if
    call inverse_matrix(s1,s1i,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular benchmark covariance')
      return
    end if
    a=quadratic_form(mu,si); a1=quadratic_form(mu1,s1i)
    res%stat=real(df2,dp)/real(m,dp)*(a-a1)/(1.0_dp+a1)
    res%pval=f_upper_tail(max(0.0_dp,res%stat),real(m,dp),real(df2,dp))
    res%h0=h0; res%status=span_ok
  end function span_f1

  function span_f2(r1,r2) result(res)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    type(span_result) :: res
    real(dp), allocatable :: r(:,:), mu(:), mu1(:), s(:,:), s1(:,:), si(:,:), s1i(:,:), one(:), one1(:)
    real(dp) :: a,b,c,d,a1,b1,c1,d1
    integer :: n,k,m,p,df2,info
    character(len=*), parameter :: h0 = 'delta = 0'
    if (.not. valid_pair(r1,r2)) then
      call set_failure(res,h0,span_invalid_input,'incompatible or empty return matrices')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2); p=k+m; df2=n-k-m+1
    if (df2 < 1) then
      call set_failure(res,h0,span_insufficient_df,'T-K-N+1 must be positive')
      return
    end if
    allocate(r(n,p),mu(p),mu1(k),s(p,p),s1(k,k),si(p,p),s1i(k,k),one(p),one1(k))
    r(:,1:k)=r1; r(:,k+1:p)=r2; one=1.0_dp; one1=1.0_dp
    mu=column_means(r); mu1=column_means(r1)
    s=sample_covariance(r); s1=sample_covariance(r1)
    call inverse_matrix(s,si,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular full covariance')
      return
    end if
    call inverse_matrix(s1,s1i,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular benchmark covariance')
      return
    end if
    a=quadratic_form(mu,si); b=dot_product(mu,matmul(si,one)); c=quadratic_form(one,si); d=a*c-b*b
    a1=quadratic_form(mu1,s1i); b1=dot_product(mu1,matmul(s1i,one1)); c1=quadratic_form(one1,s1i); d1=a1*c1-b1*b1
    res%stat=real(df2,dp)/real(m,dp)*(((c+d)/(c1+d1))*((1.0_dp+a1)/(1.0_dp+a))-1.0_dp)
    res%pval=f_upper_tail(max(0.0_dp,res%stat),real(m,dp),real(df2,dp))
    res%h0=h0; res%status=span_ok
  end function span_f2

  function span_hk(r1,r2) result(res)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    type(span_result) :: res
    real(dp), allocatable :: r(:,:), mu(:), mu1(:), s(:,:), s1(:,:), si(:,:), s1i(:,:), one(:), one1(:)
    real(dp) :: a,b,c,d,a1,b1,c1,d1,u
    integer :: n,k,m,p,df2,info
    character(len=*), parameter :: h0 = 'alpha = 0 and delta = 0'
    if (.not. valid_pair(r1,r2)) then
      call set_failure(res,h0,span_invalid_input,'incompatible or empty return matrices')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2); p=k+m; df2=2*(n-k-m)
    if (df2 < 1) then
      call set_failure(res,h0,span_insufficient_df,'T-K-N must be positive')
      return
    end if
    allocate(r(n,p),mu(p),mu1(k),s(p,p),s1(k,k),si(p,p),s1i(k,k),one(p),one1(k))
    r(:,1:k)=r1; r(:,k+1:p)=r2; one=1.0_dp; one1=1.0_dp
    mu=column_means(r); mu1=column_means(r1)
    s=sample_covariance(r); s1=sample_covariance(r1)
    call inverse_matrix(s,si,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular full covariance')
      return
    end if
    call inverse_matrix(s1,s1i,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular benchmark covariance')
      return
    end if
    a=quadratic_form(mu,si); b=dot_product(mu,matmul(si,one)); c=quadratic_form(one,si); d=a*c-b*b
    a1=quadratic_form(mu1,s1i); b1=dot_product(mu1,matmul(s1i,one1)); c1=quadratic_form(one1,s1i); d1=a1*c1-b1*b1
    u=(c1+d1)/(c+d)
    if (u <= 0.0_dp) then
      call set_failure(res,h0,span_singular,'nonpositive frontier ratio')
      return
    end if
    res%stat=real(n-k-m,dp)/real(m,dp)*(1.0_dp/sqrt(u)-1.0_dp)
    res%pval=f_upper_tail(max(0.0_dp,res%stat),real(2*m,dp),real(df2,dp))
    res%h0=h0; res%status=span_ok
  end function span_hk

  function span_km(r1,r2) result(res)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    type(span_result) :: res
    real(dp), allocatable :: r(:,:), x(:,:), y(:), beta(:), e(:), xtx(:,:), xtxi(:,:), cmat(:,:), middle(:,:), theta(:)
    real(dp) :: sigma2
    integer :: n,k,m,p,j,info
    character(len=*), parameter :: h0 = 'delta = 0'
    if (.not. valid_pair(r1,r2)) then
      call set_failure(res,h0,span_invalid_input,'incompatible or empty return matrices')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2); p=k+m
    if (n-k-m < 1) then
      call set_failure(res,h0,span_insufficient_df,'T-K-N must be positive')
      return
    end if
    allocate(r(n,p),x(n,p),y(n),beta(p),e(n),xtx(p,p),xtxi(p,p),cmat(m,p),middle(m,m),theta(m))
    r(:,1:k)=r1; r(:,k+1:p)=r2; y=r(:,1); x(:,1)=1.0_dp
    do j=2,p
      x(:,j)=r(:,1)-r(:,j)
    end do
    call ols_fit_vector(x,y,beta,e,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular Kempf-Memmel design')
      return
    end if
    xtx=matmul(transpose(x),x)
    call inverse_matrix(xtx,xtxi,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular Kempf-Memmel cross-product')
      return
    end if
    sigma2=dot_product(e,e)/real(n-p,dp)
    cmat=0.0_dp
    do j=1,m
      cmat(j,k+j)=1.0_dp
    end do
    call inverse_matrix(matmul(matmul(cmat,xtxi),transpose(cmat)),middle,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular restriction covariance')
      return
    end if
    theta=matmul(cmat,beta)
    res%stat=quadratic_form(theta,middle)/(real(m,dp)*sigma2)
    res%pval=f_upper_tail(max(0.0_dp,res%stat),real(m,dp),real(n-p,dp))
    res%h0=h0; res%status=span_ok
  end function span_km

  function span_py(r1,r2) result(res)
    real(dp), intent(in) :: r1(:,:), r2(:,:)
    type(span_result) :: res
    real(dp), allocatable :: xx(:,:), beta(:,:), e(:,:), sigma(:,:), xtx(:,:), coef(:), mxone(:), one(:), t2(:)
    real(dp) :: v, num_scalar, pn, theta_n, rhobar, rho, jalpha2, den
    integer :: n,k,m,i,j,info
    character(len=*), parameter :: h0 = 'alpha = 0'
    if (.not. valid_pair(r1,r2)) then
      call set_failure(res,h0,span_invalid_input,'incompatible or empty return matrices')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2)
    if (m < 2 .or. n-k-m < 1 .or. n-k-1 <= 4) then
      call set_failure(res,h0,span_insufficient_df,'PY requires N>=2, T-K-N>=1, and T-K-1>4')
      return
    end if
    allocate(xx(n,k+1),beta(k+1,m),e(n,m),sigma(m,m),xtx(k,k),coef(k),mxone(n),one(n),t2(m))
    xx(:,1)=1.0_dp; xx(:,2:k+1)=r1; one=1.0_dp
    call ols_fit_matrix(xx,r2,beta,e,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular benchmark design')
      return
    end if
    sigma=matmul(transpose(e),e)/real(n,dp)
    xtx=matmul(transpose(r1),r1)
    call solve_vector(xtx,matmul(transpose(r1),one),coef,info)
    if (info /= 0) then
      call set_failure(res,h0,span_singular,'singular raw benchmark cross-product')
      return
    end if
    mxone=one-matmul(r1,coef)
    v=real(n-k-1,dp)
    num_scalar=sum(mxone)*v
    do i=1,m
      if (sigma(i,i) <= 0.0_dp) then
        call set_failure(res,h0,span_singular,'nonpositive residual variance')
        return
      end if
      t2(i)=beta(1,i)**2*num_scalar/(real(n,dp)*sigma(i,i))
    end do
    pn=0.05_dp/real(m-1,dp)
    theta_n=normal_quantile(1.0_dp-pn/2.0_dp)**2
    rhobar=0.0_dp
    do i=2,m
      do j=1,i-1
        rho=sigma(i,j)/sqrt(sigma(i,i)*sigma(j,j))
        if (v*rho*rho >= theta_n) rhobar=rhobar+rho*rho
      end do
    end do
    rhobar=2.0_dp*rhobar/real(m*(m-1),dp)
    jalpha2=sum(t2-v/(v-2.0_dp))/sqrt(real(m,dp))
    den=(v/(v-2.0_dp))*sqrt(2.0_dp*(v-1.0_dp)*(1.0_dp+real(m-1,dp)*rhobar)/(v-4.0_dp))
    res%stat=jalpha2/den
    res%pval=1.0_dp-normal_cdf(res%stat)
    res%h0=h0; res%status=span_ok
  end function span_py

end module spantest_classical
