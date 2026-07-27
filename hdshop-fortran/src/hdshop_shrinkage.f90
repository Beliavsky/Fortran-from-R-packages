! SPDX-License-Identifier: GPL-3.0-only
! Derived from HDShOP 0.1.7.
module hdshop_shrinkage
  use hdshop_kinds, only: dp
  use hdshop_linalg, only: inverse_matrix, symmetric_eigen, trace_matrix, &
                           trace_product, quadratic_form, symmetrize
  use hdshop_stats, only: row_means, sample_covariance
  implicit none
  private

  type, public :: matrix_shrink_result
    real(dp), allocatable :: matrix(:,:)
    real(dp) :: alpha = 0.0_dp
    real(dp) :: beta = 0.0_dp
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type matrix_shrink_result

  type, public :: mean_shrink_result
    real(dp), allocatable :: means(:)
    real(dp) :: alpha = 0.0_dp
    real(dp) :: beta = 0.0_dp
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type mean_shrink_result

  public :: sigma_sample_estimator, cov_shrink_bgp14, inv_cov_shrink_bgp16
  public :: nonlin_shrink_lw, mean_bs, mean_js, mean_bop19
  public :: alpha_star_hat_bop19, beta_star_hat_bop19

contains

  function sigma_sample_estimator(x) result(s)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: s(:,:)
    s = sample_covariance(x)
  end function sigma_sample_estimator

  function cov_shrink_bgp14(n, target, sample_cov) result(res)
    integer, intent(in) :: n
    real(dp), intent(in) :: target(:,:), sample_cov(:,:)
    type(matrix_shrink_result) :: res
    real(dp) :: a1, a2, tt, st
    integer :: p
    p = size(sample_cov,1)
    allocate(res%matrix(p,p)); res%matrix = 0.0_dp
    if (n <= 0 .or. size(sample_cov,2)/=p .or. any(shape(target)/=[p,p])) then
      res%message='invalid dimensions'; return
    end if
    tt = trace_product(target,target)
    st = trace_product(sample_cov,target)
    a1 = trace_matrix(target)*0.0_dp ! initialized for compilers
    a1 = tt*trace_matrix(sample_cov)**2/real(n,dp)
    a2 = trace_product(sample_cov,sample_cov)*tt-st*st
    if (abs(a2) <= epsilon(1.0_dp)*max(1.0_dp,abs(a1))) then
      res%message='degenerate shrinkage denominator'; return
    end if
    res%alpha = 1.0_dp-a1/a2
    res%beta = st*(1.0_dp-res%alpha)/tt
    res%matrix = res%alpha*sample_cov+res%beta*target
    call symmetrize(res%matrix)
    res%ok=.true.
  end function cov_shrink_bgp14

  function inv_cov_shrink_bgp16(n, target, inverse_sample_cov) result(res)
    integer, intent(in) :: n
    real(dp), intent(in) :: target(:,:), inverse_sample_cov(:,:)
    type(matrix_shrink_result) :: res
    real(dp) :: a1, a2, tt, st, c
    integer :: p
    p=size(inverse_sample_cov,1)
    allocate(res%matrix(p,p)); res%matrix=0.0_dp
    if (n <= 0 .or. size(inverse_sample_cov,2)/=p .or. any(shape(target)/=[p,p])) then
      res%message='invalid dimensions'; return
    end if
    c=real(p,dp)/real(n,dp)
    tt=trace_product(target,target)
    st=trace_product(inverse_sample_cov,target)
    a1=tt*trace_matrix(inverse_sample_cov)**2/real(n,dp)
    a2=trace_product(inverse_sample_cov,inverse_sample_cov)*tt-st*st
    if (abs(a2)<=epsilon(1.0_dp)*max(1.0_dp,abs(a1))) then
      res%message='degenerate shrinkage denominator'; return
    end if
    res%alpha=1.0_dp-c-a1/a2
    res%beta=st*(1.0_dp-c-res%alpha)/tt
    res%matrix=res%alpha*inverse_sample_cov+res%beta*target
    call symmetrize(res%matrix)
    res%ok=.true.
  end function inv_cov_shrink_bgp16

  function nonlin_shrink_lw(x) result(shrunk)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: shrunk(:,:)
    real(dp), allocatable :: s(:,:), lambda(:), u(:,:), lam(:), ftilde(:), &
      hftilde(:), dtilde(:), lmat(:,:), hmat(:,:), xx(:,:), hftemp(:,:)
    real(dp) :: h, ratio, root5, pi, h0, d0, z, term
    logical :: ok
    integer :: p,n,m,i,j,start
    p=size(x,1); n=size(x,2); pi=acos(-1.0_dp); root5=sqrt(5.0_dp)
    s=sample_covariance(x)
    call symmetric_eigen(s,lambda,u,ok)
    allocate(shrunk(p,p)); shrunk=0.0_dp
    if (.not.ok .or. n<2) return
    h=real(n,dp)**(-1.0_dp/3.0_dp); ratio=real(p,dp)/real(n,dp)
    if (p<=n) then
      start=max(1,p-n+1); m=p-start+1
      allocate(lam(m)); lam=lambda(start:p)
      allocate(lmat(m,m),hmat(m,m),xx(m,m),hftemp(m,m),ftilde(m),hftilde(m),dtilde(p))
      do j=1,m; lmat(:,j)=lam; end do
      do i=1,m; hmat(i,:)=h*lam; end do
      xx=(lmat-transpose(lmat))/hmat
      do i=1,m
        ftilde(i)=0.0_dp; hftilde(i)=0.0_dp
        do j=1,m
          ftilde(i)=ftilde(i)+(3.0_dp/(4.0_dp*root5))*max(1.0_dp-xx(i,j)**2/5.0_dp,0.0_dp)/hmat(i,j)
          z=xx(i,j)
          if (abs(abs(z)-root5)<=100.0_dp*epsilon(1.0_dp)) then
            term=(-3.0_dp/(10.0_dp*pi))*z
          else
            term=(-3.0_dp/(10.0_dp*pi))*z+(3.0_dp/(4.0_dp*root5*pi))* &
              (1.0_dp-z*z/5.0_dp)*log(abs((root5-z)/(root5+z)))
          end if
          hftilde(i)=hftilde(i)+term/hmat(i,j)
        end do
        ftilde(i)=ftilde(i)/real(m,dp); hftilde(i)=hftilde(i)/real(m,dp)
      end do
      dtilde=0.0_dp
      do i=1,m
        dtilde(start+i-1)=lam(i)/((pi*ratio*lam(i)*ftilde(i))**2+ &
          (1.0_dp-ratio-pi*ratio*lam(i)*hftilde(i))**2)
      end do
    else
      start=max(1,p-n+2); m=p-start+1
      allocate(lam(m)); lam=lambda(start:p)
      allocate(lmat(m,m),hmat(m,m),xx(m,m),hftemp(m,m),ftilde(m),hftilde(m),dtilde(p))
      do j=1,m; lmat(:,j)=lam; end do
      do i=1,m; hmat(i,:)=h*lam; end do
      xx=(lmat-transpose(lmat))/hmat
      do i=1,m
        ftilde(i)=0.0_dp; hftilde(i)=0.0_dp
        do j=1,m
          ftilde(i)=ftilde(i)+(3.0_dp/(4.0_dp*root5))*max(1.0_dp-xx(i,j)**2/5.0_dp,0.0_dp)/hmat(i,j)
          z=xx(i,j)
          if (abs(abs(z)-root5)<=100.0_dp*epsilon(1.0_dp)) then
            term=(-3.0_dp/(10.0_dp*pi))*z
          else
            term=(-3.0_dp/(10.0_dp*pi))*z+(3.0_dp/(4.0_dp*root5*pi))* &
              (1.0_dp-z*z/5.0_dp)*log(abs((root5-z)/(root5+z)))
          end if
          hftilde(i)=hftilde(i)+term/hmat(i,j)
        end do
        ftilde(i)=ftilde(i)/real(m,dp); hftilde(i)=hftilde(i)/real(m,dp)
      end do
      h0=(1.0_dp/pi)*(3.0_dp/(10.0_dp*h*h)+3.0_dp/(4.0_dp*root5*h)* &
        (1.0_dp-1.0_dp/(5.0_dp*h*h))*log((1.0_dp+root5*h)/(1.0_dp-root5*h)))* &
        sum(1.0_dp/lam)/real(m,dp)
      d0=1.0_dp/(pi*real(p-n,dp)/real(n,dp)*h0)
      dtilde(1:p-n+1)=d0
      do i=1,m
        dtilde(p-n+1+i)=lam(i)/(pi*pi*lam(i)*lam(i)*(ftilde(i)**2+hftilde(i)**2))
      end do
    end if
    shrunk=matmul(u*spread(dtilde,1,p),transpose(u))
    call symmetrize(shrunk)
  end function nonlin_shrink_lw

  function mean_bs(x) result(res)
    real(dp), intent(in) :: x(:,:)
    type(mean_shrink_result) :: res
    real(dp), allocatable :: cov(:,:), inv(:,:), means(:), ones(:), diff(:)
    real(dp) :: mu0, denom
    logical :: ok
    integer :: p,n
    p=size(x,1); n=size(x,2); allocate(res%means(p)); res%means=0.0_dp
    if (n-p-2<=0) then; res%message='mean_bs requires n > p + 2'; return; end if
    cov=sample_covariance(x)*real(n-1,dp)/real(n-p-2,dp)
    call inverse_matrix(cov,inv,ok); if(.not.ok)then;res%message='singular covariance';return;endif
    means=row_means(x); allocate(ones(p)); ones=1.0_dp
    mu0=dot_product(ones,matmul(inv,means))/dot_product(ones,matmul(inv,ones))
    diff=means-mu0
    denom=real(p+2,dp)+real(n,dp)*quadratic_form(diff,inv)
    res%alpha=real(p+2,dp)/denom
    res%means=(1.0_dp-res%alpha)*means+res%alpha*mu0
    res%ok=.true.
  end function mean_bs

  function mean_js(x, y0) result(res)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: y0
    type(mean_shrink_result) :: res
    real(dp), allocatable :: cov(:,:), inv(:,:), means(:), diff(:)
    real(dp) :: target, denom
    logical :: ok
    integer :: p,n
    p=size(x,1); n=size(x,2); allocate(res%means(p)); res%means=0.0_dp
    target=1.0_dp; if(present(y0))target=y0
    cov=sample_covariance(x); call inverse_matrix(cov,inv,ok)
    if(.not.ok)then;res%message='singular covariance';return;endif
    means=row_means(x); diff=means-target
    denom=real(n,dp)*quadratic_form(diff,inv)
    if(denom<=0.0_dp)then;res%alpha=1.0_dp;else;res%alpha=min(1.0_dp,real(p-2,dp)/denom);endif
    res%means=(1.0_dp-res%alpha)*means+res%alpha*target
    res%ok=.true.
  end function mean_js

  real(dp) function alpha_star_hat_bop19(n,p,y,inv,mu0) result(alpha)
    integer,intent(in)::n,p
    real(dp),intent(in)::y(:),inv(:,:),mu0(:)
    real(dp)::i1,i2,i3,den
    i1=quadratic_form(y,inv)-real(p,dp)/real(n-p,dp)
    i2=quadratic_form(mu0,inv); i3=dot_product(y,matmul(inv,mu0))
    den=quadratic_form(y,inv)*i2-i3*i3
    alpha=(i1*i2-i3*i3)/den
  end function alpha_star_hat_bop19

  real(dp) function beta_star_hat_bop19(alpha,y,inv,mu0) result(beta)
    real(dp),intent(in)::alpha,y(:),inv(:,:),mu0(:)
    beta=(1.0_dp-alpha)*dot_product(y,matmul(inv,mu0))/quadratic_form(mu0,inv)
  end function beta_star_hat_bop19

  function mean_bop19(x,mu0) result(res)
    real(dp),intent(in)::x(:,:),mu0(:)
    type(mean_shrink_result)::res
    real(dp),allocatable::cov(:,:),inv(:,:),means(:)
    logical::ok
    integer::p,n
    p=size(x,1);n=size(x,2);allocate(res%means(p));res%means=0.0_dp
    if(size(mu0)/=p .or. n<=p)then;res%message='mean_bop19 requires matching target and n > p';return;endif
    cov=sample_covariance(x);call inverse_matrix(cov,inv,ok)
    if(.not.ok)then;res%message='singular covariance';return;endif
    means=row_means(x)
    res%alpha=alpha_star_hat_bop19(n,p,means,inv,mu0)
    res%beta=beta_star_hat_bop19(res%alpha,means,inv,mu0)
    res%means=res%alpha*means+res%beta*mu0
    res%ok=.true.
  end function mean_bop19

end module hdshop_shrinkage
