! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_distributions
  use compositions_kinds, only: dp, pi
  use compositions_geometry, only: closure, clr, clr_inv, ilr, ilr_inv, clrvar_to_ilr
  use compositions_linalg, only: invert_matrix, determinant_spd, chol_lower, symmetric_eigen
  use bayesm_rng, only: rand_gamma, randn, rand_uniform
  implicit none
  private
  public :: dirichlet_logpdf, dirichlet_pdf, rdirichlet
  public :: logistic_normal_logpdf, logistic_normal_pdf, rlogistic_normal
  public :: aitchison_integrals_result, aitchison_integrals, aitchison_logpdf, aitchison_pdf, raitchison

  type :: aitchison_integrals_result
    real(dp) :: exp_kappa = 0.0_dp
    real(dp) :: loggx_mean = 0.0_dp
    real(dp), allocatable :: clr_mean(:)
    real(dp), allocatable :: clr_second(:,:)
    real(dp), allocatable :: clr_cov(:,:)
  end type aitchison_integrals_result
contains
  real(dp) function dirichlet_logpdf(x,alpha,aitchison_measure) result(ld)
    real(dp), intent(in) :: x(:),alpha(:)
    logical, intent(in), optional :: aitchison_measure
    integer :: i,k
    logical :: am
    if(size(x)/=size(alpha)) error stop 'dirichlet_logpdf: dimension mismatch'
    if(any(x<0.0_dp) .or. abs(sum(x)-1.0_dp)>1.0e-10_dp .or. any(alpha<=0.0_dp)) then
      ld=-huge(1.0_dp); return
    end if
    am=.false.; if(present(aitchison_measure)) am=aitchison_measure
    k=1; if(am) k=0
    ld=log_gamma(sum(alpha))-sum(log_gamma(alpha))
    do i=1,size(x)
      if(x(i)==0.0_dp) then
        if(alpha(i)-real(k,dp)>0.0_dp) then
          ld=-huge(1.0_dp); return
        else if(alpha(i)-real(k,dp)<0.0_dp) then
          ld=huge(1.0_dp); return
        end if
      else
        ld=ld+(alpha(i)-real(k,dp))*log(x(i))
      end if
    end do
  end function dirichlet_logpdf

  real(dp) function dirichlet_pdf(x,alpha,aitchison_measure) result(d)
    real(dp), intent(in) :: x(:),alpha(:)
    logical, intent(in), optional :: aitchison_measure
    real(dp) :: ld
    if(present(aitchison_measure)) then
      ld=dirichlet_logpdf(x,alpha,aitchison_measure)
    else
      ld=dirichlet_logpdf(x,alpha)
    end if
    if(ld <= -0.5_dp*huge(1.0_dp)) then; d=0.0_dp
    else; d=exp(ld); end if
  end function dirichlet_pdf

  function rdirichlet(n,alpha) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha(:)
    real(dp) :: x(n,size(alpha))
    real(dp) :: row(size(alpha)),s
    integer :: i,j
    if(any(alpha<=0.0_dp)) error stop 'rdirichlet: alpha must be positive'
    do i=1,n
      do j=1,size(alpha); row(j)=rand_gamma(alpha(j),1.0_dp); end do
      s=sum(row); x(i,:)=row/s
    end do
  end function rdirichlet

  real(dp) function logistic_normal_logpdf(x,mean_comp,clr_cov,with_jacobian,source_normalization) result(ld)
    real(dp), intent(in) :: x(:),mean_comp(:),clr_cov(:,:)
    logical, intent(in), optional :: with_jacobian,source_normalization
    real(dp), allocatable :: w(:),s(:,:),sinv(:,:)
    real(dp) :: q,detv,logc
    logical :: jac,src
    if(size(x)/=size(mean_comp)) error stop 'logistic_normal_logpdf: dimension mismatch'
    w=ilr(closure(x))-ilr(closure(mean_comp)); s=clrvar_to_ilr(clr_cov)
    call invert_matrix(s,sinv); q=dot_product(w,matmul(sinv,w)); detv=determinant_spd(s)
    jac=.false.; if(present(with_jacobian)) jac=with_jacobian
    src=.true.; if(present(source_normalization)) src=source_normalization
    if(src) then
      ! Source dnorm.acomp uses sqrt(2*pi*det) regardless of dimension.
      logc=0.5_dp*log(2.0_dp*pi*detv)
    else
      logc=0.5_dp*real(size(w),dp)*log(2.0_dp*pi)+0.5_dp*log(detv)
    end if
    ld=-0.5_dp*q-logc
    if(jac) ld=ld-sum(log(x))
  end function logistic_normal_logpdf

  real(dp) function logistic_normal_pdf(x,mean_comp,clr_cov,with_jacobian,source_normalization) result(d)
    real(dp), intent(in) :: x(:),mean_comp(:),clr_cov(:,:)
    logical, intent(in), optional :: with_jacobian,source_normalization
    if(present(with_jacobian).and.present(source_normalization)) then
      d=exp(logistic_normal_logpdf(x,mean_comp,clr_cov,with_jacobian,source_normalization))
    else if(present(with_jacobian)) then
      d=exp(logistic_normal_logpdf(x,mean_comp,clr_cov,with_jacobian=with_jacobian))
    else if(present(source_normalization)) then
      d=exp(logistic_normal_logpdf(x,mean_comp,clr_cov,source_normalization=source_normalization))
    else
      d=exp(logistic_normal_logpdf(x,mean_comp,clr_cov))
    end if
  end function logistic_normal_pdf

  function rlogistic_normal(n,mean_comp,clr_cov) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: mean_comp(:),clr_cov(:,:)
    real(dp), allocatable :: x(:,:),s(:,:),l(:,:),z(:),mu(:)
    integer :: i,j,k
    mu=ilr(closure(mean_comp)); s=clrvar_to_ilr(clr_cov); call chol_lower(s,l)
    allocate(x(n,size(mean_comp)),z(size(mu)))
    do i=1,n
      do j=1,size(z); z(j)=randn(); end do
      z=mu+matmul(l,z); x(i,:)=ilr_inv(z)
    end do
  end function rlogistic_normal

  function aitchison_integrals(theta,beta,grid,mode) result(res)
    real(dp), intent(in) :: theta(:),beta(:,:)
    integer, intent(in), optional :: grid,mode
    type(aitchison_integrals_result) :: res
    integer :: d,g,md,grid_plus,p,i,j,carry,grid_size
    integer, allocatable :: gc(:)
    real(dp), allocatable :: gv(:),cm(:,:),sq(:,:)
    real(dp) :: logmean,phi,dens,oneint,kappaint
    logical :: stepped
    d=size(theta); if(size(beta,1)/=d.or.size(beta,2)/=d) error stop 'aitchison_integrals: dimension mismatch'
    if(maxval(abs(beta-transpose(beta)))>1.0e-6_dp) error stop 'aitchison_integrals: beta must be symmetric'
    if(maxval(abs(sum(beta,dim=2)))>1.0e-8_dp) error stop 'aitchison_integrals: beta must be a clr matrix'
    g=30; if(present(grid)) g=grid; md=3; if(present(mode)) md=mode
    allocate(res%clr_mean(d),res%clr_second(d,d),res%clr_cov(d,d)); res%clr_mean=0.0_dp
    res%clr_second=0.0_dp; res%clr_cov=0.0_dp
    if(md<0) return
    allocate(gc(d),gv(d)); gc=0; gc(1)=g; grid_plus=g+d
    do p=1,d; gv(p)=log(real(gc(p)+1,dp)/real(grid_plus,dp)); end do
    oneint=0.0_dp; kappaint=0.0_dp; grid_size=0
    do
      stepped=.false.
      do p=1,d-1
        if(gc(p)>0) then
          carry=gc(p)-1; gc(p+1)=gc(p+1)+1; gc(p)=0; gc(1)=carry
          gv(p+1)=log(real(gc(p+1)+1,dp)/real(grid_plus,dp))
          gv(p)=log(real(gc(p)+1,dp)/real(grid_plus,dp))
          gv(1)=log(real(gc(1)+1,dp)/real(grid_plus,dp)); stepped=.true.; exit
        end if
      end do
      if(.not.stepped) exit
      logmean=sum(gv)/real(d,dp); phi=dot_product(gv,theta-1.0_dp)+dot_product(gv,matmul(beta,gv))
      dens=exp(phi); oneint=oneint+dens; kappaint=kappaint+dens*logmean; grid_size=grid_size+1
      if(md>=1) res%clr_mean=res%clr_mean+dens*(gv-logmean)
      if(md>=2) then
        do j=1,d; do i=1,d
          res%clr_second(i,j)=res%clr_second(i,j)+dens*(gv(i)-logmean)*(gv(j)-logmean)
        end do; end do
      end if
    end do
    if(oneint<=0.0_dp.or.grid_size==0) error stop 'aitchison_integrals: empty numerical integral'
    if(md>=1) res%clr_mean=res%clr_mean/oneint
    if(md>=2) res%clr_second=res%clr_second/oneint
    res%clr_cov=res%clr_second
    if(md>=3) then
      do j=1,d; do i=1,d; res%clr_cov(i,j)=res%clr_cov(i,j)-res%clr_mean(i)*res%clr_mean(j); end do; end do
    end if
    res%loggx_mean=kappaint/oneint; res%exp_kappa=oneint/real(grid_size,dp)
  end function aitchison_integrals

  real(dp) function aitchison_logpdf(x,theta,beta,exp_kappa,realdensity) result(ld)
    real(dp), intent(in) :: x(:),theta(:),beta(:,:),exp_kappa
    logical, intent(in), optional :: realdensity
    real(dp) :: c(size(x)),lx(size(x)),cf
    if(any(x<=0.0_dp)) then; ld=-huge(1.0_dp); return; end if
    c=clr(x); lx=log(closure(x)); cf=0.0_dp
    if(present(realdensity)) then; if(realdensity) cf=1.0_dp; end if
    ld=dot_product(c,matmul(beta,c))+dot_product(lx,theta-cf)-log(exp_kappa)
  end function aitchison_logpdf

  real(dp) function aitchison_pdf(x,theta,beta,exp_kappa,realdensity) result(dens)
    real(dp), intent(in) :: x(:),theta(:),beta(:,:),exp_kappa
    logical, intent(in), optional :: realdensity
    if(present(realdensity)) then
      dens=exp(aitchison_logpdf(x,theta,beta,exp_kappa,realdensity))
    else
      dens=exp(aitchison_logpdf(x,theta,beta,exp_kappa))
    end if
  end function aitchison_pdf

  function raitchison(n,theta,sigma) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: theta(:),sigma(:,:)
    real(dp), allocatable :: x(:,:),eig(:),vec(:,:),sqrt_sigma(:,:),z(:),c(:)
    real(dp) :: alpha_d,logd,maxdir,kappa,dir,u
    integer :: d,i,j,k,info
    if(any(theta<=0.0_dp)) error stop 'raitchison: source rejection sampler requires positive theta'
    d=size(theta); if(size(sigma,1)/=d.or.size(sigma,2)/=d) error stop 'raitchison: dimension mismatch'
    call symmetric_eigen(sigma,eig,vec,info); eig=max(eig,0.0_dp)
    allocate(sqrt_sigma(d,d)); sqrt_sigma=0.0_dp
    do i=1,d; sqrt_sigma=sqrt_sigma+sqrt(eig(i))*outer(vec(:,i),vec(:,i)); end do
    alpha_d=sum(theta); logd=sum(theta*(log(theta)-log(alpha_d))); maxdir=exp(logd)
    allocate(x(n,d),z(d),c(d)); k=0
    do while(k<n)
      z=0.0_dp
      do i=1,d
        u=randn(); z=z+sqrt_sigma(:,i)*u
      end do
      kappa=log(sum(exp(z))); logd=dot_product(theta,z)-alpha_d*kappa; dir=exp(logd)
      if(dir/maxdir>=rand_uniform()) then; k=k+1; x(k,:)=clr_inv(z); end if
    end do
  contains
    function outer(a,b) result(cmat)
      real(dp), intent(in) :: a(:),b(:)
      real(dp) :: cmat(size(a),size(b)); integer :: ii,jj
      do jj=1,size(b); do ii=1,size(a); cmat(ii,jj)=a(ii)*b(jj); end do; end do
    end function outer
  end function raitchison
end module compositions_distributions
