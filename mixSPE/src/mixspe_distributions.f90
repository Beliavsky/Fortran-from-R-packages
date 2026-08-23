module mixspe_distributions
  use mvtnorm_kinds, only : dp, pi
  use mvtnorm_distributions, only : rmvnorm, dmvnorm_one
  use mvtnorm_random, only : random_normals
  use mixspe_special, only : normal_cdf, normal_logcdf
  use mixspe_linalg, only : sym_power, mahalanobis2, determinant_sym, sym_eigen
  implicit none
  private
  public :: dpe, log_dpe, rpe, dspe, log_dspe, rspe, cov_pe
contains
  function cov_pe(scale,beta) result(cov)
    real(dp),intent(in)::scale(:,:),beta
    real(dp),allocatable::cov(:,:)
    integer::p
    real(dp)::c
    p=size(scale,1)
    c=2.0_dp**(1.0_dp/beta)*gamma((real(p,dp)+2.0_dp)/(2.0_dp*beta)) / &
      (real(p,dp)*gamma(real(p,dp)/(2.0_dp*beta)))
    cov=c*scale
  end function

  real(dp) function log_dpe(y,mu,sigma,beta) result(v)
    real(dp),intent(in)::y(:),mu(:),sigma(:,:),beta
    real(dp)::q,det
    integer::p
    logical::ok
    p=size(y); q=mahalanobis2(y,mu,sigma,ok); det=determinant_sym(sigma,ok)
    if(.not.ok .or. det<=0.0_dp .or. beta<=0.0_dp) then; v=-huge(1.0_dp); return; end if
    v=log(beta)+log_gamma(real(p,dp)/2.0_dp)-real(p,dp)/2.0_dp*log(pi) &
      -0.5_dp*q**beta-real(p,dp)/(2.0_dp*beta)*log(2.0_dp) &
      -log_gamma(real(p,dp)/(2.0_dp*beta))-0.5_dp*log(det)
  end function
  real(dp) function dpe(y,mu,sigma,beta) result(v)
    real(dp),intent(in)::y(:),mu(:),sigma(:,:),beta
    v=exp(log_dpe(y,mu,sigma,beta))
  end function

  real(dp) function log_dspe(y,mu,sigma,psi,beta) result(v)
    real(dp),intent(in)::y(:),mu(:),sigma(:,:),psi(:),beta
    real(dp),allocatable::imhalf(:,:)
    real(dp)::z
    logical::ok
    imhalf=sym_power(sigma,-0.5_dp,ok)
    if(.not.ok) then; v=-huge(1.0_dp); return; end if
    z=dot_product(psi,matmul(imhalf,y-mu))
    v=log(2.0_dp)+log_dpe(y,mu,sigma,beta)+normal_logcdf(z)
  end function
  real(dp) function dspe(y,mu,sigma,psi,beta) result(v)
    real(dp),intent(in)::y(:),mu(:),sigma(:,:),psi(:),beta
    v=exp(log_dspe(y,mu,sigma,psi,beta))
  end function

  function rpe(n,beta,mean,scale) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::beta,mean(:),scale(:,:)
    real(dp),allocatable::x(:,:),vals(:),vecs(:,:),root(:,:),z(:)
    real(dp)::radius,g
    integer::i,j,p
    logical::ok
    p=size(mean); allocate(x(n,p),z(p),root(p,p)); root=0.0_dp
    call sym_eigen(scale,vals,vecs,ok)
    if(.not.ok) error stop "rpe: scale eigendecomposition failed"
    do j=1,p
      root=root+sqrt(max(vals(j),0.0_dp))*spread(vecs(:,j),2,p)*spread(vecs(:,j),1,p)
    end do
    do i=1,n
      call random_normals(z); z=z/sqrt(sum(z*z))
      g=gamma_rng(real(p,dp)/(2.0_dp*beta),2.0_dp)
      radius=g**(1.0_dp/(2.0_dp*beta))
      x(i,:)=mean+radius*matmul(root,z)
    end do
  end function

  function rspe(n,location,scale,beta,psi,burnin) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::location(:),scale(:,:),beta,psi(:)
    integer,intent(in),optional::burnin
    real(dp),allocatable::x(:,:),cov(:,:),proposal_cov(:,:),cur(:),cand(:),tmp(:,:)
    real(dp)::logalpha,u,lc,ln,hc,hn
    integer::i,acc,b
    b=2000; if(present(burnin)) b=burnin
    cov=cov_pe(scale,beta); proposal_cov=2.0_dp*cov
    allocate(x(n,size(location)),cur(size(location)),cand(size(location)))
    cur=location-0.2_dp
    lc=log_dspe(cur,location,scale,psi,beta)
    hc=dmvnorm_one(cur,location,proposal_cov,.true.)
    acc=0; i=0
    do while(acc<n+b)
      tmp=rmvnorm(1,location,proposal_cov); cand=tmp(1,:)
      ln=log_dspe(cand,location,scale,psi,beta)
      hn=dmvnorm_one(cand,location,proposal_cov,.true.)
      logalpha=min(0.0_dp,ln+hc-lc-hn)
      call random_number(u)
      if(log(max(u,tiny(1.0_dp)))<logalpha) then
        cur=cand; lc=ln; hc=hn; acc=acc+1
        if(acc>b) then; i=i+1; x(i,:)=cur; end if
      end if
    end do
  end function

  recursive real(dp) function gamma_rng(shape,scale) result(g)
    real(dp),intent(in)::shape,scale
    real(dp)::d,c,z,u
    if(shape<=0.0_dp) then; g=0.0_dp; return; end if
    if(shape<1.0_dp) then
      call random_number(u)
      g=gamma_rng(shape+1.0_dp,scale)*u**(1.0_dp/shape); return
    end if
    d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
    do
      call random_normals_scalar(z)
      if(1.0_dp+c*z<=0.0_dp) cycle
      z=(1.0_dp+c*z)**3
      call random_number(u)
      if(u<1.0_dp-0.0331_dp*((z**(1.0_dp/3.0_dp)-1.0_dp)/c)**4) exit
      if(log(u)<0.5_dp*((z**(1.0_dp/3.0_dp)-1.0_dp)/c)**2+d*(1.0_dp-z+log(z))) exit
    end do
    g=scale*d*z
  end function
  subroutine random_normals_scalar(z)
    real(dp),intent(out)::z
    real(dp)::a(1)
    call random_normals(a); z=a(1)
  end subroutine
end module
