! SPDX-License-Identifier: GPL-3.0-only
module esback_esreg
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use esback_kinds, only: dp, huge_penalty, pi
  use esback_types
  use esback_math
  use esback_optimizer
  implicit none
  private
  public :: esreg_fit, esreg_covariance, esr_loss, quantile_regression
  public :: g1_value, g1_prime, g1_second, g2_curly, g2_value, g2_prime, g2_second
  public :: conditional_mean_sigma, conditional_truncated_variance
  public :: density_quantile_function, cdf_at_quantile

  type :: quantile_context
    real(dp), allocatable :: x(:,:), y(:)
    real(dp) :: tau
  end type quantile_context

  type :: esr_context
    real(dp), allocatable :: xq(:,:), xe(:,:), y(:)
    real(dp) :: alpha
  end type esr_context

  type :: locscale_context
    real(dp), allocatable :: x(:,:), y(:)
  end type locscale_context
contains
  pure elemental real(dp) function g1_value(z,g1) result(v)
    real(dp),intent(in)::z;integer,intent(in)::g1
    if(g1==1)then;v=z;else;v=0.0_dp;end if
  end function g1_value
  pure elemental real(dp) function g1_prime(z,g1) result(v)
    real(dp),intent(in)::z;integer,intent(in)::g1
    if (g1 == 1) then
      v = 1.0_dp + 0.0_dp*z
    else
      v = 0.0_dp*z
    end if
  end function g1_prime
  pure elemental real(dp) function g1_second(z,g1) result(v)
    real(dp),intent(in)::z;integer,intent(in)::g1
    v = 0.0_dp*z + 0.0_dp*real(g1,dp)
  end function g1_second
  pure elemental real(dp) function g2_curly(z,g2) result(v)
    real(dp),intent(in)::z;integer,intent(in)::g2
    select case(g2)
    case(1); if(z<0.0_dp)then;v=-log(-z);else;v=huge_penalty;end if
    case(2); if(z<0.0_dp)then;v=-sqrt(-z);else;v=huge_penalty;end if
    case(3); if(z<0.0_dp)then;v=-1.0_dp/z;else;v=huge_penalty;end if
    case(4); if(z>35.0_dp)then;v=z;else;v=log(1.0_dp+exp(z));end if
    case(5); v=exp(min(z,700.0_dp))
    case default;v=huge_penalty
    end select
  end function g2_curly
  pure elemental real(dp) function g2_value(z,g2) result(v)
    real(dp),intent(in)::z;integer,intent(in)::g2
    select case(g2)
    case(1);if(z<0.0_dp)then;v=-1.0_dp/z;else;v=huge_penalty;end if
    case(2);if(z<0.0_dp)then;v=0.5_dp/sqrt(-z);else;v=huge_penalty;end if
    case(3);if(z<0.0_dp)then;v=1.0_dp/(z*z);else;v=huge_penalty;end if
    case(4)
      if (z >= 0.0_dp) then
        v = 1.0_dp/(1.0_dp + exp(-min(z,700.0_dp)))
      else
        v = exp(max(z,-700.0_dp))/(1.0_dp + exp(max(z,-700.0_dp)))
      end if
    case(5);v=exp(min(z,700.0_dp))
    case default;v=huge_penalty
    end select
  end function g2_value
  pure elemental real(dp) function g2_prime(z,g2) result(v)
    real(dp),intent(in)::z;integer,intent(in)::g2
    real(dp)::ez
    select case(g2)
    case(1);if(z<0.0_dp)then;v=1.0_dp/(z*z);else;v=huge_penalty;end if
    case(2);if(z<0.0_dp)then;v=0.25_dp/((-z)**1.5_dp);else;v=huge_penalty;end if
    case(3);if(z<0.0_dp)then;v=-2.0_dp/(z**3);else;v=huge_penalty;end if
    case(4);ez=exp(max(-700.0_dp,min(700.0_dp,z)));v=ez/(1.0_dp+ez)**2
    case(5);v=exp(min(z,700.0_dp))
    case default;v=huge_penalty
    end select
  end function g2_prime
  pure elemental real(dp) function g2_second(z,g2) result(v)
    real(dp),intent(in)::z;integer,intent(in)::g2
    real(dp)::ez
    select case(g2)
    case(1)
      ! Preserve the executable esreg source, which uses -1/z^3 here.
      if(z<0.0_dp)then;v=-1.0_dp/(z**3);else;v=huge_penalty;end if
    case(2);if(z<0.0_dp)then;v=0.375_dp/((-z)**2.5_dp);else;v=huge_penalty;end if
    case(3);if(z<0.0_dp)then;v=6.0_dp/(z**4);else;v=huge_penalty;end if
    case(4);ez=exp(max(-700.0_dp,min(700.0_dp,z)));v=-(ez*(ez-1.0_dp))/(ez+1.0_dp)**3
    case(5);v=exp(min(z,700.0_dp))
    case default;v=huge_penalty
    end select
  end function g2_second

  pure real(dp) function esr_loss(r,q,e,alpha,g1,g2,return_mean) result(loss)
    real(dp),intent(in)::r(:),q(:),e(:),alpha
    integer,intent(in),optional::g1,g2
    logical,intent(in),optional::return_mean
    integer::i,ig1,ig2
    logical::avg,hit
    real(dp)::li
    ig1=2;if(present(g1))ig1=g1;ig2=1;if(present(g2))ig2=g2;avg=.true.;if(present(return_mean))avg=return_mean
    loss=0.0_dp
    do i=1,size(r)
      hit=r(i)<=q(i)
      li=(merge(1.0_dp,0.0_dp,hit)-alpha)*g1_value(q(i),ig1)-merge(g1_value(r(i),ig1),0.0_dp,hit)+ &
         g2_value(e(i),ig2)*(e(i)-q(i)+merge(q(i)-r(i),0.0_dp,hit)/alpha)-g2_curly(e(i),ig2)
      if(.not.ieee_is_finite(li) .or. abs(li)>=huge_penalty)then;loss=huge_penalty;return;end if
      loss=loss+li
    end do
    if(avg)loss=loss/real(size(r),dp)
  end function esr_loss

  real(dp) function quantile_objective(b,context) result(f)
    real(dp),intent(in)::b(:);class(*),intent(in)::context
    real(dp),allocatable::u(:)
    select type(c=>context)
    type is(quantile_context)
      u=c%y-matmul(c%x,b)
      f=sum(merge(c%tau*u,(c%tau-1.0_dp)*u,u>=0.0_dp))
    class default;f=huge_penalty
    end select
  end function quantile_objective

  subroutine quantile_regression(x,y,tau,b,status,maxit,tol)
    real(dp),intent(in)::x(:,:),y(:),tau
    real(dp),allocatable,intent(out)::b(:)
    integer,intent(out)::status
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    type(quantile_context)::ctx
    real(dp),allocatable::b0(:),step(:)
    real(dp)::f
    integer::st,it,mx
    mx=3000;if(present(maxit))mx=maxit
    allocate(ctx%x(size(x,1),size(x,2)), ctx%y(size(y)))
    ctx%x = x
    ctx%y = y
    ctx%tau = tau
    call least_squares(x,y,b0,st,1.0e-8_dp)
    if(st/=esback_ok)then;allocate(b0(size(x,2)));b0=0.0_dp;end if
    allocate(step(size(b0)));step=0.15_dp*max(abs(b0),1.0_dp)
    if(present(tol))then
      call nelder_mead(quantile_objective,ctx,b0,step,mx,tol,b,f,status,it)
    else
      call nelder_mead(quantile_objective,ctx,b0,step,mx,1.0e-9_dp,b,f,status,it)
    end if
  end subroutine quantile_regression

  real(dp) function esr_objective(b,context) result(f)
    real(dp),intent(in)::b(:);class(*),intent(in)::context
    real(dp),allocatable::q(:),e(:)
    integer::kq
    select type(c=>context)
    type is(esr_context)
      kq=size(c%xq,2);q=matmul(c%xq,b(:kq));e=matmul(c%xe,b(kq+1:))
      if(any(e>=-0.01_dp))then;f=huge_penalty;return;end if
      f=esr_loss(c%y,q,e,c%alpha,2,1,.true.)
    class default;f=huge_penalty
    end select
  end function esr_objective

  subroutine esreg_fit(xq,xe,y,alpha,fit,options,compute_covariance)
    real(dp),intent(in)::xq(:,:),xe(:,:),y(:),alpha
    type(esreg_fit_result),intent(out)::fit
    type(esreg_options),intent(in),optional::options
    logical,intent(in),optional::compute_covariance
    type(esreg_options)::opt
    type(esr_context)::ctx
    type(rng_state)::rng
    real(dp),allocatable::yt(:),bq(:),be(:),b0(:),bt(:),step(:),best(:),cand(:)
    real(dp)::maxy,alpha_tilde,z,esn,bestf,f,maxe
    integer::kq,ke,st,it,j,total_it
    logical::do_cov
    opt=esreg_options();if(present(options))opt=options;do_cov=.true.;if(present(compute_covariance))do_cov=compute_covariance
    if(size(y)<5 .or. size(xq,1)/=size(y) .or. size(xe,1)/=size(y) .or. alpha<=0.0_dp .or. alpha>=1.0_dp)then
      fit%status=esback_invalid_input;return
    end if
    kq=size(xq,2);ke=size(xe,2);allocate(yt(size(y)));maxy=maxval(y);yt=y-maxy
    z=normal_quantile(alpha);esn=-normal_pdf(z)/alpha;alpha_tilde=max(1.0e-7_dp,min(alpha-1.0e-7_dp,normal_cdf(esn)))
    call quantile_regression(xq,yt,alpha,bq,st,opt%max_iterations,opt%tolerance)
    if(st/=esback_ok)then;fit%status=st;return;end if
    call quantile_regression(xe,yt,alpha_tilde,be,st,opt%max_iterations,opt%tolerance)
    if(st/=esback_ok)then;fit%status=st;return;end if
    maxe=maxval(matmul(xe,be));if(maxe>=-0.1_dp)be(1)=be(1)-(maxe+0.1_dp)
    allocate(b0(kq+ke),step(kq+ke));b0=[bq,be];step=0.10_dp*max(abs(b0),0.25_dp)
    allocate(ctx%xq(size(xq,1),size(xq,2)))
    allocate(ctx%xe(size(xe,1),size(xe,2)))
    allocate(ctx%y(size(yt)))
    ctx%xq = xq
    ctx%xe = xe
    ctx%y = yt
    ctx%alpha = alpha
    call nelder_mead(esr_objective,ctx,b0,step,opt%max_iterations,opt%tolerance,best,bestf,st,it)
    if(st/=esback_ok)then;fit%status=st;return;end if
    total_it=it;call rng_seed(rng,opt%seed)
    do j=1,max(0,opt%multistarts-1)
      allocate(bt(size(b0)));bt=best
      bt=bt+step*[(rng_normal(rng),it=1,size(bt))]
      maxe=maxval(matmul(xe,bt(kq+1:)));if(maxe>=-0.01_dp)bt(kq+1)=bt(kq+1)-(maxe+0.1_dp)
      call nelder_mead(esr_objective,ctx,bt,step,opt%max_iterations,opt%tolerance,cand,f,st,it)
      total_it=total_it+it
      if(st==esback_ok .and. f<bestf)then;best=cand;bestf=f;end if
      deallocate(bt)
    end do
    best(1)=best(1)+maxy;best(kq+1)=best(kq+1)+maxy
    allocate(fit%coefficients(kq+ke),fit%coefficients_q(kq),fit%coefficients_e(ke),fit%fitted_q(size(y)),fit%fitted_e(size(y)))
    fit%coefficients=best;fit%coefficients_q=best(:kq);fit%coefficients_e=best(kq+1:)
    fit%fitted_q=matmul(xq,fit%coefficients_q);fit%fitted_e=matmul(xe,fit%coefficients_e)
    fit%loss=bestf;fit%iterations=total_it;fit%status=esback_ok
    if(do_cov)then
      call esreg_covariance(xq,xe,y,alpha,fit%coefficients_q,fit%coefficients_e,opt,fit%covariance,st)
      fit%covariance_available=st==esback_ok
    end if
  end subroutine esreg_fit

  real(dp) function locscale_objective(b,context) result(f)
    real(dp),intent(in)::b(:);class(*),intent(in)::context
    integer::k
    real(dp),allocatable::mu(:),sig(:),z(:)
    select type(c=>context)
    type is(locscale_context)
      k=size(c%x,2);mu=matmul(c%x,b(:k));sig=matmul(c%x,b(k+1:))
      if(any(sig<=1.0e-8_dp))then;f=huge_penalty;return;end if
      z=(c%y-mu)/sig;f=sum(log(sig)+0.5_dp*z*z+0.5_dp*log(2.0_dp*pi))
    class default;f=huge_penalty
    end select
  end function locscale_objective

  subroutine conditional_mean_sigma(y,x,mu,sigma,status)
    real(dp),intent(in)::y(:),x(:,:)
    real(dp),allocatable,intent(out)::mu(:),sigma(:)
    integer,intent(out)::status
    real(dp),allocatable::b1(:),b2(:),b0(:),step(:),b(:),res(:)
    real(dp)::f,minfit
    integer::st,it,k
    type(locscale_context)::ctx
    k=size(x,2);call least_squares(x,y,b1,st,1.0e-8_dp)
    if(st/=esback_ok)then;status=st;return;end if
    res=y-matmul(x,b1);call least_squares(x,abs(res),b2,st,1.0e-8_dp)
    if(st/=esback_ok)then;status=st;return;end if
    minfit=minval(matmul(x,b2));b2(1)=b2(1)-min(0.001_dp,minfit)
    if(minval(matmul(x,b2))<=1.0e-6_dp)b2(1)=b2(1)+(1.0e-3_dp-minval(matmul(x,b2)))
    allocate(b0(2*k),step(2*k));b0=[b1,b2];step=0.08_dp*max(abs(b0),0.1_dp)
    allocate(ctx%x(size(x,1),size(x,2)), ctx%y(size(y)))
    ctx%x = x
    ctx%y = y
    call nelder_mead(locscale_objective,ctx,b0,step,5000,1.0e-9_dp,b,f,st,it)
    if(st/=esback_ok)then;b=b0;end if
    mu=matmul(x,b(:k));sigma=matmul(x,b(k+1:))
    if(any(sigma<=0.0_dp))then;status=esback_optimization_failed;else;status=esback_ok;end if
  end subroutine conditional_mean_sigma

  subroutine density_quantile_function(y,x,u,alpha,sparsity,bandwidth_estimator,density,status)
    real(dp),intent(in)::y(:),x(:,:),u(:),alpha
    integer,intent(in)::sparsity,bandwidth_estimator
    real(dp),allocatable,intent(out)::density(:)
    integer,intent(out)::status
    real(dp)::bandwidth,z,tau,eps,slope
    real(dp),allocatable::bu(:),bl(:),den(:),sel(:),xt(:)
    integer,allocatable::idx(:)
    integer::n,k,h,m,i,j,tmp,st
    n=size(y);k=size(x,2);tau=0.05_dp;z=normal_quantile(alpha);eps=epsilon(1.0_dp)**(2.0_dp/3.0_dp)
    select case(bandwidth_estimator)
    case(bandwidth_bofinger)
      bandwidth=real(n,dp)**(-0.2_dp)*((4.5_dp*normal_pdf(z)**4)/(2.0_dp*z*z+1.0_dp)**2)**0.2_dp
    case(bandwidth_chamberlain)
      bandwidth=normal_quantile(1.0_dp-alpha/2.0_dp)*sqrt(tau*(1.0_dp-tau)/real(n,dp))
    case default
      bandwidth=real(n,dp)**(-1.0_dp/3.0_dp)*normal_quantile(1.0_dp-tau/2.0_dp)**(2.0_dp/3.0_dp)* &
        ((1.5_dp*normal_pdf(z)**2)/(2.0_dp*z*z+1.0_dp))**(1.0_dp/3.0_dp)
    end select
    allocate(density(n))
    if(sparsity==sparsity_nid)then
      call quantile_regression(x, y, min(0.999999_dp, alpha + bandwidth), &
        bu, st, 4000, 1.0e-8_dp)
      if (st /= esback_ok) then
        status = st
        return
      end if
      call quantile_regression(x, y, max(0.000001_dp, alpha - bandwidth), &
        bl, st, 4000, 1.0e-8_dp)
      if (st /= esback_ok) then
        status = st
        return
      end if
      den=matmul(x,bu-bl)-eps;density=max(0.0_dp,2.0_dp*bandwidth/den)
    else
      allocate(idx(n));idx=[(i,i=1,n)]
      do i=2,n;tmp=idx(i);j=i-1;do while(j>=1 .and. abs(u(idx(j)))>abs(u(tmp)));idx(j+1)=idx(j);j=j-1;end do;idx(j+1)=tmp;end do
      h=max(k+1,ceiling(real(n,dp)*bandwidth));m=min(n-k,max(2,h+1));allocate(sel(m),xt(m))
      do i=1,m;sel(i)=u(idx(min(n,k+i)));xt(i)=real(k+i,dp)/real(n-k,dp);end do
      call sort_real(sel)
      slope=sum((xt-sum(xt)/m)*(sel-sum(sel)/m))/max(sum((xt-sum(xt)/m)**2),tiny(1.0_dp))
      density=1.0_dp/max(abs(slope),1.0e-10_dp)
    end if
    if(any(.not.ieee_is_finite(density)))then;status=esback_optimization_failed;else;status=esback_ok;end if
  end subroutine density_quantile_function

  subroutine cdf_at_quantile(y,x,q,cdf,status)
    real(dp),intent(in)::y(:),x(:,:),q(:)
    real(dp),allocatable,intent(out)::cdf(:)
    integer,intent(out)::status
    real(dp),allocatable::mu(:),sigma(:),zs(:),zq(:)
    call conditional_mean_sigma(y,x,mu,sigma,status);if(status/=esback_ok)return
    zs=(y-mu)/sigma;zq=(q-mu)/sigma;allocate(cdf(size(q)));call empirical_cdf_values(zs,zq,cdf)
  end subroutine cdf_at_quantile

  subroutine conditional_truncated_variance(y,x,approach,cv,status)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::approach
    real(dp),allocatable,intent(out)::cv(:)
    integer,intent(out)::status
    real(dp),allocatable::neg(:),mu(:),sigma(:),beta(:),zstd(:),grid(:),pdf(:),cb(:),m1(:),m2(:),mid(:),cvar(:)
    real(dp)::v,bw,lo,hi,h,phi,p
    integer::i,n,ng,st
    n=size(y);if(count(y<=0.0_dp)<=2)then;status=esback_insufficient_data;return;end if
    allocate(neg(count(y<=0.0_dp)));neg=pack(y,y<=0.0_dp);v=sample_variance(neg)
    allocate(cv(n))
    if(approach==sigma_ind)then;cv=v;status=esback_ok;return;end if
    call conditional_mean_sigma(y,x,mu,sigma,st)
    if(st/=esback_ok)then;cv=v;status=esback_ok;return;end if
    beta=-mu/sigma
    if(approach==sigma_scl_n)then
      do i=1,n
        beta(i)=max(beta(i),-30.0_dp);phi=normal_pdf(beta(i));p=max(normal_cdf(beta(i)),1.0e-300_dp)
        cv(i)=sigma(i)**2*max(0.0_dp,1.0_dp-beta(i)*phi/p-(phi/p)**2)
      end do
    else
      zstd=(y-mu)/sigma
      bw = 0.9_dp*min(sample_sd(zstd), &
        (empirical_quantile_type7(zstd,0.75_dp) - &
         empirical_quantile_type7(zstd,0.25_dp))/1.34_dp)* &
        real(n,dp)**(-0.2_dp)
      if(.not.ieee_is_finite(bw) .or. bw<=1.0e-6_dp)bw=1.06_dp*max(sample_sd(zstd),1.0e-3_dp)*real(n,dp)**(-0.2_dp)
      ng=1001;lo=minval(zstd)-4.0_dp*bw;hi=max(maxval(beta),maxval(zstd))+bw;h=(hi-lo)/real(ng-1,dp)
      allocate(grid(ng),pdf(ng),cb(ng-1),m1(ng-1),m2(ng-1),mid(ng-1),cvar(ng-1))
      do i=1,ng
        grid(i)=lo+real(i-1,dp)*h;pdf(i)=sum(exp(-0.5_dp*((grid(i)-zstd)/bw)**2))/(real(n,dp)*bw*sqrt(2.0_dp*pi))
      end do
      cb=0.0_dp;m1=0.0_dp;m2=0.0_dp
      do i=1,ng-1
        mid(i)=grid(i)+0.5_dp*h
        if(i==1)then
          cb(i)=0.5_dp*h*(pdf(i)+pdf(i+1));m1(i)=0.5_dp*h*(grid(i)*pdf(i)+grid(i+1)*pdf(i+1)); &
          m2(i)=0.5_dp*h*(grid(i)**2*pdf(i)+grid(i+1)**2*pdf(i+1))
        else
          cb(i)=cb(i-1)+0.5_dp*h*(pdf(i)+pdf(i+1));m1(i)=m1(i-1)+0.5_dp*h*(grid(i)*pdf(i)+grid(i+1)*pdf(i+1)); &
          m2(i)=m2(i-1)+0.5_dp*h*(grid(i)**2*pdf(i)+grid(i+1)**2*pdf(i+1))
        end if
        cvar(i)=max(0.0_dp,m2(i)/max(cb(i),1.0e-12_dp)-(m1(i)/max(cb(i),1.0e-12_dp))**2)
      end do
      do i=1,n;cv(i)=sigma(i)**2*interpolate_linear(mid,cvar,beta(i));end do
    end if
    if(any(.not.ieee_is_finite(cv)) .or. any(cv<0.0_dp))cv=v
    status=esback_ok
  end subroutine conditional_truncated_variance

  subroutine esreg_covariance(xq,xe,y,alpha,bq,be,opt,cov,status)
    real(dp),intent(in)::xq(:,:),xe(:,:),y(:),alpha,bq(:),be(:)
    type(esreg_options),intent(in)::opt
    real(dp),allocatable,intent(out)::cov(:,:)
    integer,intent(out)::status
    real(dp),allocatable::yt(:),bqt(:),bet(:),q(:),e(:),u(:),dens(:),cdf(:),cv(:),lambda(:,:),sigma(:,:),linv(:,:)
    real(dp)::maxy,a1,a2,a3,qdiff,cdfd
    integer::n,kq,ke,i,st
    n=size(y);kq=size(xq,2);ke=size(xe,2);yt=y;bqt=bq;bet=be;maxy=maxval(y);yt=yt-maxy;bqt(1)=bqt(1)-maxy;bet(1)=bet(1)-maxy
    q=matmul(xq,bqt);e=matmul(xe,bet);u=yt-q
    call density_quantile_function(yt,xq,u,alpha,merge(sparsity_iid,opt%sparsity,kq==1.and.ke==1),opt%bandwidth_estimator,dens,st)
    if(st/=esback_ok)then;status=st;return;end if
    call cdf_at_quantile(yt,xq,q,cdf,st);if(st/=esback_ok)then;status=st;return;end if
    call conditional_truncated_variance(u, xq, &
      merge(sigma_ind, opt%sigma_est, kq == 1 .and. ke == 1), cv, st)
    if (st /= esback_ok) then
      status = st
      return
    end if
    allocate(lambda(kq+ke,kq+ke),sigma(kq+ke,kq+ke));lambda=0.0_dp;sigma=0.0_dp
    do i=1,n
      a1=g1_prime(q(i),2)+g2_value(e(i),1)/alpha
      a2=g2_prime(e(i),1);a3=g2_second(e(i),1);cdfd=cdf(i)-alpha;qdiff=q(i)-e(i)
      lambda(:kq,:kq)=lambda(:kq,:kq)+outer_product(xq(i,:),xq(i,:))*(a1*dens(i)+merge(g1_second(q(i),2)*cdfd,0.0_dp,opt%misspec))
      if(opt%misspec)then
        lambda(:kq,kq+1:)=lambda(:kq,kq+1:)+outer_product(xq(i,:),xe(i,:))*a2*cdfd/alpha
        lambda(kq+1:,:kq)=lambda(kq+1:,:kq)+outer_product(xe(i,:),xq(i,:))*a2*cdfd/alpha
        lambda(kq+1:,kq+1:)=lambda(kq+1:,kq+1:)+outer_product(xe(i,:),xe(i,:))*(a2+a3*q(i)*cdfd/alpha)
      else
        lambda(kq+1:,kq+1:)=lambda(kq+1:,kq+1:)+outer_product(xe(i,:),xe(i,:))*a2
      end if
      a1=alpha*g1_prime(q(i),2)+g2_value(e(i),1)
      if(opt%misspec)then
        sigma(:kq,:kq) = sigma(:kq,:kq) + &
          outer_product(xq(i,:),xq(i,:))*a1*a1* &
          ((1.0_dp-alpha)/alpha + (1.0_dp-2.0_dp*alpha)*cdfd/alpha**2)
        sigma(kq+1:,:kq)=sigma(kq+1:,:kq)+outer_product(xe(i,:),xq(i,:))*a1*a2*((1.0_dp-alpha)/alpha*qdiff+ &
          (1.0_dp-alpha)/alpha*q(i)*cdfd/alpha-cdfd/alpha*qdiff)
        sigma(:kq,kq+1:)=transpose(sigma(kq+1:,:kq))
        sigma(kq+1:,kq+1:)=sigma(kq+1:,kq+1:)+outer_product(xe(i,:),xe(i,:))*a2*a2*(cv(i)/alpha+(1.0_dp-alpha)/alpha*qdiff*qdiff+ &
          2.0_dp*qdiff*q(i)*(alpha-cdf(i))/alpha)
      else
        sigma(:kq,:kq)=sigma(:kq,:kq)+outer_product(xq(i,:),xq(i,:))*a1*a1*(1.0_dp-alpha)/alpha
        sigma(kq+1:,:kq)=sigma(kq+1:,:kq)+outer_product(xe(i,:),xq(i,:))*a1*a2*(1.0_dp-alpha)/alpha*qdiff
        sigma(:kq,kq+1:)=transpose(sigma(kq+1:,:kq))
        sigma(kq+1:,kq+1:)=sigma(kq+1:,kq+1:)+outer_product(xe(i,:),xe(i,:))*a2*a2*(cv(i)/alpha+(1.0_dp-alpha)/alpha*qdiff*qdiff)
      end if
    end do
    lambda=lambda/real(n,dp);sigma=sigma/real(n,dp)
    call invert_matrix(lambda,linv,st);if(st/=esback_ok)then;status=st;return;end if
    cov=matmul(linv,matmul(sigma,linv))/real(n,dp)
    cov=0.5_dp*(cov+transpose(cov));status=esback_ok
  end subroutine esreg_covariance
end module esback_esreg
