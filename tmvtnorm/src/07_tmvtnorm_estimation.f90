! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm_estimation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mvtnorm_kinds, only : dp
  use mvtnorm_probabilities, only : mvn_prob => pmvnorm
  use mvtnorm_types, only : probability_result
  use mvtnorm_distributions, only : dmvnorm_one
  use mvtnorm_linalg, only : inverse_spd
  use mvtnorm_special, only : chi_square_cdf
  use tmvtnorm_utils, only : vech, invech, pack_cholesky, unpack_cholesky, covariance_ok
  use tmvtnorm_marginals, only : dtmvnorm_marginal
  use tmvtnorm_moments, only : tmvnorm_moments_t, mtmvnorm
  use tmvtnorm_optimize, only : nm_result_t, nelder_mead
  implicit none
  private
  public :: tmvnorm_fit_t, tmvnorm_gmm_fit_t, mle_tmvnorm, gmm_tmvnorm
  public :: gmm_moments_manjunath_wilhelm, gmm_moments_lee
  public :: gmm_mw, gmm_lee
  integer,parameter :: gmm_mw=1, gmm_lee=2

  type :: tmvnorm_fit_t
    real(dp),allocatable :: mean(:),sigma(:,:),theta(:)
    real(dp)::negloglik=huge(1.0_dp)
    integer::iterations=0,evaluations=0,convergence=1
  end type tmvnorm_fit_t

  type :: tmvnorm_gmm_fit_t
    real(dp),allocatable :: mean(:),sigma(:,:),theta(:),weight(:,:),moments(:,:)
    real(dp)::objective=huge(1.0_dp),j_stat=0.0_dp,j_pvalue=1.0_dp
    integer::iterations=0,evaluations=0,convergence=1,df=0
  end type tmvnorm_gmm_fit_t

  real(dp),allocatable,save :: sx(:,:),sl(:),su(:),base_theta(:),sw(:,:)
  logical,allocatable,save :: free_mask_saved(:)
  logical,save :: scholesky=.false.
  integer,save :: sn=0,sgmm_method=gmm_mw,slmax=1

contains

  subroutine mle_tmvnorm(x,lower,upper,mean0,sigma0,res,cholesky,free_mask,maxit,tol)
    real(dp),intent(in)::x(:,:),lower(:),upper(:),mean0(:),sigma0(:,:)
    type(tmvnorm_fit_t),intent(out)::res
    logical,intent(in),optional::cholesky,free_mask(:)
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    real(dp),allocatable::full0(:),red0(:)
    integer,allocatable::idx(:)
    integer::p,mi
    real(dp)::eps
    type(nm_result_t)::op
    logical::chol
    p=size(mean0)
    chol=.false.
    if(present(cholesky)) chol=cholesky
    mi=2000
    if(present(maxit)) mi=maxit
    eps=1.0e-7_dp
    if(present(tol)) eps=tol
    sx=x
    sl=lower
    su=upper
    sn=p
    scholesky=chol
    if(chol) then
    full0=[mean0,pack_cholesky(sigma0)]
    else
    full0=[mean0,vech(sigma0)]
    end if
    call setup_free(full0,free_mask,idx,red0)
    op=nelder_mead(mle_objective_reduced,red0,maxit=mi,reltol=eps,step=0.03_dp)
    res%theta=expand_free(op%par)
    call unpack_theta(res%theta,p,chol,res%mean,res%sigma)
    res%negloglik=op%value
    res%iterations=op%iterations
    res%evaluations=op%evaluations
    res%convergence=op%convergence
    call clear_saved()
  end subroutine mle_tmvnorm

  function mle_objective_reduced(red) result(v)
    real(dp),intent(in)::red(:)
    real(dp)::v
    real(dp),allocatable::theta(:),mu(:),sig(:,:)
    type(probability_result)::pr
    integer::i
    theta=expand_free(red)
    call unpack_theta(theta,sn,scholesky,mu,sig)
    if(.not.covariance_ok(sig)) then
    v=huge(1.0_dp)/100.0_dp
    return
    end if
    pr=mvn_prob(sl,su,mu,sig)
    if(pr%value<=tiny(1.0_dp)) then
    v=huge(1.0_dp)/100.0_dp
    return
    end if
    v=real(size(sx,1),dp)*log(pr%value)
    do i=1,size(sx,1)
      if(any(sx(i,:)<sl .or. sx(i,:)>su)) then
      v=huge(1.0_dp)/100.0_dp
      return
      end if
      v=v-dmvnorm_one(sx(i,:),mu,sig,.true.)
    end do
  end function mle_objective_reduced

  subroutine gmm_tmvnorm(x,lower,upper,mean0,sigma0,res,method,cholesky,lmax,free_mask,maxit,tol)
    real(dp),intent(in)::x(:,:),lower(:),upper(:),mean0(:),sigma0(:,:)
    type(tmvnorm_gmm_fit_t),intent(out)::res
    integer,intent(in),optional::method,lmax,maxit
    logical,intent(in),optional::cholesky,free_mask(:)
    real(dp),intent(in),optional::tol
    real(dp),allocatable::full0(:),red0(:),theta1(:),g(:,:),gc(:,:),s(:,:),sinv(:,:)
    integer,allocatable::idx(:)
    integer::p,q,mi,i
    real(dp)::eps,ridge
    logical::chol,ok
    character(len=256)::msg
    type(nm_result_t)::op1,op2
    p=size(mean0)
    chol=.false.
    if(present(cholesky)) chol=cholesky
    sgmm_method=gmm_mw
    if(present(method)) sgmm_method=method
    slmax=(p+2)/2
    if(present(lmax)) slmax=lmax
    mi=1500
    if(present(maxit)) mi=maxit
    eps=1.0e-6_dp
    if(present(tol)) eps=tol
    sx=x
    sl=lower
    su=upper
    sn=p
    scholesky=chol
    if(chol) then
    full0=[mean0,pack_cholesky(sigma0)]
    else
    full0=[mean0,vech(sigma0)]
    end if
    call setup_free(full0,free_mask,idx,red0)
    if(sgmm_method==gmm_mw) then
    q=p+p*(p+1)/2
    else
    q=(slmax+1)*p
    end if
    allocate(sw(q,q))
    sw=0.0_dp
    do i=1,q
    sw(i,i)=1.0_dp
    end do
    op1=nelder_mead(gmm_objective_reduced,red0,maxit=mi,reltol=eps,step=0.03_dp)
    theta1=expand_free(op1%par)
    g=moment_matrix(theta1)
    gc=g-spread(sum(g,dim=1)/real(size(g,1),dp),1,size(g,1))
    s=matmul(transpose(gc),gc)/real(size(g,1),dp)
    ridge=max(1.0e-10_dp,1.0e-8_dp*sum([(s(i,i),i=1,q)])/real(max(1,q),dp))
    do i=1,q
    s(i,i)=s(i,i)+ridge
    end do
    call inverse_spd(s,sinv,ok,msg)
    if(ok) sw=sinv
    op2=nelder_mead(gmm_objective_reduced,op1%par,maxit=mi,reltol=eps,step=0.02_dp)
    res%theta=expand_free(op2%par)
    call unpack_theta(res%theta,p,chol,res%mean,res%sigma)
    res%moments=moment_matrix(res%theta)
    res%weight=sw
    res%objective=op2%value
    res%df=q-size(op2%par)
    res%j_stat=real(size(x,1),dp)*res%objective
    if(res%df>0) res%j_pvalue=1.0_dp-chi_square_cdf(res%j_stat,real(res%df,dp))
    res%iterations=op1%iterations+op2%iterations
    res%evaluations=op1%evaluations+op2%evaluations
    res%convergence=op2%convergence
    call clear_saved()
  end subroutine gmm_tmvnorm

  function gmm_objective_reduced(red) result(v)
    real(dp),intent(in)::red(:)
    real(dp)::v
    real(dp),allocatable::theta(:),g(:,:),gb(:)
    theta=expand_free(red)
    g=moment_matrix(theta)
    if(any(abs(g)>=huge(1.0_dp)/1000.0_dp)) then
    v=huge(1.0_dp)/100.0_dp
    return
    end if
    gb=sum(g,dim=1)/real(size(g,1),dp)
    v=dot_product(gb,matmul(sw,gb))
  end function gmm_objective_reduced

  function moment_matrix(theta) result(g)
    real(dp),intent(in)::theta(:)
    real(dp),allocatable::g(:,:),mu(:),sig(:,:)
    call unpack_theta(theta,sn,scholesky,mu,sig)
    if(.not.covariance_ok(sig)) then
      if(sgmm_method==gmm_mw) then
      allocate(g(size(sx,1),sn+sn*(sn+1)/2))
      else
      allocate(g(size(sx,1),(slmax+1)*sn))
      end if
      g=huge(1.0_dp)/100.0_dp
      return
    end if
    if(sgmm_method==gmm_mw) then
    g=gmm_moments_manjunath_wilhelm(sx,mu,sig,sl,su)
    else
    g=gmm_moments_lee(sx,mu,sig,sl,su,slmax)
    end if
  end function moment_matrix

  function gmm_moments_manjunath_wilhelm(x,mean,sigma,lower,upper) result(g)
    real(dp),intent(in)::x(:,:),mean(:),sigma(:,:),lower(:),upper(:)
    real(dp),allocatable::g(:,:)
    type(tmvnorm_moments_t)::mom
    integer::n,t,i,j,k
    n=size(mean)
    t=size(x,1)
    allocate(g(t,n+n*(n+1)/2))
    call mtmvnorm(mean,sigma,lower,upper,mom)
    if(.not.mom%ok) then
    g=huge(1.0_dp)/100.0_dp
    return
    end if
    do i=1,n
    g(:,i)=mom%mean(i)-x(:,i)
    end do
    k=0
    do i=1,n
      do j=1,i
        k=k+1
        g(:,n+k)=(mom%mean(i)-x(:,i))*(mom%mean(j)-x(:,j))-mom%covariance(i,j)
      end do
    end do
  end function gmm_moments_manjunath_wilhelm

  function gmm_moments_lee(x,mean,sigma,lower,upper,lmax) result(g)
    real(dp),intent(in)::x(:,:),mean(:),sigma(:,:),lower(:),upper(:)
    integer,intent(in)::lmax
    real(dp),allocatable::g(:,:),h(:,:),fa(:),fb(:)
    integer::n,t,l,i,k,obs
    real(dp)::ail,bil,boundary,xi_pow,deriv,term1,term2
    logical::ok
    character(len=256)::msg
    n=size(mean)
    t=size(x,1)
    allocate(g(t,(lmax+1)*n),fa(n),fb(n))
    call inverse_spd(sigma,h,ok,msg)
    if(.not.ok) then
    g=huge(1.0_dp)/100.0_dp
    return
    end if
    do i=1,n
      fa(i)=dtmvnorm_marginal(lower(i),i,mean,sigma,lower,upper)
      fb(i)=dtmvnorm_marginal(upper(i),i,mean,sigma,lower,upper)
    end do
    k=0
    do l=0,lmax
      do i=1,n
        k=k+1
        if(ieee_is_finite(lower(i))) then
        ail=lower(i)**l
        else
        ail=0.0_dp
        end if
        if(ieee_is_finite(upper(i))) then
        bil=upper(i)**l
        else
        bil=0.0_dp
        end if
        boundary=ail*fa(i)-bil*fb(i)
        term1=dot_product(h(i,:),mean)
        do obs=1,t
          xi_pow=x(obs,i)**l
          term2=dot_product(x(obs,:),h(:,i))
          if(l==0) then
          deriv=0.0_dp
          else
          deriv=real(l,dp)*x(obs,i)**(l-1)
          end if
          g(obs,k)=term1*xi_pow-xi_pow*term2+deriv+boundary
        end do
      end do
    end do
  end function gmm_moments_lee

  subroutine unpack_theta(theta,n,chol,mean,sigma)
    real(dp),intent(in)::theta(:)
    integer,intent(in)::n
    logical,intent(in)::chol
    real(dp),allocatable,intent(out)::mean(:),sigma(:,:)
    mean=theta(1:n)
    if(chol) then
    sigma=unpack_cholesky(theta(n+1:),n)
    else
    sigma=invech(theta(n+1:))
    end if
  end subroutine unpack_theta

  subroutine setup_free(full0,mask,idx,red0)
    real(dp),intent(in)::full0(:)
    logical,intent(in),optional::mask(:)
    integer,allocatable,intent(out)::idx(:)
    real(dp),allocatable,intent(out)::red0(:)
    integer::i,k
    base_theta=full0
    allocate(free_mask_saved(size(full0)))
    free_mask_saved=.true.
    if(present(mask)) free_mask_saved=mask
    allocate(idx(count(free_mask_saved)),red0(count(free_mask_saved)))
    k=0
    do i=1,size(full0)
    if(free_mask_saved(i)) then
    k=k+1
    idx(k)=i
    red0(k)=full0(i)
    end if
    end do
  end subroutine setup_free

  function expand_free(red) result(full)
    real(dp),intent(in)::red(:)
    real(dp),allocatable::full(:)
    integer::i,k
    full=base_theta
    k=0
    do i=1,size(full)
    if(free_mask_saved(i)) then
    k=k+1
    full(i)=red(k)
    end if
    end do
  end function expand_free

  subroutine clear_saved()
    if(allocated(sx)) deallocate(sx)
    if(allocated(sl)) deallocate(sl)
    if(allocated(su)) deallocate(su)
    if(allocated(base_theta)) deallocate(base_theta)
    if(allocated(free_mask_saved)) deallocate(free_mask_saved)
    if(allocated(sw)) deallocate(sw)
  end subroutine clear_saved

end module tmvtnorm_estimation
