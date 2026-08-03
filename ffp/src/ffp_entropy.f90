module ffp_entropy
  use ffp_kinds, only : dp
  use ffp_linalg, only : solve_linear, covariance_to_correlation
  use ffp_probabilities, only : normalize_probabilities, exp_decay_probabilities
  use ffp_statistics, only : weighted_moments
  implicit none
  private
  public :: entropy_pool_equalities, entropy_pool_constrained, least_info_kernel
  public :: fit_to_moments, double_decay_covariance, double_decay_probabilities
contains
  subroutine entropy_pool_equalities(prior,aeq,beq,post,info,tol,max_iter)
    real(dp),intent(in)::prior(:),aeq(:,:),beq(:)
    real(dp),intent(out)::post(:)
    integer,intent(out)::info
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::max_iter
    real(dp),allocatable::v(:),r(:),h(:,:),delta(:),trial(:),q(:)
    real(dp)::eps,alpha,norm0,norm1
    integer::m,n,it,itmax,j,stat
    n=size(prior); m=size(beq); eps=1.0e-10_dp; if(present(tol)) eps=tol
    itmax=200; if(present(max_iter)) itmax=max_iter
    info=0
    if(size(aeq,1)/=m .or. size(aeq,2)/=n .or. size(post)/=n) then; info=-1; return; end if
    allocate(v(m),r(m),h(m,m),delta(m),trial(m),q(n)); v=0.0_dp
    do it=1,itmax
      call probabilities(v,q)
      r=matmul(aeq,q)-beq; norm0=maxval(abs(r))
      if(norm0<eps) then; post=q; call normalize_probabilities(post); return; end if
      h=0.0_dp
      do j=1,n
        h=h+q(j)*spread(aeq(:,j),2,m)*spread(aeq(:,j),1,m)
      end do
      do j=1,m; h(j,j)=h(j,j)+1.0e-12_dp; end do
      call solve_linear(h,r,delta,stat)
      if(stat/=0) then; info=10+stat; exit; end if
      alpha=1.0_dp
      do
        trial=v+alpha*delta; call probabilities(trial,q)
        norm1=maxval(abs(matmul(aeq,q)-beq))
        if(norm1<norm0 .or. alpha<1.0e-8_dp) exit
        alpha=0.5_dp*alpha
      end do
      v=trial
    end do
    call probabilities(v,post); call normalize_probabilities(post); if(info==0) info=1
  contains
    subroutine probabilities(lam,out)
      real(dp),intent(in)::lam(:); real(dp),intent(out)::out(:)
      real(dp)::zmax; real(dp),allocatable::z(:)
      allocate(z(n)); z=log(max(prior,1.0e-300_dp))-matmul(transpose(aeq),lam)
      zmax=maxval(z); out=exp(max(z-zmax,-700.0_dp)); call normalize_probabilities(out)
    end subroutine
  end subroutine entropy_pool_equalities

  subroutine entropy_pool_constrained(prior,a,b,aeq,beq,post,info,tol,max_iter)
    real(dp),intent(in)::prior(:),a(:,:),b(:),aeq(:,:),beq(:)
    real(dp),intent(out)::post(:); integer,intent(out)::info
    real(dp),intent(in),optional::tol; integer,intent(in),optional::max_iter
    real(dp),allocatable::g(:),q(:),viol(:),eqr(:)
    real(dp)::eps,eta,rho,zmax
    integer::it,itmax,n
    n=size(prior); eps=1.0e-8_dp; if(present(tol))eps=tol
    itmax=20000; if(present(max_iter))itmax=max_iter
    allocate(g(n),q(n),viol(size(b)),eqr(size(beq))); q=prior; call normalize_probabilities(q)
    rho=1000.0_dp; eta=0.02_dp; info=1
    do it=1,itmax
      viol=max(matmul(a,q)-b,0.0_dp); eqr=matmul(aeq,q)-beq
      if(max(maxval(viol),maxval(abs(eqr)))<eps) then; info=0; exit; end if
      g=log(max(q,1.0e-300_dp)/max(prior,1.0e-300_dp))+1.0_dp
      if(size(b)>0) g=g+rho*matmul(transpose(a),viol)
      if(size(beq)>0) g=g+rho*matmul(transpose(aeq),eqr)
      g=-eta*g; zmax=maxval(g); q=q*exp(max(g-zmax,-700.0_dp)); call normalize_probabilities(q)
      if(mod(it,2000)==0) eta=0.7_dp*eta
    end do
    post=q
  end subroutine entropy_pool_constrained

  subroutine least_info_kernel(y,target,h2,p,info)
    real(dp),intent(in)::y(:,:),target(:)
    real(dp),intent(in),optional::h2(:,:)
    real(dp),intent(out)::p(:); integer,intent(out)::info
    real(dp),allocatable::aeq(:,:),beq(:),prior(:)
    integer::n,k,m,i,j,row
    n=size(y,1); k=size(y,2); m=1+k
    if(present(h2)) m=m+k*(k+1)/2
    allocate(aeq(m,n),beq(m),prior(n)); aeq(1,:)=1.0_dp; beq(1)=1.0_dp
    aeq(2:1+k,:)=transpose(y); beq(2:1+k)=target; row=1+k
    if(present(h2)) then
      do i=1,k; do j=i,k; row=row+1; aeq(row,:)=y(:,i)*y(:,j); beq(row)=h2(i,j)+target(i)*target(j); end do; end do
    end if
    prior=1.0_dp/real(n,dp); call entropy_pool_equalities(prior,aeq,beq,p,info)
  end subroutine

  subroutine fit_to_moments(x,mean,sigma,p,info)
    real(dp),intent(in)::x(:,:),mean(:),sigma(:,:); real(dp),intent(out)::p(:); integer,intent(out)::info
    call least_info_kernel(x,mean,sigma,p,info)
  end subroutine

  subroutine double_decay_covariance(x,decay_low,decay_high,mean,sigma)
    real(dp),intent(in)::x(:,:),decay_low,decay_high
    real(dp),intent(out)::mean(:),sigma(:,:)
    real(dp),allocatable::pl(:),ph(:),ml(:),mh(:),sl(:,:),sh(:,:),cor(:,:),sd(:)
    integer::n,k,i,j
    n=size(x,1); k=size(x,2); allocate(pl(n),ph(n),ml(k),mh(k),sl(k,k),sh(k,k),cor(k,k),sd(k))
    call exp_decay_probabilities(n,decay_low,pl); call exp_decay_probabilities(n,decay_high,ph)
    call weighted_moments(x,pl,ml,sl); call weighted_moments(x,ph,mh,sh)
    call covariance_to_correlation(sl,cor,sd); do i=1,k; sd(i)=sqrt(max(sh(i,i),0.0_dp)); end do
    do j=1,k; do i=1,k; sigma(i,j)=sd(i)*cor(i,j)*sd(j); end do; end do
    mean=0.0_dp
  end subroutine

  subroutine double_decay_probabilities(x,decay_low,decay_high,p,info)
    real(dp),intent(in)::x(:,:),decay_low,decay_high; real(dp),intent(out)::p(:); integer,intent(out)::info
    real(dp),allocatable::m(:),s(:,:); integer::k
    k=size(x,2); allocate(m(k),s(k,k)); call double_decay_covariance(x,decay_low,decay_high,m,s)
    call fit_to_moments(x,m,s,p,info)
  end subroutine
end module ffp_entropy
