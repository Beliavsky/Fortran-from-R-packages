module mev_sampling
  use mev_kinds, only: dp
  use mev_math, only: rng_gamma, rng_normal, mvnormal_sample, chol_lower, normal_pdf
  use mev_spatial, only: lambda2cov
  implicit none
  private
  public :: rdir, rmnorm, mvrt_sample, dmvnorm
  public :: rlogspec, rneglogspec, rbilogspec, rdirmixspec, rhrspec, rdirspec
  public :: rexstudspec, rbrspec
  public :: rmev_spectral, rmev, rparp, rgparp
contains
  real(dp) function runif_open() result(u)
    call random_number(u)
    u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u))
  end function runif_open

  integer function sample_uniform_index(n) result(j)
    integer,intent(in)::n
    j=1+int(runif_open()*real(n,dp))
    j=min(n,max(1,j))
  end function sample_uniform_index

  integer function sample_weighted(w) result(j)
    real(dp),intent(in)::w(:)
    real(dp)::u,s,tot
    integer::i
    tot=sum(max(w,0.0_dp))
    if(tot<=0.0_dp)then;j=1;return;end if
    u=runif_open()*tot;s=0.0_dp
    do i=1,size(w)
      s=s+max(w(i),0.0_dp)
      if(u<=s)then;j=i;return;end if
    end do
    j=size(w)
  end function sample_weighted

  subroutine rdir(n,alpha,sample,normalize,info)
    integer,intent(in)::n
    real(dp),intent(in)::alpha(:)
    real(dp),intent(out)::sample(n,size(alpha))
    logical,intent(in),optional::normalize
    integer,intent(out),optional::info
    logical::norm
    integer::i,j
    real(dp)::s
    norm=.true.;if(present(normalize))norm=normalize
    if(present(info))info=0
    if(n<1.or.any(alpha<=0.0_dp))then
      if(present(info))info=1
      sample=0.0_dp;return
    end if
    do i=1,n
      do j=1,size(alpha)
        sample(i,j)=rng_gamma(alpha(j),1.0_dp)
      end do
      if(norm)then
        s=sum(sample(i,:));sample(i,:)=sample(i,:)/s
      end if
    end do
  end subroutine rdir

  subroutine rmnorm(n,mu,sigma,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),sigma(:,:)
    real(dp),intent(out)::sample(n,size(mu))
    integer,intent(out),optional::info
    integer::ier
    call mvnormal_sample(n,mu,sigma,sample,ier)
    if(present(info))info=ier
  end subroutine rmnorm

  subroutine mvrt_sample(n,scale,dof,loc,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::scale(:,:),dof,loc(:)
    real(dp),intent(out)::sample(n,size(loc))
    integer,intent(out),optional::info
    real(dp)::z(n,size(loc)),chi
    integer::i,ier
    if(present(info))info=0
    if(dof<=0.0_dp)then;if(present(info))info=1;sample=0.0_dp;return;end if
    call mvnormal_sample(n,0.0_dp*loc,scale,z,ier)
    if(ier/=0)then;if(present(info))info=ier;sample=0.0_dp;return;end if
    do i=1,n
      chi=2.0_dp*rng_gamma(0.5_dp*dof,1.0_dp)
      sample(i,:)=loc+z(i,:)*sqrt(dof/chi)
    end do
  end subroutine mvrt_sample

  real(dp) function dmvnorm(x,mu,sigma,log_density) result(v)
    use mev_math, only: inverse_matrix, logdet_spd
    real(dp),intent(in)::x(:),mu(:),sigma(:,:)
    logical,intent(in),optional::log_density
    real(dp)::si(size(x),size(x)),d(size(x)),ld,q
    integer::ier
    logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    call inverse_matrix(sigma,si,ier);if(ier/=0)then;v=0.0_dp;return;end if
    call logdet_spd(sigma,ld,ier);if(ier/=0)then;v=0.0_dp;return;end if
    d=x-mu;q=dot_product(d,matmul(si,d))
    v=-0.5_dp*(real(size(x),dp)*log(2.0_dp*acos(-1.0_dp))+ld+q)
    if(.not.lg)v=exp(v)
  end function dmvnorm

  subroutine rlogspec(n,d,theta,sample,info)
    integer,intent(in)::n,d
    real(dp),intent(in)::theta
    real(dp),intent(out)::sample(n,d)
    integer,intent(out),optional::info
    integer::i,k,j
    real(dp)::f0,s
    if(present(info))info=0
    if(theta<=1.0_dp)then;if(present(info))info=1;sample=0.0_dp;return;end if
    do i=1,n
      j=sample_uniform_index(d)
      f0=exp(-log(rng_gamma(1.0_dp-1.0_dp/theta,1.0_dp))/theta)
      do k=1,d
        sample(i,k)=(-log(runif_open()))**(-1.0_dp/theta)/f0
      end do
      sample(i,j)=1.0_dp;s=sum(sample(i,:));sample(i,:)=sample(i,:)/s
    end do
  end subroutine rlogspec

  subroutine rneglogspec(n,d,theta,sample,info)
    integer,intent(in)::n,d
    real(dp),intent(in)::theta
    real(dp),intent(out)::sample(n,d)
    integer,intent(out),optional::info
    integer::i,k,j
    real(dp)::sc,s
    if(present(info))info=0
    if(theta<=0.0_dp)then;if(present(info))info=1;sample=0.0_dp;return;end if
    sc=1.0_dp/gamma(1.0_dp+1.0_dp/theta)
    do i=1,n
      j=sample_uniform_index(d)
      do k=1,d
        sample(i,k)=sc*(-log(runif_open()))**(1.0_dp/theta)
      end do
      sample(i,j)=sc*rng_gamma(1.0_dp+1.0_dp/theta,1.0_dp)**(1.0_dp/theta)
      s=sum(sample(i,:));sample(i,:)=sample(i,:)/s
    end do
  end subroutine rneglogspec

  subroutine rbilogspec(n,alpha,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::alpha(:)
    real(dp),intent(out)::sample(n,size(alpha))
    integer,intent(out),optional::info
    real(dp)::astar(size(alpha)),tmp(1,size(alpha)),s
    integer::i,k,j,ier,d
    d=size(alpha);if(present(info))info=0
    if(any(alpha<0.0_dp).or.any(alpha>=1.0_dp))then;if(present(info))info=1;sample=0.0_dp;return;end if
    do i=1,n
      astar=1.0_dp;j=sample_uniform_index(d);astar(j)=1.0_dp-alpha(j)
      call rdir(1,astar,tmp,.true.,ier)
      do k=1,d
        sample(i,k)=exp(-alpha(k)*log(tmp(1,k))+log_gamma(real(d,dp)-alpha(k))-log_gamma(1.0_dp-alpha(k)))
      end do
      s=sum(sample(i,:));sample(i,:)=sample(i,:)/s
    end do
  end subroutine rbilogspec

  subroutine rdirmixspec(n,alpha,weight,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::alpha(:,:),weight(:)
    real(dp),intent(out)::sample(n,size(alpha,1))
    integer,intent(out),optional::info
    integer::i,j,k,m,d,nmix
    real(dp)::w(size(weight)),astar(size(alpha,1)),tmp(1,size(alpha,1)),s
    d=size(alpha,1);nmix=size(alpha,2);if(present(info))info=0
    if(nmix/=size(weight).or.any(alpha<=0.0_dp).or.any(weight<0.0_dp).or.sum(weight)<=0.0_dp)then
      if(present(info))info=1;sample=0.0_dp;return
    end if
    do i=1,n
      j=sample_uniform_index(d)
      do k=1,nmix
        w(k)=weight(k)*alpha(j,k)/sum(alpha(:,k))
      end do
      m=sample_weighted(w);astar=alpha(:,m);astar(j)=astar(j)+1.0_dp
      call rdir(1,astar,tmp,.false.)
      sample(i,:)=tmp(1,:);s=sum(sample(i,:));sample(i,:)=sample(i,:)/s
    end do
  end subroutine rdirmixspec

  subroutine rhrspec(n,lambda,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:,:)
    real(dp),intent(out)::sample(n,size(lambda,1))
    integer,intent(out),optional::info
    integer::i,j,k,d,ier,kk
    integer,allocatable::sub(:)
    real(dp),allocatable::cov(:,:),mu(:),z(:,:)
    real(dp)::s
    d=size(lambda,1);if(present(info))info=0
    if(size(lambda,2)/=d.or.any(abs([(lambda(k,k),k=1,d)])>1.0e-10_dp))then
      if(present(info))info=1;sample=0.0_dp;return
    end if
    allocate(sub(d-1),cov(d-1,d-1),mu(d-1),z(1,d-1))
    do i=1,n
      j=sample_uniform_index(d);kk=0
      do k=1,d;if(k/=j)then;kk=kk+1;sub(kk)=k;end if;end do
      call lambda2cov(lambda,j,sub,sub,cov,ier)
      if(ier/=0)then;if(present(info))info=ier;sample=0.0_dp;return;end if
      do k=1,d-1;mu(k)=-2.0_dp*lambda(sub(k),j);end do
      call mvnormal_sample(1,mu,cov,z,ier)
      if(ier/=0)then;if(present(info))info=ier;sample=0.0_dp;return;end if
      sample(i,j)=1.0_dp;kk=0
      do k=1,d
        if(k/=j)then;kk=kk+1;sample(i,k)=exp(z(1,kk));end if
      end do
      s=sum(sample(i,:));sample(i,:)=sample(i,:)/s
    end do
  end subroutine rhrspec

  subroutine rdirspec(n,alpha,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::alpha(:)
    real(dp),intent(out)::sample(:,:)
    integer,intent(out),optional::info
    integer::d,i,j,k
    logical::irv
    real(dp),allocatable::astar(:),tmp(:,:),m(:)
    real(dp)::s,a0
    d=size(sample,2);irv=(size(alpha)==d+1);if(present(info))info=0
    if(size(sample,1)/=n.or.(size(alpha)/=d.and..not.irv).or.any(alpha(1:d)<=0.0_dp))then
      if(present(info))info=1;sample=0.0_dp;return
    end if
    allocate(astar(d),tmp(1,d),m(d));astar=alpha(1:d)
    if(irv)then
      a0=alpha(d+1)
      if(a0<=-minval(alpha(1:d)).or.a0>1.0_dp.or.abs(a0)<1.0e-12_dp)then
        if(present(info))info=2;sample=0.0_dp;return
      end if
      do k=1,d;m(k)=-log_gamma(alpha(k))+log_gamma(alpha(k)+a0);end do
    end if
    do i=1,n
      j=sample_uniform_index(d);astar=alpha(1:d)
      if(irv)then
        astar(j)=astar(j)+a0;call rdir(1,astar,tmp,.false.)
        do k=1,d;sample(i,k)=exp(a0*log(tmp(1,k))-m(k));end do
      else
        astar(j)=astar(j)+1.0_dp;call rdir(1,astar,tmp,.false.)
        do k=1,d;sample(i,k)=tmp(1,k)/alpha(k);end do
      end if
      s=sum(sample(i,:));sample(i,:)=sample(i,:)/s
    end do
  end subroutine rdirspec

  subroutine rexstudspec(n,sigma,alpha,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::sigma(:,:),alpha
    real(dp),intent(out)::sample(n,size(sigma,1))
    integer,intent(out),optional::info
    real(dp),allocatable::cor(:,:),cov(:,:),loc(:),z(:,:),sd(:)
    integer,allocatable::idx(:)
    integer::d,i,j,k,kk,ier
    real(dp)::sm
    d=size(sigma,1);if(present(info))info=0
    if(n<1.or.d<2.or.size(sigma,2)/=d.or.alpha<=0.0_dp)then
      sample=0.0_dp;if(present(info))info=1;return
    end if
    allocate(cor(d,d),cov(d-1,d-1),loc(d-1),z(1,d-1),sd(d),idx(d-1))
    do j=1,d
      if(sigma(j,j)<=0.0_dp)then;sample=0.0_dp;if(present(info))info=2;return;end if
      sd(j)=sqrt(sigma(j,j))
    end do
    do j=1,d;do k=1,d;cor(j,k)=sigma(j,k)/(sd(j)*sd(k));end do;end do
    do i=1,n
      j=sample_uniform_index(d);kk=0
      do k=1,d;if(k/=j)then;kk=kk+1;idx(kk)=k;end if;end do
      do k=1,d-1
        loc(k)=cor(idx(k),j);do kk=1,d-1
          cov(k,kk)=(cor(idx(k),idx(kk))-cor(idx(k),j)*cor(j,idx(kk)))/(alpha+1.0_dp)
        end do
      end do
      call mvrt_sample(1,cov,alpha+1.0_dp,loc,z,ier)
      if(ier/=0)then;sample=0.0_dp;if(present(info))info=ier;return;end if
      sample(i,j)=1.0_dp;kk=0
      do k=1,d
        if(k/=j)then
          kk=kk+1;sample(i,k)=max(z(1,kk),0.0_dp)**alpha
        end if
      end do
      sm=sum(sample(i,:));if(sm<=0.0_dp)then;sample=0.0_dp;if(present(info))info=3;return;end if
      sample(i,:)=sample(i,:)/sm
    end do
  end subroutine rexstudspec

  subroutine rbrspec(n,sigma,sample,info)
    integer,intent(in)::n
    real(dp),intent(in)::sigma(:,:)
    real(dp),intent(out)::sample(n,size(sigma,1))
    integer,intent(out),optional::info
    real(dp),allocatable::z(:,:),mu(:)
    integer::d,i,j,k,ier
    real(dp)::sm,v
    d=size(sigma,1);if(present(info))info=0
    if(n<1.or.d<2.or.size(sigma,2)/=d)then
      sample=0.0_dp;if(present(info))info=1;return
    end if
    allocate(z(n,d),mu(d));mu=0.0_dp;call mvnormal_sample(n,mu,sigma,z,ier)
    if(ier/=0)then;sample=0.0_dp;if(present(info))info=ier;return;end if
    do i=1,n
      j=sample_uniform_index(d)
      do k=1,d
        v=sigma(k,k)+sigma(j,j)-2.0_dp*sigma(k,j)
        sample(i,k)=exp(z(i,k)-z(i,j)-0.5_dp*v)
      end do
      sample(i,j)=1.0_dp;sm=sum(sample(i,:));sample(i,:)=sample(i,:)/sm
    end do
  end subroutine rbrspec

  subroutine rmev_spectral(n,model,par,sample,lambda,alpha_mix,weights,info,sigma)
    integer,intent(in)::n
    character(len=*),intent(in)::model
    real(dp),intent(in),optional::par(:),lambda(:,:),alpha_mix(:,:),weights(:),sigma(:,:)
    real(dp),intent(out)::sample(:,:)
    integer,intent(out),optional::info
    integer::d,ier
    d=size(sample,2);ier=0
    select case(trim(model))
    case('log')
      if(.not.present(par))then
        ier=1
      else if(par(1)<1.0_dp.and.par(1)>0.0_dp)then
        call rlogspec(n,d,1.0_dp/par(1),sample,ier)
      else
        call rlogspec(n,d,par(1),sample,ier)
      end if
    case('neglog')
      if(.not.present(par))then;ier=1;else;call rneglogspec(n,d,par(1),sample,ier);end if
    case('bilog')
      if(.not.present(par))then
        ier=1
      else if(all(par>=1.0_dp))then
        call rbilogspec(n,1.0_dp/par,sample,ier)
      else
        call rbilogspec(n,par,sample,ier)
      end if
    case('xstud','schlather')
      if(.not.present(sigma))then
        ier=1
      else if(trim(model)=='schlather')then
        call rexstudspec(n,sigma,1.0_dp,sample,ier)
      else if(.not.present(par))then
        ier=1
      else
        call rexstudspec(n,sigma,par(1),sample,ier)
      end if
    case('br')
      if(.not.present(sigma))then;ier=1;else;call rbrspec(n,sigma,sample,ier);end if
    case('ct','sdir')
      if(.not.present(par))then;ier=1;else;call rdirspec(n,par,sample,ier);end if
    case('hr')
      if(.not.present(lambda))then;ier=1;else;call rhrspec(n,lambda,sample,ier);end if
    case('dirmix')
      if(.not.present(alpha_mix).or..not.present(weights))then
        ier=1
      else
        call rdirmixspec(n,alpha_mix,weights,sample,ier)
      end if
    case default
      ier=2;sample=0.0_dp
    end select
    if(present(info))info=ier
  end subroutine rmev_spectral

  subroutine rmev(n,model,par,sample,lambda,alpha_mix,weights,info,sigma)
    integer,intent(in)::n
    character(len=*),intent(in)::model
    real(dp),intent(in),optional::par(:),lambda(:,:),alpha_mix(:,:),weights(:),sigma(:,:)
    real(dp),intent(out)::sample(:,:)
    integer,intent(out),optional::info
    real(dp),allocatable::y(:,:)
    real(dp)::zeta
    integer::i,d,ier
    if(present(info))info=0
    d=size(sample,2)
    if(size(sample,1)/=n.or.n<1.or.d<1)then
      sample=0.0_dp
      if(present(info))info=1
      return
    end if
    allocate(y(1,d))
    sample=0.0_dp
    do i=1,n
      zeta=-log(runif_open())/real(d,dp)
      do while(1.0_dp/zeta>minval(sample(i,:)))
        call rmev_spectral(1,model,par,y,lambda,alpha_mix,weights,ier,sigma)
        if(ier/=0)then
          sample=0.0_dp
          if(present(info))info=ier
          return
        end if
        sample(i,:)=max(sample(i,:),y(1,:)/zeta)
        zeta=zeta-log(runif_open())/real(d,dp)
      end do
    end do
  end subroutine rmev

  subroutine rparp(n,shape_tail,risk,model,sample,par,lambda,alpha_mix,weights,sigma,siteindex,accept_rate,info,max_trials)
    use mev_distributions, only: qgp
    integer,intent(in)::n
    real(dp),intent(in)::shape_tail
    character(len=*),intent(in)::risk,model
    real(dp),intent(out)::sample(:,:)
    real(dp),intent(in),optional::par(:),lambda(:,:),alpha_mix(:,:),weights(:),sigma(:,:)
    integer,intent(in),optional::siteindex,max_trials
    real(dp),intent(out),optional::accept_rate
    integer,intent(out),optional::info
    real(dp),allocatable::ang(:,:),cand(:)
    real(dp)::u,rad
    integer::d,got,trials,limit,ier,site
    logical::acc
    d=size(sample,2);sample=0.0_dp;if(present(info))info=0;if(present(accept_rate))accept_rate=0.0_dp
    if(n<1.or.size(sample,1)/=n.or.d<2.or.shape_tail<=0.0_dp)then;if(present(info))info=1;return;end if
    site=1;if(present(siteindex))site=siteindex
    if(trim(risk)=='site'.and.(site<1.or.site>d))then;if(present(info))info=2;return;end if
    limit=max(100000,10000*n);if(present(max_trials))limit=max(n,max_trials)
    allocate(ang(1,d),cand(d));got=0;trials=0
    do while(got<n.and.trials<limit)
      trials=trials+1
      call rmev_spectral(1,model,par,ang,lambda,alpha_mix,weights,ier,sigma)
      if(ier/=0)then;if(present(info))info=10+ier;return;end if
      u=runif_open();rad=qgp(u,loc=1.0_dp,scale=1.0_dp,shape=shape_tail);cand=rad*ang(1,:)
      acc=risk_exceeded(cand,risk,1.0_dp,site)
      if(acc)then;got=got+1;sample(got,:)=cand;end if
    end do
    if(trials>0.and.present(accept_rate))accept_rate=real(got,dp)/real(trials,dp)
    if(got<n.and.present(info))info=3
  end subroutine rparp

  subroutine rgparp(n,shape,thresh,risk,model,loc,scale,sample,par,lambda,alpha_mix,weights,sigma, &
      siteindex,accept_rate,info,max_trials)
    integer,intent(in)::n
    real(dp),intent(in)::shape(:),thresh,loc(:),scale(:)
    character(len=*),intent(in)::risk,model
    real(dp),intent(out)::sample(:,:)
    real(dp),intent(in),optional::par(:),lambda(:,:),alpha_mix(:,:),weights(:),sigma(:,:)
    integer,intent(in),optional::siteindex,max_trials
    real(dp),intent(out),optional::accept_rate
    integer,intent(out),optional::info
    real(dp),allocatable::ang(:,:),cand(:),xi(:)
    real(dp)::rad
    integer::d,j,got,trials,limit,ier,site
    logical::acc
    d=size(sample,2);sample=0.0_dp;if(present(info))info=0;if(present(accept_rate))accept_rate=0.0_dp
    if(n<1.or.size(sample,1)/=n.or.d<2.or.size(loc)/=d.or.size(scale)/=d.or.any(scale<=0.0_dp))then
      if(present(info))info=1;return
    end if
    if(size(shape)/=1.and.size(shape)/=d)then;if(present(info))info=2;return;end if
    allocate(xi(d));if(size(shape)==1)then;xi=shape(1);else;xi=shape;end if
    site=1;if(present(siteindex))site=siteindex
    if(trim(risk)=='site'.and.(site<1.or.site>d))then;if(present(info))info=3;return;end if
    limit=max(100000,10000*n);if(present(max_trials))limit=max(n,max_trials)
    allocate(ang(1,d),cand(d));got=0;trials=0
    do while(got<n.and.trials<limit)
      trials=trials+1;call rmev_spectral(1,model,par,ang,lambda,alpha_mix,weights,ier,sigma)
      if(ier/=0)then;if(present(info))info=10+ier;return;end if
      rad=1.0_dp/runif_open();cand=rad*ang(1,:)
      do j=1,d
        if(abs(xi(j))<1.0e-10_dp)then
          cand(j)=loc(j)+scale(j)*log(cand(j))
        else
          cand(j)=loc(j)+scale(j)*(cand(j)**xi(j)-1.0_dp)/xi(j)
        end if
      end do
      acc=risk_exceeded(cand,risk,thresh,site)
      if(acc)then;got=got+1;sample(got,:)=cand;end if
    end do
    if(trials>0.and.present(accept_rate))accept_rate=real(got,dp)/real(trials,dp)
    if(got<n.and.present(info))info=4
  end subroutine rgparp

  pure logical function risk_exceeded(x,risk,thresh,site) result(ok)
    real(dp),intent(in)::x(:),thresh
    character(len=*),intent(in)::risk
    integer,intent(in)::site
    select case(trim(risk))
    case('sum');ok=sum(x)>thresh
    case('mean');ok=sum(x)/real(size(x),dp)>thresh
    case('max');ok=maxval(x)>thresh
    case('min');ok=minval(x)>thresh
    case('l2');ok=sum(x*x)>thresh
    case('site');ok=x(site)>thresh
    case default;ok=.false.
    end select
  end function risk_exceeded
end module mev_sampling
