module bayesm_blp
  use bayesm_kinds, only: dp, log2pi
  use bayesm_linalg, only: chol_upper, inverse_upper, inverse_spd, solve_spd, logdet_spd
  use bayesm_rng, only: randn, rand_uniform, rand_chisq, rwishart_draw
  use bayesm_densities, only: lnd_mvn
  use bayesm_regression, only: breg
  use bayesm_types, only: blp_result
  implicit none
  private
  public :: r2sigma, share2mu, log_jacob, rbayes_blp
contains
  pure function r2sigma(rvec,k) result(sigma)
    real(dp), intent(in) :: rvec(:)
    integer, intent(in) :: k
    real(dp) :: sigma(k,k),l(k,k)
    integer :: i,j,q
    if (size(rvec)/=k*(k+1)/2) error stop "r2sigma: wrong r length"
    l=0.0_dp
    do i=1,k
      l(i,i)=exp(rvec(i))
    end do
    q=k
    do i=1,k-1
      do j=i+1,k
        q=q+1; l(j,i)=rvec(q)
      end do
    end do
    sigma=matmul(l,transpose(l))
  end function r2sigma

  subroutine share2mu(sigma,x,v,share,jalt,tol,mu,choice_prob,info,maxiter)
    real(dp), intent(in) :: sigma(:,:),x(:,:),v(:,:),share(:),tol
    integer, intent(in) :: jalt
    real(dp), allocatable, intent(out) :: mu(:),choice_prob(:,:)
    integer, intent(out) :: info
    integer, intent(in), optional :: maxiter
    integer :: n,k,h,t,j,it,ilim,stat
    real(dp), allocatable :: r(:,:),u(:,:),expv(:,:),m0(:),m1(:),shat(:),den(:)
    real(dp) :: mx,rel
    n=size(x,1); k=size(x,2); h=size(v,2)
    if (mod(n,jalt)/=0 .or. size(share)/=n .or. size(v,1)/=k) then
      info=-1; allocate(mu(0),choice_prob(0,0)); return
    end if
    ilim=10000; if (present(maxiter)) ilim=maxiter
    allocate(r(k,k),u(n,h),expv(n,h),m0(n),m1(n),shat(n),den(h),mu(n),choice_prob(n,h))
    call chol_upper(sigma,r,stat)
    if (stat/=0) then; info=stat; mu=0.0_dp; choice_prob=0.0_dp; return; end if
    u=matmul(x,matmul(transpose(r),v)); m0=1.0_dp; m1=0.5_dp; it=0; rel=huge(1.0_dp)
    do while(rel>tol .and. it<ilim)
      m0=m1
      do t=1,n/jalt
        den=1.0_dp
        do j=1,jalt
          expv((t-1)*jalt+j,:)=exp(min(700.0_dp,u((t-1)*jalt+j,:)+m0((t-1)*jalt+j)))
          den=den+expv((t-1)*jalt+j,:)
        end do
        do j=1,jalt
          choice_prob((t-1)*jalt+j,:)=expv((t-1)*jalt+j,:)/den
        end do
      end do
      shat=sum(choice_prob,dim=2)/real(h,dp)
      m1=m0+log(max(share,tiny(1.0_dp))/max(shat,tiny(1.0_dp)))
      mx=0.0_dp
      do j=1,n
        mx=max(mx,abs(m1(j)-m0(j))/max(abs(m0(j)),1.0e-8_dp))
      end do
      rel=mx; it=it+1
    end do
    mu=m1; info=merge(0,1,it<ilim)
  end subroutine share2mu

  function log_jacob(choice_prob,jalt,info) result(val)
    real(dp), intent(in) :: choice_prob(:,:)
    integer, intent(in) :: jalt
    integer, intent(out), optional :: info
    real(dp) :: val
    integer :: n,h,t,a,b,stat
    real(dp), allocatable :: block(:,:),pbar(:)
    real(dp) :: ld
    n=size(choice_prob,1); h=size(choice_prob,2); val=0.0_dp; stat=0
    allocate(block(jalt,jalt),pbar(jalt))
    do t=1,n/jalt
      pbar=sum(choice_prob((t-1)*jalt+1:t*jalt,:),dim=2)/real(h,dp)
      do a=1,jalt
        do b=1,jalt
          block(a,b)=-sum(choice_prob((t-1)*jalt+a,:)*choice_prob((t-1)*jalt+b,:))/real(h,dp)
          if (a==b) block(a,b)=sum(choice_prob((t-1)*jalt+a,:)* &
            (1.0_dp-choice_prob((t-1)*jalt+a,:)))/real(h,dp)
        end do
      end do
      ld=logdet_spd(block,stat)
      if (stat/=0) then; val=-huge(1.0_dp); exit; end if
      val=val-ld
    end do
    if (present(info)) info=stat
  end function log_jacob

  subroutine iv_draw_once(mu,xend,z,xexo,theta_hat,a,deltabar,ad,v,nu,delta,omega,theta,info)
    real(dp), intent(in) :: mu(:),xend(:),z(:,:),xexo(:,:),theta_hat(:),a(:,:),deltabar(:),ad(:,:),v(:,:),nu
    real(dp), intent(inout) :: delta(:),omega(:,:)
    real(dp), intent(out) :: theta(:)
    integer, intent(out) :: info
    integer :: n,dimd,dimg,i,j,stat
    real(dp), allocatable :: e1(:),ee2(:),xt(:,:),bg(:),gamma(:),u(:),cm(:,:),bm(:,:),ru(:,:),lm(:,:)
    real(dp), allocatable :: yt(:),xtd(:,:),pair(:),sol(:),tmp(:),res(:,:),ss(:,:),vinv(:,:),wd(:,:),iw(:,:),cc(:,:),ci(:,:)
    real(dp) :: beta,sc
    n=size(mu); dimd=size(z,2); dimg=size(xexo,2)
    allocate(e1(n),ee2(n),xt(n,dimg+1),bg(dimg+1),gamma(dimg),u(n),cm(2,2),bm(2,2),ru(2,2),lm(2,2))
    allocate(yt(2*n),xtd(2*n,dimd),pair(2),sol(2),tmp(2),res(n,2),ss(2,2),vinv(2,2),wd(2,2),iw(2,2),cc(2,2),ci(2,2))
    e1=xend-matmul(z,delta); ee2=(omega(1,2)/omega(1,1))*e1
    sc=sqrt(max(tiny(1.0_dp),omega(2,2)-omega(1,2)*omega(1,2)/omega(1,1)))
    xt(:,1)=xend/sc; if (dimg>0) xt(:,2:dimg+1)=xexo/sc
    bg=breg((mu-ee2)/sc,xt,theta_hat,a,stat); beta=bg(1); if (dimg>0) gamma=bg(2:dimg+1)
    cm=0.0_dp; cm(1,1)=1.0_dp; cm(2,1)=beta; cm(2,2)=1.0_dp
    bm=matmul(matmul(cm,omega),transpose(cm)); call chol_upper(bm,ru,stat); lm=transpose(ru)
    if (dimg>0) then; u=mu-matmul(xexo,gamma); else; u=mu; end if
    do i=1,n
      pair=[xend(i),u(i)]; call solve_lower_local(lm,pair,sol,stat); yt(2*i-1:2*i)=sol
      do j=1,dimd
        tmp=[z(i,j),beta*z(i,j)]; call solve_lower_local(lm,tmp,sol,stat)
        xtd(2*i-1,j)=sol(1); xtd(2*i,j)=sol(2)
      end do
    end do
    delta=breg(yt,xtd,deltabar,ad,stat)
    res(:,1)=xend-matmul(z,delta)
    if (dimg>0) then; res(:,2)=mu-beta*xend-matmul(xexo,gamma); else; res(:,2)=mu-beta*xend; end if
    ss=matmul(transpose(res),res); call inverse_spd(v+ss,vinv,stat)
    call rwishart_draw(nu+real(n,dp),vinv,wd,iw,cc,ci,stat); omega=iw
    if (dimg>0) theta(1:dimg)=gamma
    theta(dimg+1)=beta; info=stat
  end subroutine iv_draw_once

  pure subroutine solve_lower_local(l,b,x,info)
    real(dp), intent(in) :: l(:,:),b(:)
    real(dp), intent(out) :: x(size(b))
    integer, intent(out) :: info
    integer :: i
    x=b; info=0
    do i=1,size(b)
      if (abs(l(i,i))<=tiny(1.0_dp)) then; info=i; return; end if
      if (i>1) x(i)=x(i)-dot_product(l(i,1:i-1),x(1:i-1))
      x(i)=x(i)/l(i,i)
    end do
  end subroutine solve_lower_local

  function rbayes_blp(x,share,jalt,vdraw,sigmasqr,a,theta_hat,nrep,keep,ssq,cand_cov,tol, &
      theta0,r0,tau0,iv,z,deltabar,ad,nu0,s0sq,vomega,delta0,omega0) result(out)
    real(dp), intent(in) :: x(:,:),share(:),vdraw(:,:),sigmasqr(:),a(:,:),theta_hat(:),ssq,cand_cov(:,:),tol
    integer, intent(in) :: jalt,nrep,keep
    real(dp), intent(in) :: theta0(:),r0(:),tau0
    logical, intent(in), optional :: iv
    real(dp), intent(in), optional :: z(:,:),deltabar(:),ad(:,:),nu0,s0sq,vomega(:,:),delta0(:),omega0(:,:)
    type(blp_result) :: out
    integer :: k,nr,ns,rep,mkeep,i,stat,nobs
    logical :: doiv
    real(dp), allocatable :: theta(:),rold(:),rnew(:),sold(:,:),snew(:,:),muold(:),munew(:),cp(:,:),cpnew(:,:),proot(:,:)
    real(dp), allocatable :: rc(:,:),zn(:),etaold(:),etanew(:),prec(:,:),rhs(:),mean(:),cholp(:,:),ri(:,:)
    real(dp), allocatable :: delta(:),omega(:,:),xexo(:,:),xend(:),pair(:),rootio(:,:)
    real(dp) :: tau,llold,llnew,jold,jnew,prold,prnew,alpha,nu1,s1,acceptn
    k=size(x,2); nobs=size(x,1); nr=k*(k+1)/2; ns=nrep/keep; doiv=.false.; if (present(iv)) doiv=iv
    allocate(theta(k),rold(nr),rnew(nr),sold(k,k),snew(k,k),rc(nr,nr),proot(nr,nr),zn(nr),etaold(nobs),etanew(nobs))
    allocate(prec(k,k),rhs(k),mean(k),cholp(k,k),ri(k,k)); theta=theta0; rold=r0; tau=tau0; rc=cand_cov
    sold=r2sigma(rold,k); call share2mu(sold,x,vdraw,share,jalt,tol,muold,cp,stat); jold=log_jacob(cp,jalt,stat)
    if (doiv) then
      if (.not.(present(z).and.present(deltabar).and.present(ad).and.present(vomega).and.present(delta0).and.present(omega0))) &
        error stop "rbayes_blp: missing IV arguments"
      allocate(delta(size(delta0)),omega(2,2),xexo(nobs,k-1),xend(nobs),pair(2),rootio(2,2))
      delta=delta0; omega=omega0; xexo=x(:,1:k-1); xend=x(:,k)
    end if
    allocate(out%thetadraw(ns,k),out%rdraw(ns,nr),out%tausqdraw(ns),out%omegadraw(ns,2,2))
    if (doiv) then; allocate(out%deltadraw(ns,size(delta))); else; allocate(out%deltadraw(ns,1)); end if
    allocate(out%llike(ns)); out%omegadraw=0.0_dp; out%deltadraw=0.0_dp; acceptn=0.0_dp; mkeep=0
    do rep=1,nrep
      call chol_upper(ssq*rc,proot,stat)
      do i=1,nr; zn(i)=randn(); end do
      rnew=rold+matmul(transpose(proot),zn); snew=r2sigma(rnew,k)
      call share2mu(snew,x,vdraw,share,jalt,tol,munew,cpnew,stat); jnew=log_jacob(cpnew,jalt,stat)
      etaold=muold-matmul(x,theta); etanew=munew-matmul(x,theta)
      if (doiv) then
        call chol_upper(omega,cholp(1:2,1:2),stat); call inverse_upper(cholp(1:2,1:2),rootio,stat)
        llold=jold; llnew=jnew
        do i=1,nobs
          pair=[xend(i)-dot_product(z(i,:),delta),etaold(i)]; llold=llold+lnd_mvn(pair,[0.0_dp,0.0_dp],rootio)
          pair(2)=etanew(i); llnew=llnew+lnd_mvn(pair,[0.0_dp,0.0_dp],rootio)
        end do
      else
        llold=jold-0.5_dp*real(nobs,dp)*(log2pi+log(tau))-0.5_dp*sum(etaold*etaold)/tau
        llnew=jnew-0.5_dp*real(nobs,dp)*(log2pi+log(tau))-0.5_dp*sum(etanew*etanew)/tau
      end if
      prold=-0.5_dp*sum(rold*rold/sigmasqr)-0.5_dp*sum(log(sigmasqr))
      prnew=-0.5_dp*sum(rnew*rnew/sigmasqr)-0.5_dp*sum(log(sigmasqr))
      alpha=min(1.0_dp,exp(min(0.0_dp,llnew+prnew-llold-prold)))
      if (rand_uniform()<=alpha) then
        rold=rnew; sold=snew; muold=munew; cp=cpnew; jold=jnew; llold=llnew; acceptn=acceptn+1.0_dp
      end if
      if (doiv) then
        call iv_draw_once(muold,xend,z,xexo,theta_hat,a,deltabar,ad,vomega,nu0,delta,omega,theta,stat)
      else
        prec=matmul(transpose(x),x)/tau+a; rhs=matmul(transpose(x),muold)/tau+matmul(a,theta_hat)
        call solve_spd(prec,rhs,mean,stat); call chol_upper(prec,cholp(1:k,1:k),stat)
        call inverse_upper(cholp(1:k,1:k),ri,stat)
        do i=1,k; zn(i)=randn(); end do
        theta=mean+matmul(ri,zn(1:k)); etaold=muold-matmul(x,theta)
        nu1=nu0+real(nobs,dp); s1=(nu0*s0sq+sum(etaold*etaold))/nu1; tau=nu1*s1/rand_chisq(nu1)
      end if
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%thetadraw(mkeep,:)=theta; out%rdraw(mkeep,:)=rold; out%tausqdraw(mkeep)=tau
        if (doiv) then; out%omegadraw(mkeep,:,:)=omega; out%deltadraw(mkeep,:)=delta; end if
        out%llike(mkeep)=llold
      end if
    end do
    out%accept=acceptn/real(nrep,dp)
  end function rbayes_blp
end module bayesm_blp
