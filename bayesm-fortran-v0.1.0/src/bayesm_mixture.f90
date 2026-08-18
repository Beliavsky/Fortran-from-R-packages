module bayesm_mixture
  use bayesm_kinds, only: dp, pi
  use bayesm_linalg, only: chol_upper, inverse_upper, inverse_spd, correlation_matrix, outer_product
  use bayesm_rng, only: randn, rand_categorical, rand_dirichlet, rwishart_draw
  use bayesm_densities, only: lnd_mvn
  use bayesm_regression, only: rmultireg
  use bayesm_types, only: normal_component, normal_mixture, mixture_step_result, nmix_result, multireg_draw
  implicit none
  private
  public :: rmixture, mix_den, mix_den_bi, mixture_moments, mom_mix, e_mix_marg_den
  public :: rmix_gibbs, rnmix_gibbs, cluster_mix
contains
  subroutine component_covariance(comp,sigma,info)
    type(normal_component), intent(in) :: comp
    real(dp), intent(out) :: sigma(size(comp%mu),size(comp%mu))
    integer, intent(out) :: info
    real(dp) :: root(size(comp%mu),size(comp%mu))
    call inverse_upper(comp%rooti,root,info)
    if (info==0) sigma=matmul(transpose(root),root)
  end subroutine component_covariance

  subroutine rmixture(n,mix,x,z,info)
    integer, intent(in) :: n
    type(normal_mixture), intent(in) :: mix
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, allocatable, intent(out) :: z(:)
    integer, intent(out) :: info
    integer :: i,j,d,k,stat
    real(dp), allocatable :: root(:,:),zn(:)
    d=size(mix%comp(1)%mu)
    allocate(x(n,d),z(n),root(d,d),zn(d)); info=0
    do i=1,n
      k=rand_categorical(mix%p); z(i)=k
      call inverse_upper(mix%comp(k)%rooti,root,stat)
      if (stat/=0) then; info=stat; return; end if
      do j=1,d; zn(j)=randn(); end do
      x(i,:)=mix%comp(k)%mu+matmul(transpose(root),zn)
    end do
  end subroutine rmixture

  function mix_den(x,mix) result(den)
    real(dp), intent(in) :: x(:,:)
    type(normal_mixture), intent(in) :: mix
    real(dp) :: den(size(x,1),size(x,2))
    real(dp) :: sigma(size(x,2),size(x,2)),sd
    integer :: i,j,k,stat
    den=0.0_dp
    do k=1,size(mix%comp)
      call component_covariance(mix%comp(k),sigma,stat)
      if (stat/=0) cycle
      do j=1,size(x,2)
        sd=sqrt(max(tiny(1.0_dp),sigma(j,j)))
        do i=1,size(x,1)
          den(i,j)=den(i,j)+mix%p(k)*exp(-0.5_dp*((x(i,j)-mix%comp(k)%mu(j))/sd)**2)/ &
            (sqrt(2.0_dp*pi)*sd)
        end do
      end do
    end do
  end function mix_den

  function mix_den_bi(iv,jv,xi,xj,mix) result(den)
    integer, intent(in) :: iv,jv
    real(dp), intent(in) :: xi(:),xj(:)
    type(normal_mixture), intent(in) :: mix
    real(dp) :: den(size(xi),size(xj))
    real(dp) :: sigma(size(mix%comp(1)%mu),size(mix%comp(1)%mu)),sub(2,2),rooti(2,2),r(2,2)
    real(dp) :: z(2),mu(2),ld
    integer :: a,b,k,stat
    den=0.0_dp
    do k=1,size(mix%comp)
      call component_covariance(mix%comp(k),sigma,stat); if (stat/=0) cycle
      sub(1,1)=sigma(iv,iv); sub(1,2)=sigma(iv,jv); sub(2,1)=sigma(jv,iv); sub(2,2)=sigma(jv,jv)
      call chol_upper(sub,r,stat); if (stat/=0) cycle
      call inverse_upper(r,rooti,stat); if (stat/=0) cycle
      mu=[mix%comp(k)%mu(iv),mix%comp(k)%mu(jv)]
      do b=1,size(xj)
        do a=1,size(xi)
          z=[xi(a),xj(b)]; ld=lnd_mvn(z,mu,rooti)
          den(a,b)=den(a,b)+mix%p(k)*exp(ld)
        end do
      end do
    end do
  end function mix_den_bi

  subroutine mixture_moments(mix,mu,sigma,sd,corr,info)
    type(normal_mixture), intent(in) :: mix
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:),sd(:),corr(:,:)
    integer, intent(out) :: info
    integer :: d,k,stat
    real(dp), allocatable :: sk(:,:),delta(:)
    d=size(mix%comp(1)%mu)
    allocate(mu(d),sigma(d,d),sd(d),corr(d,d),sk(d,d),delta(d)); mu=0.0_dp; sigma=0.0_dp
    do k=1,size(mix%comp)
      mu=mu+mix%p(k)*mix%comp(k)%mu
    end do
    info=0
    do k=1,size(mix%comp)
      call component_covariance(mix%comp(k),sk,stat)
      if (stat/=0) then; info=stat; return; end if
      delta=mix%comp(k)%mu-mu
      sigma=sigma+mix%p(k)*(sk+outer_product(delta,delta))
    end do
    sd=sqrt(max(0.0_dp,[(sigma(k,k),k=1,d)]))
    corr=correlation_matrix(sigma)
  end subroutine mixture_moments

  subroutine mom_mix(mixdraw,mu,sigma,sd,corr,info)
    type(normal_mixture), intent(in) :: mixdraw(:)
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:),sd(:),corr(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: m(:),s(:,:),st(:),c(:,:)
    integer :: d,i,stat
    d=size(mixdraw(1)%comp(1)%mu)
    allocate(mu(d),sigma(d,d),sd(d),corr(d,d)); mu=0.0_dp; sigma=0.0_dp; sd=0.0_dp; corr=0.0_dp
    info=0
    do i=1,size(mixdraw)
      call mixture_moments(mixdraw(i),m,s,st,c,stat)
      if (stat/=0) then; info=stat; return; end if
      mu=mu+m; sigma=sigma+s; sd=sd+st; corr=corr+c
    end do
    mu=mu/real(size(mixdraw),dp); sigma=sigma/real(size(mixdraw),dp)
    sd=sd/real(size(mixdraw),dp); corr=corr/real(size(mixdraw),dp)
  end subroutine mom_mix

  function e_mix_marg_den(grid,mixdraw) result(den)
    real(dp), intent(in) :: grid(:,:)
    type(normal_mixture), intent(in) :: mixdraw(:)
    real(dp) :: den(size(grid,1),size(grid,2))
    integer :: i
    den=0.0_dp
    do i=1,size(mixdraw)
      den=den+mix_den(grid,mixdraw(i))
    end do
    den=den/real(size(mixdraw),dp)
  end function e_mix_marg_den

  function rmix_gibbs(y,bbar,a,nu,v,p0,z0) result(out)
    real(dp), intent(in) :: y(:,:),bbar(:,:),a(:,:),nu,v(:,:),p0(:)
    integer, intent(in) :: z0(:)
    type(mixture_step_result) :: out
    integer :: n,d,ncomp,k,i,j,nk,stat
    real(dp), allocatable :: yk(:,:),xk(:,:),ainv(:,:),w(:,:),iw(:,:),c(:,:),ci(:,:),mu(:),r(:,:),rooti(:,:)
    real(dp), allocatable :: logp(:),prob(:),alpha(:)
    type(multireg_draw) :: mr
    n=size(y,1); d=size(y,2); ncomp=size(p0)
    allocate(out%p(ncomp),out%z(n),out%comp(ncomp),alpha(ncomp))
    do k=1,ncomp
      nk=count(z0==k)
      if (nk>0) then
        allocate(yk(nk,d),xk(nk,1)); j=0
        do i=1,n
          if (z0(i)==k) then; j=j+1; yk(j,:)=y(i,:); end if
        end do
        xk=1.0_dp
        mr=rmultireg(yk,xk,bbar,a,nu,v,stat)
        allocate(out%comp(k)%mu(d),out%comp(k)%rooti(d,d)); out%comp(k)%mu=mr%b(1,:)
        allocate(r(d,d)); call chol_upper(mr%sigma,r,stat); call inverse_upper(r,out%comp(k)%rooti,stat)
        deallocate(yk,xk,r)
      else
        allocate(ainv(d,d),w(d,d),iw(d,d),c(d,d),ci(d,d),mu(d),r(d,d),rooti(d,d))
        call inverse_spd(v,ainv,stat); call rwishart_draw(nu,ainv,w,iw,c,ci,stat)
        do i=1,d; mu(i)=randn(); end do
        mu=bbar(1,:)+matmul(ci,mu)/sqrt(a(1,1))
        call chol_upper(iw,r,stat); call inverse_upper(r,rooti,stat)
        allocate(out%comp(k)%mu(d),out%comp(k)%rooti(d,d)); out%comp(k)%mu=mu; out%comp(k)%rooti=rooti
        deallocate(ainv,w,iw,c,ci,mu,r,rooti)
      end if
    end do
    allocate(logp(ncomp),prob(ncomp))
    do i=1,n
      do k=1,ncomp
        logp(k)=log(max(tiny(1.0_dp),p0(k)))+lnd_mvn(y(i,:),out%comp(k)%mu,out%comp(k)%rooti)
      end do
      logp=logp-maxval(logp); prob=exp(logp); out%z(i)=rand_categorical(prob)
    end do
    alpha=0.0_dp
    do k=1,ncomp; alpha(k)=real(count(out%z==k),dp); end do
    alpha=alpha+1.0_dp
    ! Preserve user prior concentration approximately when p0 is only an initial probability vector.
    call rand_dirichlet(alpha,out%p)
  end function rmix_gibbs

  function rnmix_gibbs(y,bbar,a,nu,v,adir,nrep,keep,pinit,zinit) result(out)
    real(dp), intent(in) :: y(:,:),bbar(:,:),a(:,:),nu,v(:,:),adir(:)
    integer, intent(in) :: nrep,keep
    real(dp), intent(in), optional :: pinit(:)
    integer, intent(in), optional :: zinit(:)
    type(nmix_result) :: out
    integer :: ns,n,ncomp,rep,mkeep,i,k
    real(dp), allocatable :: p(:),alpha(:)
    integer, allocatable :: z(:)
    type(mixture_step_result) :: step
    ns=nrep/keep; n=size(y,1); ncomp=size(adir)
    allocate(p(ncomp),z(n),alpha(ncomp)); alpha=adir
    if (present(pinit)) then; p=pinit; else; p=1.0_dp/real(ncomp,dp); end if
    if (present(zinit)) then
      z=zinit
    else
      do i=1,n; z(i)=mod(i-1,ncomp)+1; end do
    end if
    allocate(out%probdraw(ns,ncomp),out%zdraw(ns,n),out%mixdraw(ns)); mkeep=0
    do rep=1,nrep
      step=rmix_gibbs_prior(y,bbar,a,nu,v,alpha,p,z)
      p=step%p; z=step%z
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%probdraw(mkeep,:)=p; out%zdraw(mkeep,:)=z
        allocate(out%mixdraw(mkeep)%p(ncomp),out%mixdraw(mkeep)%comp(ncomp)); out%mixdraw(mkeep)%p=p
        do k=1,ncomp
          allocate(out%mixdraw(mkeep)%comp(k)%mu(size(step%comp(k)%mu)))
          allocate(out%mixdraw(mkeep)%comp(k)%rooti(size(step%comp(k)%rooti,1),size(step%comp(k)%rooti,2)))
          out%mixdraw(mkeep)%comp(k)%mu=step%comp(k)%mu
          out%mixdraw(mkeep)%comp(k)%rooti=step%comp(k)%rooti
        end do
      end if
    end do
  end function rnmix_gibbs

  function rmix_gibbs_prior(y,bbar,a,nu,v,adir,p0,z0) result(out)
    real(dp), intent(in) :: y(:,:),bbar(:,:),a(:,:),nu,v(:,:),adir(:),p0(:)
    integer, intent(in) :: z0(:)
    type(mixture_step_result) :: out
    integer :: n,d,ncomp,k,i,j,nk,stat
    real(dp), allocatable :: yk(:,:),xk(:,:),vinv(:,:),w(:,:),iw(:,:),c(:,:),ci(:,:),mu(:),r(:,:),rooti(:,:)
    real(dp), allocatable :: logp(:),prob(:),alpha(:)
    type(multireg_draw) :: mr
    n=size(y,1); d=size(y,2); ncomp=size(p0)
    allocate(out%p(ncomp),out%z(n),out%comp(ncomp),alpha(ncomp))
    do k=1,ncomp
      nk=count(z0==k)
      if (nk>0) then
        allocate(yk(nk,d),xk(nk,1)); j=0
        do i=1,n
          if (z0(i)==k) then; j=j+1; yk(j,:)=y(i,:); end if
        end do
        xk=1.0_dp; mr=rmultireg(yk,xk,bbar,a,nu,v,stat)
        allocate(out%comp(k)%mu(d),out%comp(k)%rooti(d,d)); out%comp(k)%mu=mr%b(1,:)
        allocate(r(d,d)); call chol_upper(mr%sigma,r,stat); call inverse_upper(r,out%comp(k)%rooti,stat)
        deallocate(yk,xk,r)
      else
        allocate(vinv(d,d),w(d,d),iw(d,d),c(d,d),ci(d,d),mu(d),r(d,d),rooti(d,d))
        call inverse_spd(v,vinv,stat); call rwishart_draw(nu,vinv,w,iw,c,ci,stat)
        do i=1,d; mu(i)=randn(); end do
        mu=bbar(1,:)+matmul(ci,mu)/sqrt(a(1,1))
        call chol_upper(iw,r,stat); call inverse_upper(r,rooti,stat)
        allocate(out%comp(k)%mu(d),out%comp(k)%rooti(d,d)); out%comp(k)%mu=mu; out%comp(k)%rooti=rooti
        deallocate(vinv,w,iw,c,ci,mu,r,rooti)
      end if
    end do
    allocate(logp(ncomp),prob(ncomp))
    do i=1,n
      do k=1,ncomp
        logp(k)=log(max(tiny(1.0_dp),p0(k)))+lnd_mvn(y(i,:),out%comp(k)%mu,out%comp(k)%rooti)
      end do
      logp=logp-maxval(logp); prob=exp(logp); out%z(i)=rand_categorical(prob)
    end do
    alpha=adir
    do k=1,ncomp; alpha(k)=alpha(k)+real(count(out%z==k),dp); end do
    call rand_dirichlet(alpha,out%p)
  end function rmix_gibbs_prior

  subroutine cluster_mix(zdraw,cutoff,clustera,clusterb,pmean)
    integer, intent(in) :: zdraw(:,:)
    real(dp), intent(in) :: cutoff
    integer, allocatable, intent(out) :: clustera(:),clusterb(:)
    real(dp), allocatable, intent(out), optional :: pmean(:,:)
    integer :: nr,n,i,j,r,best,groupn,countn
    real(dp), allocatable :: pm(:,:),sim(:,:),loss(:)
    logical, allocatable :: assigned(:)
    nr=size(zdraw,1); n=size(zdraw,2)
    allocate(pm(n,n),sim(n,n),loss(nr),clustera(n),clusterb(n),assigned(n)); pm=0.0_dp
    do r=1,nr
      do j=1,n; do i=1,n
        if (zdraw(r,i)==zdraw(r,j)) pm(i,j)=pm(i,j)+1.0_dp
      end do; end do
    end do
    pm=pm/real(nr,dp)
    do r=1,nr
      loss(r)=0.0_dp
      do j=1,n; do i=1,n
        if (zdraw(r,i)==zdraw(r,j)) then
          loss(r)=loss(r)+abs(pm(i,j)-1.0_dp)
        else
          loss(r)=loss(r)+abs(pm(i,j))
        end if
      end do; end do
    end do
    best=minloc(loss,dim=1); clustera=zdraw(best,:)
    sim=0.0_dp; where(pm>=cutoff) sim=1.0_dp
    clusterb=0; assigned=.false.; groupn=1
    do j=1,n
      countn=0
      do i=1,n
        if (.not.assigned(i) .and. sim(i,j)>=0.5_dp) then
          clusterb(i)=groupn; assigned(i)=.true.; countn=countn+1
        end if
      end do
      if (countn>0) groupn=groupn+1
    end do
    if (present(pmean)) then; allocate(pmean(n,n)); pmean=pm; end if
  end subroutine cluster_mix
end module bayesm_mixture
