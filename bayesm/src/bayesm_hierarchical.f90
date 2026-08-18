module bayesm_hierarchical
  use bayesm_kinds, only: dp
  use bayesm_linalg, only: inverse_spd, chol_upper, inverse_upper, solve_spd
  use bayesm_rng, only: randn, rand_uniform, rand_chisq
  use bayesm_densities, only: lnd_mvn
  use bayesm_regression, only: rmultireg
  use bayesm_mixture, only: rnmix_gibbs
  use bayesm_mnl, only: llmnl, binlogit_loglike
  use bayesm_types, only: reg_data, mnl_data, hier_linear_result, hier_mnl_result, &
    multireg_draw, nmix_result, normal_mixture
  implicit none
  private
  public :: rhier_linear_model, rhier_bin_logit, rhier_mnl_rw_mixture
  public :: rhier_linear_mixture, rhier_mnl_dp
contains
  subroutine delta_vec_to_matrix(dvec,nz,nvar,delta)
    real(dp), intent(in) :: dvec(:)
    integer, intent(in) :: nz,nvar
    real(dp), intent(out) :: delta(nz,nvar)
    integer :: q,j
    if (size(dvec)/=nz*nvar) error stop "delta_vec_to_matrix: size mismatch"
    q=0
    do j=1,nvar
      delta(:,j)=dvec(q+1:q+nz); q=q+nz
    end do
  end subroutine delta_vec_to_matrix


  function draw_delta(zmat,betas,zz,mix,deltabar,ad,info) result(delta)
    real(dp), intent(in) :: zmat(:,:),betas(:,:),deltabar(:),ad(:,:)
    integer, intent(in) :: zz(:)
    type(normal_mixture), intent(in) :: mix
    integer, intent(out), optional :: info
    real(dp) :: delta(size(zmat,2),size(betas,2))
    integer :: n,nz,p,i,a,b,r,s,ia,ib,k,stat,q
    real(dp), allocatable :: prec(:,:),rhs(:),sigi(:,:),ym(:),meanv(:),root(:,:),ri(:,:),zn(:),dvec(:)
    n=size(betas,1); nz=size(zmat,2); p=size(betas,2)
    allocate(prec(nz*p,nz*p),rhs(nz*p),sigi(p,p),ym(p),meanv(nz*p))
    allocate(root(nz*p,nz*p),ri(nz*p,nz*p),zn(nz*p),dvec(nz*p))
    prec=ad; rhs=matmul(ad,deltabar)
    do i=1,n
      k=zz(i); sigi=matmul(mix%comp(k)%rooti,transpose(mix%comp(k)%rooti)); ym=betas(i,:)-mix%comp(k)%mu
      do a=1,p
        do r=1,nz
          ia=(a-1)*nz+r
          do b=1,p
            rhs(ia)=rhs(ia)+sigi(a,b)*ym(b)*zmat(i,r)
            do s=1,nz
              ib=(b-1)*nz+s
              prec(ia,ib)=prec(ia,ib)+sigi(a,b)*zmat(i,r)*zmat(i,s)
            end do
          end do
        end do
      end do
    end do
    call solve_spd(prec,rhs,meanv,stat); call chol_upper(prec,root,stat)
    if (stat==0) call inverse_upper(root,ri,stat)
    if (stat==0) then
      do q=1,nz*p; zn(q)=randn(); end do
      dvec=meanv+matmul(ri,zn)
    else
      dvec=meanv
    end if
    call delta_vec_to_matrix(dvec,nz,p,delta)
    if (present(info)) info=stat
  end function draw_delta

  function rhier_linear_model(data,z,deltabar,a,nu,v,nu_e,ssq,nrep,keep) result(out)
    type(reg_data), intent(in) :: data(:)
    real(dp), intent(in) :: z(:,:),deltabar(:,:),a(:,:),nu,v(:,:),nu_e,ssq(:)
    integer, intent(in) :: nrep,keep
    type(hier_linear_result) :: out
    integer :: nreg,nvar,nz,ns,rep,mkeep,i,j,stat
    real(dp), allocatable :: betas(:,:),tau(:),delta(:,:),vbeta(:,:),abeta(:,:),betabar(:,:)
    real(dp), allocatable :: pmat(:,:),rhs(:),meanb(:),root(:,:),rooti(:,:),zn(:),res(:)
    type(multireg_draw) :: mr
    nreg=size(data); nvar=size(data(1)%x,2); nz=size(z,2); ns=nrep/keep
    allocate(betas(nreg,nvar),tau(nreg),delta(nz,nvar),vbeta(nvar,nvar),abeta(nvar,nvar))
    allocate(betabar(nreg,nvar),pmat(nvar,nvar),rhs(nvar),meanb(nvar),root(nvar,nvar),rooti(nvar,nvar))
    allocate(zn(nvar)); betas=0.0_dp; delta=deltabar; vbeta=0.0_dp
    do i=1,nvar; vbeta(i,i)=1.0_dp; end do
    do i=1,nreg; tau(i)=max(ssq(i),tiny(1.0_dp)); end do
    allocate(out%betadraw(nreg,nvar,ns),out%deltadraw(ns,nz,nvar))
    allocate(out%vbetadraw(ns,nvar,nvar),out%taudraw(ns,nreg)); mkeep=0
    do rep=1,nrep
      call inverse_spd(vbeta,abeta,stat); betabar=matmul(z,delta)
      do i=1,nreg
        pmat=matmul(transpose(data(i)%x),data(i)%x)/tau(i)+abeta
        rhs=matmul(transpose(data(i)%x),data(i)%y)/tau(i)+matmul(abeta,betabar(i,:))
        call solve_spd(pmat,rhs,meanb,stat); call chol_upper(pmat,root,stat); call inverse_upper(root,rooti,stat)
        do j=1,nvar; zn(j)=randn(); end do; betas(i,:)=meanb+matmul(rooti,zn)
        allocate(res(size(data(i)%y))); res=data(i)%y-matmul(data(i)%x,betas(i,:))
        tau(i)=(nu_e*ssq(i)+sum(res*res))/rand_chisq(nu_e+real(size(res),dp)); deallocate(res)
      end do
      mr=rmultireg(betas,z,deltabar,a,nu,v,stat); delta=mr%b; vbeta=mr%sigma
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(:,:,mkeep)=betas; out%deltadraw(mkeep,:,:)=delta
        out%vbetadraw(mkeep,:,:)=vbeta; out%taudraw(mkeep,:)=tau
      end if
    end do
  end function rhier_linear_model

  function rhier_linear_mixture(data,z,deltabar,ad,mubar,amu,nu,v,nu_e,ssq,adir,nrep,keep,drawdelta) result(out)
    type(reg_data), intent(in) :: data(:)
    real(dp), intent(in) :: z(:,:),deltabar(:),ad(:,:),mubar(:),amu,nu,v(:,:),nu_e,ssq(:),adir(:)
    integer, intent(in) :: nrep,keep
    logical, intent(in), optional :: drawdelta
    type(hier_linear_result) :: out
    integer :: nreg,nvar,nz,ncomp,ns,rep,mkeep,i,j,k,stat
    logical :: dd
    real(dp), allocatable :: betas(:,:),tau(:),delta(:,:),resid(:,:),p(:),pmat(:,:),rhs(:),meanb(:)
    real(dp), allocatable :: root(:,:),ri(:,:),zn(:),res(:),bbar(:,:),amat(:,:)
    integer, allocatable :: zz(:)
    type(nmix_result) :: nm
    type(normal_mixture) :: mix
    dd=.true.; if (present(drawdelta)) dd=drawdelta
    nreg=size(data); nvar=size(data(1)%x,2); nz=size(z,2); ncomp=size(adir); ns=nrep/keep
    allocate(betas(nreg,nvar),tau(nreg),delta(nz,nvar),resid(nreg,nvar),p(ncomp),zz(nreg))
    allocate(pmat(nvar,nvar),rhs(nvar),meanb(nvar),root(nvar,nvar),ri(nvar,nvar),zn(nvar))
    allocate(bbar(1,nvar),amat(1,1)); bbar(1,:)=mubar; amat(1,1)=amu
    betas=0.0_dp; call delta_vec_to_matrix(deltabar,nz,nvar,delta); p=1.0_dp/real(ncomp,dp)
    do i=1,nreg; tau(i)=max(ssq(i),tiny(1.0_dp)); zz(i)=mod(i-1,ncomp)+1; end do
    allocate(out%betadraw(nreg,nvar,ns),out%deltadraw(ns,nz,nvar),out%vbetadraw(ns,nvar,nvar))
    allocate(out%taudraw(ns,nreg),out%probdraw(ns,ncomp),out%zdraw(ns,nreg),out%mixdraw(ns)); out%vbetadraw=0.0_dp
    mkeep=0
    do rep=1,nrep
      if (dd) then; resid=betas-matmul(z,delta); else; resid=betas; end if
      nm=rnmix_gibbs(resid,bbar,amat,nu,v,adir,1,1,p,zz); p=nm%probdraw(1,:); zz=nm%zdraw(1,:); mix=nm%mixdraw(1)
      if (dd) delta=draw_delta(z,betas,zz,mix,deltabar,ad,stat)
      do i=1,nreg
        k=zz(i); pmat=matmul(transpose(data(i)%x),data(i)%x)/tau(i)+ &
          matmul(mix%comp(k)%rooti,transpose(mix%comp(k)%rooti))
        meanb=mix%comp(k)%mu; if (dd) meanb=meanb+matmul(z(i,:),delta)
        rhs=matmul(transpose(data(i)%x),data(i)%y)/tau(i)+ &
          matmul(matmul(mix%comp(k)%rooti,transpose(mix%comp(k)%rooti)),meanb)
        call solve_spd(pmat,rhs,meanb,stat); call chol_upper(pmat,root,stat); call inverse_upper(root,ri,stat)
        do j=1,nvar; zn(j)=randn(); end do; betas(i,:)=meanb+matmul(ri,zn)
        allocate(res(size(data(i)%y))); res=data(i)%y-matmul(data(i)%x,betas(i,:))
        tau(i)=(nu_e*ssq(i)+sum(res*res))/rand_chisq(nu_e+real(size(res),dp)); deallocate(res)
      end do
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(:,:,mkeep)=betas; out%deltadraw(mkeep,:,:)=delta; out%taudraw(mkeep,:)=tau
        out%probdraw(mkeep,:)=p; out%zdraw(mkeep,:)=zz; out%mixdraw(mkeep)=mix
      end if
    end do
  end function rhier_linear_mixture

  function rhier_bin_logit(data,z,deltabar,a,nu,v,sbeta,nrep,keep) result(out)
    type(reg_data), intent(in) :: data(:)
    real(dp), intent(in) :: z(:,:),deltabar(:,:),a(:,:),nu,v(:,:),sbeta
    integer, intent(in) :: nrep,keep
    type(hier_mnl_result) :: out
    integer :: nreg,nvar,nz,ns,rep,mkeep,i,j,stat,reject
    integer, allocatable :: yi(:)
    real(dp), allocatable :: betas(:,:),delta(:,:),vbeta(:,:),vbi(:,:),root(:,:),inc(:,:),zn(:),cand(:),meanb(:)
    real(dp) :: oldll,newll,oldprior,newprior,alpha,logl
    type(multireg_draw) :: mr
    nreg=size(data); nvar=size(data(1)%x,2); nz=size(z,2); ns=nrep/keep
    allocate(betas(nreg,nvar),delta(nz,nvar),vbeta(nvar,nvar),vbi(nvar,nvar),root(nvar,nvar),inc(nvar,nvar))
    allocate(zn(nvar),cand(nvar),meanb(nvar)); betas=0.0_dp; delta=deltabar; vbeta=0.0_dp
    do i=1,nvar; vbeta(i,i)=1.0_dp; end do
    allocate(out%betadraw(nreg,nvar,ns),out%deltadraw(ns,nz,nvar),out%vbetadraw(ns,nvar,nvar))
    allocate(out%probdraw(ns,1),out%zdraw(ns,nreg),out%llike(ns),out%reject(ns)); mkeep=0
    do rep=1,nrep
      call inverse_spd(vbeta,vbi,stat); call chol_upper(vbeta,root,stat); inc=transpose(root); reject=0; logl=0.0_dp
      do i=1,nreg
        allocate(yi(size(data(i)%y))); yi=nint(data(i)%y); meanb=matmul(z(i,:),delta)
        oldll=binlogit_loglike(yi,data(i)%x,betas(i,:)); do j=1,nvar; zn(j)=randn(); end do
        cand=betas(i,:)+sbeta*matmul(inc,zn); newll=binlogit_loglike(yi,data(i)%x,cand)
        oldprior=-0.5_dp*dot_product(betas(i,:)-meanb,matmul(vbi,betas(i,:)-meanb))
        newprior=-0.5_dp*dot_product(cand-meanb,matmul(vbi,cand-meanb))
        alpha=min(1.0_dp,exp(min(0.0_dp,newll+newprior-oldll-oldprior)))
        if (rand_uniform()<=alpha) then; betas(i,:)=cand; logl=logl+newll; else; reject=reject+1; logl=logl+oldll; end if
        deallocate(yi)
      end do
      mr=rmultireg(betas,z,deltabar,a,nu,v,stat); delta=mr%b; vbeta=mr%sigma
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(:,:,mkeep)=betas; out%deltadraw(mkeep,:,:)=delta; out%vbetadraw(mkeep,:,:)=vbeta
        out%llike(mkeep)=logl; out%reject(mkeep)=real(reject,dp)/real(nreg,dp); out%probdraw(mkeep,1)=1.0_dp
        out%zdraw(mkeep,:)=1
      end if
    end do
  end function rhier_bin_logit

  function rhier_mnl_rw_mixture(data,z,deltabar,ad,mubar,amu,nu,v,adir,sbeta,nrep,keep,drawdelta) result(out)
    type(mnl_data), intent(in) :: data(:)
    real(dp), intent(in) :: z(:,:),deltabar(:),ad(:,:),mubar(:),amu,nu,v(:,:),adir(:),sbeta
    integer, intent(in) :: nrep,keep
    logical, intent(in), optional :: drawdelta
    type(hier_mnl_result) :: out
    integer :: nreg,nvar,nz,ncomp,ns,rep,mkeep,i,j,k,stat,reject
    logical :: dd
    real(dp), allocatable :: betas(:,:),delta(:,:),resid(:,:),p(:),cand(:),zn(:),root(:,:),meanb(:),bbar(:,:),amat(:,:)
    integer, allocatable :: zz(:)
    real(dp) :: oldll,newll,oldpr,newpr,alpha,logl
    type(nmix_result) :: nm
    type(normal_mixture) :: mix
    dd=.true.; if (present(drawdelta)) dd=drawdelta
    nreg=size(data); nvar=size(data(1)%x,2); nz=size(z,2); ncomp=size(adir); ns=nrep/keep
    allocate(betas(nreg,nvar),delta(nz,nvar),resid(nreg,nvar),p(ncomp),zz(nreg),cand(nvar),zn(nvar),root(nvar,nvar))
    allocate(meanb(nvar),bbar(1,nvar),amat(1,1)); bbar(1,:)=mubar; amat(1,1)=amu
    betas=0.0_dp; call delta_vec_to_matrix(deltabar,nz,nvar,delta); p=1.0_dp/real(ncomp,dp)
    do i=1,nreg; zz(i)=mod(i-1,ncomp)+1; end do
    allocate(out%betadraw(nreg,nvar,ns),out%deltadraw(ns,nz,nvar),out%vbetadraw(ns,nvar,nvar))
    allocate(out%probdraw(ns,ncomp),out%zdraw(ns,nreg),out%llike(ns),out%reject(ns)); out%vbetadraw=0.0_dp; mkeep=0
    do rep=1,nrep
      if (dd) then; resid=betas-matmul(z,delta); else; resid=betas; end if
      nm=rnmix_gibbs(resid,bbar,amat,nu,v,adir,1,1,p,zz); p=nm%probdraw(1,:); zz=nm%zdraw(1,:); mix=nm%mixdraw(1)
      if (dd) delta=draw_delta(z,betas,zz,mix,deltabar,ad,stat)
      reject=0; logl=0.0_dp
      do i=1,nreg
        k=zz(i); call inverse_upper(mix%comp(k)%rooti,root,stat); meanb=mix%comp(k)%mu
        if (dd) meanb=meanb+matmul(z(i,:),delta)
        do j=1,nvar; zn(j)=randn(); end do; cand=betas(i,:)+sbeta*matmul(transpose(root),zn)
        oldll=llmnl(betas(i,:),data(i)%y,data(i)%x); newll=llmnl(cand,data(i)%y,data(i)%x)
        oldpr=lnd_mvn(betas(i,:),meanb,mix%comp(k)%rooti); newpr=lnd_mvn(cand,meanb,mix%comp(k)%rooti)
        alpha=min(1.0_dp,exp(min(0.0_dp,newll+newpr-oldll-oldpr)))
        if (rand_uniform()<=alpha) then; betas(i,:)=cand; logl=logl+newll; else; reject=reject+1; logl=logl+oldll; end if
      end do
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(:,:,mkeep)=betas; out%deltadraw(mkeep,:,:)=delta; out%probdraw(mkeep,:)=p
        out%zdraw(mkeep,:)=zz; out%llike(mkeep)=logl; out%reject(mkeep)=real(reject,dp)/real(nreg,dp)
      end if
    end do
  end function rhier_mnl_rw_mixture

  function rhier_mnl_dp(data,z,deltabar,ad,mubar,amu,nu,v,alpha,maxcomp,sbeta,nrep,keep) result(out)
    type(mnl_data), intent(in) :: data(:)
    real(dp), intent(in) :: z(:,:),deltabar(:),ad(:,:),mubar(:),amu,nu,v(:,:),alpha,sbeta
    integer, intent(in) :: maxcomp,nrep,keep
    type(hier_mnl_result) :: out
    real(dp), allocatable :: adir(:)
    integer :: k
    if (alpha<=0.0_dp .or. maxcomp<1) error stop "rhier_mnl_dp: invalid DP approximation"
    allocate(adir(maxcomp)); adir=alpha/real(maxcomp,dp)
    out=rhier_mnl_rw_mixture(data,z,deltabar,ad,mubar,amu,nu,v,adir,sbeta,nrep,keep,.true.)
    ! This is the standard finite weak-limit representation of a DP prior; as maxcomp grows,
    ! Dirichlet(alpha/maxcomp,...,alpha/maxcomp) converges to the Dirichlet process.
    k=size(adir); if (k<1) error stop "rhier_mnl_dp: internal error"
  end function rhier_mnl_dp
end module bayesm_hierarchical
