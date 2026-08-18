module bayesm_negbin
  use bayesm_kinds, only: dp
  use bayesm_linalg, only: chol_upper, inverse_spd, matrix_sqrt_spd
  use bayesm_rng, only: randn, rand_uniform
  use bayesm_densities, only: lnd_mvn
  use bayesm_regression, only: rmultireg
  use bayesm_types, only: negbin_result, hier_negbin_result, reg_data, multireg_draw
  implicit none
  private
  public :: llnegbin, rnegbin_rw, rhier_negbin_rw
contains
  pure real(dp) function llnegbin(y,lambda,alpha,constant) result(ll)
    real(dp), intent(in) :: y(:),lambda(:),alpha
    logical, intent(in), optional :: constant
    logical :: cn
    integer :: i
    real(dp) :: prob
    cn=.true.; if (present(constant)) cn=constant; ll=0.0_dp
    do i=1,size(y)
      prob=alpha/(alpha+lambda(i))
      prob=min(1.0_dp-epsilon(1.0_dp),max(tiny(1.0_dp),prob))
      if (cn) then
        ll=ll+log_gamma(y(i)+alpha)-log_gamma(alpha)-log_gamma(y(i)+1.0_dp) &
          +alpha*log(prob)+y(i)*log(1.0_dp-prob)
      else
        ll=ll+alpha*log(prob)+y(i)*log(1.0_dp-prob)
      end if
    end do
  end function llnegbin

  real(dp) function lpost_beta(alpha,beta,x,y,betabar,roota) result(v)
    real(dp), intent(in) :: alpha,beta(:),x(:,:),y(:),betabar(:),roota(:,:)
    real(dp), allocatable :: lambda(:)
    real(dp) :: z(size(beta))
    allocate(lambda(size(y)))
    lambda=exp(min(700.0_dp,matmul(x,beta))); z=matmul(roota,beta-betabar)
    v=llnegbin(y,lambda,alpha,.false.)-0.5_dp*dot_product(z,z)
  end function lpost_beta

  real(dp) function lpost_alpha(alpha,beta,x,y,a,b) result(v)
    real(dp), intent(in) :: alpha,beta(:),x(:,:),y(:),a,b
    real(dp), allocatable :: lambda(:)
    allocate(lambda(size(y)))
    lambda=exp(min(700.0_dp,matmul(x,beta)))
    v=llnegbin(y,lambda,alpha,.true.)+(a-1.0_dp)*log(alpha)-b*alpha
  end function lpost_alpha

  function rnegbin_rw(y,x,betabar,a_prec,a_shape,b_rate,beta0,alpha0,betaroot,alphacroot,nrep,keep,fixalpha) result(out)
    real(dp), intent(in) :: y(:),x(:,:),betabar(:),a_prec(:,:),a_shape,b_rate,beta0(:),alpha0,betaroot(:,:)
    real(dp), intent(in) :: alphacroot
    integer, intent(in) :: nrep,keep
    logical, intent(in), optional :: fixalpha
    type(negbin_result) :: out
    integer :: ns,rep,mkeep,i,k,stat,naccb,nacca
    real(dp) :: beta(size(beta0)),bc(size(beta0)),alpha,ac,roota(size(a_prec,1),size(a_prec,2)),zn(size(beta0))
    real(dp) :: oldp,newp,acc
    logical :: fixa
    fixa=.false.; if (present(fixalpha)) fixa=fixalpha
    ns=nrep/keep; k=size(beta0); allocate(out%betadraw(ns,k),out%alphadraw(ns),out%llike(ns))
    beta=beta0; alpha=alpha0; call chol_upper(a_prec,roota,stat); mkeep=0; naccb=0; nacca=0
    do rep=1,nrep
      do i=1,k; zn(i)=randn(); end do; bc=beta+matmul(betaroot,zn)
      oldp=lpost_beta(alpha,beta,x,y,betabar,roota); newp=lpost_beta(alpha,bc,x,y,betabar,roota)
      acc=min(1.0_dp,exp(min(0.0_dp,newp-oldp)))
      if (rand_uniform()<=acc) then; beta=bc; naccb=naccb+1; end if
      if (.not.fixa) then
        ac=exp(log(alpha)+alphacroot*randn())
        oldp=lpost_alpha(alpha,beta,x,y,a_shape,b_rate); newp=lpost_alpha(ac,beta,x,y,a_shape,b_rate)
        acc=min(1.0_dp,exp(min(0.0_dp,newp-oldp)))
        if (rand_uniform()<=acc) then; alpha=ac; nacca=nacca+1; end if
      end if
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%alphadraw(mkeep)=alpha
        out%llike(mkeep)=llnegbin(y,exp(min(700.0_dp,matmul(x,beta))),alpha,.true.)
      end if
    end do
    out%beta_accept=real(naccb,dp)/real(nrep,dp); out%alpha_accept=real(nacca,dp)/real(nrep,dp)
  end function rnegbin_rw

  function rhier_negbin_rw(data,z,deltabar,adelta,nu,v,a_shape,b_rate,nrep,keep,sbeta,salpha, &
    beta0,delta0,vbeta0,alpha0,fixalpha) result(out)
    type(reg_data), intent(in) :: data(:)
    real(dp), intent(in) :: z(:,:),deltabar(:,:),adelta(:,:),nu,v(:,:),a_shape,b_rate,sbeta,salpha
    integer, intent(in) :: nrep,keep
    real(dp), intent(in), optional :: beta0(:,:),delta0(:,:),vbeta0(:,:),alpha0
    logical, intent(in), optional :: fixalpha
    type(hier_negbin_result) :: out
    integer :: nreg,nvar,nz,ns,rep,mkeep,i,j,stat,naccb,nacca
    real(dp), allocatable :: beta(:,:),delta(:,:),vbeta(:,:),vbi(:,:),roota(:,:),betabar(:,:),bc(:),zn(:)
    real(dp), allocatable :: hess(:,:),prop_cov(:,:),prop_root(:,:),lambda(:)
    real(dp) :: alpha,ac,oldp,newp,acc,pooled_old,pooled_new
    type(multireg_draw) :: mr
    logical :: fixa
    nreg=size(data); nvar=size(data(1)%x,2); nz=size(z,2); ns=nrep/keep
    allocate(beta(nreg,nvar),delta(nz,nvar),vbeta(nvar,nvar),vbi(nvar,nvar),roota(nvar,nvar),betabar(nreg,nvar))
    allocate(bc(nvar),zn(nvar),hess(nvar,nvar),prop_cov(nvar,nvar),prop_root(nvar,nvar))
    beta=0.0_dp; if(present(beta0)) beta=beta0; delta=0.0_dp; if(present(delta0)) delta=delta0
    vbeta=0.0_dp; do i=1,nvar; vbeta(i,i)=1.0_dp; end do; if(present(vbeta0)) vbeta=vbeta0
    alpha=1.0_dp; if(present(alpha0)) alpha=alpha0; fixa=.false.; if(present(fixalpha)) fixa=fixalpha
    allocate(out%betadraw(ns,nreg,nvar),out%alphadraw(ns),out%vbetadraw(ns,nvar,nvar),out%deltadraw(ns,nz,nvar),out%llike(ns))
    mkeep=0; naccb=0; nacca=0
    do rep=1,nrep
      call inverse_spd(vbeta,vbi,stat); call chol_upper(vbi,roota,stat); betabar=matmul(z,delta)
      do i=1,nreg
        allocate(lambda(size(data(i)%y))); lambda=exp(min(700.0_dp,matmul(data(i)%x,beta(i,:))))
        hess=0.0_dp
        do j=1,size(data(i)%y)
          hess=hess+alpha*lambda(j)/(alpha+lambda(j))*spread(data(i)%x(j,:),2,nvar)*spread(data(i)%x(j,:),1,nvar)
        end do
        deallocate(lambda)
        call inverse_spd(hess+vbi,prop_cov,stat); prop_cov=sbeta*prop_cov
        call matrix_sqrt_spd(prop_cov,prop_root,stat)
        do j=1,nvar; zn(j)=randn(); end do; bc=beta(i,:)+matmul(prop_root,zn)
        oldp=lpost_beta(alpha,beta(i,:),data(i)%x,data(i)%y,betabar(i,:),roota)
        newp=lpost_beta(alpha,bc,data(i)%x,data(i)%y,betabar(i,:),roota)
        acc=min(1.0_dp,exp(min(0.0_dp,newp-oldp)))
        if(rand_uniform()<=acc) then; beta(i,:)=bc; naccb=naccb+1; end if
      end do
      if(.not.fixa) then
        ac=exp(log(alpha)+salpha*randn()); pooled_old=0.0_dp; pooled_new=0.0_dp
        do i=1,nreg
          pooled_old=pooled_old+lpost_alpha(alpha,beta(i,:),data(i)%x,data(i)%y,a_shape,b_rate)/real(nreg,dp)
          pooled_new=pooled_new+lpost_alpha(ac,beta(i,:),data(i)%x,data(i)%y,a_shape,b_rate)/real(nreg,dp)
        end do
        acc=min(1.0_dp,exp(min(0.0_dp,pooled_new-pooled_old)))
        if(rand_uniform()<=acc) then; alpha=ac; nacca=nacca+1; end if
      end if
      mr=rmultireg(beta,z,deltabar,adelta,nu,v,stat); delta=mr%b; vbeta=mr%sigma
      if(mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(mkeep,:,:)=beta; out%alphadraw(mkeep)=alpha
        out%vbetadraw(mkeep,:,:)=vbeta; out%deltadraw(mkeep,:,:)=delta; out%llike(mkeep)=0.0_dp
        do i=1,nreg
          out%llike(mkeep)=out%llike(mkeep)+llnegbin(data(i)%y,exp(min(700.0_dp,matmul(data(i)%x,beta(i,:)))),alpha,.true.)
        end do
      end if
    end do
    out%beta_accept=real(naccb,dp)/real(nrep*nreg,dp); out%alpha_accept=real(nacca,dp)/real(nrep,dp)
  end function rhier_negbin_rw
end module bayesm_negbin
