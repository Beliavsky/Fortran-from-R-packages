module bayesm_dp
  use bayesm_kinds, only: dp
  use bayesm_linalg, only: chol_upper, inverse_upper, inverse_spd
  use bayesm_rng, only: rand_uniform, rand_beta, rand_gamma, rand_categorical, randn, rwishart_draw
  use bayesm_densities, only: lnd_mvn
  use bayesm_types, only: dp_mixture_result, normal_component
  implicit none
  private
  public :: rdp_gibbs
contains
  function rdp_gibbs(y,m0,kappa0,nu0,s0,alpha0,nrep,keep,alpha_a,alpha_b) result(out)
    real(dp), intent(in) :: y(:,:),m0(:),kappa0,nu0,s0(:,:),alpha0
    integer, intent(in) :: nrep,keep
    real(dp), intent(in), optional :: alpha_a,alpha_b
    type(dp_mixture_result) :: out
    integer :: n,d,ns,rep,mkeep,i,k,kmax,nk,stat,j,oldk,newk
    integer, allocatable :: z(:),counts(:)
    real(dp), allocatable :: mu(:,:),rooti(:,:,:),logp(:),prob(:),mean(:),scatter(:,:),vpost(:,:),vinv(:,:)
    real(dp), allocatable :: w(:,:),iw(:,:),c(:,:),ci(:,:),r(:,:),rwork(:,:),zn(:),predcov(:,:),predr(:,:),predri(:,:)
    real(dp) :: alpha,eta,aa,bb,pi_eta
    n=size(y,1); d=size(y,2); ns=nrep/keep
    if (size(m0)/=d .or. size(s0,1)/=d .or. size(s0,2)/=d) error stop "rdp_gibbs: dimension mismatch"
    allocate(z(n),counts(n),mu(n,d),rooti(n,d,d),logp(n+1),prob(n+1),mean(d),scatter(d,d))
    allocate(vpost(d,d),vinv(d,d),w(d,d),iw(d,d),c(d,d),ci(d,d),r(d,d),rwork(d,d),zn(d))
    allocate(predcov(d,d),predr(d,d),predri(d,d)); z=1; counts=0; counts(1)=n; kmax=1; alpha=alpha0
    mu=0.0_dp; rooti=0.0_dp
    aa=1.0_dp; bb=1.0_dp
    if (present(alpha_a)) aa=alpha_a
    if (present(alpha_b)) bb=alpha_b
    allocate(out%alphadraw(ns),out%zdraw(ns,n),out%ncompdraw(ns),out%mixdraw(ns)); mkeep=0
    do rep=1,nrep
      ! Draw component parameters from the normal-inverse-Wishart posterior.
      do k=1,kmax
        nk=count(z==k)
        if (nk==0) cycle
        mean=0.0_dp
        do i=1,n
          if (z(i)==k) mean=mean+y(i,:)
        end do
        mean=mean/real(nk,dp); scatter=0.0_dp
        do i=1,n
          if (z(i)==k) then
            scatter=scatter+spread(y(i,:)-mean,2,d)*spread(y(i,:)-mean,1,d)
          end if
        end do
        vpost=s0+scatter+(kappa0*real(nk,dp)/(kappa0+real(nk,dp)))* &
          spread(mean-m0,2,d)*spread(mean-m0,1,d)
        call inverse_spd(vpost,vinv,stat)
        call rwishart_draw(nu0+real(nk,dp),vinv,w,iw,c,ci,stat)
        call chol_upper(iw,r,stat); call inverse_upper(r,rwork,stat); rooti(k,:,:)=rwork
        mean=(kappa0*m0+real(nk,dp)*mean)/(kappa0+real(nk,dp))
        do j=1,d; zn(j)=randn(); end do
        mu(k,:)=mean+matmul(transpose(r),zn)/sqrt(kappa0+real(nk,dp))
      end do
      ! Collapsed-style indicator sweep with a fresh base component option.
      do i=1,n
        oldk=z(i); counts=0
        do j=1,n
          if (j/=i) counts(z(j))=counts(z(j))+1
        end do
        do k=1,kmax
          if (counts(k)>0) then
            logp(k)=log(real(counts(k),dp))+lnd_mvn(y(i,:),mu(k,:),rooti(k,:,:))
          else
            logp(k)=-huge(1.0_dp)
          end if
        end do
        predcov=s0*(1.0_dp+1.0_dp/kappa0)/max(nu0-real(d,dp)-1.0_dp,1.0_dp)
        call chol_upper(predcov,predr,stat); call inverse_upper(predr,predri,stat)
        logp(kmax+1)=log(alpha)+lnd_mvn(y(i,:),m0,predri)
        logp(1:kmax+1)=logp(1:kmax+1)-maxval(logp(1:kmax+1))
        prob(1:kmax+1)=exp(logp(1:kmax+1)); newk=rand_categorical(prob(1:kmax+1))
        if (newk==kmax+1) then
          kmax=kmax+1; newk=kmax
          call inverse_spd(s0,vinv,stat); call rwishart_draw(nu0,vinv,w,iw,c,ci,stat)
          call chol_upper(iw,r,stat); call inverse_upper(r,rwork,stat); rooti(newk,:,:)=rwork
          do j=1,d; zn(j)=randn(); end do
          mu(newk,:)=m0+matmul(transpose(r),zn)/sqrt(kappa0)
        end if
        z(i)=newk
        if (oldk<=kmax .and. count(z==oldk)==0 .and. oldk<kmax) then
          do j=1,n
            if (z(j)==kmax) z(j)=oldk
          end do
          mu(oldk,:)=mu(kmax,:); rooti(oldk,:,:)=rooti(kmax,:,:); kmax=kmax-1
        else if (oldk==kmax .and. count(z==oldk)==0) then
          kmax=kmax-1
        end if
      end do
      ! Escobar-West concentration update under Gamma(aa,bb) prior.
      eta=rand_beta(alpha+1.0_dp,real(n,dp))
      pi_eta=(aa+real(kmax,dp)-1.0_dp)/(real(n,dp)*(bb-log(eta))+aa+real(kmax,dp)-1.0_dp)
      if (rand_uniform()<pi_eta) then
        alpha=rand_gamma(aa+real(kmax,dp),1.0_dp/(bb-log(eta)))
      else
        alpha=rand_gamma(aa+real(kmax,dp)-1.0_dp,1.0_dp/(bb-log(eta)))
      end if
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%alphadraw(mkeep)=alpha; out%zdraw(mkeep,:)=z; out%ncompdraw(mkeep)=kmax
        allocate(out%mixdraw(mkeep)%p(kmax),out%mixdraw(mkeep)%comp(kmax))
        do k=1,kmax
          out%mixdraw(mkeep)%p(k)=real(count(z==k),dp)/real(n,dp)
          allocate(out%mixdraw(mkeep)%comp(k)%mu(d),out%mixdraw(mkeep)%comp(k)%rooti(d,d))
          out%mixdraw(mkeep)%comp(k)%mu=mu(k,:); out%mixdraw(mkeep)%comp(k)%rooti=rooti(k,:,:)
        end do
      end if
    end do
  end function rdp_gibbs
end module bayesm_dp
