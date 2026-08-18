module bayesm_regression
  use bayesm_kinds, only: dp
  use bayesm_linalg, only: chol_upper, inverse_upper, inverse_spd, solve_spd, identity_matrix
  use bayesm_rng, only: randn, rand_chisq, rwishart_draw
  use bayesm_types, only: unireg_result, multireg_draw, reg_data, sur_result
  implicit none
  private
  public :: breg, rmultireg, runireg, runireg_gibbs, rsur_gibbs
contains
  function breg(y,x,betabar,a,info) result(beta)
    real(dp), intent(in) :: y(:),x(:,:),betabar(:),a(:,:)
    integer, intent(out), optional :: info
    real(dp) :: beta(size(betabar)),p(size(a,1),size(a,2)),r(size(a,1),size(a,2))
    real(dp) :: ri(size(a,1),size(a,2)),rhs(size(betabar)),meanb(size(betabar)),z(size(betabar))
    integer :: stat,i
    p=matmul(transpose(x),x)+a
    rhs=matmul(transpose(x),y)+matmul(a,betabar)
    call solve_spd(p,rhs,meanb,stat)
    if (stat==0) call chol_upper(p,r,stat)
    if (stat==0) call inverse_upper(r,ri,stat)
    if (stat==0) then
      do i=1,size(z); z(i)=randn(); end do
      beta=meanb+matmul(ri,z)
    else
      beta=meanb
    end if
    if (present(info)) info=stat
  end function breg

  function rmultireg(y,x,bbar,a,nu,v,info) result(out)
    real(dp), intent(in) :: y(:,:),x(:,:),bbar(:,:),a(:,:),nu,v(:,:)
    integer, intent(out), optional :: info
    type(multireg_draw) :: out
    integer :: n,m,k,stat,i,j
    real(dp), allocatable :: w(:,:),zaug(:,:),p(:,:),r(:,:),ri(:,:),btilde(:,:),e(:,:),s(:,:)
    real(dp), allocatable :: vsinv(:,:),wdraw(:,:),iw(:,:),c(:,:),ci(:,:),zrand(:,:)
    n=size(y,1); m=size(y,2); k=size(x,2)
    allocate(w(n+k,k),zaug(n+k,m),p(k,k),r(k,k),ri(k,k),btilde(k,m),e(n+k,m),s(m,m))
    call chol_upper(a,r,stat)
    if (stat/=0) then
      allocate(out%b(k,m),out%sigma(m,m)); out%b=0.0_dp; out%sigma=0.0_dp
      if (present(info)) info=stat; return
    end if
    w(1:n,:)=x; w(n+1:n+k,:)=r
    zaug(1:n,:)=y; zaug(n+1:n+k,:)=matmul(r,bbar)
    p=matmul(transpose(w),w)
    call chol_upper(p,r,stat); if (stat==0) call inverse_upper(r,ri,stat)
    if (stat/=0) then
      allocate(out%b(k,m),out%sigma(m,m)); out%b=0.0_dp; out%sigma=0.0_dp
      if (present(info)) info=stat; return
    end if
    btilde=matmul(matmul(ri,transpose(ri)),matmul(transpose(w),zaug))
    e=zaug-matmul(w,btilde); s=matmul(transpose(e),e)
    allocate(vsinv(m,m),wdraw(m,m),iw(m,m),c(m,m),ci(m,m),zrand(k,m))
    call inverse_spd(v+s,vsinv,stat)
    if (stat==0) call rwishart_draw(nu+real(n,dp),vsinv,wdraw,iw,c,ci,stat)
    allocate(out%b(k,m),out%sigma(m,m))
    if (stat==0) then
      do j=1,m; do i=1,k; zrand(i,j)=randn(); end do; end do
      out%b=btilde+matmul(matmul(ri,zrand),transpose(ci))
      out%sigma=iw
    else
      out%b=btilde; out%sigma=0.0_dp
    end if
    if (present(info)) info=stat
  end function rmultireg

  function runireg(y,x,betabar,a,nu,ssq,nrep,keep) result(out)
    real(dp), intent(in) :: y(:),x(:,:),betabar(:),a(:,:),nu,ssq
    integer, intent(in) :: nrep,keep
    type(unireg_result) :: out
    integer :: n,k,ns,rep,mkeep,stat,i
    real(dp) :: p(size(a,1),size(a,2)),r(size(a,1),size(a,2)),ri(size(a,1),size(a,2))
    real(dp) :: rhs(size(betabar)),btilde(size(betabar)),beta(size(betabar)),z(size(betabar))
    real(dp) :: sigmasq,s
    n=size(y); k=size(x,2); ns=nrep/keep
    allocate(out%betadraw(ns,k),out%sigmasqdraw(ns))
    p=matmul(transpose(x),x)+a
    rhs=matmul(transpose(x),y)+matmul(a,betabar)
    call solve_spd(p,rhs,btilde,stat); call chol_upper(p,r,stat); call inverse_upper(r,ri,stat)
    s=sum((y-matmul(x,btilde))**2)+dot_product(btilde-betabar,matmul(a,btilde-betabar))
    mkeep=0
    do rep=1,nrep
      sigmasq=(nu*ssq+s)/rand_chisq(nu+real(n,dp))
      do i=1,k; z(i)=randn(); end do
      beta=btilde+sqrt(sigmasq)*matmul(ri,z)
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%sigmasqdraw(mkeep)=sigmasq
      end if
    end do
  end function runireg

  function runireg_gibbs(y,x,betabar,a,nu,ssq,sigmasq0,nrep,keep) result(out)
    real(dp), intent(in) :: y(:),x(:,:),betabar(:),a(:,:),nu,ssq,sigmasq0
    integer, intent(in) :: nrep,keep
    type(unireg_result) :: out
    integer :: n,k,ns,rep,mkeep,stat,i
    real(dp) :: p(size(a,1),size(a,2)),r(size(a,1),size(a,2)),ri(size(a,1),size(a,2))
    real(dp) :: rhs(size(betabar)),btilde(size(betabar)),beta(size(betabar)),z(size(betabar))
    real(dp) :: sigmasq,s
    n=size(y); k=size(x,2); ns=nrep/keep; sigmasq=sigmasq0
    allocate(out%betadraw(ns,k),out%sigmasqdraw(ns)); mkeep=0
    do rep=1,nrep
      p=matmul(transpose(x),x)/sigmasq+a
      rhs=matmul(transpose(x),y)/sigmasq+matmul(a,betabar)
      call solve_spd(p,rhs,btilde,stat); call chol_upper(p,r,stat); call inverse_upper(r,ri,stat)
      do i=1,k; z(i)=randn(); end do
      beta=btilde+matmul(ri,z)
      s=sum((y-matmul(x,beta))**2)
      sigmasq=(nu*ssq+s)/rand_chisq(nu+real(n,dp))
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%sigmasqdraw(mkeep)=sigmasq
      end if
    end do
  end function runireg_gibbs

  function rsur_gibbs(data,betabar,a,nu,v,sigma0,nrep,keep) result(out)
    type(reg_data), intent(in) :: data(:)
    real(dp), intent(in) :: betabar(:),a(:,:),nu,v(:,:),sigma0(:,:)
    integer, intent(in) :: nrep,keep
    type(sur_result) :: out
    integer :: neq,nobs,nvar,ns,rep,mkeep,i,j,ii,jj,ki,kj,offi,offj,stat,q
    integer, allocatable :: nk(:),off(:)
    real(dp), allocatable :: siginv(:,:),prec(:,:),rhs(:),meanb(:),r(:,:),ri(:,:),beta(:),zn(:)
    real(dp), allocatable :: e(:,:),s(:,:),vsinv(:,:),w(:,:),iw(:,:),c(:,:),ci(:,:),yd(:,:)
    neq=size(data); nobs=size(data(1)%y); nvar=size(betabar); ns=nrep/keep
    allocate(nk(neq),off(neq+1)); off(1)=1
    do i=1,neq; nk(i)=size(data(i)%x,2); off(i+1)=off(i)+nk(i); end do
    allocate(out%betadraw(ns,nvar),out%sigmadraw(ns,neq,neq))
    allocate(siginv(neq,neq),prec(nvar,nvar),rhs(nvar),meanb(nvar),r(nvar,nvar),ri(nvar,nvar))
    allocate(beta(nvar),zn(nvar),e(nobs,neq),s(neq,neq),vsinv(neq,neq),w(neq,neq),iw(neq,neq),c(neq,neq),ci(neq,neq))
    allocate(yd(nobs,neq)); siginv=0.0_dp
    call inverse_spd(sigma0,siginv,stat); mkeep=0
    do rep=1,nrep
      prec=a; rhs=matmul(a,betabar)
      do i=1,neq
        offi=off(i); ki=nk(i)
        do j=1,neq
          offj=off(j); kj=nk(j)
          do ii=1,ki
            do jj=1,kj
              prec(offi+ii-1,offj+jj-1)=prec(offi+ii-1,offj+jj-1)+siginv(i,j)* &
                dot_product(data(i)%x(:,ii),data(j)%x(:,jj))
            end do
          end do
          rhs(offi:offi+ki-1)=rhs(offi:offi+ki-1)+siginv(i,j)*matmul(transpose(data(i)%x),data(j)%y)
        end do
      end do
      call solve_spd(prec,rhs,meanb,stat); call chol_upper(prec,r,stat); call inverse_upper(r,ri,stat)
      do q=1,nvar; zn(q)=randn(); end do
      beta=meanb+matmul(ri,zn)
      do i=1,neq
        e(:,i)=data(i)%y-matmul(data(i)%x,beta(off(i):off(i+1)-1))
      end do
      s=matmul(transpose(e),e)+v
      call inverse_spd(s,vsinv,stat)
      call rwishart_draw(nu+real(nobs,dp),vsinv,w,iw,c,ci,stat)
      siginv=w
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%sigmadraw(mkeep,:,:)=iw
      end if
    end do
  end function rsur_gibbs
end module bayesm_regression
