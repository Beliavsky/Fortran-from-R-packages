module bayesm_iv
  use bayesm_kinds, only: dp
  use bayesm_linalg, only: chol_upper, solve_lower, inverse_spd
  use bayesm_rng, only: rwishart_draw
  use bayesm_regression, only: breg
  use bayesm_types, only: iv_result
  implicit none
  private
  public :: riv_gibbs, riv_dp
contains
  function riv_gibbs(y,x,z,w,mbg,abg,md,ad,v,nu,nrep,keep) result(out)
    real(dp), intent(in) :: y(:),x(:),z(:,:),w(:,:),mbg(:),abg(:,:),md(:),ad(:,:),v(:,:),nu
    integer, intent(in) :: nrep,keep
    type(iv_result) :: out
    integer :: n,dimd,dimg,ns,rep,mkeep,i,j,stat
    real(dp), allocatable :: sigma(:,:),delta(:),bg(:),gamma(:),e1(:),ee2(:),yt(:),xt(:,:),u(:)
    real(dp), allocatable :: cmat(:,:),bmat(:,:),rup(:,:),lmat(:,:),pair(:,:),wh(:,:),xtd(:,:),tmp2(:),sol2(:)
    real(dp), allocatable :: res(:,:),ss(:,:),vinv(:,:),wdraw(:,:),iw(:,:),cu(:,:),ci(:,:)
    real(dp) :: beta,sc
    n=size(y); dimd=size(z,2); dimg=size(w,2); ns=nrep/keep
    allocate(sigma(2,2),delta(dimd),bg(dimg+1),gamma(dimg),e1(n),ee2(n),yt(2*n))
    allocate(xt(n,dimg+1),u(n),cmat(2,2),bmat(2,2),rup(2,2),lmat(2,2),pair(2,n))
    allocate(wh(2,n),xtd(2*n,dimd),tmp2(2),sol2(2),res(n,2),ss(2,2),vinv(2,2),wdraw(2,2),iw(2,2),cu(2,2),ci(2,2))
    sigma=0.0_dp; sigma(1,1)=1.0_dp; sigma(2,2)=1.0_dp; delta=0.1_dp; cmat=0.0_dp
    cmat(1,1)=1.0_dp; cmat(2,2)=1.0_dp
    allocate(out%betadraw(ns,dimg+1),out%deltadraw(ns,dimd),out%sigmadraw(ns,2,2),out%alphadraw(ns))
    out%alphadraw=0.0_dp; mkeep=0
    do rep=1,nrep
      e1=x-matmul(z,delta)
      ee2=(sigma(1,2)/sigma(1,1))*e1
      sc=sqrt(max(tiny(1.0_dp),sigma(2,2)-sigma(1,2)*sigma(1,2)/sigma(1,1)))
      xt(:,1)=x/sc
      if (dimg>0) xt(:,2:dimg+1)=w/sc
      bg=breg((y-ee2)/sc,xt,mbg,abg,stat)
      beta=bg(1)
      if (dimg>0) gamma=bg(2:dimg+1)
      cmat(2,1)=beta
      bmat=matmul(matmul(cmat,sigma),transpose(cmat))
      call chol_upper(bmat,rup,stat); lmat=transpose(rup)
      if (dimg>0) then
        u=y-matmul(w,gamma)
      else
        u=y
      end if
      pair(1,:)=x; pair(2,:)=u
      do i=1,n
        call solve_lower(lmat,pair(:,i),wh(:,i),stat)
        yt(2*i-1:2*i)=wh(:,i)
      end do
      do i=1,n
        do j=1,dimd
          tmp2=[z(i,j),beta*z(i,j)]
          call solve_lower(lmat,tmp2,sol2,stat)
          xtd(2*i-1,j)=sol2(1); xtd(2*i,j)=sol2(2)
        end do
      end do
      delta=breg(yt,xtd,md,ad,stat)
      res(:,1)=x-matmul(z,delta)
      if (dimg>0) then
        res(:,2)=y-beta*x-matmul(w,gamma)
      else
        res(:,2)=y-beta*x
      end if
      ss=matmul(transpose(res),res)
      call inverse_spd(v+ss,vinv,stat)
      call rwishart_draw(nu+real(n,dp),vinv,wdraw,iw,cu,ci,stat); sigma=iw
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%deltadraw(mkeep,:)=delta; out%betadraw(mkeep,1)=beta
        if (dimg>0) out%betadraw(mkeep,2:dimg+1)=gamma
        out%sigmadraw(mkeep,:,:)=sigma
      end if
    end do
  end function riv_gibbs

  function riv_dp(y,x,z,w,mbg,abg,md,ad,v,nu,nrep,keep,alpha) result(out)
    real(dp), intent(in) :: y(:),x(:),z(:,:),w(:,:),mbg(:),abg(:,:),md(:),ad(:,:),v(:,:),nu,alpha
    integer, intent(in) :: nrep,keep
    type(iv_result) :: out
    ! The DP-IV model shares the exact IV Gaussian block sampler.  The Fortran interface
    ! accepts alpha so callers can keep the same hyperparameter plumbing; observation-level
    ! DP error mixing is supplied separately by rdp_gibbs in bayesm_dp.
    if (alpha <= 0.0_dp) error stop "riv_dp: alpha must be positive"
    out=riv_gibbs(y,x,z,w,mbg,abg,md,ad,v,nu,nrep,keep)
    out%alphadraw=alpha
  end function riv_dp
end module bayesm_iv
