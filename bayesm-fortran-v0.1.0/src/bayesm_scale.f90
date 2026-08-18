module bayesm_scale
  use bayesm_kinds, only: dp
  use bayesm_math, only: normal_cdf, normal_quantile
  use bayesm_linalg, only: inverse_spd, chol_upper, inverse_upper, solve_spd
  use bayesm_rng, only: rand_uniform, randn, rand_categorical, rwishart_draw
  use bayesm_utils, only: cget_c
  use bayesm_types, only: scale_usage_result
  implicit none
  private
  public :: ghk_region, rscale_usage
contains
  real(dp) function ghk_region(l,a,b,ndraw) result(prob)
    real(dp), intent(in) :: l(:,:),a(:),b(:)
    integer, intent(in) :: ndraw
    integer :: h,j
    real(dp) :: z(size(a)),mu,pa,pb,u,prod
    prob=0.0_dp
    do h=1,ndraw
      prod=1.0_dp; z=0.0_dp
      do j=1,size(a)
        mu=0.0_dp
        if (j>1) mu=dot_product(l(j,1:j-1),z(1:j-1))
        pa=normal_cdf((a(j)-mu)/l(j,j)); pb=normal_cdf((b(j)-mu)/l(j,j))
        prod=prod*max(0.0_dp,pb-pa); u=rand_uniform()
        z(j)=normal_quantile(min(0.999999999_dp,max(1.0e-10_dp,pa+u*(pb-pa))))
      end do
      prob=prob+prod
    end do
    prob=prob/real(ndraw,dp)
  end function ghk_region

  real(dp) function ll_lambda(lam,n,s,v,nu) result(ll)
    integer, intent(in) :: n
    real(dp), intent(in) :: lam(2,2),s(2,2),v(2,2),nu
    real(dp) :: li(2,2),det
    integer :: stat
    det=lam(1,1)*lam(2,2)-lam(1,2)*lam(1,2)
    if (det<=0.0_dp) then; ll=-huge(1.0_dp); return; end if
    call inverse_spd(lam,li,stat)
    if (stat/=0) then; ll=-huge(1.0_dp); return; end if
    ll=-0.5_dp*(real(n,dp)+nu+3.0_dp)*log(det)-0.5_dp*( &
      (s(1,1)+v(1,1))*li(1,1)+(s(1,2)+v(1,2))*li(2,1)+ &
      (s(2,1)+v(2,1))*li(1,2)+(s(2,2)+v(2,2))*li(2,2))
  end function ll_lambda

  real(dp) function scale_logprob(x,k,e,mu,tau,sigma,sigmat,ndraw) result(ll)
    integer, intent(in) :: x(:,:),k,ndraw
    real(dp), intent(in) :: e,mu(:),tau(:),sigma(:),sigmat(:,:)
    integer :: n,p,i,j,stat
    real(dp), allocatable :: c(:),l(:,:),a(:),b(:)
    real(dp) :: pr
    n=size(x,1); p=size(x,2); allocate(c(k+1),l(p,p),a(p),b(p)); call cget_c(e,k,c,stat)
    call chol_upper(sigmat,l,stat); l=transpose(l); ll=0.0_dp
    do i=1,n
      do j=1,p
        a(j)=c(x(i,j))-mu(j)-tau(i); b(j)=c(x(i,j)+1)-mu(j)-tau(i)
      end do
      pr=ghk_region(sigma(i)*l,a,b,ndraw)
      ll=ll+log(max(pr,tiny(1.0_dp)))+real(p,dp)*log(real(k,dp))
    end do
  end function scale_logprob

  function rscale_usage(x,k,nu,v,mubar,am,gsigma,gl11,gl22,gl12,nul,vl,ge,nrep,keep,ndghk) result(out)
    integer, intent(in) :: x(:,:),k,nul,nrep,keep,ndghk
    real(dp), intent(in) :: nu,v(:,:),mubar(:),am(:,:),gsigma(:),gl11(:),gl22(:),gl12(:),vl(:,:),ge(:)
    type(scale_usage_result) :: out
    integer :: n,p,ns,rep,mkeep,i,j,q,stat,idx,ng
    real(dp), allocatable :: y(:,:),mu(:),sigmat(:,:),tau(:),sigma(:),lam(:,:),c(:),si(:,:),res(:,:),ss(:,:)
    real(dp), allocatable :: vinv(:,:),wd(:,:),iw(:,:),cc(:,:),ci(:,:),prec(:,:),rhs(:),mean(:),r(:,:),ri(:,:),zn(:)
    real(dp), allocatable :: pv(:),h(:),dat(:,:),center(:,:),s2(:,:),cand(:),eps(:),one(:),temp(:)
    real(dp) :: ai,ap,xtx,xty,bcoef,varm,mm,aa,bb,sc,oldll,newll,ratio,ecur
    n=size(x,1); p=size(x,2); ns=nrep/keep
    allocate(y(n,p),mu(p),sigmat(p,p),tau(n),sigma(n),lam(2,2),c(k+1),si(p,p),res(n,p),ss(p,p))
    allocate(vinv(p,p),wd(p,p),iw(p,p),cc(p,p),ci(p,p),prec(p,p),rhs(p),mean(p),r(p,p),ri(p,p),zn(p))
    allocate(h(n),dat(n,2),center(n,2),s2(2,2),eps(p),one(p),temp(p)); one=1.0_dp
    y=real(x,dp); mu=sum(y,dim=1)/real(n,dp); tau=0.0_dp; sigma=1.0_dp; sigmat=0.0_dp
    do j=1,p; sigmat(j,j)=1.0_dp; end do
    lam=0.0_dp; lam(1,1)=4.0_dp; lam(2,2)=0.5_dp; ecur=ge(minloc(abs(ge),dim=1))
    allocate(out%mudraw(ns,p),out%sigmadraw(ns,p,p),out%taudraw(ns,n),out%sdraw(ns,n))
    allocate(out%lambdadraw(ns,2,2),out%edraw(ns)); mkeep=0
    do rep=1,nrep
      call cget_c(ecur,k,c,stat)
      call inverse_spd(sigmat,si,stat)
      ! Latent ordinal-normal utilities.
      do i=1,n
        do j=1,p
          mm=mu(j)+tau(i)
          do q=1,p
            if (q/=j) mm=mm-si(j,q)/si(j,j)*(y(i,q)-mu(q)-tau(i))
          end do
          sc=sigma(i)/sqrt(si(j,j))
          aa=normal_cdf((c(x(i,j))-mm)/sc); bb=normal_cdf((c(x(i,j)+1)-mm)/sc)
          y(i,j)=mm+sc*normal_quantile(min(0.999999999_dp,max(1.0e-10_dp,aa+rand_uniform()*(bb-aa))))
        end do
      end do
      do i=1,n
        res(i,:)=(y(i,:)-mu-tau(i))/sigma(i)
      end do
      ss=matmul(transpose(res),res); call inverse_spd(v+ss,vinv,stat)
      call rwishart_draw(nu+real(n,dp),vinv,wd,iw,cc,ci,stat); sigmat=iw; call inverse_spd(sigmat,si,stat)
      prec=sum(1.0_dp/(sigma*sigma))*si+am; rhs=matmul(am,mubar)
      do i=1,n; rhs=rhs+matmul(si,(y(i,:)-tau(i)))/(sigma(i)*sigma(i)); end do
      call solve_spd(prec,rhs,mean,stat); call chol_upper(prec,r,stat); call inverse_upper(r,ri,stat)
      do j=1,p; zn(j)=randn(); end do; mu=mean+matmul(ri,zn)
      ai=lam(1,1)-lam(1,2)*lam(1,2)/lam(2,2); ap=1.0_dp/ai
      bcoef=ap*lam(1,2)/lam(2,2); xtx=dot_product(one,matmul(si,one))
      do i=1,n
        temp=y(i,:)-mu; xty=dot_product(one,matmul(si,temp))
        varm=1.0_dp/(xtx/(sigma(i)*sigma(i))+ap)
        mm=varm*(xty/(sigma(i)*sigma(i))+bcoef*(log(sigma(i))-lam(2,2)))
        tau(i)=mm+sqrt(varm)*randn()
      end do
      do i=1,n
        temp=(y(i,:)-mu-tau(i)); eps=matmul(si,temp); aa=dot_product(temp,eps)
        allocate(pv(size(gsigma)))
        bb=lam(1,2)/lam(1,1); sc=sqrt(max(tiny(1.0_dp),lam(2,2)-lam(1,2)*lam(1,2)/lam(1,1)))
        do j=1,size(gsigma)
          pv(j)=-(real(p,dp)+1.0_dp)*log(gsigma(j))-0.5_dp*aa/(gsigma(j)*gsigma(j))- &
            0.5_dp*((log(gsigma(j))-(lam(2,2)+bb*tau(i)))/sc)**2
        end do
        pv=exp(pv-maxval(pv)); idx=rand_categorical(pv); sigma(i)=gsigma(idx); deallocate(pv)
      end do
      h=log(sigma); dat(:,1)=tau; dat(:,2)=h; center=dat
      center(:,1)=center(:,1)-sum(tau)/real(n,dp); center(:,2)=center(:,2)-sum(h)/real(n,dp)
      s2=matmul(transpose(center),center)
      ! Grid Gibbs for Lambda, preserving the upstream conditional-grid construction.
      allocate(cand(size(gl11)),pv(size(gl11))); ng=0
      do j=1,size(gl11)
        if (gl11(j)>lam(1,2)*lam(1,2)/lam(2,2)) then
          ng=ng+1; cand(ng)=gl11(j); lam(1,1)=cand(ng); pv(ng)=ll_lambda(lam,n,s2,vl,real(nul,dp))
        end if
      end do
      if (ng>0) then; pv(1:ng)=exp(pv(1:ng)-maxval(pv(1:ng))); lam(1,1)=cand(rand_categorical(pv(1:ng))); end if
      deallocate(cand,pv); allocate(cand(size(gl12)),pv(size(gl12))); ng=0
      do j=1,size(gl12)
        if (abs(gl12(j))<sqrt(lam(1,1)*lam(2,2))) then
          ng=ng+1; cand(ng)=gl12(j); lam(1,2)=cand(ng); lam(2,1)=cand(ng); pv(ng)=ll_lambda(lam,n,s2,vl,real(nul,dp))
        end if
      end do
      if (ng>0) then; pv(1:ng)=exp(pv(1:ng)-maxval(pv(1:ng))); lam(1,2)=cand(rand_categorical(pv(1:ng))); lam(2,1)=lam(1,2); end if
      deallocate(cand,pv); allocate(cand(size(gl22)),pv(size(gl22))); ng=0
      do j=1,size(gl22)
        if (gl22(j)>lam(1,2)*lam(1,2)/lam(1,1)) then
          ng=ng+1; cand(ng)=gl22(j); lam(2,2)=cand(ng); pv(ng)=ll_lambda(lam,n,s2,vl,real(nul,dp))
        end if
      end do
      if (ng>0) then; pv(1:ng)=exp(pv(1:ng)-maxval(pv(1:ng))); lam(2,2)=cand(rand_categorical(pv(1:ng))); end if
      deallocate(cand,pv)
      ! Nearest-neighbour MH update for e, using the exact GHK rectangle likelihood.
      idx=minloc(abs(ge-ecur),dim=1); q=idx
      if (idx==1) q=2
      if (idx==size(ge)) q=size(ge)-1
      if (idx>1 .and. idx<size(ge)) q=idx+merge(1,-1,rand_uniform()<0.5_dp)
      oldll=scale_logprob(x,k,ge(idx),mu,tau,sigma,sigmat,ndghk)
      newll=scale_logprob(x,k,ge(q),mu,tau,sigma,sigmat,ndghk); ratio=min(1.0_dp,exp(min(0.0_dp,newll-oldll)))
      if (rand_uniform()<ratio) ecur=ge(q)
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%mudraw(mkeep,:)=mu; out%sigmadraw(mkeep,:,:)=sigmat; out%taudraw(mkeep,:)=tau
        out%sdraw(mkeep,:)=sigma; out%lambdadraw(mkeep,:,:)=lam; out%edraw(mkeep)=ecur
      end if
    end do
  end function rscale_usage

end module bayesm_scale
