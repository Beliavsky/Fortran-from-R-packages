! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_stats
  use compositions_kinds, only: dp
  use compositions_geometry, only: closure, clr_rows, clr_inv, cpt_rows, ilr_rows, ilr_inv_rows, ilr_base, &
    ilrvar_to_clr, clrvar_to_variation, acomp_mean, build_ilr_base, balance_coordinate
  use compositions_linalg, only: covariance_matrix, symmetric_eigen, invert_matrix, determinant_spd, solve_least_squares
  use robustbase_detmcd, only: detmcd_result, cov_detmcd
  use robustbase_probability, only: chi_square_cdf
  implicit none
  private
  public :: fit_dirichlet_result, fit_dirichlet, compositional_covariance, robust_compositional_covariance
  public :: pca_result, compositional_pca, normal_location_result, acomp_normal_location_one_sample
  public :: acomp_normal_location_two_sample_equal, compositional_lm_result, compositional_lm_fit
  public :: compositional_lm_predict, principal_balance_maxvar, principal_balance_hclust, &
    principal_balance_angprox, digamma_approx, trigamma_approx
  public :: mahalanobis_distances

  type :: fit_dirichlet_result
    real(dp), allocatable :: alpha(:)
    real(dp) :: loglikelihood = 0.0_dp
    integer :: df = 0
    integer :: iterations = 0
    logical :: converged = .false.
  end type

  type :: pca_result
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: ilr_loadings(:,:)
    real(dp), allocatable :: clr_loadings(:,:)
    real(dp), allocatable :: scores(:,:)
  end type

  type :: normal_location_result
    real(dp) :: statistic=0.0_dp
    integer :: df=0
    real(dp) :: p_value=1.0_dp
  end type

  type :: compositional_lm_result
    real(dp), allocatable :: coefficients(:,:)
    real(dp), allocatable :: residuals(:,:)
    real(dp), allocatable :: fitted_ilr(:,:)
    real(dp), allocatable :: residual_covariance(:,:)
    real(dp), allocatable :: basis(:,:)
    integer :: df_residual=0
    logical :: ok=.false.
  end type
contains
  real(dp) function digamma_approx(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z,r
    if(x<=0.0_dp) error stop 'digamma_approx: x must be positive'
    z=x; y=0.0_dp
    do while(z<8.0_dp); y=y-1.0_dp/z; z=z+1.0_dp; end do
    r=1.0_dp/z
    y=y+log(z)-0.5_dp*r-r*r*(1.0_dp/12.0_dp-r*r*(1.0_dp/120.0_dp-r*r/252.0_dp))
  end function digamma_approx

  real(dp) function trigamma_approx(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z,r,r2
    if(x<=0.0_dp) error stop 'trigamma_approx: x must be positive'
    z=x; y=0.0_dp
    do while(z<8.0_dp); y=y+1.0_dp/(z*z); z=z+1.0_dp; end do
    r=1.0_dp/z; r2=r*r
    y=y+r+0.5_dp*r2+r*r2/6.0_dp-r*r2*r2/30.0_dp+r*r2*r2*r2/42.0_dp
  end function trigamma_approx

  function fit_dirichlet(x,alpha0,max_iter,tol) result(res)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: alpha0(:),tol
    integer, intent(in), optional :: max_iter
    type(fit_dirichlet_result) :: res
    real(dp), allocatable :: alpha(:),elog(:),e(:),v(:,:),vinv(:,:),update(:),trial(:)
    real(dp) :: delta,eps,step,ll
    integer :: d,n,it,mi,i,j
    n=size(x,1); d=size(x,2)
    if(any(x<=0.0_dp)) error stop 'fit_dirichlet: positive compositions required'
    elog=sum(log(x),dim=1)/real(n,dp)
    allocate(alpha(d)); alpha=1.0_dp
    if(present(alpha0)) then; if(size(alpha0)/=d) error stop 'fit_dirichlet: alpha0 mismatch'; alpha=alpha0; end if
    mi=100; if(present(max_iter)) mi=max_iter; eps=1.0e-9_dp; if(present(tol)) eps=tol
    allocate(e(d),v(d,d),update(d),trial(d))
    do it=1,mi
      do i=1,d; e(i)=digamma_approx(alpha(i))-digamma_approx(sum(alpha)); end do
      v=-trigamma_approx(sum(alpha))
      do i=1,d; v(i,i)=v(i,i)+trigamma_approx(alpha(i)); end do
      call invert_matrix(v,vinv); update=matmul(vinv,elog-e); delta=sqrt(sum((elog-e)**2))
      step=1.0_dp
      do
        trial=alpha+step*update
        if(all(trial>0.0_dp)) exit
        step=0.5_dp*step
        if(step<1.0e-12_dp) exit
      end do
      alpha=trial
      if(delta<eps) exit
    end do
    ll=-real(n,dp)*(sum(log_gamma(alpha))-log_gamma(sum(alpha))+sum(elog*(alpha-1.0_dp)))
    res%alpha=alpha; res%loglikelihood=ll; res%df=n*(d-1)-d; res%iterations=it; res%converged=(delta<eps)
  end function fit_dirichlet

  subroutine compositional_covariance(x,cov,center)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:)
    real(dp), allocatable, intent(out), optional :: center(:)
    real(dp), allocatable :: z(:,:),mu(:)
    z=clr_rows(x); call covariance_matrix(z,cov,mu)
    if(present(center)) center=acomp_mean(x)
  end subroutine compositional_covariance

  subroutine robust_compositional_covariance(x,cov,center,alpha)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:),center(:)
    real(dp), intent(in), optional :: alpha
    real(dp), allocatable :: z(:,:),v(:,:),ctmp(:,:)
    type(detmcd_result) :: r
    real(dp) :: a
    z=ilr_rows(x); a=0.5_dp; if(present(alpha)) a=alpha
    call cov_detmcd(z,r,alpha=a)
    v=ilr_base(size(x,2)); cov=ilrvar_to_clr(r%estimate%covariance,v)
    ctmp=ilr_inv_rows(reshape(r%estimate%center,[1,size(r%estimate%center)]),v)
    center=ctmp(1,:)
  end subroutine robust_compositional_covariance

  function compositional_pca(x,robust) result(res)
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: robust
    type(pca_result) :: res
    real(dp), allocatable :: z(:,:),cov(:,:),mu(:),vals(:),vecs(:,:),v(:,:),tmp(:,:),ctmp(:,:)
    type(detmcd_result) :: rr
    logical :: rb
    integer :: info,n,p,i
    z=ilr_rows(x); n=size(z,1); p=size(z,2); rb=.false.; if(present(robust)) rb=robust
    if(rb) then
      call cov_detmcd(z,rr); cov=rr%estimate%covariance; mu=rr%estimate%center
    else
      call covariance_matrix(z,cov,mu)
    end if
    call symmetric_eigen(cov,vals,vecs,info)
    ! dsyev returns ascending order; reverse for PCA.
    allocate(res%eigenvalues(p),res%ilr_loadings(p,p),res%scores(n,p))
    do i=1,p
      res%eigenvalues(i)=vals(p-i+1); res%ilr_loadings(:,i)=vecs(:,p-i+1)
    end do
    v=ilr_base(size(x,2)); res%clr_loadings=matmul(v,res%ilr_loadings)
    ctmp=ilr_inv_rows(reshape(mu,[1,p]),v)
    res%center=ctmp(1,:)
    res%scores=matmul(z-spread(mu,1,n),res%ilr_loadings)
  end function compositional_pca

  function acomp_normal_location_one_sample(x) result(res)
    real(dp), intent(in) :: x(:,:)
    type(normal_location_result) :: res
    real(dp), allocatable :: w(:,:),v(:,:),tss(:,:),rss(:,:),mu(:)
    real(dp) :: dt,dr
    integer :: n,m
    w=ilr_rows(x); n=size(w,1); m=size(w,2); mu=sum(w,dim=1)/real(n,dp)
    v=w-spread(mu,1,n); tss=matmul(transpose(w),w); rss=matmul(transpose(v),v)
    dt=determinant_spd(tss/real(n,dp)); dr=determinant_spd(rss/real(n,dp))
    if(dt<=0.0_dp.or.dr<=0.0_dp) then; res%statistic=0.0_dp; res%p_value=1.0_dp
    else
      res%statistic=real(n,dp)*(log(dt)-log(dr)); res%p_value=1.0_dp-chi_square_cdf(res%statistic,real(m,dp))
    end if
    res%df=m
  end function acomp_normal_location_one_sample

  function acomp_normal_location_two_sample_equal(x,y) result(res)
    real(dp), intent(in) :: x(:,:),y(:,:)
    type(normal_location_result) :: res
    real(dp), allocatable :: zx(:,:),zy(:,:),allz(:,:),tss(:,:),rss(:,:),mx(:),my(:),ma(:)
    real(dp) :: dt,dr
    integer :: nx,ny,n,m
    zx=ilr_rows(x); zy=ilr_rows(y); nx=size(zx,1); ny=size(zy,1); n=nx+ny; m=size(zx,2)
    allocate(allz(n,m)); allz(1:nx,:)=zx; allz(nx+1:n,:)=zy
    mx=sum(zx,dim=1)/real(nx,dp); my=sum(zy,dim=1)/real(ny,dp); ma=sum(allz,dim=1)/real(n,dp)
    tss=matmul(transpose(allz-spread(ma,1,n)),allz-spread(ma,1,n))
    rss=matmul(transpose(zx-spread(mx,1,nx)),zx-spread(mx,1,nx))+ &
        matmul(transpose(zy-spread(my,1,ny)),zy-spread(my,1,ny))
    dt=determinant_spd(tss/real(n,dp)); dr=determinant_spd(rss/real(n,dp))
    res%statistic=real(n,dp)*(log(max(dt,tiny(1.0_dp)))-log(max(dr,tiny(1.0_dp))))
    res%df=m; res%p_value=1.0_dp-chi_square_cdf(max(res%statistic,0.0_dp),real(m,dp))
  end function acomp_normal_location_two_sample_equal

  function compositional_lm_fit(x,y,basis) result(res)
    real(dp), intent(in) :: x(:,:),y(:,:)
    real(dp), intent(in), optional :: basis(:,:)
    type(compositional_lm_result) :: res
    real(dp), allocatable :: z(:,:),v(:,:),beta(:,:),r(:,:),cov(:,:)
    integer :: info,n,p
    if(size(x,1)/=size(y,1)) error stop 'compositional_lm_fit: row mismatch'
    if(present(basis)) then; v=basis; z=ilr_rows(y,v); else; v=ilr_base(size(y,2)); z=ilr_rows(y,v); end if
    call solve_least_squares(x,z,beta,r,info); res%coefficients=beta; res%residuals=r
    res%fitted_ilr=matmul(x,beta); res%basis=v; n=size(x,1); p=size(x,2); res%df_residual=max(0,n-p)
    if(res%df_residual>0) then
      allocate(cov(size(z,2),size(z,2))); cov=matmul(transpose(r),r)/real(res%df_residual,dp); res%residual_covariance=cov
    else
      allocate(res%residual_covariance(size(z,2),size(z,2))); res%residual_covariance=0.0_dp
    end if
    res%ok=(info==0)
  end function compositional_lm_fit

  function compositional_lm_predict(model,xnew) result(yhat)
    type(compositional_lm_result), intent(in) :: model
    real(dp), intent(in) :: xnew(:,:)
    real(dp), allocatable :: yhat(:,:),z(:,:)
    z=matmul(xnew,model%coefficients); yhat=ilr_inv_rows(z,model%basis)
  end function compositional_lm_predict

  function mahalanobis_distances(x,center,cov) result(d)
    real(dp), intent(in) :: x(:,:),center(:),cov(:,:)
    real(dp) :: d(size(x,1))
    real(dp), allocatable :: inv(:,:),u(:)
    integer :: i
    call invert_matrix(cov,inv)
    do i=1,size(x,1); u=x(i,:)-center; d(i)=sqrt(max(0.0_dp,dot_product(u,matmul(inv,u)))); end do
  end function mahalanobis_distances

  function principal_balance_maxvar(x) result(v)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: v(:,:)
    real(dp), allocatable :: xwork(:,:),xorig(:,:),cov(:,:),vals(:),vecs(:,:),vb(:,:),proj(:,:)
    integer, allocatable :: code(:,:),status(:,:),active(:),suggest(:),trial(:),newstatus(:)
    real(dp), allocatable :: newvar(:)
    integer :: d,iorder,i,j,r,s,info,bestj,iter,lev
    real(dp) :: prevvar,bestvar

    d=size(x,2)
    if(d<2) then
      allocate(v(d,0))
      return
    end if

    ! gsi.PrinBal(method="PBmaxvar") passes cpt(x) into gsi.mvOPTIMAL.
    xorig=cpt_rows(x); xwork=xorig
    allocate(code(d-1,d),status(d,d),active(d),suggest(d),trial(d),newvar(d),newstatus(d))
    code=0; status=0; active=1

    do iorder=1,d-1
      ! gsi.mvSVDSUGG: leading right singular vector, then mask inactive parts.
      call covariance_matrix(xwork,cov)
      call symmetric_eigen(cov,vals,vecs,info)
      suggest=0
      do j=1,d
        if(active(j)/=0) then
          if(vecs(j,size(vecs,2))>=0.0_dp) then
            suggest(j)=1
          else
            suggest(j)=-1
          end if
        end if
      end do
      if(count(suggest>0)==0.or.count(suggest<0)==0) then
        call fallback_split(active,suggest)
      end if
      code(iorder,:)=suggest

      ! gsi.mvREFINESBP: greedily flip one sign when variance increases.
      trial=code(iorder,:); r=count(trial==1); s=count(trial==-1)
      prevvar=balance_variance_source(xorig,trial); iter=0
      do while(iter<50)
        newvar=0.0_dp
        do j=1,d
          if(trial(j)==1.and.r>1) then
            suggest=trial; suggest(j)=-1; newvar(j)=balance_variance_source(xorig,suggest)
          else if(trial(j)==-1.and.s>1) then
            suggest=trial; suggest(j)=1; newvar(j)=balance_variance_source(xorig,suggest)
          end if
        end do
        bestj=maxloc(newvar,dim=1); bestvar=newvar(bestj)
        if(bestvar<=prevvar) exit
        if(trial(bestj)==1) then
          trial(bestj)=-1; r=r-1; s=s+1
        else if(trial(bestj)==-1) then
          trial(bestj)=1; r=r+1; s=s-1
        else
          exit
        end if
        prevvar=bestvar; iter=iter+1
      end do
      code(iorder,:)=trial

      if(iorder<d-1) then
        ! gsi.mvPROJclr: remove the subspace of completed balances.
        vb=build_ilr_base(transpose(code(1:iorder,:)))
        proj=matmul(vb,transpose(vb))
        xwork=xwork-matmul(xwork,proj)

        ! gsi.mvNEXTACTIVE state update.
        r=count(code(iorder,:)==1); s=count(code(iorder,:)==-1)
        do j=1,d
          if(code(iorder,j)==0) then
            if(status(iorder,j)/=d+1) then
              newstatus(j)=status(iorder,j)+1
            else
              newstatus(j)=d+1
            end if
          else if(code(iorder,j)==-1.and.status(iorder,j)/=d+1) then
            if(s==1) then
              newstatus(j)=d+1
            else
              newstatus(j)=status(iorder,j)+1
            end if
          else if(code(iorder,j)==1.and.status(iorder,j)/=d+1) then
            if(r==1) then
              newstatus(j)=d+1
            else
              newstatus(j)=status(iorder,j)
            end if
          else
            newstatus(j)=d+1
          end if
        end do
        status(iorder+1,:)=newstatus
        lev=minval(status(iorder+1,:)); active=0
        where(status(iorder+1,:)==lev) active=1
      end if
    end do
    v=build_ilr_base(transpose(code))
  contains
    subroutine fallback_split(act,sgn)
      integer, intent(in) :: act(:)
      integer, intent(inout) :: sgn(:)
      integer :: nact,k,cut
      nact=count(act/=0); cut=max(1,nact/2); k=0; sgn=0
      do j=1,size(act)
        if(act(j)/=0) then
          k=k+1
          if(k<=cut) then; sgn(j)=1; else; sgn(j)=-1; end if
        end if
      end do
    end subroutine fallback_split

    real(dp) function balance_variance_source(a,sgn) result(vr)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: sgn(:)
      real(dp) :: bal(size(a,1)),mu,cp,cm
      integer :: ii,jj,rn,sn
      rn=count(sgn==1); sn=count(sgn==-1); bal=0.0_dp
      if(rn==0.or.sn==0) then; vr=0.0_dp; return; end if
      cp=sqrt(real(sn,dp)/real(rn*(rn+sn),dp))
      cm=sqrt(real(rn,dp)/real(sn*(rn+sn),dp))
      do ii=1,size(a,1)
        do jj=1,size(a,2)
          if(sgn(jj)==1) bal(ii)=bal(ii)+cp*a(ii,jj)
          if(sgn(jj)==-1) bal(ii)=bal(ii)-cm*a(ii,jj)
        end do
      end do
      mu=sum(bal)/real(size(bal),dp)
      vr=sum((bal-mu)**2)/real(max(1,size(bal)-1),dp)
    end function balance_variance_source
  end function principal_balance_maxvar

  function principal_balance_hclust(x) result(v)
    !! gsi.PrinBal(method="PBhclust"): Ward clustering of the variation distances.
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: v(:,:)
    real(dp), allocatable :: lx(:,:),dist2(:,:)
    integer, allocatable :: members(:,:),nmem(:),active(:),code(:,:)
    integer :: n,d,i,j,k,a,b,nc,step,best_a,best_b,newidx
    real(dp) :: vr,best,duv,ds,dt,den
    d=size(x,2); n=size(x,1)
    if(d<2) then; allocate(v(d,0)); return; end if
    if(any(x<=0.0_dp)) error stop 'principal_balance_hclust: positive compositions required'
    allocate(lx(n,d)); lx=log(x)
    do j=1,d; lx(:,j)=lx(:,j)-sum(lx(:,j))/real(n,dp); end do
    allocate(dist2(2*d-1,2*d-1)); dist2=huge(1.0_dp)
    do i=1,d; dist2(i,i)=0.0_dp; end do
    do i=1,d-1; do j=i+1,d
      vr=sum((lx(:,i)-lx(:,j))**2)/real(max(1,n-1),dp)
      dist2(i,j)=vr; dist2(j,i)=vr
    end do; end do
    allocate(members(2*d-1,d),nmem(2*d-1),active(2*d-1),code(d-1,d))
    members=0; nmem=0; active=0; code=0
    do i=1,d; members(i,1)=i; nmem(i)=1; active(i)=1; end do
    nc=d
    do step=1,d-1
      best=huge(1.0_dp); best_a=0; best_b=0
      do a=1,nc-1; if(active(a)==0) cycle
        do b=a+1,nc; if(active(b)==0) cycle
          if(dist2(a,b)<best) then; best=dist2(a,b); best_a=a; best_b=b; end if
        end do
      end do
      if(best_a==0) error stop 'principal_balance_hclust: clustering failed'
      do k=1,nmem(best_a); code(step,members(best_a,k))=-1; end do
      do k=1,nmem(best_b); code(step,members(best_b,k))=1; end do
      nc=nc+1; newidx=nc; nmem(newidx)=nmem(best_a)+nmem(best_b)
      members(newidx,1:nmem(best_a))=members(best_a,1:nmem(best_a))
      members(newidx,nmem(best_a)+1:nmem(newidx))=members(best_b,1:nmem(best_b))
      active(best_a)=0; active(best_b)=0; active(newidx)=1; dist2(newidx,newidx)=0.0_dp
      do k=1,newidx-1; if(active(k)==0) cycle
        den=real(nmem(best_a)+nmem(best_b)+nmem(k),dp)
        ds=dist2(best_a,k); dt=dist2(best_b,k); duv=dist2(best_a,best_b)
        dist2(newidx,k)=(real(nmem(best_a)+nmem(k),dp)*ds &
          +real(nmem(best_b)+nmem(k),dp)*dt-real(nmem(k),dp)*duv)/den
        dist2(k,newidx)=dist2(newidx,k)
      end do
    end do
    v=build_ilr_base(transpose(code))
  end function principal_balance_hclust

  function principal_balance_angprox(x) result(v)
    !! gsi.PrinBal(method="PBangprox"): exhaustive angular proximity recursion.
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: v(:,:)
    real(dp), allocatable :: cptx(:,:),cov(:,:),vals(:),vecs(:,:)
    integer, allocatable :: code(:,:),parts(:)
    integer :: d,info,pos
    d=size(x,2)
    if(d<2) then; allocate(v(d,0)); return; end if
    cptx=cpt_rows(x); call covariance_matrix(cptx,cov); call symmetric_eigen(cov,vals,vecs,info)
    allocate(code(d-1,d),parts(d)); code=0; parts=[(pos,pos=1,d)]; pos=0
    call recurse(parts,d)
    v=build_ilr_base(transpose(code))
  contains
    recursive subroutine recurse(sub,nsub)
      integer, intent(in) :: sub(:),nsub
      integer, allocatable :: left(:),right(:),sgn(:)
      integer :: nl,nr,q
      if(nsub<=1) return
      call best_partition(sub,nsub,sgn)
      pos=pos+1; code(pos,:)=0
      do q=1,nsub; code(pos,sub(q))=sgn(q); end do
      nl=count(sgn>0); nr=count(sgn<0); allocate(left(nl),right(nr)); nl=0; nr=0
      do q=1,nsub
        if(sgn(q)>0) then; nl=nl+1; left(nl)=sub(q); else; nr=nr+1; right(nr)=sub(q); end if
      end do
      if(nl>1) call recurse(left,nl)
      if(nr>1) call recurse(right,nr)
    end subroutine recurse

    subroutine best_partition(sub,nsub,sgn)
      integer, intent(in) :: sub(:),nsub
      integer, allocatable, intent(out) :: sgn(:)
      integer :: mask,q,r,np,nn,maxmask,pc
      real(dp) :: normv,dotv,best,cp,cm
      real(dp), allocatable :: cand(:)
      if(nsub>20) error stop 'principal_balance_angprox: exhaustive method limited to 20 parts'
      maxmask=2**nsub-2; allocate(sgn(nsub),cand(nsub)); best=-huge(1.0_dp); sgn=0
      do mask=1,maxmask
        np=popcnt(mask); nn=nsub-np
        if(np==0.or.nn==0) cycle
        cp=real(nn,dp); cm=-real(np,dp); cand=cm
        do q=1,nsub; if(btest(mask,q-1)) cand(q)=cp; end do
        normv=sqrt(sum(cand*cand)); cand=cand/normv
        do pc=1,d-1
          dotv=0.0_dp
          do q=1,nsub; dotv=dotv+cand(q)*vecs(sub(q),d-pc+1); end do
          if(dotv>best) then
            best=dotv
            do r=1,nsub; sgn(r)=merge(1,-1,cand(r)>0.0_dp); end do
          end if
        end do
      end do
    end subroutine best_partition
  end function principal_balance_angprox

end module compositions_stats
