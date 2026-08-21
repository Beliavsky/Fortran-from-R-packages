! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_imputation
  use compositions_kinds, only: dp, pi
  use compositions_geometry, only: closure
  use compositions_linalg, only: pseudoinverse, solve_least_squares, covariance_matrix, chol_lower
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none
  private
  integer, parameter, public :: mt_observed=0, mt_bdl=1, mt_mar=2, mt_sz=3, mt_error=4, mt_mnar=5
  public :: classify_missingness, missing_pattern_indices, conditional_alr_moments
  public :: impute_acomp_conditional, fit_acomp_projection, fit_acomp_em

  type, public :: acomp_imputation_result
    real(dp), allocatable :: composition(:,:)
    real(dp), allocatable :: clr(:,:)
    real(dp), allocatable :: beta(:,:)
    real(dp), allocatable :: clr_cov(:,:)
    integer, allocatable :: accepted(:)
    logical :: ok=.false.
  end type
contains
  subroutine classify_missingness(comp, mt, dl, detection_limit)
    real(dp), intent(in) :: comp(:,:)
    integer, allocatable, intent(out) :: mt(:,:)
    real(dp), allocatable, intent(out) :: dl(:,:)
    real(dp), intent(in), optional :: detection_limit
    integer :: i,j
    real(dp) :: dgen,v
    dgen=0.0_dp; if(present(detection_limit)) dgen=detection_limit
    allocate(mt(size(comp,1),size(comp,2)),dl(size(comp,1),size(comp,2))); dl=0.0_dp
    do j=1,size(comp,2); do i=1,size(comp,1)
      v=comp(i,j)
      if(ieee_is_finite(v)) then
        if(v>0.0_dp) then
          mt(i,j)=mt_observed
        else
          mt(i,j)=mt_bdl
          if(v<0.0_dp) then; dl(i,j)=-v; else; dl(i,j)=dgen; end if
        end if
      else if(ieee_is_nan(v)) then
        mt(i,j)=mt_mar
      else if(v>0.0_dp) then
        mt(i,j)=mt_error
      else
        mt(i,j)=mt_sz
      end if
    end do; end do
  end subroutine classify_missingness

  subroutine missing_pattern_indices(mt, pattern_id, representative, index_order, nmissing)
    integer, intent(in) :: mt(:,:)
    integer, allocatable, intent(out) :: pattern_id(:),representative(:),index_order(:,:),nmissing(:)
    integer :: n,d,i,j,k,np
    logical :: same
    n=size(mt,1); d=size(mt,2)
    allocate(pattern_id(n),representative(n),index_order(n,d),nmissing(n))
    pattern_id=0; representative=0; index_order=0; nmissing=0; np=0
    do i=1,n
      k=np+1
      do k=1,np
        same=.true.
        do j=1,d
          if((mt(i,j)==mt_observed).neqv.(mt(representative(k),j)==mt_observed)) then
            same=.false.; exit
          end if
        end do
        if(same) exit
      end do
      if(k>np.or.np==0) then
        np=np+1; representative(np)=i; k=np
      end if
      pattern_id(i)=k
    end do
    do k=1,np
      i=representative(k); nmissing(k)=count(mt(i,:)/=mt_observed); j=0
      call fill_order(mt(i,:),index_order(k,:),j)
    end do
    representative=representative(1:np); index_order=index_order(1:np,:); nmissing=nmissing(1:np)
  contains
    subroutine fill_order(row,ord,pos)
      integer, intent(in) :: row(:)
      integer, intent(out) :: ord(:)
      integer, intent(inout) :: pos
      integer :: q
      pos=0
      do q=1,size(row); if(row(q)/=mt_observed) then; pos=pos+1; ord(pos)=q; end if; end do
      do q=1,size(row); if(row(q)==mt_observed) then; pos=pos+1; ord(pos)=q; end if; end do
    end subroutine fill_order
  end subroutine missing_pattern_indices

  subroutine conditional_alr_moments(clr_cov, order, nmiss, lambda, residual_cov)
    real(dp), intent(in) :: clr_cov(:,:)
    integer, intent(in) :: order(:),nmiss
    real(dp), allocatable, intent(out) :: lambda(:,:),residual_cov(:,:)
    real(dp), allocatable :: a(:,:),vgg(:,:),vgm(:,:),vmg(:,:),vmm(:,:),vggi(:,:)
    integer :: d,i,j,ref,ng
    d=size(clr_cov,1); ng=d-nmiss; ref=order(d)
    if(size(clr_cov,2)/=d.or.size(order)/=d) error stop 'conditional_alr_moments: dimension mismatch'
    allocate(a(d,d)); a=0.0_dp
    do i=1,d; do j=1,d
      a(i,j)=clr_cov(order(i),order(j))+clr_cov(ref,ref) &
        -clr_cov(ref,order(j))-clr_cov(order(i),ref)
    end do; end do
    if(nmiss==0) then
      allocate(lambda(0,ng),residual_cov(0,0)); return
    end if
    allocate(vmm(nmiss,nmiss)); vmm=a(1:nmiss,1:nmiss)
    if(ng<2) then
      allocate(lambda(nmiss,max(0,ng))); lambda=0.0_dp; residual_cov=vmm; return
    end if
    allocate(vgg(ng,ng),vgm(ng,nmiss),vmg(nmiss,ng))
    vgg=a(nmiss+1:d,nmiss+1:d); vgm=a(nmiss+1:d,1:nmiss); vmg=transpose(vgm)
    call pseudoinverse(vgg,vggi,1.0e-8_dp)
    lambda=matmul(vmg,vggi)
    residual_cov=vmm-matmul(lambda,vgm)
    residual_cov=0.5_dp*(residual_cov+transpose(residual_cov))
  end subroutine conditional_alr_moments

  function impute_acomp_conditional(comp,pred_clr,clr_cov,mt,dl,nsim,seed) result(res)
    real(dp), intent(in) :: comp(:,:),pred_clr(:,:),clr_cov(:,:)
    integer, intent(in) :: mt(:,:)
    real(dp), intent(in), optional :: dl(:,:)
    integer, intent(in), optional :: nsim,seed
    type(acomp_imputation_result) :: res
    integer :: n,d,i,j,k,nm,ng,ref,ns,nacc,t,info,nseed
    integer, allocatable :: ord(:),put(:)
    real(dp), allocatable :: lambda(:,:),rcov(:,:),lchol(:,:),obs(:),pm(:),pg(:),mis(:),draw(:),alr(:)
    real(dp) :: lref,pref,lim
    logical :: accept
    n=size(comp,1); d=size(comp,2)
    if(any(shape(pred_clr)/=[n,d]).or.any(shape(mt)/=[n,d])) error stop 'impute_acomp_conditional: shape mismatch'
    if(any(shape(clr_cov)/=[d,d])) error stop 'impute_acomp_conditional: covariance mismatch'
    ns=0; if(present(nsim)) ns=max(0,nsim)
    if(present(seed)) then
      call random_seed(size=nseed); allocate(put(nseed));
      do k=1,nseed; put(k)=mod(abs(seed)+104729*k,huge(1)-1)+1; end do
      call random_seed(put=put)
    end if
    allocate(res%composition(n,d),res%clr(n,d),res%accepted(n),ord(d),alr(d)); res%accepted=0
    do i=1,n
      nm=count(mt(i,:)/=mt_observed); ng=d-nm; k=0
      do j=1,d; if(mt(i,j)/=mt_observed) then; k=k+1; ord(k)=j; end if; end do
      do j=1,d; if(mt(i,j)==mt_observed) then; k=k+1; ord(k)=j; end if; end do
      if(ng>0) then; ref=ord(d); else; ref=d; end if
      if(nm==0) then
        res%composition(i,:)=closure(comp(i,:)); res%clr(i,:)=log(res%composition(i,:))
        res%clr(i,:)=res%clr(i,:)-sum(res%clr(i,:))/real(d,dp); cycle
      end if
      if(ng==0) then
        res%clr(i,:)=pred_clr(i,:); res%composition(i,:)=closure(exp(pred_clr(i,:)-maxval(pred_clr(i,:))))
        cycle
      end if
      call conditional_alr_moments(clr_cov,ord,nm,lambda,rcov)
      allocate(obs(ng),pg(ng),pm(nm),mis(nm),draw(nm))
      lref=log(comp(i,ref)); pref=pred_clr(i,ref)
      do j=1,ng
        obs(j)=log(comp(i,ord(nm+j)))-lref
        pg(j)=pred_clr(i,ord(nm+j))-pref
      end do
      do j=1,nm; pm(j)=pred_clr(i,ord(j))-pref; end do
      mis=pm
      if(ng>=2) mis=pm+matmul(lambda,obs-pg)
      nacc=0
      if(ns>0.and.nm>0) then
        call safe_chol(rcov,lchol,info)
        if(info==0) then
          draw=0.0_dp
          do t=1,ns
            call normal_vector(nm,alr(1:nm))
            alr(1:nm)=mis+matmul(lchol,alr(1:nm)); accept=.true.
            if(present(dl)) then
              do j=1,nm
                if(mt(i,ord(j))==mt_bdl.and.dl(i,ord(j))>0.0_dp) then
                  lim=log(dl(i,ord(j)))-lref
                  if(alr(j)>=lim) then; accept=.false.; exit; end if
                end if
              end do
            end if
            if(accept) then; draw=draw+alr(1:nm); nacc=nacc+1; end if
          end do
          if(nacc>0) mis=draw/real(nacc,dp)
        end if
      end if
      if(nacc==0.and.present(dl)) then
        do j=1,nm
          if(mt(i,ord(j))==mt_bdl.and.dl(i,ord(j))>0.0_dp) then
            lim=log(dl(i,ord(j)))-lref; mis(j)=min(mis(j),lim)
          end if
        end do
      end if
      res%accepted(i)=nacc; alr=0.0_dp
      do j=1,nm; alr(ord(j))=mis(j); end do
      do j=1,ng; alr(ord(nm+j))=obs(j); end do
      alr(ref)=0.0_dp
      res%clr(i,:)=alr-sum(alr)/real(d,dp)
      res%composition(i,:)=closure(exp(res%clr(i,:)-maxval(res%clr(i,:))))
      deallocate(obs,pg,pm,mis,draw,lambda,rcov)
      if(allocated(lchol)) deallocate(lchol)
    end do
    call covariance_matrix(res%clr,res%clr_cov)
    res%ok=.true.
  contains
    subroutine safe_chol(a,l,istat)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      integer, intent(out) :: istat
      real(dp), allocatable :: aa(:,:)
      real(dp) :: ridge
      integer :: tries,q
      ridge=0.0_dp
      do tries=1,7
        aa=a
        if(ridge>0.0_dp) then; do q=1,size(aa,1); aa(q,q)=aa(q,q)+ridge; end do; end if
        call chol_try(aa,l,istat)
        if(istat==0) return
        ridge=merge(1.0e-12_dp,10.0_dp*ridge,ridge==0.0_dp)
      end do
    end subroutine safe_chol
    subroutine chol_try(a,l,istat)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      integer, intent(out) :: istat
      integer :: q,r
      real(dp) :: v
      allocate(l(size(a,1),size(a,2))); l=0.0_dp; istat=0
      do q=1,size(a,1)
        v=a(q,q)-sum(l(q,1:q-1)**2)
        if(v<=0.0_dp) then; istat=1; return; end if
        l(q,q)=sqrt(v)
        do r=q+1,size(a,1)
          l(r,q)=(a(r,q)-sum(l(r,1:q-1)*l(q,1:q-1)))/l(q,q)
        end do
      end do
    end subroutine chol_try
    subroutine normal_vector(m,z)
      integer, intent(in) :: m
      real(dp), intent(out) :: z(m)
      integer :: q
      real(dp) :: a,b
      q=1
      do while(q<=m)
        call random_number(a); call random_number(b); a=max(a,tiny(1.0_dp))
        z(q)=sqrt(-2.0_dp*log(a))*cos(2.0_dp*pi*b)
        if(q+1<=m) z(q+1)=sqrt(-2.0_dp*log(a))*sin(2.0_dp*pi*b)
        q=q+2
      end do
    end subroutine normal_vector
  end function impute_acomp_conditional

  function fit_acomp_projection(comp,x,mt,source_compatible) result(res)
    real(dp), intent(in) :: comp(:,:),x(:,:)
    integer, intent(in) :: mt(:,:)
    logical, intent(in), optional :: source_compatible
    type(acomp_imputation_result) :: res
    real(dp), allocatable :: xb(:,:),yb(:,:),bvec(:,:),rr(:,:),proj(:,:),target(:)
    integer :: n,d,p,c,i,j,l,row,col,no
    real(dp) :: obnm,mu,xv
    logical :: src
    n=size(comp,1); d=size(comp,2); p=size(x,2)
    if(size(x,1)/=n.or.any(shape(mt)/=[n,d])) error stop 'fit_acomp_projection: shape mismatch'
    src=.false.; if(present(source_compatible)) src=source_compatible
    allocate(xb(n*d,p*d),yb(n*d,1),proj(d,d),target(d)); xb=0.0_dp; yb=0.0_dp
    do c=1,n
      no=count(mt(c,:)==mt_observed); proj=0.0_dp; target=0.0_dp
      if(no>0) then
        obnm=-1.0_dp/real(no,dp)
        do i=1,d; do j=1,d
          if(mt(c,i)==mt_observed.and.mt(c,j)==mt_observed) proj(i,j)=obnm+merge(1.0_dp,0.0_dp,i==j)
        end do; end do
        mu=0.0_dp
        do i=1,d; if(mt(c,i)==mt_observed) mu=mu+log(comp(c,i)); end do
        if(.not.src) mu=mu/real(no,dp)
        do i=1,d; if(mt(c,i)==mt_observed) target(i)=log(comp(c,i))-mu; end do
      end if
      do i=1,d
        row=(c-1)*d+i; yb(row,1)=target(i)
        do j=1,d; do l=1,p
          col=(j-1)*p+l; xv=merge(1.0_dp,x(c,l),src)
          xb(row,col)=xv*proj(i,j)
        end do; end do
      end do
    end do
    call solve_least_squares(xb,yb,bvec,rr)
    allocate(res%beta(p,d));
    do j=1,d; res%beta(:,j)=bvec((j-1)*p+1:j*p,1); end do
    res%clr=matmul(x,res%beta); allocate(res%composition(n,d))
    do c=1,n; res%composition(c,:)=closure(exp(res%clr(c,:)-maxval(res%clr(c,:)))); end do
    call covariance_matrix(res%clr,res%clr_cov); allocate(res%accepted(n)); res%accepted=0; res%ok=.true.
  end function fit_acomp_projection

  function fit_acomp_em(comp,x,mt,dl,steps,seed) result(res)
    real(dp), intent(in) :: comp(:,:),x(:,:)
    integer, intent(in) :: mt(:,:)
    real(dp), intent(in), optional :: dl(:,:)
    integer, intent(in), optional :: steps,seed
    type(acomp_imputation_result) :: res,tmp,imp
    real(dp), allocatable :: resid(:,:),cov(:,:),b(:,:),rr(:,:),y(:,:)
    integer :: it,nit,c,j
    nit=10; if(present(steps)) nit=max(1,steps)
    tmp=fit_acomp_projection(comp,x,mt,source_compatible=.false.)
    do it=1,nit
      resid=0.0_dp*tmp%clr
      do c=1,size(comp,1)
        do j=1,size(comp,2)
          if(mt(c,j)==mt_observed) resid(c,j)=log(comp(c,j))
        end do
        if(count(mt(c,:)==mt_observed)>0) then
          resid(c,:)=resid(c,:)-sum(resid(c,:),mask=mt(c,:)==mt_observed) &
            /real(count(mt(c,:)==mt_observed),dp)
        end if
      end do
      call covariance_matrix(resid,cov)
      if(present(dl)) then
        imp=impute_acomp_conditional(comp,tmp%clr,cov,mt,dl,nsim=200,seed=seed)
      else
        imp=impute_acomp_conditional(comp,tmp%clr,cov,mt,nsim=0,seed=seed)
      end if
      y=imp%clr; call solve_least_squares(x,y,b,rr); tmp%beta=b; tmp%clr=matmul(x,b)
      do c=1,size(comp,1); tmp%composition(c,:)=closure(exp(tmp%clr(c,:)-maxval(tmp%clr(c,:)))); end do
    end do
    res=imp; res%beta=tmp%beta; res%clr=tmp%clr; res%composition=tmp%composition
    call covariance_matrix(y-matmul(x,res%beta),res%clr_cov); res%ok=.true.
  end function fit_acomp_em
end module compositions_imputation
