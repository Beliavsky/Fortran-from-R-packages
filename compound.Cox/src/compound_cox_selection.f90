! SPDX-License-Identifier: GPL-2.0-only
module compound_cox_selection
  use compound_cox_kinds, only : dp
  use compound_cox_types, only : univariate_result, selection_result
  use compound_cox_math, only : chisq_cdf, shuffle_int
  use compound_cox_regression, only : cox_loglik_breslow
  use survival_cox, only : coxph_fit
  use survival_types, only : coxph_result, concordance_result
  use survival_stats, only : concordance_right
  implicit none
  private
  public :: uni_score, uni_wald, uni_selection
contains
  subroutine uni_score(time,status,x,res,d0)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    type(univariate_result),intent(out)::res
    real(dp),intent(in),optional::d0
    integer::n,p,i,j,k
    real(dp)::dd,s0,s1,s2,s,v
    n=size(time)
    p=size(x,2)
    dd=0
    if(present(d0))dd=d0
    allocate(res%beta(p),res%z(p),res%p(p),res%se(p))
    res%se=0
    do j=1,p
    s=0
    v=0
    do i=1,n
    if(status(i)==0)cycle
    s0=0
    s1=0
    s2=0
    do k=1,n
    if(time(k)>=time(i))then
    s0=s0+1
    s1=s1+x(k,j)
    s2=s2+x(k,j)**2
    end if
    end do
    s=s+x(i,j)-s1/s0
    v=v+s2/s0-(s1/s0)**2
    end do
    res%z(j)=s/(sqrt(max(v,0.0_dp))+dd)
    res%p(j)=1-chisq_cdf(res%z(j)**2,1.0_dp)
    res%beta(j)=s/v
    if(v>0)res%se(j)=1/sqrt(v)
    end do
  end subroutine uni_score

  subroutine uni_wald(time,status,x,res)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    type(univariate_result),intent(out)::res
    integer::p,j
    real(dp)::xx(size(time),1)
    type(coxph_result)::cr
    p=size(x,2)
    allocate(res%beta(p),res%se(p),res%z(p),res%p(p))
    do j=1,p
    xx(:,1)=x(:,j)
    call coxph_fit(time,status,xx,cr,method='efron')
    res%beta(j)=cr%coef(1)
    res%se(j)=sqrt(max(cr%var(1,1),0.0_dp))
    res%z(j)=res%beta(j)/res%se(j)
    res%p(j)=1-chisq_cdf(res%z(j)**2,1.0_dp)
    end do
  end subroutine uni_wald

  subroutine fit_scalar_cox(time,status,risk,beta,ll)
    real(dp),intent(in)::time(:),risk(:)
    integer,intent(in)::status(:)
    real(dp),intent(out)::beta,ll
    real(dp)::xx(size(time),1)
    type(coxph_result)::cr
    xx(:,1)=risk
    call coxph_fit(time,status,xx,cr,method='efron')
    beta=cr%coef(1)
    ll=cr%loglik
  end subroutine fit_scalar_cox

  subroutine uni_selection(time,status,x,res,p_value,kfold,use_score,d0,permutation,nperm)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    type(selection_result),intent(out)::res
    real(dp),intent(in),optional::p_value,d0
    integer,intent(in),optional::kfold,nperm
    logical,intent(in),optional::use_score,permutation
    type(univariate_result)::ur,urk
    real(dp)::pth,dd,beta,lltrain,llfull
    integer::kfold0,mperm,n,p,q,j,k,lo,hi,ns,i,kk
    logical::score0,perm0
    integer,allocatable::sel(:),selk(:),permidx(:)
    logical,allocatable::keep(:)
    real(dp),allocatable::cc(:),cccv(:),cctest(:),risktrain(:),w(:),tt(:),xx(:,:),xcut(:,:),ccall(:)
    integer,allocatable::ddv(:)
    type(concordance_result)::conc
    real(dp)::fcount
    pth=0.001_dp
    if(present(p_value))pth=p_value
    kfold0=10
    if(present(kfold))kfold0=kfold
    score0=.true.
    if(present(use_score))score0=use_score
    dd=0
    if(present(d0))dd=d0
    perm0=.false.
    if(present(permutation))perm0=permutation
    mperm=200
    if(present(nperm))mperm=nperm
    n=size(time)
    p=size(x,2)
    if(score0)then
    call uni_score(time,status,x,ur,dd)
    else
    call uni_wald(time,status,x,ur)
    end if
    q=count(ur%p<pth)
    if(q==0)error stop 'uni_selection: no feature selected'
    allocate(sel(q))
    kk=0
    do j=1,p
    if(ur%p(j)<pth)then
    kk=kk+1
    sel(kk)=j
    end if
    end do
    allocate(cc(n),cccv(n),cctest(n),xcut(n,q),w(q))
    xcut=x(:,sel)
    if(score0)then
    w=ur%z(sel)
    else
    w=ur%beta(sel)
    end if
    cc=matmul(xcut,w)
    call concordance_right(time,status,cc,conc)
    res%c_index_no_cv=conc%cindex
    cccv=0
    cctest=0
    res%cvl=0
    res%rcvl1=0
    res%rcvl2=0
    allocate(keep(n))
    do k=1,kfold0
      lo=(k-1)*n/kfold0+1
      hi=k*n/kfold0
      keep=.true.
      keep(lo:hi)=.false.
      call subset(time,status,x,keep,tt,ddv,xx)
      if(score0)then
      call uni_score(tt,ddv,xx,urk,dd)
      else
      call uni_wald(tt,ddv,xx,urk)
      end if
      ns=count(urk%p<pth)
      allocate(ccall(n))
      ccall=0
      if(ns>0)then
      allocate(selk(ns))
      kk=0
      do j=1,p
      if(urk%p(j)<pth)then
      kk=kk+1
      selk(kk)=j
      end if
      end do
      if(score0)then
      ccall=matmul(x(:,selk),urk%z(selk))
      else
      ccall=matmul(x(:,selk),urk%beta(selk))
      end if
      deallocate(selk)
      end if
      cctest(lo:hi)=ccall(lo:hi)
      allocate(risktrain(count(keep)))
      risktrain=pack(cc,keep)
      call fit_scalar_cox(tt,ddv,risktrain,beta,lltrain)
      llfull=cox_scalar_full(time,status,cc,beta)
      res%rcvl1=res%rcvl1+llfull-lltrain
      deallocate(risktrain)
      if(score0)then
      call uni_score(tt,ddv,pack_matrix_rows(xcut,keep),urk,dd)
      w=urk%z
      else
      call uni_wald(tt,ddv,pack_matrix_rows(xcut,keep),urk)
      w=urk%beta
      end if
      cccv(lo:hi)=matmul(xcut(lo:hi,:),w)
      risktrain=matmul(pack_matrix_rows(xcut,keep),w)
      call fit_scalar_cox(tt,ddv,risktrain,beta,lltrain)
      llfull=cox_scalar_full(time,status,cc,beta)
      res%rcvl2=res%rcvl2+llfull-lltrain
      deallocate(risktrain)
      risktrain=pack(ccall,keep)
      call fit_scalar_cox(tt,ddv,risktrain,beta,lltrain)
      llfull=cox_scalar_full(time,status,ccall,beta)
      res%cvl=res%cvl+llfull-lltrain
      deallocate(risktrain,ccall)
    end do
    call concordance_right(time,status,cccv,conc)
    res%c_index_incomplete=conc%cindex
    call concordance_right(time,status,cctest,conc)
    res%c_index_full=conc%cindex
    res%n_genes=p
    res%n_selected=q
    allocate(res%beta(q),res%z(q),res%p(q),res%selected(q))
    res%beta=ur%beta(sel)
    res%z=ur%z(sel)
    res%p=ur%p(sel)
    res%selected=sel
    call sort_selected(res%p,res%beta,res%z,res%selected)
    res%fdr_formula=pth*real(p,dp)/real(q,dp)
    if(perm0)then
    allocate(permidx(n))
    permidx=[(i,i=1,n)]
    fcount=0
    do i=1,mperm
    call shuffle_int(permidx)
    if(score0)then
    call uni_score(time,status,x(permidx,:),urk,dd)
    else
    call uni_wald(time,status,x(permidx,:),urk)
    end if
    fcount=fcount+count(urk%p<pth)
    end do
    res%false_selected=fcount/real(mperm,dp)
    res%fdr_permutation=res%false_selected/real(q,dp)
    end if
  end subroutine uni_selection

  subroutine subset(time,status,x,keep,t2,d2,x2)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    logical,intent(in)::keep(:)
    real(dp),allocatable,intent(out)::t2(:),x2(:,:)
    integer,allocatable,intent(out)::d2(:)
    integer::m,i,k
    m=count(keep)
    allocate(t2(m),d2(m),x2(m,size(x,2)))
    k=0
    do i=1,size(time)
    if(keep(i))then
    k=k+1
    t2(k)=time(i)
    d2(k)=status(i)
    x2(k,:)=x(i,:)
    end if
    end do
  end subroutine subset
  function pack_matrix_rows(x,keep) result(y)
    real(dp),intent(in)::x(:,:)
    logical,intent(in)::keep(:)
    real(dp),allocatable::y(:,:)
    integer::i,k
    allocate(y(count(keep),size(x,2)))
    k=0
    do i=1,size(x,1)
    if(keep(i))then
    k=k+1
    y(k,:)=x(i,:)
    end if
    end do
  end function pack_matrix_rows
  real(dp) function cox_scalar_full(time,status,risk,beta) result(ll)
    real(dp),intent(in)::time(:),risk(:),beta
    integer,intent(in)::status(:)
    real(dp)::xx(size(time),1),bb(1)
    xx(:,1)=risk
    bb(1)=beta
    ll=cox_loglik_breslow(time,status,xx,bb)
  end function cox_scalar_full
  subroutine sort_selected(p,b,z,idx)
    real(dp),intent(inout)::p(:),b(:),z(:)
    integer,intent(inout)::idx(:)
    integer::i,j,iv
    real(dp)::pv,bv,zv
    do i=2,size(p)
    pv=p(i)
    bv=b(i)
    zv=z(i)
    iv=idx(i)
    j=i-1
    do while(j>=1)
    if(p(j)<=pv)exit
    p(j+1)=p(j)
    b(j+1)=b(j)
    z(j+1)=z(j)
    idx(j+1)=idx(j)
    j=j-1
    end do
    p(j+1)=pv
    b(j+1)=bv
    z(j+1)=zv
    idx(j+1)=iv
    end do
  end subroutine sort_selected
end module compound_cox_selection
