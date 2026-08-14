! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn_learning
  use, intrinsic :: iso_fortran_env, only : int64
  use fnn_kinds, only : dp
  use fnn_types, only : knn_result, classification_result, regression_result, ownn_result
  use fnn_neighbors, only : get_knn, get_knnx
  implicit none
  private
  public :: knn_classify, knn_cv, knn_reg, ownn
contains
  function knn_classify(train,test,cl,k,algorithm) result(out)
    real(dp), target, intent(in) :: train(:,:),test(:,:)
    integer, intent(in) :: cl(:),k
    character(len=*), intent(in), optional :: algorithm
    type(classification_result) :: out
    type(knn_result) :: z
    integer :: i,keff
    if(size(train,1)/=size(cl)) error stop "knn_classify: class length mismatch"
    if(size(train,2)/=size(test,2)) error stop "knn_classify: dimensions differ"
    keff=min(k,size(train,1)); if(keff<1) error stop "knn_classify: invalid k"
    z=get_knnx(train,test,keff,algorithm)
    allocate(out%class(size(test,1)),out%probability(size(test,1)))
    allocate(out%nn_index(size(z%index,1),size(z%index,2)))
    allocate(out%nn_distance(size(z%distance,1),size(z%distance,2)))
    out%nn_index=z%index
    out%nn_distance=z%distance
    do i=1,size(test,1)
      call weighted_vote(cl(z%index(i,:)),spread(1.0_dp,1,keff),out%class(i),out%probability(i))
    end do
  end function knn_classify

  function knn_cv(train,cl,k,algorithm) result(out)
    real(dp), target, intent(in) :: train(:,:)
    integer, intent(in) :: cl(:),k
    character(len=*), intent(in), optional :: algorithm
    type(classification_result) :: out
    type(knn_result) :: z
    integer :: i,keff
    if(size(train,1)/=size(cl)) error stop "knn_cv: class length mismatch"
    keff=min(k,size(train,1)-1); if(keff<1) error stop "knn_cv: invalid k"
    z=get_knn(train,keff,algorithm)
    allocate(out%class(size(train,1)),out%probability(size(train,1)))
    allocate(out%nn_index(size(z%index,1),size(z%index,2)))
    allocate(out%nn_distance(size(z%distance,1),size(z%distance,2)))
    out%nn_index=z%index
    out%nn_distance=z%distance
    do i=1,size(train,1)
      call weighted_vote(cl(z%index(i,:)),spread(1.0_dp,1,keff),out%class(i),out%probability(i))
    end do
  end function knn_cv

  function knn_reg(train,y,k,test,algorithm) result(out)
    real(dp), target, intent(in) :: train(:,:)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: k
    real(dp), target, intent(in), optional :: test(:,:)
    character(len=*), intent(in), optional :: algorithm
    type(regression_result) :: out
    type(knn_result) :: z
    integer :: i
    real(dp) :: denom
    if(size(train,1)/=size(y)) error stop "knn_reg: y length mismatch"
    if(present(test)) then
      z=get_knnx(train,test,min(k,size(train,1)),algorithm)
      allocate(out%prediction(size(test,1)))
    else
      if(k>=size(train,1)) error stop "knn_reg: k must be < n for CV"
      z=get_knn(train,k,algorithm)
      allocate(out%prediction(size(train,1)))
      out%cross_validated=.true.
    end if
    do i=1,size(out%prediction)
      out%prediction(i)=sum(y(z%index(i,:)))/real(size(z%index,2),dp)
    end do
    if(.not.present(test)) then
      allocate(out%residuals(size(y)))
      out%residuals=y-out%prediction
      out%press=sum(out%residuals**2)
      denom=sum((y-sum(y)/real(size(y),dp))**2)
      if(denom>0.0_dp) then
        out%r2_predict=1.0_dp-out%press/denom
      else
        out%r2_predict=0.0_dp
      end if
    end if
  end function knn_reg

  function ownn(train,test,cl,k,algorithm,seed,testcl) result(out)
    real(dp), target, intent(in) :: train(:,:),test(:,:)
    integer, intent(in) :: cl(:)
    integer, intent(in), optional :: k,seed
    integer, intent(in), optional :: testcl(:)
    character(len=*), intent(in), optional :: algorithm
    type(ownn_result) :: out
    type(knn_result) :: z
    integer :: kval,n,d,i,kstar
    real(dp), allocatable :: wknn(:),wownn(:),wbnn(:)
    real(dp) :: q,alpha
    n=size(train,1); d=size(train,2)
    if(size(cl)/=n) error stop "ownn: class length mismatch"
    if(size(test,2)/=d) error stop "ownn: dimensions differ"
    if(present(k)) then; kval=max(1,min(k,n))
    else; kval=choose_k(train,cl,algorithm,seed); end if
    out%k=kval
    z=get_knnx(train,test,n,algorithm)
    allocate(wknn(n),wownn(n),wbnn(n)); wknn=0.0_dp; wownn=0.0_dp; wbnn=0.0_dp
    wknn(1:kval)=1.0_dp/real(kval,dp)
    kstar=floor((2.0_dp*real(d+4,dp)/real(d+2,dp))**(real(d,dp)/real(d+4,dp))*real(kval,dp))
    kstar=max(1,min(kstar,n))
    do i=1,kstar
      alpha=real(i,dp)**(1.0_dp+2.0_dp/real(d,dp))-real(i-1,dp)**(1.0_dp+2.0_dp/real(d,dp))
      wownn(i)=(1.0_dp+0.5_dp*real(d,dp)-real(d,dp)*alpha/(2.0_dp*real(kstar,dp)**(2.0_dp/real(d,dp))))/real(kstar,dp)
    end do
    q=2.0_dp**(real(d,dp)/real(d+4,dp))*gamma(2.0_dp+2.0_dp/real(d,dp))**(2.0_dp*real(d,dp)/real(d+4,dp))/real(kval,dp)
    q=min(max(q,epsilon(1.0_dp)),1.0_dp)
    do i=1,n
      wbnn(i)=q*(1.0_dp-q)**real(i-1,dp)/(1.0_dp-(1.0_dp-q)**real(n,dp))
    end do
    allocate(out%knn_class(size(test,1)),out%ownn_class(size(test,1)),out%bnn_class(size(test,1)))
    allocate(out%knn_probability(size(test,1)),out%ownn_probability(size(test,1)),out%bnn_probability(size(test,1)))
    do i=1,size(test,1)
      call weighted_vote(cl(z%index(i,:)),wknn,out%knn_class(i),out%knn_probability(i))
      call weighted_vote(cl(z%index(i,:)),wownn,out%ownn_class(i),out%ownn_probability(i))
      call weighted_vote(cl(z%index(i,:)),wbnn,out%bnn_class(i),out%bnn_probability(i))
    end do
    if(present(testcl)) then
      if(size(testcl)/=size(test,1)) error stop "ownn: test class length mismatch"
      out%accuracy(1)=real(count(out%knn_class==testcl),dp)/real(size(testcl),dp)
      out%accuracy(2)=real(count(out%ownn_class==testcl),dp)/real(size(testcl),dp)
      out%accuracy(3)=real(count(out%bnn_class==testcl),dp)/real(size(testcl),dp)
      out%has_accuracy=.true.
    end if
  end function ownn

  subroutine weighted_vote(labels,weights,winner,prob)
    integer, intent(in) :: labels(:)
    real(dp), intent(in) :: weights(:)
    integer, intent(out) :: winner
    real(dp), intent(out) :: prob
    integer, allocatable :: u(:)
    real(dp), allocatable :: score(:)
    integer :: i,j,nu,best
    allocate(u(size(labels)),score(size(labels))); nu=0; score=0.0_dp
    do i=1,size(labels)
      j=0
      if(nu>0) then
        do best=1,nu
          if(u(best)==labels(i)) then; j=best; exit; end if
        end do
      end if
      if(j==0) then
        nu=nu+1; u(nu)=labels(i); j=nu
      end if
      score(j)=score(j)+weights(i)
    end do
    best=1
    do i=2,nu
      if(score(i)>score(best)) then
        best=i
      else if(score(i)<score(best)) then
        cycle
      else if(u(i)<u(best)) then
        best=i
      end if
    end do
    winner=u(best)
    if(sum(score(1:nu))>0.0_dp) then; prob=score(best)/sum(score(1:nu)); else; prob=0.0_dp; end if
  end subroutine weighted_vote

  integer function choose_k(train,cl,algorithm,seed) result(kchosen)
    real(dp), target, intent(in) :: train(:,:)
    integer, intent(in) :: cl(:)
    character(len=*), intent(in), optional :: algorithm
    integer, intent(in), optional :: seed
    integer :: n,d,g,i,jj,kcand,bestk,besterr,err,lo,hi,s
    integer, allocatable :: fold(:),train_idx(:),test_idx(:)
    type(classification_result) :: cv
    real(dp), allocatable, target :: tr(:,:),te(:,:)
    integer, allocatable :: yc(:)
    n=size(train,1); d=size(train,2)
    if(n<6) then; kchosen=max(1,n/2); return; end if
    s=13579; if(present(seed)) s=seed
    allocate(fold(n)); call make_folds(fold,s)
    lo=min(5,max(1,n/10)); hi=max(lo,n/2); besterr=huge(1); bestk=lo
    do g=0,20
      kcand=nint(real(20-g,dp)*real(lo,dp)/20.0_dp+real(g,dp)*real(hi,dp)/20.0_dp)
      err=0
      do i=1,5
        train_idx=pack([(jj,jj=1,n)],fold/=i); test_idx=pack([(jj,jj=1,n)],fold==i)
        if(size(test_idx)==0 .or. size(train_idx)==0) cycle
        allocate(tr(size(train_idx),d),te(size(test_idx),d),yc(size(train_idx)))
        tr=train(train_idx,:); te=train(test_idx,:); yc=cl(train_idx)
        cv=knn_classify(tr,te,yc,min(kcand,size(train_idx)),algorithm)
        err=err+count(cv%class/=cl(test_idx))
        deallocate(tr,te,yc)
      end do
      if(err<besterr) then; besterr=err; bestk=kcand; end if
    end do
    kchosen=floor(real(bestk,dp)*(5.0_dp/4.0_dp)**(4.0_dp/real(d+4,dp)))
    kchosen=max(1,min(kchosen,n))
  end function choose_k

  subroutine make_folds(fold,seed)
    integer, intent(out) :: fold(:)
    integer, intent(in) :: seed
    integer, allocatable :: perm(:)
    integer :: i,j,tmp,n
    integer(int64) :: state
    n=size(fold); allocate(perm(n)); perm=[(i,i=1,n)]; state=abs(int(seed,int64))+1_int64
    do i=n,2,-1
      state=mod(1103515245_int64*state+12345_int64,2147483647_int64)
      j=1+int(mod(state,int(i,int64)))
      tmp=perm(i); perm(i)=perm(j); perm(j)=tmp
    end do
    do i=1,n; fold(perm(i))=1+mod(i-1,5); end do
  end subroutine make_folds
end module fnn_learning
