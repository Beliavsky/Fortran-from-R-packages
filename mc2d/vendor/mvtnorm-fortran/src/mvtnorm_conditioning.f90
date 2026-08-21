! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_conditioning
  use mvtnorm_kinds, only : dp
  use mvtnorm_types, only : conditional_result, mvnormal_model
  use mvtnorm_linalg, only : inverse_spd, solve_spd, cholesky_lower, symmetrize
  use mvtnorm_distributions, only : rmvnorm
  implicit none
  private
  public :: make_mvnormal, marginal_mvnormal, conditional_mvnormal
  public :: simulate_mvnormal, permute_mvnormal

contains

  function make_mvnormal(mean,covariance) result(model)
    real(dp),intent(in)::mean(:),covariance(:,:)
    type(mvnormal_model)::model
    logical::ok
    character(len=256)::message
    allocate(model%mean(size(mean))); model%mean=mean
    if(size(covariance,1)/=size(mean) .or. size(covariance,2)/=size(mean)) then
      model%message='non-conforming dimensions'; return
    end if
    call cholesky_lower(covariance,model%chol,ok,message)
    model%ok=ok; model%message=message
  end function make_mvnormal

  function marginal_mvnormal(mean,covariance,which) result(out)
    real(dp),intent(in)::mean(:),covariance(:,:)
    integer,intent(in)::which(:)
    type(conditional_result)::out
    integer::i,j,n
    n=size(which); allocate(out%mean(n),out%covariance(n,n))
    if(any(which<1) .or. any(which>size(mean))) then
      out%message='marginal indices out of range'; return
    end if
    do i=1,n
      out%mean(i)=mean(which(i))
      do j=1,n
        out%covariance(i,j)=covariance(which(i),which(j))
      end do
    end do
    call symmetrize(out%covariance); out%ok=.true.
  end function marginal_mvnormal

  function conditional_mvnormal(mean,covariance,which_given,given) result(out)
    real(dp),intent(in)::mean(:),covariance(:,:),given(:)
    integer,intent(in)::which_given(:)
    type(conditional_result)::out
    integer,allocatable::which_free(:)
    real(dp),allocatable::s11(:,:),s12(:,:),s21(:,:),s22(:,:),rhs(:,:),sol(:,:),gain(:,:)
    logical::ok
    character(len=256)::message
    integer::p,ng,nf,i,j,k

    p=size(mean); ng=size(which_given); nf=p-ng
    if(size(given)/=ng .or. any(which_given<1) .or. any(which_given>p)) then
      out%message='invalid conditioning dimensions or indices'; return
    end if
    allocate(which_free(nf)); k=0
    do i=1,p
      if(.not.any(which_given==i)) then; k=k+1; which_free(k)=i; end if
    end do
    allocate(s11(ng,ng),s12(ng,nf),s21(nf,ng),s22(nf,nf))
    do i=1,ng
      do j=1,ng; s11(i,j)=covariance(which_given(i),which_given(j)); end do
      do j=1,nf; s12(i,j)=covariance(which_given(i),which_free(j)); end do
    end do
    do i=1,nf
      do j=1,ng; s21(i,j)=covariance(which_free(i),which_given(j)); end do
      do j=1,nf; s22(i,j)=covariance(which_free(i),which_free(j)); end do
    end do
    call solve_spd(s11,s12,sol,ok,message)
    if(.not.ok) then; out%message='conditioning covariance is singular'; return; end if
    gain=transpose(sol)
    allocate(rhs(ng,1)); rhs(:,1)=given-mean(which_given)
    allocate(out%mean(nf)); out%mean=mean(which_free)+matmul(gain,rhs(:,1))
    allocate(out%covariance(nf,nf)); out%covariance=s22-matmul(gain,s12)
    call symmetrize(out%covariance); out%ok=.true.
  end function conditional_mvnormal

  function simulate_mvnormal(model,n,seed,standardize) result(x)
    type(mvnormal_model),intent(in)::model
    integer,intent(in)::n
    integer,intent(in),optional::seed
    logical,intent(in),optional::standardize
    real(dp),allocatable::x(:,:),cov(:,:)
    logical::std
    integer::i,p
    real(dp),allocatable::sd(:)
    std=.false.; if(present(standardize)) std=standardize
    p=size(model%mean); cov=matmul(model%chol,transpose(model%chol))
    if(std) then
      allocate(sd(p)); do i=1,p; sd(i)=sqrt(cov(i,i)); end do
      do i=1,p
        cov(i,:)=cov(i,:)/sd(i)
        cov(:,i)=cov(:,i)/sd(i)
      end do
      x=rmvnorm(n,spread(0.0_dp,1,p),cov,seed)
    else
      x=rmvnorm(n,model%mean,cov,seed)
    end if
  end function simulate_mvnormal

  function permute_mvnormal(model,perm) result(out)
    type(mvnormal_model),intent(in)::model
    integer,intent(in)::perm(:)
    type(mvnormal_model)::out
    real(dp),allocatable::cov(:,:),pcov(:,:)
    integer::i,j,n
    logical::ok
    character(len=256)::message
    n=size(perm); allocate(out%mean(n),pcov(n,n)); cov=matmul(model%chol,transpose(model%chol))
    do i=1,n
      out%mean(i)=model%mean(perm(i))
      do j=1,n; pcov(i,j)=cov(perm(i),perm(j)); end do
    end do
    call cholesky_lower(pcov,out%chol,ok,message); out%ok=ok; out%message=message
  end function permute_mvnormal
end module mvtnorm_conditioning
