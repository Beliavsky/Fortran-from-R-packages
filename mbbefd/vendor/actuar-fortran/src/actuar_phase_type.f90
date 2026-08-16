! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_phase_type
  use actuar_kinds, only : dp
  use actuar_special, only : nan_dp
  use actuar_rng, only : runif, rexp
  implicit none
  private
  public :: dphtype, pphtype, rphtype, mphtype, mgfphtype
  public :: matrix_exponential
contains

  function matrix_exponential(a) result(ea)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable :: ea(:,:)
    real(dp), allocatable :: b(:,:),term(:,:),tmp(:,:),ident(:,:)
    real(dp) :: norm1
    integer :: n,i,k,s
    n=size(a,1)
    allocate(ea(n,n),b(n,n),term(n,n),tmp(n,n),ident(n,n))
    ident=0.0_dp
    do i=1,n; ident(i,i)=1.0_dp; end do
    norm1=maxval(sum(abs(a),dim=1))
    s=max(0,ceiling(log(max(1.0_dp,norm1))/log(2.0_dp)))
    b=a/(2.0_dp**s)
    ea=ident; term=ident
    do k=1,80
      tmp=matmul(term,b)/real(k,dp)
      term=tmp; ea=ea+term
      if(maxval(abs(term))<1.0e-16_dp*max(1.0_dp,maxval(abs(ea)))) exit
    end do
    do k=1,s
      ea=matmul(ea,ea)
    end do
  end function matrix_exponential

  pure function valid_phase(pi0,tmat) result(ok)
    real(dp), intent(in) :: pi0(:),tmat(:,:)
    logical :: ok
    integer :: i,j,n
    n=size(pi0)
    ok=size(tmat,1)==n .and. size(tmat,2)==n .and. &
       all(pi0>=0.0_dp) .and. sum(pi0)<=1.0_dp+1.0e-12_dp
    if(.not.ok) return
    do i=1,n
      if(tmat(i,i)>=0.0_dp) ok=.false.
      do j=1,n
        if(i/=j .and. tmat(i,j)<0.0_dp) ok=.false.
      end do
      if(sum(tmat(i,:))>1.0e-12_dp) ok=.false.
    end do
  end function valid_phase

  function dphtype(x,pi0,tmat) result(d)
    real(dp), intent(in) :: x,pi0(:),tmat(:,:)
    real(dp) :: d
    real(dp), allocatable :: e(:,:),exitv(:)
    if(.not.valid_phase(pi0,tmat)) then
      d=nan_dp(); return
    end if
    if(x<0.0_dp) then
      d=0.0_dp; return
    end if
    e=matrix_exponential(tmat*x)
    allocate(exitv(size(pi0))); exitv=-sum(tmat,dim=2)
    d=dot_product(pi0,matmul(e,exitv))
    d=max(0.0_dp,d)
  end function dphtype

  function pphtype(x,pi0,tmat) result(p)
    real(dp), intent(in) :: x,pi0(:),tmat(:,:)
    real(dp) :: p
    real(dp), allocatable :: e(:,:),ones(:)
    if(.not.valid_phase(pi0,tmat)) then
      p=nan_dp(); return
    end if
    if(x<0.0_dp) then
      p=0.0_dp; return
    end if
    e=matrix_exponential(tmat*x)
    allocate(ones(size(pi0))); ones=1.0_dp
    p=1.0_dp-dot_product(pi0,matmul(e,ones))
    p=max(0.0_dp,min(1.0_dp,p))
  end function pphtype

  function rphtype(pi0,tmat) result(x)
    real(dp), intent(in) :: pi0(:),tmat(:,:)
    real(dp) :: x,u,cum,rate
    integer :: state,j,n
    real(dp), allocatable :: probs(:)
    if(.not.valid_phase(pi0,tmat)) then
      x=nan_dp(); return
    end if
    n=size(pi0); allocate(probs(n+1))
    probs(:n)=pi0; probs(n+1)=max(0.0_dp,1.0_dp-sum(pi0))
    u=runif(); cum=0.0_dp; state=n+1
    do j=1,n+1
      cum=cum+probs(j)
      if(u<=cum) then; state=j; exit; end if
    end do
    x=0.0_dp
    do while(state<=n)
      rate=-tmat(state,state)
      x=x+rexp(rate)
      probs=0.0_dp
      do j=1,n
        if(j/=state) probs(j)=tmat(state,j)/rate
      end do
      probs(n+1)=max(0.0_dp,-sum(tmat(state,:))/rate)
      u=runif(); cum=0.0_dp
      do j=1,n+1
        cum=cum+probs(j)
        if(u<=cum) then; state=j; exit; end if
      end do
    end do
  end function rphtype

  function solve_matrix(a,b,ok) result(x)
    real(dp), intent(in) :: a(:,:),b(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: aug(:,:),row(:)
    real(dp) :: pivot,factor
    integer :: n,m,i,k,p
    n=size(a,1); m=size(b,2)
    allocate(aug(n,n+m),x(n,m),row(n+m))
    aug(:,:n)=a; aug(:,n+1:)=b; ok=.true.
    do k=1,n
      p=k
      do i=k+1,n
        if(abs(aug(i,k))>abs(aug(p,k))) p=i
      end do
      if(abs(aug(p,k))<1.0e-14_dp) then
        ok=.false.; x=0.0_dp; return
      end if
      if(p/=k) then; row=aug(k,:); aug(k,:)=aug(p,:); aug(p,:)=row; end if
      pivot=aug(k,k); aug(k,:)=aug(k,:)/pivot
      do i=1,n
        if(i==k) cycle
        factor=aug(i,k); aug(i,:)=aug(i,:)-factor*aug(k,:)
      end do
    end do
    x=aug(:,n+1:)
  end function solve_matrix

  function mphtype(order,pi0,tmat) result(m)
    integer, intent(in) :: order
    real(dp), intent(in) :: pi0(:),tmat(:,:)
    real(dp) :: m
    real(dp), allocatable :: rhs(:,:),sol(:,:),ones(:)
    logical :: ok
    integer :: j,n
    if(order<0 .or. .not.valid_phase(pi0,tmat)) then
      m=nan_dp(); return
    end if
    if(order==0) then; m=1.0_dp; return; end if
    n=size(pi0); allocate(ones(n)); ones=1.0_dp
    allocate(rhs(n,1)); rhs(:,1)=ones
    sol=rhs
    do j=1,order
      sol=solve_matrix(-tmat,sol,ok)
      if(.not.ok) then; m=nan_dp(); return; end if
    end do
    m=gamma(real(order+1,dp))*dot_product(pi0,sol(:,1))
  end function mphtype

  function mgfphtype(s,pi0,tmat) result(m)
    real(dp), intent(in) :: s,pi0(:),tmat(:,:)
    real(dp) :: m
    real(dp), allocatable :: a(:,:),rhs(:,:),sol(:,:),exitv(:)
    logical :: ok
    integer :: n,i
    if(.not.valid_phase(pi0,tmat)) then
      m=nan_dp(); return
    end if
    n=size(pi0); allocate(a(n,n),rhs(n,1),exitv(n))
    a=-tmat
    do i=1,n; a(i,i)=a(i,i)-s; end do
    exitv=-sum(tmat,dim=2); rhs(:,1)=exitv
    sol=solve_matrix(a,rhs,ok)
    if(.not.ok) then
      m=nan_dp()
    else
      m=dot_product(pi0,sol(:,1))+max(0.0_dp,1.0_dp-sum(pi0))
    end if
  end function mgfphtype

end module actuar_phase_type
