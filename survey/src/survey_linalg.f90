! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_linalg
  use survey_kinds, only : dp
  implicit none
  private
  public :: solve_linear, inverse_matrix, sym_pinv, symmetric_eigen, matrix_rank_sym, outer_product
contains

  pure function outer_product(x,y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    a = spread(x,2,size(y))*spread(y,1,size(x))
  end function outer_product

  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:), bb(:), rowtmp(:)
    real(dp) :: piv, fac, tmp
    integer :: n,i,j,k,p
    n=size(b); info=0
    if(size(a,1)/=n .or. size(a,2)/=n .or. size(x)/=n) then
      info=-1; x=0; return
    end if
    allocate(aa(n,n),bb(n),rowtmp(n)); aa=a; bb=b
    do k=1,n-1
      p=k
      do i=k+1,n
        if(abs(aa(i,k))>abs(aa(p,k))) p=i
      end do
      if(abs(aa(p,k)) <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa)))) then
        info=k; x=0; return
      end if
      if(p/=k) then
        rowtmp=aa(k,:); aa(k,:)=aa(p,:); aa(p,:)=rowtmp
        tmp=bb(k); bb(k)=bb(p); bb(p)=tmp
      end if
      do i=k+1,n
        fac=aa(i,k)/aa(k,k)
        aa(i,k)=0.0_dp
        do j=k+1,n
          aa(i,j)=aa(i,j)-fac*aa(k,j)
        end do
        bb(i)=bb(i)-fac*bb(k)
      end do
    end do
    if(abs(aa(n,n)) <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa)))) then
      info=n; x=0; return
    end if
    x(n)=bb(n)/aa(n,n)
    do i=n-1,1,-1
      piv=bb(i)-dot_product(aa(i,i+1:n),x(i+1:n))
      x(i)=piv/aa(i,i)
    end do
  end subroutine solve_linear

  subroutine inverse_matrix(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: e(:), x(:)
    integer :: n,j,istat
    n=size(a,1); info=0
    if(size(a,2)/=n .or. size(ainv,1)/=n .or. size(ainv,2)/=n) then
      info=-1; ainv=0; return
    end if
    allocate(e(n),x(n)); ainv=0
    do j=1,n
      e=0; e(j)=1
      call solve_linear(a,e,x,istat)
      if(istat/=0) then
        info=istat; ainv=0; return
      end if
      ainv(:,j)=x
    end do
  end subroutine inverse_matrix

  subroutine symmetric_eigen(a, eval, evec, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: eval(:), evec(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: b(:,:)
    real(dp) :: app,aqq,apq,tau,t,c,s,bkp,bkq,vkp,vkq,tol,off
    integer :: n,p,q,k,iter,maxit,i,j,imax
    real(dp) :: tmp
    real(dp), allocatable :: coltmp(:)
    n=size(a,1); info=0
    if(size(a,2)/=n .or. size(eval)/=n .or. size(evec,1)/=n .or. size(evec,2)/=n) then
      info=-1; return
    end if
    allocate(b(n,n),coltmp(n)); b=0.5_dp*(a+transpose(a)); evec=0
    do i=1,n; evec(i,i)=1.0_dp; end do
    tol=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(b)))
    maxit=max(50,50*n*n)
    do iter=1,maxit
      off=0; p=1; q=min(2,n)
      do i=1,n-1
        do j=i+1,n
          if(abs(b(i,j))>off) then
            off=abs(b(i,j)); p=i; q=j
          end if
        end do
      end do
      if(n==1 .or. off<=tol) exit
      app=b(p,p); aqq=b(q,q); apq=b(p,q)
      tau=(aqq-app)/(2.0_dp*apq)
      if(tau>=0) then
        t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
      do k=1,n
        if(k/=p .and. k/=q) then
          bkp=b(k,p); bkq=b(k,q)
          b(k,p)=c*bkp-s*bkq; b(p,k)=b(k,p)
          b(k,q)=s*bkp+c*bkq; b(q,k)=b(k,q)
        end if
      end do
      b(p,p)=c*c*app-2*c*s*apq+s*s*aqq
      b(q,q)=s*s*app+2*c*s*apq+c*c*aqq
      b(p,q)=0; b(q,p)=0
      do k=1,n
        vkp=evec(k,p); vkq=evec(k,q)
        evec(k,p)=c*vkp-s*vkq
        evec(k,q)=s*vkp+c*vkq
      end do
    end do
    if(iter>maxit) info=1
    do i=1,n; eval(i)=b(i,i); end do
    ! descending order
    do i=1,n-1
      imax=i
      do j=i+1,n
        if(eval(j)>eval(imax)) imax=j
      end do
      if(imax/=i) then
        tmp=eval(i); eval(i)=eval(imax); eval(imax)=tmp
        coltmp=evec(:,i); evec(:,i)=evec(:,imax); evec(:,imax)=coltmp
      end if
    end do
  end subroutine symmetric_eigen

  subroutine sym_pinv(a,ainv,rank,tol,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: rank, info
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: eval(:), evec(:,:)
    real(dp) :: cutoff
    integer :: n,j
    n=size(a,1); allocate(eval(n),evec(n,n)); ainv=0; rank=0
    call symmetric_eigen(a,eval,evec,info)
    if(info<0) return
    if(present(tol)) then
      cutoff=tol
    else
      cutoff=256.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(eval)))
    end if
    do j=1,n
      if(abs(eval(j))>cutoff) then
        ainv=ainv + outer_product(evec(:,j),evec(:,j))/eval(j)
        rank=rank+1
      end if
    end do
  end subroutine sym_pinv

  integer function matrix_rank_sym(a,tol) result(r)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: eval(:), evec(:,:)
    real(dp) :: cutoff
    integer :: n,info
    n=size(a,1); allocate(eval(n),evec(n,n))
    call symmetric_eigen(a,eval,evec,info)
    if(present(tol)) then
      cutoff=tol
    else
      cutoff=256.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(eval)))
    end if
    r=count(abs(eval)>cutoff)
  end function matrix_rank_sym
end module survey_linalg
