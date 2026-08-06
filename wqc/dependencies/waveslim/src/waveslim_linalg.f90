! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_linalg
  use waveslim_kinds, only : dp
  implicit none
  private
  public :: symmetric_eigen, solve_linear, inverse_matrix
contains
  subroutine symmetric_eigen(a,values,vectors,info)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::values(:),vectors(:,:)
    integer,intent(out),optional::info
    real(dp),allocatable::b(:,:)
    real(dp)::app,aqq,apq,tau,t,c,s,tmp,maxoff
    integer::n,p,q,i,j,iter,maxiter,k
    n=size(a,1)
    allocate(b(n,n),vectors(n,n),values(n))
    b=0.5_dp*(a+transpose(a))
    vectors=0.0_dp
    do i=1,n
    vectors(i,i)=1.0_dp
    end do
    maxiter=max(100,50*n*n)
    do iter=1,maxiter
      maxoff=0.0_dp
      p=1
      q=min(2,n)
      do i=1,n-1
      do j=i+1,n
      if(abs(b(i,j))>maxoff)then
      maxoff=abs(b(i,j))
      p=i
      q=j
      end if
      end do
      end do
      if(maxoff<1e-13_dp*max(1.0_dp,maxval(abs(b))))exit
      app=b(p,p)
      aqq=b(q,q)
      apq=b(p,q)
      tau=(aqq-app)/(2.0_dp*apq)
      t=sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau))
      c=1.0_dp/sqrt(1.0_dp+t*t)
      s=t*c
      do k=1,n
        if(k/=p.and.k/=q)then
        tmp=b(k,p)
        b(k,p)=c*tmp-s*b(k,q)
        b(p,k)=b(k,p)
        b(k,q)=s*tmp+c*b(k,q)
        b(q,k)=b(k,q)
        end if
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
      b(p,q)=0.0_dp
      b(q,p)=0.0_dp
      do k=1,n
      tmp=vectors(k,p)
      vectors(k,p)=c*tmp-s*vectors(k,q)
      vectors(k,q)=s*tmp+c*vectors(k,q)
      end do
    end do
    do i=1,n
    values(i)=b(i,i)
    end do
    call sort_eigen(values,vectors)
    if(present(info))info=merge(0,1,iter<=maxiter)
  end subroutine symmetric_eigen

  subroutine sort_eigen(values,vectors)
    real(dp),intent(inout)::values(:),vectors(:,:)
    integer::i,j,k,n
    real(dp)::tmp
    real(dp),allocatable::col(:)
    n=size(values)
    allocate(col(n))
    do i=1,n-1
    k=i
    do j=i+1,n
    if(values(j)>values(k))k=j
    end do
    if(k/=i)then
    tmp=values(i)
    values(i)=values(k)
    values(k)=tmp
    col=vectors(:,i)
    vectors(:,i)=vectors(:,k)
    vectors(:,k)=col
    end if
    end do
  end subroutine sort_eigen

  subroutine solve_linear(a,b,x,info)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),allocatable,intent(out)::x(:)
    integer,intent(out),optional::info
    real(dp),allocatable::m(:,:),rhs(:),row(:)
    real(dp)::factor,piv
    integer::n,i,k,p
    n=size(b)
    allocate(m(n,n),rhs(n),x(n),row(n))
    m=a
    rhs=b
    do k=1,n
      p=k
      do i=k+1,n
      if(abs(m(i,k))>abs(m(p,k)))p=i
      end do
      if(abs(m(p,k))<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(m))))then
      x=0.0_dp
      if(present(info))info=1
      return
      end if
      if(p/=k)then
      row=m(k,:)
      m(k,:)=m(p,:)
      m(p,:)=row
      piv=rhs(k)
      rhs(k)=rhs(p)
      rhs(p)=piv
      end if
      do i=k+1,n
      factor=m(i,k)/m(k,k)
      m(i,k:n)=m(i,k:n)-factor*m(k,k:n)
      rhs(i)=rhs(i)-factor*rhs(k)
      end do
    end do
    do i=n,1,-1
    x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
    end do
    if(present(info))info=0
  end subroutine solve_linear

  subroutine inverse_matrix(a,ainv,info)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::ainv(:,:)
    integer,intent(out),optional::info
    real(dp),allocatable::e(:),x(:)
    integer::n,j,st
    n=size(a,1)
    allocate(ainv(n,n),e(n))
    ainv=0.0_dp
    st=0
    do j=1,n
    e=0.0_dp
    e(j)=1.0_dp
    call solve_linear(a,e,x,st)
    if(st/=0)exit
    ainv(:,j)=x
    end do
    if(present(info))info=st
  end subroutine inverse_matrix
end module waveslim_linalg
