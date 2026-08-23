module mixspe_linalg
  use mvtnorm_kinds, only : dp
  use mvtnorm_linalg, only : jacobi_eigen, cholesky_lower, solve_lower
  implicit none
  private
  public :: sym_eigen, sym_power, mahalanobis2, covariance_weighted, inverse_sym, determinant_sym
contains
  subroutine sym_eigen(a, vals, vecs, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: vals(:), vecs(:,:)
    logical, intent(out) :: ok
    integer :: i,j,n,k
    real(dp) :: tv
    real(dp), allocatable :: tc(:)
    call jacobi_eigen(a,vals,vecs,ok)
    if (.not.ok) return
    n=size(vals); allocate(tc(n))
    do i=1,n-1
      k=i
      do j=i+1,n
        if (vals(j)>vals(k)) k=j
      end do
      if (k/=i) then
        tv=vals(i); vals(i)=vals(k); vals(k)=tv
        tc=vecs(:,i); vecs(:,i)=vecs(:,k); vecs(:,k)=tc
      end if
    end do
  end subroutine

  function sym_power(a,pow,ok) result(b)
    real(dp), intent(in) :: a(:,:), pow
    logical, intent(out), optional :: ok
    real(dp), allocatable :: b(:,:), vals(:), vecs(:,:)
    logical :: lok
    integer :: i,n
    call sym_eigen(a,vals,vecs,lok)
    n=size(a,1); allocate(b(n,n)); b=0.0_dp
    if (.not.lok) then
      if(present(ok)) ok=.false.; return
    end if
    do i=1,n
      if (vals(i)<=0.0_dp .and. pow<0.0_dp) then
        if(present(ok)) ok=.false.; return
      end if
      b=b+max(abs(vals(i)),tiny(1.0_dp))**pow * spread(vecs(:,i),2,n)*spread(vecs(:,i),1,n)
    end do
    if(present(ok)) ok=.true.
  end function

  real(dp) function mahalanobis2(x,mu,sigma,ok) result(q)
    real(dp),intent(in)::x(:),mu(:),sigma(:,:)
    logical,intent(out),optional::ok
    real(dp),allocatable::l(:,:),rhs(:,:),sol(:,:)
    logical::lok
    character(len=256)::msg
    integer::p
    p=size(x); allocate(rhs(p,1)); rhs(:,1)=x-mu
    call cholesky_lower(sigma,l,lok,msg)
    if(.not.lok) then; q=huge(1.0_dp); if(present(ok))ok=.false.; return; end if
    call solve_lower(l,rhs,sol,lok)
    q=sum(sol(:,1)**2); if(present(ok))ok=lok
  end function

  function inverse_sym(a,ok) result(inv)
    real(dp),intent(in)::a(:,:)
    logical,intent(out),optional::ok
    real(dp),allocatable::inv(:,:),vals(:),vecs(:,:)
    logical::lok
    integer::i,n
    call sym_eigen(a,vals,vecs,lok); n=size(a,1); allocate(inv(n,n)); inv=0.0_dp
    if(lok) then
      do i=1,n
        if(vals(i)<=tiny(1.0_dp)) then; lok=.false.; exit; end if
        inv=inv+(1.0_dp/vals(i))*spread(vecs(:,i),2,n)*spread(vecs(:,i),1,n)
      end do
    end if
    if(present(ok))ok=lok
  end function

  real(dp) function determinant_sym(a,ok) result(d)
    real(dp),intent(in)::a(:,:)
    logical,intent(out),optional::ok
    real(dp),allocatable::vals(:),vecs(:,:)
    logical::lok
    call sym_eigen(a,vals,vecs,lok)
    if(lok) then; d=product(vals); else; d=0.0_dp; end if
    if(present(ok))ok=lok
  end function

  function covariance_weighted(x,w,mu) result(s)
    real(dp),intent(in)::x(:,:),w(:),mu(:)
    real(dp),allocatable::s(:,:)
    real(dp)::sw
    integer::i,p
    p=size(x,2); allocate(s(p,p)); s=0.0_dp; sw=sum(w)
    do i=1,size(x,1)
      s=s+w(i)*spread(x(i,:)-mu,2,p)*spread(x(i,:)-mu,1,p)
    end do
    s=s/max(sw,tiny(1.0_dp))
  end function
end module
