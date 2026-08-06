! SPDX-License-Identifier: GPL-3.0-only
module stochvoltmb_linalg
  use stochvoltmb_kinds, only : dp, tiny_dp
  use stochvoltmb_status, only : sv_ok, sv_singular
  implicit none
  private
  public :: solve_tridiagonal, tridiagonal_logdet, tridiagonal_inverse_diag
  public :: solve_linear, invert_matrix, cholesky_lower, mvn_draws

contains

  subroutine solve_tridiagonal(diag, off, rhs, x, info)
    real(dp), intent(in) :: diag(:), off(:), rhs(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: d(:), r(:)
    real(dp) :: mult
    integer :: n, i
    n=size(diag); info=sv_ok
    if (size(rhs)/=n .or. size(x)/=n .or. size(off)/=max(0,n-1)) then
      info=sv_singular; return
    end if
    allocate(d(n),r(n)); d=diag; r=rhs
    do i=2,n
      if (abs(d(i-1)) <= tiny_dp) then
        info=sv_singular; return
      end if
      mult=off(i-1)/d(i-1)
      d(i)=d(i)-mult*off(i-1)
      r(i)=r(i)-mult*r(i-1)
    end do
    if (abs(d(n)) <= tiny_dp) then
      info=sv_singular; return
    end if
    x(n)=r(n)/d(n)
    do i=n-1,1,-1
      if (abs(d(i)) <= tiny_dp) then
        info=sv_singular; return
      end if
      x(i)=(r(i)-off(i)*x(i+1))/d(i)
    end do
  end subroutine solve_tridiagonal

  real(dp) function tridiagonal_logdet(diag, off, info) result(ld)
    real(dp), intent(in) :: diag(:), off(:)
    integer, intent(out) :: info
    real(dp), allocatable :: d(:)
    integer :: n, i
    n=size(diag); info=sv_ok; ld=0.0_dp
    if (size(off)/=max(0,n-1)) then
      info=sv_singular; ld=huge(1.0_dp); return
    end if
    allocate(d(n)); d=diag
    if (d(1)<=tiny_dp) then
      info=sv_singular; ld=huge(1.0_dp); return
    end if
    ld=log(d(1))
    do i=2,n
      d(i)=d(i)-off(i-1)*off(i-1)/d(i-1)
      if (d(i)<=tiny_dp) then
        info=sv_singular; ld=huge(1.0_dp); return
      end if
      ld=ld+log(d(i))
    end do
  end function tridiagonal_logdet

  subroutine tridiagonal_inverse_diag(diag, off, invdiag, info)
    real(dp), intent(in) :: diag(:), off(:)
    real(dp), intent(out) :: invdiag(:)
    integer, intent(out) :: info
    real(dp), allocatable :: ldiag(:), lsub(:)
    integer :: n, i
    n=size(diag); info=sv_ok
    if (size(invdiag)/=n .or. size(off)/=max(0,n-1)) then
      info=sv_singular; return
    end if
    allocate(ldiag(n),lsub(n)); lsub=0.0_dp
    if (diag(1)<=tiny_dp) then
      info=sv_singular; return
    end if
    ldiag(1)=sqrt(diag(1))
    do i=2,n
      lsub(i)=off(i-1)/ldiag(i-1)
      if (diag(i)-lsub(i)*lsub(i)<=tiny_dp) then
        info=sv_singular; return
      end if
      ldiag(i)=sqrt(diag(i)-lsub(i)*lsub(i))
    end do
    invdiag(n)=1.0_dp/(ldiag(n)*ldiag(n))
    do i=n-1,1,-1
      invdiag(i)=1.0_dp/(ldiag(i)*ldiag(i)) + &
        (lsub(i+1)/ldiag(i))**2*invdiag(i+1)
    end do
  end subroutine tridiagonal_inverse_diag

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: pivot, factor, tmp
    integer :: n, i, j, k, p
    n=size(b); info=sv_ok
    if (size(a,1)/=n .or. size(a,2)/=n .or. size(x)/=n) then
      info=sv_singular; return
    end if
    allocate(aa(n,n),bb(n)); aa=a; bb=b
    do k=1,n-1
      p=k
      do i=k+1,n
        if (abs(aa(i,k))>abs(aa(p,k))) p=i
      end do
      if (abs(aa(p,k))<=tiny_dp) then
        info=sv_singular; return
      end if
      if (p/=k) then
        do j=k,n
          tmp=aa(k,j); aa(k,j)=aa(p,j); aa(p,j)=tmp
        end do
        tmp=bb(k); bb(k)=bb(p); bb(p)=tmp
      end if
      pivot=aa(k,k)
      do i=k+1,n
        factor=aa(i,k)/pivot
        aa(i,k)=0.0_dp
        aa(i,k+1:n)=aa(i,k+1:n)-factor*aa(k,k+1:n)
        bb(i)=bb(i)-factor*bb(k)
      end do
    end do
    if (abs(aa(n,n))<=tiny_dp) then
      info=sv_singular; return
    end if
    x(n)=bb(n)/aa(n,n)
    do i=n-1,1,-1
      x(i)=(bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
    end do
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: e(:), x(:)
    integer :: n, j, stat
    n=size(a,1); info=sv_ok
    if (size(a,2)/=n .or. any(shape(ainv)/=[n,n])) then
      info=sv_singular; return
    end if
    allocate(e(n),x(n))
    do j=1,n
      e=0.0_dp; e(j)=1.0_dp
      call solve_linear(a,e,x,stat)
      if (stat/=sv_ok) then
        info=stat; return
      end if
      ainv(:,j)=x
    end do
    ainv=0.5_dp*(ainv+transpose(ainv))
  end subroutine invert_matrix

  subroutine cholesky_lower(a, l, info, jitter)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: jitter
    real(dp) :: s, jit
    integer :: n, i, j, k
    n=size(a,1); info=sv_ok; l=0.0_dp; jit=0.0_dp
    if (present(jitter)) jit=jitter
    if (size(a,2)/=n .or. any(shape(l)/=[n,n])) then
      info=sv_singular; return
    end if
    do i=1,n
      do j=1,i
        s=a(i,j)
        if (i==j) s=s+jit
        do k=1,j-1
          s=s-l(i,k)*l(j,k)
        end do
        if (i==j) then
          if (s<=tiny_dp) then
            info=sv_singular; return
          end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine mvn_draws(mean, cov, z, draws, info)
    real(dp), intent(in) :: mean(:), cov(:,:), z(:,:)
    real(dp), intent(out) :: draws(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:)
    integer :: n, ns, j, stat
    n=size(mean); ns=size(z,2); info=sv_ok
    if (any(shape(cov)/=[n,n]) .or. size(z,1)/=n .or. any(shape(draws)/=[n,ns])) then
      info=sv_singular; return
    end if
    allocate(l(n,n))
    call cholesky_lower(cov,l,stat,jitter=1.0e-12_dp)
    if (stat/=sv_ok) then
      info=stat; return
    end if
    do j=1,ns
      draws(:,j)=mean+matmul(l,z(:,j))
    end do
  end subroutine mvn_draws

end module stochvoltmb_linalg
