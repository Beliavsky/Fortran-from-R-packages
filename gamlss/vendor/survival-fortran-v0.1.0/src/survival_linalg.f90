! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_linalg
  use survival_kinds, only : dp
  implicit none
  private
  public :: solve_sym, invert_matrix, matrix_rank
contains
  subroutine solve_sym(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: m(:,:), rhs(:), row(:)
    real(dp) :: piv, fac
    integer :: n, i, k, p
    n=size(b); allocate(m(n,n),rhs(n),row(n)); m=a; rhs=b; ok=.true.
    do k=1,n
       p=k
       do i=k+1,n
          if (abs(m(i,k)) > abs(m(p,k))) p=i
       end do
       if (abs(m(p,k)) <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(m)))) then
          ok=.false.; x=0.0_dp; return
       end if
       if (p/=k) then
          row=m(k,:); m(k,:)=m(p,:); m(p,:)=row
          piv=rhs(k); rhs(k)=rhs(p); rhs(p)=piv
       end if
       piv=m(k,k)
       do i=k+1,n
          fac=m(i,k)/piv
          m(i,k)=0.0_dp
          m(i,k+1:n)=m(i,k+1:n)-fac*m(k,k+1:n)
          rhs(i)=rhs(i)-fac*rhs(k)
       end do
    end do
    x=rhs
    do i=n,1,-1
       if (i<n) x(i)=x(i)-dot_product(m(i,i+1:n),x(i+1:n))
       x(i)=x(i)/m(i,i)
    end do
  end subroutine

  subroutine invert_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: e(:), x(:)
    integer :: n,j
    n=size(a,1); allocate(e(n),x(n)); ainv=0.0_dp; ok=.true.
    do j=1,n
       e=0.0_dp; e(j)=1.0_dp
       call solve_sym(a,e,x,ok)
       if (.not.ok) return
       ainv(:,j)=x
    end do
  end subroutine

  integer function matrix_rank(a, tol) result(r)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: m(:,:), row(:)
    real(dp) :: t, fac
    integer :: nr,nc,i,j,k,p
    nr=size(a,1); nc=size(a,2); allocate(m(nr,nc),row(nc)); m=a
    t=sqrt(epsilon(1.0_dp))*max(1.0_dp,maxval(abs(a))); if(present(tol)) t=tol
    r=0; i=1
    do j=1,nc
       if(i>nr) exit
       p=i
       do k=i+1,nr
          if(abs(m(k,j))>abs(m(p,j))) p=k
       end do
       if(abs(m(p,j))<=t) cycle
       if(p/=i) then; row=m(i,:);m(i,:)=m(p,:);m(p,:)=row;end if
       do k=i+1,nr
          fac=m(k,j)/m(i,j); m(k,j:nc)=m(k,j:nc)-fac*m(i,j:nc)
       end do
       r=r+1; i=i+1
    end do
  end function
end module survival_linalg
