module kriginv_linalg
  use kriginv_kinds, only : dp
  implicit none
  private
  public :: cholesky_lower, solve_lower, solve_upper_from_lower, solve_chol, invert_spd, symmetrize
contains
  subroutine cholesky_lower(a,l,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    logical, intent(out) :: ok
    integer :: n,i,j,k
    real(dp) :: s
    n=size(a,1)
    allocate(l(n,n)); l=0.0_dp; ok=.false.
    if(size(a,2)/=n) return
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1
          s=s-l(i,k)*l(j,k)
        end do
        if(i==j) then
          if(s<=0.0_dp) return
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
    ok=.true.
  end subroutine cholesky_lower

  subroutine solve_lower(l,b,x)
    real(dp), intent(in) :: l(:,:),b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer :: n,m,i,j,k
    n=size(l,1); m=size(b,2); allocate(x(n,m)); x=0.0_dp
    do j=1,m
      do i=1,n
        x(i,j)=b(i,j)
        do k=1,i-1
          x(i,j)=x(i,j)-l(i,k)*x(k,j)
        end do
        x(i,j)=x(i,j)/l(i,i)
      end do
    end do
  end subroutine solve_lower

  subroutine solve_upper_from_lower(l,b,x)
    real(dp), intent(in) :: l(:,:),b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer :: n,m,i,j,k
    n=size(l,1); m=size(b,2); allocate(x(n,m)); x=0.0_dp
    do j=1,m
      do i=n,1,-1
        x(i,j)=b(i,j)
        do k=i+1,n
          x(i,j)=x(i,j)-l(k,i)*x(k,j)
        end do
        x(i,j)=x(i,j)/l(i,i)
      end do
    end do
  end subroutine solve_upper_from_lower

  subroutine solve_chol(l,b,x)
    real(dp), intent(in) :: l(:,:),b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    real(dp), allocatable :: y(:,:)
    call solve_lower(l,b,y)
    call solve_upper_from_lower(l,y,x)
  end subroutine solve_chol

  subroutine invert_spd(a,ainv,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: l(:,:),eye(:,:)
    integer :: n,i
    n=size(a,1)
    call cholesky_lower(a,l,ok)
    if(.not.ok) then
      allocate(ainv(n,n)); ainv=0.0_dp; return
    end if
    allocate(eye(n,n)); eye=0.0_dp
    do i=1,n; eye(i,i)=1.0_dp; end do
    call solve_chol(l,eye,ainv)
    call symmetrize(ainv)
  end subroutine invert_spd

  subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    integer :: i,j
    do j=1,size(a,2)
      do i=j+1,size(a,1)
        a(i,j)=0.5_dp*(a(i,j)+a(j,i))
        a(j,i)=a(i,j)
      end do
    end do
  end subroutine symmetrize
end module kriginv_linalg
