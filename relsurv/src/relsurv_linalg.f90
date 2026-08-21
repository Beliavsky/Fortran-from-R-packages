module relsurv_linalg
  use relsurv_kinds, only : dp
  implicit none
  private
  public :: solve_linear, inverse_matrix, outer_product
contains
  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: pivot, factor, tmp
    integer :: n, i, j, k, p
    n = size(b)
    ok = .false.
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
    allocate(aa(n,n), bb(n)); aa=a; bb=b
    do k=1,n
      p = k
      do i=k+1,n
        if (abs(aa(i,k)) > abs(aa(p,k))) p=i
      end do
      if (abs(aa(p,k)) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa)))) return
      if (p /= k) then
        do j=1,n
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
    x=0.0_dp
    do i=n,1,-1
      x(i)=(bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
    end do
    ok=.true.
  end subroutine solve_linear

  subroutine inverse_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: e(:), x(:)
    integer :: n, j
    n=size(a,1); ok=.false.
    if (size(a,2)/=n .or. size(ainv,1)/=n .or. size(ainv,2)/=n) return
    allocate(e(n),x(n)); ainv=0.0_dp
    do j=1,n
      e=0.0_dp; e(j)=1.0_dp
      call solve_linear(a,e,x,ok)
      if (.not.ok) return
      ainv(:,j)=x
    end do
  end subroutine inverse_matrix

  pure function outer_product(x,y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i=1,size(x)
      a(i,:)=x(i)*y
    end do
  end function outer_product
end module relsurv_linalg
