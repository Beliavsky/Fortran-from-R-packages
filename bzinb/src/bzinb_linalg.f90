module bzinb_linalg
  use bzinb_kinds, only : dp
  implicit none
  private
  public :: invert_matrix, symmetrize
contains
  subroutine invert_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:)
    real(dp) :: piv, fac, tmp
    integer :: n, i, j, k, ip
    n = size(a,1); ok = .false.
    if (size(a,2) /= n) return
    allocate(aug(n,2*n)); aug = 0.0_dp
    aug(:,1:n) = a
    do i=1,n; aug(i,n+i)=1.0_dp; end do
    do i=1,n
      ip = i
      do k=i+1,n
        if (abs(aug(k,i)) > abs(aug(ip,i))) ip = k
      end do
      if (abs(aug(ip,i)) < 1.0e-12_dp) return
      if (ip /= i) then
        do j=1,2*n
          tmp=aug(i,j); aug(i,j)=aug(ip,j); aug(ip,j)=tmp
        end do
      end if
      piv=aug(i,i); aug(i,:)=aug(i,:)/piv
      do k=1,n
        if (k==i) cycle
        fac=aug(k,i); aug(k,:)=aug(k,:)-fac*aug(i,:)
      end do
    end do
    ainv=aug(:,n+1:2*n); ok=.true.
  end subroutine invert_matrix

  pure subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    integer :: i,j
    do j=1,size(a,2)
      do i=j+1,size(a,1)
        a(i,j)=0.5_dp*(a(i,j)+a(j,i)); a(j,i)=a(i,j)
      end do
    end do
  end subroutine symmetrize
end module bzinb_linalg
