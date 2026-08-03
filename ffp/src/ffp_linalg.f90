module ffp_linalg
  use ffp_kinds, only : dp
  implicit none
  private
  public :: solve_linear, covariance_to_correlation
contains
  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: m(:, :), rhs(:), row(:)
    real(dp) :: factor, pivot
    integer :: n, i, k, ip
    n = size(b); info = 0
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
      info = -1; return
    end if
    allocate(m(n,n), rhs(n), row(n)); m=a; rhs=b
    do k=1,n-1
      ip=k
      do i=k+1,n
        if (abs(m(i,k)) > abs(m(ip,k))) ip=i
      end do
      if (abs(m(ip,k)) <= 100.0_dp*epsilon(1.0_dp)) then
        info=k; return
      end if
      if (ip /= k) then
        row=m(k,:); m(k,:)=m(ip,:); m(ip,:)=row
        pivot=rhs(k); rhs(k)=rhs(ip); rhs(ip)=pivot
      end if
      do i=k+1,n
        factor=m(i,k)/m(k,k)
        m(i,k:n)=m(i,k:n)-factor*m(k,k:n)
        rhs(i)=rhs(i)-factor*rhs(k)
      end do
    end do
    if (abs(m(n,n)) <= 100.0_dp*epsilon(1.0_dp)) then
      info=n; return
    end if
    x(n)=rhs(n)/m(n,n)
    do i=n-1,1,-1
      x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
    end do
  end subroutine solve_linear

  subroutine covariance_to_correlation(sigma, cor, sd)
    real(dp), intent(in) :: sigma(:, :)
    real(dp), intent(out) :: cor(:, :), sd(:)
    integer :: i,j,n
    n=size(sigma,1)
    do i=1,n
      sd(i)=sqrt(max(sigma(i,i),0.0_dp))
    end do
    do j=1,n; do i=1,n
      if (sd(i)>0.0_dp .and. sd(j)>0.0_dp) then
        cor(i,j)=sigma(i,j)/(sd(i)*sd(j))
      else
        cor(i,j)=0.0_dp
      end if
    end do; end do
  end subroutine covariance_to_correlation
end module ffp_linalg
