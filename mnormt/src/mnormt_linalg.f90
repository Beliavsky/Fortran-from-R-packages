module mnormt_linalg
  use mnormt_special, only: dp
  implicit none
  private
  public :: cholesky_upper, pd_solve, solve_linear, covariance_to_correlation
contains
  subroutine cholesky_upper(a, r, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: r(size(a,1),size(a,2))
    integer, intent(out) :: info
    integer :: n,i,j,k
    real(dp) :: s
    n=size(a,1); r=0.0_dp; info=0
    if (size(a,2)/=n) then; info=-1; return; end if
    do j=1,n
      do i=1,j
        s=a(i,j)
        do k=1,i-1
          s=s-r(k,i)*r(k,j)
        end do
        if (i==j) then
          if (s<=0.0_dp) then; info=i; return; end if
          r(i,i)=sqrt(s)
        else
          r(i,j)=s/r(i,i)
        end if
      end do
    end do
  end subroutine cholesky_upper

  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), intent(out) :: x(size(b))
    integer, intent(out) :: info
    real(dp), allocatable :: aug(:,:)
    real(dp) :: tmp,fac
    integer :: n,i,j,k,p
    n=size(b); info=0
    if (size(a,1)/=n .or. size(a,2)/=n) then; info=-1; return; end if
    allocate(aug(n,n+1)); aug(:,1:n)=a; aug(:,n+1)=b
    do k=1,n
      p=k
      do i=k+1,n
        if (abs(aug(i,k))>abs(aug(p,k))) p=i
      end do
      if (abs(aug(p,k))<=tiny(1.0_dp)) then; info=k; x=0.0_dp; return; end if
      if (p/=k) then
        do j=k,n+1
          tmp=aug(k,j); aug(k,j)=aug(p,j); aug(p,j)=tmp
        end do
      end if
      do i=k+1,n
        fac=aug(i,k)/aug(k,k); aug(i,k:n+1)=aug(i,k:n+1)-fac*aug(k,k:n+1)
      end do
    end do
    do i=n,1,-1
      x(i)=(aug(i,n+1)-dot_product(aug(i,i+1:n),x(i+1:n)))/aug(i,i)
    end do
  end subroutine solve_linear

  subroutine pd_solve(a, ainv, logdet, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    real(dp), intent(out), optional :: logdet
    integer, intent(out) :: info
    real(dp), allocatable :: r(:,:), y(:), x(:)
    integer :: n,i,j
    n=size(a,1); allocate(r(n,n),y(n),x(n)); call cholesky_upper(a,r,info)
    if (info/=0) then; ainv=0.0_dp; if(present(logdet)) logdet=huge(1.0_dp); return; end if
    if(present(logdet)) logdet=2.0_dp*sum(log([(r(i,i),i=1,n)]))
    ainv=0.0_dp
    do j=1,n
      y=0.0_dp
      do i=1,n
        y(i)=merge(1.0_dp,0.0_dp,i==j)
        if(i>1) y(i)=y(i)-dot_product(r(1:i-1,i),y(1:i-1))
        y(i)=y(i)/r(i,i)
      end do
      x=0.0_dp
      do i=n,1,-1
        x(i)=y(i)
        if(i<n) x(i)=x(i)-dot_product(r(i,i+1:n),x(i+1:n))
        x(i)=x(i)/r(i,i)
      end do
      ainv(:,j)=x
    end do
    ainv=0.5_dp*(ainv+transpose(ainv))
  end subroutine pd_solve

  subroutine covariance_to_correlation(s,rho,sd,info)
    real(dp), intent(in) :: s(:,:)
    real(dp), intent(out) :: rho(size(s,1),size(s,2)), sd(size(s,1))
    integer, intent(out) :: info
    integer :: i,j,n
    n=size(s,1); info=0
    do i=1,n
      if(s(i,i)<=0.0_dp) then; info=i; return; end if
      sd(i)=sqrt(s(i,i))
    end do
    do j=1,n; do i=1,n; rho(i,j)=s(i,j)/(sd(i)*sd(j)); end do; end do
  end subroutine covariance_to_correlation
end module mnormt_linalg
