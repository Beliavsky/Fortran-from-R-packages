! SPDX-License-Identifier: GPL-3.0-only
module imputefin_linalg
  use imputefin_kinds, only : dp
  use imputefin_rng, only : rng_state, rng_normal
  implicit none
  private
  public :: solve_linear, inverse_matrix, cholesky_lower, solve_spd, mvn_sample, logdet_spd, symmetrize
contains
  subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    a = 0.5_dp*(a + transpose(a))
  end subroutine symmetrize

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(out) :: x(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: aug(:,:), rowtmp(:)
    real(dp) :: pivot, factor
    integer :: n, nrhs, i, k, p
    n = size(a,1)
    nrhs = size(b,2)
    info = 0
    if (size(a,2) /= n .or. size(b,1) /= n .or. size(x,1) /= n .or. size(x,2) /= nrhs) then
      info = 1
      return
    end if
    allocate(aug(n,n+nrhs), rowtmp(n+nrhs))
    aug(:,1:n) = a
    aug(:,n+1:n+nrhs) = b
    do k=1,n
      p = k
      do i=k+1,n
        if (abs(aug(i,k)) > abs(aug(p,k))) p=i
      end do
      if (abs(aug(p,k)) <= 100.0_dp*epsilon(1.0_dp)) then
        info = k
        return
      end if
      if (p /= k) then
        rowtmp = aug(k,:)
        aug(k,:) = aug(p,:)
        aug(p,:) = rowtmp
      end if
      pivot = aug(k,k)
      aug(k,:) = aug(k,:)/pivot
      do i=1,n
        if (i == k) cycle
        factor = aug(i,k)
        if (abs(factor) > tiny(1.0_dp)) aug(i,:) = aug(i,:) - factor*aug(k,:)
      end do
    end do
    x = aug(:,n+1:n+nrhs)
  end subroutine solve_linear

  subroutine inverse_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: eye(:,:)
    integer :: i, n
    n=size(a,1)
    allocate(eye(n,n)); eye=0.0_dp
    do i=1,n; eye(i,i)=1.0_dp; end do
    call solve_linear(a, eye, ainv, info)
  end subroutine inverse_matrix

  subroutine cholesky_lower(a, l, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out) :: info
    real(dp) :: s
    integer :: n, i, j, k
    n=size(a,1); l=0.0_dp; info=0
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1; s=s-l(i,k)*l(j,k); end do
        if (i==j) then
          if (s <= 0.0_dp) then; info=i; return; end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine solve_spd(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(out) :: x(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:), y(:,:)
    integer :: n, nrhs, i, j, k
    n=size(a,1); nrhs=size(b,2)
    allocate(l(n,n), y(n,nrhs))
    call cholesky_lower(a,l,info)
    if (info/=0) return
    y=0.0_dp
    do j=1,nrhs
      do i=1,n
        y(i,j)=b(i,j)
        do k=1,i-1; y(i,j)=y(i,j)-l(i,k)*y(k,j); end do
        y(i,j)=y(i,j)/l(i,i)
      end do
      do i=n,1,-1
        x(i,j)=y(i,j)
        do k=i+1,n; x(i,j)=x(i,j)-l(k,i)*x(k,j); end do
        x(i,j)=x(i,j)/l(i,i)
      end do
    end do
  end subroutine solve_spd

  function logdet_spd(a, info) result(v)
    real(dp), intent(in) :: a(:,:)
    integer, intent(out) :: info
    real(dp) :: v
    real(dp), allocatable :: l(:,:)
    integer :: i,n
    n=size(a,1); allocate(l(n,n))
    call cholesky_lower(a,l,info)
    if(info/=0) then; v=huge(1.0_dp); return; end if
    v=0.0_dp
    do i=1,n; v=v+2.0_dp*log(l(i,i)); end do
  end function logdet_spd

  subroutine mvn_sample(mean, cov, state, x, info)
    real(dp), intent(in) :: mean(:), cov(:,:)
    type(rng_state), intent(inout) :: state
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:), z(:), a(:,:)
    integer :: n, i
    n=size(mean); allocate(l(n,n),z(n),a(n,n)); a=cov
    call symmetrize(a)
    do i=1,n; a(i,i)=a(i,i)+1.0e-12_dp*max(1.0_dp,abs(a(i,i))); end do
    call cholesky_lower(a,l,info)
    if(info/=0) return
    do i=1,n; z(i)=rng_normal(state); end do
    x=mean+matmul(l,z)
  end subroutine mvn_sample
end module imputefin_linalg
