! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm_utils
  use mvtnorm_kinds, only : dp
  use mvtnorm_linalg, only : cholesky_lower, inverse_spd
  implicit none
  private
  public :: vech, invech, pack_cholesky, unpack_cholesky, in_box, finite_or_zero
  public :: covariance_ok, identity_matrix_local

contains

  function vech(a) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable :: v(:)
    integer :: n,i,j,k
    n=size(a,1)
    allocate(v(n*(n+1)/2))
    k=0
    do i=1,n
      do j=i,n
        k=k+1
        v(k)=a(i,j)
      end do
    end do
  end function vech

  function invech(v) result(a)
    real(dp), intent(in) :: v(:)
    real(dp), allocatable :: a(:,:)
    integer :: n,i,j,k
    n=int((sqrt(1.0_dp+8.0_dp*real(size(v),dp))-1.0_dp)/2.0_dp+0.5_dp)
    allocate(a(n,n))
    a=0.0_dp
    k=0
    do i=1,n
      do j=i,n
        k=k+1
        a(i,j)=v(k)
        a(j,i)=v(k)
      end do
    end do
  end function invech

  function pack_cholesky(sigma) result(theta)
    real(dp), intent(in) :: sigma(:,:)
    real(dp), allocatable :: theta(:), l(:,:)
    logical :: ok
    character(len=256) :: msg
    integer :: n,i,j,k
    n=size(sigma,1)
    allocate(theta(n*(n+1)/2))
    call cholesky_lower(sigma,l,ok,msg)
    if (.not.ok) then
      theta=0.0_dp
      return
    end if
    k=0
    do i=1,n
      do j=1,i
        k=k+1
        if (i==j) then
          theta(k)=log(max(l(i,j),tiny(1.0_dp)))
        else
          theta(k)=l(i,j)
        end if
      end do
    end do
  end function pack_cholesky

  function unpack_cholesky(theta,n) result(sigma)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: n
    real(dp), allocatable :: sigma(:,:), l(:,:)
    integer :: i,j,k
    allocate(l(n,n))
    l=0.0_dp
    k=0
    do i=1,n
      do j=1,i
        k=k+1
        if (i==j) then
          l(i,j)=exp(theta(k))
        else
          l(i,j)=theta(k)
        end if
      end do
    end do
    sigma=matmul(l,transpose(l))
  end function unpack_cholesky

  logical function in_box(x,lower,upper) result(ok)
    real(dp), intent(in) :: x(:),lower(:),upper(:)
    ok=all(x>=lower .and. x<=upper) .and. all(abs(x)<huge(1.0_dp))
  end function in_box

  elemental real(dp) function finite_or_zero(x) result(v)
    real(dp), intent(in) :: x
    if (abs(x)>=huge(1.0_dp)/2.0_dp) then
      v=0.0_dp
    else
      v=x
    end if
  end function finite_or_zero

  logical function covariance_ok(sigma) result(ok)
    real(dp), intent(in) :: sigma(:,:)
    real(dp), allocatable :: l(:,:)
    character(len=256) :: msg
    if (size(sigma,1)/=size(sigma,2)) then
    ok=.false.
    return
    end if
    call cholesky_lower(sigma,l,ok,msg)
  end function covariance_ok

  function identity_matrix_local(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a=0.0_dp
    do i=1,n
    a(i,i)=1.0_dp
    end do
  end function identity_matrix_local

end module tmvtnorm_utils
