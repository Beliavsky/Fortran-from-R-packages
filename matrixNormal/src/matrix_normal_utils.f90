! SPDX-License-Identifier: GPL-3.0-only
module matrix_normal_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use mvtnorm_kinds, only : dp
  use mvtnorm_linalg, only : jacobi_eigen
  implicit none
  private
  public :: identity_matrix, ones_matrix, tr, vec, vech, kronecker_product
  public :: is_square_matrix, is_symmetric_matrix
  public :: is_positive_definite, is_positive_semidefinite

contains

  function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp), allocatable :: a(:,:)
    integer :: k
    allocate(a(n,n))
    a = 0.0_dp
    do k = 1, n
      a(k,k) = 1.0_dp
    end do
  end function identity_matrix

  function ones_matrix(n,m) result(a)
    integer, intent(in) :: n
    integer, intent(in), optional :: m
    real(dp), allocatable :: a(:,:)
    integer :: mm
    mm = n
    if (present(m)) mm = m
    allocate(a(n,mm))
    a = 1.0_dp
  end function ones_matrix

  pure logical function is_square_matrix(a) result(ok)
    real(dp), intent(in) :: a(:,:)
    ok = size(a,1) == size(a,2)
  end function is_square_matrix

  logical function is_symmetric_matrix(a,tol) result(ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp) :: eps
    if (.not. is_square_matrix(a)) then
      ok = .false.
      return
    end if
    if (any(ieee_is_nan(a))) then
      ok = .false.
      return
    end if
    eps = sqrt(epsilon(1.0_dp))
    if (present(tol)) eps = tol
    ok = sum(abs(a-transpose(a))) < eps
  end function is_symmetric_matrix

  logical function is_positive_semidefinite(a,tol) result(ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: values(:), vectors(:,:)
    real(dp) :: eps
    logical :: eig_ok
    integer :: k
    eps = sqrt(epsilon(1.0_dp))
    if (present(tol)) eps = tol
    if (.not. is_symmetric_matrix(a,eps)) then
      ok = .false.
      return
    end if
    call jacobi_eigen(a,values,vectors,eig_ok)
    if (.not. eig_ok) then
      ok = .false.
      return
    end if
    do k=1,size(values)
      if (abs(values(k)) < eps) values(k)=0.0_dp
    end do
    ok = minval(values) >= 0.0_dp
  end function is_positive_semidefinite

  logical function is_positive_definite(a,tol) result(ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: values(:), vectors(:,:)
    real(dp) :: eps
    logical :: eig_ok
    integer :: k
    eps = sqrt(epsilon(1.0_dp))
    if (present(tol)) eps = tol
    if (.not. is_symmetric_matrix(a,eps)) then
      ok = .false.
      return
    end if
    call jacobi_eigen(a,values,vectors,eig_ok)
    if (.not. eig_ok) then
      ok = .false.
      return
    end if
    do k=1,size(values)
      if (abs(values(k)) < eps) values(k)=0.0_dp
    end do
    ok = minval(values) > 0.0_dp
  end function is_positive_definite

  pure real(dp) function tr(a) result(v)
    real(dp), intent(in) :: a(:,:)
    integer :: k,n
    v = 0.0_dp
    n = min(size(a,1),size(a,2))
    do k=1,n
      v = v + a(k,k)
    end do
  end function tr

  pure function vec(a) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable :: v(:)
    allocate(v(size(a)))
    v = reshape(a,[size(a)])
  end function vec

  function vech(a,tol,ok) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    logical, intent(out), optional :: ok
    real(dp), allocatable :: v(:)
    integer :: n,i,j,k
    logical :: sym
    sym = is_symmetric_matrix(a,tol)
    if (.not.sym) then
      allocate(v(0))
      if (present(ok)) ok=.false.
      return
    end if
    n=size(a,1)
    allocate(v(n*(n+1)/2))
    k=0
    do j=1,n
      do i=j,n
        k=k+1
        v(k)=a(i,j)
      end do
    end do
    if (present(ok)) ok=.true.
  end function vech

  pure function kronecker_product(a,b) result(c)
    real(dp), intent(in) :: a(:,:),b(:,:)
    real(dp), allocatable :: c(:,:)
    integer :: i,j,mb,nb
    mb=size(b,1)
    nb=size(b,2)
    allocate(c(size(a,1)*mb,size(a,2)*nb))
    do j=1,size(a,2)
      do i=1,size(a,1)
        c((i-1)*mb+1:i*mb,(j-1)*nb+1:j*nb)=a(i,j)*b
      end do
    end do
  end function kronecker_product

end module matrix_normal_utils
