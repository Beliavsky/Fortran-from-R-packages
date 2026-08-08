! SPDX-License-Identifier: CECILL-2.0
! Derived from the R package neldermead 1.0-13 and its Scilab lineage.
! See LICENSE and UPSTREAM_PROVENANCE.md.

module neldermead_simplex
  use neldermead_kinds, only : dp
  use neldermead_types, only : nm_simplex, nm_objective
  implicit none
  private
  public :: simplex_sort, simplex_center, simplex_xbar, simplex_size
  public :: simplex_delta_f, simplex_fmean, simplex_fvariance, simplex_fstd
  public :: simplex_build_axes, simplex_build_spendley, simplex_build_pfeffer
  public :: simplex_build_randbounds, simplex_build_given, simplex_oriented_restart
  public :: simplex_shrink, simplex_gradient, set_random_seed

contains

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i, huge(1)-1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine set_random_seed

  subroutine simplex_sort(s)
    type(nm_simplex), intent(inout) :: s
    integer :: i, j, k
    real(dp) :: tf
    real(dp), allocatable :: tx(:)
    allocate(tx(size(s%x,1)))
    do i = 1, size(s%f)-1
      k = i
      do j = i+1, size(s%f)
        if (s%f(j) < s%f(k)) k = j
      end do
      if (k /= i) then
        tf = s%f(i); s%f(i) = s%f(k); s%f(k) = tf
        tx = s%x(:,i); s%x(:,i) = s%x(:,k); s%x(:,k) = tx
      end if
    end do
  end subroutine simplex_sort

  function simplex_center(s) result(c)
    type(nm_simplex), intent(in) :: s
    real(dp), allocatable :: c(:)
    integer :: j
    allocate(c(size(s%x,1)))
    c = 0.0_dp
    do j = 1, size(s%x,2)
      c = c + s%x(:,j)
    end do
    c = c / real(size(s%x,2), dp)
  end function simplex_center

  function simplex_xbar(s, exclude) result(c)
    type(nm_simplex), intent(in) :: s
    integer, intent(in), optional :: exclude(:)
    real(dp), allocatable :: c(:)
    integer :: j, k, m
    logical :: skip
    allocate(c(size(s%x,1)))
    c = 0.0_dp; m = 0
    do j = 1, size(s%x,2)
      skip = .false.
      if (present(exclude)) then
        do k = 1, size(exclude)
          if (j == exclude(k)) then
            skip = .true.; exit
          end if
        end do
      end if
      if (.not. skip) then
        c = c + s%x(:,j); m = m + 1
      end if
    end do
    if (m > 0) c = c / real(m,dp)
  end function simplex_xbar

  real(dp) function simplex_size(s) result(v)
    type(nm_simplex), intent(in) :: s
    integer :: j
    v = 0.0_dp
    do j = 2, size(s%x,2)
      v = max(v, sqrt(sum((s%x(:,j)-s%x(:,1))**2)))
    end do
  end function simplex_size

  real(dp) function simplex_delta_f(s) result(v)
    type(nm_simplex), intent(in) :: s
    v = maxval(s%f) - minval(s%f)
  end function simplex_delta_f

  real(dp) function simplex_fmean(s) result(v)
    type(nm_simplex), intent(in) :: s
    v = sum(s%f)/real(size(s%f),dp)
  end function simplex_fmean

  real(dp) function simplex_fvariance(s) result(v)
    type(nm_simplex), intent(in) :: s
    real(dp) :: m
    m = simplex_fmean(s)
    v = sum((s%f-m)**2)/real(size(s%f),dp)
  end function simplex_fvariance

  real(dp) function simplex_fstd(s) result(v)
    type(nm_simplex), intent(in) :: s
    v = sqrt(max(0.0_dp,simplex_fvariance(s)))
  end function simplex_fstd

  subroutine simplex_build_axes(fn, x0, len, s, fevals)
    procedure(nm_objective) :: fn
    real(dp), intent(in) :: x0(:), len
    type(nm_simplex), intent(out) :: s
    integer, intent(inout) :: fevals
    integer :: n, i
    n=size(x0); allocate(s%x(n,n+1),s%f(n+1))
    s%x(:,1)=x0
    do i=1,n
      s%x(:,i+1)=x0; s%x(i,i+1)=s%x(i,i+1)+len
    end do
    do i=1,n+1
      s%f(i)=fn(s%x(:,i)); fevals=fevals+1
    end do
  end subroutine simplex_build_axes

  subroutine simplex_build_spendley(fn, x0, len, s, fevals)
    procedure(nm_objective) :: fn
    real(dp), intent(in) :: x0(:), len
    type(nm_simplex), intent(out) :: s
    integer, intent(inout) :: fevals
    integer :: n, i, j
    real(dp) :: p, q
    n=size(x0); allocate(s%x(n,n+1),s%f(n+1))
    p = len*(sqrt(real(n+1,dp)) - 1.0_dp + real(n,dp))/(real(n,dp)*sqrt(2.0_dp))
    q = len*(sqrt(real(n+1,dp)) - 1.0_dp)/(real(n,dp)*sqrt(2.0_dp))
    s%x(:,1)=x0
    do j=1,n
      s%x(:,j+1)=x0
      do i=1,n
        if (i==j) then
          s%x(i,j+1)=s%x(i,j+1)+p
        else
          s%x(i,j+1)=s%x(i,j+1)+q
        end if
      end do
    end do
    do i=1,n+1
      s%f(i)=fn(s%x(:,i)); fevals=fevals+1
    end do
  end subroutine simplex_build_spendley

  subroutine simplex_build_pfeffer(fn, x0, deltausual, deltazero, s, fevals)
    procedure(nm_objective) :: fn
    real(dp), intent(in) :: x0(:), deltausual, deltazero
    type(nm_simplex), intent(out) :: s
    integer, intent(inout) :: fevals
    integer :: n, i
    n=size(x0); allocate(s%x(n,n+1),s%f(n+1)); s%x(:,1)=x0
    do i=1,n
      s%x(:,i+1)=x0
      if (abs(x0(i)) > tiny(1.0_dp)) then
        s%x(i,i+1)=(1.0_dp+deltausual)*x0(i)
      else
        s%x(i,i+1)=deltazero
      end if
    end do
    do i=1,n+1
      s%f(i)=fn(s%x(:,i)); fevals=fevals+1
    end do
  end subroutine simplex_build_pfeffer

  subroutine simplex_build_randbounds(fn, x0, lower, upper, npoints, s, fevals)
    procedure(nm_objective) :: fn
    real(dp), intent(in) :: x0(:), lower(:), upper(:)
    integer, intent(in) :: npoints
    type(nm_simplex), intent(out) :: s
    integer, intent(inout) :: fevals
    integer :: n, i
    real(dp), allocatable :: u(:)
    n=size(x0); allocate(s%x(n,npoints),s%f(npoints),u(n))
    s%x(:,1)=x0
    do i=2,npoints
      call random_number(u)
      s%x(:,i)=lower+u*(upper-lower)
    end do
    do i=1,npoints
      s%f(i)=fn(s%x(:,i)); fevals=fevals+1
    end do
  end subroutine simplex_build_randbounds

  subroutine simplex_build_given(fn, coords, s, fevals)
    procedure(nm_objective) :: fn
    real(dp), intent(in) :: coords(:,:)
    type(nm_simplex), intent(out) :: s
    integer, intent(inout) :: fevals
    integer :: i
    allocate(s%x(size(coords,1),size(coords,2)),s%f(size(coords,2)))
    s%x=coords
    do i=1,size(coords,2)
      s%f(i)=fn(s%x(:,i)); fevals=fevals+1
    end do
  end subroutine simplex_build_given

  subroutine simplex_oriented_restart(fn, old, s, fevals)
    procedure(nm_objective) :: fn
    type(nm_simplex), intent(in) :: old
    type(nm_simplex), intent(out) :: s
    integer, intent(inout) :: fevals
    integer :: n, i
    real(dp), allocatable :: best(:)
    n=size(old%x,1); allocate(s%x(n,n+1),s%f(n+1),best(n)); best=old%x(:,1)
    s%x(:,1)=best
    do i=1,n
      s%x(:,i+1)=best + 0.5_dp*(old%x(:,min(i+1,size(old%x,2)))-best)
    end do
    do i=1,n+1
      s%f(i)=fn(s%x(:,i)); fevals=fevals+1
    end do
  end subroutine simplex_oriented_restart

  subroutine simplex_shrink(fn, s, sigma, fevals)
    procedure(nm_objective) :: fn
    type(nm_simplex), intent(inout) :: s
    real(dp), intent(in) :: sigma
    integer, intent(inout) :: fevals
    integer :: j
    do j=2,size(s%x,2)
      s%x(:,j)=s%x(:,1)+sigma*(s%x(:,j)-s%x(:,1))
      s%f(j)=fn(s%x(:,j)); fevals=fevals+1
    end do
  end subroutine simplex_shrink

  subroutine simplex_gradient(s, g, ok)
    type(nm_simplex), intent(in) :: s
    real(dp), intent(out) :: g(:)
    logical, intent(out) :: ok
    integer :: n, m, i, j, k, pivot
    real(dp), allocatable :: a(:,:), b(:), row(:)
    real(dp) :: fac, mx, tmp
    n=size(s%x,1); m=size(s%x,2)-1
    g=0.0_dp; ok=.false.
    if (m < n) return
    allocate(a(n,n),b(n),row(n))
    do i=1,n
      a(i,:)=s%x(:,i+1)-s%x(:,1)
      b(i)=s%f(i+1)-s%f(1)
    end do
    do k=1,n
      pivot=k; mx=abs(a(k,k))
      do i=k+1,n
        if (abs(a(i,k))>mx) then; pivot=i; mx=abs(a(i,k)); end if
      end do
      if (mx <= 100.0_dp*epsilon(1.0_dp)) return
      if (pivot/=k) then
        row=a(k,:); a(k,:)=a(pivot,:); a(pivot,:)=row
        tmp=b(k); b(k)=b(pivot); b(pivot)=tmp
      end if
      do i=k+1,n
        fac=a(i,k)/a(k,k)
        a(i,k:n)=a(i,k:n)-fac*a(k,k:n); b(i)=b(i)-fac*b(k)
      end do
    end do
    do i=n,1,-1
      tmp=b(i)
      do j=i+1,n; tmp=tmp-a(i,j)*g(j); end do
      g(i)=tmp/a(i,i)
    end do
    ok=.true.
  end subroutine simplex_gradient

end module neldermead_simplex
