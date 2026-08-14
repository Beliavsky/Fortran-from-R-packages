! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! Vendored in KrigInv-fortran under the GPL-3 option; see
! licenses/DiceKriging/LICENSE-GPL-3 and UPSTREAM.md.
module dk_linalg
  use dk_kinds, only : dp, pi_dp
  implicit none
  private
  public :: chol_lower, solve_lower, solve_upper, solve_chol, invert_spd
  public :: least_squares_normal, logdet_from_chol, diag_aba, normal_fill
  public :: symmetrize
contains

  subroutine chol_lower(a, l, info, jitter)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: jitter
    integer :: n, i, j, k
    real(dp) :: s, jit

    n = size(a,1)
    allocate(l(n,n))
    l = 0.0_dp
    info = 0
    jit = 0.0_dp
    if (present(jitter)) jit = jitter
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        if (i == j) s = s + jit
        do k = 1, j-1
          s = s - l(i,k)*l(j,k)
        end do
        if (i == j) then
          if (s <= max(1.0e-15_dp, epsilon(1.0_dp)*max(1.0_dp, abs(a(i,i))))) then
            info = i
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s/l(j,j)
        end if
      end do
    end do
  end subroutine chol_lower

  subroutine solve_lower(l, b, x)
    real(dp), intent(in) :: l(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer :: n, nrhs, i, j, k
    real(dp) :: s
    n = size(l,1); nrhs = size(b,2)
    allocate(x(n,nrhs)); x = 0.0_dp
    do j = 1, nrhs
      do i = 1, n
        s = b(i,j)
        do k = 1, i-1
          s = s - l(i,k)*x(k,j)
        end do
        x(i,j) = s/l(i,i)
      end do
    end do
  end subroutine solve_lower

  subroutine solve_upper(u, b, x)
    real(dp), intent(in) :: u(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer :: n, nrhs, i, j, k
    real(dp) :: s
    n = size(u,1); nrhs = size(b,2)
    allocate(x(n,nrhs)); x = 0.0_dp
    do j = 1, nrhs
      do i = n, 1, -1
        s = b(i,j)
        do k = i+1, n
          s = s - u(i,k)*x(k,j)
        end do
        x(i,j) = s/u(i,i)
      end do
    end do
  end subroutine solve_upper

  subroutine solve_chol(l, b, x)
    real(dp), intent(in) :: l(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    real(dp), allocatable :: y(:,:), u(:,:)
    call solve_lower(l, b, y)
    u = transpose(l)
    call solve_upper(u, y, x)
  end subroutine solve_chol

  subroutine invert_spd(l, ainv)
    real(dp), intent(in) :: l(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    real(dp), allocatable :: eye(:,:)
    integer :: n, i
    n = size(l,1)
    allocate(eye(n,n)); eye = 0.0_dp
    do i = 1, n
      eye(i,i) = 1.0_dp
    end do
    call solve_chol(l, eye, ainv)
    call symmetrize(ainv)
  end subroutine invert_spd

  subroutine least_squares_normal(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: ata(:,:), atb(:,:), l(:,:)
    ata = matmul(transpose(a), a)
    atb = matmul(transpose(a), b)
    call chol_lower(ata, l, info, jitter=1.0e-14_dp*max(1.0_dp,maxval(abs(ata))))
    if (info /= 0) then
      allocate(x(size(a,2),size(b,2))); x = 0.0_dp
      return
    end if
    call solve_chol(l, atb, x)
  end subroutine least_squares_normal

  pure function logdet_from_chol(l) result(v)
    real(dp), intent(in) :: l(:,:)
    real(dp) :: v
    integer :: i
    v = 0.0_dp
    do i = 1, size(l,1)
      v = v + 2.0_dp*log(l(i,i))
    end do
  end function logdet_from_chol

  function diag_aba(a, b) result(v)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable :: v(:)
    integer :: n, i, j, k
    real(dp) :: s
    n = size(a,1)
    allocate(v(n)); v = 0.0_dp
    do i = 1, n
      s = 0.0_dp
      do j = 1, size(b,1)
        do k = 1, size(b,2)
          s = s + a(i,j)*b(j,k)*a(i,k)
        end do
      end do
      v(i) = s
    end do
  end function diag_aba

  subroutine normal_fill(z)
    real(dp), intent(out) :: z(:,:)
    integer :: i, j
    real(dp) :: u1, u2, r, th
    do j = 1, size(z,2)
      i = 1
      do while (i <= size(z,1))
        call random_number(u1); call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        r = sqrt(-2.0_dp*log(u1)); th = 2.0_dp*pi_dp*u2
        z(i,j) = r*cos(th)
        if (i+1 <= size(z,1)) z(i+1,j) = r*sin(th)
        i = i + 2
      end do
    end do
  end subroutine normal_fill

  subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    integer :: i, j
    do j = 1, size(a,2)
      do i = j+1, size(a,1)
        a(i,j) = 0.5_dp*(a(i,j)+a(j,i))
        a(j,i) = a(i,j)
      end do
    end do
  end subroutine symmetrize

end module dk_linalg
