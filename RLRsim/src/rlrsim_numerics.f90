module rlrsim_numerics
  use, intrinsic :: iso_fortran_env, only : int64
  use rlrsim_kinds, only : dp
  implicit none
  private

  public :: seed_random_number
  public :: chisq_random
  public :: singular_values_squared
  public :: residualize_against
  public :: sort_descending

contains

  subroutine seed_random_number(seed)
    integer, intent(in) :: seed !! User seed used to initialize Fortran's intrinsic pseudorandom-number generator.
    integer :: n
    integer :: i
    integer(int64) :: state
    integer(int64) :: modulus
    integer, allocatable :: put(:)

    call random_seed(size = n)
    allocate(put(n))
    modulus = int(huge(1), int64)
    state = modulo(int(seed, int64), modulus)
    do i = 1, n
      state = modulo(1103515245_int64 * state + 12345_int64 + 104729_int64 * int(i, int64), modulus)
      put(i) = int(max(1_int64, state))
    end do
    call random_seed(put = put)
  end subroutine seed_random_number

  real(dp) function normal_random() result(z)
    real(dp) :: u1
    real(dp) :: u2
    real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)

    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp * log(u1)) * cos(twopi * u2)
  end function normal_random

  recursive real(dp) function gamma_random(shape) result(x)
    real(dp), intent(in) :: shape !! Gamma shape parameter; must be strictly positive, with unit scale.
    real(dp) :: c
    real(dp) :: d
    real(dp) :: u
    real(dp) :: v
    real(dp) :: z

    if (shape <= 0.0_dp) error stop "gamma_random: shape must be positive"
    if (shape < 1.0_dp) then
      call random_number(u)
      u = max(u, tiny(1.0_dp))
      x = gamma_random(shape + 1.0_dp) * u ** (1.0_dp / shape)
      return
    end if

    d = shape - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      z = normal_random()
      v = 1.0_dp + c * z
      if (v <= 0.0_dp) cycle
      v = v ** 3
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp * z ** 4) exit
      if (log(max(u, tiny(1.0_dp))) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
    end do
    x = d * v
  end function gamma_random

  real(dp) function chisq_random(df) result(x)
    real(dp), intent(in) :: df !! Chi-square degrees of freedom; zero returns exactly zero and positive values are sampled.

    if (df < 0.0_dp) error stop "chisq_random: df must be nonnegative"
    if (df == 0.0_dp) then
      x = 0.0_dp
    else
      x = 2.0_dp * gamma_random(0.5_dp * df)
    end if
  end function chisq_random

  subroutine residualize_against(x, z, zr)
    real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix, shape (n,p), whose column space is projected out.
    real(dp), intent(in) :: z(:, :) !! Random-effect design matrix, shape (n,k), to be residualized against x.
    real(dp), allocatable, intent(out) :: zr(:, :) !! Residualized z, shape (n,k), for full-rank x.
    real(dp), allocatable :: q(:, :)
    real(dp), allocatable :: v(:)
    real(dp) :: nrm
    real(dp) :: tol
    integer :: i
    integer :: j
    integer :: rank
    integer :: n
    integer :: p

    n = size(x, 1)
    p = size(x, 2)
    if (size(z, 1) /= n) error stop "residualize_against: incompatible row counts"
    allocate(q(n, p), source = 0.0_dp)
    allocate(v(n))
    rank = 0
    tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(x))) * sqrt(real(max(1, n), dp))

    do j = 1, p
      v = x(:, j)
      do i = 1, rank
        v = v - dot_product(q(:, i), v) * q(:, i)
      end do
      do i = 1, rank
        v = v - dot_product(q(:, i), v) * q(:, i)
      end do
      nrm = sqrt(max(0.0_dp, dot_product(v, v)))
      if (nrm > tol) then
        rank = rank + 1
        q(:, rank) = v / nrm
      end if
    end do

    allocate(zr(size(z, 1), size(z, 2)), source = z)
    do j = 1, size(zr, 2)
      do i = 1, rank
        zr(:, j) = zr(:, j) - dot_product(q(:, i), zr(:, j)) * q(:, i)
      end do
    end do
  end subroutine residualize_against

  subroutine singular_values_squared(a, values)
    real(dp), intent(in) :: a(:, :) !! Real matrix whose squared singular values are requested.
    real(dp), allocatable, intent(out) :: values(:) !! Squared singular values in descending order, length min(size(a,1),size(a,2)).
    real(dp), allocatable :: gram(:, :)
    real(dp), allocatable :: eig(:)
    integer :: m
    integer :: n
    integer :: r

    m = size(a, 1)
    n = size(a, 2)
    r = min(m, n)
    if (m <= n) then
      gram = matmul(a, transpose(a))
    else
      gram = matmul(transpose(a), a)
    end if
    call jacobi_eigenvalues(gram, eig)
    eig = max(eig, 0.0_dp)
    call sort_descending(eig)
    allocate(values(r))
    values = eig(1:r)
  end subroutine singular_values_squared

  subroutine jacobi_eigenvalues(a_in, eig)
    real(dp), intent(in) :: a_in(:, :) !! Real symmetric square matrix whose eigenvalues are requested.
    real(dp), allocatable, intent(out) :: eig(:) !! Eigenvalues of a_in in unspecified order.
    real(dp), allocatable :: a(:, :)
    real(dp) :: app
    real(dp) :: aqq
    real(dp) :: apq
    real(dp) :: c
    real(dp) :: s
    real(dp) :: tau
    real(dp) :: t
    real(dp) :: aik
    real(dp) :: akq
    real(dp) :: off
    real(dp) :: scale
    integer :: i
    integer :: k
    integer :: p
    integer :: q
    integer :: sweep
    integer :: n

    n = size(a_in, 1)
    if (size(a_in, 2) /= n) error stop "jacobi_eigenvalues: matrix must be square"
    allocate(a(n, n), source = a_in)
    if (n == 0) then
      allocate(eig(0))
      return
    end if
    scale = max(1.0_dp, maxval(abs(a)))
    do sweep = 1, max(20, 12 * n * n)
      off = 0.0_dp
      p = 1
      q = min(2, n)
      do i = 1, n - 1
        do k = i + 1, n
          if (abs(a(i, k)) > off) then
            off = abs(a(i, k))
            p = i
            q = k
          end if
        end do
      end do
      if (n == 1 .or. off <= 100.0_dp * epsilon(1.0_dp) * scale) exit
      apq = a(p, q)
      app = a(p, p)
      aqq = a(q, q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t * t)
      s = t * c
      do k = 1, n
        if (k == p .or. k == q) cycle
        aik = a(k, p)
        akq = a(k, q)
        a(k, p) = c * aik - s * akq
        a(p, k) = a(k, p)
        a(k, q) = s * aik + c * akq
        a(q, k) = a(k, q)
      end do
      a(p, p) = c * c * app - 2.0_dp * s * c * apq + s * s * aqq
      a(q, q) = s * s * app + 2.0_dp * s * c * apq + c * c * aqq
      a(p, q) = 0.0_dp
      a(q, p) = 0.0_dp
    end do
    allocate(eig(n))
    do i = 1, n
      eig(i) = a(i, i)
    end do
  end subroutine jacobi_eigenvalues

  pure subroutine sort_descending(x)
    real(dp), intent(inout) :: x(:) !! Real vector sorted in place from largest to smallest.
    real(dp) :: key
    integer :: i
    integer :: j

    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) >= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_descending

end module rlrsim_numerics
