module vamc_math
  use vamc_kinds, only: dp, i8
  use vamc_status, only: status_type, vamc_success, vamc_dimension_error, vamc_numerical_error
  implicit none
  private
  public :: cholesky_upper, solve_linear, seed_rng, uniform_random, normal_random
  public :: sample_integer, random_permutation, is_positive_semidefinite

  integer(i8), save :: rng_state(4) = [123456789_i8, 362436069_i8, 521288629_i8, 88675123_i8]
contains
  subroutine cholesky_upper(a, r, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: r(:,:)
    type(status_type), intent(inout), optional :: status
    integer :: n, i, j, k
    real(dp) :: s
    if (present(status)) call status%clear()
    n = size(a, 1)
    if (size(a, 2) /= n) then
      allocate(r(0,0))
      if (present(status)) call status%set(vamc_dimension_error, 'Cholesky input must be square.')
      return
    end if
    allocate(r(n,n), source=0.0_dp)
    do j = 1, n
      do i = 1, j
        s = a(i,j)
        do k = 1, i - 1
          s = s - r(k,i) * r(k,j)
        end do
        if (i == j) then
          if (s <= epsilon(1.0_dp) * max(1.0_dp, abs(a(i,i)))) then
            if (present(status)) call status%set(vamc_numerical_error, 'Matrix is not positive definite.')
            return
          end if
          r(i,j) = sqrt(s)
        else
          r(i,j) = s / r(i,i)
        end if
      end do
    end do
  end subroutine cholesky_upper

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    type(status_type), intent(inout), optional :: status
    real(dp), allocatable :: aug(:,:)
    real(dp) :: pivot, factor
    integer :: n, i, j, k, p
    if (present(status)) call status%clear()
    n = size(b)
    if (size(a,1) /= n .or. size(a,2) /= n) then
      allocate(x(0))
      if (present(status)) call status%set(vamc_dimension_error, 'Linear system dimensions do not conform.')
      return
    end if
    allocate(aug(n,n+1))
    aug(:,1:n) = a
    aug(:,n+1) = b
    do k = 1, n
      p = k
      do i = k + 1, n
        if (abs(aug(i,k)) > abs(aug(p,k))) p = i
      end do
      if (abs(aug(p,k)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aug(:,k))))) then
        allocate(x(0))
        if (present(status)) call status%set(vamc_numerical_error, 'Singular linear system.')
        return
      end if
      if (p /= k) then
        do j = k, n + 1
          pivot = aug(k,j)
          aug(k,j) = aug(p,j)
          aug(p,j) = pivot
        end do
      end if
      pivot = aug(k,k)
      aug(k,k:n+1) = aug(k,k:n+1) / pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        aug(i,k:n+1) = aug(i,k:n+1) - factor * aug(k,k:n+1)
      end do
    end do
    allocate(x(n))
    x = aug(:,n+1)
  end subroutine solve_linear

  pure integer(i8) function rotl(x, k)
    integer(i8), intent(in) :: x
    integer, intent(in) :: k
    rotl = ior(shiftl(x, k), shiftr(x, 64-k))
  end function rotl

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer(i8) :: z
    integer :: i
    z = int(seed, i8)
    if (z == 0_i8) z = 1_i8
    do i = 1, 4
      z = z + int(z'9E3779B97F4A7C15', i8)
      z = ieor(z, shiftr(z, 30)) * int(z'BF58476D1CE4E5B9', i8)
      z = ieor(z, shiftr(z, 27)) * int(z'94D049BB133111EB', i8)
      rng_state(i) = ieor(z, shiftr(z, 31))
    end do
  end subroutine seed_rng

  integer(i8) function next_u64()
    integer(i8) :: result, t
    result = rotl(rng_state(2) * 5_i8, 7) * 9_i8
    t = shiftl(rng_state(2), 17)
    rng_state(3) = ieor(rng_state(3), rng_state(1))
    rng_state(4) = ieor(rng_state(4), rng_state(2))
    rng_state(2) = ieor(rng_state(2), rng_state(3))
    rng_state(1) = ieor(rng_state(1), rng_state(4))
    rng_state(3) = ieor(rng_state(3), t)
    rng_state(4) = rotl(rng_state(4), 45)
    next_u64 = result
  end function next_u64

  real(dp) function uniform_random()
    integer(i8) :: x
    x = shiftr(next_u64(), 11)
    uniform_random = (real(iand(x, int(z'001FFFFFFFFFFFFF', i8)), dp) + 0.5_dp) / 9007199254740992.0_dp
  end function uniform_random

  real(dp) function normal_random()
    real(dp) :: u1, u2
    real(dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp
    u1 = max(uniform_random(), tiny(1.0_dp))
    u2 = uniform_random()
    normal_random = sqrt(-2.0_dp * log(u1)) * cos(twopi * u2)
  end function normal_random

  integer function sample_integer(low, high)
    integer, intent(in) :: low, high
    if (high < low) then
      sample_integer = low
    else
      sample_integer = low + min(high-low, int(uniform_random() * real(high-low+1,dp)))
    end if
  end function sample_integer

  subroutine random_permutation(values)
    integer, intent(inout) :: values(:)
    integer :: i, j, temp
    do i = size(values), 2, -1
      j = sample_integer(1, i)
      temp = values(i)
      values(i) = values(j)
      values(j) = temp
    end do
  end subroutine random_permutation

  logical function is_positive_semidefinite(a, tolerance)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: b(:,:), r(:,:)
    real(dp) :: tol
    type(status_type) :: status
    integer :: n, i
    tol = 1.0e-12_dp
    if (present(tolerance)) tol = tolerance
    n = size(a,1)
    if (size(a,2) /= n) then
      is_positive_semidefinite = .false.
      return
    end if
    allocate(b(n,n))
    b = 0.5_dp * (a + transpose(a))
    do i = 1, n
      b(i,i) = b(i,i) + tol
    end do
    call cholesky_upper(b, r, status)
    is_positive_semidefinite = status%code == vamc_success
  end function is_positive_semidefinite
end module vamc_math
