! SPDX-License-Identifier: GPL-2.0-or-later
module mcga_bytes
  use, intrinsic :: iso_c_binding, only : c_int, c_long
  use mcga_kinds, only : dp, i8, i32
  use mcga_random, only : runif_scalar, randint
  implicit none
  private

  public :: max_double, size_of_double, size_of_int, size_of_long
  public :: double_to_bytes, double_vector_to_bytes, bytes_to_double, byte_vector_to_doubles
  public :: one_point_crossover, one_point_crossover_doubles
  public :: two_point_crossover, two_point_crossover_doubles
  public :: uniform_crossover, uniform_crossover_doubles
  public :: byte_code_mutation, byte_code_mutation_doubles, byte_code_mutation_doubles_random
  public :: ensure_bounds

contains
  pure function max_double() result(x)
    real(dp) :: x
    x = huge(0.0_dp)
  end function max_double

  pure function size_of_double() result(n)
    integer :: n
    n = storage_size(0.0_dp) / 8
  end function size_of_double

  pure function size_of_int() result(n)
    integer :: n
    integer(c_int) :: x
    n = storage_size(x) / 8
  end function size_of_int

  pure function size_of_long() result(n)
    integer :: n
    integer(c_long) :: x
    n = storage_size(x) / 8
  end function size_of_long

  pure function double_to_bytes(x) result(bytes)
    real(dp), intent(in) :: x
    integer(i32) :: bytes(storage_size(x) / 8)
    integer(i8) :: raw(storage_size(x) / 8)
    integer :: i

    raw = transfer(x, raw)
    do i = 1, size(raw)
      bytes(i) = iand(int(raw(i), i32), int(z'FF', i32))
    end do
  end function double_to_bytes

  pure function bytes_to_double(bytes) result(x)
    integer(i32), intent(in) :: bytes(:)
    real(dp) :: x
    integer(i8) :: raw(storage_size(x) / 8)
    integer :: i, nb

    nb = size(raw)
    if (size(bytes) /= nb) error stop "bytes_to_double: wrong byte count"
    do i = 1, nb
      raw(i) = int(merge(bytes(i), bytes(i) - 256_i32, bytes(i) <= 127_i32), i8)
    end do
    x = transfer(raw, x)
  end function bytes_to_double

  pure function double_vector_to_bytes(x) result(bytes)
    real(dp), intent(in) :: x(:)
    integer(i32), allocatable :: bytes(:)
    integer(i32), allocatable :: b(:)
    integer :: i, nb, p

    nb = size_of_double()
    allocate(bytes(nb * size(x)))
    p = 1
    do i = 1, size(x)
      b = double_to_bytes(x(i))
      bytes(p:p + nb - 1) = b
      p = p + nb
    end do
  end function double_vector_to_bytes

  pure function byte_vector_to_doubles(bytes) result(x)
    integer(i32), intent(in) :: bytes(:)
    real(dp), allocatable :: x(:)
    integer :: nb, n, i, p

    nb = size_of_double()
    if (mod(size(bytes), nb) /= 0) error stop "byte_vector_to_doubles: byte count not divisible by sizeof(double)"
    n = size(bytes) / nb
    allocate(x(n))
    p = 1
    do i = 1, n
      x(i) = bytes_to_double(bytes(p:p + nb - 1))
      p = p + nb
    end do
  end function byte_vector_to_doubles

  subroutine one_point_crossover(bytes1, bytes2, cutpoint, child1, child2)
    integer(i32), intent(in) :: bytes1(:), bytes2(:)
    integer, intent(in) :: cutpoint
    integer(i32), allocatable, intent(out) :: child1(:), child2(:)
    integer :: cp, n

    n = size(bytes1)
    if (size(bytes2) /= n) error stop "one_point_crossover: shape mismatch"
    cp = max(0, min(cutpoint, n))
    allocate(child1(n), child2(n))
    if (cp > 0) then
      child1(:cp) = bytes1(:cp)
      child2(:cp) = bytes2(:cp)
    end if
    if (cp < n) then
      child1(cp + 1:) = bytes2(cp + 1:)
      child2(cp + 1:) = bytes1(cp + 1:)
    end if
  end subroutine one_point_crossover

  subroutine one_point_crossover_doubles(x1, x2, cutpoint, child1, child2)
    real(dp), intent(in) :: x1(:), x2(:)
    integer, intent(in) :: cutpoint
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    integer(i32), allocatable :: b1(:), b2(:), c1(:), c2(:)

    b1 = double_vector_to_bytes(x1)
    b2 = double_vector_to_bytes(x2)
    call one_point_crossover(b1, b2, cutpoint, c1, c2)
    child1 = byte_vector_to_doubles(c1)
    child2 = byte_vector_to_doubles(c2)
  end subroutine one_point_crossover_doubles

  subroutine two_point_crossover(bytes1, bytes2, cutpoint1, cutpoint2, child1, child2)
    integer(i32), intent(in) :: bytes1(:), bytes2(:)
    integer, intent(in) :: cutpoint1, cutpoint2
    integer(i32), allocatable, intent(out) :: child1(:), child2(:)
    integer :: i, lo, hi, n

    n = size(bytes1)
    if (size(bytes2) /= n) error stop "two_point_crossover: shape mismatch"
    lo = max(0, min(min(cutpoint1, cutpoint2), n))
    hi = max(0, min(max(cutpoint1, cutpoint2), n))
    allocate(child1(n), child2(n))
    do i = 1, n
      ! Rcpp uses zero-based i and swaps for cutpoint1 <= i <= cutpoint2.
      if ((i - 1 < lo) .or. (i - 1 > hi)) then
        child1(i) = bytes1(i)
        child2(i) = bytes2(i)
      else
        child1(i) = bytes2(i)
        child2(i) = bytes1(i)
      end if
    end do
  end subroutine two_point_crossover

  subroutine two_point_crossover_doubles(x1, x2, cutpoint1, cutpoint2, child1, child2)
    real(dp), intent(in) :: x1(:), x2(:)
    integer, intent(in) :: cutpoint1, cutpoint2
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    integer(i32), allocatable :: b1(:), b2(:), c1(:), c2(:)

    b1 = double_vector_to_bytes(x1)
    b2 = double_vector_to_bytes(x2)
    call two_point_crossover(b1, b2, cutpoint1, cutpoint2, c1, c2)
    child1 = byte_vector_to_doubles(c1)
    child2 = byte_vector_to_doubles(c2)
  end subroutine two_point_crossover_doubles

  subroutine uniform_crossover(bytes1, bytes2, child1, child2)
    integer(i32), intent(in) :: bytes1(:), bytes2(:)
    integer(i32), allocatable, intent(out) :: child1(:), child2(:)
    integer :: i, n

    n = size(bytes1)
    if (size(bytes2) /= n) error stop "uniform_crossover: shape mismatch"
    allocate(child1(n), child2(n))
    do i = 1, n
      if (runif_scalar(0.0_dp, 1.0_dp) < 0.5_dp) then
        child1(i) = bytes1(i)
        child2(i) = bytes2(i)
      else
        child1(i) = bytes2(i)
        child2(i) = bytes1(i)
      end if
    end do
  end subroutine uniform_crossover

  subroutine uniform_crossover_doubles(x1, x2, child1, child2)
    real(dp), intent(in) :: x1(:), x2(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    integer(i32), allocatable :: b1(:), b2(:), c1(:), c2(:)

    b1 = double_vector_to_bytes(x1)
    b2 = double_vector_to_bytes(x2)
    call uniform_crossover(b1, b2, c1, c2)
    child1 = byte_vector_to_doubles(c1)
    child2 = byte_vector_to_doubles(c2)
  end subroutine uniform_crossover_doubles

  subroutine byte_code_mutation(bytes, pmutation)
    integer(i32), intent(inout) :: bytes(:)
    real(dp), intent(in) :: pmutation
    integer :: i

    do i = 1, size(bytes)
      if (runif_scalar(0.0_dp, 1.0_dp) < pmutation) then
        if (runif_scalar(0.0_dp, 1.0_dp) < 0.5_dp) then
          bytes(i) = modulo(bytes(i) + 1_i32, 256_i32)
        else
          bytes(i) = modulo(bytes(i) - 1_i32, 256_i32)
        end if
      end if
    end do
  end subroutine byte_code_mutation

  subroutine byte_code_mutation_doubles(x, pmutation)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: pmutation
    integer(i32), allocatable :: bytes(:)

    bytes = double_vector_to_bytes(x)
    call byte_code_mutation(bytes, pmutation)
    x = byte_vector_to_doubles(bytes)
  end subroutine byte_code_mutation_doubles

  subroutine byte_code_mutation_doubles_random(x, pmutation)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: pmutation
    integer(i32), allocatable :: bytes(:)
    integer :: i

    bytes = double_vector_to_bytes(x)
    do i = 1, size(bytes)
      if (runif_scalar(0.0_dp, 1.0_dp) < pmutation) bytes(i) = int(randint(0, 255), i32)
    end do
    x = byte_vector_to_doubles(bytes)
  end subroutine byte_code_mutation_doubles_random

  subroutine ensure_bounds(x, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:), upper(:)
    integer :: i

    if (size(lower) /= size(x) .or. size(upper) /= size(x)) error stop "ensure_bounds: shape mismatch"
    do i = 1, size(x)
      if (.not. (x(i) >= lower(i) .and. x(i) <= upper(i))) then
        x(i) = runif_scalar(lower(i), upper(i))
      end if
    end do
  end subroutine ensure_bounds
end module mcga_bytes
