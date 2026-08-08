program test_bytes
  use, intrinsic :: iso_fortran_env, only : int64
  use mcga, only : dp, i32, double_to_bytes, bytes_to_double, double_vector_to_bytes, byte_vector_to_doubles, &
                   one_point_crossover, two_point_crossover, byte_code_mutation, set_random_seed
  implicit none
  real(dp) :: x, y
  real(dp), allocatable :: xv(:), yv(:)
  integer(i32) :: bs(storage_size(0.0_dp)/8)
  integer(i32), allocatable :: b(:), b1(:), b2(:), c1(:), c2(:)
  integer :: i

  x = -12345.6789012345_dp
  bs = double_to_bytes(x)
  y = bytes_to_double(bs)
  if (transfer(x, 0_int64) /= transfer(y, 0_int64)) error stop "double byte roundtrip failed"

  xv = [1.25_dp, -2.5_dp, 1.0e100_dp]
  allocate(b(size(xv) * storage_size(0.0_dp)/8))
  b = double_vector_to_bytes(xv)
  yv = byte_vector_to_doubles(b)
  if (size(yv) /= size(xv)) error stop "vector byte size failed"
  do i = 1, size(xv)
    if (transfer(xv(i), 0_int64) /= transfer(yv(i), 0_int64)) error stop "vector byte roundtrip failed"
  end do

  b1 = [1_i32, 2_i32, 3_i32, 4_i32]
  b2 = [11_i32, 12_i32, 13_i32, 14_i32]
  call one_point_crossover(b1, b2, 2, c1, c2)
  if (any(c1 /= [1_i32,2_i32,13_i32,14_i32])) error stop "one point crossover failed"
  if (any(c2 /= [11_i32,12_i32,3_i32,4_i32])) error stop "one point crossover failed"
  call two_point_crossover(b1, b2, 1, 2, c1, c2)
  if (any(c1 /= [1_i32,12_i32,13_i32,4_i32])) error stop "two point crossover failed"

  call set_random_seed(123)
  deallocate(b)
  allocate(b(2))
  b = [0_i32, 255_i32]
  call byte_code_mutation(b, 1.0_dp)
  if (any(b < 0_i32) .or. any(b > 255_i32)) error stop "mutation byte range failed"

  print *, "test_bytes: PASS"
end program test_bytes
