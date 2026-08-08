program test_operators
  use mcga, only : dp, set_random_seed, byte_crossover, byte_mutation, sbx_crossover, flat_crossover, &
                   arithmetic_crossover, blx_crossover, unfair_average_crossover, linear_crossover
  implicit none
  real(dp) :: p1(3), p2(3), lo(3), hi(3)
  real(dp), allocatable :: c1(:), c2(:)

  p1 = [1.0_dp, 2.0_dp, 3.0_dp]
  p2 = [4.0_dp, 5.0_dp, 6.0_dp]
  lo = [-10.0_dp, -10.0_dp, -10.0_dp]
  hi = [ 10.0_dp,  10.0_dp,  10.0_dp]
  call set_random_seed(77)

  call arithmetic_crossover(p1, p2, c1, c2)
  if (maxval(abs(c1 + c2 - p1 - p2)) > 1.0e-12_dp) error stop "arithmetic crossover invariant failed"
  call flat_crossover(p1, p2, c1, c2)
  if (any(c1 < p1) .or. any(c1 > p2) .or. any(c2 < p1) .or. any(c2 > p2)) error stop "flat crossover bounds failed"
  call blx_crossover(p1, p2, c1, c2)
  if (size(c1) /= 3 .or. size(c2) /= 3) error stop "BLX shape failed"
  call sbx_crossover(p1, p2, c1, c2)
  if (maxval(abs(c1 + c2 - p1 - p2)) > 1.0e-10_dp) error stop "SBX symmetry failed"
  call unfair_average_crossover(p1, p2, c1, c2, split_index=2, alpha_value=0.25_dp)
  if (maxval(abs(c1 + c2 - p1 - p2)) > 1.0e-12_dp) error stop "unfair average invariant failed"
  call linear_crossover(p1, p2, fitness_max, c1, c2)
  if (size(c1) /= 3 .or. size(c2) /= 3) error stop "linear crossover shape failed"
  call byte_crossover(p1, p2, lo, hi, c1, c2)
  if (any(c1 < lo) .or. any(c1 > hi) .or. any(c2 < lo) .or. any(c2 > hi)) error stop "byte crossover bounds failed"
  c1 = p1
  call byte_mutation(c1, 1.0_dp, lo, hi)
  if (any(c1 < lo) .or. any(c1 > hi)) error stop "byte mutation bounds failed"

  print *, "test_operators: PASS"
contains
  function fitness_max(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = -sum((x - 4.0_dp)**2)
  end function fitness_max
end program test_operators
