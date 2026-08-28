program test_randomization
use qrng, only: dp, set_seed_int, korobov, ghalton, sobol
implicit none
real(dp), allocatable :: a(:,:), b(:,:)
real(dp) :: shift(2)

call set_seed_int(271)
a = korobov(17,2,3,randomize=.true.)
call set_seed_int(271)
b = korobov(17,2,3,randomize=.true.)
if (any(a /= b)) error stop "korobov randomization reproducibility"
if (any(a < 0.0_dp) .or. any(a > 1.0_dp)) error stop "korobov randomization range"

call set_seed_int(271)
a = ghalton(20,4)
call set_seed_int(271)
b = ghalton(20,4)
if (any(a /= b)) error stop "ghalton randomization reproducibility"
if (any(a < 0.0_dp) .or. any(a >= 1.0_dp)) error stop "ghalton range"

shift = [0.1_dp,0.2_dp]
a = sobol(16,2,randomize=.true.,digital_shift=shift)
b = sobol(16,2,randomize=.true.,digital_shift=shift)
if (any(a /= b)) error stop "sobol fixed digital shift reproducibility"
if (any(a < 0.0_dp) .or. any(a >= 1.0_dp)) error stop "sobol randomized range"

print '(a)', 'test_randomization: PASS'
end program test_randomization
