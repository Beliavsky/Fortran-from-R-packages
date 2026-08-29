! Umbrella module for the computational qrng translation.
module qrng
use r_compat, only: dp, set_seed_int
use qrng_korobov_mod, only: korobov
use qrng_ghalton_mod, only: ghalton
use qrng_sobol_mod, only: sobol
use qrng_utils_mod, only: to_array_matrix, to_array_3d
use qrng_test_functions_mod, only: sum_of_squares, sobol_g, exceedance_indicator, &
   exceedance_individual_given_sum, exceedance_sum_given_sum
implicit none
public
end module qrng
