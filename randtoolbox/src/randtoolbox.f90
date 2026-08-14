module randtoolbox
   use randtoolbox_base, only : dp
   use randtoolbox_bits, only : int2bit, bit2int, bit2unitreal, bit_xor
   use randtoolbox_primes, only : get_primes, nth_prime
   use randtoolbox_quasi, only : torus, halton, sobol, halton_radical_inverse, sobol_directions
   use randtoolbox_pseudo, only : set_seed, congru_rand, sfmt_generate, well_generate, knuth_taocp, mt19937_generate
   use randtoolbox_congruential, only : congru_rng
   use randtoolbox_sfmt, only : sfmt_rng, reset_sfmt_parameter_sets, sfmt_supported
   use randtoolbox_mt19937, only : mt19937_rng
   use randtoolbox_knuth, only : knuth_rng
   use rngwell, only : well_rng, well_from_options, well_state_size, well_variant_supported
   use randtoolbox_tests, only : rng_test_result, gap_test, frequency_test, serial_test, poker_test, order_test, &
      collision_count, collision_test_counts, stirling_second, stirling_divided_by_k, permutations
   implicit none
   public
end module randtoolbox
