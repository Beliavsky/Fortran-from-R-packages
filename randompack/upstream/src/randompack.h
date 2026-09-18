// randompack.h – see https://github.com/jonasson2/randompack

#ifndef RANDOMPACK_H
#define RANDOMPACK_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define RANDOMPACK_VERSION "0.1.10"

typedef struct randompack_rng randompack_rng;

#ifdef __cplusplus
extern "C" {
#endif

//================================== CREATION AND SETUP ==================================

randompack_rng *randompack_create( // Create randomized RNG with given engine, error → 0
  const char *engine    // in      Engine identifier, one of:
  );                    //           x256++simd   x256++   xoro++    philox   ranlux++
//								//				 x256**simd	  x256**	  pcg64     sfc64 	chacha20
//								//				 sfc64simd 	  x128+ 	  squares   cwg128
//                      //         randompack_create(0) gives he default, x256++simd.

bool randompack_seed( // Seed RNG deterministically, false on error
  int seed,             // in      Any integer seed; expanded with a hash to fill state
  uint32_t *spawn_key,  // in      Optional spawn key array (may be 0 if n_key==0)
  int n_key,            // in      Number of spawn_key entries
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_randomize( // Randomize RNG state from system entropy, false on error
  randompack_rng *rng   // in/out  Random number generator
);

void randompack_free( // Free an RNG created with randompack_create
  randompack_rng *rng   // in      Random number generator
);

randompack_rng *randompack_duplicate( // Clone an RNG (identical engine+state), error → 0
  randompack_rng *rng   // in      RNG to duplicate
);

//================================== ENGINE INFORMATION ==================================

bool randompack_engines( // Return supported RNG engine names and their descriptions
  char *engines,      // out     n×eng_maxlen buffer of engine names (0 → query mode)
  char *descriptions, // out     n×desc_maxlen buffer of descriptions (or 0 in query mode)
  int *nengines,      // in/out  number of engines (= n, query → fill, non-query → check)
  int *eng_maxlen,    // in/out  maxi engine-name length including trailing 0
  int *desc_maxlen    // in/out  maxi description length including trailing 0
);

char *randompack_last_error( // Get last error string, or 0 if none
  randompack_rng *rng   // in      Random number generator
);

//=================================== STREAM SELECTION ===================================

bool randompack_jump( // Jump pcg64 / ranlux / xor-family rng by 2^p steps, false on error
  int p,               // in      Jump exponent, ≤128 for pcg, otherwise 32/64/96/128/192
  randompack_rng *rng  // in/out  Random number generator
);

bool randompack_advance( // Advance pcg64 by arbitrary 128-bit delta
  uint64_t delta[2],     // in      128-bit advance delta {low, high}
  randompack_rng *rng    // in/out  Random number generator
);

bool randompack_pcg64_set_inc( // Set PCG increment (state unchanged)
  uint64_t inc[2],     // in      128-bit increment {low, high}; low must be odd
  randompack_rng *rng  // in/out  target RNG
);

bool randompack_cwg128_set_weyl( // Set CWG128 Weyl increment and skip 96 states
  uint64_t weyl[2],    // in      128-bit Weyl increment {low, high}; low must be odd
  randompack_rng *rng  // in/out  target RNG
);

bool randompack_chacha_set_nonce( // Set ChaCha20 nonce (state otherwise unchanged)
  uint32_t nonce[3],   // in      96-bit nonce {w0, w1, w2}
  randompack_rng *rng  // in/out  target RNG
);

bool randompack_philox_set_key( // Set philox key state directly
  uint64_t key[2],     // in      key state
  randompack_rng *rng  // in/out  target RNG
);

bool randompack_squares_set_key( // Set squares64 key state directly
  uint64_t key,        // in      key state
  randompack_rng *rng  // in/out  target RNG
);

bool randompack_sfc64_set_abc( // Set the sfc64 a,b,c state words and warm up
  uint64_t abc[3],     // in      state words a, b, c
  randompack_rng *rng  // in/out  target RNG
);

//==================================== CONFIGURATION =====================================

bool randompack_full_mantissa( // Use full mantissa (default false)
  randompack_rng *rng,  // in/out  Random number generator
  bool enable           // in      true → 24/53-bit mantissa, false → 23/52-bit
);

bool randompack_bitexact( // Use bitexact log/exp (default false)
  randompack_rng *rng,  // in/out  Random number generator
  bool enable           // in      true → bitexact, false → fast, vectorized
);

bool randompack_set_state( // Set state of general engine directly
  uint64_t state[],     // in      state words (must be nonzero for xor-family)
  int nstate,           // in      number of state words provided, depends on engine
  randompack_rng *rng   // in/out  target RNG
);

//==================================== SERIALIZATION =====================================

bool randompack_serialize( // Serialize an RNG to an opaque byte buffer
  uint8_t *buf,          // out     buffer for serialization (may be 0 if *len==0)
  int *len,              // in/out  0 → query size; otherwise buffer length
  randompack_rng *rng    // in      RNG to serialize
);

bool randompack_deserialize( // Restore an RNG from an opaque byte buffer
  const uint8_t *buf,    // in      buffer with serialization
  int len,               // in      buffer length
  randompack_rng *rng    // out     target RNG (must be allocated)
);

//==================================== RAW BITSTREAMS ====================================

bool randompack_raw( // Generate random bytes (raw bitstream), false on error
  void *out,            // out     vector of bytes (compatible with any type)
  size_t nbytes,        // in      number of bytes requested
  randompack_rng *rng   // in/out  random number generator
);

bool randompack_uint8( // Generate uint8 in [0, bound), false on error
  uint8_t x[],         // out     vector of bytes
  size_t len,          // in      number requested
  uint8_t bound,       // in      exclusive upper bound, or 0 for unbounded
  randompack_rng *rng  // in/out  random number generator
);

bool randompack_uint16( // Generate uint16 in [0, bound), false on error
  uint16_t x[],        // out     vector of short integers
  size_t len,          // in      number requested
  uint16_t bound,      // in      exclusive upper bound, or 0 for unbounded
  randompack_rng *rng  // in/out  random number generator
);

bool randompack_uint32( // Generate uint32 in [0, bound), false on error
  uint32_t x[],        // out     vector of integers
  size_t len,          // in      number requested
  uint32_t bound,      // in      exclusive upper bound, or 0 for unbounded
  randompack_rng *rng  // in/out  random number generator
);

bool randompack_uint64( // Generate uint64 in [0, bound), false on error
  uint64_t x[],        // out     vector of uint64_t integers
  size_t len,          // in      number requested
  uint64_t bound,      // in      exclusive upper bound, or 0 for unbounded
  randompack_rng *rng  // in/out  random number generator
);

//================================ DISCRETE DISTRIBUTIONS ================================

bool randompack_int( // Generate uniform integers in [m, n], false on error
  int x[],              // out     vector of integers
  size_t len,           // in      Number of integers requested
  int m,                // in      Inclusive minimum
  int n,                // in      Inclusive maximum
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_long_long( // Generate uniform long long in [m, n], false on error
  long long x[],        // out     vector of integers
  size_t len,           // in      Number of integers requested
  long long m,          // in      Inclusive minimum
  long long n,          // in      Inclusive maximum
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_perm( // Generate a random permutation of 0..n-1, false on error
  int x[],              // out     vector containing the permutation
  int n,                // in      Permutation size
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_sample( // Sample without replacement from 0..n-1, false on error
  int x[],              // out     k-vector of sampled indices
  int n,                // in      Population size
  int k,                // in      Sample size (0 <= k <= n)
  randompack_rng *rng   // in/out  Random number generator
);

//=============================== CONTINUOUS DISTRIBUTIONS ===============================

bool randompack_u01( // Generate uniform random numbers in [0,1), false on error
  double x[],           // out     vector: uniform random numbers in [0,1)
  size_t len,           // in      Number of variates
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_unif( // Generate uniform random numbers in [a,b], false on error
  double x[],           // out     vector: uniform random numbers in [a,b]
  size_t len,           // in      Number of variates
  double a,             // in      Minimum (a < b)
  double b,             // in      Maximum (a < b)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_norm( // Generate standard normal random numbers N(0,1), false on error
  double x[],           // out     vector: standard normal random numbers
  size_t len,           // in      Number of variates
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_normal( // Generate normal random numbers N(mu,sigma), false on error
  double x[],           // out     vector: normal random numbers
  size_t len,           // in      Number of variates
  double mu,            // in      Mean
  double sigma,         // in      Standard deviation (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_exp( // Generate exponential random numbers, false on error
  double x[],           // out     vector: exponential random numbers
  size_t len,           // in      Number of variates
  double scale,         // in      Scale parameter (1.0 → standard exponential)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_lognormal( // Generate lognormal random numbers, false on error
  double x[],           // out     vector: lognormal random numbers
  size_t len,           // in      Number of variates
  double mu,            // in      Mean of underlying normal
  double sigma,         // in      Std dev of underlying normal (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_gamma( // Generate gamma random numbers, false on error
  double x[],           // out     vector: gamma random numbers
  size_t len,           // in      Number of variates
  double shape,         // in      Shape parameter (> 0)
  double scale,         // in      Scale parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_beta( // Generate beta random numbers, false on error
  double x[],           // out     vector: beta random numbers
  size_t len,           // in      Number of variates
  double a,             // in      Shape parameter a (> 0)
  double b,             // in      Shape parameter b (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_chi2( // Generate chi-square random numbers, false on error
  double x[],           // out     vector: chi-square random numbers
  size_t len,           // in      Number of variates
  double nu,            // in      Degrees of freedom (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_t( // Generate t random numbers, false on error
  double x[],           // out     vector: t random numbers
  size_t len,           // in      Number of variates
  double nu,            // in      Degrees of freedom (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_f( // Generate F random numbers, false on error
  double x[],           // out     vector: F random numbers
  size_t len,           // in      Number of variates
  double nu1,           // in      Numerator degrees of freedom (> 0)
  double nu2,           // in      Denominator degrees of freedom (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_gumbel( // Generate gumbel random numbers, false on error
  double x[],           // out     vector: gumbel random numbers
  size_t len,           // in      Number of variates
  double mu,            // in      Location parameter
  double beta,          // in      Scale parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_pareto( // Generate pareto random numbers, false on error
  double x[],           // out     vector: pareto random numbers
  size_t len,           // in      Number of variates
  double xm,            // in      Scale parameter (> 0)
  double alpha,         // in      Shape parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_weibull( // Generate weibull random numbers, false on error
  double x[],           // out     vector: weibull random numbers
  size_t len,           // in      Number of variates
  double shape,         // in      Shape parameter (> 0)
  double scale,         // in      Scale parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_skew_normal( // Generate skew-normal random numbers, false on error
  double x[],           // out     vector: skew-normal random numbers
  size_t len,           // in      Number of variates
  double mu,            // in      Location parameter
  double sigma,         // in      Scale parameter (> 0)
  double alpha,         // in      Skew parameter
  randompack_rng *rng   // in/out  Random number generator
);

//================================= MULTIVARIATE NORMAL ==================================

bool randompack_mvn( // Generate multivariate normal randoms N(mu,Sig), false on error
  char *transp,         // in      "N" to get n×d X, "T" to get d×n X
  double mu[],          // in      d-vector: mean (0 → zero-mean)
  double Sig[],         // in      d×d covariance matrix (0 → use L as-is)
  int d,                // in      Dimension of each vector
  size_t len,           // in      Number of replicates
  double X[],           // out     n×d or d×n matrix of generated vectors
  int ldx,              // in      Leading dimension of X
  double L[],           // in/out  d×d lower Cholesky factor of Sig (or 0)
  randompack_rng *rng   // in/out  Random number generator
);

//=========================== CONTINUOUS DISTRIBUTIONS (FLOAT) ===========================

bool randompack_u01f( // Generate uniform random floats in [0,1), false on error
  float x[],            // out     vector: uniform random numbers in [0,1)
  size_t len,           // in      Number of variates
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_uniff( // Generate uniform random floats U(a,b), false on error
  float x[],            // out     vector: uniform random numbers in [a,b)
  size_t len,           // in      Number of variates
  float a,              // in      Lower bound
  float b,              // in      Upper bound (> a)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_normf( // Generate normal random floats N(0,1), false on error
  float x[],            // out     vector: normal random numbers
  size_t len,           // in      Number of variates
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_normalf( // Generate normal random floats N(mu,sigma), false on error
  float x[],            // out     vector: normal random numbers
  size_t len,           // in      Number of variates
  float mu,             // in      Mean
  float sigma,          // in      Standard deviation (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_expf( // Generate exponential random floats, false on error
  float x[],            // out     vector: exponential random numbers
  size_t len,           // in      Number of variates
  float scale,          // in      Scale parameter (1.0f → standard exponential)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_lognormalf( // Generate lognormal random floats, false on error
  float x[],            // out     vector: lognormal random numbers
  size_t len,           // in      Number of variates
  float mu,             // in      Mean of underlying normal
  float sigma,          // in      Std dev of underlying normal (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_gammaf( // Generate gamma random floats, false on error
  float x[],            // out     vector: gamma random numbers
  size_t len,           // in      Number of variates
  float shape,          // in      Shape parameter (> 0)
  float scale,          // in      Scale parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_betaf( // Generate beta random floats, false on error
  float x[],            // out     vector: beta random numbers
  size_t len,           // in      Number of variates
  float a,              // in      First shape parameter (> 0)
  float b,              // in      Second shape parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_chi2f( // Generate chi-square random floats, false on error
  float x[],            // out     vector: chi-square random numbers
  size_t len,           // in      Number of variates
  float nu,             // in      Degrees of freedom (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_tf( // Generate t random floats, false on error
  float x[],            // out     vector: t random numbers
  size_t len,           // in      Number of variates
  float nu,             // in      Degrees of freedom (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_ff( // Generate F random floats, false on error
  float x[],            // out     vector: F random numbers
  size_t len,           // in      Number of variates
  float nu1,            // in      Numerator degrees of freedom (> 0)
  float nu2,            // in      Denominator degrees of freedom (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_gumbelf( // Generate gumbel random floats, false on error
  float x[],            // out     vector: gumbel random numbers
  size_t len,           // in      Number of variates
  float mu,             // in      Location parameter
  float beta,           // in      Scale parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_paretof( // Generate pareto random floats, false on error
  float x[],            // out     vector: pareto random numbers
  size_t len,           // in      Number of variates
  float xm,             // in      Scale parameter (> 0)
  float alpha,          // in      Shape parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_weibullf( // Generate weibull random floats, false on error
  float x[],            // out     vector: weibull random numbers
  size_t len,           // in      Number of variates
  float shape,          // in      Shape parameter (> 0)
  float scale,          // in      Scale parameter (> 0)
  randompack_rng *rng   // in/out  Random number generator
);

bool randompack_skew_normalf( // Generate skew-normal random floats, false on error
  float x[],            // out     vector: skew-normal random numbers
  size_t len,           // in      Number of variates
  float mu,             // in      Location parameter
  float sigma,          // in      Scale parameter (> 0)
  float alpha,          // in      Skew parameter
  randompack_rng *rng   // in/out  Random number generator
);

#ifdef __cplusplus
}
#endif
#endif /* RANDOMPACK_H */
