program basic_usage
  use rmpfr
  implicit none

  integer, parameter :: prec_bits = 256
  type(mpfr_real) :: pi_value, probability, quantile
  type(mpfr_root_result) :: qres

  pi_value = mpfr_const_pi(prec_bits)
  probability = mpfr_from_string('0.975', prec_bits)
  qres = mpfr_qnorm(probability, tol=mpfr_from_string('1e-60', prec_bits))
  quantile = qres%root

  print '(a)', 'MPFR version: ' // trim(mpfr_version())
  print '(a,i0)', 'precision bits: ', mpfr_precision(pi_value)
  print '(a)', 'pi: ' // trim(mpfr_to_string(pi_value, 70))
  print '(a)', 'qnorm(0.975): ' // trim(mpfr_to_string(quantile, 55))
  print '(a)', 'Bernoulli(10): ' // trim(mpfr_to_string(mpfr_bernoulli(10, prec_bits), 40))
end program basic_usage
