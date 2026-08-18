! SPDX-License-Identifier: GPL-3.0-only
program example_poisson_binomial
  use poisson_binomial, only : dp, gpb_table, dpbinom, qpbinom, dgpbinom
  implicit none
  real(dp), parameter :: probs(4) = [0.1_dp, 0.3_dp, 0.6_dp, 0.8_dp]
  integer, parameter :: val_p(3) = [3, -1, 5]
  integer, parameter :: val_q(3) = [0, 2, 1]
  real(dp), parameter :: gp(3) = [0.2_dp, 0.7_dp, 0.4_dp]
  real(dp), allocatable :: pmf(:)
  type(gpb_table) :: gtab
  integer :: k

  pmf = dpbinom(probs, "DivideFFT")
  print '(a)', 'Ordinary Poisson-binomial PMF'
  do k = 0, size(pmf)-1
    print '(i3,2x,f12.8)', k, pmf(k+1)
  end do
  print '(a,i0)', 'Median: ', qpbinom(0.5_dp, probs, "Convolve")

  gtab = dgpbinom(gp, val_p, val_q, "Characteristic")
  print '(/,a)', 'Generalized Poisson-binomial PMF (nonzero entries)'
  do k = gtab%lower, gtab%upper
    if (gtab%values(k-gtab%lower+1) > 1.0e-14_dp) then
      print '(i4,2x,f12.8)', k, gtab%values(k-gtab%lower+1)
    end if
  end do
end program example_poisson_binomial
