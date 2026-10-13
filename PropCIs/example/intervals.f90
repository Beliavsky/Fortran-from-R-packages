! SPDX-License-Identifier: GPL-2.0-only
program intervals
  use propcis
  implicit none
  type(ci_result) :: ci
  ci = exactci(5, 20, 0.95_dp)
  call show('Exact binomial', ci)
  ci = diffscoreci(10, 40, 18, 50, 0.95_dp)
  call show('Independent difference score', ci)
  ci = riskscoreci(10, 100, 90, 100, 0.95_dp)
  call show('Risk-ratio score', ci)
  ci = orscoreci(10, 40, 18, 50, 0.95_dp)
  call show('Odds-ratio score', ci)
  ci = diffci_bayes(5, 20, 10, 30, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call show('Bayesian difference', ci)
  ci = rrci_bayes(5, 20, 10, 30, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call show('Bayesian risk ratio', ci)
  ci = orci_bayes(5, 20, 30, 30, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call show('Bayesian odds ratio, all successes in sample 2', ci)
  ci = diffci_bayes_hpd(5, 20, 10, 30, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call show('Shortest connected posterior difference', ci)
contains
  subroutine show(label, value)
    character(len=*), intent(in) :: label !! Human-readable method name.
    type(ci_result), intent(in) :: value !! Computed interval, required to have successful numerical status.
    if (value%status /= ci_ok) error stop 'interval computation failed'
    print '(a,2f15.8)', label, value%lower, value%upper
  end subroutine show
end program intervals
