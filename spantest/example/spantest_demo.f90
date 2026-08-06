! SPDX-License-Identifier: GPL-3.0-only
program spantest_demo
  use spantest
  implicit none
  type(simulation_result) :: sim
  type(span_result) :: grs,hk,km
  type(gl_result) :: gl
  type(as_result) :: robust
  integer :: i

  sim=span_simulate(n=250,k=3,n_test=8,ncp=0.20_dp,dgp=7,seed=2026)
  if (sim%status/=span_ok) error stop trim(sim%message)

  grs=span_grs(sim%r1,sim%r2)
  hk=span_hk(sim%r1,sim%r2)
  km=span_km(sim%r1,sim%r2)
  gl=span_gl_ad(sim%r1,sim%r2,totsim=199,seed=2026)
  robust=span_as(sim%r1,sim%r2,b_draws=20,seed=2026)

  write(*,'(a)') 'spantest modern Fortran demonstration'
  write(*,'(a,es13.5,a,es13.5)') 'GRS statistic = ',grs%stat,', p = ',grs%pval
  write(*,'(a,es13.5,a,es13.5)') 'HK statistic  = ',hk%stat,', p = ',hk%pval
  write(*,'(a,es13.5,a,es13.5)') 'KM statistic  = ',km%stat,', p = ',km%pval
  write(*,'(a,es13.5,a,es13.5)') 'GL joint LMC p = ',gl%pval_lmc,', BMC p = ',gl%pval_bmc
  write(*,'(a,a)') 'GL decision: ',trim(gl%decision_string)
  do i=1,size(robust%pvalues)
    write(*,'(a,1x,es13.5)') trim(robust%names(i)),robust%pvalues(i)
  end do
end program spantest_demo
