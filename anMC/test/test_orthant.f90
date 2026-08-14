program test_orthant
  use anmc
  implicit none
  integer,parameter::d=20
  real(dp)::mu(d),sigma(d,d),e(d,1),target
  integer::i
  type(probability_estimate)::pmc,panmc
  type(simulation_control)::ctl
  mu=0.0_dp; sigma=0.5_dp
  do i=1,d
    sigma(i,i)=1.0_dp
    e(i,1)=real(i-1,dp)/real(d-1,dp)
  end do
  target=real(d,dp)/real(d+1,dp)
  ctl%max_outer=20000; ctl%max_inner=100; ctl%max_rejection_batch=200000; ctl%enforce_budget=.false.
  pmc=proba_max(0.03_dp,0.0_dp,mu,sigma,e,q=5,method=0,algo='MC',sim_control=ctl, &
                prob_control=genz_bretz(maxpts=200000,abseps=1e-5_dp),seed=123)
  if(.not.pmc%ok) then
    write(*,*) 'MC failure: ',trim(pmc%message); error stop 1
  end if
  if(abs(pmc%probability-target)>0.04_dp) then
    write(*,*) 'MC orthant mismatch ',pmc%probability,target; error stop 1
  end if
  panmc=proba_max(0.03_dp,0.0_dp,mu,sigma,e,q=5,method=0,algo='ANMC',sim_control=ctl, &
                  prob_control=genz_bretz(maxpts=200000,abseps=1e-5_dp),seed=321)
  if(.not.panmc%ok) then
    write(*,*) 'ANMC failure: ',trim(panmc%message); error stop 1
  end if
  if(abs(panmc%probability-target)>0.04_dp) then
    write(*,*) 'ANMC orthant mismatch ',panmc%probability,target; error stop 1
  end if
  write(*,'(a,f10.6,a,f10.6,a,f10.6)') 'MC=',pmc%probability,' ANMC=',panmc%probability,' target=',target
  print *, 'test_orthant: PASS'
end program test_orthant
