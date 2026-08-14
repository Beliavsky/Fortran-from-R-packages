program test_core
  use kofnga, only : dp, i64, kofnga_control, kofnga_result, kofn_ga, mutation_probability
  implicit none
  real(dp), parameter :: vals(20) = [ &
    0.01_dp,0.02_dp,0.03_dp,0.04_dp,0.5_dp,0.6_dp,0.7_dp,0.8_dp,0.9_dp,1.0_dp, &
    1.1_dp,1.2_dp,1.3_dp,1.4_dp,1.5_dp,1.6_dp,1.7_dp,1.8_dp,1.9_dp,2.0_dp]
  type(kofnga_control) :: ctl
  type(kofnga_result) :: res
  integer :: i,j
  ctl%popsize=80; ctl%keepbest=8; ctl%ngen=80; ctl%tourneysize=8
  ctl%mutprob=0.03_dp; ctl%seed=123456_i64
  call kofn_ga(20,4,obj,res,ctl)
  if (any(res%bestsol /= [1,2,3,4])) error stop 'best subset mismatch'
  if (abs(res%bestobj-0.10_dp)>1.0e-12_dp) error stop 'best objective mismatch'
  do i=1,size(res%pop,1)
    if(any(res%pop(i,:)<1) .or. any(res%pop(i,:)>20)) error stop 'range failure'
    do j=2,size(res%pop,2)
      if(res%pop(i,j)<=res%pop(i,j-1)) error stop 'row not sorted/unique'
    end do
  end do
  do i=2,size(res%obj_history)
    if(res%obj_history(i)>res%obj_history(i-1)+1.0e-14_dp) error stop 'elitist best worsened'
  end do
  if(abs(mutation_probability(4,0.5_dp)-(1.0_dp-0.5_dp**0.25_dp))>1.0e-15_dp) error stop 'mutfrac mapping'
  print *, 'test_core: PASS'
contains
  function obj(s) result(v)
    integer,intent(in)::s(:)
    real(dp)::v
    v=sum(vals(s))
  end function obj
end program test_core
