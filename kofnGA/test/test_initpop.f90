program test_initpop
  use kofnga, only : dp, i64, kofnga_control, kofnga_result, kofn_ga, summarize_result, kofnga_summary
  implicit none
  type(kofnga_control)::ctl
  type(kofnga_result)::res
  type(kofnga_summary)::s
  integer :: init(6,2)
  init = reshape([1,2, 1,3, 1,4, 2,3, 2,4, 3,4],[6,2],order=[2,1])
  ctl%popsize=6; ctl%keepbest=2; ctl%ngen=10; ctl%tourneysize=2; ctl%mutprob=0.0_dp; ctl%seed=7_i64
  call kofn_ga(4,2,obj,res,ctl,init)
  if(abs(res%bestobj-3.0_dp)>1.0e-14_dp) error stop 'initpop best objective'
  if(any(res%bestsol/=[1,2])) error stop 'initpop best solution'
  s=summarize_result(res)
  if(s%generations/=10) error stop 'summary generations'
  if(s%unique_final<1 .or. s%unique_final>6) error stop 'summary unique'
  print *, 'test_initpop: PASS'
contains
  function obj(x) result(v)
    integer,intent(in)::x(:)
    real(dp)::v
    v=real(sum(x),dp)
  end function obj
end program test_initpop
