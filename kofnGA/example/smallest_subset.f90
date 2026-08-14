program smallest_subset
  use kofnga, only : dp, i64, kofnga_control, kofnga_result, kofn_ga
  implicit none
  real(dp), parameter :: x(10)=[0.05_dp,0.10_dp,0.15_dp,0.20_dp,0.8_dp,0.9_dp,1.0_dp,1.1_dp,1.2_dp,1.3_dp]
  type(kofnga_control)::ctl
  type(kofnga_result)::res
  ctl%popsize=50; ctl%ngen=40; ctl%keepbest=5; ctl%tourneysize=5; ctl%mutprob=0.03_dp; ctl%seed=42_i64
  call kofn_ga(10,3,obj,res,ctl)
  write(*,'(a,*(i0,1x))') 'best subset: ',res%bestsol
  write(*,'(a,f10.6)') 'best objective: ',res%bestobj
contains
  function obj(s) result(v)
    integer,intent(in)::s(:)
    real(dp)::v
    v=sum(x(s))
  end function obj
end program smallest_subset
