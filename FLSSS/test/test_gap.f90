program test_gap
  use flsss_api
  implicit none
  type(gap_result) :: r, rmin, rg
  real(dp) :: cost(2,3), profit(2,3), budget(2)
  cost=reshape([2.0_dp,4.0_dp, 3.0_dp,2.0_dp, 4.0_dp,1.0_dp],[2,3])
  profit=reshape([8.0_dp,7.0_dp, 6.0_dp,9.0_dp, 7.0_dp,5.0_dp],[2,3])
  budget=[6.0_dp,5.0_dp]
  r=aux_gap_bb(cost,profit,budget,optim="max")
  if(.not.r%feasible .or. abs(r%total_profit_or_loss-24.0_dp)>1e-12_dp) error stop "GAP max"
  if(any(r%assignment/=[1,2,1])) error stop "GAP assignment"
  rmin=aux_gap_bbdp(cost,profit,budget,optim="min")
  if(.not.rmin%feasible .or. abs(rmin%total_profit_or_loss-18.0_dp)>1e-12_dp) error stop "GAP min"
  r=gap_solve(cost,profit,budget)
  if(abs(r%total_profit_or_loss-24.0_dp)>1e-12_dp) error stop "GAP high-level"
  rg=aux_gap_ga(cost,profit,budget,trials=2,population_size=30,generations=60,seed=42_i8)
  if(.not.rg%feasible) error stop "GAP GA feasibility"
  print *, "test_gap: PASS"
end program test_gap
