program test_knapsack
  use flsss_api
  implicit none
  type(knapsack_result) :: r, ri
  type(knapsack_multi_result) :: d, b
  real(dp) :: profits(5), costs(5,2), caps(2)

  profits=[10.0_dp,8.0_dp,7.0_dp,6.0_dp,5.0_dp]
  costs=reshape([2.0_dp,3.0_dp,4.0_dp,5.0_dp,1.0_dp, &
                 3.0_dp,2.0_dp,4.0_dp,1.0_dp,5.0_dp],[5,2])
  caps=[7.0_dp,6.0_dp]
  r=mm_knapsack(2,profits,costs,caps)
  if(.not.r%feasible) error stop "mmKnapsack infeasible"
  if(abs(r%selection_profit-18.0_dp)>1.0e-12_dp) error stop "mmKnapsack objective"
  if(any(r%solution/=[1,2])) error stop "mmKnapsack selection"

  ri=mm_knapsack_integerized(2,profits,costs,caps,precision_level=[16,16])
  if(abs(ri%selection_profit-18.0_dp)>1.0e-12_dp) error stop "mmKnapsackIntegerized"

  d=aux_knapsack01dp([2,3,4,5],[3.0_dp,4.0_dp,5.0_dp,6.0_dp],[5,7])
  if(maxval(abs(d%max_value-[7.0_dp,9.0_dp]))>1.0e-12_dp) error stop "aux DP"
  b=aux_knapsack01bb([2.0_dp,3.0_dp,4.0_dp,5.0_dp], &
                     [3.0_dp,4.0_dp,5.0_dp,6.0_dp],[5.0_dp,7.0_dp])
  if(maxval(abs(b%max_value-[7.0_dp,9.0_dp]))>1.0e-12_dp) error stop "aux BB"
  print *, "test_knapsack: PASS"
end program test_knapsack
