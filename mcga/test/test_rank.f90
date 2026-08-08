program test_rank
  use mcga, only : dp, calculate_rank_scores
  implicit none
  real(dp) :: costs(2,3), ranks(3)

  costs(:,1) = [1.0_dp, 5.0_dp]
  costs(:,2) = [2.0_dp, 2.0_dp]
  costs(:,3) = [5.0_dp, 1.0_dp]
  call calculate_rank_scores(costs, ranks)
  ! Upstream score: each point is better than each other point in at least
  ! one objective, so all receive two points.
  if (any(abs(ranks - 2.0_dp) > 0.0_dp)) error stop "rank score mismatch"
  print *, "test_rank: PASS"
end program test_rank
