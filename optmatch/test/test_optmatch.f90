program test_optmatch
   use optmatch_kinds, only : dp
   use optmatch_types, only : distance_spec, match_result
   use optmatch_distance, only : dense_distance, score_distance, euclidean_distance
   use optmatch_distance, only : mahalanobis_distance, rank_mahalanobis_distance
   use optmatch_distance, only : caliper_distance, exact_match_distance, add_distances
   use optmatch_distance, only : to_sparse, from_sparse, num_eligible_matches
   use optmatch_matching, only : fullmatch, pairmatch
   use optmatch_utilities, only : effective_sample_size, standardization_scale, integer_digits
   use optmatch_feasibility, only : caliper_size, max_controls_cap, min_controls_cap
   implicit none
   type(distance_spec) :: d, c, e, combined, d2
   type(match_result) :: m
   real(dp) :: values(2,2), scores(4), data(4,2)
   logical :: z(4)
   integer :: block(4)
   real(dp) :: cap

   values = reshape([1.0_dp, 5.0_dp, 4.0_dp, 1.0_dp], [2,2])
   d = dense_distance(values)
   m = pairmatch(d)
   call assert_true(m%feasible, 'pairmatch feasible')
   call assert_true(m%n_selected == 2, 'pairmatch selects two edges')
   call assert_close(m%objective, 2.0_dp, 1.0e-10_dp, 'pairmatch objective')
   call assert_true(all(m%treatment_group > 0), 'all treatments matched')
   call assert_true(all(m%control_group > 0), 'all controls matched')

   m = fullmatch(d)
   call assert_true(m%feasible, 'fullmatch feasible')
   call assert_close(m%objective, 2.0_dp, 1.0e-10_dp, 'fullmatch objective')

   z = [.true., .true., .false., .false.]
   scores = [0.0_dp, 10.0_dp, 1.0_dp, 9.0_dp]
   d = score_distance(scores, z)
   call assert_close(d%value(1,1), 1.0_dp, 1.0e-12_dp, 'score distance')
   c = caliper_distance(d, 2.0_dp)
   call assert_true(count(c%allowed) == 2, 'caliper keeps two edges')

   block = [1,2,1,2]
   e = exact_match_distance(block, z)
   combined = add_distances(d, e)
   call assert_true(count(combined%allowed) == 2, 'exact matching restricts edges')
   m = pairmatch(combined)
   call assert_true(m%feasible, 'exact pairmatch feasible')
   call assert_close(m%objective, 2.0_dp, 1.0e-10_dp, 'exact pairmatch objective')

   d2 = from_sparse(to_sparse(combined))
   call assert_true(num_eligible_matches(d2) == 2, 'sparse round trip')

   data = reshape([0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp, &
                   0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [4,2])
   d = euclidean_distance(data, z)
   call assert_close(d%value(1,1), 3.0_dp, 1.0e-10_dp, 'euclidean distance')

   data = reshape([0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp, &
                   0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [4,2])
   d = mahalanobis_distance(data, z)
   call assert_true(all(d%value >= 0.0_dp), 'Mahalanobis nonnegative')
   d = rank_mahalanobis_distance(data, z)
   call assert_true(all(d%value >= 0.0_dp), 'rank Mahalanobis nonnegative')

   values = reshape([1.0_dp, 5.0_dp, 4.0_dp, 1.0_dp], [2,2])
   d = dense_distance(values)
   m = pairmatch(d)
   call assert_close(effective_sample_size(m), 2.0_dp, 1.0e-10_dp, 'effective sample size')
   call assert_true(caliper_size(scores, z, 2.0_dp) == 2, 'caliper size')
   call assert_true(integer_digits(-12345) == 6, 'integer digits')
   call assert_true(standardization_scale(scores,z,use_mad=.false.) > 0.0_dp, 'standardization scale')
   cap = max_controls_cap(d)
   call assert_close(cap, 1.0_dp, 1.0e-10_dp, 'max controls cap')
   cap = min_controls_cap(d)
   call assert_close(cap, 1.0_dp, 1.0e-10_dp, 'min controls cap')

   d = dense_distance(reshape([1.0_dp,2.0_dp,3.0_dp], [1,3]))
   m = fullmatch(d)
   call assert_true(m%feasible .and. m%n_selected == 3, 'one-to-many full match')
   cap = max_controls_cap(d)
   call assert_close(cap, 3.0_dp, 1.0e-10_dp, 'one-to-many max cap')
   cap = min_controls_cap(d)
   call assert_close(cap, 3.0_dp, 1.0e-10_dp, 'one-to-many min cap')

   d = dense_distance(reshape([1.0_dp,2.0_dp], [2,1]))
   m = fullmatch(d)
   call assert_true(m%feasible .and. m%n_selected == 2, 'many-to-one full match')

   print '(a)', 'All optmatch-fortran tests passed.'

contains

subroutine assert_true(condition, message)
   logical, intent(in) :: condition
   character(len=*), intent(in) :: message
   if (.not. condition) then
      write(*, '(a)') 'FAIL: '//message
      error stop 1
   end if
end subroutine assert_true

subroutine assert_close(x, y, tol, message)
   real(dp), intent(in) :: x, y, tol
   character(len=*), intent(in) :: message
   call assert_true(abs(x-y) <= tol, message)
end subroutine assert_close

end program test_optmatch
