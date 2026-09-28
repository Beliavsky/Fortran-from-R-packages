program test_pdfcluster_clustering
   use pdfcluster, only : dbs_from_labels, dbs_from_pdfcluster, dbs_result, dp
   use pdfcluster, only : groups_at_stage, pdfcluster_fit, pdfcluster_regraph_pairs, pdfcluster_result
   implicit none

   type(dbs_result) :: dres
   type(pdfcluster_result) :: pair_fit
   type(pdfcluster_result) :: pair_regraph
   type(pdfcluster_result) :: result
   real(dp) :: x1(10, 1)
   real(dp) :: x2(6, 2)
   integer, allocatable :: core(:)
   integer, allocatable :: groups(:)
   integer :: labels(10)
   integer :: i
   logical :: ok

   x1(:, 1) = [-5.0_dp, -4.9_dp, -4.8_dp, -4.7_dp, -4.6_dp, &
                4.6_dp, 4.7_dp, 4.8_dp, 4.9_dp, 5.0_dp]
   labels = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]

   call pdfcluster_fit(x1, result, ok, n_grid=10, n_stage=0)
   call assert_true(ok, "unidimensional pdfCluster fit")
   call assert_true(result%noc == 2, "two cluster cores are detected")
   call groups_at_stage(result, stage=0, groups=core, success=ok)
   call assert_true(ok, "groups stage zero succeeds")
   call assert_true(all(core(2:4) == core(2)) .and. core(2) > 0, "left core is coherent")
   call assert_true(all(core(6:9) == core(6)) .and. core(6) > 0, "right core is coherent")
   call assert_true(core(2) /= core(6), "cluster cores remain distinct")

   call pdfcluster_fit(x1, result, ok, n_grid=10, n_stage=5)
   call assert_true(ok, "staged pdfClassification succeeds")
   call groups_at_stage(result, groups=groups, success=ok)
   call assert_true(ok, "groups default stage succeeds")
   call assert_true(all(groups(1:5) == groups(1)), "left observations share a final class")
   call assert_true(all(groups(6:10) == groups(6)), "right observations share a final class")
   call assert_true(groups(1) /= groups(6), "final classes remain distinct")

   call dbs_from_labels(x1, labels, dres, ok)
   call assert_true(ok, "dbs from explicit labels succeeds")
   call assert_true(all(dres%values > 0.99_dp), "well-separated groups have high density silhouettes")

   call dbs_from_pdfcluster(result, dres, ok)
   call assert_true(ok, "dbs from pdfCluster result succeeds")
   call assert_true(count(dres%clusters > 0) == 10, "dbs retains all classified observations")

   x2(1, :) = [-3.0_dp, -3.0_dp]
   x2(2, :) = [-2.8_dp, -3.1_dp]
   x2(3, :) = [-3.1_dp, -2.8_dp]
   x2(4, :) = [3.0_dp, 3.0_dp]
   x2(5, :) = [2.8_dp, 3.1_dp]
   x2(6, :) = [3.1_dp, 2.8_dp]

   call pdfcluster_fit(x2, pair_fit, ok, graphtype="pairs", lambda=0.10_dp, grid_pairs=10, n_grid=6, n_stage=0)
   call assert_true(ok, "pair-graph pdfCluster fit")
   call assert_true(pair_fit%noc == 2, "pair graph separates two compact groups")
   call pdfcluster_regraph_pairs(pair_fit, pair_regraph, ok, lambda=1.0_dp, n_grid=6, n_stage=0)
   call assert_true(ok, "pair-graph pdfCluster refit from stored valley areas")
   call assert_true(pair_regraph%noc == 1, "large pair threshold merges the two groups")
   call assert_true(all(pair_regraph%adjacency), "lambda one retains every pair edge in this fixture")

   call pdfcluster_fit(x2, result, ok, graphtype="delaunay", n_grid=6, n_stage=0)
   call assert_true(ok, "Delaunay pdfCluster path succeeds through sibling geometry package")
   call assert_true(all([(result%adjacency(i, i), i=1, size(x2, 1))]), &
                    "Delaunay graph includes every self edge")

   print *, "All pdfCluster clustering tests passed."

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Boolean condition that must evaluate true.
      character(len=*), intent(in) :: message !! Assertion description printed on failure.

      if (.not. condition) then
         print *, "FAIL:", trim(message)
         error stop 1
      end if
   end subroutine assert_true

end program test_pdfcluster_clustering
