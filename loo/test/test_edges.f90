! SPDX-License-Identifier: GPL-3.0-or-later
program test_edges
   use loo_api
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_is_nan, ieee_is_finite, ieee_negative_inf
   implicit none
   integer :: perm(10), cats(10), groups(10), folds(10), i
   integer, allocatable :: ids(:)
   real(dp) :: x(100), grid(100), tab(3, 3), tiny_eff(1)
   real(dp), allocatable :: w(:, :)
   type(importance_sampling_result) :: object
   type(expectation_result) :: out
   type(loo_result) :: result
   perm = [10, 1, 3, 7, 2, 9, 5, 4, 8, 6]
   cats = [-1, -1, 8, 8, -1, 3, 8, 3, 3, 3]
   groups = [8, -1, 2, 8, 10, -1, 2, 1, 10, 1]
   folds = kfold_split_random(3, perm)
   call fold_reference('random', folds)
   folds = kfold_split_stratified(3, cats, perm)
   call fold_reference('stratified', folds)
   folds = kfold_split_grouped(3, groups, [5, 2, 1, 4, 3])
   call fold_reference('grouped', folds)
   do i = 1, 10
      call require(all(pack(folds, groups == groups(i)) == folds(i)), 'group stays together')
   end do
   folds = kfold_split_grouped(5, groups, [1])
   call require(all(folds == [4, 1, 3, 4, 5, 1, 3, 2, 5, 2]), 'one group per fold')
   folds = kfold_split_random(1, perm)
   call require(all(folds == 0), 'invalid fold count')
   folds = kfold_split_random(3, [1, 1, 2, 3, 4, 5, 6, 7, 8, 9])
   call require(all(folds == 0), 'invalid permutation')
   x = [(real(i, dp), i = 1,100)]
   grid = 0.0_dp
   object = psis(grid)
   out = E_loo(x, object, 'quantile', [0.1_dp, 0.5_dp, 0.9_dp])
   call require(maxval(abs(out%value(:, 1) - [10.9_dp, 50.5_dp, 90.1_dp])) < 1.0e-12_dp, 'equal weight type 7 quantiles')
   out = E_loo(x, object, 'mean')
   call require(abs(out%value(1, 1) - 50.5_dp) < 1.0e-12_dp, 'Eloo vector mean')
   out = E_loo(x, object, 'variance')
   call require(abs(out%value(1, 1) - 841.6666666666667_dp) < 1.0e-10_dp, 'Eloo unbiased uniform variance')
   out = E_loo(x, object, 'quantile', [0.0_dp])
   call require(out%status == 1, 'endpoint quantile rejected')
   out = E_loo(x, object, 'unknown')
   call require(out%status == 1, 'unknown expectation rejected')
   object = psis(log(x))
   tab = pareto_k_table(object)
   call require(sum(tab(:, 1)) == 1.0_dp, 'diagnostic total count')
   call require(abs(sum(tab(:, 2)) - 1.0_dp) < 1.0e-14_dp, 'diagnostic proportion sum')
   ids = pareto_k_ids(object, threshold = -10.0_dp)
   call require(size(ids) == 1, 'k threshold selection')
   ids = pareto_k_ids(object, threshold = 10.0_dp)
   call require(size(ids) == 0, 'no bad k')
   tiny_eff = [tiny(1.0_dp) / 1024.0_dp]
   object = psis(log(x), tiny_eff)
   call require(object%status == 0 .and. object%tail_len(1) == 20, 'subnormal efficiency tail cap')
   call require(ieee_is_finite(object%n_eff(1)), 'subnormal efficiency ESS finite')
   object = psis([0.0_dp])
   call require(object%status == 0, 'single draw psis')
   call require(pareto_threshold(1) == ieee_value(0.0_dp, ieee_negative_inf), 'single draw reliability threshold')
   ids = pareto_k_ids(object)
   call require(size(ids) == 1, 'single draw invalid reliability selected')
   object = psis([(log(real(1 + mod(i, 7), dp)), i = 1,100)])
   w = importance_weights(object, logarithm = .false.)
   call require(object%status == 0 .and. abs(sum(w) - 1.0_dp) < 1.0e-12_dp, 'discrete tied-tail normalization')
   object = psis(log(x), [1.0_dp, 2.0_dp])
   call require(object%status == 1, 'invalid efficiency length')
   result = waic(reshape([1.0_dp], [1,1]))
   call require(result%status == 1, 'WAIC requires two draws')
   result = loo(reshape(log(x), [100,1]))
   call require(all(ieee_is_nan(result%estimates(:, 2))), 'one observation standard errors undefined')
   print *, 'Fold oracles, ties, quantiles, diagnostics, and degenerate-dimension tests passed.'
contains
   subroutine require(condition, label)
      logical, intent(in) :: condition !! Required property of a numerical edge case.
      character(*), intent(in) :: label !! Human-readable description of the failed edge case.
      if (.not. condition) then
         print *, 'FAILED: ', label
         error stop 1
      end if
   end subroutine require
   subroutine fold_reference(name, actual)
      character(*), intent(in) :: name !! Name of the actual-R deterministic-permutation fold fixture.
      integer, intent(in) :: actual(:) !! Returned fold assignments in observation order.
      integer :: expected(size(actual)), unit, k
      open(newunit = unit, file = 'test/reference/data/fold_' // name // '.txt', status = 'old', action = 'read')
      do k = 1, size(actual)
         read(unit, *) expected(k)
      end do
      close(unit)
      call require(all(actual == expected), name // ' fold oracle')
   end subroutine fold_reference
end program test_edges
