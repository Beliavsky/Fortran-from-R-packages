! SPDX-License-Identifier: BSD-3-Clause
program test_classification
    use metrics, only : dp, ce, accuracy, score_quadratic_weighted_kappa, mean_quadratic_weighted_kappa, &
                        scorequadraticweightedkappa, meanquadraticweightedkappa
    use test_support
    implicit none

    call check_close('integer ce', ce([1, 1, 1, 0, 0, 0], [1, 1, 1, 1, 0, 0]), 1.0_dp / 6.0_dp)
    call check_close('real accuracy', accuracy([1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 0.0_dp, 3.0_dp]), 2.0_dp / 3.0_dp)
    call check_close('string ce', ce(['cat ', 'dog ', 'bird'], ['cat ', 'dog ', 'fish']), 1.0_dp / 3.0_dp)
    call check_close('qwk 1', score_quadratic_weighted_kappa([1, 2, 1], [1, 2, 2]), 0.4_dp)
    call check_close('direct qwk alias', scorequadraticweightedkappa([1, 2, 1], [1, 2, 2]), 0.4_dp)
    call check_close('qwk 2', score_quadratic_weighted_kappa([1, 2, 3, 1, 2, 3], [1, 2, 3, 1, 3, 2]), 0.75_dp)
    call check_close('qwk exact', score_quadratic_weighted_kappa([1, 2, 3], [1, 2, 3]), 1.0_dp)
    call check_close('qwk gaps', score_quadratic_weighted_kappa([1, 3, 5], [2, 4, 6]), 0.8421052631578947_dp)
    call check_close('qwk explicit range', score_quadratic_weighted_kappa([1, 3, 3, 5], [2, 4, 5, 6], 1, 6), &
                     0.6956521739130435_dp)
    call check_close('mean qwk positive', mean_quadratic_weighted_kappa([1.0_dp, 1.0_dp]), 0.999_dp)
    call check_close('direct mean qwk alias', meanquadraticweightedkappa([1.0_dp, 1.0_dp]), 0.999_dp)
    call check_close('mean qwk cancelling', mean_quadratic_weighted_kappa([1.0_dp, -1.0_dp]), 0.0_dp, 1.0e-15_dp)
    call check_close('mean qwk weighted', mean_quadratic_weighted_kappa([0.5_dp, 0.8_dp], [1.0_dp, 0.5_dp]), &
                     0.624536446425734_dp, 1.0e-14_dp)

    call finish_tests('test_classification')
end program test_classification
