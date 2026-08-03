! SPDX-License-Identifier: BSD-3-Clause
program classification_metrics
    use metrics, only : accuracy, ce, score_quadratic_weighted_kappa
    implicit none

    character(len=4), parameter :: actual(5) = ['a   ', 'a   ', 'c   ', 'b   ', 'c   ']
    character(len=4), parameter :: predicted(5) = ['a   ', 'b   ', 'c   ', 'b   ', 'a   ']

    write(*, '(a, f10.6)') 'Accuracy: ', accuracy(actual, predicted)
    write(*, '(a, f10.6)') 'Error:    ', ce(actual, predicted)
    write(*, '(a, f10.6)') 'QWK:      ', score_quadratic_weighted_kappa([1, 2, 3, 1], [1, 2, 2, 1])
end program classification_metrics
