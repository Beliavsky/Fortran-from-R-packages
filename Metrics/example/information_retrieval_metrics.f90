! SPDX-License-Identifier: BSD-3-Clause
program information_retrieval_metrics
    use metrics, only : f1, apk, mapk, integer_vector
    implicit none

    type(integer_vector) :: actual(2), predicted(2)

    actual(1)%values = [1, 3, 4]
    actual(2)%values = [1, 2, 4]
    predicted(1)%values = [1, 2, 3, 4, 5]
    predicted(2)%values = [1, 2, 3, 4, 5]

    write(*, '(a, f10.6)') 'Set F1: ', f1([3, 4, 5], [3, 4])
    write(*, '(a, f10.6)') 'AP@3:   ', apk(3, actual(1)%values, predicted(1)%values)
    write(*, '(a, f10.6)') 'MAP@3:  ', mapk(3, actual, predicted)
end program information_retrieval_metrics
