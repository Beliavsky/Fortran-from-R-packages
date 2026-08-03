! SPDX-License-Identifier: BSD-3-Clause
program test_information_retrieval
    use metrics, only : dp, f1, apk, mapk, integer_vector, string_vector
    use test_support
    implicit none

    type(integer_vector) :: actual_i(3), predicted_i(3)
    type(string_vector) :: actual_s(2), predicted_s(2)

    call check_close('f1 integer', f1([3, 4, 5], [3, 4]), 0.8_dp)
    call check_close('f1 no match', f1([7], [1, 1]), 0.0_dp)
    call check_close('f1 string', f1(['a', 'c', 'd'], ['d', 'e', 'x']), 1.0_dp / 3.0_dp)
    call check_close('apk 1', apk(2, [1, 2, 3, 4, 5], [6, 4, 7, 1, 2]), 0.25_dp)
    call check_close('apk duplicate', apk(5, [1, 2, 3, 4, 5], [1, 1, 1, 1, 1]), 0.2_dp)
    call check_close('apk limit', apk(3, [1, 3], [1, 2, 3, 4, 5]), 5.0_dp / 6.0_dp)

    actual_i(1)%values = [1, 3, 4]
    actual_i(2)%values = [1, 2, 4]
    actual_i(3)%values = [1, 3]
    predicted_i(1)%values = [1, 2, 3, 4, 5]
    predicted_i(2)%values = [1, 2, 3, 4, 5]
    predicted_i(3)%values = [1, 2, 3, 4, 5]
    call check_close('mapk integer', mapk(3, actual_i, predicted_i), 0.685185185185185_dp, 1.0e-14_dp)

    allocate(character(len=1) :: actual_s(1)%values(2), predicted_s(1)%values(3))
    actual_s(1)%values = ['a', 'b']
    predicted_s(1)%values = ['a', 'c', 'b']
    allocate(character(len=1) :: actual_s(2)%values(1), predicted_s(2)%values(3))
    actual_s(2)%values = ['x']
    predicted_s(2)%values = ['y', 'x', 'z']
    call check_close('mapk string', mapk(3, actual_s, predicted_s), &
                     0.5_dp * ((1.0_dp + 2.0_dp / 3.0_dp) / 2.0_dp + 0.5_dp))

    call finish_tests('test_information_retrieval')
end program test_information_retrieval
