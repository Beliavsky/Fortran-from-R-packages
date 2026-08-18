! SPDX-License-Identifier: GPL-2.0-or-later
program tensor_example
    use tensora
    implicit none

    type(tensor_t) :: a, b, c
    real(dp) :: am(2,3), bm(3,2)

    am = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp], [2,3])
    bm = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp], [3,2])

    a = tensor(reshape(am,[6]), [2,3], ['i','j'])
    b = tensor(reshape(bm,[6]), [3,2], ['j','k'])
    c = einstein_pair(a,b)

    print '(a)', 'A * B using named Einstein contraction:'
    print '(2f10.4)', reshape(real(c%data,dp), [2,2])
end program tensor_example
