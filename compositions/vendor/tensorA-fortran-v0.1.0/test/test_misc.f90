! SPDX-License-Identifier: GPL-2.0-or-later
program test_misc
    use tensora
    implicit none

    type(tensor_t) :: x, z, d, dg, op, g, gi, v, up, down
    complex(dp) :: m1(2,2), m2(2,2)
    real(dp), parameter :: tol = 1.0e-10_dp

    x = tensor([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, &
                3.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], &
               [2,2,2], [character(len=5) :: 'row','col','batch'])

    z = to_matrix_tensor(x, ['row'], ['col'], ['batch'])
    if (any(z%shape /= [2,2,2])) error stop 'to_matrix shape'
    if (maxval(abs(z%data-x%data)) > tol) error stop 'to_matrix data'

    op = opnorm_by_tensor(x, ['row'], ['batch'], ['col'])
    if (any(op%shape /= [2])) error stop 'opnorm_by shape'
    if (abs(real(op%data(1),dp)-2.0_dp) > tol) error stop 'opnorm batch 1'
    if (abs(real(op%data(2),dp)-4.0_dp) > tol) error stop 'opnorm batch 2'
    if (abs(opnorm_tensor(x, ['row'], ['col'], ['batch'])-4.0_dp) > tol) then
        error stop 'opnorm global'
    end if

    d = delta_tensor([2,3], [character(len=5) :: 'i','batch'], by=[2])
    if (any(d%shape /= [2,2,3])) error stop 'delta by shape'
    if (count(abs(d%data-cmplx(1.0_dp,0.0_dp,dp)) < tol) /= 6) then
        error stop 'delta by count'
    end if

    x = tensor([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp], &
               [2,3], [character(len=5) :: 'i','batch'])
    dg = diag_tensor(x, by=[2])
    if (any(dg%shape /= [2,2,3])) error stop 'diag by shape'
    if (abs(sum(real(dg%data,dp))-21.0_dp) > tol) error stop 'diag by values'

    g = tensor([2.0_dp,0.5_dp,0.5_dp,1.0_dp], [2,2], ['i','j'])
    v = tensor([3.0_dp,5.0_dp], [2], ['i'])
    up = drag_tensor(v, g, ['i'])
    m1 = reshape(g%data, [2,2])
    gi = inv_tensor(g, ['i'])
    m2 = reshape(gi%data, [2,2])
    if (maxval(abs(up%data-matmul(m2,v%data))) > tol) then
        error stop 'drag nontrivial raise'
    end if
    down = drag_tensor(up, g, ['^i'])
    if (maxval(abs(down%data-v%data)) > tol) error stop 'drag nontrivial lower'

    print '(a)', 'test_misc: PASS'
end program test_misc
