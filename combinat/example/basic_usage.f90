program basic_usage
    use combinat
    implicit none

    integer, allocatable :: combinations(:, :)
    integer, allocatable :: simplex(:, :)
    real(dp), allocatable :: lattice(:, :)
    integer :: info

    call combn_indices(5, 3, combinations, info)
    if (info /= combinat_success) error stop 'combn_indices failed'
    print '(a, i0, a, i0)', 'combn(5,3) shape: ', size(combinations, 1), ' x ', size(combinations, 2)
    print '(a, *(1x, i0))', 'first combination:', combinations(:, 1)
    print '(a, *(1x, i0))', 'last combination:', combinations(:, size(combinations, 2))

    call xsimplex(3, 4, simplex, info)
    if (info /= combinat_success) error stop 'xsimplex failed'
    print '(a, i0)', 'number of {3,4} simplex points: ', size(simplex, 2)

    call hcube([2, 3], lattice, scale=[1.0_dp, 0.5_dp], translation=[-1.0_dp, 0.0_dp], info=info)
    if (info /= combinat_success) error stop 'hcube failed'
    print '(a) ', 'scaled 2 x 3 lattice:'
    print '(2f8.3)', transpose(lattice)

    print '(a, f10.6)', 'dmnom([1,1],[0.25,0.75]) = ', dmnom([1.0_dp, 1.0_dp], [0.25_dp, 0.75_dp])
end program basic_usage
