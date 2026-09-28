program test_ng_reference
    use sparsem, only: dp, csr_matrix, csr_cholesky, csr_from_dense, csr_chol
    implicit none

    integer, parameter :: n = 16
    integer, parameter :: width = 4
    integer, parameter :: expected_perm(n) = [ &
        16, 12, 15, 11, 4, 3, 8, 13, 9, 1, 5, 6, 2, 7, 10, 14]
    integer, parameter :: expected_col_count(n) = [ &
        3, 4, 4, 5, 3, 4, 5, 3, 4, 3, 5, 5, 4, 3, 2, 1]
    integer, parameter :: expected_xsuper(13) = [ &
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 17]
    integer, parameter :: expected_xlindx(13) = [ &
        1, 4, 8, 12, 17, 20, 24, 29, 32, 36, 39, 44, 49]
    integer, parameter :: expected_lindx(48) = [ &
        1, 2, 3, 2, 3, 4, 7, 3, 4, 7, 16, 4, 7, 14, 15, 16, &
        5, 6, 7, 6, 7, 13, 14, 7, 13, 14, 15, 16, 8, 9, 16, 9, &
        11, 15, 16, 10, 11, 13, 11, 12, 13, 15, 16, 12, 13, 14, 15, 16]
    integer, parameter :: expected_xlnz(17) = [ &
        1, 4, 8, 12, 17, 20, 24, 29, 32, 36, 39, 44, 49, 53, 56, 58, 59]
    integer, parameter :: expected_split(n) = [ &
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 0, 0, 0]
    real(dp), parameter :: expected_lnz(58) = [ &
        2.2360679774997898_dp, -0.44721359549995793_dp, -0.44721359549995793_dp, &
        2.1908902300206643_dp, -0.091287092917527679_dp, -0.45643546458763845_dp, &
        -0.45643546458763845_dp, 2.1889875894272826_dp, -0.47586686726680072_dp, &
        -0.019034674690672028_dp, -0.45683219257612867_dp, 2.1366369348357592_dp, &
        -0.10174461594455995_dp, -0.46802523334497576_dp, -0.46802523334497576_dp, &
        -0.10174461594455998_dp, 2.2360679774997898_dp, -0.44721359549995793_dp, &
        -0.44721359549995793_dp, 2.1908902300206643_dp, -0.091287092917527679_dp, &
        -0.45643546458763845_dp, -0.45643546458763845_dp, 2.1383683143039338_dp, &
        -0.019485261911126706_dp, -0.50940041853374118_dp, -0.022268870755573378_dp, &
        -0.0089075483022293539_dp, 2.2360679774997898_dp, -0.44721359549995793_dp, &
        -0.44721359549995793_dp, 2.1908902300206643_dp, -0.45643546458763845_dp, &
        -0.45643546458763845_dp, -0.091287092917527679_dp, 2.2360679774997898_dp, &
        -0.44721359549995793_dp, -0.44721359549995793_dp, 2.1428174599500225_dp, &
        -0.46667530888203762_dp, -0.093335061776407519_dp, -0.097224022683757835_dp, &
        -0.019444804536751566_dp, 2.1868274179916107_dp, -0.47720142897350137_dp, &
        -0.45728345628591177_dp, -0.47803134450215090_dp, -0.0041495776432478385_dp, &
        2.0868287791919502_dp, -0.20915729025966812_dp, -0.11386922013076668_dp, &
        -0.0019017527287762926_dp, 2.0150124887342860_dp, -0.23464062715623543_dp, &
        -0.027023085668735566_dp, 2.0654622351656573_dp, -0.53252751455210467_dp, &
        2.0707063930472267_dp]

    real(dp) :: dense(n, n)
    type(csr_matrix) :: a
    type(csr_cholesky) :: factor
    integer :: column
    integer :: info
    integer :: neighbor
    integer :: node
    integer :: row

    dense = 0.0_dp
    do row = 1, width
        do column = 1, width
            node = (row - 1) * width + column
            dense(node, node) = 5.0_dp
            if (column < width) then
                neighbor = node + 1
                dense(node, neighbor) = -1.0_dp
                dense(neighbor, node) = -1.0_dp
            end if
            if (row < width) then
                neighbor = node + width
                dense(node, neighbor) = -1.0_dp
                dense(neighbor, node) = -1.0_dp
            end if
        end do
    end do

    a = csr_from_dense(dense)
    factor = csr_chol(a, info, cache_kb=32, level=8, nsubmax=1000, nnzlmax=1000, tmpmax=1000)
    if (info /= 0) error stop 'reference factorization status'
    if (factor%nsub /= 48 .or. factor%nnzl /= 58 .or. factor%nsuper /= 12) error stop 'reference counts'
    if (factor%tmp_size /= 10) error stop 'reference workspace size'
    if (any(factor%perm /= expected_perm)) error stop 'reference permutation'
    if (any(factor%col_count /= expected_col_count)) error stop 'reference column counts'
    if (any(factor%xsuper(1:13) /= expected_xsuper)) error stop 'reference supernodes'
    if (any(factor%xlindx /= expected_xlindx)) error stop 'reference compressed pointers'
    if (any(factor%lindx /= expected_lindx)) error stop 'reference compressed indices'
    if (any(factor%xlnz /= expected_xlnz)) error stop 'reference numerical pointers'
    if (any(factor%split /= expected_split)) error stop 'reference cache split'
    if (maxval(abs(factor%lnz - expected_lnz)) > 5.0e-14_dp) error stop 'reference numerical factor'

    print '(a)', 'test_ng_reference: PASS'
end program test_ng_reference
