program test_utils
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use e1071
    implicit none

    real(dp), allocatable :: d(:)
    real(dp), allocatable :: p(:)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: dense(:, :)
    real(dp), allocatable :: dense2(:, :)
    real(dp), allocatable :: labels(:)
    real(dp) :: cost(4, 4)
    integer, allocatable :: path(:)
    type(shortest_paths_result) :: paths
    type(matrix_csr) :: csr
    type(matrix_csr) :: csr2
    type(stft_result) :: spectrum
    type(probplot_result) :: pp
    logical :: has_y
    integer :: i

    d = ddiscrete([1.0_dp, 2.0_dp, 3.0_dp], [0.25_dp, 0.75_dp], [1.0_dp, 2.0_dp])
    if (maxval(abs(d - [0.25_dp, 0.75_dp, 0.0_dp])) > 1.0e-12_dp) error stop "ddiscrete failed"
    p = pdiscrete([0.5_dp, 1.0_dp, 2.0_dp], [0.25_dp, 0.75_dp], [1.0_dp, 2.0_dp])
    if (maxval(abs(p - [0.0_dp, 0.25_dp, 1.0_dp])) > 1.0e-12_dp) error stop "pdiscrete failed"
    q = qdiscrete([0.0_dp, 0.25_dp, 0.26_dp, 1.0_dp], [0.25_dp, 0.75_dp], [1.0_dp, 2.0_dp])
    if (maxval(abs(q - [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp])) > 1.0e-12_dp) error stop "qdiscrete failed"

    if (abs(e1071_moment([1.0_dp, 2.0_dp, 3.0_dp], 1) - 2.0_dp) > 1.0e-12_dp) error stop "moment failed"
    if (size(bincombinations(3), 1) /= 8) error stop "bincombinations failed"
    if (size(permutations(4), 1) /= 24) error stop "permutations failed"
    if (array_linear_index([2, 3, 4], [2, 2, 3]) /= 16) error stop "element linear index failed"
    if (hamming_distance_vector([1, 0, 1], [1, 1, 0]) /= 2) error stop "hamming distance failed"

    cost = reshape([0.0_dp, 1.0_dp, 9.0_dp, 9.0_dp, &
                    1.0_dp, 0.0_dp, 2.0_dp, 9.0_dp, &
                    9.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, &
                    9.0_dp, 9.0_dp, 1.0_dp, 0.0_dp], [4, 4])
    call all_shortest_paths(cost, paths)
    if (abs(paths%distance(1, 4) - 4.0_dp) > 1.0e-12_dp) error stop "Floyd-Warshall failed"
    path = extract_path(paths, 1, 4)
    if (any(path /= [1, 2, 3, 4])) error stop "shortest path extraction failed"

    dense = reshape([1.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], [2, 3])
    call dense_to_csr(dense, csr)
    call csr_to_dense(csr, dense2)
    if (maxval(abs(dense - dense2)) > 0.0_dp) error stop "CSR dense roundtrip failed"
    call write_matrix_csr(csr, "/tmp/e1071-csr.dat", [1.0_dp, -1.0_dp])
    call read_matrix_csr("/tmp/e1071-csr.dat", csr2, labels, has_y)
    call csr_to_dense(csr2, dense2)
    if (.not. has_y .or. maxval(abs(labels - [1.0_dp, -1.0_dp])) > 0.0_dp) error stop "CSR labels failed"
    if (maxval(abs(dense - dense2)) > 0.0_dp) error stop "CSR file roundtrip failed"

    call stft([(sin(2.0_dp * e1071_pi * real(i, dp) / 8.0_dp), i = 1, 32)], &
              spectrum, win=16, increment=8, coef=8)
    if (size(spectrum%values, 2) /= 8 .or. .not. all(ieee_is_finite(spectrum%values))) error stop "STFT failed"

    call probplot_normal([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], pp)
    if (.not. ieee_is_finite(pp%slope) .or. pp%slope <= 0.0_dp) error stop "probplot failed"

    print '(a)', "test_utils: PASS"
end program test_utils
