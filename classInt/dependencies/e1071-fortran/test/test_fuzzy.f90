program test_fuzzy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use e1071
    implicit none

    real(dp) :: x(8, 2)
    real(dp) :: initial(2, 2)
    real(dp) :: shell_initial(2, 2)
    real(dp) :: radii(2)
    type(fuzzy_cluster_model) :: cfit
    type(fuzzy_cluster_model) :: sfit
    type(fclust_indices_result) :: index_result
    type(rng_state) :: rng
    integer :: i

    x = reshape([0.0_dp, 0.0_dp, 0.2_dp, 0.1_dp, -0.1_dp, 0.2_dp, 0.1_dp, -0.2_dp, &
                 5.0_dp, 5.0_dp, 5.2_dp, 5.1_dp, 4.9_dp, 5.2_dp, 5.1_dp, 4.8_dp], [8, 2], order=[2, 1])
    initial(1, :) = [0.0_dp, 0.0_dp]
    initial(2, :) = [5.0_dp, 5.0_dp]
    call cmeans_fit(x, initial, cfit, m=2.0_dp, iter_max=100)
    if (cfit%iterations /= 3) error stop "cmeans upstream iteration fixture failed"
    if (abs(cfit%within_error - 0.0343428119711938_dp) > 1.0e-14_dp) error stop "cmeans upstream objective fixture failed"
    if (maxval(abs(cfit%centers - reshape([0.050014989247376726_dp, 5.0500097344945685_dp, &
                                           0.02500397838521997_dp, 5.0250087500686869_dp], [2, 2]))) > 1.0e-13_dp) then
        error stop "cmeans upstream center fixture failed"
    end if
    if (maxval(abs(sum(cfit%membership, dim=2) - 1.0_dp)) > 1.0e-10_dp) error stop "cmeans membership normalization failed"
    if (sqrt(sum(cfit%centers(1, :)**2)) > 0.5_dp) error stop "cmeans first center failed"
    if (sqrt(sum((cfit%centers(2, :) - 5.0_dp)**2)) > 0.5_dp) error stop "cmeans second center failed"
    if (.not. ieee_is_finite(cfit%within_error)) error stop "cmeans objective failed"

    call rng_seed(rng, 1234)
    call cmeans_fit_k(x, 2, rng, cfit)
    if (maxval(abs(sum(cfit%membership, dim=2) - 1.0_dp)) > 1.0e-10_dp) error stop "cmeans scalar-k API failed"

    initial(1, :) = [0.0_dp, 0.0_dp]
    initial(2, :) = [5.0_dp, 5.0_dp]
    call cmeans_fit(x, initial, cfit, m=2.0_dp, distance=fuzzy_manhattan, iter_max=100)
    if (cfit%iterations /= 1) error stop "Manhattan cmeans upstream iteration fixture failed"
    if (abs(cfit%within_error - 0.21844433535388025_dp) > 1.0e-14_dp) then
        error stop "Manhattan cmeans upstream objective fixture failed"
    end if
    if (maxval(abs(cfit%centers - initial)) > 1.0e-14_dp) then
        error stop "Manhattan cmeans upstream weighted-median center fixture failed"
    end if

    index_result = fclust_indices(cfit, x)
    if (.not. ieee_is_finite(index_result%partition_coefficient)) error stop "fclustIndex failed"
    if (index_result%partition_coefficient <= 0.0_dp) error stop "fclustIndex partition coefficient failed"

    do i = 1, 8
        x(i, 1) = cos(2.0_dp * e1071_pi * real(i - 1, dp) / 8.0_dp)
        x(i, 2) = sin(2.0_dp * e1071_pi * real(i - 1, dp) / 8.0_dp)
    end do
    shell_initial(1, :) = [0.0_dp, 0.0_dp]
    shell_initial(2, :) = [0.2_dp, 0.0_dp]
    radii = [0.9_dp, 1.1_dp]
    call cshell_fit(x, shell_initial, sfit, radius=radii, m=2.0_dp, iter_max=8)
    if (maxval(abs(sum(sfit%membership, dim=2) - 1.0_dp)) > 1.0e-8_dp) error stop "cshell normalization failed"
    if (.not. all(ieee_is_finite(sfit%radius))) error stop "cshell radii failed"

    print '(a)', "test_fuzzy: PASS"
end program test_fuzzy
