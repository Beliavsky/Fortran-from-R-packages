program test_ampute
    use iso_fortran_env, only : int64
    use mice, only : dp, mice_ok, mice_rng_state, rng_seed, ampute_mcar, ampute_continuous, ampute_right
    implicit none
    integer :: assignments(4), patterns(2, 3), types(2), info
    real(dp) :: scores(4, 2)
    logical, allocatable :: missing(:, :)
    type(mice_rng_state) :: rng

    assignments = [1, 1, 2, 2]
    patterns(1, :) = [0, 1, 1]
    patterns(2, :) = [1, 0, 1]
    call rng_seed(rng, 1_int64)
    call ampute_mcar(assignments, patterns, 1.0_dp, rng, missing, info)
    if (info /= mice_ok) error stop "MCAR status"
    if (.not. all(missing(1:2, 1)) .or. any(missing(1:2, 2:3))) error stop "MCAR pattern 1"
    if (.not. all(missing(3:4, 2)) .or. any(missing(3:4, [1,3]))) error stop "MCAR pattern 2"
    scores(:, 1) = [-1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    scores(:, 2) = [0.0_dp, 0.0_dp, -1.0_dp, 1.0_dp]
    types = ampute_right
    call rng_seed(rng, 44_int64)
    call ampute_continuous(assignments, patterns, scores, 0.5_dp, types, rng, missing, info)
    if (info /= mice_ok .or. any(shape(missing) /= [4, 3])) error stop "continuous amputation"
    print *, "test_ampute: PASS"
end program test_ampute
