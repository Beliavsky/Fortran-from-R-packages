program test_binary_cache_compat
    use genalg, only : dp, rbga_bin_control, rbga_bin_result, rbga_bin
    implicit none
    type(rbga_bin_control) :: legacy_ctl, fixed_ctl
    type(rbga_bin_result) :: legacy, fixed
    integer :: i, mismatch_legacy, mismatch_fixed

    legacy_ctl%pop_size = 80
    legacy_ctl%iters = 4
    legacy_ctl%elitism = 8
    legacy_ctl%mutation_chance = 1.0_dp
    legacy_ctl%zero_to_one_ratio = 1.0_dp
    legacy_ctl%seed = 818181
    legacy_ctl%legacy_binary_eval_cache = .true.
    fixed_ctl = legacy_ctl
    fixed_ctl%legacy_binary_eval_cache = .false.

    call rbga_bin(8, objective, legacy, legacy_ctl)
    call rbga_bin(8, objective, fixed, fixed_ctl)

    mismatch_legacy = 0
    mismatch_fixed = 0
    do i = 1, legacy_ctl%pop_size
        if (abs(legacy%evaluations(i) - objective(legacy%population(i,:))) > 1.0e-12_dp) then
            mismatch_legacy = mismatch_legacy + 1
        end if
        if (abs(fixed%evaluations(i) - objective(fixed%population(i,:))) > 1.0e-12_dp) then
            mismatch_fixed = mismatch_fixed + 1
        end if
    end do
    if (mismatch_legacy == 0) error stop "legacy cache behavior was not reproduced"
    if (mismatch_fixed /= 0) error stop "corrected cache mode has stale evaluations"
    if (legacy%nfe >= fixed%nfe) error stop "legacy cache should evaluate fewer chromosomes"
    print *, "test_binary_cache_compat: PASS", mismatch_legacy, legacy%nfe, fixed%nfe

contains

    function objective(x) result(f)
        integer, intent(in) :: x(:)
        real(dp) :: f
        integer :: j
        f = 0.0_dp
        do j = 1, size(x)
            f = f + real(j*x(j), dp)
        end do
    end function objective

end program test_binary_cache_compat
