program test_state
    use iso_fortran_env, only : int64
    use zigg, only : dp, ziggurat_rng, zsetseed, zrunif, zgetstate, zsetstate
    implicit none

    type(ziggurat_rng) :: a, b
    integer(int64) :: state(4), expected_state(4)
    real(dp) :: x1(12), x2(12), y1(8), y2(8)
    real(dp), allocatable :: tmp(:)

    expected_state = [1229274879_int64, 2002354486_int64, 453738772_int64, 3452246052_int64]

    call zsetseed(12345)
    tmp = zrunif(20)
    call zgetstate(state)
    if (any(state /= expected_state)) error stop 'upstream state mismatch'

    call zsetseed(24680)
    tmp = zrunif(7)
    call zgetstate(state)
    x1 = zrunif(12)
    call zsetstate(state)
    x2 = zrunif(12)
    if (maxval(abs(x1 - x2)) > 0.0_dp) error stop 'module state restore failed'

    call a%set_seed(777)
    call b%set_seed(777)
    call a%fill_normal(y1)
    call b%fill_normal(y2)
    if (maxval(abs(y1 - y2)) > 0.0_dp) error stop 'independent generator mismatch'

    call a%set_seed(1)
    call b%set_seed(2)
    call a%fill_exponential(y1)
    call b%fill_exponential(y2)
    if (maxval(abs(y1 - y2)) < 1.0e-15_dp) error stop 'independent generators not independent'

    print '(a)', 'test_state: PASS'
end program test_state
