program test_tables_sampling
    use multiplicative_multinomial
    implicit none
    type(gunter_type) :: g
    type(mb_type) :: mb
    type(gunter_mb_type) :: gm
    type(glm_fit_type) :: mbfit
    type(paras_type) :: par
    integer, allocatable :: draws(:,:)
    integer :: obs(5,3), bobs(8,2), i

    obs(1,:) = [1, 0, 1]
    obs(2,:) = [1, 0, 1]
    obs(3,:) = [0, 2, 0]
    obs(4,:) = [2, 0, 0]
    obs(5,:) = [1, 1, 0]
    g = gunter(obs)
    if (sum(g%d) /= 5) error stop "gunter count total"
    call check_support_count(g, [1,0,1], 2)
    call check_support_count(g, [0,2,0], 1)

    mb = make_mb(obs(:,1:2), [2,2])
    gm = gunter_mb(mb)
    if (size(gm%tbl,1) /= 9 .or. sum(gm%d) /= 5) error stop "gunter_mb support"

    bobs(1,:) = [0, 0]
    bobs(2,:) = [0, 1]
    bobs(3,:) = [0, 2]
    bobs(4,:) = [1, 0]
    bobs(5,:) = [1, 1]
    bobs(6,:) = [1, 2]
    bobs(7,:) = [2, 0]
    bobs(8,:) = [2, 1]
    mb = make_mb(bobs, [2,2])
    call lindsey_mb(mb, mbfit)
    if (.not. allocated(mbfit%coefficients)) error stop "Lindsey_MB no coefficients"
    if (size(mbfit%coefficients) /= 6) error stop "Lindsey_MB wrong coefficient count"

    par = paras(3)
    call rmm(100, 7, par, draws, burnin=40, every=10)
    if (size(draws,1) /= 100 .or. size(draws,2) /= 3) error stop "rMM shape"
    do i = 1, size(draws,1)
        if (sum(draws(i,:)) /= 7 .or. any(draws(i,:) < 0)) error stop "rMM invalid state"
    end do

    print '(a)', 'test_tables_sampling: PASS'

contains

    subroutine check_support_count(x, state, ref)
        type(gunter_type), intent(in) :: x
        integer, intent(in) :: state(:), ref
        integer :: s
        do s = 1, size(x%tbl,1)
            if (all(x%tbl(s,:) == state)) then
                if (x%d(s) /= ref) error stop "gunter support count mismatch"
                return
            end if
        end do
        error stop "gunter support state missing"
    end subroutine check_support_count

end program test_tables_sampling
