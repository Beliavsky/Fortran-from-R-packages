! Computational translation of the R package MM 1.7-0.
! Upstream license: GPL-2. This translation is GPL-2.
module mm_sampling
    use mm_kinds, only : dp
    use mm_types, only : paras_type
    use mm_parameters, only : paras_dimension, p
    use mm_core, only : mm_single_log
    implicit none
    private

    public :: rmm

contains

    subroutine rmm(n, y_total, par, out, burnin, every, start)
        integer, intent(in) :: n, y_total
        type(paras_type), intent(in) :: par
        integer, allocatable, intent(out) :: out(:,:)
        integer, intent(in), optional :: burnin, every
        integer, intent(in), optional :: start(:)
        integer, allocatable :: state(:)
        integer :: k, nburn, nevery, i, j

        if (n < 0 .or. y_total < 0) error stop "rmm: n and y_total must be nonnegative"
        k = paras_dimension(par)
        allocate(state(k), out(n, k))
        nburn = 4 * y_total
        nevery = 4 * y_total
        if (present(burnin)) nburn = burnin
        if (present(every)) nevery = every
        if (nburn < 0 .or. nevery < 0) error stop "rmm: burnin/every must be nonnegative"

        if (present(start)) then
            if (size(start) /= k .or. any(start < 0) .or. sum(start) /= y_total) then
                error stop "rmm: invalid start"
            end if
            state = start
        else
            call multinomial_start(y_total, p(par), state)
        end if

        do i = 1, nburn
            call mh_update(state, par)
        end do
        do i = 1, n
            do j = 1, nevery
                call mh_update(state, par)
            end do
            out(i, :) = state
        end do
    end subroutine rmm

    subroutine multinomial_start(y_total, prob, state)
        integer, intent(in) :: y_total
        real(dp), intent(in) :: prob(:)
        integer, intent(out) :: state(:)
        real(dp) :: u, cs
        integer :: t, j

        state = 0
        do t = 1, y_total
            call random_number(u)
            cs = 0.0_dp
            do j = 1, size(prob)
                cs = cs + prob(j)
                if (u <= cs .or. j == size(prob)) then
                    state(j) = state(j) + 1
                    exit
                end if
            end do
        end do
    end subroutine multinomial_start

    subroutine mh_update(state, par)
        integer, intent(inout) :: state(:)
        type(paras_type), intent(in) :: par
        integer, allocatable :: proposed(:)
        integer :: i, j, k
        real(dp) :: u, lognum, logden, alpha

        k = size(state)
        if (k < 2) return
        call random_number(u)
        i = 1 + int(u * real(k, dp))
        i = min(i, k)
        do
            call random_number(u)
            j = 1 + int(u * real(k, dp))
            j = min(j, k)
            if (j /= i) exit
        end do
        proposed = state
        proposed(i) = proposed(i) - 1
        proposed(j) = proposed(j) + 1
        if (proposed(i) < 0) return

        lognum = mm_single_log(proposed, par)
        logden = mm_single_log(state, par)
        alpha = min(1.0_dp, exp(min(0.0_dp, lognum - logden)))
        call random_number(u)
        if (u < alpha) state = proposed
    end subroutine mh_update

end module mm_sampling
