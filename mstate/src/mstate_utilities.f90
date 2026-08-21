module mstate_utilities
    use mstate_kinds, only : dp
    use mstate_types, only : transition_map, msdata_type, etmdata_type
    implicit none
    private
    public :: msdata_to_etm, etm_to_msdata, bootstrap_msdata
    public :: optimal_weights_multiple, optimal_weights_matrix

contains

    subroutine msdata_to_etm(ms, etm)
        type(msdata_type), intent(in) :: ms
        type(etmdata_type), intent(out) :: etm
        integer :: i, j, n, tostate
        logical :: same
        n = 0
        i = 1
        do while (i <= ms%n)
            n = n + 1
            j = i + 1
            do while (j <= ms%n)
                same = ms%id(j) == ms%id(i) .and. ms%from(j) == ms%from(i) .and. &
                       ms%tstart(j) == ms%tstart(i) .and. ms%tstop(j) == ms%tstop(i)
                if (.not. same) exit
                j = j + 1
            end do
            i = j
        end do
        etm%n = n
        allocate(etm%id(n), etm%from(n), etm%to(n), etm%entry(n), etm%exit(n))
        n = 0
        i = 1
        do while (i <= ms%n)
            n = n + 1
            j = i
            tostate = 99
            do while (j <= ms%n)
                same = ms%id(j) == ms%id(i) .and. ms%from(j) == ms%from(i) .and. &
                       ms%tstart(j) == ms%tstart(i) .and. ms%tstop(j) == ms%tstop(i)
                if (.not. same) exit
                if (ms%status(j) == 1) tostate = ms%to(j)
                j = j + 1
            end do
            etm%id(n) = ms%id(i)
            etm%entry(n) = ms%tstart(i)
            etm%exit(n) = ms%tstop(i)
            etm%from(n) = ms%from(i)
            etm%to(n) = tostate
            i = j
        end do
    end subroutine msdata_to_etm

    subroutine etm_to_msdata(etm, tr, ms, info)
        type(etmdata_type), intent(in) :: etm
        type(transition_map), intent(in) :: tr
        type(msdata_type), intent(out) :: ms
        integer, intent(out), optional :: info
        integer :: i, j, n, row
        if (present(info)) info = 0
        n = 0
        do i = 1, etm%n
            if (etm%from(i) < 1 .or. etm%from(i) > tr%nstate) then
                if (present(info)) info = 1
                return
            end if
            n = n + count(tr%trans(etm%from(i), :) > 0)
        end do
        ms%n = n
        allocate(ms%id(n), ms%from(n), ms%to(n), ms%trans(n), ms%status(n))
        allocate(ms%tstart(n), ms%tstop(n), ms%time(n))
        row = 0
        do i = 1, etm%n
            do j = 1, tr%nstate
                if (tr%trans(etm%from(i), j) <= 0) cycle
                row = row + 1
                ms%id(row) = etm%id(i)
                ms%from(row) = etm%from(i)
                ms%to(row) = j
                ms%trans(row) = tr%trans(etm%from(i), j)
                ms%tstart(row) = etm%entry(i)
                ms%tstop(row) = etm%exit(i)
                ms%time(row) = etm%exit(i) - etm%entry(i)
                ms%status(row) = merge(1, 0, etm%to(i) == j)
            end do
        end do
    end subroutine etm_to_msdata

    subroutine bootstrap_msdata(ms, out)
        type(msdata_type), intent(in) :: ms
        type(msdata_type), intent(out) :: out
        integer, allocatable :: uid(:), pick(:), counts(:)
        integer :: i, j, nuid, nrow, row
        real(dp) :: u
        call unique_ids(ms%id, uid)
        nuid = size(uid)
        allocate(pick(nuid), counts(nuid)); counts = 0
        do i = 1, nuid
            call random_number(u)
            pick(i) = uid(min(nuid, 1 + int(u * real(nuid, dp))))
            counts(i) = count(ms%id == pick(i))
        end do
        nrow = sum(counts)
        out%n = nrow
        allocate(out%id(nrow), out%from(nrow), out%to(nrow), out%trans(nrow), out%status(nrow))
        allocate(out%tstart(nrow), out%tstop(nrow), out%time(nrow))
        row = 0
        do i = 1, nuid
            do j = 1, ms%n
                if (ms%id(j) /= pick(i)) cycle
                row = row + 1
                out%id(row) = i
                out%from(row) = ms%from(j)
                out%to(row) = ms%to(j)
                out%trans(row) = ms%trans(j)
                out%status(row) = ms%status(j)
                out%tstart(row) = ms%tstart(j)
                out%tstop(row) = ms%tstop(j)
                out%time(row) = ms%time(j)
            end do
        end do
    end subroutine bootstrap_msdata

    subroutine optimal_weights_multiple(ms, tr, grid, transition, weights, min_time, info)
        type(msdata_type), intent(in) :: ms
        type(transition_map), intent(in) :: tr
        real(dp), intent(in) :: grid(:)
        integer, intent(in) :: transition
        real(dp), allocatable, intent(out) :: weights(:)
        real(dp), intent(in), optional :: min_time
        integer, intent(out), optional :: info
        real(dp), allocatable :: wm(:, :)
        integer :: ierr
        ierr = 0
        call compute_weight_matrix(ms, tr, grid, transition, wm, .false., min_time, ierr)
        if (present(info)) info = ierr
        if (ierr /= 0) then
            allocate(weights(0))
            return
        end if
        allocate(weights(size(grid)))
        weights = maxval(wm, dim=2)
    end subroutine optimal_weights_multiple

    subroutine optimal_weights_matrix(ms, tr, grid, transition, weights, min_time, info)
        type(msdata_type), intent(in) :: ms
        type(transition_map), intent(in) :: tr
        real(dp), intent(in) :: grid(:)
        integer, intent(in) :: transition
        real(dp), allocatable, intent(out) :: weights(:, :)
        real(dp), intent(in), optional :: min_time
        integer, intent(out), optional :: info
        call compute_weight_matrix(ms, tr, grid, transition, weights, .true., min_time, info)
    end subroutine optimal_weights_matrix

    subroutine compute_weight_matrix(ms, tr, grid, transition, weights, root_form, min_time, info)
        type(msdata_type), intent(in) :: ms
        type(transition_map), intent(in) :: tr
        real(dp), intent(in) :: grid(:)
        integer, intent(in) :: transition
        real(dp), allocatable, intent(out) :: weights(:, :)
        logical, intent(in) :: root_form
        real(dp), intent(in), optional :: min_time
        integer, intent(out), optional :: info
        type(etmdata_type) :: etm
        real(dp), allocatable :: nums(:, :), subevent(:)
        integer :: g, s, i, froms, tos
        real(dp) :: val, total, mt, dt, prev_grid
        if (present(info)) info = 0
        if (transition < 1 .or. transition > tr%ntrans) then
            if (present(info)) info = 1
            allocate(weights(0, 0))
            return
        end if
        froms = tr%from(transition); tos = tr%to(transition)
        call msdata_to_etm(ms, etm)
        allocate(nums(tr%nstate, size(grid)), subevent(size(grid)))
        nums = 0.0_dp; subevent = 0.0_dp
        do g = 1, size(grid)
            do i = 1, etm%n
                if (etm%entry(i) <= grid(g) .and. etm%exit(i) > grid(g)) then
                    nums(etm%from(i), g) = nums(etm%from(i), g) + 1.0_dp
                end if
                if (etm%from(i) == froms .and. etm%to(i) == tos .and. etm%exit(i) > grid(g)) then
                    subevent(g) = subevent(g) + 1.0_dp
                end if
            end do
        end do
        allocate(weights(size(grid), tr%nstate)); weights = 0.0_dp
        mt = 0.0_dp
        if (present(min_time)) mt = min_time
        prev_grid = mt
        do g = 1, size(grid)
            total = sum(nums(:, g))
            dt = grid(g) - prev_grid
            prev_grid = grid(g)
            if (total <= 0.0_dp) cycle
            do s = 1, tr%nstate
                val = subevent(g) * nums(s, g) * (total - nums(s, g)) / (total * total)
                if (root_form) val = sqrt(max(0.0_dp, val))
                weights(g, s) = val * dt
            end do
        end do
    end subroutine compute_weight_matrix

    subroutine unique_ids(ids, u)
        integer, intent(in) :: ids(:)
        integer, allocatable, intent(out) :: u(:)
        integer, allocatable :: tmp(:)
        integer :: i, n
        allocate(tmp(size(ids))); n = 0
        do i = 1, size(ids)
            if (n == 0) then
                n = 1; tmp(1) = ids(i)
            else if (.not. any(tmp(1:n) == ids(i))) then
                n = n + 1; tmp(n) = ids(i)
            end if
        end do
        allocate(u(n)); if (n > 0) u = tmp(1:n)
    end subroutine unique_ids

end module mstate_utilities
