! SPDX-License-Identifier: GPL-2.0-or-later
module evir_data
    use evir_kinds, only : dp
    use evir_types, only : decluster_result, evir_ok, evir_invalid_input, evir_domain_error
    use evir_math, only : sort_descending
    implicit none
    private
    public :: find_threshold, block_maxima, block_maxima_groups, decluster

contains

    real(dp) function find_threshold(data, ne, status) result(threshold)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: ne
        integer, intent(out), optional :: status
        real(dp), allocatable :: sorted(:)
        real(dp) :: target
        integer :: i
        if (present(status)) status = evir_ok
        if (size(data) == 0 .or. ne < 1 .or. ne > size(data)) then
            threshold = 0.0_dp
            if (present(status)) status = evir_invalid_input
            return
        end if
        sorted = data
        call sort_descending(sorted)
        target = sorted(ne)
        threshold = sorted(size(sorted))
        do i = ne+1, size(sorted)
            if (sorted(i) < target) then
                threshold = sorted(i)
                return
            end if
        end do
    end function find_threshold

    function block_maxima(data, block_size) result(maxima)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: block_size
        real(dp), allocatable :: maxima(:)
        integer :: nblocks, b, lo, hi
        if (block_size <= 0 .or. size(data) == 0) then
            allocate(maxima(0))
            return
        end if
        nblocks = (size(data)+block_size-1)/block_size
        allocate(maxima(nblocks))
        do b = 1, nblocks
            lo = (b-1)*block_size+1
            hi = min(b*block_size, size(data))
            maxima(b) = maxval(data(lo:hi))
        end do
    end function block_maxima

    function block_maxima_groups(data, group) result(maxima)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: group(:)
        real(dp), allocatable :: maxima(:)
        integer, allocatable :: labels(:)
        integer :: i, j, ng
        logical :: found
        if (size(data) /= size(group) .or. size(data) == 0) then
            allocate(maxima(0))
            return
        end if
        allocate(labels(size(group)))
        ng = 0
        do i = 1, size(group)
            found = .false.
            do j = 1, ng
                if (labels(j) == group(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                ng = ng+1
                labels(ng) = group(i)
            end if
        end do
        allocate(maxima(ng))
        do j = 1, ng
            maxima(j) = -huge(1.0_dp)
            do i = 1, size(data)
                if (group(i) == labels(j)) maxima(j) = max(maxima(j), data(i))
            end do
        end do
    end function block_maxima_groups

    function decluster(series, times, run) result(out)
        real(dp), intent(in) :: series(:), times(:), run
        type(decluster_result) :: out
        integer, allocatable :: starts(:), ends(:)
        integer :: n, i, nc, j, idx
        real(dp) :: best
        n = size(series)
        if (n == 0 .or. size(times) /= n .or. run < 0.0_dp) then
            out%status = evir_invalid_input
            allocate(out%values(0), out%times(0))
            return
        end if
        if (n == 1) then
            allocate(out%values(1), out%times(1))
            out%values(1) = series(1)
            out%times(1) = times(1)
            return
        end if
        allocate(starts(n), ends(n))
        nc = 1
        starts(1) = 1
        do i = 1, n-1
            if (times(i+1)-times(i) > run) then
                ends(nc) = i
                nc = nc+1
                starts(nc) = i+1
            end if
        end do
        ends(nc) = n
        if (nc <= 2) then
            out%status = evir_domain_error
            allocate(out%values(0), out%times(0))
            return
        end if
        allocate(out%values(nc), out%times(nc))
        do j = 1, nc
            best = series(starts(j))
            idx = starts(j)
            do i = starts(j)+1, ends(j)
                if (series(i) > best) then
                    best = series(i)
                    idx = i
                end if
            end do
            out%values(j) = best
            out%times(j) = times(idx)
        end do
    end function decluster

end module evir_data
