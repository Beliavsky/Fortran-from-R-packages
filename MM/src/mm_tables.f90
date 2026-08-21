! Computational translation of the R package MM 1.7-0.
! Upstream license: GPL-2. This translation is GPL-2.
module mm_tables
    use mm_types, only : mb_type, gunter_type, gunter_mb_type
    use partitions, only : compositions
    implicit none
    private

    public :: gunter, gunter_mb

contains

    function gunter(obs) result(out)
        integer, intent(in) :: obs(:,:)
        type(gunter_type) :: out
        integer, allocatable :: comp(:,:)
        integer :: nr, k, y_total, s, r

        nr = size(obs, 1)
        k = size(obs, 2)
        if (nr < 1 .or. k < 1) error stop "gunter: empty observation matrix"
        if (any(obs < 0)) error stop "gunter: negative observations"
        if (.not. all(sum(obs, dim=2) == sum(obs(1, :)))) error stop "gunter: row sums differ"
        y_total = sum(obs(1, :))
        if (y_total == 0) then
            allocate(out%tbl(1, k), out%d(1))
            out%tbl = 0
            out%d(1) = nr
            return
        end if
        comp = compositions(y_total, k)
        allocate(out%tbl(size(comp, 2), k), out%d(size(comp, 2)))
        out%tbl = transpose(comp)
        out%d = 0
        do r = 1, nr
            do s = 1, size(out%tbl, 1)
                if (all(obs(r, :) == out%tbl(s, :))) then
                    out%d(s) = out%d(s) + 1
                    exit
                end if
            end do
        end do
    end function gunter

    function gunter_mb(obs) result(out)
        type(mb_type), intent(in) :: obs
        type(gunter_mb_type) :: out
        integer :: k, ns, s, j, r, code, radix

        k = size(obs%m)
        ns = 1
        do j = 1, k
            if (obs%m(j) > huge(ns) / ns) error stop "gunter_mb: support too large"
            ns = ns * (obs%m(j) + 1)
        end do
        allocate(out%tbl(ns, k), out%d(ns), out%m(k))
        out%m = obs%m
        out%d = 0

        do s = 1, ns
            code = s - 1
            do j = 1, k
                radix = obs%m(j) + 1
                out%tbl(s, j) = mod(code, radix)
                code = code / radix
            end do
        end do
        do r = 1, size(obs%counts, 1)
            do s = 1, ns
                if (all(obs%counts(r, :) == out%tbl(s, :))) then
                    out%d(s) = out%d(s) + 1
                    exit
                end if
            end do
        end do
    end function gunter_mb

end module mm_tables
