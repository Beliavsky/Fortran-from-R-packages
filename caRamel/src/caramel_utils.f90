module caramel_utils
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use caramel_kinds, only: dp
    use caramel_linalg, only: determinant
    use caramel_pareto, only: dominate, dominated
    use caramel_random, only: random_uniform, random_index
    implicit none
    private
    public :: val2rank, boxes, rselect, matvcov, vol_splx, dimprove
    public :: downsize, decrease_pop, column_mean, column_sd

contains

    subroutine val2rank(x, opt, rank)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: opt
        real(dp), intent(out) :: rank(size(x))
        integer, allocatable :: ord(:)
        integer :: n, i, j, key, first, last, unique_rank
        real(dp) :: value

        n = size(x)
        if (n == 0) return
        allocate(ord(n))
        ord = [(i, i=1,n)]
        do i = 2, n
            key = ord(i)
            value = x(key)
            j = i - 1
            do while (j >= 1)
                if (x(ord(j)) <= value) exit
                ord(j+1) = ord(j)
                j = j - 1
            end do
            ord(j+1) = key
        end do

        unique_rank = 0
        first = 1
        do while (first <= n)
            last = first
            do while (last < n)
                if (abs(x(ord(last+1)) - x(ord(first))) > 0.0_dp) exit
                last = last + 1
            end do
            unique_rank = unique_rank + 1
            select case (opt)
            case (1)
                value = 0.5_dp * real(first + last, dp)
            case (2)
                value = real(unique_rank, dp)
            case default
                value = real(last, dp)
            end select
            do i = first, last
                rank(ord(i)) = value
            end do
            first = last + 1
        end do
    end subroutine val2rank

    subroutine boxes(points, prec, box_number)
        real(dp), intent(in) :: points(:,:), prec(:)
        integer(int64), intent(out) :: box_number(size(points,1))
        real(dp), allocatable :: ranks(:), levels(:,:)
        real(dp) :: vmin
        integer(int64), allocatable :: base(:)
        integer :: n, d, i, j

        n = size(points,1)
        d = size(points,2)
        if (size(prec) /= d) error stop "boxes: inconsistent precision dimension"
        if (any(prec <= 0.0_dp)) error stop "boxes: precision must be positive"
        allocate(levels(n,d), ranks(n), base(d))
        do j = 1, d
            vmin = minval(points(:,j))
            levels(:,j) = floor((points(:,j) - vmin) / prec(j))
            call val2rank(levels(:,j), 2, ranks)
            levels(:,j) = ranks
        end do
        base(1) = 1_int64
        do j = 2, d
            base(j) = base(j-1) * (1_int64 + int(maxval(levels(:,j-1)), int64))
        end do
        box_number = 0_int64
        do i = 1, n
            do j = 1, d
                box_number(i) = box_number(i) + int(levels(i,j),int64) * base(j)
            end do
        end do
    end subroutine boxes

    subroutine rselect(n, facp, ix)
        integer, intent(in) :: n
        real(dp), intent(in) :: facp(:)
        integer, allocatable, intent(out) :: ix(:)
        real(dp), allocatable :: fc(:)
        real(dp) :: total, f
        integer :: i, j

        total = sum(facp)
        if (n <= 0 .or. total <= 0.0_dp) then
            allocate(ix(0))
            return
        end if
        if (any(facp < 0.0_dp)) error stop "rselect: weights must be nonnegative"
        allocate(ix(n), fc(size(facp)))
        fc(1) = facp(1) / total
        do i = 2, size(facp)
            fc(i) = fc(i-1) + facp(i) / total
        end do
        do i = 1, n
            f = random_uniform()
            ix(i) = size(facp)
            do j = 1, size(facp)
                if (fc(j) > f) then
                    ix(i) = j
                    exit
                end if
            end do
        end do
    end subroutine rselect

    function column_mean(x) result(mean_x)
        real(dp), intent(in) :: x(:,:)
        real(dp) :: mean_x(size(x,2))
        if (size(x,1) == 0) then
            mean_x = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            mean_x = sum(x,dim=1) / real(size(x,1),dp)
        end if
    end function column_mean

    function column_sd(x) result(sd_x)
        real(dp), intent(in) :: x(:,:)
        real(dp) :: sd_x(size(x,2))
        real(dp) :: m
        integer :: j, n

        n = size(x,1)
        if (n <= 1) then
            sd_x = 0.0_dp
            return
        end if
        do j = 1, size(x,2)
            m = sum(x(:,j)) / real(n,dp)
            sd_x(j) = sqrt(sum((x(:,j)-m)**2) / real(n-1,dp))
        end do
    end function column_sd

    subroutine matvcov(x, g, rr)
        real(dp), intent(in) :: x(:,:), g(:)
        real(dp), intent(out) :: rr(size(x,2),size(x,2))
        real(dp), allocatable :: xred(:,:), sd(:)
        real(dp) :: nanv
        integer :: n, p, i, j

        n = size(x,1)
        p = size(x,2)
        if (size(g) /= p) error stop "matvcov: inconsistent center dimension"
        nanv = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(xred(n,p), sd(p))
        sd = column_sd(x)
        do j = 1, p
            if (sd(j) > 0.0_dp) then
                xred(:,j) = (x(:,j) - g(j)) / sd(j)
            else
                xred(:,j) = nanv
            end if
        end do
        rr = 0.0_dp
        do i = 1, p
            rr(i,i) = 1.0_dp
            do j = i + 1, p
                rr(i,j) = sum(xred(:,i) * xred(:,j)) / real(max(1,n),dp)
                rr(j,i) = rr(i,j)
            end do
        end do
    end subroutine matvcov

    real(dp) function vol_splx(s) result(volume)
        real(dp), intent(in) :: s(:,:)
        real(dp), allocatable :: a(:,:)
        integer :: d, i

        d = size(s,2)
        if (size(s,1) /= d + 1) error stop "vol_splx: simplex must have d+1 vertices"
        allocate(a(d,d))
        do i = 1, d
            a(i,:) = s(i+1,:) - s(1,:)
        end do
        volume = abs(determinant(a)) / factorial_real(d)
    end function vol_splx

    subroutine dimprove(o_splx, f_splx, oriented_edges, edge_length)
        real(dp), intent(in) :: o_splx(:,:)
        integer, intent(in) :: f_splx(:)
        integer, allocatable, intent(out) :: oriented_edges(:,:)
        real(dp), allocatable, intent(out) :: edge_length(:)
        logical, allocatable :: isdom(:)
        integer, allocatable :: tmp_edges(:,:)
        real(dp), allocatable :: tmp_length(:)
        integer :: n, d, i, j, m

        n = size(o_splx,1)
        d = size(o_splx,2)
        if (size(f_splx) /= n) error stop "dimprove: inconsistent dimensions"
        allocate(tmp_edges(n*n,2), tmp_length(n*n), isdom(n))
        m = 0
        do i = 1, n
            if (f_splx(i) /= 1) cycle
            call dominated(o_splx(i,:), o_splx, isdom)
            do j = 1, n
                if (isdom(j)) then
                    m = m + 1
                    tmp_edges(m,:) = [j,i]
                    tmp_length(m) = sqrt(sum((o_splx(i,:) - o_splx(j,:))**2))
                end if
            end do
        end do
        allocate(oriented_edges(m,2), edge_length(m))
        if (m > 0) then
            oriented_edges = tmp_edges(1:m,:)
            edge_length = tmp_length(1:m)
        end if
    end subroutine dimprove

    subroutine downsize(points, front_rank, prec, indices)
        real(dp), intent(in) :: points(:,:), prec(:)
        integer, intent(in) :: front_rank(:)
        integer, allocatable, intent(out) :: indices(:)
        integer(int64), allocatable :: box_id(:), unique_box(:)
        integer, allocatable :: tmp(:), candidates(:), maxima(:)
        integer :: n, d, i, j, k, nb, m, nc, fmin, chosen
        logical :: found

        n = size(points,1)
        d = size(points,2)
        if (size(front_rank) /= n) error stop "downsize: inconsistent front ranks"
        allocate(box_id(n), maxima(d), tmp(n), unique_box(n))
        call boxes(points, prec, box_id)
        do j = 1, d
            maxima(j) = maxloc(points(:,j),dim=1)
        end do

        nb = 0
        do i = 1, n
            found = .false.
            do j = 1, nb
                if (unique_box(j) == box_id(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                nb = nb + 1
                unique_box(nb) = box_id(i)
            end if
        end do
        call sort_int64(unique_box(1:nb))

        m = 0
        do k = 1, nb
            nc = count(box_id == unique_box(k))
            allocate(candidates(nc))
            candidates = pack([(i,i=1,n)], box_id == unique_box(k))
            fmin = minval(front_rank(candidates))
            candidates = pack(candidates, front_rank(candidates) == fmin)
            chosen = 0
            if (fmin == 1) then
                outer: do i = 1, size(candidates)
                    do j = 1, d
                        if (candidates(i) == maxima(j)) then
                            chosen = candidates(i)
                            exit outer
                        end if
                    end do
                end do outer
            end if
            if (chosen == 0) chosen = candidates(random_index(size(candidates)))
            m = m + 1
            tmp(m) = chosen
            deallocate(candidates)
        end do
        allocate(indices(m))
        indices = tmp(1:m)
    end subroutine downsize

    subroutine decrease_pop(matobj, minmax, prec, archsize, popsize, ind_arch, ind_pop)
        real(dp), intent(in) :: matobj(:,:), prec(:)
        logical, intent(in) :: minmax(:)
        integer, intent(in) :: archsize, popsize
        integer, allocatable, intent(out) :: ind_arch(:), ind_pop(:)
        real(dp), allocatable :: work(:,:), arch(:,:), arch_prec(:)
        integer, allocatable :: rank(:), indices(:), original(:), ord(:), rest_rank(:)
        integer, allocatable :: candidates(:), selected(:), tmp(:)
        integer :: n, d, j, na, np, fmax, nlow, need, iter

        n = size(matobj,1)
        d = size(matobj,2)
        if (size(minmax) /= d .or. size(prec) /= d) error stop "decrease_pop: inconsistent dimensions"
        allocate(work(n,d), rank(n), original(n))
        work = matobj
        do j = 1, d
            if (.not. minmax(j)) work(:,j) = -work(:,j)
        end do
        call dominate(work, rank)
        call downsize(work, rank, prec, indices)
        work = work(indices,:)
        original = 0
        original(1:size(indices)) = indices
        original = original(1:size(indices))

        deallocate(rank)
        allocate(rank(size(indices)), ord(size(indices)))
        call dominate(work, rank)
        call order_int(rank, ord)
        work = work(ord,:)
        rank = rank(ord)
        original = original(ord)

        na = count(rank == 1)
        np = size(rank) - na
        allocate(ind_arch(na), ind_pop(np))
        if (na > 0) ind_arch = pack(original, rank == 1)
        if (np > 0) ind_pop = pack(original, rank > 1)

        if (np > popsize) then
            rest_rank = pack(rank, rank > 1)
            fmax = rest_rank(popsize)
            nlow = count(rest_rank < fmax)
            allocate(tmp(popsize))
            if (nlow > 0) tmp(1:nlow) = pack(ind_pop, rest_rank < fmax)
            need = popsize - nlow
            if (need > 0) then
                candidates = pack(ind_pop, rest_rank == fmax)
                call choose_without_replacement(candidates, need, selected)
                tmp(nlow+1:popsize) = selected
            end if
            deallocate(ind_pop)
            allocate(ind_pop(popsize))
            ind_pop = tmp
        end if

        if (allocated(rank)) deallocate(rank)
        if (na > archsize) then
            allocate(arch(na,d), arch_prec(d))
            arch = matobj(ind_arch,:)
            do j = 1, d
                if (.not. minmax(j)) arch(:,j) = -arch(:,j)
            end do
            arch_prec = 2.0_dp * prec
            do iter = 1, 60
                allocate(rank(na))
                rank = 1
                call downsize(arch, rank, arch_prec, indices)
                deallocate(rank)
                if (size(indices) <= archsize) exit
                arch_prec = 2.0_dp * arch_prec
                deallocate(indices)
            end do
            if (size(indices) > archsize) indices = indices(1:archsize)
            tmp = ind_arch(indices)
            deallocate(ind_arch)
            allocate(ind_arch(size(tmp)))
            ind_arch = tmp
        end if
    end subroutine decrease_pop

    real(dp) function factorial_real(n) result(f)
        integer, intent(in) :: n
        integer :: i
        f = 1.0_dp
        do i = 2, n
            f = f * real(i,dp)
        end do
    end function factorial_real

    subroutine order_int(x, ord)
        integer, intent(in) :: x(:)
        integer, intent(out) :: ord(size(x))
        integer :: i, j, key
        ord = [(i,i=1,size(x))]
        do i = 2, size(x)
            key = ord(i)
            j = i - 1
            do while (j >= 1)
                if (x(ord(j)) <= x(key)) exit
                ord(j+1) = ord(j)
                j = j - 1
            end do
            ord(j+1) = key
        end do
    end subroutine order_int

    subroutine sort_int64(x)
        integer(int64), intent(inout) :: x(:)
        integer(int64) :: key
        integer :: i, j
        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j+1) = x(j)
                j = j - 1
            end do
            x(j+1) = key
        end do
    end subroutine sort_int64

    subroutine choose_without_replacement(candidates, n, selected)
        integer, intent(in) :: candidates(:), n
        integer, allocatable, intent(out) :: selected(:)
        integer, allocatable :: pool(:)
        integer :: i, j, tmp

        if (n < 0 .or. n > size(candidates)) error stop "choose_without_replacement: invalid size"
        pool = candidates
        allocate(selected(n))
        do i = 1, n
            j = i - 1 + random_index(size(pool) - i + 1)
            tmp = pool(i)
            pool(i) = pool(j)
            pool(j) = tmp
            selected(i) = pool(i)
        end do
    end subroutine choose_without_replacement

end module caramel_utils
