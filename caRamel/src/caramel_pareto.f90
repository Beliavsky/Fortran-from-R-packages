module caramel_pareto
    use caramel_kinds, only: dp
    implicit none
    private
    public :: pareto, pareto_2d, pareto_3d, dominate, dominated

contains

    subroutine pareto(x, front)
        real(dp), intent(in) :: x(:,:)
        integer, intent(out) :: front(size(x,1))
        integer, allocatable :: ord(:), pf(:)
        integer :: n, d, i, j, k, npf, w, idx, pj
        logical :: is_dominated, dominates

        n = size(x,1)
        d = size(x,2)
        front = 0
        if (n == 0) return
        allocate(ord(n), pf(n))
        call lex_order_desc(x, ord)
        npf = 0

        do i = 1, n
            idx = ord(i)
            is_dominated = .false.
            do k = 1, npf
                pj = pf(k)
                if (all(x(pj,:) >= x(idx,:))) then
                    if (any(x(pj,:) > x(idx,:)) .or. maxval(abs(x(pj,:) - x(idx,:))) <= 0.0_dp) then
                        is_dominated = .true.
                        exit
                    end if
                end if
            end do
            if (is_dominated) cycle

            w = 0
            do k = 1, npf
                pj = pf(k)
                dominates = all(x(idx,:) >= x(pj,:)) .and. any(x(idx,:) > x(pj,:))
                if (.not. dominates) then
                    w = w + 1
                    pf(w) = pj
                end if
            end do
            npf = w + 1
            pf(npf) = idx
        end do

        do j = 1, npf
            front(pf(j)) = 1
        end do
    end subroutine pareto

    subroutine pareto_2d(x, front)
        real(dp), intent(in) :: x(:,:)
        integer, intent(out) :: front(size(x,1))
        if (size(x,2) /= 2) error stop "pareto_2d: second dimension must be 2"
        call pareto(x, front)
    end subroutine pareto_2d

    subroutine pareto_3d(x, front)
        real(dp), intent(in) :: x(:,:)
        integer, intent(out) :: front(size(x,1))
        if (size(x,2) /= 3) error stop "pareto_3d: second dimension must be 3"
        call pareto(x, front)
    end subroutine pareto_3d

    subroutine dominate(matobj, rank)
        real(dp), intent(in) :: matobj(:,:)
        integer, intent(out) :: rank(size(matobj,1))
        logical, allocatable :: left(:)
        integer, allocatable :: idx(:), ft(:)
        integer :: n, i, m, level

        n = size(matobj,1)
        rank = 0
        if (n == 0) return
        allocate(left(n))
        left = .true.
        level = 1
        do while (any(left))
            m = count(left)
            allocate(idx(m), ft(m))
            idx = pack([(i, i=1,n)], left)
            call pareto(matobj(idx,:), ft)
            do i = 1, m
                if (ft(i) == 1) rank(idx(i)) = level
            end do
            left = rank == 0
            deallocate(idx, ft)
            level = level + 1
        end do
    end subroutine dominate

    subroutine dominated(x, y, is_dominated)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:,:)
        logical, intent(out) :: is_dominated(size(y,1))
        integer :: i

        if (size(y,2) /= size(x)) error stop "dominated: inconsistent dimensions"
        do i = 1, size(y,1)
            is_dominated(i) = all(x >= y(i,:)) .and. any(x > y(i,:))
        end do
    end subroutine dominated

    subroutine lex_order_desc(x, ord)
        real(dp), intent(in) :: x(:,:)
        integer, intent(out) :: ord(size(x,1))
        integer :: i, j, key

        ord = [(i, i=1,size(x,1))]
        do i = 2, size(ord)
            key = ord(i)
            j = i - 1
            do while (j >= 1)
                if (.not. lex_before(x, key, ord(j))) exit
                ord(j+1) = ord(j)
                j = j - 1
            end do
            ord(j+1) = key
        end do
    end subroutine lex_order_desc

    logical function lex_before(x, ia, ib) result(before)
        real(dp), intent(in) :: x(:,:)
        integer, intent(in) :: ia, ib
        integer :: k

        before = .false.
        do k = 1, size(x,2)
            if (x(ia,k) > x(ib,k)) then
                before = .true.
                return
            else if (x(ia,k) < x(ib,k)) then
                return
            end if
        end do
        before = ia < ib
    end function lex_before

end module caramel_pareto
