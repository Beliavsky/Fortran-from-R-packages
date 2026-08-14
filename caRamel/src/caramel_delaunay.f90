module caramel_delaunay
    use caramel_kinds, only: dp
    use caramel_linalg, only: solve_linear
    implicit none
    private
    public :: delaunay_nd

contains

    subroutine delaunay_nd(points, simplices, ok)
        real(dp), intent(in) :: points(:,:)
        integer, allocatable, intent(out) :: simplices(:,:)
        logical, intent(out) :: ok
        real(dp), allocatable :: work(:,:), center(:)
        integer, allocatable :: current(:,:), next(:,:), facets(:,:), boundary(:,:)
        logical, allocatable :: bad(:), used(:)
        integer :: n, d, ns, i, j, k, f, t, nb, ngood, nfac, count_same
        integer :: facet_pos, ns_new, nkeep, row
        real(dp) :: radius, m, scale, delta

        n = size(points,1)
        d = size(points,2)
        ok = .false.
        if (d < 1 .or. n < d + 1) then
            allocate(simplices(0, d + 1))
            return
        end if

        allocate(work(n+d+1,d), center(d))
        center = 0.5_dp * (minval(points, dim=1) + maxval(points, dim=1))
        radius = maxval(abs(points - spread(center,1,n)))
        if (radius <= 100.0_dp * epsilon(1.0_dp)) then
            allocate(simplices(0,d+1))
            return
        end if

        scale = max(1.0_dp, maxval(abs(points)))
        do i = 1, n
            do j = 1, d
                delta = 1.0e-12_dp * scale * real(mod(37*i + 101*j + 13*i*j, 97) - 48, dp) / 97.0_dp
                work(i,j) = points(i,j) + delta
            end do
        end do

        ! A simplex containing the full bounding box. Relative to center, its
        ! feasible region is y_j >= -m and sum(y_j) <= m.
        m = max(1.0_dp, 1000.0_dp * real(d,dp) * radius)
        work(n+1,:) = center - m
        do j = 1, d
            work(n+1+j,:) = center - m
            work(n+1+j,j) = center(j) + real(d,dp) * m
        end do

        allocate(current(1,d+1))
        current(1,:) = [(n+i, i=1,d+1)]
        ns = 1

        do i = 1, n
            allocate(bad(ns))
            do j = 1, ns
                bad(j) = point_in_circumsphere(work(i,:), work(current(j,:),:))
            end do
            nb = count(bad)
            if (nb == 0) then
                deallocate(bad)
                allocate(simplices(0,d+1))
                return
            end if

            nfac = nb * (d + 1)
            allocate(facets(nfac,d), used(nfac))
            facet_pos = 0
            do j = 1, ns
                if (.not. bad(j)) cycle
                do k = 1, d + 1
                    facet_pos = facet_pos + 1
                    f = 0
                    do t = 1, d + 1
                        if (t == k) cycle
                        f = f + 1
                        facets(facet_pos,f) = current(j,t)
                    end do
                    call sort_int(facets(facet_pos,:))
                end do
            end do

            used = .false.
            nb = 0
            allocate(boundary(nfac,d))
            do f = 1, nfac
                if (used(f)) cycle
                count_same = 1
                used(f) = .true.
                do t = f + 1, nfac
                    if (.not. used(t)) then
                        if (all(facets(t,:) == facets(f,:))) then
                            used(t) = .true.
                            count_same = count_same + 1
                        end if
                    end if
                end do
                if (count_same == 1) then
                    nb = nb + 1
                    boundary(nb,:) = facets(f,:)
                end if
            end do

            ngood = count(.not. bad)
            ns_new = ngood + nb
            allocate(next(ns_new,d+1))
            row = 0
            do j = 1, ns
                if (.not. bad(j)) then
                    row = row + 1
                    next(row,:) = current(j,:)
                end if
            end do
            do j = 1, nb
                row = row + 1
                next(row,1:d) = boundary(j,:)
                next(row,d+1) = i
            end do

            call move_alloc(next, current)
            ns = ns_new
            deallocate(bad, facets, used, boundary)
        end do

        nkeep = 0
        do i = 1, ns
            if (all(current(i,:) <= n)) nkeep = nkeep + 1
        end do
        allocate(simplices(nkeep,d+1))
        row = 0
        do i = 1, ns
            if (all(current(i,:) <= n)) then
                row = row + 1
                simplices(row,:) = current(i,:)
                call sort_int(simplices(row,:))
            end if
        end do
        call unique_simplex_rows(simplices)
        ok = size(simplices,1) > 0
    end subroutine delaunay_nd

    logical function point_in_circumsphere(point, vertices) result(inside)
        real(dp), intent(in) :: point(:)
        real(dp), intent(in) :: vertices(:,:)
        real(dp), allocatable :: a(:,:), b(:), center(:)
        real(dp) :: r2, d2, tol, sc
        logical :: solved
        integer :: d, i

        d = size(point)
        allocate(a(d,d), b(d), center(d))
        do i = 1, d
            a(i,:) = 2.0_dp * (vertices(i+1,:) - vertices(1,:))
            b(i) = dot_product(vertices(i+1,:), vertices(i+1,:)) - &
                   dot_product(vertices(1,:), vertices(1,:))
        end do
        call solve_linear(a, b, center, solved)
        if (.not. solved) then
            inside = .false.
            return
        end if
        r2 = sum((vertices(1,:) - center)**2)
        d2 = sum((point - center)**2)
        sc = max(1.0_dp, r2, d2)
        tol = 1.0e3_dp * epsilon(1.0_dp) * sc
        inside = d2 <= r2 + tol
    end function point_in_circumsphere

    subroutine sort_int(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, key
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
    end subroutine sort_int

    subroutine unique_simplex_rows(a)
        integer, allocatable, intent(inout) :: a(:,:)
        integer, allocatable :: tmp(:,:)
        logical :: duplicate
        integer :: i, j, nuniq, d

        d = size(a,2)
        if (size(a,1) <= 1) return
        allocate(tmp(size(a,1),d))
        nuniq = 0
        do i = 1, size(a,1)
            duplicate = .false.
            do j = 1, nuniq
                if (all(a(i,:) == tmp(j,:))) then
                    duplicate = .true.
                    exit
                end if
            end do
            if (.not. duplicate) then
                nuniq = nuniq + 1
                tmp(nuniq,:) = a(i,:)
            end if
        end do
        deallocate(a)
        allocate(a(nuniq,d))
        a = tmp(1:nuniq,:)
    end subroutine unique_simplex_rows

end module caramel_delaunay
