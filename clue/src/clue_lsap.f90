! SPDX-License-Identifier: GPL-2.0-only
module clue_lsap
    use clue_kinds, only: dp
    implicit none
    private
    public :: solve_lsap, assignment_cost
contains
    subroutine solve_lsap(cost, assignment, maximum, status)
        real(dp), intent(in) :: cost(:,:)
        integer, allocatable, intent(out) :: assignment(:)
        logical, intent(in), optional :: maximum
        integer, intent(out), optional :: status
        integer :: n, m, i, j, j0, j1, i0
        real(dp) :: cur, delta, maxc
        real(dp), allocatable :: a(:,:), u(:), v(:), minv(:)
        integer, allocatable :: p(:), way(:)
        logical, allocatable :: used(:)
        logical :: want_max

        n = size(cost,1)
        m = size(cost,2)
        if (n <= 0 .or. m <= 0 .or. n > m) then
            allocate(assignment(0))
            if (present(status)) status = 1
            return
        end if
        want_max = .false.
        if (present(maximum)) want_max = maximum
        allocate(a(n,m))
        a = cost
        if (want_max) then
            maxc = maxval(a)
            a = maxc - a
        end if
        allocate(u(0:n), v(0:m), p(0:m), way(0:m), minv(0:m), used(0:m))
        u=0.0_dp
        v=0.0_dp
        p=0
        way=0
        do i=1,n
            p(0)=i
            j0=0
            minv=huge(1.0_dp)
            used=.false.
            do
                used(j0)=.true.
                i0=p(j0)
                delta=huge(1.0_dp)
                j1=0
                do j=1,m
                    if (.not. used(j)) then
                        cur=a(i0,j)-u(i0)-v(j)
                        if (cur < minv(j)) then
                            minv(j)=cur
                            way(j)=j0
                        end if
                        if (minv(j) < delta) then
                            delta=minv(j)
                            j1=j
                        end if
                    end if
                end do
                do j=0,m
                    if (used(j)) then
                        u(p(j))=u(p(j))+delta
                        v(j)=v(j)-delta
                    else if (j>0) then
                        minv(j)=minv(j)-delta
                    end if
                end do
                j0=j1
                if (p(j0)==0) exit
            end do
            do
                j1=way(j0)
                p(j0)=p(j1)
                j0=j1
                if (j0==0) exit
            end do
        end do
        allocate(assignment(n))
        assignment=0
        do j=1,m
            if (p(j)>=1 .and. p(j)<=n) assignment(p(j))=j
        end do
        if (present(status)) status = merge(0,2,all(assignment>0))
    end subroutine solve_lsap

    function assignment_cost(cost, assignment) result(value)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: assignment(:)
        real(dp) :: value
        integer :: i
        value=0.0_dp
        do i=1,min(size(cost,1),size(assignment))
            if (assignment(i)>=1 .and. assignment(i)<=size(cost,2)) value=value+cost(i,assignment(i))
        end do
    end function assignment_cost
end module clue_lsap
