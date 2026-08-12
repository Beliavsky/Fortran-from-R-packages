! SPDX-License-Identifier: LGPL-2.0-only
module lpsolve_special
    use lpsolve_types
    use lpsolve_core, only: solve_lp
    implicit none
    private

    public :: lp_assign, lp_transport, make_q8, q8_triplets, lp_status_message

    type, public :: q8_triplets
        integer, allocatable :: constraint(:)
        integer, allocatable :: variable(:)
        real(dp), allocatable :: value(:)
    end type q8_triplets

contains

    subroutine lp_assign(cost, result, direction, control, assignment)
        real(dp), intent(in) :: cost(:,:)
        type(lp_result), intent(out) :: result
        integer, intent(in), optional :: direction
        type(lp_control), intent(in), optional :: control
        real(dp), intent(out), optional :: assignment(:,:)

        integer :: nr, nc, n, m, i, j, k, dir
        real(dp), allocatable :: c(:), a(:,:), rhs(:)
        integer, allocatable :: sense(:), ints(:)

        nr = size(cost,1)
        nc = size(cost,2)
        n = nr * nc
        m = nr + nc
        dir = LP_MIN
        if (present(direction)) dir = direction
        allocate(c(n), a(m,n), rhs(m), sense(m), ints(n))
        a = 0.0_dp
        rhs = 1.0_dp
        sense = LP_EQ
        ints = [(i, i=1,n)]

        k = 0
        do i = 1, nr
            do j = 1, nc
                k = k + 1
                c(k) = cost(i,j)
                a(i,k) = 1.0_dp
                a(nr+j,k) = 1.0_dp
            end do
        end do

        if (present(control)) then
            call solve_lp(dir, c, a, sense, rhs, result, control, integer_variables=ints)
        else
            call solve_lp(dir, c, a, sense, rhs, result, integer_variables=ints)
        end if
        if (present(assignment)) then
            if (size(assignment,1) == nr .and. size(assignment,2) == nc .and. &
                allocated(result%solution)) then
                k = 0
                do i = 1, nr
                    do j = 1, nc
                        k = k + 1
                        assignment(i,j) = result%solution(k)
                    end do
                end do
            end if
        end if
    end subroutine lp_assign


    subroutine lp_transport(cost, row_sense, row_rhs, col_sense, col_rhs, result, &
        direction, integer_variables, control, flow)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: row_sense(:), col_sense(:)
        real(dp), intent(in) :: row_rhs(:), col_rhs(:)
        type(lp_result), intent(out) :: result
        integer, intent(in), optional :: direction
        integer, intent(in), optional :: integer_variables(:)
        type(lp_control), intent(in), optional :: control
        real(dp), intent(out), optional :: flow(:,:)

        integer :: nr, nc, n, m, i, j, k, dir
        real(dp), allocatable :: c(:), a(:,:), rhs(:)
        integer, allocatable :: sense(:), ints(:)

        nr = size(cost,1)
        nc = size(cost,2)
        n = nr * nc
        m = nr + nc
        dir = LP_MIN
        if (present(direction)) dir = direction
        if (size(row_sense) /= nr .or. size(row_rhs) /= nr .or. &
            size(col_sense) /= nc .or. size(col_rhs) /= nc) then
            result%status = LP_NUMFAILURE
            return
        end if

        allocate(c(n), a(m,n), rhs(m), sense(m))
        a = 0.0_dp
        rhs(1:nr) = row_rhs
        rhs(nr+1:m) = col_rhs
        sense(1:nr) = row_sense
        sense(nr+1:m) = col_sense
        k = 0
        do i = 1, nr
            do j = 1, nc
                k = k + 1
                c(k) = cost(i,j)
                a(i,k) = 1.0_dp
                a(nr+j,k) = 1.0_dp
            end do
        end do

        if (present(integer_variables)) then
            if (present(control)) then
                call solve_lp(dir, c, a, sense, rhs, result, control, &
                    integer_variables=integer_variables)
            else
                call solve_lp(dir, c, a, sense, rhs, result, &
                    integer_variables=integer_variables)
            end if
        else
            allocate(ints(n))
            ints = [(i, i=1,n)]
            if (present(control)) then
                call solve_lp(dir, c, a, sense, rhs, result, control, integer_variables=ints)
            else
                call solve_lp(dir, c, a, sense, rhs, result, integer_variables=ints)
            end if
        end if
        if (present(flow)) then
            if (size(flow,1) == nr .and. size(flow,2) == nc .and. allocated(result%solution)) then
                k = 0
                do i = 1, nr
                    do j = 1, nc
                        k = k + 1
                        flow(i,j) = result%solution(k)
                    end do
                end do
            end if
        end if
    end subroutine lp_transport


    function lp_status_message(status) result(message)
        integer, intent(in) :: status
        character(len=32) :: message
        select case (status)
        case (LP_OPTIMAL)
            message = 'optimal'
        case (LP_SUBOPTIMAL)
            message = 'suboptimal'
        case (LP_INFEASIBLE)
            message = 'infeasible'
        case (LP_UNBOUNDED)
            message = 'unbounded'
        case (LP_NUMFAILURE)
            message = 'numerical failure'
        case (LP_TIMEOUT)
            message = 'timeout'
        case default
            message = 'unknown status'
        end select
    end function lp_status_message


    subroutine make_q8(q8)
        type(q8_triplets), intent(out) :: q8
        integer :: board(8,8), i, j, k, c, nent, r, cc, s
        integer, allocatable :: tr(:), tv(:)
        real(dp), allocatable :: tx(:)

        k = 0
        do i = 1, 8
            do j = 1, 8
                k = k + 1
                board(i,j) = k
            end do
        end do

        allocate(tr(512), tv(512), tx(512))
        nent = 0
        c = 0

        do i = 1, 8
            c = c + 1
            do j = 1, 8
                nent = nent + 1
                tr(nent) = c
                tv(nent) = board(i,j)
                tx(nent) = 1.0_dp
            end do
        end do
        do j = 1, 8
            c = c + 1
            do i = 1, 8
                nent = nent + 1
                tr(nent) = c
                tv(nent) = board(i,j)
                tx(nent) = 1.0_dp
            end do
        end do

        ! Diagonals with i+j = 3,...,15 (R indexing).
        do s = 3, 15
            c = c + 1
            do i = 1, 8
                j = s - i
                if (j >= 1 .and. j <= 8) then
                    nent = nent + 1
                    tr(nent) = c
                    tv(nent) = board(i,j)
                    tx(nent) = 1.0_dp
                end if
            end do
        end do

        ! Other diagonal direction, retaining only diagonals of length >= 2.
        do r = -6, 6
            c = c + 1
            do i = 1, 8
                cc = i + r
                if (cc >= 1 .and. cc <= 8) then
                    nent = nent + 1
                    tr(nent) = c
                    tv(nent) = board(i,cc)
                    tx(nent) = 1.0_dp
                end if
            end do
        end do

        allocate(q8%constraint(nent), q8%variable(nent), q8%value(nent))
        q8%constraint = tr(1:nent)
        q8%variable = tv(1:nent)
        q8%value = tx(1:nent)
    end subroutine make_q8

end module lpsolve_special
