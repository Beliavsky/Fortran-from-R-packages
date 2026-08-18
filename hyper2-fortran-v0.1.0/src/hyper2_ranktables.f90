! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_ranktables
    use hyper2_kinds, only : dp, name_len
    use hyper2_types, only : hyper2_model, hyper3_model, hyper2_create, hyper3_create, operator(+)
    use hyper2_models, only : rankvec_likelihood, ordervec2supp3
    implicit none
    private

    public :: ranktable_to_ordertable, ordertable_to_ranktable
    public :: ordertable2supp, ordertable2supp3

contains

    function ranktable_to_ordertable(rank_ids) result(order)
        integer, intent(in) :: rank_ids(:,:)
        integer, allocatable :: order(:,:)
        integer :: r, j, pos, id, n

        n = size(rank_ids,2)
        allocate(order(n,size(rank_ids,1)))
        order = 0
        do r = 1, size(rank_ids,1)
            do pos = 1, n
                id = rank_ids(r,pos)
                if (id < 1 .or. id > n) error stop "ranktable_to_ordertable: invalid player id"
                if (order(id,r) /= 0) error stop "ranktable_to_ordertable: duplicate player id"
                order(id,r) = pos
            end do
        end do
        do j = 1, size(order,2)
            if (any(order(:,j) == 0)) error stop "ranktable_to_ordertable: incomplete ranking"
        end do
    end function ranktable_to_ordertable

    function ordertable_to_ranktable(order) result(rank_ids)
        integer, intent(in) :: order(:,:)
        integer, allocatable :: rank_ids(:,:)
        integer :: n, j, player, pos

        n = size(order,1)
        allocate(rank_ids(size(order,2),n))
        rank_ids = 0
        do j = 1, size(order,2)
            do player = 1, n
                pos = order(player,j)
                if (pos < 1 .or. pos > n) error stop "ordertable_to_ranktable: invalid rank"
                if (rank_ids(j,pos) /= 0) error stop "ordertable_to_ranktable: duplicate rank"
                rank_ids(j,pos) = player
            end do
        end do
    end function ordertable_to_ranktable

    function ordertable2supp(order, names, times) result(h)
        integer, intent(in) :: order(:,:)
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in), optional :: times(:)
        type(hyper2_model) :: h, tmp
        character(len=name_len), allocatable :: v(:)
        real(dp), allocatable :: tt(:)
        integer :: j, player, pos, n

        n = size(order,1)
        if (size(names) /= n) error stop "ordertable2supp: names size mismatch"
        allocate(v(n), tt(size(order,2)))
        tt = 1.0_dp
        if (present(times)) then
            if (size(times) /= size(order,2)) error stop "ordertable2supp: times size mismatch"
            tt = times
        end if
        h = hyper2_create(names)
        do j = 1, size(order,2)
            do player = 1, n
                pos = order(player,j)
                if (pos < 1 .or. pos > n) error stop "ordertable2supp: invalid rank"
                v(pos) = names(player)
            end do
            tmp = rankvec_likelihood(v)
            call tmp%scale(tt(j))
            h = h + tmp
        end do
    end function ordertable2supp

    function ordertable2supp3(order, names, times) result(h)
        integer, intent(in) :: order(:,:)
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in), optional :: times(:)
        type(hyper3_model) :: h, tmp
        character(len=name_len), allocatable :: v(:)
        real(dp), allocatable :: tt(:)
        integer :: j, player, pos, n

        n = size(order,1)
        if (size(names) /= n) error stop "ordertable2supp3: names size mismatch"
        allocate(v(n), tt(size(order,2)))
        tt = 1.0_dp
        if (present(times)) then
            if (size(times) /= size(order,2)) error stop "ordertable2supp3: times size mismatch"
            tt = times
        end if
        h = hyper3_create(names)
        do j = 1, size(order,2)
            do player = 1, n
                pos = order(player,j)
                if (pos < 1 .or. pos > n) error stop "ordertable2supp3: invalid rank"
                v(pos) = names(player)
            end do
            tmp = ordervec2supp3(v)
            call tmp%scale(tt(j))
            h = h + tmp
        end do
    end function ordertable2supp3

end module hyper2_ranktables
