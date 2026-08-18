! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_types
    use hyper2_kinds, only : dp, name_len
    implicit none
    private

    real(dp), parameter :: zero_tol = 0.0_dp

    type, public :: hyper2_term
        integer, allocatable :: ids(:)
        real(dp) :: power = 0.0_dp
    end type hyper2_term

    type, public :: hyper3_term
        integer, allocatable :: ids(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: power = 0.0_dp
    end type hyper3_term

    type, public :: hyper2_model
        character(len=name_len), allocatable :: pnames(:)
        type(hyper2_term), allocatable :: terms(:)
    contains
        procedure :: size => h2_size
        procedure :: nterms => h2_nterms
        procedure :: add_player => h2_add_player
        procedure :: add_term => h2_add_term_names
        procedure :: add_term_ids => h2_add_term_ids
        procedure :: set_term => h2_set_term_names
        procedure :: get_power => h2_get_power_names
        procedure :: scale => h2_scale
        procedure :: tidy => h2_tidy
    end type hyper2_model

    type, public :: hyper3_model
        character(len=name_len), allocatable :: pnames(:)
        type(hyper3_term), allocatable :: terms(:)
    contains
        procedure :: size => h3_size
        procedure :: nterms => h3_nterms
        procedure :: add_player => h3_add_player
        procedure :: add_term => h3_add_term_names
        procedure :: add_term_ids => h3_add_term_ids
        procedure :: set_term => h3_set_term_names
        procedure :: get_power => h3_get_power_names
        procedure :: scale => h3_scale
    end type hyper3_model

    interface operator(+)
        module procedure h2_plus_h2
        module procedure h3_plus_h3
    end interface
    interface operator(-)
        module procedure h2_minus_h2
        module procedure h3_minus_h3
    end interface
    interface operator(*)
        module procedure h2_times_real
        module procedure real_times_h2
        module procedure h3_times_real
        module procedure real_times_h3
    end interface
    interface operator(==)
        module procedure h2_equal
        module procedure h3_equal
    end interface

    public :: hyper2_create, hyper3_create, player_index, names_to_ids
    public :: operator(+), operator(-), operator(*), operator(==)

contains

    function hyper2_create(pnames) result(h)
        character(len=*), intent(in), optional :: pnames(:)
        type(hyper2_model) :: h
        if (present(pnames)) then
            allocate(h%pnames(size(pnames)))
            h%pnames = pnames
        else
            allocate(h%pnames(0))
        end if
        allocate(h%terms(0))
    end function hyper2_create

    function hyper3_create(pnames) result(h)
        character(len=*), intent(in), optional :: pnames(:)
        type(hyper3_model) :: h
        if (present(pnames)) then
            allocate(h%pnames(size(pnames)))
            h%pnames = pnames
        else
            allocate(h%pnames(0))
        end if
        allocate(h%terms(0))
    end function hyper3_create

    integer function h2_size(self)
        class(hyper2_model), intent(in) :: self
        h2_size = size(self%pnames)
    end function h2_size

    integer function h3_size(self)
        class(hyper3_model), intent(in) :: self
        h3_size = size(self%pnames)
    end function h3_size

    integer function h2_nterms(self)
        class(hyper2_model), intent(in) :: self
        h2_nterms = size(self%terms)
    end function h2_nterms

    integer function h3_nterms(self)
        class(hyper3_model), intent(in) :: self
        h3_nterms = size(self%terms)
    end function h3_nterms

    integer function player_index(pnames, name)
        character(len=*), intent(in) :: pnames(:)
        character(len=*), intent(in) :: name
        integer :: i
        player_index = 0
        do i = 1, size(pnames)
            if (trim(pnames(i)) == trim(name)) then
                player_index = i
                return
            end if
        end do
    end function player_index

    subroutine h2_add_player(self, name)
        class(hyper2_model), intent(inout) :: self
        character(len=*), intent(in) :: name
        character(len=name_len), allocatable :: tmp(:)
        integer :: n
        if (player_index(self%pnames, name) > 0) return
        n = size(self%pnames)
        allocate(tmp(n + 1))
        if (n > 0) tmp(1:n) = self%pnames
        tmp(n + 1) = name
        call move_alloc(tmp, self%pnames)
    end subroutine h2_add_player

    subroutine h3_add_player(self, name)
        class(hyper3_model), intent(inout) :: self
        character(len=*), intent(in) :: name
        character(len=name_len), allocatable :: tmp(:)
        integer :: n
        if (player_index(self%pnames, name) > 0) return
        n = size(self%pnames)
        allocate(tmp(n + 1))
        if (n > 0) tmp(1:n) = self%pnames
        tmp(n + 1) = name
        call move_alloc(tmp, self%pnames)
    end subroutine h3_add_player

    function names_to_ids(pnames, names, add_missing, h2, h3) result(ids)
        character(len=*), intent(in) :: pnames(:)
        character(len=*), intent(in) :: names(:)
        logical, intent(in), optional :: add_missing
        class(hyper2_model), intent(inout), optional :: h2
        class(hyper3_model), intent(inout), optional :: h3
        integer, allocatable :: ids(:)
        integer :: i, id
        logical :: add
        add = .false.
        if (present(add_missing)) add = add_missing
        allocate(ids(size(names)))
        do i = 1, size(names)
            id = player_index(pnames, names(i))
            if (id == 0 .and. add) then
                if (present(h2)) then
                    call h2%add_player(names(i))
                    id = player_index(h2%pnames, names(i))
                else if (present(h3)) then
                    call h3%add_player(names(i))
                    id = player_index(h3%pnames, names(i))
                end if
            end if
            ids(i) = id
        end do
    end function names_to_ids

    subroutine sort_unique_int(x, ok)
        integer, allocatable, intent(inout) :: x(:)
        logical, intent(out) :: ok
        integer :: i, j, t
        ok = .true.
        do i = 2, size(x)
            t = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= t) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = t
        end do
        if (size(x) > 1) then
            do i = 2, size(x)
                if (x(i) == x(i - 1)) then
                    ok = .false.
                    return
                end if
            end do
        end if
    end subroutine sort_unique_int

    logical function same_ints(a, b)
        integer, intent(in) :: a(:), b(:)
        if (size(a) /= size(b)) then
            same_ints = .false.
        else if (size(a) == 0) then
            same_ints = .true.
        else
            same_ints = all(a == b)
        end if
    end function same_ints

    logical function same_real_bits(a, b)
        use iso_fortran_env, only : int64
        real(dp), intent(in) :: a, b
        same_real_bits = transfer(a, 0_int64) == transfer(b, 0_int64)
    end function same_real_bits

    logical function same_h3_key(a_ids, a_w, b_ids, b_w)
        integer, intent(in) :: a_ids(:), b_ids(:)
        real(dp), intent(in) :: a_w(:), b_w(:)
        integer :: i
        same_h3_key = same_ints(a_ids, b_ids) .and. size(a_w) == size(b_w)
        if (.not. same_h3_key) return
        do i = 1, size(a_w)
            if (.not. same_real_bits(a_w(i), b_w(i))) then
                same_h3_key = .false.
                return
            end if
        end do
    end function same_h3_key

    subroutine remove_h2_term(self, k)
        class(hyper2_model), intent(inout) :: self
        integer, intent(in) :: k
        type(hyper2_term), allocatable :: tmp(:)
        integer :: n
        n = size(self%terms)
        allocate(tmp(n - 1))
        if (k > 1) tmp(1:k - 1) = self%terms(1:k - 1)
        if (k < n) tmp(k:n - 1) = self%terms(k + 1:n)
        call move_alloc(tmp, self%terms)
    end subroutine remove_h2_term

    subroutine remove_h3_term(self, k)
        class(hyper3_model), intent(inout) :: self
        integer, intent(in) :: k
        type(hyper3_term), allocatable :: tmp(:)
        integer :: n
        n = size(self%terms)
        allocate(tmp(n - 1))
        if (k > 1) tmp(1:k - 1) = self%terms(1:k - 1)
        if (k < n) tmp(k:n - 1) = self%terms(k + 1:n)
        call move_alloc(tmp, self%terms)
    end subroutine remove_h3_term

    subroutine h2_add_term_names(self, names, power)
        class(hyper2_model), intent(inout) :: self
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: power
        integer, allocatable :: ids(:)
        integer :: i
        do i = 1, size(names)
            call self%add_player(names(i))
        end do
        allocate(ids(size(names)))
        do i = 1, size(names)
            ids(i) = player_index(self%pnames, names(i))
        end do
        call self%add_term_ids(ids, power)
    end subroutine h2_add_term_names

    subroutine h2_add_term_ids(self, ids_in, power)
        class(hyper2_model), intent(inout) :: self
        integer, intent(in) :: ids_in(:)
        real(dp), intent(in) :: power
        integer, allocatable :: ids(:)
        type(hyper2_term), allocatable :: tmp(:)
        integer :: i, n
        logical :: ok
        if (abs(power) <= zero_tol) return
        ids = ids_in
        if (size(ids) == 0 .or. any(ids < 1) .or. any(ids > size(self%pnames))) error stop "invalid hyper2 bracket"
        call sort_unique_int(ids, ok)
        if (.not. ok) error stop "hyper2 bracket contains repeated player"
        do i = 1, size(self%terms)
            if (same_ints(self%terms(i)%ids, ids)) then
                self%terms(i)%power = self%terms(i)%power + power
                if (abs(self%terms(i)%power) <= zero_tol) call remove_h2_term(self, i)
                return
            end if
        end do
        n = size(self%terms)
        allocate(tmp(n + 1))
        if (n > 0) tmp(1:n) = self%terms
        tmp(n + 1)%ids = ids
        tmp(n + 1)%power = power
        call move_alloc(tmp, self%terms)
    end subroutine h2_add_term_ids

    subroutine h2_set_term_names(self, names, power)
        class(hyper2_model), intent(inout) :: self
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: power
        integer, allocatable :: ids(:)
        integer :: i, j
        logical :: ok
        do i = 1, size(names)
            call self%add_player(names(i))
        end do
        allocate(ids(size(names)))
        do i = 1, size(names)
            ids(i) = player_index(self%pnames, names(i))
        end do
        call sort_unique_int(ids, ok)
        if (.not. ok) error stop "hyper2 bracket contains repeated player"
        do j = 1, size(self%terms)
            if (same_ints(self%terms(j)%ids, ids)) then
                if (abs(power) <= zero_tol) then
                    call remove_h2_term(self, j)
                else
                    self%terms(j)%power = power
                end if
                return
            end if
        end do
        if (abs(power) > zero_tol) call self%add_term_ids(ids, power)
    end subroutine h2_set_term_names

    real(dp) function h2_get_power_names(self, names)
        class(hyper2_model), intent(in) :: self
        character(len=*), intent(in) :: names(:)
        integer, allocatable :: ids(:)
        integer :: i
        logical :: ok
        allocate(ids(size(names)))
        do i = 1, size(names)
            ids(i) = player_index(self%pnames, names(i))
            if (ids(i) == 0) then
                h2_get_power_names = 0.0_dp
                return
            end if
        end do
        call sort_unique_int(ids, ok)
        if (.not. ok) then
            h2_get_power_names = 0.0_dp
            return
        end if
        h2_get_power_names = 0.0_dp
        do i = 1, size(self%terms)
            if (same_ints(self%terms(i)%ids, ids)) then
                h2_get_power_names = self%terms(i)%power
                return
            end if
        end do
    end function h2_get_power_names

    subroutine h2_scale(self, a)
        class(hyper2_model), intent(inout) :: self
        real(dp), intent(in) :: a
        integer :: i
        if (abs(a) <= zero_tol) then
            deallocate(self%terms)
            allocate(self%terms(0))
            return
        end if
        do i = 1, size(self%terms)
            self%terms(i)%power = a * self%terms(i)%power
        end do
    end subroutine h2_scale

    subroutine h2_tidy(self)
        class(hyper2_model), intent(inout) :: self
        ! Terms are canonicalized as they are inserted; pnames are intentionally retained.
        if (self%nterms() < 0) error stop "unreachable"
    end subroutine h2_tidy

    subroutine canonical_h3(ids, weights, ok)
        integer, allocatable, intent(inout) :: ids(:)
        real(dp), allocatable, intent(inout) :: weights(:)
        logical, intent(out) :: ok
        integer :: i, j, ti, m, k
        real(dp) :: tw
        integer, allocatable :: ni(:)
        real(dp), allocatable :: nw(:)
        ok = size(ids) == size(weights)
        if (.not. ok) return
        do i = 2, size(ids)
            ti = ids(i)
            tw = weights(i)
            j = i - 1
            do while (j >= 1)
                if (ids(j) <= ti) exit
                ids(j + 1) = ids(j)
                weights(j + 1) = weights(j)
                j = j - 1
            end do
            ids(j + 1) = ti
            weights(j + 1) = tw
        end do
        m = 0
        do i = 1, size(ids)
            if (weights(i) < 0.0_dp) then
                ok = .false.
                return
            end if
            if (weights(i) > 0.0_dp) then
                if (m > 0) then
                    if (ids(i) == ids(m)) then
                        weights(m) = weights(m) + weights(i)
                    else
                        m = m + 1
                        ids(m) = ids(i)
                        weights(m) = weights(i)
                    end if
                else
                    m = 1
                    ids(m) = ids(i)
                    weights(m) = weights(i)
                end if
            end if
        end do
        allocate(ni(m), nw(m))
        do k = 1, m
            ni(k) = ids(k)
            nw(k) = weights(k)
        end do
        call move_alloc(ni, ids)
        call move_alloc(nw, weights)
    end subroutine canonical_h3

    subroutine h3_add_term_names(self, names, weights_in, power)
        class(hyper3_model), intent(inout) :: self
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: weights_in(:)
        real(dp), intent(in) :: power
        integer, allocatable :: ids(:)
        integer :: i
        if (size(names) /= size(weights_in)) error stop "hyper3 names/weights size mismatch"
        do i = 1, size(names)
            call self%add_player(names(i))
        end do
        allocate(ids(size(names)))
        do i = 1, size(names)
            ids(i) = player_index(self%pnames, names(i))
        end do
        call self%add_term_ids(ids, weights_in, power)
    end subroutine h3_add_term_names

    subroutine h3_add_term_ids(self, ids_in, weights_in, power)
        class(hyper3_model), intent(inout) :: self
        integer, intent(in) :: ids_in(:)
        real(dp), intent(in) :: weights_in(:)
        real(dp), intent(in) :: power
        integer, allocatable :: ids(:)
        real(dp), allocatable :: w(:)
        type(hyper3_term), allocatable :: tmp(:)
        logical :: ok
        integer :: i, n
        if (abs(power) <= zero_tol) return
        ids = ids_in
        w = weights_in
        if (any(ids < 1) .or. any(ids > size(self%pnames))) error stop "invalid hyper3 player id"
        call canonical_h3(ids, w, ok)
        if (.not. ok .or. size(ids) == 0) error stop "invalid hyper3 weighted bracket"
        do i = 1, size(self%terms)
            if (same_h3_key(self%terms(i)%ids, self%terms(i)%weights, ids, w)) then
                self%terms(i)%power = self%terms(i)%power + power
                if (abs(self%terms(i)%power) <= zero_tol) call remove_h3_term(self, i)
                return
            end if
        end do
        n = size(self%terms)
        allocate(tmp(n + 1))
        if (n > 0) tmp(1:n) = self%terms
        tmp(n + 1)%ids = ids
        tmp(n + 1)%weights = w
        tmp(n + 1)%power = power
        call move_alloc(tmp, self%terms)
    end subroutine h3_add_term_ids

    subroutine h3_set_term_names(self, names, weights_in, power)
        class(hyper3_model), intent(inout) :: self
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: weights_in(:)
        real(dp), intent(in) :: power
        integer, allocatable :: ids(:)
        real(dp), allocatable :: w(:)
        integer :: i, j
        logical :: ok
        if (size(names) /= size(weights_in)) error stop "hyper3 names/weights size mismatch"
        do i = 1, size(names)
            call self%add_player(names(i))
        end do
        allocate(ids(size(names)), w(size(weights_in)))
        do i = 1, size(names)
            ids(i) = player_index(self%pnames, names(i))
        end do
        w = weights_in
        call canonical_h3(ids, w, ok)
        if (.not. ok) error stop "invalid hyper3 weighted bracket"
        do j = 1, size(self%terms)
            if (same_h3_key(self%terms(j)%ids, self%terms(j)%weights, ids, w)) then
                if (abs(power) <= zero_tol) then
                    call remove_h3_term(self, j)
                else
                    self%terms(j)%power = power
                end if
                return
            end if
        end do
        if (abs(power) > zero_tol) call self%add_term_ids(ids, w, power)
    end subroutine h3_set_term_names

    real(dp) function h3_get_power_names(self, names, weights_in)
        class(hyper3_model), intent(in) :: self
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: weights_in(:)
        integer, allocatable :: ids(:)
        real(dp), allocatable :: w(:)
        integer :: i
        logical :: ok
        if (size(names) /= size(weights_in)) then
            h3_get_power_names = 0.0_dp
            return
        end if
        allocate(ids(size(names)), w(size(weights_in)))
        do i = 1, size(names)
            ids(i) = player_index(self%pnames, names(i))
            if (ids(i) == 0) then
                h3_get_power_names = 0.0_dp
                return
            end if
        end do
        w = weights_in
        call canonical_h3(ids, w, ok)
        h3_get_power_names = 0.0_dp
        if (.not. ok) return
        do i = 1, size(self%terms)
            if (same_h3_key(self%terms(i)%ids, self%terms(i)%weights, ids, w)) then
                h3_get_power_names = self%terms(i)%power
                return
            end if
        end do
    end function h3_get_power_names

    subroutine h3_scale(self, a)
        class(hyper3_model), intent(inout) :: self
        real(dp), intent(in) :: a
        integer :: i
        if (abs(a) <= zero_tol) then
            deallocate(self%terms)
            allocate(self%terms(0))
            return
        end if
        do i = 1, size(self%terms)
            self%terms(i)%power = a * self%terms(i)%power
        end do
    end subroutine h3_scale

    function remap_h2(src, names) result(out)
        type(hyper2_model), intent(in) :: src
        character(len=name_len), intent(in) :: names(:)
        type(hyper2_model) :: out
        character(len=name_len), allocatable :: bn(:)
        integer :: i, j
        out = hyper2_create(names)
        do i = 1, size(src%terms)
            allocate(bn(size(src%terms(i)%ids)))
            do j = 1, size(bn)
                bn(j) = src%pnames(src%terms(i)%ids(j))
            end do
            call out%add_term(bn, src%terms(i)%power)
            deallocate(bn)
        end do
    end function remap_h2

    function remap_h3(src, names) result(out)
        type(hyper3_model), intent(in) :: src
        character(len=name_len), intent(in) :: names(:)
        type(hyper3_model) :: out
        character(len=name_len), allocatable :: bn(:)
        integer :: i, j
        out = hyper3_create(names)
        do i = 1, size(src%terms)
            allocate(bn(size(src%terms(i)%ids)))
            do j = 1, size(bn)
                bn(j) = src%pnames(src%terms(i)%ids(j))
            end do
            call out%add_term(bn, src%terms(i)%weights, src%terms(i)%power)
            deallocate(bn)
        end do
    end function remap_h3

    function union_names(a, b) result(names)
        character(len=name_len), intent(in) :: a(:), b(:)
        character(len=name_len), allocatable :: names(:), tmp(:)
        integer :: i, n
        names = a
        do i = 1, size(b)
            if (player_index(names, b(i)) == 0) then
                n = size(names)
                allocate(tmp(n + 1))
                if (n > 0) tmp(1:n) = names
                tmp(n + 1) = b(i)
                call move_alloc(tmp, names)
            end if
        end do
    end function union_names

    function h2_plus_h2(a, b) result(c)
        type(hyper2_model), intent(in) :: a, b
        type(hyper2_model) :: c, aa, bb
        character(len=name_len), allocatable :: names(:)
        integer :: i
        names = union_names(a%pnames, b%pnames)
        aa = remap_h2(a, names)
        bb = remap_h2(b, names)
        c = aa
        do i = 1, size(bb%terms)
            call c%add_term_ids(bb%terms(i)%ids, bb%terms(i)%power)
        end do
    end function h2_plus_h2

    function h2_minus_h2(a, b) result(c)
        type(hyper2_model), intent(in) :: a, b
        type(hyper2_model) :: c, bm
        bm = b
        call bm%scale(-1.0_dp)
        c = a + bm
    end function h2_minus_h2

    function h2_times_real(a, x) result(c)
        type(hyper2_model), intent(in) :: a
        real(dp), intent(in) :: x
        type(hyper2_model) :: c
        c = a
        call c%scale(x)
    end function h2_times_real

    function real_times_h2(x, a) result(c)
        real(dp), intent(in) :: x
        type(hyper2_model), intent(in) :: a
        type(hyper2_model) :: c
        c = a * x
    end function real_times_h2

    function h3_plus_h3(a, b) result(c)
        type(hyper3_model), intent(in) :: a, b
        type(hyper3_model) :: c, aa, bb
        character(len=name_len), allocatable :: names(:)
        integer :: i
        names = union_names(a%pnames, b%pnames)
        aa = remap_h3(a, names)
        bb = remap_h3(b, names)
        c = aa
        do i = 1, size(bb%terms)
            call c%add_term_ids(bb%terms(i)%ids, bb%terms(i)%weights, bb%terms(i)%power)
        end do
    end function h3_plus_h3

    function h3_minus_h3(a, b) result(c)
        type(hyper3_model), intent(in) :: a, b
        type(hyper3_model) :: c, bm
        bm = b
        call bm%scale(-1.0_dp)
        c = a + bm
    end function h3_minus_h3

    function h3_times_real(a, x) result(c)
        type(hyper3_model), intent(in) :: a
        real(dp), intent(in) :: x
        type(hyper3_model) :: c
        c = a
        call c%scale(x)
    end function h3_times_real

    function real_times_h3(x, a) result(c)
        real(dp), intent(in) :: x
        type(hyper3_model), intent(in) :: a
        type(hyper3_model) :: c
        c = a * x
    end function real_times_h3

    logical function h2_equal(a, b)
        type(hyper2_model), intent(in) :: a, b
        type(hyper2_model) :: bb
        integer :: i, j
        if (size(a%terms) /= size(b%terms)) then
            h2_equal = .false.
            return
        end if
        bb = remap_h2(b, union_names(a%pnames, b%pnames))
        h2_equal = .true.
        do i = 1, size(a%terms)
            block
                character(len=name_len), allocatable :: bn(:)
                integer, allocatable :: ids(:)
                logical :: found
                allocate(bn(size(a%terms(i)%ids)), ids(size(a%terms(i)%ids)))
                do j = 1, size(bn)
                    bn(j) = a%pnames(a%terms(i)%ids(j))
                    ids(j) = player_index(bb%pnames, bn(j))
                end do
                call sort_unique_int(ids, found)
                found = .false.
                do j = 1, size(bb%terms)
                    if (same_ints(ids, bb%terms(j)%ids)) then
                        if (same_real_bits(a%terms(i)%power, bb%terms(j)%power)) found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    h2_equal = .false.
                    return
                end if
            end block
        end do
    end function h2_equal

    logical function h3_equal(a, b)
        type(hyper3_model), intent(in) :: a, b
        type(hyper3_model) :: aa, bb
        character(len=name_len), allocatable :: names(:)
        integer :: i, j
        if (size(a%terms) /= size(b%terms)) then
            h3_equal = .false.
            return
        end if
        names = union_names(a%pnames, b%pnames)
        aa = remap_h3(a, names)
        bb = remap_h3(b, names)
        h3_equal = .true.
        do i = 1, size(aa%terms)
            block
                logical :: found
                found = .false.
                do j = 1, size(bb%terms)
                    if (same_h3_key(aa%terms(i)%ids, aa%terms(i)%weights, &
                                    bb%terms(j)%ids, bb%terms(j)%weights)) then
                        if (same_real_bits(aa%terms(i)%power, bb%terms(j)%power)) found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    h3_equal = .false.
                    return
                end if
            end block
        end do
    end function h3_equal

end module hyper2_types
