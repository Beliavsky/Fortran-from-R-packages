! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_sums
    use hyper2_kinds, only : dp, name_len
    use hyper2_types, only : hyper2_model, hyper2_create, operator(+)
    use hyper2_likelihood, only : loglik_h2
    use hyper2_optimize, only : fit_result
    use hyper2_models, only : rankvec_likelihood
    use partitions, only : perms
    implicit none
    private

    type, public :: name_group
        character(len=name_len), allocatable :: items(:)
    end type name_group

    type, public :: hyper2_list
        type(hyper2_model), allocatable :: models(:)
    contains
        procedure :: size => suplist_size
    end type hyper2_list

    type, public :: lsl_component
        type(hyper2_list) :: alternatives
        real(dp) :: power = 1.0_dp
    end type lsl_component

    type, public :: lsl_model
        type(lsl_component), allocatable :: components(:)
    end type lsl_model

    public :: as_suplist, general_grouped_rank_likelihood
    public :: like_single_list, loglik_suplist, suplist_add, suplist_scale
    public :: lsl_create, lsl_add_component, loglik_lsl, maxp_lsl

    type :: int_matrix_box
        integer, allocatable :: a(:,:)
    end type int_matrix_box

contains

    integer function suplist_size(self) result(n)
        class(hyper2_list), intent(in) :: self
        if (allocated(self%models)) then
            n = size(self%models)
        else
            n = 0
        end if
    end function suplist_size

    function as_suplist(h) result(s)
        type(hyper2_model), intent(in) :: h
        type(hyper2_list) :: s
        allocate(s%models(1))
        s%models(1) = h
    end function as_suplist

    real(dp) function like_single_list(p, s) result(v)
        real(dp), intent(in) :: p(:)
        type(hyper2_list), intent(in) :: s
        real(dp), allocatable :: ll(:)
        real(dp) :: m
        integer :: i

        if (.not. allocated(s%models) .or. size(s%models) == 0) then
            v = 0.0_dp
            return
        end if
        allocate(ll(size(s%models)))
        do i = 1, size(s%models)
            ll(i) = loglik_h2(p, s%models(i))
        end do
        m = maxval(ll)
        if (m <= -huge(1.0_dp)/4.0_dp) then
            v = 0.0_dp
        else
            v = exp(m)*sum(exp(ll-m))
        end if
    end function like_single_list

    real(dp) function loglik_suplist(p, s) result(v)
        real(dp), intent(in) :: p(:)
        type(hyper2_list), intent(in) :: s
        real(dp) :: z
        z = like_single_list(p, s)
        if (z > 0.0_dp) then
            v = log(z)
        else
            v = -huge(1.0_dp)
        end if
    end function loglik_suplist

    function suplist_add(a, b) result(c)
        type(hyper2_list), intent(in) :: a, b
        type(hyper2_list) :: c
        integer :: i, j, k, na, nb

        na = a%size()
        nb = b%size()
        allocate(c%models(na*nb))
        k = 0
        do j = 1, nb
            do i = 1, na
                k = k + 1
                c%models(k) = a%models(i) + b%models(j)
            end do
        end do
    end function suplist_add

    function suplist_scale(a, x) result(c)
        type(hyper2_list), intent(in) :: a
        integer, intent(in) :: x
        type(hyper2_list) :: c, tmp
        type(hyper2_model) :: zero
        integer :: i

        if (x < 0) error stop "suplist_scale: multiplier must be nonnegative"
        if (a%size() == 0) then
            allocate(c%models(0))
        else if (x == 0) then
            zero = a%models(1)
            call zero%scale(0.0_dp)
            allocate(c%models(1))
            c%models(1) = zero
        else
            c = a
            do i = 2, x
                tmp = suplist_add(c, a)
                c = tmp
            end do
        end if
    end function suplist_scale

    function general_grouped_rank_likelihood(groups, nonfinishers) result(s)
        type(name_group), intent(in) :: groups(:)
        character(len=*), intent(in), optional :: nonfinishers(:)
        type(hyper2_list) :: s
        type(int_matrix_box), allocatable :: pm(:)
        character(len=name_len), allocatable :: v(:), common_names(:)
        integer, allocatable :: radix(:), which(:)
        integer :: g, i, j, k, col, norder, nplayer, offset, q, nextra

        allocate(pm(size(groups)), radix(size(groups)), which(size(groups)))
        nplayer = 0
        norder = 1
        do g = 1, size(groups)
            q = size(groups(g)%items)
            if (q < 1) error stop "general_grouped_rank_likelihood: empty group"
            pm(g)%a = perms(q)
            radix(g) = size(pm(g)%a,2)
            if (radix(g) > 0 .and. norder > huge(norder)/radix(g)) then
                error stop "general_grouped_rank_likelihood: too many permutations"
            end if
            norder = norder*radix(g)
            nplayer = nplayer + q
        end do

        nextra = 0
        if (present(nonfinishers)) nextra = size(nonfinishers)
        allocate(s%models(norder), v(nplayer), common_names(nplayer+nextra))
        offset = 0
        do g = 1, size(groups)
            common_names(offset+1:offset+size(groups(g)%items)) = groups(g)%items
            offset = offset + size(groups(g)%items)
        end do
        if (present(nonfinishers)) common_names(nplayer+1:nplayer+nextra) = nonfinishers
        do k = 0, norder-1
            j = k
            do g = 1, size(groups)
                which(g) = modulo(j, radix(g)) + 1
                j = j/radix(g)
            end do
            offset = 0
            do g = 1, size(groups)
                col = which(g)
                do i = 1, size(groups(g)%items)
                    v(offset+i) = groups(g)%items(pm(g)%a(i,col))
                end do
                offset = offset + size(groups(g)%items)
            end do
            if (present(nonfinishers)) then
                s%models(k+1) = canonical_players(rankvec_likelihood(v, nonfinishers), common_names)
            else
                s%models(k+1) = canonical_players(rankvec_likelihood(v), common_names)
            end if
        end do
    end function general_grouped_rank_likelihood

    function lsl_create() result(x)
        type(lsl_model) :: x
        allocate(x%components(0))
    end function lsl_create

    subroutine lsl_add_component(self, alternatives, power)
        type(lsl_model), intent(inout) :: self
        type(hyper2_list), intent(in) :: alternatives
        real(dp), intent(in), optional :: power
        type(lsl_component), allocatable :: tmp(:)
        integer :: n

        if (.not. allocated(self%components)) allocate(self%components(0))
        n = size(self%components)
        allocate(tmp(n+1))
        if (n > 0) tmp(1:n) = self%components
        tmp(n+1)%alternatives = alternatives
        tmp(n+1)%power = 1.0_dp
        if (present(power)) tmp(n+1)%power = power
        call move_alloc(tmp, self%components)
    end subroutine lsl_add_component

    real(dp) function loglik_lsl(p, x) result(v)
        real(dp), intent(in) :: p(:)
        type(lsl_model), intent(in) :: x
        real(dp) :: z
        integer :: i

        v = 0.0_dp
        if (.not. allocated(x%components)) return
        do i = 1, size(x%components)
            z = like_single_list(p, x%components(i)%alternatives)
            if (z <= 0.0_dp) then
                v = -huge(1.0_dp)
                return
            end if
            v = v + x%components(i)%power*log(z)
        end do
    end function loglik_lsl


    function canonical_players(h, names) result(out)
        type(hyper2_model), intent(in) :: h
        character(len=*), intent(in) :: names(:)
        type(hyper2_model) :: out
        character(len=name_len), allocatable :: bracket(:)
        integer :: i, j

        out = hyper2_create(names)
        do i = 1, size(h%terms)
            allocate(bracket(size(h%terms(i)%ids)))
            do j = 1, size(bracket)
                bracket(j) = h%pnames(h%terms(i)%ids(j))
            end do
            call out%add_term(bracket, h%terms(i)%power)
            deallocate(bracket)
        end do
    end function canonical_players


    function softmax_last(eta) result(p)
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable :: p(:), z(:)
        real(dp) :: m, den
        integer :: n

        n = size(eta) + 1
        allocate(p(n), z(size(eta)))
        if (size(eta) == 0) then
            p = 1.0_dp
            return
        end if
        m = max(0.0_dp, maxval(eta))
        z = exp(eta-m)
        den = exp(-m) + sum(z)
        p(1:n-1) = z/den
        p(n) = exp(-m)/den
    end function softmax_last

    function eta_from_p(p) result(eta)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: eta(:)
        integer :: n
        n = size(p)
        allocate(eta(max(0,n-1)))
        if (n > 1) eta = log(max(p(1:n-1),tiny(1.0_dp))/max(p(n),tiny(1.0_dp)))
    end function eta_from_p

    function numeric_eta_gradient(x, eta) result(g)
        type(lsl_model), intent(in) :: x
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable :: g(:), ep(:), em(:), pp(:), pm(:)
        real(dp) :: h
        integer :: i

        allocate(g(size(eta)), ep(size(eta)), em(size(eta)))
        do i = 1, size(eta)
            h = 2.0e-6_dp*max(1.0_dp,abs(eta(i)))
            ep = eta
            em = eta
            ep(i) = ep(i) + h
            em(i) = em(i) - h
            pp = softmax_last(ep)
            pm = softmax_last(em)
            g(i) = (loglik_lsl(pp,x)-loglik_lsl(pm,x))/(2.0_dp*h)
        end do
    end function numeric_eta_gradient

    function maxp_lsl(x, startp, tol, max_iter) result(res)
        type(lsl_model), intent(in) :: x
        real(dp), intent(in), optional :: startp(:), tol
        integer, intent(in), optional :: max_iter
        type(fit_result) :: res
        real(dp), allocatable :: p(:), eta(:), en(:), g(:)
        real(dp) :: ll, lln, alpha, tt
        integer :: n, it, mi
        logical :: accepted

        if (.not. allocated(x%components) .or. size(x%components) == 0) then
            res%status = 2
            allocate(res%p(0))
            return
        end if
        if (x%components(1)%alternatives%size() == 0) then
            res%status = 2
            allocate(res%p(0))
            return
        end if
        n = x%components(1)%alternatives%models(1)%size()
        allocate(p(n))
        p = 1.0_dp/real(n,dp)
        if (present(startp)) then
            if (size(startp) /= n) error stop "maxp_lsl: startp size mismatch"
            if (any(startp <= 0.0_dp) .or. abs(sum(startp)-1.0_dp) > 1.0e-8_dp) then
                error stop "maxp_lsl: invalid startp"
            end if
            p = startp
        end if
        eta = eta_from_p(p)
        ll = loglik_lsl(p,x)
        tt = 1.0e-8_dp
        if (present(tol)) tt = tol
        mi = 1000
        if (present(max_iter)) mi = max_iter

        do it = 1, mi
            g = numeric_eta_gradient(x,eta)
            if (size(g) == 0 .or. maxval(abs(g)) <= tt) then
                res%converged = .true.
                exit
            end if
            g = g/max(1.0_dp,maxval(abs(g)))
            alpha = 1.0_dp
            accepted = .false.
            do while (alpha >= 1.0e-12_dp)
                en = eta + alpha*g
                p = softmax_last(en)
                lln = loglik_lsl(p,x)
                if (lln > ll) then
                    accepted = .true.
                    exit
                end if
                alpha = alpha*0.5_dp
            end do
            if (.not. accepted) exit
            if (abs(lln-ll) <= tt) then
                eta = en
                ll = lln
                res%converged = .true.
                exit
            end if
            eta = en
            ll = lln
        end do
        res%p = softmax_last(eta)
        res%log_likelihood = loglik_lsl(res%p,x)
        res%iterations = min(it,mi)
        if (.not. res%converged) res%status = 1
    end function maxp_lsl

end module hyper2_sums
