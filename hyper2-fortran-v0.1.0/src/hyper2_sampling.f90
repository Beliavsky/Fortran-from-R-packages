! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_sampling
    use hyper2_kinds, only : dp, name_len
    use hyper2_types, only : hyper2_model, hyper3_model, hyper2_create, hyper3_create, operator(+)
    use hyper2_likelihood, only : loglik_h2, fillup, equalp
    use hyper2_models, only : rankvec_likelihood, ordervec2supp3, zipf, trial, dirichlet
    implicit none
    private

    public :: seed_rng, rdirichlet, rp, rrank_single, rorder_single, rrace
    public :: rwinner3, rrace3, rpair3, rracehyper3, rhyper2

contains

    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729*i + 8191*i*i, huge(1)-1)
            if (put(i) <= 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine seed_rng

    real(dp) function randn() result(z)
        real(dp) :: u1, u2
        real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        z = sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
    end function randn

    recursive real(dp) function rgamma1(shape) result(x)
        real(dp), intent(in) :: shape
        real(dp) :: d, c, z, u, v

        if (shape <= 0.0_dp) then
            x = 0.0_dp
            return
        end if
        if (shape < 1.0_dp) then
            call random_number(u)
            u = max(u, tiny(1.0_dp))
            x = rgamma1(shape + 1.0_dp)*u**(1.0_dp/shape)
            return
        end if

        d = shape - 1.0_dp/3.0_dp
        c = 1.0_dp/sqrt(9.0_dp*d)
        do
            z = randn()
            v = 1.0_dp + c*z
            if (v <= 0.0_dp) cycle
            v = v*v*v
            call random_number(u)
            if (u < 1.0_dp - 0.0331_dp*z**4) exit
            if (log(max(u, tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
        end do
        x = d*v
    end function rgamma1

    function rdirichlet(n, alpha) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha(:)
        real(dp), allocatable :: x(:,:)
        real(dp) :: sx
        integer :: i, j

        if (any(alpha <= 0.0_dp)) error stop "rdirichlet: alpha must be positive"
        allocate(x(n, size(alpha)))
        do i = 1, n
            do j = 1, size(alpha)
                x(i,j) = rgamma1(alpha(j))
            end do
            sx = sum(x(i,:))
            if (sx > 0.0_dp) x(i,:) = x(i,:)/sx
        end do
    end function rdirichlet

    integer function sample_weighted(w) result(idx)
        real(dp), intent(in) :: w(:)
        real(dp) :: u, cum, s
        integer :: i

        s = sum(max(w, 0.0_dp))
        if (s <= 0.0_dp) then
            idx = 0
            return
        end if
        call random_number(u)
        u = u*s
        cum = 0.0_dp
        do i = 1, size(w)
            cum = cum + max(w(i), 0.0_dp)
            if (u <= cum) then
                idx = i
                return
            end if
        end do
        idx = size(w)
    end function sample_weighted

    function rrank_single(p) result(order)
        real(dp), intent(in) :: p(:)
        integer, allocatable :: order(:)
        real(dp), allocatable :: w(:)
        integer :: i, j

        allocate(order(size(p)))
        w = max(p, 0.0_dp)
        do i = 1, size(p)
            j = sample_weighted(w)
            if (j == 0) then
                order(i) = 0
            else
                order(i) = j
                w(j) = 0.0_dp
            end if
        end do
    end function rrank_single

    function rorder_single(p) result(ord)
        real(dp), intent(in) :: p(:)
        integer, allocatable :: ord(:), r(:)
        integer :: i

        r = rrank_single(p)
        allocate(ord(size(r)))
        ord = 0
        do i = 1, size(r)
            if (r(i) > 0) ord(r(i)) = i
        end do
    end function rorder_single

    function rrace(strengths) result(order)
        real(dp), intent(in) :: strengths(:)
        integer, allocatable :: order(:)
        order = rrank_single(strengths)
    end function rrace

    function rp(n, h, startp, sigma, small) result(samples)
        integer, intent(in) :: n
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in), optional :: startp(:), sigma, small
        real(dp), allocatable :: samples(:,:), x(:), xn(:), p(:), pn(:)
        real(dp) :: sd, sm, ll, lln, u
        integer :: i, j, m

        m = h%size() - 1
        if (m < 1) error stop "rp: model needs at least two players"
        allocate(samples(n, h%size()), x(m), xn(m))
        sd = 0.1_dp
        if (present(sigma)) sd = sigma
        sm = 1.0e-6_dp
        if (present(small)) sm = small

        x = 1.0_dp/real(h%size(), dp)
        if (present(startp)) then
            if (size(startp) == h%size()) then
                x = startp(1:m)
            else if (size(startp) == m) then
                x = startp
            else
                error stop "rp: startp has wrong size"
            end if
        end if
        p = fillup(x)
        ll = loglik_h2(p, h)

        do i = 1, n
            xn = x
            do j = 1, m
                xn(j) = x(j) + sd*randn()
            end do
            if (all(xn >= sm) .and. 1.0_dp-sum(xn) >= sm) then
                pn = fillup(xn)
                lln = loglik_h2(pn, h)
                call random_number(u)
                if (log(max(u, tiny(1.0_dp))) <= lln-ll) then
                    x = xn
                    p = pn
                    ll = lln
                end if
            end if
            samples(i,:) = p
        end do
    end function rp

    integer function rwinner3(counts, strengths) result(idx)
        integer, intent(in) :: counts(:)
        real(dp), intent(in) :: strengths(:)
        real(dp), allocatable :: w(:)

        if (size(counts) /= size(strengths)) error stop "rwinner3: size mismatch"
        allocate(w(size(counts)))
        w = real(counts, dp)*strengths
        idx = sample_weighted(w)
    end function rwinner3

    function rrace3(counts_in, strengths) result(order)
        integer, intent(in) :: counts_in(:)
        real(dp), intent(in) :: strengths(:)
        integer, allocatable :: order(:), counts(:)
        integer :: i, j

        if (any(counts_in < 0)) error stop "rrace3: negative count"
        counts = counts_in
        allocate(order(sum(counts)))
        do i = 1, size(order)
            j = rwinner3(counts, strengths)
            if (j < 1) error stop "rrace3: no available competitor"
            order(i) = j
            counts(j) = counts(j) - 1
        end do
    end function rrace3

    function rpair3(n, s, lambda) result(h)
        integer, intent(in) :: n, s
        real(dp), intent(in) :: lambda
        type(hyper3_model) :: h
        character(len=name_len), allocatable :: names(:), one(:), pair(:)
        real(dp), allocatable :: strength(:)
        real(dp) :: u, pw
        integer :: i, a, b

        if (n < 2 .or. s < 0 .or. lambda <= 0.0_dp) error stop "rpair3: invalid argument"
        allocate(names(n), one(1), pair(2))
        do i = 1, n
            write(names(i), '("p",i0)') i
        end do
        strength = zipf(n)
        h = hyper3_create(names)
        do i = 1, s
            call random_number(u)
            a = 1 + int(u*real(n,dp))
            do
                call random_number(u)
                b = 1 + int(u*real(n,dp))
                if (b /= a) exit
            end do
            pw = strength(a)*lambda/(strength(a)*lambda + strength(b))
            call random_number(u)
            if (u < pw) then
                one(1) = names(a)
                pair = [names(a), names(b)]
                call h%add_term(one, [lambda], 1.0_dp)
                call h%add_term(pair, [lambda, 1.0_dp], -1.0_dp)
            else
                one(1) = names(b)
                pair = [names(b), names(a)]
                call h%add_term(one, [lambda], 1.0_dp)
                call h%add_term(pair, [lambda, 1.0_dp], -1.0_dp)
            end if
        end do
    end function rpair3

    function rracehyper3(n, total_size, strengths, races) result(h)
        integer, intent(in) :: n, total_size, races
        real(dp), intent(in), optional :: strengths(:)
        type(hyper3_model) :: h, tmp
        real(dp), allocatable :: ps(:)
        integer, allocatable :: counts(:), ord(:)
        character(len=name_len), allocatable :: names(:), v(:)
        real(dp) :: u
        integer :: i, j

        if (n < 1 .or. total_size < 1 .or. races < 0) error stop "rracehyper3: invalid argument"
        allocate(names(n), counts(n))
        do i = 1, n
            write(names(i), '("p",i0)') i
        end do
        if (present(strengths)) then
            if (size(strengths) /= n) error stop "rracehyper3: strengths size mismatch"
            ps = strengths
        else
            ps = zipf(n)
        end if
        counts = 0
        do i = 1, total_size
            call random_number(u)
            j = min(n, 1 + int(u*real(n,dp)))
            counts(j) = counts(j) + 1
        end do
        h = hyper3_create(names)
        do i = 1, races
            ord = rrace3(counts, ps)
            allocate(v(size(ord)))
            do j = 1, size(ord)
                v(j) = names(ord(j))
            end do
            tmp = ordervec2supp3(v)
            h = h + tmp
            deallocate(v)
        end do
    end function rracehyper3

    function rhyper2(n, s, pairs, teams, race_flag) result(h)
        integer, intent(in) :: n, s
        logical, intent(in), optional :: pairs, teams, race_flag
        type(hyper2_model) :: h, tmp
        character(len=name_len), allocatable :: names(:), v(:), winners(:), one(:), pair(:)
        real(dp), allocatable :: ps(:)
        logical :: lp, lt, lr
        real(dp) :: u
        integer, allocatable :: ord(:)
        integer :: nn, i, j, a, b

        nn = n - modulo(n, 2)
        if (nn < 2 .or. s < 0) error stop "rhyper2: invalid argument"
        allocate(names(nn), one(1), pair(2))
        do i = 1, nn
            write(names(i), '("p",i0)') i
        end do
        lp = .true.
        if (present(pairs)) lp = pairs
        lt = .true.
        if (present(teams)) lt = teams
        lr = .true.
        if (present(race_flag)) lr = race_flag
        h = hyper2_create(names)

        if (lp) then
            do i = 1, s
                call random_number(u)
                a = min(nn, 1 + int(u*real(nn,dp)))
                do
                    call random_number(u)
                    b = min(nn, 1 + int(u*real(nn,dp)))
                    if (b /= a) exit
                end do
                one(1) = names(a)
                pair = [names(a), names(b)]
                tmp = trial(one, pair)
                h = h + tmp
            end do
        end if

        if (lt) then
            allocate(winners(nn/2))
            ps = equalp(nn)
            do i = 1, s
                ord = rrank_single(ps)
                do j = 1, nn/2
                    winners(j) = names(ord(j))
                end do
                tmp = trial(winners, names)
                h = h + tmp
            end do
        end if

        if (lr .and. s > 0) then
            ps = zipf(nn)
            ord = rrank_single(ps)
            allocate(v(nn))
            do j = 1, nn
                v(j) = names(ord(j))
            end do
            tmp = rankvec_likelihood(v)
            h = h + tmp
        end if
    end function rhyper2

end module hyper2_sampling
