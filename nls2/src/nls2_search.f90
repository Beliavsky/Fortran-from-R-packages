! SPDX-License-Identifier: GPL-2.0-only
module nls2_search
    use nls2_kinds, only : dp
    use nls2_types, only : nls_control, nls_result, nls2_search_result, nls_model, nls_jacobian, &
        plinear_basis, nls2_ok, nls2_bad_input, nls2_no_finite_start
    use nls2_core, only : fit_nls, fit_plinear, evaluate_model, evaluate_plinear
    use nls2_random, only : random_uniform, latin_hypercube
    implicit none
    private
    public :: nls2_fit, nls2_fit_plinear, make_grid, cpoptim_compat

contains

    subroutine nls2_fit(model, x, y, start, algorithm, result, control, weights, jacobian, lower, upper)
        procedure(nls_model) :: model
        real(dp), intent(in) :: x(:,:), y(:), start(:,:)
        character(len=*), intent(in) :: algorithm
        type(nls2_search_result), intent(out) :: result
        type(nls_control), intent(in), optional :: control
        real(dp), intent(in), optional :: weights(:), lower(:), upper(:)
        procedure(nls_jacobian), optional :: jacobian
        type(nls_control) :: ctl
        character(len=:), allocatable :: alg
        logical :: do_opt

        ctl = nls_control()
        if (present(control)) ctl = control
        alg = lower_string(trim(algorithm))
        if (alg == 'grid-search') alg = 'brute-force'
        do_opt = .not. (alg == 'brute-force' .or. alg == 'random-search' .or. alg == 'lhs')
        if (alg == 'cpoptim') then
            call cpoptim_compat(model, x, y, start, result, ctl, weights)
            return
        end if
        call prepare_starts(start, alg, ctl%maxiter, result%starts, result%status)
        if (result%status /= nls2_ok) return
        call run_model_starts(model, x, y, result%starts, do_opt, result, ctl, weights, jacobian, lower, upper)
    end subroutine nls2_fit

    subroutine nls2_fit_plinear(basis_fn, x, y, start, n_linear, algorithm, result, control, weights, lower, upper)
        procedure(plinear_basis) :: basis_fn
        real(dp), intent(in) :: x(:,:), y(:), start(:,:)
        integer, intent(in) :: n_linear
        character(len=*), intent(in) :: algorithm
        type(nls2_search_result), intent(out) :: result
        type(nls_control), intent(in), optional :: control
        real(dp), intent(in), optional :: weights(:), lower(:), upper(:)
        type(nls_control) :: ctl
        character(len=:), allocatable :: alg, base_alg
        logical :: do_opt

        ctl = nls_control()
        if (present(control)) ctl = control
        alg = lower_string(trim(algorithm))
        select case (alg)
        case ('plinear-brute-force', 'plinear-brute')
            base_alg = 'brute-force'
            do_opt = .false.
        case ('plinear-random')
            base_alg = 'random-search'
            do_opt = .false.
        case ('plinear-lhs')
            base_alg = 'lhs'
            do_opt = .false.
        case ('plinear')
            base_alg = 'default'
            do_opt = .true.
        case default
            result%status = nls2_bad_input
            return
        end select
        call prepare_starts(start, base_alg, ctl%maxiter, result%starts, result%status)
        if (result%status /= nls2_ok) return
        call run_plinear_starts(basis_fn, x, y, result%starts, n_linear, do_opt, result, ctl, weights, lower, upper)
    end subroutine nls2_fit_plinear

    subroutine prepare_starts(start, alg, maxiter, starts, status)
        real(dp), intent(in) :: start(:,:)
        character(len=*), intent(in) :: alg
        integer, intent(in) :: maxiter
        real(dp), allocatable, intent(out) :: starts(:,:)
        integer, intent(out) :: status
        integer :: p
        p = size(start,2)
        status = nls2_ok
        if (size(start,1) < 1 .or. p < 1) then
            status = nls2_bad_input
            return
        end if
        if (size(start,1) == 2) then
            if (any(start(1,:) > start(2,:))) then
                status = nls2_bad_input
                return
            end if
            select case (trim(alg))
            case ('brute-force')
                call make_grid(start(1,:), start(2,:), maxiter, starts)
            case ('lhs')
                allocate(starts(maxiter,p))
                call latin_hypercube(start(1,:), start(2,:), starts)
            case default
                ! This matches nls2.R: default and random-search both create maxiter uniform points.
                allocate(starts(maxiter,p))
                call random_uniform(start(1,:), start(2,:), starts)
            end select
        else
            ! nls2.R evaluates the supplied rows verbatim; its help page historically claimed
            ! random-search sampled them, but the current source does not do so.
            allocate(starts(size(start,1),p))
            starts = start
        end if
    end subroutine prepare_starts

    subroutine make_grid(lower, upper, maxiter, grid)
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: maxiter
        real(dp), allocatable, intent(out) :: grid(:,:)
        integer :: p, k1, k, ngrid, i, j, idx, q
        integer, allocatable :: digit(:)
        real(dp) :: frac

        p = size(lower)
        k1 = max(count(upper > lower), 1)
        k = max(ceiling(real(maxiter,dp) ** (1.0_dp / real(k1,dp))), 1)
        ngrid = 1
        do j = 1, p
            if (upper(j) > lower(j)) ngrid = ngrid * k
        end do
        allocate(grid(ngrid,p), digit(p))
        do i = 1, ngrid
            idx = i - 1
            digit = 0
            do j = 1, p
                if (upper(j) > lower(j)) then
                    q = modulo(idx, k)
                    idx = idx / k
                    digit(j) = q
                    if (k == 1) then
                        frac = 0.0_dp
                    else
                        frac = real(q,dp) / real(k-1,dp)
                    end if
                    grid(i,j) = lower(j) + frac * (upper(j)-lower(j))
                else
                    grid(i,j) = lower(j)
                end if
            end do
        end do
    end subroutine make_grid

    subroutine run_model_starts(model, x, y, starts, do_opt, result, ctl, weights, jacobian, lower, upper)
        procedure(nls_model) :: model
        real(dp), intent(in) :: x(:,:), y(:), starts(:,:)
        logical, intent(in) :: do_opt
        type(nls2_search_result), intent(inout) :: result
        type(nls_control), intent(in) :: ctl
        real(dp), intent(in), optional :: weights(:), lower(:), upper(:)
        procedure(nls_jacobian), optional :: jacobian
        type(nls_result) :: one
        integer :: i, best_i
        real(dp) :: best_rss

        allocate(result%fits(size(starts,1)), result%start_rss(size(starts,1)))
        result%n_candidates = size(starts,1)
        result%start_rss = huge(1.0_dp)
        best_rss = huge(1.0_dp)
        best_i = 0
        do i = 1, size(starts,1)
            if (do_opt) then
                if (present(jacobian)) then
                    call fit_nls(model, x, y, starts(i,:), one, ctl, weights, jacobian, lower, upper)
                else
                    call fit_nls(model, x, y, starts(i,:), one, ctl, weights=weights, lower=lower, upper=upper)
                end if
            else
                call evaluate_model(model, x, y, starts(i,:), one, weights)
            end if
            one%start_index = i
            result%fits(i) = one
            result%start_rss(i) = one%rss
            if (one%rss < best_rss .and. one%status /= nls2_bad_input) then
                best_rss = one%rss
                best_i = i
            end if
        end do
        if (best_i == 0) then
            result%status = nls2_no_finite_start
        else
            result%best = result%fits(best_i)
            result%status = nls2_ok
        end if
    end subroutine run_model_starts

    subroutine run_plinear_starts(basis_fn, x, y, starts, n_linear, do_opt, result, ctl, weights, lower, upper)
        procedure(plinear_basis) :: basis_fn
        real(dp), intent(in) :: x(:,:), y(:), starts(:,:)
        integer, intent(in) :: n_linear
        logical, intent(in) :: do_opt
        type(nls2_search_result), intent(inout) :: result
        type(nls_control), intent(in) :: ctl
        real(dp), intent(in), optional :: weights(:), lower(:), upper(:)
        type(nls_result) :: one
        integer :: i, best_i
        real(dp) :: best_rss

        allocate(result%fits(size(starts,1)), result%start_rss(size(starts,1)))
        result%n_candidates = size(starts,1)
        result%start_rss = huge(1.0_dp)
        best_rss = huge(1.0_dp)
        best_i = 0
        do i = 1, size(starts,1)
            if (do_opt) then
                call fit_plinear(basis_fn, x, y, starts(i,:), n_linear, one, ctl, weights, lower, upper)
            else
                call evaluate_plinear(basis_fn, x, y, starts(i,:), one, n_linear, weights)
            end if
            one%start_index = i
            result%fits(i) = one
            result%start_rss(i) = one%rss
            if (one%rss < best_rss .and. one%status /= nls2_bad_input) then
                best_rss = one%rss
                best_i = i
            end if
        end do
        if (best_i == 0) then
            result%status = nls2_no_finite_start
        else
            result%best = result%fits(best_i)
            result%status = nls2_ok
        end if
    end subroutine run_plinear_starts

    subroutine cpoptim_compat(model, x, y, bounds, result, ctl, weights)
        procedure(nls_model) :: model
        real(dp), intent(in) :: x(:,:), y(:), bounds(:,:)
        type(nls2_search_result), intent(out) :: result
        type(nls_control), intent(in) :: ctl
        real(dp), intent(in), optional :: weights(:)
        integer :: sample_size
        real(dp), allocatable :: starts(:,:)

        if (size(bounds,1) /= 2 .or. any(bounds(1,:) > bounds(2,:))) then
            result%status = nls2_bad_input
            return
        end if
        ! CPoptim is an external suggested R package, not code contained in nls2.
        ! This compatibility routine provides bounded space-filling sampling only;
        ! it is deliberately not advertised as an exact CPoptim translation.
        sample_size = max(1000, ctl%maxiter)
        allocate(starts(sample_size,size(bounds,2)))
        call latin_hypercube(bounds(1,:), bounds(2,:), starts)
        allocate(result%starts(sample_size,size(bounds,2)))
        result%starts = starts
        call run_model_starts(model, x, y, result%starts, .false., result, ctl, weights)
    end subroutine cpoptim_compat

    pure function lower_string(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, c
        out = s
        do i = 1, len(s)
            c = iachar(out(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
        end do
    end function lower_string

end module nls2_search
