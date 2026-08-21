module leaps

use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use leaps_lsq, only: dp, startup, includ, tolset, sing, ss, regcf, cov, reordr, &
    nobs_qr => nobs, ncol_qr => ncol, d_qr => d, rhs_qr => rhs, r_qr => r, &
    tol_qr => tol, rss_qr => rss, sserr_qr => sserr, vorder_qr => vorder, &
    tol_set_qr => tol_set, rss_set_qr => rss_set
use leaps_find_subsets, only: init_subsets, xhaust, bakwrd, forwrd, seqrep, &
    nbest_search => nbest, max_size_search => max_size, bound_best => bound, &
    ress_best => ress, lopt_best => lopt

implicit none
private

integer, parameter, public :: method_exhaustive = 1
integer, parameter, public :: method_backward = 2
integer, parameter, public :: method_forward = 3
integer, parameter, public :: method_seqrep = 4

public :: dp
public :: regsubsets_result
public :: regsubsets_fit
public :: get_model
public :: model_coefficients
public :: method_name

type :: regsubsets_result
    integer :: status = 0
    integer :: method = method_exhaustive
    integer :: nobs = 0
    integer :: nvar = 0
    integer :: ncol = 0
    integer :: nvmax = 0
    integer :: nbest = 0
    integer :: first = 0
    integer :: last = 0
    logical :: intercept = .true.
    real(dp) :: nullrss = 0.0_dp
    real(dp) :: sserr = 0.0_dp
    real(dp) :: sigma2 = 0.0_dp
    logical, allocatable :: force_in(:)
    logical, allocatable :: force_out(:)
    logical, allocatable :: dependent(:)
    integer, allocatable :: order(:)
    logical, allocatable :: valid(:,:)
    integer, allocatable :: model(:,:,:)
    real(dp), allocatable :: rss(:,:)
    real(dp), allocatable :: rsq(:,:)
    real(dp), allocatable :: adjr2(:,:)
    real(dp), allocatable :: cp(:,:)
    real(dp), allocatable :: bic(:,:)
    real(dp), allocatable :: qr_d(:)
    real(dp), allocatable :: qr_rhs(:)
    real(dp), allocatable :: qr_r(:)
    real(dp), allocatable :: qr_tol(:)
    real(dp), allocatable :: qr_rss(:)
    integer, allocatable :: qr_vorder(:)
end type regsubsets_result

contains

subroutine regsubsets_fit(x, y, result, weights, nvmax, nbest, method, &
                          intercept, force_in, force_out, ier)

    real(dp), intent(in) :: x(:,:), y(:)
    type(regsubsets_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:)
    integer, intent(in), optional :: nvmax, nbest
    character(len=*), intent(in), optional :: method
    logical, intent(in), optional :: intercept
    logical, intent(in), optional :: force_in(:), force_out(:)
    integer, intent(out), optional :: ier

    integer :: n, p, requested_nvmax, requested_nbest, meth
    integer :: i, j, pos, ierr, first, last, lastsafe, max_internal
    integer :: intercept_i, n_force_in, n_force_out
    integer, allocatable :: ord(:), neword(:)
    logical :: fit_const, reordered, nested
    logical, allocatable :: fi(:), fo(:), lindep_pos(:), right(:), dep_orig(:)
    real(dp), allocatable :: w(:)

    call clear_result(result)
    ierr = 0
    n = size(x, 1)
    p = size(x, 2)

    if (size(y) /= n) then
        ierr = 1
        call finish_error(result, ierr, ier)
        return
    end if
    if (n < 1 .or. p < 1) then
        ierr = 2
        call finish_error(result, ierr, ier)
        return
    end if

    allocate(w(n))
    w = 1.0_dp
    if (present(weights)) then
        if (size(weights) /= n) then
            ierr = 3
            call finish_error(result, ierr, ier)
            return
        end if
        w = weights
    end if
    if (any(w < 0.0_dp)) then
        ierr = 4
        call finish_error(result, ierr, ier)
        return
    end if

    fit_const = .true.
    if (present(intercept)) fit_const = intercept
    intercept_i = merge(1, 0, fit_const)

    requested_nvmax = min(8, p)
    if (present(nvmax)) requested_nvmax = nvmax
    if (requested_nvmax < 0) then
        ierr = 5
        call finish_error(result, ierr, ier)
        return
    end if
    requested_nvmax = min(requested_nvmax, p)

    requested_nbest = 1
    if (present(nbest)) requested_nbest = nbest
    if (requested_nbest < 1) then
        ierr = 6
        call finish_error(result, ierr, ier)
        return
    end if

    meth = method_exhaustive
    if (present(method)) then
        meth = parse_method(method)
        if (meth == 0) then
            ierr = 7
            call finish_error(result, ierr, ier)
            return
        end if
    end if

    allocate(fi(p), fo(p), dep_orig(p))
    fi = .false.
    fo = .false.
    dep_orig = .false.
    if (present(force_in)) then
        if (size(force_in) /= p) then
            ierr = 8
            call finish_error(result, ierr, ier)
            return
        end if
        fi = force_in
    end if
    if (present(force_out)) then
        if (size(force_out) /= p) then
            ierr = 9
            call finish_error(result, ierr, ier)
            return
        end if
        fo = force_out
    end if
    if (any(fi .and. fo)) then
        ierr = 10
        call finish_error(result, ierr, ier)
        return
    end if

    n_force_in = count(fi)
    n_force_out = count(fo)
    requested_nvmax = max(requested_nvmax, n_force_in)

    allocate(ord(p), neword(p))
    pos = 0
    do i = 1, p
        if (fi(i)) then
            pos = pos + 1
            ord(pos) = i
        end if
    end do
    do i = 1, p
        if (.not. fi(i) .and. .not. fo(i)) then
            pos = pos + 1
            ord(pos) = i
        end if
    end do
    do i = 1, p
        if (fo(i)) then
            pos = pos + 1
            ord(pos) = i
        end if
    end do

    call fit_qr(x, y, w, fit_const, ord, lindep_pos, ierr)
    if (ierr /= 0) then
        call finish_error(result, 20 + abs(ierr), ier)
        return
    end if

    do j = 1, p
        pos = j + intercept_i
        if (lindep_pos(pos)) dep_orig(ord(j)) = .true.
    end do
    if (any(dep_orig .and. fi)) then
        ierr = 11
        call finish_error(result, ierr, ier)
        return
    end if

    allocate(right(p + intercept_i))
    call build_right_mask(ord, fo, dep_orig, fit_const, right)
    reordered = has_true_before_false(right)

    if (reordered) then
        pos = 0
        do j = 1, p
            i = ord(j)
            if (.not. dep_orig(i) .and. .not. fo(i)) then
                pos = pos + 1
                neword(pos) = i
            end if
        end do
        do j = 1, p
            i = ord(j)
            if (dep_orig(i) .or. fo(i)) then
                pos = pos + 1
                neword(pos) = i
            end if
        end do
        ord = neword
        call fit_qr(x, y, w, fit_const, ord, lindep_pos, ierr)
        if (ierr /= 0) then
            call finish_error(result, 40 + abs(ierr), ier)
            return
        end if
        dep_orig = .false.
        do j = 1, p
            pos = j + intercept_i
            if (lindep_pos(pos)) dep_orig(ord(j)) = .true.
        end do
        if (any(dep_orig .and. fi)) then
            ierr = 11
            call finish_error(result, ierr, ier)
            return
        end if
        call build_right_mask(ord, fo, dep_orig, fit_const, right)
    end if

    first = 1 + n_force_in + intercept_i
    last = p - n_force_out + intercept_i
    lastsafe = 0
    do i = 1, size(right)
        if (.not. right(i)) lastsafe = i
    end do
    max_internal = min(requested_nvmax + intercept_i, last, lastsafe)

    result%nobs = n
    result%nvar = p
    result%ncol = p + intercept_i
    result%nvmax = max(0, max_internal - intercept_i)
    result%nbest = requested_nbest
    result%first = first
    result%last = last
    result%intercept = fit_const
    result%method = meth
    allocate(result%force_in(p), result%force_out(p), result%dependent(p))
    result%force_in = fi
    result%force_out = fo
    result%dependent = dep_orig

    if (max_internal < 1) then
        ierr = 12
        call finish_error(result, ierr, ier)
        return
    end if

    nbest_search = requested_nbest
    call init_subsets(max_internal - intercept_i, fit_const)

    result%sserr = sserr_qr
    if (fit_const) then
        result%nullrss = rss_qr(1)
    else
        result%nullrss = sum(y * y)
    end if
    if (n - last > 0) then
        result%sigma2 = result%sserr / real(n - last, dp)
    else
        result%sigma2 = ieee_value(0.0_dp, ieee_quiet_nan)
    end if

    call save_qr_state(result)

    if (max_internal >= first .and. first <= last) then
        nested = requested_nbest == 1 .and. &
                 (meth == method_backward .or. meth == method_forward)
        if (nested) nbest_search = 0

        select case (meth)
        case (method_exhaustive)
            call xhaust(first, last, ierr)
        case (method_backward)
            call bakwrd(first, last, ierr)
        case (method_forward)
            call forwrd(first, last, ierr)
        case (method_seqrep)
            call seqrep(first, last, ierr)
        end select

        if (ierr /= 0) then
            result%status = ierr
            if (present(ier)) ier = ierr
            return
        end if

        if (nested) then
            nbest_search = 1
            call rebuild_nested_results(first, max_internal)
        end if
    end if

    call copy_best_results(result, first, max_internal)
    call compute_metrics(result)
    result%status = 0
    if (present(ier)) ier = 0
end subroutine regsubsets_fit


subroutine get_model(result, size_predictors, rank, ids, rss, ier)

    type(regsubsets_result), intent(in) :: result
    integer, intent(in) :: size_predictors, rank
    integer, allocatable, intent(out) :: ids(:)
    real(dp), intent(out), optional :: rss
    integer, intent(out), optional :: ier

    integer :: p, ierr

    ierr = 0
    p = size_predictors + merge(1, 0, result%intercept)
    if (.not. allocated(result%valid)) then
        ierr = 1
    else if (p < 1 .or. p > size(result%valid, 1)) then
        ierr = 2
    else if (rank < 1 .or. rank > size(result%valid, 2)) then
        ierr = 3
    else if (.not. result%valid(p, rank)) then
        ierr = 4
    end if

    if (ierr /= 0) then
        allocate(ids(0))
        if (present(rss)) rss = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(ier)) ier = ierr
        return
    end if

    allocate(ids(p))
    ids = result%model(1:p, p, rank)
    if (present(rss)) rss = result%rss(p, rank)
    if (present(ier)) ier = 0
end subroutine get_model


subroutine model_coefficients(result, size_predictors, rank, beta, ids, ier, vcov)

    type(regsubsets_result), intent(in) :: result
    integer, intent(in) :: size_predictors, rank
    real(dp), allocatable, intent(out) :: beta(:)
    integer, allocatable, intent(out) :: ids(:)
    integer, intent(out), optional :: ier
    real(dp), allocatable, intent(out), optional :: vcov(:,:)

    integer :: p, ierr, ifault, dimcov, pos, row, col
    real(dp) :: var
    real(dp), allocatable :: packed(:), sterr(:)

    call get_model(result, size_predictors, rank, ids, ier=ierr)
    if (ierr /= 0) then
        allocate(beta(0))
        if (present(vcov)) allocate(vcov(0,0))
        if (present(ier)) ier = ierr
        return
    end if

    p = size(ids)
    call restore_qr_state(result)
    call reordr(ids, p, 1, ifault)
    if (ifault /= 0) then
        allocate(beta(0))
        if (present(vcov)) allocate(vcov(0,0))
        if (present(ier)) ier = 10 + abs(ifault)
        return
    end if

    allocate(beta(p))
    call regcf(beta, p, ifault)
    if (ifault /= 0) then
        deallocate(beta)
        allocate(beta(0))
        if (present(vcov)) allocate(vcov(0,0))
        if (present(ier)) ier = 20 + abs(ifault)
        return
    end if

    if (present(vcov)) then
        dimcov = p * (p + 1) / 2
        allocate(packed(dimcov), sterr(p), vcov(p,p))
        call cov(p, var, packed, dimcov, sterr, ifault)
        if (ifault /= 0) then
            deallocate(vcov)
            allocate(vcov(0,0))
            if (present(ier)) ier = 30 + abs(ifault)
            return
        end if
        vcov = 0.0_dp
        pos = 0
        do row = 1, p
            do col = row, p
                pos = pos + 1
                vcov(row,col) = packed(pos)
                vcov(col,row) = packed(pos)
            end do
        end do
    end if

    if (present(ier)) ier = 0
end subroutine model_coefficients


function method_name(method) result(name)

    integer, intent(in) :: method
    character(len=24) :: name

    select case (method)
    case (method_exhaustive)
        name = 'exhaustive'
    case (method_backward)
        name = 'backward'
    case (method_forward)
        name = 'forward'
    case (method_seqrep)
        name = 'seqrep'
    case default
        name = 'unknown'
    end select
end function method_name


subroutine fit_qr(x, y, w, fit_const, ord, lindep, ierr)

    real(dp), intent(in) :: x(:,:), y(:), w(:)
    logical, intent(in) :: fit_const
    integer, intent(in) :: ord(:)
    logical, allocatable, intent(out) :: lindep(:)
    integer, intent(out) :: ierr

    integer :: i, j, offset
    real(dp), allocatable :: row(:)

    call startup(size(x,2), fit_const)
    allocate(row(ncol_qr), lindep(ncol_qr))
    offset = merge(1, 0, fit_const)
    if (fit_const) vorder_qr(1) = 0
    do j = 1, size(ord)
        vorder_qr(j + offset) = ord(j)
    end do

    do i = 1, size(x,1)
        if (fit_const) row(1) = 1.0_dp
        do j = 1, size(ord)
            row(j + offset) = x(i,ord(j))
        end do
        call includ(w(i), row, y(i))
    end do

    call tolset(1.0e-14_dp)
    call sing(lindep, ierr)
    if (ierr < 0) ierr = 0
    call ss()
end subroutine fit_qr


subroutine build_right_mask(ord, fo, dep_orig, fit_const, right)

    integer, intent(in) :: ord(:)
    logical, intent(in) :: fo(:), dep_orig(:), fit_const
    logical, intent(out) :: right(:)

    integer :: j, offset

    right = .false.
    offset = merge(1, 0, fit_const)
    do j = 1, size(ord)
        right(j + offset) = dep_orig(ord(j)) .or. fo(ord(j))
    end do
end subroutine build_right_mask


logical function has_true_before_false(mask)

    logical, intent(in) :: mask(:)
    integer :: i
    logical :: seen_true

    seen_true = .false.
    has_true_before_false = .false.
    do i = 1, size(mask)
        if (mask(i)) then
            seen_true = .true.
        else if (seen_true) then
            has_true_before_false = .true.
            return
        end if
    end do
end function has_true_before_false


subroutine save_qr_state(result)

    type(regsubsets_result), intent(inout) :: result

    allocate(result%qr_d(size(d_qr)))
    allocate(result%qr_rhs(size(rhs_qr)))
    allocate(result%qr_r(size(r_qr)))
    allocate(result%qr_tol(size(tol_qr)))
    allocate(result%qr_rss(size(rss_qr)))
    allocate(result%qr_vorder(size(vorder_qr)))
    result%qr_d = d_qr
    result%qr_rhs = rhs_qr
    result%qr_r = r_qr
    result%qr_tol = tol_qr
    result%qr_rss = rss_qr
    result%qr_vorder = vorder_qr
    allocate(result%order(size(vorder_qr)))
    result%order = vorder_qr
end subroutine save_qr_state


subroutine restore_qr_state(result)

    type(regsubsets_result), intent(in) :: result

    call startup(result%nvar, result%intercept)
    nobs_qr = result%nobs
    d_qr = result%qr_d
    rhs_qr = result%qr_rhs
    r_qr = result%qr_r
    tol_qr = result%qr_tol
    rss_qr = result%qr_rss
    vorder_qr = result%qr_vorder
    sserr_qr = result%sserr
    tol_set_qr = .true.
    rss_set_qr = .true.
end subroutine restore_qr_state


subroutine rebuild_nested_results(first, max_internal)

    integer, intent(in) :: first, max_internal
    integer :: p, pos1
    integer, allocatable :: ids(:)

    ress_best = huge(1.0_dp)
    lopt_best = 0
    bound_best = huge(1.0_dp)

    do p = first, max_internal
        allocate(ids(p))
        ids = vorder_qr(1:p)
        call sort_int(ids)
        pos1 = p * (p - 1) / 2 + 1
        ress_best(p,1) = rss_qr(p)
        lopt_best(pos1:pos1+p-1,1) = ids
        bound_best(p) = rss_qr(p)
        deallocate(ids)
    end do
end subroutine rebuild_nested_results


subroutine copy_best_results(result, first, max_internal)

    type(regsubsets_result), intent(inout) :: result
    integer, intent(in) :: first, max_internal

    integer :: p, rank, pos1, nbest_local
    real(dp) :: big

    nbest_local = result%nbest
    allocate(result%valid(max_internal, nbest_local))
    allocate(result%model(max_internal, max_internal, nbest_local))
    allocate(result%rss(max_internal, nbest_local))
    allocate(result%rsq(max_internal, nbest_local))
    allocate(result%adjr2(max_internal, nbest_local))
    allocate(result%cp(max_internal, nbest_local))
    allocate(result%bic(max_internal, nbest_local))
    result%valid = .false.
    result%model = 0
    big = huge(1.0_dp)
    result%rss = big
    result%rsq = ieee_value(0.0_dp, ieee_quiet_nan)
    result%adjr2 = ieee_value(0.0_dp, ieee_quiet_nan)
    result%cp = ieee_value(0.0_dp, ieee_quiet_nan)
    result%bic = ieee_value(0.0_dp, ieee_quiet_nan)

    do p = max(1, first), max_internal
        pos1 = p * (p - 1) / 2 + 1
        do rank = 1, nbest_local
            if (ress_best(p,rank) < big * 0.5_dp) then
                result%valid(p,rank) = .true.
                result%rss(p,rank) = ress_best(p,rank)
                result%model(1:p,p,rank) = lopt_best(pos1:pos1+p-1,rank)
            end if
        end do
    end do
end subroutine copy_best_results


subroutine compute_metrics(result)

    type(regsubsets_result), intent(inout) :: result

    integer :: p, rank, n1, intercept_i
    real(dp) :: vr, nn, nanv

    nanv = ieee_value(0.0_dp, ieee_quiet_nan)
    intercept_i = merge(1, 0, result%intercept)
    n1 = result%nobs - intercept_i
    nn = real(result%nobs, dp)

    do p = 1, size(result%valid,1)
        do rank = 1, size(result%valid,2)
            if (.not. result%valid(p,rank)) cycle
            if (result%nullrss > 0.0_dp) then
                vr = result%rss(p,rank) / result%nullrss
                result%rsq(p,rank) = 1.0_dp - vr
                if (n1 + intercept_i - p > 0) then
                    result%adjr2(p,rank) = 1.0_dp - vr * real(n1,dp) / &
                                           real(n1 + intercept_i - p,dp)
                else
                    result%adjr2(p,rank) = nanv
                end if
                if (vr > 0.0_dp) then
                    result%bic(p,rank) = nn * log(vr) + real(p,dp) * log(nn)
                else
                    result%bic(p,rank) = -huge(1.0_dp)
                end if
            else
                result%rsq(p,rank) = nanv
                result%adjr2(p,rank) = nanv
                result%bic(p,rank) = nanv
            end if
            if (result%sigma2 > 0.0_dp) then
                result%cp(p,rank) = result%rss(p,rank) / result%sigma2 - &
                                    real(result%nobs - 2*p, dp)
            else
                result%cp(p,rank) = nanv
            end if
        end do
    end do
end subroutine compute_metrics


integer function parse_method(text)

    character(len=*), intent(in) :: text
    character(len=:), allocatable :: word

    word = trim(lower_ascii(adjustl(text)))
    select case (word)
    case ('exhaustive')
        parse_method = method_exhaustive
    case ('backward')
        parse_method = method_backward
    case ('forward')
        parse_method = method_forward
    case ('seqrep', 'sequential replacement', 'sequential_replacement')
        parse_method = method_seqrep
    case default
        parse_method = 0
    end select
end function parse_method


function lower_ascii(text) result(out)

    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code

    out = text
    do i = 1, len(text)
        code = iachar(out(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
            out(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
    end do
end function lower_ascii


subroutine sort_int(a)

    integer, intent(inout) :: a(:)
    integer :: i, j, key

    do i = 2, size(a)
        key = a(i)
        j = i - 1
        do while (j >= 1)
            if (a(j) <= key) exit
            a(j+1) = a(j)
            j = j - 1
        end do
        a(j+1) = key
    end do
end subroutine sort_int


subroutine clear_result(result)

    type(regsubsets_result), intent(out) :: result

    result%status = 0
end subroutine clear_result


subroutine finish_error(result, ierr, ier)

    type(regsubsets_result), intent(inout) :: result
    integer, intent(in) :: ierr
    integer, intent(out), optional :: ier

    result%status = ierr
    if (present(ier)) ier = ierr
end subroutine finish_error

end module leaps
