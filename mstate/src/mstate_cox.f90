module mstate_cox
    use mstate_kinds, only : dp
    use mstate_types, only : hazard_type
    use survival_types, only : coxph_result
    use survival_linalg, only : solve_sym, invert_matrix, matrix_rank
    use mstate_msfit, only : msfit_from_cox_arrays
    implicit none
    private
    public :: coxph_fit_stratified_counting, msfit_from_survival_cox, msfit_cox

contains

    subroutine coxph_fit_stratified_counting(start, stop, status, x, stratum, result, method, &
                                              weights, offset, maxiter, eps)
        real(dp), intent(in) :: start(:), stop(:), x(:, :)
        integer, intent(in) :: status(:), stratum(:)
        type(coxph_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        real(dp), intent(in), optional :: weights(:), offset(:)
        integer, intent(in), optional :: maxiter
        real(dp), intent(in), optional :: eps
        integer :: n, p, iter, max_it
        real(dp) :: tolerance, ll, newll, alpha
        real(dp), allocatable :: beta(:), score(:), info(:, :), step(:), w(:), off(:), means(:), vinv(:, :)
        logical :: ok
        character(len=12) :: tie_method

        n = size(stop)
        p = size(x, 2)
        max_it = 20
        if (present(maxiter)) max_it = maxiter
        tolerance = 1.0e-9_dp
        if (present(eps)) tolerance = eps
        tie_method = 'efron'
        if (present(method)) tie_method = adjustl(method)

        allocate(beta(p), score(p), info(p,p), step(p), w(n), off(n), means(p), vinv(p,p))
        beta = 0.0_dp
        if (present(weights)) then
            w = weights
        else
            w = 1.0_dp
        end if
        if (present(offset)) then
            off = offset
        else
            off = 0.0_dp
        end if
        if (p > 0) then
            means = matmul(transpose(x), w) / max(sum(w), tiny(1.0_dp))
        end if

        call stratified_cox_stats(start, stop, status, x, stratum, beta, w, off, tie_method, ll, score, info)
        result%loglik_initial = ll
        if (p > 0) then
            call solve_sym(info, score, step, ok)
            if (ok) result%score_test = dot_product(score, step)
        else
            ok = .true.
            result%score_test = 0.0_dp
        end if

        result%converged = p == 0
        iter = 0
        if (p > 0) then
            do iter = 1, max_it
                call solve_sym(info, score, step, ok)
                if (.not. ok) exit
                alpha = 1.0_dp
                do
                    call stratified_cox_stats(start, stop, status, x, stratum, beta + alpha*step, w, off, &
                                              tie_method, newll, score, info)
                    if (newll >= ll .or. alpha < 1.0e-8_dp) exit
                    alpha = alpha / 2.0_dp
                end do
                beta = beta + alpha*step
                if (abs(newll - ll) <= tolerance*(1.0_dp + abs(newll))) then
                    result%converged = .true.
                    ll = newll
                    exit
                end if
                ll = newll
                call stratified_cox_stats(start, stop, status, x, stratum, beta, w, off, tie_method, &
                                          ll, score, info)
            end do
        end if

        call stratified_cox_stats(start, stop, status, x, stratum, beta, w, off, tie_method, ll, score, info)
        if (p > 0) then
            call invert_matrix(info, vinv, ok)
            if (.not. ok) vinv = 0.0_dp
        end if
        allocate(result%coef(p), result%var(p,p), result%score(p), result%means(p))
        result%coef = beta
        result%score = score
        result%means = means
        if (p > 0) result%var = vinv
        result%loglik = ll
        result%iterations = min(iter, max_it)
        if (p > 0) then
            result%rank = matrix_rank(info)
        else
            result%rank = 0
        end if
    end subroutine coxph_fit_stratified_counting

    subroutine stratified_cox_stats(start, stop, status, x, stratum, beta, w, off, method, loglik, score, info)
        real(dp), intent(in) :: start(:), stop(:), x(:, :), beta(:), w(:), off(:)
        integer, intent(in) :: status(:), stratum(:)
        character(len=*), intent(in) :: method
        real(dp), intent(out) :: loglik, score(:), info(:, :)
        integer :: n, p, i, j, k, s, ns, ndeath, ntime
        integer, allocatable :: strata_values(:)
        real(dp) :: t, denom, event_denom, event_weight, frac, eta, rr, part
        real(dp), allocatable :: event_times(:), a(:), a2(:, :), d1(:), d2(:, :), meanx(:)

        n = size(stop)
        p = size(beta)
        call unique_int(stratum, strata_values)
        ns = size(strata_values)
        allocate(event_times(n), a(p), a2(p,p), d1(p), d2(p,p), meanx(p))
        loglik = 0.0_dp
        score = 0.0_dp
        info = 0.0_dp
        do s = 1, ns
            ntime = 0
            do i = 1, n
                if (stratum(i) /= strata_values(s) .or. status(i) == 0) cycle
                if (.not. contains_time(event_times, ntime, stop(i))) then
                    ntime = ntime + 1
                    event_times(ntime) = stop(i)
                end if
            end do
            if (ntime == 0) cycle
            call sort_real(event_times(1:ntime))
            do j = 1, ntime
                t = event_times(j)
                denom = 0.0_dp
                a = 0.0_dp
                a2 = 0.0_dp
                event_weight = 0.0_dp
                event_denom = 0.0_dp
                d1 = 0.0_dp
                d2 = 0.0_dp
                ndeath = 0
                do i = 1, n
                    if (stratum(i) /= strata_values(s)) cycle
                    eta = dot_product(x(i,:), beta) + off(i)
                    rr = w(i)*exp(min(eta, 700.0_dp))
                    if (start(i) < t .and. stop(i) >= t) then
                        denom = denom + rr
                        if (p > 0) then
                            a = a + rr*x(i,:)
                            a2 = a2 + rr*outer(x(i,:), x(i,:))
                        end if
                    end if
                    if (status(i) /= 0 .and. same_time(stop(i), t)) then
                        ndeath = ndeath + 1
                        event_weight = event_weight + w(i)
                        loglik = loglik + w(i)*eta
                        if (p > 0) score = score + w(i)*x(i,:)
                        event_denom = event_denom + rr
                        if (p > 0) then
                            d1 = d1 + rr*x(i,:)
                            d2 = d2 + rr*outer(x(i,:), x(i,:))
                        end if
                    end if
                end do
                if (ndeath == 0 .or. denom <= 0.0_dp) cycle
                if (index(method, 'efron') == 1 .and. ndeath > 1) then
                    do k = 0, ndeath - 1
                        frac = real(k,dp)/real(ndeath,dp)
                        part = denom - frac*event_denom
                        loglik = loglik - (event_weight/real(ndeath,dp))*log(part)
                        if (p > 0) then
                            meanx = (a - frac*d1)/part
                            score = score - (event_weight/real(ndeath,dp))*meanx
                            info = info + (event_weight/real(ndeath,dp))* &
                                   ((a2 - frac*d2)/part - outer(meanx,meanx))
                        end if
                    end do
                else
                    loglik = loglik - event_weight*log(denom)
                    if (p > 0) then
                        meanx = a/denom
                        score = score - event_weight*meanx
                        info = info + event_weight*(a2/denom - outer(meanx,meanx))
                    end if
                end if
            end do
        end do
    end subroutine stratified_cox_stats

    subroutine msfit_from_survival_cox(fit, start, stop, event, xmat, strata, kstrata, times, newx, &
                                       method, variance, hazout, info, offset, newoffset)
        type(coxph_result), intent(in) :: fit
        real(dp), intent(in) :: start(:), stop(:), xmat(:, :), times(:), newx(:, :)
        integer, intent(in) :: event(:), strata(:), kstrata(:), method
        logical, intent(in) :: variance
        type(hazard_type), intent(out) :: hazout
        integer, intent(out), optional :: info
        real(dp), intent(in), optional :: offset(:), newoffset(:)
        call msfit_from_cox_arrays(start, stop, event, xmat, fit%coef, fit%var, strata, kstrata, times, &
                                   newx, method, variance, hazout, info, offset, newoffset)
    end subroutine msfit_from_survival_cox

    subroutine msfit_cox(start, stop, event, xmat, stratum, newx, kstrata, times, method, variance, &
                         fit, hazout, info, offset, newoffset)
        real(dp), intent(in) :: start(:), stop(:), xmat(:, :), newx(:, :), times(:)
        integer, intent(in) :: event(:), stratum(:), kstrata(:), method
        logical, intent(in) :: variance
        type(coxph_result), intent(out) :: fit
        type(hazard_type), intent(out) :: hazout
        integer, intent(out), optional :: info
        real(dp), intent(in), optional :: offset(:), newoffset(:)
        integer, allocatable :: perm(:), sbounds(:), estate(:), stratum_sorted(:)
        real(dp), allocatable :: sstart(:), sstop(:), sx(:, :), soff(:)
        integer :: n, i, h, ierr, conv_info
        character(len=12) :: mth

        if (present(info)) info = 0
        ierr = 0
        conv_info = 0
        n = size(stop)
        mth = 'breslow'
        if (method == 2) mth = 'efron'
        if (present(offset)) then
            call coxph_fit_stratified_counting(start, stop, event, xmat, stratum, fit, mth, offset=offset)
        else
            call coxph_fit_stratified_counting(start, stop, event, xmat, stratum, fit, mth)
        end if
        if (.not. fit%converged .and. size(fit%coef) > 0) conv_info = 10
        if (minval(stratum) < 1) then
            if (present(info)) info = 11
            return
        end if
        call order_by_int(stratum, perm)
        allocate(sstart(n), sstop(n), estate(n), stratum_sorted(n), sx(n,size(xmat,2)))
        sstart = start(perm); sstop = stop(perm); estate = event(perm); stratum_sorted = stratum(perm)
        if (size(xmat,2) > 0) sx = xmat(perm,:)
        if (present(offset)) then
            allocate(soff(n)); soff = offset(perm)
        end if
        h = maxval(stratum)
        allocate(sbounds(h+1)); sbounds = 1
        do i = 1, h
            sbounds(i) = 1 + count(stratum_sorted < i)
        end do
        sbounds(h+1) = n + 1
        if (present(offset)) then
            call msfit_from_survival_cox(fit, sstart, sstop, estate, sx, sbounds, kstrata, times, newx, &
                                         method, variance, hazout, ierr, soff, newoffset)
        else
            call msfit_from_survival_cox(fit, sstart, sstop, estate, sx, sbounds, kstrata, times, newx, &
                                         method, variance, hazout, ierr, newoffset=newoffset)
        end if
        if (present(info)) then
            if (ierr /= 0) then
                info = ierr
            else
                info = conv_info
            end if
        end if
    end subroutine msfit_cox

    pure function outer(a, b) result(c)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: c(size(a),size(b))
        integer :: i
        do i = 1, size(a)
            c(i,:) = a(i)*b
        end do
    end function outer

    pure logical function same_time(a, b)
        real(dp), intent(in) :: a, b
        same_time = abs(a-b) <= 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(a),abs(b))
    end function same_time

    logical function contains_time(x, n, value)
        real(dp), intent(in) :: x(:), value
        integer, intent(in) :: n
        integer :: i
        contains_time = .false.
        do i = 1, n
            if (same_time(x(i), value)) then
                contains_time = .true.
                return
            end if
        end do
    end function contains_time

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: key
        do i = 2, size(x)
            key = x(i); j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j+1) = x(j); j = j - 1
            end do
            x(j+1) = key
        end do
    end subroutine sort_real

    subroutine unique_int(x, u)
        integer, intent(in) :: x(:)
        integer, allocatable, intent(out) :: u(:)
        integer, allocatable :: temp(:)
        integer :: i, n
        allocate(temp(size(x))); n = 0
        do i = 1, size(x)
            if (n == 0 .or. .not. any(temp(1:n) == x(i))) then
                n = n + 1; temp(n) = x(i)
            end if
        end do
        allocate(u(n)); if (n > 0) u = temp(1:n)
    end subroutine unique_int

    subroutine order_by_int(x, perm)
        integer, intent(in) :: x(:)
        integer, allocatable, intent(out) :: perm(:)
        integer :: i, j, key
        allocate(perm(size(x)))
        perm = [(i, i=1,size(x))]
        do i = 2, size(x)
            key = perm(i); j = i - 1
            do while (j >= 1)
                if (x(perm(j)) <= x(key)) exit
                perm(j+1) = perm(j); j = j - 1
            end do
            perm(j+1) = key
        end do
    end subroutine order_by_int

end module mstate_cox
