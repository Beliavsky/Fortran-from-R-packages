! Computational translation of the R package MM 1.7-0.
! Upstream license: GPL-2. This translation is GPL-2.
module mm_core
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf
    use mm_kinds, only : dp
    use mm_types, only : paras_type, suffstats_type
    use mm_parameters, only : paras_dimension, p, theta, valid_paras
    use partitions, only : compositions
    use quadform, only : quad_tdiag
    implicit none
    private

    public :: lmultinomial, multinomial
    public :: mm_single_log, mm_single, normc_log, normc, dmm
    public :: mm_loglik, mm_allsamesum, mm_differsums
    public :: mm_allsamesum_a, mm_differsums_a, mm_support
    public :: suffstats, expected_suffstats

contains

    pure real(dp) function lmultinomial(x) result(ans)
        integer, intent(in) :: x(:)
        integer :: i, n

        if (any(x < 0)) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
            return
        end if
        n = sum(x)
        ans = log_gamma(real(n + 1, dp))
        do i = 1, size(x)
            ans = ans - log_gamma(real(x(i) + 1, dp))
        end do
    end function lmultinomial

    pure real(dp) function multinomial(x) result(ans)
        integer, intent(in) :: x(:)
        ans = exp(lmultinomial(x))
    end function multinomial

    real(dp) function mm_single_log(yrow, par) result(ans)
        integer, intent(in) :: yrow(:)
        type(paras_type), intent(in) :: par
        real(dp), allocatable :: prob(:), th(:,:), lm(:,:), yr(:,:), qv(:)
        real(dp) :: lp
        integer :: k, i, j

        k = paras_dimension(par)
        if (size(yrow) /= k) error stop "mm_single_log: yrow has wrong length"
        if (any(yrow < 0)) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
            return
        end if
        prob = p(par)
        th = theta(par)
        if (any(prob < 0.0_dp) .or. any(th <= 0.0_dp)) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
            return
        end if

        lp = 0.0_dp
        do i = 1, k
            if (yrow(i) > 0) then
                if (prob(i) <= 0.0_dp) then
                    ans = ieee_value(0.0_dp, ieee_negative_inf)
                    return
                end if
                lp = lp + real(yrow(i), dp) * log(prob(i))
            end if
        end do

        allocate(lm(k, k), yr(1, k))
        lm = 0.0_dp
        do j = 2, k
            do i = 1, j - 1
                lm(i, j) = log(th(i, j))
            end do
        end do
        yr(1, :) = real(yrow, dp)
        qv = quad_tdiag(lm, yr)
        ans = lmultinomial(yrow) + lp + qv(1)
    end function mm_single_log

    real(dp) function mm_single(yrow, par) result(ans)
        integer, intent(in) :: yrow(:)
        type(paras_type), intent(in) :: par
        ans = exp(mm_single_log(yrow, par))
    end function mm_single

    real(dp) function normc_log(y_total, par) result(ans)
        integer, intent(in) :: y_total
        type(paras_type), intent(in) :: par
        integer, allocatable :: comp(:,:)
        real(dp), allocatable :: terms(:)
        real(dp) :: mx
        integer :: k, j

        if (y_total < 0) error stop "normc_log: y_total must be nonnegative"
        if (.not. valid_paras(par)) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
            return
        end if
        if (y_total == 0) then
            ans = 0.0_dp
            return
        end if
        k = paras_dimension(par)
        comp = compositions(y_total, k)
        allocate(terms(size(comp, 2)))
        do j = 1, size(comp, 2)
            terms(j) = mm_single_log(comp(:, j), par)
        end do
        mx = maxval(terms)
        ans = mx + log(sum(exp(terms - mx)))
    end function normc_log

    real(dp) function normc(y_total, par) result(ans)
        integer, intent(in) :: y_total
        type(paras_type), intent(in) :: par
        ans = exp(normc_log(y_total, par))
    end function normc

    real(dp) function dmm(yrow, par) result(ans)
        integer, intent(in) :: yrow(:)
        type(paras_type), intent(in) :: par
        ans = exp(mm_single_log(yrow, par) - normc_log(sum(yrow), par))
    end function dmm

    real(dp) function mm_loglik(y, par, n) result(ans)
        integer, intent(in) :: y(:,:)
        type(paras_type), intent(in) :: par
        real(dp), intent(in), optional :: n(:)
        integer, allocatable :: totals(:)

        if (size(y, 2) /= paras_dimension(par)) error stop "mm_loglik: wrong number of columns"
        totals = sum(y, dim=2)
        if (all(totals == totals(1))) then
            ans = mm_allsamesum(y, par, n)
        else
            ans = mm_differsums(y, par, n)
        end if
    end function mm_loglik

    real(dp) function mm_allsamesum(y, par, n) result(ans)
        integer, intent(in) :: y(:,:)
        type(paras_type), intent(in) :: par
        real(dp), intent(in), optional :: n(:)
        real(dp), allocatable :: wt(:), prob(:), th(:,:), lm(:,:), yr(:,:), qv(:)
        real(dp) :: nsum, lp
        integer :: nr, k, i, j, y_total

        nr = size(y, 1)
        k = size(y, 2)
        if (k /= paras_dimension(par)) error stop "mm_allsamesum: wrong number of columns"
        if (.not. all(sum(y, dim=2) == sum(y(1, :)))) error stop "mm_allsamesum: row sums differ"
        call make_weights(nr, n, wt)
        y_total = sum(y(1, :))
        nsum = sum(wt)
        prob = p(par)
        th = theta(par)
        if (any(prob <= 0.0_dp) .or. any(th <= 0.0_dp)) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
            return
        end if

        ans = -normc_log(y_total, par) * nsum
        ans = ans + log_gamma(real(y_total + 1, dp)) * nsum
        do i = 1, nr
            do j = 1, k
                ans = ans - wt(i) * log_gamma(real(y(i, j) + 1, dp))
            end do
        end do
        lp = 0.0_dp
        do i = 1, nr
            lp = lp + wt(i) * sum(real(y(i, :), dp) * log(prob))
        end do
        ans = ans + lp

        allocate(lm(k, k), yr(nr, k))
        lm = 0.0_dp
        do j = 2, k
            do i = 1, j - 1
                lm(i, j) = log(th(i, j))
            end do
        end do
        yr = real(y, dp)
        qv = quad_tdiag(lm, yr)
        ans = ans + sum(wt * qv)
    end function mm_allsamesum

    real(dp) function mm_differsums(y, par, n) result(ans)
        integer, intent(in) :: y(:,:)
        type(paras_type), intent(in) :: par
        real(dp), intent(in), optional :: n(:)
        real(dp), allocatable :: wt(:)
        integer :: i, nr

        nr = size(y, 1)
        if (size(y, 2) /= paras_dimension(par)) error stop "mm_differsums: wrong number of columns"
        call make_weights(nr, n, wt)
        ans = 0.0_dp
        do i = 1, nr
            ans = ans + wt(i) * (mm_single_log(y(i, :), par) - normc_log(sum(y(i, :)), par))
        end do
    end function mm_differsums

    real(dp) function mm_allsamesum_a(y, par) result(ans)
        integer, intent(in) :: y(:,:)
        type(paras_type), intent(in) :: par
        real(dp), allocatable :: prob(:), th(:,:), lm(:,:), yr(:,:)
        integer :: nr, k, i, j, y_total

        ! Compatibility with upstream MM_allsamesum_A(), including its use
        ! of NormC rather than log(NormC) in the first term.
        nr = size(y, 1)
        k = size(y, 2)
        if (.not. all(sum(y, dim=2) == sum(y(1, :)))) error stop "mm_allsamesum_a: row sums differ"
        y_total = sum(y(1, :))
        prob = p(par)
        th = theta(par)
        ans = -real(nr, dp) * normc(y_total, par)
        ans = ans + real(nr, dp) * log_gamma(real(y_total + 1, dp))
        do i = 1, nr
            do j = 1, k
                ans = ans - log_gamma(real(y(i, j) + 1, dp))
                ans = ans + real(y(i, j), dp) * log(prob(j))
            end do
        end do
        allocate(lm(k, k), yr(nr, k))
        lm = 0.0_dp
        do j = 2, k
            do i = 1, j - 1
                lm(i, j) = log(th(i, j))
            end do
        end do
        yr = real(y, dp)
        ans = ans + sum(quad_tdiag(lm, yr))
    end function mm_allsamesum_a

    real(dp) function mm_differsums_a(y, par) result(ans)
        integer, intent(in) :: y(:,:)
        type(paras_type), intent(in) :: par
        ans = mm_differsums(y, par)
    end function mm_differsums_a

    real(dp) function mm_support(par, ss) result(ans)
        type(paras_type), intent(in) :: par
        type(suffstats_type), intent(in) :: ss
        real(dp), allocatable :: prob(:), th(:,:)
        integer :: k, i, j

        k = paras_dimension(par)
        if (size(ss%row_sums) /= k) error stop "mm_support: incompatible sufficient statistics"
        prob = p(par)
        th = theta(par)
        if (any(prob <= 0.0_dp) .or. any(th <= 0.0_dp)) then
            ans = ieee_value(0.0_dp, ieee_negative_inf)
            return
        end if
        ans = -normc_log(ss%y_total, par) * ss%nobs
        ans = ans + sum(log(prob) * ss%row_sums)
        do j = 2, k
            do i = 1, j - 1
                ans = ans + log(th(i, j)) * ss%cross_prods(i, j)
            end do
        end do
    end function mm_support

    function suffstats(y, n) result(ss)
        integer, intent(in) :: y(:,:)
        real(dp), intent(in), optional :: n(:)
        type(suffstats_type) :: ss
        real(dp), allocatable :: wt(:)
        integer :: nr, k, i, j, r

        nr = size(y, 1)
        k = size(y, 2)
        if (.not. all(sum(y, dim=2) == sum(y(1, :)))) error stop "suffstats: row sums differ"
        call make_weights(nr, n, wt)
        ss%y_total = sum(y(1, :))
        ss%nobs = sum(wt)
        allocate(ss%row_sums(k), ss%cross_prods(k, k))
        ss%row_sums = 0.0_dp
        ss%cross_prods = 0.0_dp
        do r = 1, nr
            ss%row_sums = ss%row_sums + wt(r) * real(y(r, :), dp)
            do j = 1, k
                do i = 1, k
                    ss%cross_prods(i, j) = ss%cross_prods(i, j) + &
                        wt(r) * real(y(r, i) * y(r, j), dp)
                end do
            end do
        end do
    end function suffstats

    function expected_suffstats(par, y_total) result(ss)
        type(paras_type), intent(in) :: par
        integer, intent(in) :: y_total
        type(suffstats_type) :: ss
        integer, allocatable :: comp(:,:)
        real(dp), allocatable :: lp(:), pr(:)
        real(dp) :: mx
        integer :: k, i, j, c

        k = paras_dimension(par)
        ss%y_total = y_total
        ss%nobs = 1.0_dp
        allocate(ss%row_sums(k), ss%cross_prods(k, k))
        ss%row_sums = 0.0_dp
        ss%cross_prods = 0.0_dp
        if (y_total == 0) return
        comp = compositions(y_total, k)
        allocate(lp(size(comp, 2)), pr(size(comp, 2)))
        do c = 1, size(comp, 2)
            lp(c) = mm_single_log(comp(:, c), par)
        end do
        mx = maxval(lp)
        pr = exp(lp - mx)
        pr = pr / sum(pr)
        do c = 1, size(comp, 2)
            ss%row_sums = ss%row_sums + pr(c) * real(comp(:, c), dp)
            do j = 1, k
                do i = 1, k
                    ss%cross_prods(i, j) = ss%cross_prods(i, j) + &
                        pr(c) * real(comp(i, c) * comp(j, c), dp)
                end do
            end do
        end do
    end function expected_suffstats

    subroutine make_weights(nr, n, wt)
        integer, intent(in) :: nr
        real(dp), intent(in), optional :: n(:)
        real(dp), allocatable, intent(out) :: wt(:)

        allocate(wt(nr))
        if (present(n)) then
            if (size(n) /= nr) error stop "weight vector has wrong length"
            if (any(n < 0.0_dp)) error stop "weights must be nonnegative"
            wt = n
        else
            wt = 1.0_dp
        end if
    end subroutine make_weights

end module mm_core
