! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

module dirmult_core
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use dirmult_types, only : dp, dirmult_fit_type, mom_result_type, dirmult_summary_type
    use dirmult_linalg, only : solve_linear, invert_matrix
    implicit none
    private

    public :: dirmult_loglik, multinomial_loglik, score_function
    public :: observed_fim, expected_fim, theta_fim
    public :: weir_mom, fit_dirmult, summarize_dirmult
    public :: clean_counts

contains

    pure function recip_sum(a, n) result(s)
        real(dp), intent(in) :: a
        integer, intent(in) :: n
        real(dp) :: s
        integer :: r
        s = 0.0_dp
        do r = 0, n - 1
            s = s + 1.0_dp / (a + real(r,dp))
        end do
    end function recip_sum

    pure function recip_square_sum(a, n) result(s)
        real(dp), intent(in) :: a
        integer, intent(in) :: n
        real(dp) :: s, z
        integer :: r
        s = 0.0_dp
        do r = 0, n - 1
            z = a + real(r,dp)
            s = s + 1.0_dp / (z*z)
        end do
    end function recip_square_sum

    pure function rising_log_sum(a, n) result(s)
        real(dp), intent(in) :: a
        integer, intent(in) :: n
        real(dp) :: s
        if (n <= 0) then
            s = 0.0_dp
        else
            s = log_gamma(a + real(n,dp)) - log_gamma(a)
        end if
    end function rising_log_sum

    subroutine clean_counts(data, remove_rows, remove_cols, x, row_index, col_index)
        integer, intent(in) :: data(:,:)
        logical, intent(in) :: remove_rows, remove_cols
        integer, allocatable, intent(out) :: x(:,:)
        integer, allocatable, intent(out), optional :: row_index(:), col_index(:)
        logical, allocatable :: keep_r(:), keep_c(:)
        integer :: nr, nc, i, j, ir, jc

        nr = size(data,1)
        nc = size(data,2)
        allocate(keep_r(nr), keep_c(nc))
        keep_r = .true.
        keep_c = .true.
        if (remove_rows) then
            do i = 1, nr
                keep_r(i) = sum(data(i,:)) /= 0
            end do
        end if
        if (remove_cols) then
            do j = 1, nc
                keep_c(j) = sum(data(:,j)) /= 0
            end do
        end if

        allocate(x(count(keep_r), count(keep_c)))
        if (present(row_index)) allocate(row_index(count(keep_r)))
        if (present(col_index)) allocate(col_index(count(keep_c)))
        ir = 0
        do i = 1, nr
            if (.not. keep_r(i)) cycle
            ir = ir + 1
            if (present(row_index)) row_index(ir) = i
            jc = 0
            do j = 1, nc
                if (.not. keep_c(j)) cycle
                jc = jc + 1
                if (ir == 1 .and. present(col_index)) col_index(jc) = j
                x(ir,jc) = data(i,j)
            end do
        end do
    end subroutine clean_counts

    function score_function(data, gamma) result(s)
        integer, intent(in) :: data(:,:)
        real(dp), intent(in) :: gamma(:)
        real(dp) :: s(size(gamma))
        real(dp) :: total_gamma, sn
        integer :: i, j, n, k

        k = size(gamma)
        s = 0.0_dp
        if (size(data,2) /= k) return
        total_gamma = sum(gamma)
        do j = 1, size(data,1)
            n = sum(data(j,:))
            sn = recip_sum(total_gamma, n)
            do i = 1, k
                s(i) = s(i) - sn + recip_sum(gamma(i), data(j,i))
            end do
        end do
    end function score_function

    function observed_fim(data, gamma) result(f)
        integer, intent(in) :: data(:,:)
        real(dp), intent(in) :: gamma(:)
        real(dp) :: f(size(gamma),size(gamma))
        real(dp) :: d(size(gamma)), od, total_gamma
        integer :: i, j, n, k

        k = size(gamma)
        f = 0.0_dp
        if (size(data,2) /= k) return
        d = 0.0_dp
        od = 0.0_dp
        total_gamma = sum(gamma)
        do j = 1, size(data,1)
            n = sum(data(j,:))
            od = od + recip_square_sum(total_gamma, n)
            do i = 1, k
                d(i) = d(i) + recip_square_sum(gamma(i), data(j,i))
            end do
        end do
        f = od
        do i = 1, k
            f(i,i) = f(i,i) - d(i)
        end do
    end function observed_fim

    pure function beta_binomial_pmf(x, n, a, b) result(p)
        integer, intent(in) :: x, n
        real(dp), intent(in) :: a, b
        real(dp) :: p, lp

        if (x < 0 .or. x > n .or. a <= 0.0_dp .or. b <= 0.0_dp) then
            p = 0.0_dp
            return
        end if
        lp = log_gamma(real(n+1,dp)) - log_gamma(real(x+1,dp)) - &
             log_gamma(real(n-x+1,dp))
        lp = lp + rising_log_sum(a, x) + rising_log_sum(b, n-x) - &
             rising_log_sum(a+b, n)
        p = exp(lp)
    end function beta_binomial_pmf

    function expected_fim(data, gamma) result(f)
        integer, intent(in) :: data(:,:)
        real(dp), intent(in) :: gamma(:)
        real(dp) :: f(size(gamma),size(gamma))
        real(dp) :: inner(size(gamma)), total_gamma, od, tail, z
        integer :: j, k, r, n, kk

        kk = size(gamma)
        f = 0.0_dp
        if (size(data,2) /= kk .or. kk < 2) return
        inner = 0.0_dp
        od = 0.0_dp
        total_gamma = sum(gamma)
        do j = 1, size(data,1)
            n = sum(data(j,:))
            do k = 1, kk
                tail = 0.0_dp
                do r = n, 1, -1
                    tail = tail + beta_binomial_pmf(r, n, gamma(k), total_gamma-gamma(k))
                    z = gamma(k) + real(r-1,dp)
                    inner(k) = inner(k) + tail / (z*z)
                end do
            end do
            od = od + recip_square_sum(total_gamma, n)
        end do
        f = -od
        do k = 1, kk
            f(k,k) = f(k,k) + inner(k)
        end do
    end function expected_fim

    function dirmult_loglik(data, gamma) result(value)
        integer, intent(in) :: data(:,:)
        real(dp), intent(in) :: gamma(:)
        real(dp) :: value, total_gamma
        integer :: i, j, n, k

        value = 0.0_dp
        k = size(gamma)
        if (size(data,2) /= k .or. any(gamma <= 0.0_dp)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        total_gamma = sum(gamma)
        do j = 1, size(data,1)
            n = sum(data(j,:))
            value = value - rising_log_sum(total_gamma, n)
            do i = 1, k
                value = value + rising_log_sum(gamma(i), data(j,i))
            end do
        end do
    end function dirmult_loglik

    function multinomial_loglik(data) result(value)
        integer, intent(in) :: data(:,:)
        real(dp) :: value
        integer, allocatable :: x(:,:)
        real(dp), allocatable :: p(:)
        integer :: j, k

        call clean_counts(data, .false., .true., x)
        if (size(x,2) == 0 .or. sum(x) <= 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        allocate(p(size(x,2)))
        do k = 1, size(x,2)
            p(k) = real(sum(x(:,k)),dp) / real(sum(x),dp)
        end do
        value = 0.0_dp
        do j = 1, size(x,1)
            do k = 1, size(x,2)
                if (x(j,k) > 0) value = value + real(x(j,k),dp) * log(p(k))
            end do
        end do
    end function multinomial_loglik

    function theta_fim(gamma, fim_gamma) result(f)
        real(dp), intent(in) :: gamma(:)
        real(dp), intent(in) :: fim_gamma(:,:)
        real(dp) :: f(size(gamma),size(gamma))
        real(dp) :: d(size(gamma),size(gamma)), pi(size(gamma))
        real(dp) :: theta, gsum
        integer :: i, k

        k = size(gamma)
        f = 0.0_dp
        if (size(fim_gamma,1) /= k .or. size(fim_gamma,2) /= k) return
        gsum = sum(gamma)
        theta = 1.0_dp / (gsum + 1.0_dp)
        pi = gamma / gsum
        d = 0.0_dp
        do i = 1, k
            d(i,i) = (1.0_dp-theta) / theta
        end do
        d(k,:) = -(1.0_dp-theta) / theta
        d(:,k) = -pi / (theta*theta)
        f = matmul(transpose(d), matmul(fim_gamma, d))
    end function theta_fim

    subroutine weir_mom(data, result)
        integer, intent(in) :: data(:,:)
        type(mom_result_type), intent(out) :: result
        integer, allocatable :: x(:,:)
        real(dp), allocatable :: p(:), sn(:), rowp(:)
        real(dp) :: msp, msg, nc, total
        integer :: i, j, jj, k

        result = mom_result_type()
        call clean_counts(data, .true., .true., x)
        jj = size(x,1)
        k = size(x,2)
        if (jj <= 1 .or. k == 0 .or. sum(x) <= jj) then
            result%theta = ieee_value(0.0_dp, ieee_quiet_nan)
            result%se = result%theta
            result%info = 1
            return
        end if
        allocate(p(k), sn(jj), rowp(k))
        total = real(sum(x),dp)
        do i = 1, k
            p(i) = real(sum(x(:,i)),dp) / total
        end do
        do j = 1, jj
            sn(j) = real(sum(x(j,:)),dp)
        end do

        msp = 0.0_dp
        msg = 0.0_dp
        do j = 1, jj
            rowp = real(x(j,:),dp) / sn(j)
            msp = msp + sum((rowp-p)**2) * sn(j)
            msg = msg + sum(rowp*(1.0_dp-rowp)) * sn(j)
        end do
        msp = msp / real(jj-1,dp)
        msg = msg / (total-real(jj,dp))
        nc = (sum(sn)-sum(sn*sn)/sum(sn)) / real(jj-1,dp)
        if (abs(msp+(nc-1.0_dp)*msg) <= tiny(1.0_dp)) then
            result%theta = ieee_value(0.0_dp, ieee_quiet_nan)
            result%se = result%theta
            result%info = 2
            return
        end if
        result%theta = (msp-msg) / (msp+(nc-1.0_dp)*msg)
        result%se = sqrt(2.0_dp*(1.0_dp-result%theta)**2/real(jj-1,dp) * &
                         ((1.0_dp+(nc-1.0_dp)*result%theta)/nc)**2)
    end subroutine weir_mom

    subroutine fit_dirmult(data, fit, init, initscalar, epsilon, trace, mode, maxit)
        integer, intent(in) :: data(:,:)
        type(dirmult_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: init(:), initscalar, epsilon
        logical, intent(in), optional :: trace
        character(len=*), intent(in), optional :: mode
        integer, intent(in), optional :: maxit
        integer, allocatable :: x(:,:)
        real(dp), allocatable :: gamma(:), fim(:,:), score(:), delta(:)
        type(mom_result_type) :: mr
        real(dp) :: scalar, eps, lik1, lik2
        integer :: k, i, ite, lim, stat
        logical :: tr, stop_after
        character(len=3) :: md

        fit = dirmult_fit_type()
        eps = 1.0e-4_dp
        if (present(epsilon)) eps = epsilon
        tr = .true.
        if (present(trace)) tr = trace
        lim = 1000
        if (present(maxit)) lim = maxit
        md = 'obs'
        if (present(mode)) then
            if (trim(mode) == 'exp') md = 'exp'
            if (trim(mode) == 'obs') md = 'obs'
        end if

        call clean_counts(data, .true., .true., x)
        k = size(x,2)
        if (size(x,1) == 0 .or. k == 0 .or. sum(x) <= 0) then
            fit%info = 1
            return
        end if
        allocate(gamma(k), fim(k,k), score(k), delta(k))
        if (present(initscalar)) then
            scalar = initscalar
        else
            call weir_mom(x, mr)
            scalar = mr%theta
            if (.not. (scalar > 0.0_dp)) scalar = 0.005_dp
            scalar = (1.0_dp-scalar) / scalar
        end if
        if (present(init)) then
            if (size(init) /= k) then
                fit%info = 2
                return
            end if
            gamma = init
        else
            do i = 1, k
                gamma(i) = real(sum(x(:,i)),dp) / real(sum(x),dp) * scalar
            end do
        end if
        if (any(gamma <= 0.0_dp)) then
            fit%info = 3
            return
        end if

        lik1 = 0.0_dp
        lik2 = eps * 10.0_dp
        fit%converged = .false.
        do ite = 1, lim
            stop_after = abs(lik2-lik1) < eps
            if (md == 'exp') then
                fim = expected_fim(x, gamma)
            else
                fim = -observed_fim(x, gamma)
            end if
            lik1 = dirmult_loglik(x, gamma)
            score = score_function(x, gamma)
            call solve_linear(fim, score, delta, stat)
            if (stat /= 0) then
                fit%info = 10 + stat
                exit
            end if
            gamma = gamma + delta
            where (gamma < 0.0_dp) gamma = 0.01_dp
            if (tr) write(*,'(a,i0,a,es24.15)') 'Iteration ', ite, &
                ': Log-likelihood value: ', lik1
            lik2 = dirmult_loglik(x, gamma)
            if (stop_after) then
                fit%converged = .true.
                exit
            end if
        end do

        fit%iterations = min(ite,lim)
        fit%loglik = lik1
        allocate(fit%gamma(k), fit%pi(k))
        fit%gamma = gamma
        scalar = sum(gamma)
        fit%pi = gamma / scalar
        fit%theta = 1.0_dp / (scalar+1.0_dp)
        if (.not. fit%converged .and. fit%info == 0) fit%info = 4
    end subroutine fit_dirmult

    subroutine summarize_dirmult(data, fit, summary, expected)
        integer, intent(in) :: data(:,:)
        type(dirmult_fit_type), intent(in) :: fit
        type(dirmult_summary_type), intent(out) :: summary
        logical, intent(in), optional :: expected
        integer, allocatable :: x(:,:)
        real(dp), allocatable :: fim(:,:), tf(:,:), invf(:,:), p(:), sn(:), rowp(:)
        type(mom_result_type) :: mr
        logical :: use_exp
        real(dp) :: var_last, total
        integer :: i, j, k, jj, stat

        summary = dirmult_summary_type()
        use_exp = .false.
        if (present(expected)) use_exp = expected
        call clean_counts(data, .true., .true., x)
        k = size(x,2)
        jj = size(x,1)
        if (.not. allocated(fit%gamma) .or. size(fit%gamma) /= k .or. k < 2) then
            summary%info = 1
            return
        end if
        allocate(fim(k,k), tf(k,k), invf(k,k), p(k), sn(jj), rowp(k))
        if (use_exp) then
            fim = expected_fim(x, fit%gamma)
        else
            fim = -observed_fim(x, fit%gamma)
        end if
        tf = theta_fim(fit%gamma, fim)
        call invert_matrix(tf, invf, stat)
        if (stat /= 0) then
            summary%info = 10 + stat
            return
        end if

        allocate(summary%mle(k+1), summary%se_mle(k+1))
        allocate(summary%mom(k+1), summary%se_mom(k+1))
        summary%mle(1:k) = fit%pi
        summary%mle(k+1) = fit%theta
        do i = 1, k-1
            summary%se_mle(i) = sqrt(max(0.0_dp, invf(i,i)))
        end do
        var_last = sum(invf(1:k-1,1:k-1))
        summary%se_mle(k) = sqrt(max(0.0_dp, var_last))
        summary%se_mle(k+1) = sqrt(max(0.0_dp, invf(k,k)))

        total = real(sum(x),dp)
        do i = 1, k
            p(i) = real(sum(x(:,i)),dp) / total
        end do
        summary%mom(1:k) = p
        call weir_mom(x, mr)
        summary%mom(k+1) = mr%theta
        do j = 1, jj
            sn(j) = real(sum(x(j,:)),dp)
        end do
        do i = 1, k
            summary%se_mom(i) = 0.0_dp
            do j = 1, jj
                rowp = real(x(j,:),dp) / sn(j)
                summary%se_mom(i) = summary%se_mom(i) + (rowp(i)-p(i))**2
            end do
            summary%se_mom(i) = sqrt(summary%se_mom(i)/real(jj-1,dp))
        end do
        summary%se_mom(k+1) = mr%se
    end subroutine summarize_dirmult

end module dirmult_core
