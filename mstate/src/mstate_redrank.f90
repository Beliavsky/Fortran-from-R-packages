module mstate_redrank
    use mstate_kinds, only : dp
    use mstate_types, only : redrank_result
    use mstate_cox, only : coxph_fit_stratified_counting
    use survival_types, only : coxph_result
    implicit none
    private
    public :: redrank_fit

contains

    subroutine redrank_fit(start, stop, status, transition, z, rank_r, result, full_x, strata, &
                           gamma_start, method, eps, maxiter, clock, info)
        real(dp), intent(in) :: start(:), stop(:), z(:, :)
        integer, intent(in) :: status(:), transition(:), rank_r
        type(redrank_result), intent(out) :: result
        real(dp), intent(in), optional :: full_x(:, :), gamma_start(:, :)
        integer, intent(in), optional :: strata(:), maxiter
        character(len=*), intent(in), optional :: method, clock
        real(dp), intent(in), optional :: eps
        integer, intent(out), optional :: info
        integer :: n, p, p2, ktrans, r, k, i, j, iter, max_it, n1, n2
        integer, allocatable :: strat(:)
        real(dp), allocatable :: gamma(:, :), alpha(:, :), beta2(:), alphax(:, :)
        real(dp), allocatable :: x1(:, :), x2(:, :), fstart(:), fstop(:), bmat(:, :)
        type(coxph_result) :: fit1, fit2
        real(dp) :: tolerance, previous_ll, delta
        character(len=12) :: tie_method, clock_type
        logical :: have_previous

        if (present(info)) info = 0
        n = size(stop)
        p = size(z, 2)
        p2 = 0
        if (present(full_x)) p2 = size(full_x, 2)
        if (size(start) /= n .or. size(status) /= n .or. size(transition) /= n .or. size(z,1) /= n) then
            if (present(info)) info = 1
            return
        end if
        if (present(full_x)) then
            if (size(full_x,1) /= n) then
                if (present(info)) info = 2
                return
            end if
        end if
        if (p == 0) then
            if (present(info)) info = 3
            return
        end if
        ktrans = maxval(transition)
        if (minval(transition) < 1 .or. rank_r < 1 .or. rank_r > min(p,ktrans)) then
            if (present(info)) info = 4
            return
        end if
        allocate(strat(n))
        if (present(strata)) then
            if (size(strata) /= n) then
                if (present(info)) info = 5
                return
            end if
            strat = strata
        else
            strat = transition
        end if

        tie_method = 'breslow'
        if (present(method)) tie_method = adjustl(method)
        clock_type = 'forward'
        if (present(clock)) clock_type = adjustl(clock)
        tolerance = 1.0e-5_dp
        if (present(eps)) tolerance = eps
        max_it = 100
        if (present(maxiter)) max_it = maxiter

        allocate(fstart(n), fstop(n))
        if (index(clock_type, 'reset') == 1) then
            fstart = 0.0_dp
            fstop = stop - start
        else
            fstart = start
            fstop = stop
        end if

        allocate(gamma(rank_r,ktrans), alpha(p,rank_r), beta2(p2), alphax(n,rank_r))
        if (present(gamma_start)) then
            if (size(gamma_start,1) /= rank_r .or. size(gamma_start,2) /= ktrans) then
                if (present(info)) info = 6
                return
            end if
            gamma = gamma_start
        else
            call fill_normal(gamma)
        end if
        alpha = 0.0_dp
        beta2 = 0.0_dp
        alphax = 0.0_dp

        n1 = p*rank_r + p2
        n2 = rank_r*ktrans + p2
        allocate(x1(n,n1), x2(n,n2))
        previous_ll = 0.0_dp
        have_previous = .false.
        result%converged = .false.

        do iter = 1, max_it
            x1 = 0.0_dp
            do i = 1, n
                k = transition(i)
                do r = 1, rank_r
                    do j = 1, p
                        x1(i,(r-1)*p+j) = gamma(r,k)*z(i,j)
                    end do
                end do
            end do
            if (p2 > 0) x1(:,p*rank_r+1:n1) = full_x
            call coxph_fit_stratified_counting(fstart, fstop, status, x1, strat, fit1, tie_method, &
                                               maxiter=50, eps=1.0e-9_dp)
            if (size(fit1%coef) < p*rank_r) then
                if (present(info)) info = 7
                return
            end if
            do r = 1, rank_r
                alpha(:,r) = fit1%coef((r-1)*p+1:r*p)
            end do
            if (p2 > 0) beta2 = fit1%coef(p*rank_r+1:n1)
            alphax = matmul(z, alpha)

            x2 = 0.0_dp
            do i = 1, n
                k = transition(i)
                do r = 1, rank_r
                    x2(i,(r-1)*ktrans+k) = alphax(i,r)
                end do
            end do
            if (p2 > 0) x2(:,rank_r*ktrans+1:n2) = full_x
            call coxph_fit_stratified_counting(fstart, fstop, status, x2, strat, fit2, tie_method, &
                                               maxiter=50, eps=1.0e-9_dp)
            if (size(fit2%coef) < rank_r*ktrans) then
                if (present(info)) info = 8
                return
            end if
            do r = 1, rank_r
                gamma(r,:) = fit2%coef((r-1)*ktrans+1:r*ktrans)
            end do
            if (p2 > 0) beta2 = fit2%coef(rank_r*ktrans+1:n2)

            if (have_previous) then
                delta = abs(fit1%loglik - previous_ll)
                if (delta <= tolerance) then
                    result%converged = .true.
                    exit
                end if
            end if
            previous_ll = fit1%loglik
            have_previous = .true.
        end do

        allocate(bmat(p,ktrans))
        bmat = matmul(alpha, gamma)
        call normalize_reduced_rank(bmat, rank_r, alpha, gamma)
        alphax = matmul(z, alpha)

        result%rank = rank_r
        result%niter = min(iter, max_it)
        result%df = rank_r*(p + ktrans - rank_r)
        result%loglik = fit1%loglik
        allocate(result%alpha(p,rank_r), result%gamma(rank_r,ktrans), result%beta(p,ktrans))
        allocate(result%beta2(p2), result%alphax(n,rank_r))
        result%alpha = alpha
        result%gamma = gamma
        result%beta = matmul(alpha,gamma)
        result%beta2 = beta2
        result%alphax = alphax
    end subroutine redrank_fit

    subroutine normalize_reduced_rank(b, rank_r, alpha, gamma)
        real(dp), intent(in) :: b(:, :)
        integer, intent(in) :: rank_r
        real(dp), intent(out) :: alpha(:, :), gamma(:, :)
        real(dp), allocatable :: a(:, :), v(:, :), sval(:), tmpcol(:)
        real(dp) :: app, aqq, apq, tau, t, c, s, scale, tol
        integer :: m, n, i, j, sweep, maxs, r, imax

        m = size(b,1); n = size(b,2)
        allocate(a(m,n), v(n,n), sval(n), tmpcol(max(m,n)))
        a = b; v = 0.0_dp
        do i = 1, n
            v(i,i) = 1.0_dp
        end do
        tol = 64.0_dp*epsilon(1.0_dp)
        maxs = max(20, 10*n*n)
        do sweep = 1, maxs
            scale = 0.0_dp
            do i = 1, n-1
                do j = i+1, n
                    app = dot_product(a(:,i),a(:,i))
                    aqq = dot_product(a(:,j),a(:,j))
                    apq = dot_product(a(:,i),a(:,j))
                    scale = max(scale, abs(apq))
                    if (abs(apq) <= tol*sqrt(max(app*aqq,tiny(1.0_dp)))) cycle
                    tau = (aqq-app)/(2.0_dp*apq)
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp/(tau + sqrt(1.0_dp+tau*tau))
                    else
                        t = -1.0_dp/(-tau + sqrt(1.0_dp+tau*tau))
                    end if
                    c = 1.0_dp/sqrt(1.0_dp+t*t); s = c*t
                    tmpcol(1:m) = c*a(:,i) - s*a(:,j)
                    a(:,j) = s*a(:,i) + c*a(:,j)
                    a(:,i) = tmpcol(1:m)
                    tmpcol(1:n) = c*v(:,i) - s*v(:,j)
                    v(:,j) = s*v(:,i) + c*v(:,j)
                    v(:,i) = tmpcol(1:n)
                end do
            end do
            if (scale <= tol*max(1.0_dp,maxval(abs(a)))) exit
        end do
        do j = 1, n
            sval(j) = sqrt(max(0.0_dp,dot_product(a(:,j),a(:,j))))
        end do
        do i = 1, min(rank_r,n)
            imax = i - 1 + maxloc(sval(i:n),dim=1)
            if (imax /= i) then
                t=sval(i);sval(i)=sval(imax);sval(imax)=t
                tmpcol(1:m)=a(:,i);a(:,i)=a(:,imax);a(:,imax)=tmpcol(1:m)
                tmpcol(1:n)=v(:,i);v(:,i)=v(:,imax);v(:,imax)=tmpcol(1:n)
            end if
        end do
        alpha = 0.0_dp; gamma = 0.0_dp
        do r = 1, rank_r
            if (sval(r) > sqrt(tiny(1.0_dp))) alpha(:,r) = a(:,r)/sval(r)
            gamma(r,:) = sval(r)*v(:,r)
        end do
    end subroutine normalize_reduced_rank

    subroutine fill_normal(x)
        real(dp), intent(out) :: x(:, :)
        integer :: i, j
        real(dp) :: u1, u2, pi
        pi = acos(-1.0_dp)
        do j = 1, size(x,2)
            i = 1
            do while (i <= size(x,1))
                call random_number(u1); call random_number(u2)
                u1 = max(u1,tiny(1.0_dp))
                x(i,j) = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
                if (i+1 <= size(x,1)) x(i+1,j) = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*pi*u2)
                i = i + 2
            end do
        end do
    end subroutine fill_normal

end module mstate_redrank
