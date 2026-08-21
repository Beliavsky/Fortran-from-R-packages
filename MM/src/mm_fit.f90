! Computational translation of the R package MM 1.7-0.
! Upstream license: GPL-2. This translation is GPL-2.
module mm_fit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use mm_kinds, only : dp
    use mm_types, only : paras_type, mb_type, gunter_mb_type, glm_fit_type, mm_fit_type
    use mm_parameters, only : paras, paras_from_values, paras_dimension, p
    use mm_core, only : mm_loglik, suffstats, mm_support
    use mm_tables, only : gunter, gunter_mb
    use partitions, only : compositions
    use quadform_linalg, only : solve_linear
    implicit none
    private

    public :: lindsey, lindsey_fit, lindsey_mb
    public :: optimizer, optimizer_allsamesum, optimizer_differsums
    public :: poisson_glm_fit

contains

    function lindsey(obs, n) result(par)
        integer, intent(in) :: obs(:,:)
        integer, intent(in), optional :: n(:)
        type(paras_type) :: par
        type(glm_fit_type) :: fit
        real(dp), allocatable :: wt(:)

        if (present(n)) then
            if (size(n) /= size(obs, 1) .or. any(n < 0)) error stop "lindsey: invalid n"
            wt = real(n, dp)
            call lindsey_fit(obs, par, fit, wt)
        else
            call lindsey_fit(obs, par, fit)
        end if
    end function lindsey

    subroutine lindsey_fit(obs, par, fit, n)
        integer, intent(in) :: obs(:,:)
        type(paras_type), intent(out) :: par
        type(glm_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: n(:)
        integer, allocatable :: comp(:,:)
        real(dp), allocatable :: d(:), x(:,:), offset(:), wt(:), ecoef(:), vals(:)
        integer :: nr, k, y_total, ns, q, r, s, i, j, col

        nr = size(obs, 1)
        k = size(obs, 2)
        if (nr < 1 .or. k < 2) error stop "lindsey_fit: invalid observation matrix"
        if (any(obs < 0)) error stop "lindsey_fit: negative observations"
        if (.not. all(sum(obs, dim=2) == sum(obs(1, :)))) error stop "lindsey_fit: row sums differ"
        allocate(wt(nr))
        if (present(n)) then
            if (size(n) /= nr .or. any(n < 0.0_dp)) error stop "lindsey_fit: invalid n"
            wt = n
        else
            wt = 1.0_dp
        end if

        y_total = sum(obs(1, :))
        if (y_total == 0) error stop "lindsey_fit: total count must be positive"
        comp = compositions(y_total, k)
        ns = size(comp, 2)
        q = k + k * (k - 1) / 2
        allocate(d(ns), x(ns, q), offset(ns))
        d = 0.0_dp
        do r = 1, nr
            do s = 1, ns
                if (all(obs(r, :) == comp(:, s))) then
                    d(s) = d(s) + wt(r)
                    exit
                end if
            end do
        end do

        do s = 1, ns
            x(s, 1:k) = real(comp(:, s), dp)
            col = k + 1
            do j = 2, k
                do i = 1, j - 1
                    x(s, col) = real(comp(i, s) * comp(j, s), dp)
                    col = col + 1
                end do
            end do
            offset(s) = 0.0_dp
            do i = 1, k
                offset(s) = offset(s) - log_gamma(real(comp(i, s) + 1, dp))
            end do
        end do

        call poisson_glm_fit(x, d, offset, fit)
        if (.not. allocated(fit%coefficients)) error stop "lindsey_fit: GLM fit failed"
        allocate(ecoef(q), vals(q - 1))
        ecoef = exp(fit%coefficients)
        if (sum(ecoef(1:k)) <= 0.0_dp) error stop "lindsey_fit: invalid main effects"
        ecoef(1:k) = ecoef(1:k) / sum(ecoef(1:k))
        vals(1:k - 1) = ecoef(1:k - 1)
        vals(k:q - 1) = ecoef(k + 1:q)
        par = paras_from_values(vals)
    end subroutine lindsey_fit

    subroutine lindsey_mb(a, fit)
        type(mb_type), intent(in) :: a
        type(glm_fit_type), intent(out) :: fit
        type(gunter_mb_type) :: g
        real(dp), allocatable :: x(:,:), d(:), offset(:)
        integer :: ns, x1, x2, m1, m2, s

        if (size(a%m) /= 2) error stop "lindsey_mb: only bivariate case is implemented upstream"
        g = gunter_mb(a)
        ns = size(g%d)
        m1 = a%m(1)
        m2 = a%m(2)
        allocate(x(ns, 6), d(ns), offset(ns))
        d = real(g%d, dp)
        do s = 1, ns
            x1 = g%tbl(s, 1)
            x2 = g%tbl(s, 2)
            x(s, 1) = 1.0_dp
            x(s, 2) = real(x1, dp)
            x(s, 3) = real(x2, dp)
            x(s, 4) = real(x1 * (m1 - x1), dp)
            x(s, 5) = real(x2 * (m2 - x2), dp)
            x(s, 6) = real(x1 * x2, dp)
            offset(s) = log_choose(m1, x1) + log_choose(m2, x2)
        end do
        call poisson_glm_fit(x, d, offset, fit)
    end subroutine lindsey_mb

    subroutine poisson_glm_fit(x, response, offset, fit, max_iter, tol)
        real(dp), intent(in) :: x(:,:), response(:), offset(:)
        type(glm_fit_type), intent(out) :: fit
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: beta(:), eta(:), mu(:), grad(:), h(:,:), wx(:,:)
        real(dp), allocatable :: rhs(:,:), sol(:,:), trial(:), fitted(:)
        real(dp) :: eps, old_ll, new_ll, step, ridge
        integer :: n, q, it, mit, info, j

        n = size(x, 1)
        q = size(x, 2)
        if (size(response) /= n .or. size(offset) /= n) error stop "poisson_glm_fit: incompatible arrays"
        if (any(response < 0.0_dp)) error stop "poisson_glm_fit: negative response"
        mit = 100
        if (present(max_iter)) mit = max_iter
        eps = 1.0e-10_dp
        if (present(tol)) eps = tol
        allocate(beta(q), eta(n), mu(n), grad(q), h(q, q), wx(n, q), rhs(q, 1), trial(q))
        beta = 0.0_dp
        old_ll = poisson_loglik(x, response, offset, beta)

        do it = 1, mit
            eta = offset + matmul(x, beta)
            eta = max(-700.0_dp, min(700.0_dp, eta))
            mu = exp(eta)
            grad = matmul(transpose(x), response - mu)
            wx = x * spread(mu, dim=2, ncopies=q)
            h = matmul(transpose(x), wx)
            rhs(:, 1) = grad
            sol = solve_linear(h, rhs, info)
            if (info /= 0) then
                ridge = 1.0e-10_dp * max(1.0_dp, maxval(abs(h)))
                do j = 1, q
                    h(j, j) = h(j, j) + ridge
                end do
                sol = solve_linear(h, rhs, info)
            end if
            if (info /= 0) exit

            step = 1.0_dp
            do
                trial = beta + step * sol(:, 1)
                new_ll = poisson_loglik(x, response, offset, trial)
                if (new_ll >= old_ll .or. step <= 1.0e-8_dp) exit
                step = 0.5_dp * step
            end do
            beta = trial
            if (maxval(abs(step * sol(:, 1))) <= eps * (1.0_dp + maxval(abs(beta)))) then
                fit%converged = .true.
                old_ll = new_ll
                exit
            end if
            old_ll = new_ll
        end do

        fit%iterations = min(it, mit)
        fit%loglik = old_ll
        fit%coefficients = beta
        allocate(fitted(n))
        fitted = exp(max(-700.0_dp, min(700.0_dp, offset + matmul(x, beta))))
        fit%fitted = fitted
    end subroutine poisson_glm_fit

    subroutine optimizer(y, fit, n, start, max_iter, tol, method)
        integer, intent(in) :: y(:,:)
        type(mm_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: n(:)
        type(paras_type), intent(in), optional :: start
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        character(len=*), intent(in), optional :: method
        integer, allocatable :: totals(:)

        totals = sum(y, dim=2)
        if (all(totals == totals(1))) then
            call optimizer_allsamesum(y, fit, n, start, max_iter, tol, method)
        else
            call optimizer_differsums(y, fit, n, start, max_iter, tol, method)
        end if
    end subroutine optimizer

    subroutine optimizer_allsamesum(y, fit, n, start, max_iter, tol, method)
        integer, intent(in) :: y(:,:)
        type(mm_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: n(:)
        type(paras_type), intent(in), optional :: start
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        character(len=*), intent(in), optional :: method
        type(paras_type) :: p0
        type(glm_fit_type) :: gl

        if (.not. all(sum(y, dim=2) == sum(y(1, :)))) error stop "optimizer_allsamesum: row sums differ"
        if (present(start)) then
            p0 = start
        else
            if (present(n)) then
                call lindsey_fit(y, p0, gl, n)
            else
                call lindsey_fit(y, p0, gl)
            end if
        end if
        if (present(method)) then
            if (is_nelder_method(method)) then
                call nelder_mead_mm(y, p0, fit, n, max_iter, tol)
            else
                call bfgs_mm(y, p0, fit, n, max_iter, tol)
            end if
        else
            call bfgs_mm(y, p0, fit, n, max_iter, tol)
        end if
    end subroutine optimizer_allsamesum

    subroutine optimizer_differsums(y, fit, n, start, max_iter, tol, method)
        integer, intent(in) :: y(:,:)
        type(mm_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: n(:)
        type(paras_type), intent(in), optional :: start
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        character(len=*), intent(in), optional :: method
        type(paras_type) :: p0

        if (present(start)) then
            p0 = start
        else
            p0 = paras(size(y, 2))
        end if
        if (present(method)) then
            if (is_nelder_method(method)) then
                call nelder_mead_mm(y, p0, fit, n, max_iter, tol)
            else
                call bfgs_mm(y, p0, fit, n, max_iter, tol)
            end if
        else
            call bfgs_mm(y, p0, fit, n, max_iter, tol)
        end if
    end subroutine optimizer_differsums

    subroutine nelder_mead_mm(y, p0, fit, n, max_iter, tol)
        integer, intent(in) :: y(:,:)
        type(paras_type), intent(in) :: p0
        type(mm_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: n(:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: simp(:,:), fv(:), centroid(:), xr(:), xe(:), xc(:), x0(:)
        real(dp) :: fr, fe, fc, eps, step, spreadx
        integer :: d, i, it, mit
        logical :: conv

        d = size(p0%vals)
        mit = 1000
        if (present(max_iter)) mit = max_iter
        eps = 1.0e-8_dp
        if (present(tol)) eps = tol
        allocate(simp(d, d + 1), fv(d + 1), centroid(d), xr(d), xe(d), xc(d), x0(d))
        x0 = log(max(p0%vals, tiny(1.0_dp)))
        simp(:, 1) = x0
        step = 0.08_dp
        do i = 1, d
            simp(:, i + 1) = x0
            simp(i, i + 1) = simp(i, i + 1) + step
        end do
        do i = 1, d + 1
            fv(i) = objective(simp(:, i))
        end do

        conv = .false.
        do it = 1, mit
            call sort_simplex(simp, fv)
            spreadx = maxval(abs(simp(:, 2:d + 1) - spread(simp(:, 1), 2, d)))
            if (maxval(abs(fv - fv(1))) <= eps * (1.0_dp + abs(fv(1))) .and. &
                spreadx <= sqrt(eps) * (1.0_dp + maxval(abs(simp(:, 1))))) then
                conv = .true.
                exit
            end if
            centroid = sum(simp(:, 1:d), dim=2) / real(d, dp)
            xr = centroid + (centroid - simp(:, d + 1))
            fr = objective(xr)
            if (fr < fv(1)) then
                xe = centroid + 2.0_dp * (xr - centroid)
                fe = objective(xe)
                if (fe < fr) then
                    simp(:, d + 1) = xe
                    fv(d + 1) = fe
                else
                    simp(:, d + 1) = xr
                    fv(d + 1) = fr
                end if
            else if (fr < fv(d)) then
                simp(:, d + 1) = xr
                fv(d + 1) = fr
            else
                if (fr < fv(d + 1)) then
                    xc = centroid + 0.5_dp * (xr - centroid)
                else
                    xc = centroid + 0.5_dp * (simp(:, d + 1) - centroid)
                end if
                fc = objective(xc)
                if (fc < min(fr, fv(d + 1))) then
                    simp(:, d + 1) = xc
                    fv(d + 1) = fc
                else
                    do i = 2, d + 1
                        simp(:, i) = simp(:, 1) + 0.5_dp * (simp(:, i) - simp(:, 1))
                        fv(i) = objective(simp(:, i))
                    end do
                end if
            end if
        end do
        call sort_simplex(simp, fv)
        fit%parameters = paras_from_values(exp(simp(:, 1)))
        fit%loglik = -fv(1)
        fit%iterations = min(it, mit)
        fit%converged = conv

    contains

        real(dp) function objective(x) result(val)
            real(dp), intent(in) :: x(:)
            type(paras_type) :: par
            real(dp), allocatable :: vals(:), prob(:)
            real(dp) :: ll

            if (any(x > 700.0_dp) .or. any(x < -700.0_dp)) then
                val = 1.0e200_dp
                return
            end if
            vals = exp(x)
            par = paras_from_values(vals)
            prob = p(par)
            if (any(prob <= 0.0_dp) .or. any(prob >= 1.0_dp)) then
                val = 1.0e200_dp
                return
            end if
            if (present(n)) then
                ll = mm_loglik(y, par, n)
            else
                ll = mm_loglik(y, par)
            end if
            if (ieee_is_nan(ll) .or. abs(ll) >= huge(1.0_dp) / 10.0_dp) then
                val = 1.0e200_dp
            else
                val = -ll
            end if
        end function objective

    end subroutine nelder_mead_mm

    subroutine bfgs_mm(y, p0, fit, n, max_iter, tol)
        integer, intent(in) :: y(:,:)
        type(paras_type), intent(in) :: p0
        type(mm_fit_type), intent(out) :: fit
        real(dp), intent(in), optional :: n(:)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: x(:), xnew(:), g(:), gnew(:), direction(:), hmat(:,:)
        real(dp), allocatable :: svec(:), yvec(:), ident(:,:), left(:,:), right(:,:)
        real(dp) :: f, fnew, alpha, dg, ys, rho, eps
        integer :: d, i, it, mit
        logical :: conv

        d = size(p0%vals)
        mit = 250
        if (present(max_iter)) mit = max_iter
        eps = 1.0e-8_dp
        if (present(tol)) eps = tol
        allocate(x(d), xnew(d), g(d), gnew(d), direction(d), hmat(d,d))
        allocate(svec(d), yvec(d), ident(d,d), left(d,d), right(d,d))
        x = log(max(p0%vals, tiny(1.0_dp)))
        ident = 0.0_dp
        do i = 1, d
            ident(i,i) = 1.0_dp
        end do
        hmat = ident
        f = objective(x)
        call numerical_gradient(x, f, g)
        conv = .false.

        do it = 1, mit
            if (maxval(abs(g)) <= eps * (1.0_dp + abs(f))) then
                conv = .true.
                exit
            end if
            direction = -matmul(hmat, g)
            dg = dot_product(g, direction)
            if (dg >= -sqrt(epsilon(1.0_dp)) * max(1.0_dp, dot_product(g,g))) then
                direction = -g
                hmat = ident
                dg = -dot_product(g,g)
            end if

            alpha = 1.0_dp
            do
                xnew = x + alpha * direction
                fnew = objective(xnew)
                if (fnew <= f + 1.0e-4_dp * alpha * dg) exit
                alpha = 0.5_dp * alpha
                if (alpha < 1.0e-8_dp) exit
            end do
            if (alpha < 1.0e-8_dp .and. fnew >= f) exit

            call numerical_gradient(xnew, fnew, gnew)
            svec = xnew - x
            yvec = gnew - g
            ys = dot_product(yvec, svec)
            if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, &
                sqrt(dot_product(yvec,yvec) * dot_product(svec,svec)))) then
                rho = 1.0_dp / ys
                left = ident - rho * outer_product(svec, yvec)
                right = ident - rho * outer_product(yvec, svec)
                hmat = matmul(left, matmul(hmat, right)) + rho * outer_product(svec, svec)
            else
                hmat = ident
            end if

            x = xnew
            g = gnew
            f = fnew
            if (maxval(abs(svec)) <= eps * (1.0_dp + maxval(abs(x)))) then
                conv = .true.
                exit
            end if
        end do

        fit%parameters = paras_from_values(exp(x))
        fit%loglik = -f
        fit%iterations = min(it, mit)
        fit%converged = conv

    contains

        real(dp) function objective(z) result(val)
            real(dp), intent(in) :: z(:)
            type(paras_type) :: par
            real(dp), allocatable :: vals(:), prob(:)
            real(dp) :: ll

            if (any(z > 700.0_dp) .or. any(z < -700.0_dp)) then
                val = 1.0e200_dp
                return
            end if
            vals = exp(z)
            par = paras_from_values(vals)
            prob = p(par)
            if (any(prob <= 0.0_dp) .or. any(prob >= 1.0_dp)) then
                val = 1.0e200_dp
                return
            end if
            if (present(n)) then
                ll = mm_loglik(y, par, n)
            else
                ll = mm_loglik(y, par)
            end if
            if (ieee_is_nan(ll) .or. abs(ll) >= huge(1.0_dp) / 10.0_dp) then
                val = 1.0e200_dp
            else
                val = -ll
            end if
        end function objective

        subroutine numerical_gradient(z, fz, grad)
            real(dp), intent(in) :: z(:), fz
            real(dp), intent(out) :: grad(:)
            real(dp), allocatable :: zp(:), zm(:)
            real(dp) :: hp, fp, fm
            integer :: j

            allocate(zp(size(z)), zm(size(z)))
            do j = 1, size(z)
                hp = 1.0e-5_dp * (1.0_dp + abs(z(j)))
                zp = z
                zm = z
                zp(j) = zp(j) + hp
                zm(j) = zm(j) - hp
                fp = objective(zp)
                fm = objective(zm)
                if (fp < 1.0e190_dp .and. fm < 1.0e190_dp) then
                    grad(j) = (fp - fm) / (2.0_dp * hp)
                else if (fp < 1.0e190_dp) then
                    grad(j) = (fp - fz) / hp
                else if (fm < 1.0e190_dp) then
                    grad(j) = (fz - fm) / hp
                else
                    grad(j) = 0.0_dp
                end if
            end do
        end subroutine numerical_gradient

    end subroutine bfgs_mm

    pure function outer_product(a, b) result(out)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: out(size(a), size(b))
        out = spread(a, 2, size(b)) * spread(b, 1, size(a))
    end function outer_product

    pure logical function is_nelder_method(method) result(ans)
        character(len=*), intent(in) :: method
        character(len=:), allocatable :: m
        m = trim(adjustl(method))
        ans = (m == 'Nelder' .or. m == 'nelder' .or. m == 'Nelder-Mead' .or. &
            m == 'nelder-mead' .or. m == 'NELDER')
    end function is_nelder_method

    subroutine sort_simplex(simp, fv)
        real(dp), intent(inout) :: simp(:,:), fv(:)
        real(dp), allocatable :: tmp(:)
        real(dp) :: tf
        integer :: i, j, best

        allocate(tmp(size(simp, 1)))
        do i = 1, size(fv) - 1
            best = i
            do j = i + 1, size(fv)
                if (fv(j) < fv(best)) best = j
            end do
            if (best /= i) then
                tf = fv(i)
                fv(i) = fv(best)
                fv(best) = tf
                tmp = simp(:, i)
                simp(:, i) = simp(:, best)
                simp(:, best) = tmp
            end if
        end do
    end subroutine sort_simplex

    pure real(dp) function poisson_loglik(x, response, offset, beta) result(ans)
        real(dp), intent(in) :: x(:,:), response(:), offset(:), beta(:)
        real(dp), allocatable :: eta(:)

        eta = offset + matmul(x, beta)
        if (any(eta > 700.0_dp)) then
            ans = -huge(1.0_dp)
            return
        end if
        eta = max(eta, -700.0_dp)
        ans = sum(response * eta - exp(eta))
    end function poisson_loglik

    pure real(dp) function log_choose(n, k) result(ans)
        integer, intent(in) :: n, k
        if (k < 0 .or. k > n) then
            ans = -huge(1.0_dp)
        else
            ans = log_gamma(real(n + 1, dp)) - log_gamma(real(k + 1, dp)) - &
                log_gamma(real(n - k + 1, dp))
        end if
    end function log_choose

end module mm_fit
