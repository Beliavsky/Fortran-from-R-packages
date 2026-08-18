! SPDX-License-Identifier: GPL-3.0-or-later
!
! Native Fortran implementations of the four Cuba-style integration families.
! Cuhre uses deterministic h-adaptive cubature.  Divonne uses adaptive
! stratified quasi-Monte-Carlo partitioning.  Suave combines adaptive
! subdivision with randomized low-discrepancy sampling.  Vegas uses an
! adaptive separable importance grid and inverse-variance iteration combining.
module cuba_mod
    use cubature_kinds, only : dp, i8
    use cubature_types, only : cubature_result, cubature_integrand, cubature_integrand_v, &
        cuhre_options, divonne_options, suave_options, vegas_options, CUBATURE_SUCCESS, &
        CUBATURE_MAXEVAL, CUBATURE_BADARG, ERROR_INDIVIDUAL, errors_converged
    use cubature_rng, only : rng_state, rng_seed, rng_uniform
    use cubature_utils, only : halton_point
    use hcubature_mod, only : hcubature, hcubature_v
    implicit none
    private
    public :: cuhre, cuhre_v, divonne, divonne_v, suave, suave_v, vegas, vegas_v

contains

    subroutine cuhre(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(cuhre_options), intent(in), optional :: options
        type(cuhre_options) :: opt
        real(dp) :: rt, at
        opt = cuhre_options()
        if (present(options)) opt = options
        rt = 1.0e-5_dp; at = 1.0e-12_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        call hcubature(f, lower, upper, ncomp, result, rt, at, opt%max_eval, ERROR_INDIVIDUAL)
    end subroutine cuhre

    subroutine cuhre_v(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand_v) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(cuhre_options), intent(in), optional :: options
        type(cuhre_options) :: opt
        real(dp) :: rt, at
        opt = cuhre_options()
        if (present(options)) opt = options
        rt = 1.0e-5_dp; at = 1.0e-12_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        call hcubature_v(f, lower, upper, ncomp, result, rt, at, opt%max_eval, ERROR_INDIVIDUAL)
    end subroutine cuhre_v

    subroutine divonne(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(divonne_options), intent(in), optional :: options
        type(divonne_options) :: opt
        real(dp) :: rt, at
        integer :: nsample
        opt = divonne_options()
        if (present(options)) opt = options
        rt = 1.0e-5_dp; at = 1.0e-12_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        nsample = max(128, 4 * abs(opt%key1))
        call adaptive_stratified(f=f, lower=lower, upper=upper, ncomp=ncomp, result=result, rel_tol=rt, abs_tol=at, &
            min_eval=opt%min_eval, max_eval=opt%max_eval, nsample=nsample, max_regions=max(64, opt%max_regions), &
            seed=opt%seed, randomized=.true.)
    end subroutine divonne

    subroutine divonne_v(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand_v) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(divonne_options), intent(in), optional :: options
        type(divonne_options) :: opt
        real(dp) :: rt, at
        integer :: nsample
        opt = divonne_options()
        if (present(options)) opt = options
        rt = 1.0e-5_dp; at = 1.0e-12_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        nsample = max(128, 4 * abs(opt%key1))
        call adaptive_stratified(fv=f, lower=lower, upper=upper, ncomp=ncomp, result=result, rel_tol=rt, abs_tol=at, &
            min_eval=opt%min_eval, max_eval=opt%max_eval, nsample=nsample, max_regions=max(64, opt%max_regions), &
            seed=opt%seed, randomized=.true.)
    end subroutine divonne_v

    subroutine suave(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(suave_options), intent(in), optional :: options
        type(suave_options) :: opt
        real(dp) :: rt, at
        opt = suave_options()
        if (present(options)) opt = options
        rt = 1.0e-5_dp; at = 1.0e-12_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        call adaptive_stratified(f=f, lower=lower, upper=upper, ncomp=ncomp, result=result, rel_tol=rt, abs_tol=at, &
            min_eval=opt%min_eval, max_eval=opt%max_eval, nsample=max(opt%nmin, opt%nnew), &
            max_regions=max(64, opt%max_regions), seed=opt%seed, randomized=.true.)
    end subroutine suave

    subroutine suave_v(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand_v) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(suave_options), intent(in), optional :: options
        type(suave_options) :: opt
        real(dp) :: rt, at
        opt = suave_options()
        if (present(options)) opt = options
        rt = 1.0e-5_dp; at = 1.0e-12_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        call adaptive_stratified(fv=f, lower=lower, upper=upper, ncomp=ncomp, result=result, rel_tol=rt, abs_tol=at, &
            min_eval=opt%min_eval, max_eval=opt%max_eval, nsample=max(opt%nmin, opt%nnew), &
            max_regions=max(64, opt%max_regions), seed=opt%seed, randomized=.true.)
    end subroutine suave_v

    subroutine adaptive_stratified(f, fv, lower, upper, ncomp, result, rel_tol, abs_tol, min_eval, max_eval, &
        nsample, max_regions, seed, randomized)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:), rel_tol, abs_tol
        integer, intent(in) :: ncomp, nsample, max_regions, seed
        type(cubature_result), intent(out) :: result
        integer(i8), intent(in) :: min_eval, max_eval
        logical, intent(in) :: randomized
        integer :: dim, nreg, worst, d, newidx, status
        integer(i8) :: neval, used, seq
        real(dp), allocatable :: lo(:, :), hi(:, :), val(:, :), err(:, :), quality(:)
        real(dp), allocatable :: total(:), terr(:), midlo(:), midhi(:)
        type(rng_state) :: rng

        dim = size(lower)
        call init_result(result, max(ncomp, 0))
        if (ncomp <= 0 .or. dim <= 0 .or. size(upper) /= dim .or. any(upper < lower)) then
            result%return_code = CUBATURE_BADARG
            return
        end if
        allocate(lo(dim, max_regions), hi(dim, max_regions), val(ncomp, max_regions), err(ncomp, max_regions))
        allocate(quality(max_regions), total(ncomp), terr(ncomp), midlo(dim), midhi(dim))
        call rng_seed(rng, merge(seed, 104729, seed /= 0))
        seq = 1_i8
        lo(:, 1) = lower; hi(:, 1) = upper
        call sample_region(f, fv, lo(:, 1), hi(:, 1), ncomp, nsample, seq, rng, randomized, &
            val(:, 1), err(:, 1), used, status)
        if (status /= 0) then
            result%return_code = CUBATURE_BADARG
            return
        end if
        seq = seq + used
        neval = used
        quality(1) = maxval(err(:, 1))
        nreg = 1
        do
            total = sum(val(:, 1:nreg), dim=2)
            terr = sqrt(sum(err(:, 1:nreg) ** 2, dim=2))
            if (neval >= min_eval .and. errors_converged(total, terr, abs_tol, rel_tol, ERROR_INDIVIDUAL)) exit
            if (neval + 2_i8 * nsample > max_eval .and. max_eval > 0_i8) exit
            if (nreg >= max_regions) exit
            worst = maxloc(quality(1:nreg), dim=1)
            d = maxloc(hi(:, worst) - lo(:, worst), dim=1)
            midlo = lo(:, worst); midhi = hi(:, worst)
            midlo(d) = 0.5_dp * (lo(d, worst) + hi(d, worst))
            midhi(d) = midlo(d)
            newidx = nreg + 1
            lo(:, newidx) = midlo; hi(:, newidx) = hi(:, worst)
            hi(:, worst) = midhi
            call sample_region(f, fv, lo(:, worst), hi(:, worst), ncomp, nsample, seq, rng, randomized, &
                val(:, worst), err(:, worst), used, status)
            seq = seq + used; neval = neval + used
            call sample_region(f, fv, lo(:, newidx), hi(:, newidx), ncomp, nsample, seq, rng, randomized, &
                val(:, newidx), err(:, newidx), used, status)
            seq = seq + used; neval = neval + used
            quality(worst) = maxval(err(:, worst))
            quality(newidx) = maxval(err(:, newidx))
            nreg = newidx
        end do
        result%integral = sum(val(:, 1:nreg), dim=2)
        result%error = sqrt(sum(err(:, 1:nreg) ** 2, dim=2))
        result%evaluations = neval
        result%nregions = nreg
        if (errors_converged(result%integral, result%error, abs_tol, rel_tol, ERROR_INDIVIDUAL)) then
            result%return_code = CUBATURE_SUCCESS
        else
            result%return_code = CUBATURE_MAXEVAL
        end if
        result%prob = 1.0_dp
    end subroutine adaptive_stratified

    subroutine sample_region(f, fv, lower, upper, ncomp, nsample, seq0, rng, randomized, val, err, neval, status)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp, nsample
        integer(i8), intent(in) :: seq0
        type(rng_state), intent(inout) :: rng
        logical, intent(in) :: randomized
        real(dp), intent(out) :: val(ncomp), err(ncomp)
        integer(i8), intent(out) :: neval
        integer, intent(out) :: status
        integer :: dim, i, r, nrep, nr
        real(dp) :: vol
        real(dp), allocatable :: u(:), x(:), fx(:), shift(:), repmean(:, :)

        dim = size(lower)
        nrep = min(8, max(2, nsample / 16))
        nr = max(4, nsample / nrep)
        allocate(u(dim), x(dim), fx(ncomp), shift(dim), repmean(ncomp, nrep))
        repmean = 0.0_dp
        vol = product(upper - lower)
        do r = 1, nrep
            if (randomized) then
                do i = 1, dim
                    shift(i) = rng_uniform(rng)
                end do
            else
                do i = 1, dim
                    shift(i) = modulo(real(r * (2 * i + 1), dp) * 0.6180339887498948482_dp, 1.0_dp)
                end do
            end if
            do i = 1, nr
                call halton_point(seq0 + int((r - 1) * nr + i, i8), dim, u, shift)
                x = lower + (upper - lower) * u
                call eval_point(f, fv, x, fx)
                repmean(:, r) = repmean(:, r) + fx
            end do
            repmean(:, r) = vol * repmean(:, r) / real(nr, dp)
        end do
        val = sum(repmean, dim=2) / real(nrep, dp)
        if (nrep > 1) then
            do i = 1, ncomp
                err(i) = sqrt(sum((repmean(i, :) - val(i)) ** 2) / real(nrep * (nrep - 1), dp))
            end do
        else
            err = huge(1.0_dp)
        end if
        neval = int(nrep * nr, i8)
        status = 0
    end subroutine sample_region

    subroutine vegas_core(f, fv, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(vegas_options), intent(in), optional :: options
        type(vegas_options) :: opt
        real(dp) :: rt, at
        integer :: dim, iter, nsamp, i, d, b, nb, status
        integer(i8) :: neval
        real(dp), allocatable :: grid(:, :), imp(:, :), sumy(:), sumy2(:), y(:), fx(:), x(:)
        real(dp), allocatable :: iter_val(:), iter_err(:), wsum(:), vsum(:)
        real(dp) :: u, uu, width, jac, var
        type(rng_state) :: rng

        opt = vegas_options()
        if (present(options)) opt = options
        rt = 1.0e-5_dp; at = 1.0e-12_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        dim = size(lower)
        call init_result(result, max(ncomp, 0))
        if (ncomp <= 0 .or. dim <= 0 .or. size(upper) /= dim .or. any(upper < lower)) then
            result%return_code = CUBATURE_BADARG
            return
        end if
        nb = max(4, min(128, opt%nbins))
        allocate(grid(0:nb, dim), imp(nb, dim), sumy(ncomp), sumy2(ncomp), y(ncomp), fx(ncomp), x(dim))
        allocate(iter_val(ncomp), iter_err(ncomp), wsum(ncomp), vsum(ncomp))
        do d = 1, dim
            do b = 0, nb
                grid(b, d) = real(b, dp) / real(nb, dp)
            end do
        end do
        call rng_seed(rng, merge(opt%seed, 12345, opt%seed /= 0))
        neval = 0_i8
        wsum = 0.0_dp; vsum = 0.0_dp
        result%integral = 0.0_dp; result%error = huge(1.0_dp)
        status = CUBATURE_MAXEVAL
        do iter = 1, max(1, opt%max_iter)
            nsamp = max(64, opt%nstart + (iter - 1) * opt%nincrease)
            if (opt%max_eval > 0_i8) nsamp = min(nsamp, int(max(0_i8, opt%max_eval - neval)))
            if (nsamp < 16) exit
            imp = 0.0_dp; sumy = 0.0_dp; sumy2 = 0.0_dp
            do i = 1, nsamp
                jac = 1.0_dp
                do d = 1, dim
                    ! Independent importance-grid sampling.  The grid is adapted
                    ! between iterations, as in the VEGAS family of algorithms.
                    u = rng_uniform(rng)
                    b = min(nb, int(u * real(nb, dp)) + 1)
                    uu = rng_uniform(rng)
                    width = grid(b, d) - grid(b - 1, d)
                    x(d) = lower(d) + (upper(d) - lower(d)) * (grid(b - 1, d) + uu * width)
                    jac = jac * (upper(d) - lower(d)) * real(nb, dp) * width
                end do
                call eval_point(f, fv, x, fx)
                y = jac * fx
                sumy = sumy + y
                sumy2 = sumy2 + y * y
                do d = 1, dim
                    u = (x(d) - lower(d)) / max(upper(d) - lower(d), tiny(1.0_dp))
                    b = locate_bin(grid(:, d), u, nb)
                    imp(b, d) = imp(b, d) + sum(abs(y))
                end do
            end do
            neval = neval + int(nsamp, i8)
            iter_val = sumy / real(nsamp, dp)
            do i = 1, ncomp
                var = max(0.0_dp, sumy2(i) / real(nsamp, dp) - iter_val(i) ** 2)
                iter_err(i) = sqrt(var / real(nsamp, dp))
                if (iter_err(i) <= 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(iter_val(i)))) then
                    iter_err(i) = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(iter_val(i)))
                end if
                wsum(i) = wsum(i) + 1.0_dp / (iter_err(i) * iter_err(i))
                vsum(i) = vsum(i) + iter_val(i) / (iter_err(i) * iter_err(i))
                result%integral(i) = vsum(i) / wsum(i)
                result%error(i) = sqrt(1.0_dp / wsum(i))
            end do
            call adapt_vegas_grid(grid, imp, nb, dim)
            if (neval >= opt%min_eval .and. errors_converged(result%integral, result%error, at, rt, &
                ERROR_INDIVIDUAL)) then
                status = CUBATURE_SUCCESS
                exit
            end if
            if (opt%max_eval > 0_i8 .and. neval >= opt%max_eval) exit
        end do
        result%evaluations = neval
        result%nregions = 1
        result%return_code = status
        result%prob = 1.0_dp
    end subroutine vegas_core

    subroutine vegas(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(vegas_options), intent(in), optional :: options
        call vegas_core(f=f, lower=lower, upper=upper, ncomp=ncomp, result=result, &
            rel_tol=rel_tol, abs_tol=abs_tol, options=options)
    end subroutine vegas

    subroutine vegas_v(f, lower, upper, ncomp, result, rel_tol, abs_tol, options)
        procedure(cubature_integrand_v) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: ncomp
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(vegas_options), intent(in), optional :: options
        call vegas_core(fv=f, lower=lower, upper=upper, ncomp=ncomp, result=result, &
            rel_tol=rel_tol, abs_tol=abs_tol, options=options)
    end subroutine vegas_v

    integer function locate_bin(grid, u, nb) result(b)
        real(dp), intent(in) :: grid(0:), u
        integer, intent(in) :: nb
        integer :: lo, hi, mid
        if (u <= grid(0)) then
            b = 1; return
        end if
        if (u >= grid(nb)) then
            b = nb; return
        end if
        lo = 0; hi = nb
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (u >= grid(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        b = lo + 1
    end function locate_bin

    subroutine adapt_vegas_grid(grid, imp, nb, dim)
        integer, intent(in) :: nb, dim
        real(dp), intent(inout) :: grid(0:nb, dim)
        real(dp), intent(in) :: imp(nb, dim)
        real(dp), allocatable :: smooth(:), cum(:), old(:)
        real(dp) :: total, target, frac
        integer :: d, b, j
        allocate(smooth(nb), cum(0:nb), old(0:nb))
        do d = 1, dim
            if (nb == 1) then
                smooth(1) = sqrt(max(imp(1, d), tiny(1.0_dp)))
            else
                smooth(1) = sqrt(max(imp(1, d) + 0.5_dp * imp(2, d), tiny(1.0_dp)))
                do b = 2, nb - 1
                    smooth(b) = imp(b, d) + 0.5_dp * (imp(b - 1, d) + imp(b + 1, d))
                    smooth(b) = sqrt(max(smooth(b), tiny(1.0_dp)))
                end do
                smooth(nb) = sqrt(max(imp(nb, d) + 0.5_dp * imp(nb - 1, d), tiny(1.0_dp)))
            end if
            total = sum(smooth)
            if (total <= tiny(1.0_dp)) cycle
            old = grid(:, d)
            cum(0) = 0.0_dp
            do b = 1, nb
                cum(b) = cum(b - 1) + smooth(b)
            end do
            grid(0, d) = 0.0_dp
            grid(nb, d) = 1.0_dp
            j = 1
            do b = 1, nb - 1
                target = total * real(b, dp) / real(nb, dp)
                do while (j < nb .and. cum(j) < target)
                    j = j + 1
                end do
                frac = (target - cum(j - 1)) / max(cum(j) - cum(j - 1), tiny(1.0_dp))
                grid(b, d) = old(j - 1) + frac * (old(j) - old(j - 1))
            end do
        end do
    end subroutine adapt_vegas_grid

    subroutine eval_point(f, fv, x, value)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value(:)
        real(dp) :: xx(size(x), 1), vv(size(value), 1)
        if (present(f)) then
            call f(x, value)
        else if (present(fv)) then
            xx(:, 1) = x
            call fv(xx, vv)
            value = vv(:, 1)
        else
            value = 0.0_dp
        end if
    end subroutine eval_point

    subroutine init_result(result, fdim)
        type(cubature_result), intent(out) :: result
        integer, intent(in) :: fdim
        allocate(result%integral(fdim), result%error(fdim), result%prob(fdim))
        result%integral = 0.0_dp
        result%error = huge(1.0_dp)
        result%prob = 1.0_dp
        result%evaluations = 0_i8
        result%return_code = CUBATURE_SUCCESS
        result%nregions = 0
    end subroutine init_result

end module cuba_mod
