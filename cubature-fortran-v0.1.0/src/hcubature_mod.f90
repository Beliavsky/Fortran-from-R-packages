! SPDX-License-Identifier: GPL-3.0-or-later
!
! Free-format Fortran translation of the numerical method in
! Steven G. Johnson's hcubature.c (cubature 1.0.4): globally adaptive
! h-refinement, 15-point Gauss-Kronrod in one dimension, and the
! embedded Genz-Malik degree-7/degree-5 rule in multiple dimensions.
module hcubature_mod
    use cubature_kinds, only : dp, i8
    use cubature_types, only : cubature_result, cubature_integrand, cubature_integrand_v, &
        CUBATURE_SUCCESS, CUBATURE_MAXEVAL, CUBATURE_BADARG, CUBATURE_FAILURE, &
        ERROR_INDIVIDUAL, errors_converged
    implicit none
    private
    public :: hcubature, hcubature_v, adapt_integrate, adaptIntegrate

    interface adaptIntegrate
        module procedure adapt_integrate
    end interface adaptIntegrate

contains

    subroutine hcubature_core(f, fv, lower, upper, fdim, result, rel_tol, abs_tol, max_eval, norm)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        integer, intent(in), optional :: norm

        real(dp) :: rtol, atol
        integer(i8) :: maxev, neval, npoints, max_regions_i8
        integer :: nrm, dim, nreg, max_regions, worst, d, newidx
        real(dp), allocatable :: lo(:, :), hi(:, :), vals(:, :), errs(:, :), errmax(:)
        integer, allocatable :: split_dim(:)
        real(dp), allocatable :: total_val(:), total_err(:), midlo(:), midhi(:)
        integer(i8) :: used
        integer :: status
        logical :: ok

        dim = size(lower)
        rtol = 1.0e-5_dp
        atol = 100.0_dp * epsilon(1.0_dp)
        maxev = 0_i8
        nrm = ERROR_INDIVIDUAL
        if (present(rel_tol)) rtol = rel_tol
        if (present(abs_tol)) atol = abs_tol
        if (present(max_eval)) maxev = max_eval
        if (present(norm)) nrm = norm

        call init_result(result, max(fdim, 0))
        if (fdim <= 0 .or. dim <= 0 .or. size(upper) /= dim .or. any(upper < lower) .or. &
            rtol < 0.0_dp .or. atol < 0.0_dp) then
            result%return_code = CUBATURE_BADARG
            return
        end if
        if (any(.not. is_finite(lower)) .or. any(.not. is_finite(upper))) then
            result%return_code = CUBATURE_BADARG
            return
        end if

        if (dim == 1) then
            npoints = 15_i8
        else
            if (dim >= bit_size(1) - 1) then
                result%return_code = CUBATURE_BADARG
                return
            end if
            npoints = 1_i8 + 4_i8 * dim + 2_i8 * dim * (dim - 1) + shiftl(1_i8, dim)
        end if
        if (maxev > 0_i8) then
            max_regions_i8 = max(4_i8, maxev / max(1_i8, npoints) + 4_i8)
            max_regions_i8 = min(max_regions_i8, 200000_i8)
        else
            max_regions_i8 = 200000_i8
        end if
        max_regions = int(max_regions_i8)

        allocate(lo(dim, max_regions), hi(dim, max_regions))
        allocate(vals(fdim, max_regions), errs(fdim, max_regions), errmax(max_regions))
        allocate(split_dim(max_regions), total_val(fdim), total_err(fdim))
        allocate(midlo(dim), midhi(dim))

        lo(:, 1) = lower
        hi(:, 1) = upper
        call eval_region(f, fv, lo(:, 1), hi(:, 1), fdim, vals(:, 1), errs(:, 1), split_dim(1), used, status)
        if (status /= 0) then
            result%return_code = CUBATURE_FAILURE
            return
        end if
        neval = used
        errmax(1) = maxval(errs(:, 1))
        nreg = 1

        do
            total_val = sum(vals(:, 1:nreg), dim=2)
            total_err = sum(errs(:, 1:nreg), dim=2)
            ok = errors_converged(total_val, total_err, atol, rtol, nrm)
            if (ok .and. neval >= 1_i8) exit
            if (maxev > 0_i8 .and. neval + 2_i8 * npoints > maxev) exit
            if (nreg >= max_regions) exit

            worst = maxloc(errmax(1:nreg), dim=1)
            d = split_dim(worst)
            if (d < 1 .or. d > dim) d = maxloc(hi(:, worst) - lo(:, worst), dim=1)
            midlo = lo(:, worst)
            midhi = hi(:, worst)
            midlo(d) = 0.5_dp * (lo(d, worst) + hi(d, worst))
            midhi(d) = midlo(d)

            newidx = nreg + 1
            lo(:, newidx) = midlo
            hi(:, newidx) = hi(:, worst)
            lo(:, worst) = lo(:, worst)
            hi(:, worst) = midhi

            call eval_region(f, fv, lo(:, worst), hi(:, worst), fdim, vals(:, worst), errs(:, worst), &
                split_dim(worst), used, status)
            if (status /= 0) then
                result%return_code = CUBATURE_FAILURE
                return
            end if
            neval = neval + used
            call eval_region(f, fv, lo(:, newidx), hi(:, newidx), fdim, vals(:, newidx), errs(:, newidx), &
                split_dim(newidx), used, status)
            if (status /= 0) then
                result%return_code = CUBATURE_FAILURE
                return
            end if
            neval = neval + used
            errmax(worst) = maxval(errs(:, worst))
            errmax(newidx) = maxval(errs(:, newidx))
            nreg = newidx
        end do

        result%integral = sum(vals(:, 1:nreg), dim=2)
        result%error = sum(errs(:, 1:nreg), dim=2)
        result%evaluations = neval
        result%nregions = nreg
        if (errors_converged(result%integral, result%error, atol, rtol, nrm)) then
            result%return_code = CUBATURE_SUCCESS
        else
            result%return_code = CUBATURE_MAXEVAL
        end if
        if (allocated(result%prob)) result%prob = 1.0_dp
    end subroutine hcubature_core

    subroutine hcubature(f, lower, upper, fdim, result, rel_tol, abs_tol, max_eval, norm)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        integer, intent(in), optional :: norm
        call hcubature_core(f=f, lower=lower, upper=upper, fdim=fdim, result=result, &
            rel_tol=rel_tol, abs_tol=abs_tol, max_eval=max_eval, norm=norm)
    end subroutine hcubature

    subroutine hcubature_v(f, lower, upper, fdim, result, rel_tol, abs_tol, max_eval, norm)
        procedure(cubature_integrand_v) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: rel_tol, abs_tol
        integer(i8), intent(in), optional :: max_eval
        integer, intent(in), optional :: norm
        call hcubature_core(fv=f, lower=lower, upper=upper, fdim=fdim, result=result, &
            rel_tol=rel_tol, abs_tol=abs_tol, max_eval=max_eval, norm=norm)
    end subroutine hcubature_v

    subroutine adapt_integrate(f, lower, upper, fdim, result, tol, abs_error, max_eval, norm)
        procedure(cubature_integrand) :: f
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        type(cubature_result), intent(out) :: result
        real(dp), intent(in), optional :: tol, abs_error
        integer(i8), intent(in), optional :: max_eval
        integer, intent(in), optional :: norm
        call hcubature(f, lower, upper, fdim, result, tol, abs_error, max_eval, norm)
    end subroutine adapt_integrate

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

    pure elemental logical function is_finite(x)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
        real(dp), intent(in) :: x
        is_finite = ieee_is_finite(x)
    end function is_finite

    subroutine eval_region(f, fv, lower, upper, fdim, val, err, split_dim, neval, status)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        real(dp), intent(out) :: val(fdim), err(fdim)
        integer, intent(out) :: split_dim
        integer(i8), intent(out) :: neval
        integer, intent(out) :: status
        if (size(lower) == 1) then
            call eval_gk15(f, fv, lower(1), upper(1), fdim, val, err, neval)
            split_dim = 1
            status = 0
        else
            call eval_genz_malik(f, fv, lower, upper, fdim, val, err, split_dim, neval, status)
        end if
    end subroutine eval_region

    subroutine eval_gk15(f, fv, a, b, fdim, val, err, neval)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: a, b
        integer, intent(in) :: fdim
        real(dp), intent(out) :: val(fdim), err(fdim)
        integer(i8), intent(out) :: neval
        real(dp), parameter :: xgk(8) = [ &
            0.991455371120812639206854697526329_dp, &
            0.949107912342758524526189684047851_dp, &
            0.864864423359769072789712788640926_dp, &
            0.741531185599394439863864773280788_dp, &
            0.586087235467691130294144838258730_dp, &
            0.405845151377397166906606412076961_dp, &
            0.207784955007898467600689403773245_dp, 0.0_dp ]
        real(dp), parameter :: wg(4) = [ &
            0.129484966168869693270611432679082_dp, &
            0.279705391489276667901467771423780_dp, &
            0.381830050505118944950369775488975_dp, &
            0.417959183673469387755102040816327_dp ]
        real(dp), parameter :: wgk(8) = [ &
            0.022935322010529224963732008058970_dp, &
            0.063092092629978553290700663189204_dp, &
            0.104790010322250183839876322541518_dp, &
            0.140653259715525918745189590510238_dp, &
            0.169004726639267902826583426598550_dp, &
            0.190350578064785409913256402421014_dp, &
            0.204432940075298892414161999234649_dp, &
            0.209482141084727828012999174891714_dp ]
        real(dp) :: center, halfwidth, x(1)
        real(dp), allocatable :: fc(:), f1(:), f2(:), result_k(:), result_g(:), result_abs(:), result_asc(:)
        real(dp), allocatable :: vm(:,:), vp(:,:)
        real(dp) :: absc, mean, scale, min_err
        integer :: j, k

        allocate(fc(fdim), f1(fdim), f2(fdim), result_k(fdim), result_g(fdim), result_abs(fdim), result_asc(fdim))
        allocate(vm(fdim, 7), vp(fdim, 7))
        center = 0.5_dp * (a + b)
        halfwidth = 0.5_dp * (b - a)
        x(1) = center
        call eval_point(f, fv, x, fc)
        result_g = wg(4) * fc
        result_k = wgk(8) * fc
        result_abs = wgk(8) * abs(fc)
        do j = 1, 7
            absc = halfwidth * xgk(j)
            x(1) = center - absc; call eval_point(f, fv, x, vm(:, j))
            x(1) = center + absc; call eval_point(f, fv, x, vp(:, j))
            result_k = result_k + wgk(j) * (vm(:, j) + vp(:, j))
            result_abs = result_abs + wgk(j) * (abs(vm(:, j)) + abs(vp(:, j)))
            select case (j)
            case (2)
                result_g = result_g + wg(1) * (vm(:, j) + vp(:, j))
            case (4)
                result_g = result_g + wg(2) * (vm(:, j) + vp(:, j))
            case (6)
                result_g = result_g + wg(3) * (vm(:, j) + vp(:, j))
            end select
        end do
        val = result_k * halfwidth
        result_abs = result_abs * halfwidth
        do k = 1, fdim
            mean = 0.5_dp * result_k(k)
            result_asc(k) = wgk(8) * abs(fc(k) - mean)
            do j = 1, 7
                result_asc(k) = result_asc(k) + wgk(j) * (abs(vm(k, j) - mean) + abs(vp(k, j) - mean))
            end do
            result_asc(k) = result_asc(k) * halfwidth
            err(k) = abs(result_k(k) - result_g(k)) * halfwidth
            if (result_asc(k) > 0.0_dp .and. err(k) > 0.0_dp) then
                scale = (200.0_dp * err(k) / result_asc(k)) ** 1.5_dp
                err(k) = result_asc(k) * min(1.0_dp, scale)
            end if
            if (result_abs(k) > tiny(1.0_dp) / (50.0_dp * epsilon(1.0_dp))) then
                min_err = 50.0_dp * epsilon(1.0_dp) * result_abs(k)
                err(k) = max(err(k), min_err)
            end if
        end do
        neval = 15_i8
    end subroutine eval_gk15

    subroutine eval_genz_malik(f, fv, lower, upper, fdim, val, err, split_dim, neval, status)
        procedure(cubature_integrand), optional :: f
        procedure(cubature_integrand_v), optional :: fv
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: fdim
        real(dp), intent(out) :: val(fdim), err(fdim)
        integer, intent(out) :: split_dim
        integer(i8), intent(out) :: neval
        integer, intent(out) :: status
        real(dp), parameter :: lambda2 = 0.3585685828003180919906451539079375_dp
        real(dp), parameter :: lambda4 = 0.9486832980505137995996680633298156_dp
        real(dp), parameter :: lambda5 = 0.6882472016116852977216287342936235_dp
        real(dp), parameter :: weight2 = 980.0_dp / 6561.0_dp
        real(dp), parameter :: weight4 = 200.0_dp / 19683.0_dp
        real(dp), parameter :: weighte2 = 245.0_dp / 486.0_dp
        real(dp), parameter :: weighte4 = 25.0_dp / 729.0_dp
        real(dp), parameter :: ratio = (lambda2 * lambda2) / (lambda4 * lambda4)
        integer :: dim, i, j, k, mask, ncorners
        real(dp) :: volume, weight1, weight3, weight5, weighte1, weighte3
        real(dp), allocatable :: center(:), half(:), x(:), f0(:), fm(:), fp(:), g1(:), g2(:)
        real(dp), allocatable :: sum2(:), sum3(:), sum4(:), sum5(:), res5(:), diff(:)
        real(dp) :: maxdiff, tol_diff

        dim = size(lower)
        status = 0
        if (dim < 2 .or. dim >= bit_size(1) - 1) then
            status = 1
            val = 0.0_dp; err = huge(1.0_dp); split_dim = 1; neval = 0_i8
            return
        end if
        allocate(center(dim), half(dim), x(dim), f0(fdim), fm(fdim), fp(fdim), g1(fdim), g2(fdim))
        allocate(sum2(fdim), sum3(fdim), sum4(fdim), sum5(fdim), res5(fdim), diff(dim))
        center = 0.5_dp * (lower + upper)
        half = 0.5_dp * (upper - lower)
        volume = product(upper - lower)
        weight1 = real(12824 - 9120 * dim + 400 * dim * dim, dp) / 19683.0_dp
        weight3 = real(1820 - 400 * dim, dp) / 19683.0_dp
        weight5 = 6859.0_dp / 19683.0_dp / real(shiftl(1, dim), dp)
        weighte1 = real(729 - 950 * dim + 50 * dim * dim, dp) / 729.0_dp
        weighte3 = real(265 - 100 * dim, dp) / 1458.0_dp
        sum2 = 0.0_dp; sum3 = 0.0_dp; sum4 = 0.0_dp; sum5 = 0.0_dp; diff = 0.0_dp
        x = center
        call eval_point(f, fv, x, f0)
        neval = 1_i8
        do i = 1, dim
            x = center; x(i) = center(i) - lambda2 * half(i); call eval_point(f, fv, x, fm)
            x(i) = center(i) + lambda2 * half(i); call eval_point(f, fv, x, fp)
            g1 = fm + fp
            sum2 = sum2 + g1
            x = center; x(i) = center(i) - lambda4 * half(i); call eval_point(f, fv, x, fm)
            x(i) = center(i) + lambda4 * half(i); call eval_point(f, fv, x, fp)
            g2 = fm + fp
            sum3 = sum3 + g2
            diff(i) = sum(abs(g1 - 2.0_dp * f0 - ratio * (g2 - 2.0_dp * f0)))
            neval = neval + 4_i8
        end do
        do i = 1, dim - 1
            do j = i + 1, dim
                do k = 0, 3
                    x = center
                    if (btest(k, 0)) then
                        x(i) = center(i) + lambda4 * half(i)
                    else
                        x(i) = center(i) - lambda4 * half(i)
                    end if
                    if (btest(k, 1)) then
                        x(j) = center(j) + lambda4 * half(j)
                    else
                        x(j) = center(j) - lambda4 * half(j)
                    end if
                    call eval_point(f, fv, x, g1)
                    sum4 = sum4 + g1
                    neval = neval + 1_i8
                end do
            end do
        end do
        ncorners = shiftl(1, dim)
        do mask = 0, ncorners - 1
            do i = 1, dim
                if (btest(mask, i - 1)) then
                    x(i) = center(i) + lambda5 * half(i)
                else
                    x(i) = center(i) - lambda5 * half(i)
                end if
            end do
            call eval_point(f, fv, x, g1)
            sum5 = sum5 + g1
            neval = neval + 1_i8
        end do
        val = volume * (weight1 * f0 + weight2 * sum2 + weight3 * sum3 + weight4 * sum4 + weight5 * sum5)
        res5 = volume * (weighte1 * f0 + weighte2 * sum2 + weighte3 * sum3 + weighte4 * sum4)
        err = abs(res5 - val)

        split_dim = 1
        maxdiff = diff(1)
        tol_diff = sum(err) / max(volume * 10.0_dp ** min(dim, 30), tiny(1.0_dp))
        do i = 2, dim
            if (diff(i) - maxdiff > tol_diff) then
                maxdiff = diff(i)
                split_dim = i
            else if (abs(diff(i) - maxdiff) <= tol_diff .and. half(i) > half(split_dim)) then
                split_dim = i
            end if
        end do
    end subroutine eval_genz_malik

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

end module hcubature_mod
