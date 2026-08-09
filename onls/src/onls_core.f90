! SPDX-License-Identifier: GPL-2.0-or-later
module onls_core
    use onls_kinds, only : dp
    use onls_lm, only : lm_control, lm_result, lm_solve
    use onls_minimize, only : brent_minimize
    use onls_linalg, only : invert_spd
    implicit none
    private
    public :: dp, model_fn, onls_control, onls_result, fit_onls
    public :: orthogonal_residuals, vertical_loglik, orthogonal_loglik
    public :: check_orthogonality

    abstract interface
        subroutine model_fn(x, par, y, ierr)
            import dp
            real(dp), intent(in) :: x(:), par(:)
            real(dp), intent(out) :: y(:)
            integer, intent(out) :: ierr
        end subroutine model_fn
    end interface

    type :: onls_control
        type(lm_control) :: lm = lm_control()
        integer :: window = 12
        real(dp) :: extend(2) = [0.2_dp, 0.2_dp]
        real(dp) :: projection_tol = sqrt(epsilon(1.0_dp))
        logical :: mimic_r_unsorted_weights = .true.
    end type onls_control

    type :: onls_result
        real(dp), allocatable :: par_nls(:), par_onls(:)
        real(dp), allocatable :: x(:), y(:), weights(:)
        real(dp), allocatable :: x0(:), y0(:)
        real(dp), allocatable :: fitted_nls(:), fitted_onls(:)
        real(dp), allocatable :: resid_nls(:), resid_onls(:), resid_o(:), distance_o(:)
        real(dp), allocatable :: gradient(:,:), covariance(:,:), stderr(:)
        logical, allocatable :: ortho(:)
        real(dp), allocatable :: ortho_angle(:)
        real(dp) :: rss_nls = huge(1.0_dp)
        real(dp) :: rss_vertical = huge(1.0_dp)
        real(dp) :: rss_orthogonal = huge(1.0_dp)
        real(dp) :: sigma_vertical = huge(1.0_dp)
        real(dp) :: sigma_orthogonal = huge(1.0_dp)
        integer :: niter_nls = 0, niter_onls = 0
        integer :: info_nls = 0, info_onls = 0
        logical :: converged = .false.
        character(len=128) :: message = ''
    end type onls_result

    type :: projection_context
        procedure(model_fn), pointer, nopass :: model => null()
        real(dp), allocatable :: par(:)
        real(dp) :: xobs = 0.0_dp
        real(dp) :: yobs = 0.0_dp
        integer :: ierr = 0
    end type projection_context

    type :: model_holder
        procedure(model_fn), pointer, nopass :: model => null()
        real(dp), allocatable :: x(:), y(:), wroot(:)
        real(dp), allocatable :: fixed_start(:)
        logical, allocatable :: fixed(:)
        integer :: window = 12
        real(dp) :: lower_x = 0.0_dp, upper_x = 0.0_dp
        real(dp) :: projection_tol = sqrt(epsilon(1.0_dp))
        real(dp), allocatable :: x0(:)
    end type model_holder
contains
    subroutine fit_onls(model, x, y, start, result, control, lower, upper, weights, fixed)
        procedure(model_fn) :: model
        real(dp), intent(in) :: x(:), y(:), start(:)
        type(onls_result), intent(out) :: result
        type(onls_control), intent(in), optional :: control
        real(dp), intent(in), optional :: lower(:), upper(:), weights(:)
        logical, intent(in), optional :: fixed(:)
        type(onls_control) :: ctl
        type(model_holder) :: ctx
        type(lm_result) :: nlsfit, ofit
        real(dp), allocatable :: xs(:), ys(:), ws(:), wroot(:), p0(:), lo(:), hi(:)
        integer, allocatable :: ord(:)
        logical, allocatable :: fmask(:)
        integer :: n, p, ierr, i, rdf
        real(dp) :: xrange

        ctl = onls_control()
        if (present(control)) ctl = control
        n = size(x)
        p = size(start)
        if (n < 2 .or. size(y) /= n .or. p < 1) then
            result%message = 'Invalid data or parameter dimensions'
            result%info_onls = -1
            return
        end if
        allocate(xs(n), ys(n), ws(n), wroot(n), ord(n), p0(p), lo(p), hi(p), fmask(p))
        call order_real(x, ord)
        xs = x(ord)
        ys = y(ord)
        ws = 1.0_dp
        if (present(weights)) then
            if (size(weights) /= n .or. any(weights < 0.0_dp)) then
                result%message = 'Weights must be nonnegative and match data length'
                result%info_onls = -2
                return
            end if
            if (ctl%mimic_r_unsorted_weights) then
                ws = weights
            else
                ws = weights(ord)
            end if
        end if
        wroot = sqrt(ws)
        lo = -huge(1.0_dp)
        hi = huge(1.0_dp)
        if (present(lower)) then
            if (size(lower) == 1) then
                lo = lower(1)
            else if (size(lower) == p) then
                lo = lower
            else
                result%message = 'lower must be scalar or match parameter count'
                result%info_onls = -3
                return
            end if
        end if
        if (present(upper)) then
            if (size(upper) == 1) then
                hi = upper(1)
            else if (size(upper) == p) then
                hi = upper
            else
                result%message = 'upper must be scalar or match parameter count'
                result%info_onls = -4
                return
            end if
        end if
        if (any(lo > hi)) then
            result%message = 'lower exceeds upper'
            result%info_onls = -5
            return
        end if
        fmask = .false.
        if (present(fixed)) then
            if (size(fixed) /= p) then
                result%message = 'fixed must match parameter count'
                result%info_onls = -6
                return
            end if
            fmask = fixed
        end if
        ctx%model => model
        ctx%x = xs
        ctx%y = ys
        ctx%wroot = wroot
        ctx%fixed_start = start
        ctx%fixed = fmask
        ctx%window = max(1, ctl%window)
        xrange = maxval(xs)-minval(xs)
        ctx%lower_x = minval(xs) - ctl%extend(1)*xrange
        ctx%upper_x = maxval(xs) + ctl%extend(2)*xrange
        ctx%projection_tol = ctl%projection_tol
        allocate(ctx%x0(n))
        ctx%x0 = xs

        call lm_solve(vertical_residual_dispatch, ctx, start, n, nlsfit, ctl%lm, lo, hi)
        p0 = nlsfit%par
        do i = 1, p
            if (fmask(i)) p0(i) = start(i)
        end do
        call lm_solve(orthogonal_residual_dispatch, ctx, p0, n, ofit, ctl%lm, lo, hi)

        result%par_nls = nlsfit%par
        result%par_onls = ofit%par
        do i = 1, p
            if (fmask(i)) result%par_onls(i) = start(i)
        end do
        result%x = xs
        result%y = ys
        result%weights = ws
        result%x0 = ctx%x0
        allocate(result%y0(n), result%fitted_nls(n), result%fitted_onls(n))
        call evaluate_model(model, xs, result%par_nls, result%fitted_nls, ierr)
        if (ierr /= 0) then
            result%message = 'Model evaluation failed at ordinary NLS solution'
            return
        end if
        ! Upstream stores MINPACK's weighted fvec and subtracts it from RESP.
        result%resid_nls = nlsfit%fvec
        result%fitted_nls = ys - result%resid_nls
        call evaluate_model(model, xs, result%par_onls, result%fitted_onls, ierr)
        if (ierr /= 0) then
            result%message = 'Model evaluation failed at ONLS solution'
            return
        end if
        call evaluate_model(model, result%x0, result%par_onls, result%y0, ierr)
        if (ierr /= 0) then
            result%message = 'Model evaluation failed at orthogonal foot points'
            return
        end if
        result%resid_onls = ys-result%fitted_onls
        result%distance_o = sqrt((xs-result%x0)**2 + (ys-result%y0)**2)
        result%resid_o = wroot * result%distance_o
        result%rss_nls = sum(result%resid_nls**2)
        result%rss_vertical = sum(result%resid_onls**2)
        result%rss_orthogonal = sum(result%resid_o**2)
        result%niter_nls = nlsfit%niter
        result%niter_onls = ofit%niter
        result%info_nls = nlsfit%info
        result%info_onls = ofit%info
        result%converged = ofit%info >= 1 .and. ofit%info <= 4
        result%message = ofit%message
        rdf = count(ws > 0.0_dp) - p
        if (rdf > 0) then
            result%sigma_vertical = sqrt(result%rss_vertical/real(rdf,dp))
            result%sigma_orthogonal = sqrt(result%rss_orthogonal/real(rdf,dp))
        end if
        call parameter_statistics(model, xs, result%par_onls, wroot, result%rss_vertical, rdf, &
                                  result%gradient, result%covariance, result%stderr)
        call check_orthogonality(model, xs, ys, result%x0, result%y0, result%par_onls, &
                                 result%ortho_angle, result%ortho)
    end subroutine fit_onls

    subroutine vertical_residual_dispatch(par, r, rawctx, ierr)
        real(dp), intent(in) :: par(:)
        real(dp), intent(out) :: r(:)
        class(*), intent(inout) :: rawctx
        integer, intent(out) :: ierr
        real(dp), allocatable :: pred(:)
        select type (ctx => rawctx)
        type is (model_holder)
            allocate(pred(size(ctx%x)))
            call evaluate_model(ctx%model, ctx%x, par, pred, ierr)
            if (ierr == 0) r = ctx%wroot * (ctx%y-pred)
        class default
            ierr = 99
        end select
    end subroutine vertical_residual_dispatch

    subroutine orthogonal_residual_dispatch(par_in, r, rawctx, ierr)
        real(dp), intent(in) :: par_in(:)
        real(dp), intent(out) :: r(:)
        class(*), intent(inout) :: rawctx
        integer, intent(out) :: ierr
        real(dp), allocatable :: par(:)
        integer :: i, n, ilo, ihi
        real(dp) :: a, b, xmin, fmin
        select type (ctx => rawctx)
        type is (model_holder)
            n = size(ctx%x)
            allocate(par(size(par_in)))
            par = par_in
            where (ctx%fixed) par = ctx%fixed_start
            ierr = 0
            do i = 1, n
                if (n <= 25) then
                    a = ctx%lower_x
                    b = ctx%upper_x
                else
                    ilo = max(1, i-(ctx%window-1))
                    ihi = min(n, i+(ctx%window-1))
                    if (i <= ctx%window) then
                        a = ctx%lower_x
                        b = ctx%x(ihi)
                    else if (i >= n-(ctx%window-1)) then
                        a = ctx%x(ilo)
                        b = ctx%upper_x
                    else
                        a = ctx%x(ilo)
                        b = ctx%x(ihi)
                    end if
                end if
                call project_one(ctx%model, par, ctx%x(i), ctx%y(i), a, b, &
                                 ctx%projection_tol, xmin, fmin, ierr)
                if (ierr /= 0) return
                ctx%x0(i) = xmin
                r(i) = ctx%wroot(i) * fmin
            end do
        class default
            ierr = 99
        end select
    end subroutine orthogonal_residual_dispatch

    subroutine project_one(model, par, xobs, yobs, a, b, tol, xmin, fmin, ierr)
        procedure(model_fn) :: model
        real(dp), intent(in) :: par(:), xobs, yobs, a, b, tol
        real(dp), intent(out) :: xmin, fmin
        integer, intent(out) :: ierr
        type(projection_context) :: pctx
        pctx%model => model
        pctx%par = par
        pctx%xobs = xobs
        pctx%yobs = yobs
        pctx%ierr = 0
        call brent_minimize(distance_dispatch, pctx, a, b, xmin, fmin, tol)
        ierr = pctx%ierr
        if (.not.(fmin <= huge(1.0_dp))) ierr = 1
    end subroutine project_one

    function distance_dispatch(xx, rawctx) result(dist)
        real(dp), intent(in) :: xx
        class(*), intent(inout) :: rawctx
        real(dp) :: dist
        real(dp) :: xv(1), yv(1)
        integer :: ie
        select type (ctx => rawctx)
        type is (projection_context)
            xv(1) = xx
            call ctx%model(xv, ctx%par, yv, ie)
            if (ie /= 0 .or. .not.(abs(yv(1)) <= huge(1.0_dp))) then
                dist = huge(1.0_dp)/100.0_dp
                ctx%ierr = 1
            else
                dist = sqrt((xx-ctx%xobs)**2 + (yv(1)-ctx%yobs)**2)
            end if
        class default
            dist = huge(1.0_dp)/100.0_dp
        end select
    end function distance_dispatch

    subroutine evaluate_model(model, x, par, y, ierr)
        procedure(model_fn) :: model
        real(dp), intent(in) :: x(:), par(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: ierr
        call model(x, par, y, ierr)
    end subroutine evaluate_model

    subroutine parameter_statistics(model, x, par, wroot, rss, rdf, grad, covar, stderr)
        procedure(model_fn) :: model
        real(dp), intent(in) :: x(:), par(:), wroot(:), rss
        integer, intent(in) :: rdf
        real(dp), allocatable, intent(out) :: grad(:,:), covar(:,:), stderr(:)
        real(dp), allocatable :: pp(:), pm(:), yp(:), ym(:), jtwj(:,:), inv(:,:)
        real(dp) :: h, s2
        integer :: n, p, j, ierr
        logical :: ok

        n = size(x)
        p = size(par)
        allocate(grad(n,p), pp(p), pm(p), yp(n), ym(n), jtwj(p,p), inv(p,p), stderr(p))
        do j = 1, p
            h = sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(par(j)))
            pp = par
            pm = par
            pp(j) = pp(j)+h
            pm(j) = pm(j)-h
            call evaluate_model(model,x,pp,yp,ierr)
            if (ierr /= 0) then
                grad(:,j) = 0.0_dp
                cycle
            end if
            call evaluate_model(model,x,pm,ym,ierr)
            if (ierr /= 0) then
                grad(:,j) = 0.0_dp
                cycle
            end if
            grad(:,j) = (yp-ym)/(2.0_dp*h)
        end do
        grad = spread(wroot,2,p)*grad
        jtwj = matmul(transpose(grad),grad)
        call invert_spd(jtwj,inv,ok)
        if (.not. ok .or. rdf <= 0) then
            allocate(covar(p,p))
            covar = huge(1.0_dp)
            stderr = huge(1.0_dp)
            return
        end if
        s2 = rss/real(rdf,dp)
        allocate(covar(p,p))
        covar = s2*inv
        do j=1,p
            stderr(j)=sqrt(max(0.0_dp,covar(j,j)))
        end do
    end subroutine parameter_statistics

    subroutine check_orthogonality(model, x, y, x0, y0, par, angle, ok)
        procedure(model_fn) :: model
        real(dp), intent(in) :: x(:), y(:), x0(:), y0(:), par(:)
        real(dp), allocatable, intent(out) :: angle(:)
        logical, allocatable, intent(out) :: ok(:)
        real(dp), allocatable :: up(:), down(:), yup(:), ydown(:), ms(:), ds(:), ta(:)
        real(dp) :: epsx
        integer :: n, ierr, i
        n=size(x)
        allocate(angle(n),ok(n),up(n),down(n),yup(n),ydown(n),ms(n),ds(n),ta(n))
        epsx=sqrt(epsilon(1.0_dp))
        up=x0+epsx
        down=x0-epsx
        call evaluate_model(model,up,par,yup,ierr)
        if (ierr /= 0) then
            angle=0.0_dp; ok=.false.; return
        end if
        call evaluate_model(model,down,par,ydown,ierr)
        if (ierr /= 0) then
            angle=0.0_dp; ok=.false.; return
        end if
        ms=(yup-ydown)/(up-down)
        do i=1,n
            if (abs(x(i)-x0(i)) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i)))) then
                if (abs(y(i)-y0(i)) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(y(i)))) then
                    ds(i)=huge(1.0_dp)
                else
                    ds(i)=sign(huge(1.0_dp),y(i)-y0(i))
                end if
            else
                ds(i)=(y(i)-y0(i))/(x(i)-x0(i))
            end if
            if (abs(ds(i)) > 0.1_dp*huge(1.0_dp)) then
                ta(i)=abs(1.0_dp/ms(i))
            else
                ta(i)=abs((ms(i)-ds(i))/(1.0_dp+ms(i)*ds(i)))
            end if
        end do
        angle=atan(ta)*180.0_dp/acos(-1.0_dp)
        ok=angle>89.95_dp .and. angle<90.05_dp
    end subroutine check_orthogonality

    function vertical_loglik(result) result(val)
        type(onls_result), intent(in) :: result
        real(dp) :: val
        val=gaussian_loglik(result%resid_onls,result%weights)
    end function vertical_loglik

    function orthogonal_loglik(result) result(val)
        type(onls_result), intent(in) :: result
        real(dp) :: val
        val=gaussian_loglik(result%resid_o,result%weights)
    end function orthogonal_loglik

    function gaussian_loglik(res,w) result(val)
        real(dp), intent(in) :: res(:),w(:)
        real(dp) :: val, swlog, rss
        integer :: n,i
        n=size(res)
        swlog=0.0_dp
        rss=0.0_dp
        do i=1,n
            if (w(i)>0.0_dp) swlog=swlog+log(w(i))
            rss=rss+w(i)*res(i)*res(i)
        end do
        if (rss<=0.0_dp) then
            val=huge(1.0_dp)
        else
            val=-0.5_dp*real(n,dp)*(log(2.0_dp*acos(-1.0_dp))+1.0_dp-log(real(n,dp)) &
                -swlog+log(rss))
        end if
    end function gaussian_loglik

    function orthogonal_residuals(result) result(r)
        type(onls_result), intent(in) :: result
        real(dp), allocatable :: r(:)
        r=result%resid_o
    end function orthogonal_residuals

    subroutine order_real(x,ord)
        real(dp),intent(in)::x(:)
        integer,intent(out)::ord(:)
        integer::i,j,key
        do i=1,size(x); ord(i)=i; end do
        do i=2,size(x)
            key=ord(i); j=i-1
            do while(j>=1)
                if(x(ord(j))<=x(key)) exit
                ord(j+1)=ord(j); j=j-1
            end do
            ord(j+1)=key
        end do
    end subroutine order_real
end module onls_core
