! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_arimax
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
  use tsa_kinds, only : dp
  use tsa_types, only : arimax_result, bootstrap_result, transfer_spec
  use tsa_utils, only : difference_series, seasonal_difference, &
    polynomial_convolution, convolution_filter, recursive_filter, &
    sample_with_replacement
  use tseries_linalg, only : least_squares, invert_matrix_lu, solve_linear, right_singular_vectors
  use tseries_optimize, only : bfgs, optim_hessian
  use tseries_random, only : random_normal
  implicit none
  private

  public :: arima_fit, arimax_fit, arima_sim, arima_bootstrap, arima_bootstrap_sample
  public :: transfer_filter, io_regressor

contains

  function arima_fit(x, p, d, q, include_mean, seasonal_p, seasonal_d, &
      seasonal_q, period, method, fixed, init, transform_pars, n_cond, kappa) result(fit)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: p, d, q
    logical, intent(in), optional :: include_mean, transform_pars
    integer, intent(in), optional :: seasonal_p, seasonal_d, seasonal_q, period
    integer, intent(in), optional :: n_cond
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: fixed(:), init(:), kappa
    type(arimax_result) :: fit
    real(dp), allocatable :: empty(:,:)

    allocate(empty(size(x), 0))
    fit = arimax_fit(x, p, d, q, empty, include_mean, seasonal_p, seasonal_d, &
      seasonal_q, period, method=method, fixed=fixed, init=init, &
      transform_pars=transform_pars, n_cond=n_cond, kappa=kappa)
  end function arima_fit

  function arimax_fit(x, p, d, q, xreg, include_mean, seasonal_p, seasonal_d, &
      seasonal_q, period, max_iterations, tolerance, xtransf, transfer, io, &
      fixed, init, method, transform_pars, n_cond, kappa) result(fit)
    real(dp), intent(in) :: x(:), xreg(:,:)
    integer, intent(in) :: p, d, q
    logical, intent(in), optional :: include_mean, transform_pars
    integer, intent(in), optional :: seasonal_p, seasonal_d, seasonal_q, period
    integer, intent(in), optional :: max_iterations, io(:), n_cond
    real(dp), intent(in), optional :: tolerance, xtransf(:,:), fixed(:), init(:), kappa
    type(transfer_spec), intent(in), optional :: transfer(:)
    character(len=*), intent(in), optional :: method
    type(arimax_result) :: fit

    real(dp), allocatable :: yd(:), par_work(:), par_actual(:), init_actual(:)
    real(dp), allocatable :: fixed_actual(:), free(:), free_scale(:), parscale_actual(:)
    real(dp), allocatable :: hess(:,:), hinv(:,:)
    real(dp), allocatable :: jac(:,:), covfree(:,:), covwork(:,:), residual_d(:)
    real(dp), allocatable :: fitted_d(:), beta0(:), xinit(:,:)
    real(dp), allocatable :: reg_rotation(:,:), par_output(:), cov_output(:,:)
    logical, allocatable :: free_mask(:)
    integer, allocatable :: free_index(:)
    integer :: sp, sd, sq, per, nbase, nio, nreg, ntrans, narma, ncoef
    integer :: maxit, iterations, istat, i, j, off, nfree, ndiff, css_drop, ml_nused
    real(dp) :: tol, obj, sigma2, logdet, diffuse_kappa
    logical :: inc, do_transform, has_missing, use_ml, css_ml, rotate_regression
    character(len=6) :: fit_method

    sp = 0
    sd = 0
    sq = 0
    per = 1
    if (present(seasonal_p)) sp = seasonal_p
    if (present(seasonal_d)) sd = seasonal_d
    if (present(seasonal_q)) sq = seasonal_q
    if (present(period)) per = max(1, period)

    inc = .true.
    if (present(include_mean)) inc = include_mean
    do_transform = .true.
    if (present(transform_pars)) do_transform = transform_pars
    maxit = 3000
    if (present(max_iterations)) maxit = max_iterations
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = tolerance
    diffuse_kappa = 1.0e6_dp
    if (present(kappa)) diffuse_kappa = kappa
    if (diffuse_kappa <= 0.0_dp) then
      fit%status = 11
      return
    end if

    fit_method = 'CSS'
    if (present(method)) then
      select case (method(1:1))
      case ('M','m')
        fit_method = 'ML'
      case ('C','c')
        if (len_trim(method) >= 5) then
          if (index(method, 'ML') > 0 .or. index(method, 'ml') > 0) then
            fit_method = 'CSS-ML'
          else
            fit_method = 'CSS'
          end if
        else
          fit_method = 'CSS'
        end if
      case default
        fit%status = 10
        return
      end select
    end if

    fit%p = p
    fit%d = d
    fit%q = q
    fit%seasonal_p = sp
    fit%seasonal_d = sd
    fit%seasonal_q = sq
    fit%period = per
    fit%include_mean = inc
    fit%method = fit_method
    fit%n_xreg = size(xreg,2)
    fit%n_io = 0
    if (present(io)) fit%n_io = size(io)
    fit%n_transfer = 0
    if (present(transfer)) fit%n_transfer = size(transfer)

    if (size(xreg,1) /= size(x) .or. any([p,d,q,sp,sd,sq] < 0)) then
      fit%status = 1
      return
    end if
    if (per < 1 .or. (sp + sd + sq > 0 .and. per == 1)) then
      fit%status = 2
      return
    end if
    if (present(io)) then
      if (any(io < 1) .or. any(io > size(x))) then
        fit%status = 3
        return
      end if
    end if
    if (present(transfer)) then
      if (.not. present(xtransf)) then
        fit%status = 4
        return
      end if
      if (size(xtransf,1) /= size(x) .or. size(xtransf,2) /= size(transfer)) then
        fit%status = 4
        return
      end if
      do i = 1, size(transfer)
        if (transfer(i)%ar_order < 0 .or. transfer(i)%ma_order < 0) then
          fit%status = 4
          return
        end if
      end do
    else if (present(xtransf)) then
      if (size(xtransf,2) > 0) then
        fit%status = 4
        return
      end if
    end if

    yd = difference_series(x, d)
    yd = seasonal_difference(yd, sd, per)
    ndiff = size(x) - size(yd)
    if (size(yd) < 4) then
      fit%status = 5
      return
    end if

    nbase = size(xreg,2)
    nio = 0
    if (present(io)) nio = size(io)
    nreg = nbase + nio + merge(1,0,inc .and. d + sd == 0)
    narma = p + q + sp + sq
    ntrans = 0
    if (present(transfer)) then
      do i = 1, size(transfer)
        ntrans = ntrans + transfer(i)%ar_order + transfer(i)%ma_order + 1
      end do
    end if
    ncoef = narma + nreg + ntrans

    css_drop = p + sp*per
    if (present(n_cond)) css_drop = max(css_drop, max(0,n_cond))
    fit%n_cond = css_drop + ndiff
    if (size(yd) <= css_drop + 2) then
      fit%status = 6
      return
    end if

    has_missing = any(.not. ieee_is_finite(x))
    if (size(xreg,2) > 0) has_missing = has_missing .or. any(.not. ieee_is_finite(xreg))
    if (present(xtransf)) then
      if (size(xtransf,2) > 0) has_missing = has_missing .or. any(.not. ieee_is_finite(xtransf))
    end if
    if (fit_method == 'CSS-ML' .and. has_missing) fit_method = 'ML'
    if (fit_method == 'CSS' .and. has_missing) then
      fit%status = 7
      return
    end if
    fit%method = fit_method
    if (fit_method == 'ML') fit%n_cond = 0

    allocate(init_actual(ncoef), fixed_actual(ncoef), free_mask(ncoef))
    init_actual = 0.0_dp
    fixed_actual = 0.0_dp
    free_mask = .true.
    rotate_regression = .false.

    ! R/TSA determine whether xreg may be SVD-rotated from the fixed mask
    ! before constructing regression starting values.
    if (present(fixed)) then
      if (size(fixed) /= ncoef) then
        fit%status = 9
        return
      end if
      do i = 1, ncoef
        if (.not. ieee_is_nan(fixed(i))) then
          free_mask(i) = .false.
          fixed_actual(i) = fixed(i)
        end if
      end do
    end if

    allocate(reg_rotation(max(0,nreg),max(0,nreg)))
    if (nreg > 0) then
      reg_rotation = 0.0_dp
      do i = 1, nreg
        reg_rotation(i,i) = 1.0_dp
      end do
      if (nreg > 1 .and. all(free_mask(narma+1:narma+nreg))) then
        call regression_svd_rotation(reg_rotation, istat)
        rotate_regression = istat == 0
        istat = 0
      end if
    end if

    if (nreg > 0) then
      call initial_regression_matrix(xinit, init_actual)
      if (size(xinit,2) == nreg) then
        allocate(beta0(nreg))
        call complete_case_least_squares(xinit,yd,beta0,istat)
        if (istat == 0) init_actual(narma+1:narma+nreg) = beta0
      end if
    end if

    if (present(init)) then
      if (size(init) /= ncoef) then
        fit%status = 8
        return
      end if
      do i = 1, ncoef
        if (.not. ieee_is_nan(init(i))) init_actual(i) = init(i)
      end do
    end if

    do i = 1, ncoef
      if (.not. free_mask(i)) init_actual(i) = fixed_actual(i)
    end do

    if (do_transform) then
      if (p > 0) then
        if (any(.not. free_mask(1:p))) do_transform = .false.
      end if
      if (sp > 0) then
        off = p + q
        if (any(.not. free_mask(off+1:off+sp))) do_transform = .false.
      end if
    end if
    if (do_transform .and. fit_method == 'ML') call invert_free_ma_blocks(init_actual)

    allocate(parscale_actual(ncoef))
    call parameter_scales(parscale_actual)
    par_work = work_from_actual(init_actual, do_transform)
    nfree = count(free_mask)
    allocate(free_index(nfree), free(nfree), free_scale(nfree))
    j = 0
    do i = 1, ncoef
      if (free_mask(i)) then
        j = j + 1
        free_index(j) = i
        free(j) = par_work(i)
        free_scale(j) = parscale_actual(i)
      end if
    end do

    css_ml = fit_method == 'CSS-ML'
    use_ml = fit_method == 'ML'
    iterations = 0
    istat = 0
    obj = 0.0_dp

    if (nfree > 0) then
      if (css_ml) then
        use_ml = .false.
        call bfgs(objective_free, free, obj, iterations, istat, &
          max_iterations=maxit, reltol=tol, parscale=free_scale)
        if (istat == 0) then
          call merge_free(free, par_work)
          par_actual = actual_from_work(par_work, do_transform)
          if (do_transform) call invert_free_ma_blocks(par_actual)
          par_work = work_from_actual(par_actual, do_transform)
          do j = 1, nfree
            free(j) = par_work(free_index(j))
          end do
        end if
        use_ml = .true.
      end if
      if (istat == 0) then
        call bfgs(objective_free, free, obj, iterations, istat, &
          max_iterations=maxit, reltol=tol, parscale=free_scale)
      end if
      call merge_free(free, par_work)
    else
      obj = objective_full(par_work)
    end if

    par_actual = actual_from_work(par_work, do_transform)
    if (use_ml .and. do_transform) then
      call invert_free_ma_blocks(par_actual)
      par_work = work_from_actual(par_actual,do_transform)
      do j = 1, nfree
        free(j) = par_work(free_index(j))
      end do
    end if
    call evaluate_model(par_actual, residual_d, fitted_d, sigma2, logdet, use_ml, ml_nused)

    par_output = par_actual
    if (rotate_regression) then
      par_output(narma+1:narma+nreg) = &
        matmul(reg_rotation,par_actual(narma+1:narma+nreg))
    end if

    allocate(fit%coefficients(ncoef), fit%covariance(ncoef,ncoef), &
      fit%residuals(size(x)), fit%fitted(size(x)), fit%estimated(ncoef), &
      fit%series(size(x)))
    fit%coefficients = par_output
    fit%estimated = free_mask
    fit%series = x
    fit%covariance = 0.0_dp
    if (use_ml) then
      fit%residuals = residual_d
    else
      fit%residuals = 0.0_dp
      fit%residuals(ndiff+1:) = residual_d
    end if
    fit%fitted = x - fit%residuals

    allocate(fit%ar(p), fit%ma(q), fit%sar(sp), fit%sma(sq), &
      fit%regression(nreg), fit%transfer(ntrans))
    off = 0
    if (p > 0) then
      fit%ar = par_output(1:p)
      off = p
    end if
    if (q > 0) then
      fit%ma = par_output(off+1:off+q)
      off = off + q
    end if
    if (sp > 0) then
      fit%sar = par_output(off+1:off+sp)
      off = off + sp
    end if
    if (sq > 0) then
      fit%sma = par_output(off+1:off+sq)
      off = off + sq
    end if
    if (nreg > 0) fit%regression = par_output(narma+1:narma+nreg)
    if (ntrans > 0) fit%transfer = par_output(narma+nreg+1:)

    fit%sigma2 = sigma2
    if (use_ml) then
      fit%loglik = concentrated_loglik(sigma2, logdet, ml_nused)
      fit%aic = -2.0_dp*fit%loglik + 2.0_dp*real(nfree+1,dp)
    else
      fit%loglik = css_loglik(residual_d, css_drop, sigma2)
      fit%aic = huge(1.0_dp)
    end if
    fit%iterations = iterations
    fit%status = istat

    if (nfree > 0 .and. istat == 0) then
      allocate(hess(nfree,nfree), hinv(nfree,nfree), covfree(nfree,nfree))
      call optim_hessian(objective_free, free, hess, istat, parscale=free_scale)
      if (istat == 0) call invert_matrix_lu(hess, hinv, istat)
      if (istat == 0) then
        covfree = hinv / real(max(1,merge(ml_nused,count_finite(residual_d),use_ml)),dp)
        allocate(covwork(ncoef,ncoef), jac(ncoef,ncoef))
        covwork = 0.0_dp
        do i = 1, nfree
          do j = 1, nfree
            covwork(free_index(i),free_index(j)) = covfree(i,j)
          end do
        end do
        call parameter_jacobian(par_work, do_transform, jac)
        fit%covariance = matmul(jac, matmul(covwork, transpose(jac)))
        if (rotate_regression) then
          allocate(cov_output(ncoef,ncoef))
          call rotate_covariance_to_original(fit%covariance,cov_output)
          fit%covariance = cov_output
        end if
      end if
    end if
    if (fit%status == 0 .and. istat /= 0) fit%status = istat

  contains

    subroutine regression_svd_rotation(v, status)
      ! R stats::arima/TSA rotate multiple all-free regression columns by
      ! the right singular vectors of complete-case xreg before fitting.
      real(dp), intent(out) :: v(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: raw(:,:), complete(:,:), sval(:)
      logical, allocatable :: keep(:)
      integer :: row, nr, k, col

      status = 0
      v = 0.0_dp
      do k = 1, size(v,1)
        v(k,k) = 1.0_dp
      end do
      if (nreg <= 1) return

      allocate(raw(size(x),nreg),keep(size(x)))
      raw = 0.0_dp
      col = 0
      if (inc .and. d+sd == 0) then
        col = col+1
        raw(:,col) = 1.0_dp
      end if
      do k = 1, nbase
        col = col+1
        raw(:,col) = xreg(:,k)
      end do
      do k = 1, nio
        col = col+1
        raw(io(k),col) = 1.0_dp
      end do
      if (col /= nreg) then
        status = 1
        return
      end if

      keep = .true.
      do k = 1, nreg
        keep = keep .and. ieee_is_finite(raw(:,k))
      end do
      nr = count(keep)
      if (nr <= 0) then
        status = 2
        return
      end if
      allocate(complete(nr,nreg))
      row = 0
      do k = 1, size(x)
        if (keep(k)) then
          row = row+1
          complete(row,:) = raw(k,:)
        end if
      end do
      allocate(sval(nreg))
      call right_singular_vectors(complete,v,sval,status)
    end subroutine regression_svd_rotation

    subroutine rotate_covariance_to_original(covin,covout)
      real(dp), intent(in) :: covin(:,:)
      real(dp), intent(out) :: covout(:,:)
      real(dp), allocatable :: a(:,:)
      integer :: k

      allocate(a(ncoef,ncoef))
      a = 0.0_dp
      do k = 1, ncoef
        a(k,k) = 1.0_dp
      end do
      a(narma+1:narma+nreg,narma+1:narma+nreg) = reg_rotation
      covout = matmul(a,matmul(covin,transpose(a)))
    end subroutine rotate_covariance_to_original

    subroutine apply_regression_rotation(mat)
      real(dp), intent(inout) :: mat(:,:)
      real(dp), allocatable :: tmp(:,:)
      if (.not. rotate_regression .or. size(mat,2) /= nreg) return
      tmp = matmul(mat,reg_rotation)
      mat = tmp
    end subroutine apply_regression_rotation

    subroutine invert_free_ma_blocks(par)
      real(dp), intent(inout) :: par(:)
      real(dp), allocatable :: tmp(:)
      integer :: k0

      if (q > 0) then
        k0 = p+1
        if (all(free_mask(k0:k0+q-1))) then
          tmp = ma_invert(par(k0:k0+q-1))
          par(k0:k0+q-1) = tmp
        end if
      end if
      if (sq > 0) then
        k0 = p+q+sp+1
        if (all(free_mask(k0:k0+sq-1))) then
          tmp = ma_invert(par(k0:k0+sq-1))
          par(k0:k0+sq-1) = tmp
        end if
      end if
    end subroutine invert_free_ma_blocks

    subroutine parameter_scales(scales)
      ! TSA arimax() defaults optim.control$parscale to one for ARMA
      ! parameters, 10*OLS standard errors for regression coefficients, and
      ! 10 times the simple-regression slope SE for each transfer block.
      real(dp), intent(out) :: scales(:)
      real(dp), allocatable :: mat(:,:), xx(:,:), yy(:), beta(:), cov(:,:)
      logical, allocatable :: keep(:)
      real(dp) :: se
      integer :: cols, col, k, irow, nr, st, pos, pa, qa, slope_col

      scales = 1.0_dp
      if (size(scales) /= ncoef) return

      if (nreg > 0) then
        if (nio == 0 .and. ntrans == 0) then
          ! Plain TSA arimax() delegates to stats::arima(). Match current R:
          ! compute xreg starting-value SEs on the differenced response and
          ! differenced rotated regressors, with no extra intercept column.
          call regression_matrix(init_actual,mat)
          cols = nreg
          slope_col = 1
          allocate(keep(size(yd)))
          keep = ieee_is_finite(yd)
          do k = 1, cols
            keep = keep .and. ieee_is_finite(mat(:,k))
          end do
          nr = count(keep)
          st = 1
          if (nr > cols) then
            allocate(xx(nr,cols),yy(nr),beta(cols),cov(cols,cols))
            irow = 0
            do k = 1, size(yd)
              if (keep(k)) then
                irow = irow+1
                xx(irow,:) = mat(k,:)
                yy(irow) = yd(k)
              end if
            end do
            call least_squares(xx,yy,beta,covariance=cov,status=st)
            if (st == 0) then
              do k = 1, nreg
                se = sqrt(max(0.0_dp,cov(k,k)))
                if (ieee_is_finite(se) .and. se > tiny(1.0_dp)) &
                  scales(narma+k) = 10.0_dp*se
              end do
            end if
            deallocate(xx,yy,beta,cov)
          end if
          deallocate(mat,keep)

          ! R falls back to the undifferenced regression for degenerate
          ! differenced designs. Preserve that old-code compatibility path.
          if (st /= 0) then
            allocate(mat(size(x),nreg))
            mat = 0.0_dp
            col = 0
            if (inc .and. d+sd == 0) then
              col = col+1
              mat(:,col) = 1.0_dp
            end if
            do k = 1, nbase
              col = col+1
              mat(:,col) = xreg(:,k)
            end do
            call apply_regression_rotation(mat)
            allocate(keep(size(x)))
            keep = ieee_is_finite(x)
            do k = 1, nreg
              keep = keep .and. ieee_is_finite(mat(:,k))
            end do
            nr = count(keep)
            if (nr > nreg) then
              allocate(xx(nr,nreg),yy(nr),beta(nreg),cov(nreg,nreg))
              irow = 0
              do k = 1, size(x)
                if (keep(k)) then
                  irow = irow+1
                  xx(irow,:) = mat(k,:)
                  yy(irow) = x(k)
                end if
              end do
              call least_squares(xx,yy,beta,covariance=cov,status=st)
              if (st == 0) then
                do k = 1, nreg
                  se = sqrt(max(0.0_dp,cov(k,k)))
                  if (ieee_is_finite(se) .and. se > tiny(1.0_dp)) &
                    scales(narma+k) = 10.0_dp*se
                end do
              end if
              deallocate(xx,yy,beta,cov)
            end if
            deallocate(mat,keep)
          end if
        else
          ! TSA's specialized transfer/IO branch uses the raw-series lm()
          ! construction from arimax.R for optim.control$parscale.
          if (inc .and. d+sd == 0) then
            cols = nreg
            slope_col = 1
            allocate(mat(size(x),cols))
            col = 1
            mat(:,1) = 1.0_dp
          else
            cols = nreg+1
            slope_col = 2
            allocate(mat(size(x),cols))
            mat(:,1) = 1.0_dp
            col = 1
          end if
          do k = 1, nbase
            col = col+1
            mat(:,col) = xreg(:,k)
          end do
          do k = 1, nio
            col = col+1
            mat(:,col) = 0.0_dp
            mat(io(k),col) = 1.0_dp
          end do
          if (rotate_regression) then
            if (cols == nreg) then
              call apply_regression_rotation(mat)
            else if (cols == nreg+1) then
              mat(:,2:) = matmul(mat(:,2:),reg_rotation)
            end if
          end if
          allocate(keep(size(x)))
          keep = ieee_is_finite(x)
          do k=1,cols
            keep = keep .and. ieee_is_finite(mat(:,k))
          end do
          nr = count(keep)
          if (nr > cols) then
            allocate(xx(nr,cols),yy(nr),beta(cols),cov(cols,cols))
            irow=0
            do k=1,size(x)
              if (keep(k)) then
                irow=irow+1
                xx(irow,:)=mat(k,:)
                yy(irow)=x(k)
              end if
            end do
            call least_squares(xx,yy,beta,covariance=cov,status=st)
            if (st == 0) then
              do k=1,nreg
                se=sqrt(max(0.0_dp,cov(slope_col+k-1,slope_col+k-1)))
                if (ieee_is_finite(se) .and. se > tiny(1.0_dp)) &
                  scales(narma+k)=10.0_dp*se
              end do
            end if
            deallocate(xx,yy,beta,cov)
          end if
          deallocate(mat,keep)
        end if
      end if

      if (ntrans > 0) then
        pos=narma+nreg
        do k=1,size(transfer)
          pa=transfer(k)%ar_order
          qa=transfer(k)%ma_order+1
          allocate(keep(size(x)))
          keep=ieee_is_finite(x) .and. ieee_is_finite(xtransf(:,k))
          nr=count(keep)
          se=0.0_dp
          if (nr > 2) then
            allocate(xx(nr,2),yy(nr),beta(2),cov(2,2))
            irow=0
            do col=1,size(x)
              if (keep(col)) then
                irow=irow+1
                xx(irow,1)=1.0_dp
                xx(irow,2)=xtransf(col,k)
                yy(irow)=x(col)
              end if
            end do
            call least_squares(xx,yy,beta,covariance=cov,status=st)
            if (st == 0) se=sqrt(max(0.0_dp,cov(2,2)))
            deallocate(xx,yy,beta,cov)
          end if
          if (ieee_is_finite(se) .and. se > tiny(1.0_dp)) then
            scales(pos+1:pos+pa+qa)=10.0_dp*se
          end if
          pos=pos+pa+qa
          deallocate(keep)
        end do
      end if
    end subroutine parameter_scales

    subroutine initial_regression_matrix(mat, initial)
      real(dp), allocatable, intent(out) :: mat(:,:)
      real(dp), intent(in) :: initial(:)
      real(dp), allocatable :: dummy(:)
      if (nreg == 0) then
        allocate(mat(size(yd),0))
        return
      end if
      call regression_matrix(initial, mat)
      if (ntrans > 0) then
        call transfer_effect_differenced(initial, dummy)
      end if
    end subroutine initial_regression_matrix

    subroutine complete_case_least_squares(mat,response,beta,status)
      real(dp), intent(in) :: mat(:,:), response(:)
      real(dp), intent(out) :: beta(:)
      integer, intent(out) :: status
      logical, allocatable :: keep(:)
      real(dp), allocatable :: xx(:,:), yy(:)
      integer :: i, j, nr

      status = 0
      beta = 0.0_dp
      if (size(mat,1) /= size(response) .or. size(mat,2) /= size(beta)) then
        status = 1
        return
      end if
      allocate(keep(size(response)))
      keep = ieee_is_finite(response)
      do j = 1, size(mat,2)
        keep = keep .and. ieee_is_finite(mat(:,j))
      end do
      nr = count(keep)
      if (nr <= size(mat,2)) then
        status = 2
        return
      end if
      allocate(xx(nr,size(mat,2)),yy(nr))
      j = 0
      do i = 1, size(response)
        if (keep(i)) then
          j = j+1
          xx(j,:) = mat(i,:)
          yy(j) = response(i)
        end if
      end do
      call least_squares(xx,yy,beta,status=status)
    end subroutine complete_case_least_squares

    subroutine merge_free(fv, full)
      real(dp), intent(in) :: fv(:)
      real(dp), intent(inout) :: full(:)
      integer :: k
      do k = 1, size(fv)
        full(free_index(k)) = fv(k)
      end do
      do k = 1, ncoef
        if (.not. free_mask(k)) full(k) = fixed_actual(k)
      end do
    end subroutine merge_free

    function objective_free(fv) result(v)
      real(dp), intent(in) :: fv(:)
      real(dp) :: v
      real(dp), allocatable :: full(:)
      allocate(full(ncoef))
      full = par_work
      call merge_free(fv, full)
      v = objective_full(full)
    end function objective_free

    function objective_full(work) result(v)
      real(dp), intent(in) :: work(:)
      real(dp) :: v, sig2, ldet
      real(dp), allocatable :: actual(:), rr(:)
      integer :: nobs, nused_eval

      if (size(work) > 0) then
        if (maxval(abs(work), mask=ieee_is_finite(work)) > 50.0_dp) then
          v = huge(1.0_dp)/100.0_dp
          return
        end if
      end if
      actual = actual_from_work(work, do_transform)
      if (.not. all_transfer_stable(actual)) then
        v = huge(1.0_dp)/100.0_dp
        return
      end if
      call evaluate_model(actual, rr, sigma2_out=sig2, logdet_out=ldet, ml=use_ml, nused_out=nused_eval)
      nobs = merge(nused_eval,count_finite(rr),use_ml)
      if (nobs <= 0 .or. .not. ieee_is_finite(sig2) .or. sig2 <= 0.0_dp) then
        v = huge(1.0_dp)/100.0_dp
      else if (use_ml) then
        v = 0.5_dp*(log(sig2) + ldet/real(nobs,dp))
      else
        v = 0.5_dp*log(sig2)
      end if
    end function objective_full

    subroutine evaluate_model(actual, rr, ff, sigma2_out, logdet_out, ml, nused_out)
      real(dp), intent(in) :: actual(:)
      real(dp), allocatable, intent(out) :: rr(:)
      real(dp), allocatable, intent(out), optional :: ff(:)
      real(dp), intent(out) :: sigma2_out, logdet_out
      logical, intent(in) :: ml
      integer, intent(out), optional :: nused_out
      real(dp), allocatable :: regmat(:,:), transv(:), z(:), delta(:)
      real(dp), allocatable :: fullar(:), fullma(:)
      integer :: nobs, start_css, nu

      call full_arma(actual, fullar, fullma)
      if (ml) then
        allocate(z(size(x)))
        z = x
        if (nreg > 0) then
          call regression_matrix_ml(actual, regmat)
          z = z - matmul(regmat, actual(narma+1:narma+nreg))
        end if
        if (ntrans > 0) then
          call transfer_effect_raw(actual, transv)
          z = z - transv
        end if
        delta = differencing_recursion(d,sd,per)
        call arima_diffuse_innovations(z, fullar, fullma, delta, diffuse_kappa, &
          rr, sigma2_out, logdet_out, nu)
        if (present(nused_out)) nused_out = nu
        if (present(ff)) ff = x - rr
      else
        allocate(z(size(yd)))
        z = yd
        if (nreg > 0) then
          call regression_matrix(actual, regmat)
          z = z - matmul(regmat, actual(narma+1:narma+nreg))
        end if
        if (ntrans > 0) then
          call transfer_effect_differenced(actual, transv)
          z = z - transv
        end if
        call css_innovations(z, fullar, fullma, rr)
        start_css = min(size(rr), css_drop) + 1
        nobs = count(ieee_is_finite(rr(start_css:)))
        if (nobs > 0) then
          sigma2_out = sum(pack(rr(start_css:)**2, &
            ieee_is_finite(rr(start_css:)))) / real(nobs,dp)
        else
          sigma2_out = huge(1.0_dp)
        end if
        logdet_out = 0.0_dp
        if (present(nused_out)) nused_out = nobs
        if (present(ff)) ff = yd - rr
      end if
    end subroutine evaluate_model

    subroutine regression_matrix_ml(actual, mat)
      real(dp), intent(in) :: actual(:)
      real(dp), allocatable, intent(out) :: mat(:,:)
      real(dp), allocatable :: raw(:), tmp(:), arfull(:), mafull(:), delta(:)
      integer :: col, k

      allocate(mat(size(x),nreg))
      if (nreg == 0) return
      mat = 0.0_dp
      col = 0
      if (inc .and. d + sd == 0) then
        col = col + 1
        mat(:,col) = 1.0_dp
      end if
      do k = 1, nbase
        col = col + 1
        mat(:,col) = xreg(:,k)
      end do
      if (nio > 0) then
        call full_arma(actual, arfull, mafull)
        delta = differencing_recursion(d,sd,per)
        do k = 1, nio
          col = col + 1
          allocate(raw(size(x)))
          raw = 0.0_dp
          raw(io(k)) = 1.0_dp
          call apply_arma_filter(raw, arfull, mafull, delta, tmp)
          mat(:,col) = tmp
          deallocate(raw)
        end do
      end if
      call apply_regression_rotation(mat)
    end subroutine regression_matrix_ml

    subroutine transfer_effect_raw(actual, effect)
      real(dp), intent(in) :: actual(:)
      real(dp), allocatable, intent(out) :: effect(:)
      real(dp), allocatable :: one(:), phi(:), theta(:)
      integer :: k, pos, pa, qa

      allocate(effect(size(x)))
      effect = 0.0_dp
      pos = narma + nreg
      do k = 1, size(transfer)
        pa = transfer(k)%ar_order
        qa = transfer(k)%ma_order + 1
        allocate(phi(pa), theta(qa))
        if (pa > 0) phi = actual(pos+1:pos+pa)
        if (qa > 0) theta = actual(pos+pa+1:pos+pa+qa)
        call transfer_filter(xtransf(:,k), phi, theta, one)
        effect = effect + one
        pos = pos + pa + qa
        deallocate(phi,theta,one)
      end do
    end subroutine transfer_effect_raw

    subroutine regression_matrix(actual, mat)
      real(dp), intent(in) :: actual(:)
      real(dp), allocatable, intent(out) :: mat(:,:)
      real(dp), allocatable :: raw(:), tmp(:), arfull(:), mafull(:), delta(:)
      integer :: col, k

      allocate(mat(size(yd),nreg))
      if (nreg == 0) return
      mat = 0.0_dp
      col = 0
      if (inc .and. d + sd == 0) then
        col = col + 1
        mat(:,col) = 1.0_dp
      end if
      do k = 1, nbase
        col = col + 1
        tmp = difference_series(xreg(:,k),d)
        tmp = seasonal_difference(tmp,sd,per)
        mat(:,col) = tmp
      end do
      if (nio > 0) then
        call full_arma(actual, arfull, mafull)
        delta = differencing_recursion(d,sd,per)
        do k = 1, nio
          col = col + 1
          allocate(raw(size(x)))
          raw = 0.0_dp
          raw(io(k)) = 1.0_dp
          call apply_arma_filter(raw, arfull, mafull, delta, tmp)
          raw = difference_series(tmp,d)
          raw = seasonal_difference(raw,sd,per)
          mat(:,col) = raw
          deallocate(raw)
        end do
      end if
      call apply_regression_rotation(mat)
    end subroutine regression_matrix

    subroutine transfer_effect_differenced(actual, effectd)
      real(dp), intent(in) :: actual(:)
      real(dp), allocatable, intent(out) :: effectd(:)
      real(dp), allocatable :: raw(:), one(:), tmp(:), phi(:), theta(:)
      integer :: k, pos, pa, qa

      allocate(raw(size(x)))
      raw = 0.0_dp
      pos = narma + nreg
      do k = 1, size(transfer)
        pa = transfer(k)%ar_order
        qa = transfer(k)%ma_order + 1
        allocate(phi(pa), theta(qa))
        if (pa > 0) phi = actual(pos+1:pos+pa)
        if (qa > 0) theta = actual(pos+pa+1:pos+pa+qa)
        call transfer_filter(xtransf(:,k), phi, theta, one)
        raw = raw + one
        pos = pos + pa + qa
        deallocate(phi,theta,one)
      end do
      tmp = difference_series(raw,d)
      effectd = seasonal_difference(tmp,sd,per)
    end subroutine transfer_effect_differenced

    subroutine full_arma(actual, arcoef, macoef)
      real(dp), intent(in) :: actual(:)
      real(dp), allocatable, intent(out) :: arcoef(:), macoef(:)
      real(dp), allocatable :: rp(:), spoly(:), mp(:), smpoly(:), convp(:)
      integer :: k, ofs

      allocate(rp(p+1), mp(q+1), spoly(sp*per+1), smpoly(sq*per+1))
      rp = 0.0_dp
      mp = 0.0_dp
      spoly = 0.0_dp
      smpoly = 0.0_dp
      rp(1) = 1.0_dp
      mp(1) = 1.0_dp
      spoly(1) = 1.0_dp
      smpoly(1) = 1.0_dp
      ofs = 0
      if (p > 0) then
        rp(2:) = -actual(1:p)
        ofs = p
      end if
      if (q > 0) then
        mp(2:) = actual(ofs+1:ofs+q)
        ofs = ofs + q
      end if
      if (sp > 0) then
        do k = 1, sp
          spoly(k*per+1) = -actual(ofs+k)
        end do
        ofs = ofs + sp
      end if
      if (sq > 0) then
        do k = 1, sq
          smpoly(k*per+1) = actual(ofs+k)
        end do
      end if
      convp = polynomial_convolution(rp,spoly)
      allocate(arcoef(max(0,size(convp)-1)))
      if (size(arcoef) > 0) arcoef = -convp(2:)
      convp = polynomial_convolution(mp,smpoly)
      allocate(macoef(max(0,size(convp)-1)))
      if (size(macoef) > 0) macoef = convp(2:)
    end subroutine full_arma

    function actual_from_work(work, trans) result(actual)
      real(dp), intent(in) :: work(:)
      logical, intent(in) :: trans
      real(dp), allocatable :: actual(:)
      integer :: ofs

      allocate(actual(size(work)))
      actual = work
      if (.not. trans) return
      ofs = 0
      if (p > 0) then
        actual(1:p) = levinson(work(1:p))
        ofs = p
      end if
      ofs = ofs + q
      if (sp > 0) then
        actual(ofs+1:ofs+sp) = levinson(work(ofs+1:ofs+sp))
      end if
    end function actual_from_work

    function work_from_actual(actual, trans) result(work)
      real(dp), intent(in) :: actual(:)
      logical, intent(in) :: trans
      real(dp), allocatable :: work(:)
      integer :: ofs

      allocate(work(size(actual)))
      work = actual
      if (.not. trans) return
      ofs = 0
      if (p > 0) then
        work(1:p) = ar_to_unconstrained(actual(1:p))
        ofs = p
      end if
      ofs = ofs + q
      if (sp > 0) then
        work(ofs+1:ofs+sp) = ar_to_unconstrained(actual(ofs+1:ofs+sp))
      end if
    end function work_from_actual

    subroutine parameter_jacobian(work, trans, jmat)
      ! Analytic derivative of TSA's Levinson/PACF stationarity transform.
      ! jmat(i,j) = d(actual_i)/d(work_j). TSA's difpar() stores the
      ! transpose of this matrix before applying t(A) H^-1 A.
      real(dp), intent(in) :: work(:)
      logical, intent(in) :: trans
      real(dp), intent(out) :: jmat(:,:)
      real(dp), allocatable :: block(:,:)
      integer :: k, ofs

      jmat = 0.0_dp
      do k = 1, size(work)
        jmat(k,k) = 1.0_dp
      end do
      if (.not. trans) return
      ofs = 0
      if (p > 0) then
        allocate(block(p,p))
        call levinson_jacobian(work(1:p),block)
        jmat(1:p,1:p) = block
        deallocate(block)
        ofs = p
      end if
      ofs = ofs+q
      if (sp > 0) then
        allocate(block(sp,sp))
        call levinson_jacobian(work(ofs+1:ofs+sp),block)
        jmat(ofs+1:ofs+sp,ofs+1:ofs+sp) = block
      end if
    end subroutine parameter_jacobian

    subroutine levinson_jacobian(u,jmat)
      real(dp), intent(in) :: u(:)
      real(dp), intent(out) :: jmat(:,:)
      real(dp), allocatable :: oldphi(:), newphi(:), oldj(:,:), newj(:,:)
      real(dp), allocatable :: pacf(:), dpacf(:)
      integer :: n, k, i, m

      n = size(u)
      jmat = 0.0_dp
      if (n == 0) return
      allocate(pacf(n),dpacf(n),oldphi(1),oldj(1,n))
      pacf = tanh(u)
      dpacf = 1.0_dp-pacf*pacf
      oldphi(1) = pacf(1)
      oldj = 0.0_dp
      oldj(1,1) = dpacf(1)
      if (n == 1) then
        jmat(1,1) = dpacf(1)
        return
      end if

      do k = 2, n
        allocate(newphi(k),newj(k,n))
        newphi = 0.0_dp
        newj = 0.0_dp
        newphi(k) = pacf(k)
        newj(k,k) = dpacf(k)
        do i = 1, k-1
          newphi(i) = oldphi(i)-pacf(k)*oldphi(k-i)
          do m = 1, k-1
            newj(i,m) = oldj(i,m)-pacf(k)*oldj(k-i,m)
          end do
          newj(i,k) = -dpacf(k)*oldphi(k-i)
        end do
        call move_alloc(newphi,oldphi)
        call move_alloc(newj,oldj)
      end do
      jmat = oldj
    end subroutine levinson_jacobian

    logical function all_transfer_stable(actual) result(ok)
      real(dp), intent(in) :: actual(:)
      integer :: k, pos, pa, qa
      ok = .true.
      if (.not. present(transfer)) return
      pos = narma + nreg
      do k = 1, size(transfer)
        pa = transfer(k)%ar_order
        qa = transfer(k)%ma_order + 1
        if (pa > 0) then
          if (sum(abs(actual(pos+1:pos+pa))) > 20.0_dp) then
            ok = .false.
            return
          end if
        end if
        pos = pos + pa + qa
      end do
    end function all_transfer_stable

  end function arimax_fit

  pure function levinson(u) result(phi)
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: phi(:), old(:), newv(:)
    real(dp) :: mult
    integer :: k, j

    allocate(phi(size(u)))
    if (size(u) == 0) return
    phi = tanh(u)
    if (size(u) == 1) return
    allocate(old(1))
    old(1) = phi(1)
    do k = 2, size(u)
      allocate(newv(k))
      mult = tanh(u(k))
      newv(k) = mult
      do j = 1, k-1
        newv(j) = old(j) - mult*old(k-j)
      end do
      call move_alloc(newv,old)
    end do
    phi = old
  end function levinson

  pure function ar_to_unconstrained(phi) result(u)
    real(dp), intent(in) :: phi(:)
    real(dp), allocatable :: u(:), cur(:), prev(:)
    real(dp) :: kappa, den
    integer :: m, j

    allocate(u(size(phi)),cur(size(phi)))
    if (size(phi) == 0) return
    cur = phi
    do m = size(phi), 2, -1
      kappa = max(-1.0_dp+1.0e-10_dp, min(1.0_dp-1.0e-10_dp, cur(m)))
      u(m) = atanh(kappa)
      den = max(1.0e-12_dp, 1.0_dp-kappa*kappa)
      allocate(prev(m-1))
      do j = 1, m-1
        prev(j) = (cur(j)+kappa*cur(m-j))/den
      end do
      deallocate(cur)
      allocate(cur(m-1))
      cur = prev
      deallocate(prev)
    end do
    u(1) = atanh(max(-1.0_dp+1.0e-10_dp,min(1.0_dp-1.0e-10_dp,cur(1))))
  end function ar_to_unconstrained

  subroutine transfer_filter(x, phi, theta, y)
    real(dp), intent(in) :: x(:), phi(:), theta(:)
    real(dp), allocatable, intent(out) :: y(:)
    real(dp), allocatable :: tmp(:)
    real(dp), parameter :: safe_limit = 1.0e100_dp
    real(dp) :: value, term, prev
    integer :: i, j, k
    logical :: bad

    allocate(tmp(size(x)),y(size(x)))
    bad = any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(phi)) .or. &
      any(.not. ieee_is_finite(theta))
    if (bad) then
      y = safe_limit
      return
    end if

    ! R's recursive filter may produce Inf for trial transfer parameters.
    ! Under strict IEEE traps, detect the same unusable optimization region
    ! before overflow and return a large finite regressor instead.
    tmp = 0.0_dp
    do i = 1, size(x)
      value = x(i)
      if (abs(value) > safe_limit) then
        bad = .true.
        exit
      end if
      do j = 1, size(phi)
        k = i-j
        if (k < 1) cycle
        prev = tmp(k)
        if (abs(prev) > tiny(1.0_dp)) then
          if (abs(phi(j)) > safe_limit/abs(prev)) then
            bad = .true.
            exit
          end if
        end if
        term = phi(j)*prev
        if (abs(value) > safe_limit-abs(term)) then
          bad = .true.
          exit
        end if
        value = value+term
      end do
      if (bad) exit
      tmp(i) = value
    end do
    if (bad) then
      y = safe_limit
      return
    end if

    if (size(theta) == 0) then
      y = tmp
      return
    end if
    y = 0.0_dp
    do i = 1, size(tmp)
      value = 0.0_dp
      do j = 1, min(size(theta),i)
        prev = tmp(i-j+1)
        if (abs(prev) > tiny(1.0_dp)) then
          if (abs(theta(j)) > safe_limit/abs(prev)) then
            bad = .true.
            exit
          end if
        end if
        term = theta(j)*prev
        if (abs(value) > safe_limit-abs(term)) then
          bad = .true.
          exit
        end if
        value = value+term
      end do
      if (bad) exit
      y(i)=value
    end do
    if (bad) y = safe_limit
  end subroutine transfer_filter

  subroutine io_regressor(n, point, ar, ma, d, seasonal_d, period, y, status)
    integer, intent(in) :: n, point, d, seasonal_d, period
    real(dp), intent(in) :: ar(:), ma(:)
    real(dp), allocatable, intent(out) :: y(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: pulse(:), delta(:)

    if (point < 1 .or. point > n) then
      allocate(y(0))
      if (present(status)) status = 1
      return
    end if
    allocate(pulse(n))
    pulse = 0.0_dp
    pulse(point) = 1.0_dp
    delta = differencing_recursion(d,seasonal_d,max(1,period))
    call apply_arma_filter(pulse,ar,ma,delta,y)
    if (present(status)) status = 0
  end subroutine io_regressor

  function differencing_recursion(d, seasonal_d, period) result(delta)
    integer, intent(in) :: d, seasonal_d, period
    real(dp), allocatable :: delta(:), poly(:), fac(:), tmp(:)
    integer :: k, s

    s = max(1,period)
    allocate(poly(1))
    poly = 1.0_dp
    do k = 1, d
      allocate(fac(2))
      fac = [1.0_dp,-1.0_dp]
      tmp = polynomial_convolution(poly,fac)
      call move_alloc(tmp,poly)
      deallocate(fac)
    end do
    do k = 1, seasonal_d
      allocate(fac(s+1))
      fac = 0.0_dp
      fac(1) = 1.0_dp
      fac(s+1) = -1.0_dp
      tmp = polynomial_convolution(poly,fac)
      call move_alloc(tmp,poly)
      deallocate(fac)
    end do
    allocate(delta(max(0,size(poly)-1)))
    if (size(delta) > 0) delta = -poly(2:)
  end function differencing_recursion

  subroutine apply_arma_filter(x, ar, ma, delta, y)
    real(dp), intent(in) :: x(:), ar(:), ma(:), delta(:)
    real(dp), allocatable, intent(out) :: y(:)
    real(dp), allocatable :: t1(:), t2(:), theta(:)

    allocate(t1(size(x)),t2(size(x)))
    if (size(delta) > 0) then
      call recursive_filter(x,delta,t1)
    else
      t1 = x
    end if
    if (size(ar) > 0) then
      call recursive_filter(t1,ar,t2)
    else
      t2 = t1
    end if
    allocate(theta(size(ma)+1))
    theta(1) = 1.0_dp
    if (size(ma) > 0) theta(2:) = ma
    allocate(y(size(x)))
    call convolution_filter(t2,theta,y)
  end subroutine apply_arma_filter

  subroutine css_innovations(z, ar, ma, residuals)
    real(dp), intent(in) :: z(:), ar(:), ma(:)
    real(dp), allocatable, intent(out) :: residuals(:)
    real(dp) :: pred
    integer :: i,j

    allocate(residuals(size(z)))
    residuals = 0.0_dp
    do i = 1, size(z)
      if (.not. ieee_is_finite(z(i))) then
        residuals(i) = z(i)
        cycle
      end if
      pred = 0.0_dp
      do j = 1, min(size(ar),i-1)
        if (ieee_is_finite(z(i-j))) pred = pred + ar(j)*z(i-j)
      end do
      do j = 1, min(size(ma),i-1)
        if (ieee_is_finite(residuals(i-j))) pred = pred + ma(j)*residuals(i-j)
      end do
      residuals(i) = z(i)-pred
    end do
  end subroutine css_innovations

  subroutine arima_diffuse_innovations(z, ar, ma, delta, kappa, residuals, sigma2, logdet, nused)
    ! R stats::arima-style innovations likelihood for ARIMA models.
    ! The ARMA part is stationary and the integrated components receive
    ! diffuse variance kappa. Innovations with gain >= 1e4 are updated but
    ! excluded from the concentrated likelihood, matching C_ARIMA_Like.
    real(dp), intent(in) :: z(:), ar(:), ma(:), delta(:), kappa
    real(dp), allocatable, intent(out) :: residuals(:)
    real(dp), intent(out) :: sigma2, logdet
    integer, intent(out) :: nused
    real(dp), allocatable :: tmat(:,:), zvec(:), rvec(:), a(:), anew(:)
    real(dp), allocatable :: pmat(:,:), pnew(:,:), pstat(:,:), mvec(:)
    real(dp) :: innov, gain, ssq
    integer :: p, q, r, dint, rd, i, k, info

    p = size(ar)
    q = size(ma)
    r = max(p,q+1)
    dint = size(delta)
    rd = r+dint
    allocate(tmat(rd,rd),zvec(rd),rvec(rd),a(rd),anew(rd), &
      pmat(rd,rd),pnew(rd,rd),mvec(rd),residuals(size(z)))
    tmat = 0.0_dp
    zvec = 0.0_dp
    rvec = 0.0_dp
    a = 0.0_dp
    pmat = 0.0_dp
    residuals = 0.0_dp

    zvec(1) = 1.0_dp
    if (dint > 0) zvec(r+1:rd) = delta
    if (p > 0) tmat(1:p,1) = ar
    if (r > 1) then
      do i = 2, r
        tmat(i-1,i) = 1.0_dp
      end do
    end if
    if (dint > 0) then
      tmat(r+1,:) = zvec
      if (dint > 1) then
        do i = 2, dint
          tmat(r+i,r+i-1) = 1.0_dp
        end do
      end if
    end if
    rvec(1) = 1.0_dp
    if (q > 0) rvec(2:q+1) = ma

    allocate(pstat(r,r))
    call stationary_state_covariance_as154(ar,ma,pstat,info)
    if (info /= 0) then
      ! The AS154 packed system grows as O(r^4) in storage.  For unusually
      ! large orders or an allocation failure, use the algebraically
      ! equivalent discrete Lyapunov solver rather than failing the fit.
      call stationary_state_covariance(tmat(1:r,1:r), rvec(1:r), pstat, info)
    end if
    if (info /= 0) then
      sigma2 = huge(1.0_dp)
      logdet = huge(1.0_dp)
      nused = 0
      residuals = z
      return
    end if
    pmat(1:r,1:r) = pstat
    if (dint > 0) then
      do i = 1, dint
        pmat(r+i,r+i) = kappa
      end do
    end if

    ssq = 0.0_dp
    logdet = 0.0_dp
    nused = 0
    do k = 1, size(z)
      anew = matmul(tmat,a)
      ! R starts with Pn directly at the first observation; subsequent
      ! observations receive one state prediction plus process variance.
      if (k == 1) then
        pnew = pmat
      else
        pnew = matmul(tmat,matmul(pmat,transpose(tmat))) + outer_product(rvec,rvec)
        pnew = 0.5_dp*(pnew+transpose(pnew))
      end if

      if (ieee_is_finite(z(k))) then
        innov = z(k)-dot_product(zvec,anew)
        mvec = matmul(pnew,zvec)
        gain = dot_product(zvec,mvec)
        if (.not. ieee_is_finite(gain) .or. gain <= tiny(1.0_dp)) then
          sigma2 = huge(1.0_dp)
          logdet = huge(1.0_dp)
          nused = 0
          residuals = z
          return
        end if
        residuals(k) = innov/sqrt(gain)
        if (gain < 1.0e4_dp) then
          nused = nused+1
          ssq = ssq + innov*innov/gain
          logdet = logdet + log(gain)
        end if
        a = anew + mvec*(innov/gain)
        pmat = pnew - outer_product(mvec,mvec)/gain
        pmat = 0.5_dp*(pmat+transpose(pmat))
      else
        residuals(k) = z(k)
        a = anew
        pmat = pnew
      end if
    end do
    if (nused > 0) then
      sigma2 = ssq/real(nused,dp)
    else
      sigma2 = huge(1.0_dp)
    end if
  end subroutine arima_diffuse_innovations

  subroutine stationary_state_covariance_as154(ar, ma, pmat, status)
    ! R's default Gardner1980/AS154 initialization (C_getQ0), translated
    ! from src/library/stats/src/arima.c.  A packed upper-triangular system
    ! is accumulated with the AS154 inclu2 update and then back-substituted.
    real(dp), intent(in) :: ar(:), ma(:)
    real(dp), intent(out) :: pmat(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: xnext(:), xrow(:), rbar(:), thetab(:), vpack(:), ppack(:)
    real(dp) :: vj, vi, phij, phii, ynext, bi
    integer :: p, q, r, np, nrbar, npr, npr1
    integer :: indi, indj, indn, i, j, ithisr, ind, ind1, ind2, im, jm, astat

    p = size(ar)
    q = size(ma)
    r = max(p,q+1)
    status = 0
    if (size(pmat,1) /= r .or. size(pmat,2) /= r) then
      status = 1
      return
    end if
    pmat = 0.0_dp
    if (r == 1) then
      if (p == 0) then
        pmat(1,1) = 1.0_dp
      else if (abs(ar(1)) < 1.0_dp) then
        pmat(1,1) = 1.0_dp/(1.0_dp-ar(1)*ar(1))
      else
        status = 2
      end if
      return
    end if

    ! R permits r <= 350, but nrbar is quadratic in np = r(r+1)/2.
    ! Avoid multi-hundred-MB/GB temporary allocations in the standalone
    ! Fortran port; the caller falls back to an equivalent Lyapunov solve.
    if (r > 80) then
      status = 4
      return
    end if
    np = r*(r+1)/2
    nrbar = np*(np-1)/2
    allocate(xnext(0:np-1),xrow(0:np-1),thetab(0:np-1), &
      vpack(0:np-1),ppack(0:np-1),rbar(0:nrbar-1),stat=astat)
    if (astat /= 0) then
      status = 4
      return
    end if

    ind = 0
    do j = 0, r-1
      vj = 0.0_dp
      if (j == 0) then
        vj = 1.0_dp
      else if (j <= q) then
        vj = ma(max(1,j))
      end if
      do i = j, r-1
        vi = 0.0_dp
        if (i == 0) then
          vi = 1.0_dp
        else if (i <= q) then
          vi = ma(max(1,i))
        end if
        vpack(ind) = vi*vj
        ind = ind+1
      end do
    end do

    if (p > 0) then
      rbar = 0.0_dp
      ppack = 0.0_dp
      thetab = 0.0_dp
      xnext = 0.0_dp
      ind = 0
      ind1 = -1
      npr = np-r
      npr1 = npr+1
      indj = npr
      ind2 = npr-1
      do j = 0, r-1
        phij = 0.0_dp
        if (j < p) phij = ar(j+1)
        xnext(indj) = 0.0_dp
        indj = indj+1
        indi = npr1+j
        do i = j, r-1
          ynext = vpack(ind)
          ind = ind+1
          phii = 0.0_dp
          if (i < p) phii = ar(i+1)
          if (j /= r-1) then
            xnext(indj) = -phii
            if (i /= r-1) then
              xnext(indi) = xnext(indi)-phij
              ind1 = ind1+1
              xnext(ind1) = -1.0_dp
            end if
          end if
          xnext(npr) = -phii*phij
          ind2 = ind2+1
          if (ind2 >= np) ind2 = 0
          xnext(ind2) = xnext(ind2)+1.0_dp
          call inclu2_as154(np,xnext,xrow,ynext,ppack,rbar,thetab)
          xnext(ind2) = 0.0_dp
          if (i /= r-1) then
            xnext(indi) = 0.0_dp
            indi = indi+1
            xnext(ind1) = 0.0_dp
          end if
        end do
      end do
      ithisr = nrbar-1
      im = np-1
      do i = 0, np-1
        bi = thetab(im)
        jm = np-1
        do j = 0, i-1
          bi = bi-rbar(ithisr)*ppack(jm)
          ithisr = ithisr-1
          jm = jm-1
        end do
        ppack(im) = bi
        im = im-1
      end do

      ind = npr
      do i = 0, r-1
        xnext(i) = ppack(ind)
        ind = ind+1
      end do
      ind = np-1
      ind1 = npr-1
      do i = 0, npr-1
        ppack(ind) = ppack(ind1)
        ind = ind-1
        ind1 = ind1-1
      end do
      do i = 0, r-1
        ppack(i) = xnext(i)
      end do
    else
      indn = np
      ind = np
      do i = 0, r-1
        do j = 0, i
          ind = ind-1
          ppack(ind) = vpack(ind)
          if (j /= 0) then
            indn = indn-1
            ppack(ind) = ppack(ind)+ppack(indn)
          end if
        end do
      end do
    end if

    ! C_getQ0 stores the packed solution in the first np slots of the
    ! column-major r-by-r result before unpacking the remaining triangle.
    pmat = 0.0_dp
    do ind = 0, np-1
      i = mod(ind,r)
      j = ind/r
      pmat(i+1,j+1) = ppack(ind)
    end do
    ind = np
    do i = r-1, 1, -1
      do j = r-1, i, -1
        ind = ind-1
        pmat(j+1,i+1) = ppack(ind)
      end do
    end do
    do i = 0, r-2
      do j = i+1, r-1
        pmat(i+1,j+1) = pmat(j+1,i+1)
      end do
    end do
    if (any(.not. ieee_is_finite(pmat))) status = 3
  end subroutine stationary_state_covariance_as154

  subroutine inclu2_as154(np, xnext, xrow, ynext, dpack, rbar, thetab)
    integer, intent(in) :: np
    real(dp), intent(in) :: xnext(0:), ynext
    real(dp), intent(out) :: xrow(0:)
    real(dp), intent(inout) :: dpack(0:), rbar(0:), thetab(0:)
    real(dp) :: ywork, cbar, sbar, di, xi, xk, rbthis, dpi
    integer :: i, k, ithisr

    xrow(0:np-1) = xnext(0:np-1)
    ywork = ynext
    ithisr = 0
    do i = 0, np-1
      if (abs(xrow(i)) > 0.0_dp) then
        xi = xrow(i)
        di = dpack(i)
        dpi = di+xi*xi
        dpack(i) = dpi
        cbar = di/dpi
        sbar = xi/dpi
        do k = i+1, np-1
          xk = xrow(k)
          rbthis = rbar(ithisr)
          xrow(k) = xk-xi*rbthis
          rbar(ithisr) = cbar*rbthis+sbar*xk
          ithisr = ithisr+1
        end do
        xk = ywork
        ywork = xk-xi*thetab(i)
        thetab(i) = cbar*thetab(i)+sbar*xk
        if (di <= 0.0_dp) return
      else
        ithisr = ithisr+np-i-1
      end if
    end do
  end subroutine inclu2_as154

  subroutine stationary_state_covariance(tmat, rvec, pmat, status)
    ! Solve P = T P T' + R R'.  For moderate state dimension use the
    ! vectorized Lyapunov system; fall back to fixed-point iteration for
    ! larger orders to avoid a very large dense linear system.
    real(dp), intent(in) :: tmat(:,:), rvec(:)
    real(dp), intent(out) :: pmat(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: amat(:,:), b(:), sol(:), pnext(:,:), vmat(:,:)
    real(dp) :: diff, scale
    integer :: r, i, j, k, l, eq, var, iter

    r = size(rvec)
    status = 0
    pmat = 0.0_dp
    if (size(tmat,1) /= r .or. size(tmat,2) /= r .or. &
        size(pmat,1) /= r .or. size(pmat,2) /= r) then
      status = 1
      return
    end if
    vmat = outer_product(rvec,rvec)
    if (r == 1) then
      if (abs(1.0_dp-tmat(1,1)*tmat(1,1)) <= 100.0_dp*epsilon(1.0_dp)) then
        status = 2
      else
        pmat(1,1) = vmat(1,1)/(1.0_dp-tmat(1,1)*tmat(1,1))
      end if
      return
    end if

    if (r <= 12) then
      allocate(amat(r*r,r*r),b(r*r),sol(r*r))
      amat = 0.0_dp
      do i = 1, r
        do j = 1, r
          eq = (i-1)*r+j
          amat(eq,eq) = 1.0_dp
          b(eq) = vmat(i,j)
          do k = 1, r
            do l = 1, r
              var = (k-1)*r+l
              amat(eq,var) = amat(eq,var) - tmat(i,k)*tmat(j,l)
            end do
          end do
        end do
      end do
      call solve_linear(amat,b,sol,status)
      if (status == 0) then
        do i = 1, r
          do j = 1, r
            pmat(i,j) = sol((i-1)*r+j)
          end do
        end do
        pmat = 0.5_dp*(pmat+transpose(pmat))
        return
      end if
      status = 0
    end if

    allocate(pnext(r,r))
    pmat = 0.0_dp
    do iter = 1, 100000
      pnext = matmul(tmat,matmul(pmat,transpose(tmat))) + vmat
      diff = maxval(abs(pnext-pmat))
      scale = max(1.0_dp,maxval(abs(pnext)))
      pmat = pnext
      if (diff <= 5.0e-14_dp*scale) exit
    end do
    if (iter > 100000) status = 3
    pmat = 0.5_dp*(pmat+transpose(pmat))
  end subroutine stationary_state_covariance

  pure function outer_product(a,b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a),size(b))
    integer :: i,j
    do j = 1, size(b)
      do i = 1, size(a)
        c(i,j) = a(i)*b(j)
      end do
    end do
  end function outer_product

  pure integer function count_finite(x) result(n)
    real(dp), intent(in) :: x(:)
    n = count(ieee_is_finite(x))
  end function count_finite

  pure real(dp) function concentrated_loglik(sigma2, logdet, nobs) result(v)
    real(dp), intent(in) :: sigma2, logdet
    integer, intent(in) :: nobs
    if (nobs <= 0 .or. sigma2 <= 0.0_dp) then
      v = -huge(1.0_dp)
    else
      v = -0.5_dp*(real(nobs,dp)*(log(2.0_dp*acos(-1.0_dp)) + &
        1.0_dp + log(sigma2)) + logdet)
    end if
  end function concentrated_loglik

  pure real(dp) function css_loglik(residuals, drop, sigma2) result(v)
    real(dp), intent(in) :: residuals(:), sigma2
    integer, intent(in) :: drop
    integer :: first, nobs
    first = min(size(residuals),max(0,drop)) + 1
    nobs = count(ieee_is_finite(residuals(first:)))
    if (nobs <= 0 .or. sigma2 <= 0.0_dp) then
      v = -huge(1.0_dp)
    else
      v = -0.5_dp*real(nobs,dp)*(log(2.0_dp*acos(-1.0_dp)*sigma2)+1.0_dp)
    end if
  end function css_loglik

  subroutine arima_sim(ar, ma, d, n, x, sigma, mean, ntrans, init)
    real(dp), intent(in) :: ar(:), ma(:)
    integer, intent(in) :: d, n
    real(dp), allocatable, intent(out) :: x(:)
    real(dp), intent(in), optional :: sigma, mean, init(:)
    integer, intent(in), optional :: ntrans
    real(dp) :: sig, mu, pred
    real(dp), allocatable :: work(:), innovations(:), out(:)
    integer :: nt, p, q, total, i, j, k, start

    sig = 1.0_dp
    mu = 0.0_dp
    nt = 100
    if (present(sigma)) sig = sigma
    if (present(mean)) mu = mean
    if (present(ntrans)) nt = ntrans
    p = size(ar)
    q = size(ma)
    start = max(1,p+d)
    total = n + nt + start
    allocate(work(total),innovations(total))
    work = mu
    innovations = 0.0_dp
    if (present(init)) then
      work(1:min(size(init),total)) = init(1:min(size(init),total))
    end if
    do i = 1, total
      innovations(i) = sig*random_normal()
    end do
    do i = start+1, total
      pred = mu
      do j = 1, min(p,i-1)
        pred = pred + ar(j)*(work(i-j)-mu)
      end do
      do j = 1, min(q,i-1)
        pred = pred + ma(j)*innovations(i-j)
      end do
      work(i) = pred + innovations(i)
    end do
    allocate(out(n))
    out = work(total-n+1:)
    do k = 1, d
      do i = 2, n
        out(i) = out(i) + out(i-1)
      end do
    end do
    call move_alloc(out,x)
  end subroutine arima_sim

  function ma_invert(ma) result(out)
    real(dp), intent(in) :: ma(:)
    real(dp), allocatable :: out(:)
    complex(dp), allocatable :: roots(:), poly(:), next(:)
    complex(dp) :: z, denom, delta
    real(dp) :: radius, twopi, err
    integer :: q0, i, j, iter
    logical :: converged, changed

    allocate(out(size(ma)))
    out = ma
    q0 = 0
    do i = 1, size(ma)
      if (abs(ma(i)) > epsilon(1.0_dp)) q0 = i
    end do
    if (q0 == 0) return
    if (q0 == 1) then
      if (abs(ma(1)) > 1.0_dp) out(1) = 1.0_dp/ma(1)
      return
    end if

    allocate(roots(q0))
    radius = 1.0_dp+maxval(abs(ma(:q0)))/max(abs(ma(q0)),sqrt(epsilon(1.0_dp)))
    radius = max(1.1_dp,min(radius,10.0_dp))
    twopi = 2.0_dp*acos(-1.0_dp)
    do i = 1, q0
      roots(i) = cmplx(radius*cos(twopi*(real(i,dp)-0.37_dp)/real(q0,dp)), &
        radius*sin(twopi*(real(i,dp)-0.37_dp)/real(q0,dp)),kind=dp)
    end do
    converged = .false.
    do iter = 1, 500
      err = 0.0_dp
      do i = 1, q0
        z = roots(i)
        denom = cmplx(1.0_dp,0.0_dp,kind=dp)
        do j = 1, q0
          if (j /= i) denom = denom*(z-roots(j))
        end do
        if (abs(denom) <= tiny(1.0_dp)) cycle
        delta = poly_value(ma(:q0),z)/denom
        roots(i) = z-delta
        err = max(err,abs(delta))
      end do
      if (err < 1.0e-12_dp) then
        converged = .true.
        exit
      end if
    end do
    if (.not. converged) return

    changed = .false.
    do i = 1, q0
      if (abs(roots(i)) < 1.0_dp) then
        roots(i) = 1.0_dp/roots(i)
        changed = .true.
      end if
    end do
    if (.not. changed) return

    allocate(poly(1))
    poly(1) = cmplx(1.0_dp,0.0_dp,kind=dp)
    do i = 1, q0
      allocate(next(size(poly)+1))
      next = cmplx(0.0_dp,0.0_dp,kind=dp)
      next(:size(poly)) = next(:size(poly))+poly
      next(2:) = next(2:)-poly/roots(i)
      call move_alloc(next,poly)
    end do
    out(:q0) = real(poly(2:),kind=dp)
  contains
    function poly_value(coef,x) result(value)
      real(dp), intent(in) :: coef(:)
      complex(dp), intent(in) :: x
      complex(dp) :: value
      integer :: k
      value = cmplx(coef(size(coef)),0.0_dp,kind=dp)
      do k = size(coef)-1, 1, -1
        value = value*x+cmplx(coef(k),0.0_dp,kind=dp)
      end do
      value = value*x+cmplx(1.0_dp,0.0_dp,kind=dp)
    end function poly_value
  end function ma_invert

  subroutine arima_bootstrap_sample(fit, sample, normal, cond_boot, ntrans, &
      init, status)
    type(arimax_result), intent(in) :: fit
    real(dp), allocatable, intent(out) :: sample(:)
    logical, intent(in), optional :: normal, cond_boot
    integer, intent(in), optional :: ntrans
    real(dp), intent(in), optional :: init(:)
    integer, intent(out) :: status

    logical :: norm, cond
    integer :: nt, n, pp, qq, total, i, j, idx, nkeep
    real(dp) :: intercept
    real(dp), allocatable :: arint(:), initial(:), z(:), noise(:), work(:)
    real(dp), allocatable :: finite_resid(:), bootall(:)

    norm = .true.
    cond = .false.
    if (present(normal)) norm = normal
    if (present(cond_boot)) cond = cond_boot
    nt = 100
    if (present(ntrans)) nt = max(0,ntrans)
    if (cond) nt = 0

    n = size(fit%residuals)
    arint = bootstrap_integrated_ar(fit%ar,fit%d)
    pp = size(arint)
    qq = size(fit%ma)
    if (n < 1 .or. pp > n) then
      status = 1
      allocate(sample(0))
      return
    end if

    allocate(initial(pp))
    if (pp > 0) then
      if (present(init)) then
        if (size(init) < pp) then
          status = 2
          allocate(sample(0))
          return
        end if
        initial = init(:pp)
      else if (allocated(fit%series) .and. size(fit%series) >= pp) then
        initial = fit%series(:pp)
      else
        initial = 0.0_dp
      end if
    end if

    intercept = 0.0_dp
    if (fit%include_mean .and. allocated(fit%regression)) then
      if (size(fit%regression) > 0) intercept = fit%regression(1)
    end if

    total = n+nt
    allocate(z(total+qq),noise(total),bootall(total),work(pp+total))
    if (norm) then
      do i = 1, size(z)
        z(i) = sqrt(max(fit%sigma2,0.0_dp))*random_normal()
      end do
    else
      finite_resid = pack(fit%residuals,ieee_is_finite(fit%residuals))
      if (size(finite_resid) == 0) then
        status = 3
        allocate(sample(0))
        return
      end if
      call sample_with_replacement(finite_resid,z)
    end if

    ! TSA arima.boot uses stats::filter(..., filter=ma, method='convolution')
    ! without a leading unit coefficient.  Preserve that source behavior.
    if (qq > 0) then
      noise = 0.0_dp
      do i = 1, total
        do j = 1, qq
          noise(i) = noise(i)+fit%ma(j)*z(qq+i-j+1)
        end do
      end do
    else
      noise = z(1:total)
    end if

    work = 0.0_dp
    if (pp > 0) work(:pp) = initial-intercept
    do i = 1, total
      idx = pp+i
      work(idx) = noise(i)
      do j = 1, pp
        work(idx) = work(idx)+arint(j)*work(idx-j)
      end do
      bootall(i) = work(idx)+intercept
    end do

    allocate(sample(n))
    sample = bootall(nt+1:nt+n)
    if (cond .and. pp > 0) then
      nkeep = min(pp,n)
      sample(:nkeep) = initial(:nkeep)
      if (nkeep < n) sample(nkeep+1:) = bootall(1:n-nkeep)
    end if
    status = 0
  end subroutine arima_bootstrap_sample

  subroutine arima_bootstrap(fit, b, boot, normal, cond_boot, ntrans, init, status)
    type(arimax_result), intent(in) :: fit
    integer, intent(in) :: b
    type(bootstrap_result), intent(out) :: boot
    logical, intent(in), optional :: normal, cond_boot
    integer, intent(in), optional :: ntrans
    real(dp), intent(in), optional :: init(:)
    integer, intent(out) :: status

    integer :: attempts, max_attempts, st, pcoef
    real(dp), allocatable :: sim(:)
    type(arimax_result) :: bf

    if (b < 1) then
      status = 1
      allocate(boot%coefficients(0,0),boot%sigma2(0))
      boot%successful = 0
      return
    end if
    pcoef = size(fit%coefficients)
    allocate(boot%coefficients(b,pcoef),boot%sigma2(b))
    boot%coefficients = 0.0_dp
    boot%sigma2 = 0.0_dp
    boot%successful = 0
    attempts = 0
    max_attempts = max(20,20*b)

    do while (boot%successful < b .and. attempts < max_attempts)
      attempts = attempts+1
      call arima_bootstrap_sample(fit,sim,normal,cond_boot,ntrans,init,st)
      if (st /= 0) then
        status = st
        return
      end if
      bf = arima_fit(sim,fit%p,fit%d,fit%q,fit%include_mean,method='CSS-ML')
      if (bf%status /= 0 .or. size(bf%coefficients) /= pcoef) cycle
      boot%successful = boot%successful+1
      boot%coefficients(boot%successful,:) = bf%coefficients
      boot%sigma2(boot%successful) = bf%sigma2
    end do
    status = merge(0,1,boot%successful == b)
  end subroutine arima_bootstrap

  function bootstrap_integrated_ar(ar,dord) result(out)
    real(dp), intent(in) :: ar(:)
    integer, intent(in) :: dord
    real(dp), allocatable :: out(:), poly(:), fac(:), tmp(:)
    integer :: k

    allocate(poly(size(ar)+1))
    poly(1) = 1.0_dp
    if (size(ar) > 0) poly(2:) = -ar
    do k = 1, dord
      fac = [1.0_dp,-1.0_dp]
      tmp = polynomial_convolution(poly,fac)
      call move_alloc(tmp,poly)
    end do
    allocate(out(size(poly)-1))
    out = -poly(2:)
  end function bootstrap_integrated_ar

end module tsa_arimax
