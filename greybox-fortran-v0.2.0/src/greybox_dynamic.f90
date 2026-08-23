module greybox_dynamic
    use greybox_kinds, only: dp
    use greybox_regression, only: alm_model, alm_fit, alm_predict
    implicit none
    private
    public :: alm_occurrence_model, alm_dynamic_model
    public :: alm_fit_occurrence, alm_fit_arima_errors

    type :: alm_occurrence_model
        character(len=16) :: distribution = 'dnorm'
        character(len=16) :: occurrence_link = 'plogis'
        type(alm_model) :: occurrence
        type(alm_model) :: positive
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: point_loglik(:)
        real(dp) :: loglik = 0.0_dp
        real(dp) :: aic = huge(1.0_dp)
        integer :: nobs = 0
        integer :: nparam = 0
        logical :: converged = .false.
    end type alm_occurrence_model

    type :: alm_dynamic_model
        character(len=16) :: distribution = 'dnorm'
        integer :: orders(3) = 0
        type(alm_model) :: conditional
        real(dp), allocatable :: exogenous_beta(:)
        real(dp), allocatable :: ar(:)
        real(dp), allocatable :: ma(:)
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residuals(:)
        logical :: converged = .false.
    end type alm_dynamic_model

contains

    subroutine alm_fit_occurrence(x,y,distribution,occurrence_link,model,max_iter,tol)
        real(dp), intent(in) :: x(:,:), y(:)
        character(len=*), intent(in) :: distribution, occurrence_link
        type(alm_occurrence_model), intent(out) :: model
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: occ_y(:), xp(:,:), yp(:), pocc(:), posmean(:)
        integer, allocatable :: ids(:)
        integer :: n, i, j, np

        n = size(y)
        if (size(x,1) /= n) error stop 'alm_fit_occurrence: incompatible x and y'
        allocate(occ_y(n),pocc(n),posmean(n))
        where(y > 0.0_dp)
            occ_y = 1.0_dp
        elsewhere
            occ_y = 0.0_dp
        end where
        call alm_fit(x,occ_y,occurrence_link,model%occurrence,max_iter=max_iter,tol=tol)
        np = count(y > 0.0_dp)
        if (np < max(2,size(x,2))) then
            model%converged = .false.
            return
        end if
        allocate(ids(np),xp(np,size(x,2)),yp(np))
        j = 0
        do i=1,n
            if (y(i) <= 0.0_dp) cycle
            j = j + 1
            ids(j)=i; xp(j,:)=x(i,:); yp(j)=y(i)
        end do
        call alm_fit(xp,yp,distribution,model%positive,max_iter=max_iter,tol=tol)
        if (.not.model%occurrence%converged .or. .not.model%positive%converged) then
            model%converged=.false.; return
        end if
        call alm_predict(model%occurrence,x,pocc)
        call alm_predict(model%positive,x,posmean)
        allocate(model%fitted(n),model%point_loglik(n))
        model%fitted = pocc*posmean
        model%point_loglik = log(max(1.0_dp-pocc,1.0e-15_dp))
        j=0
        do i=1,n
            if (y(i) > 0.0_dp) then
                j=j+1
                model%point_loglik(i)=log(max(pocc(i),1.0e-15_dp))+model%positive%point_loglik(j)
            end if
        end do
        model%distribution=distribution
        model%occurrence_link=occurrence_link
        model%loglik=sum(model%point_loglik)
        model%nobs=n
        model%nparam=model%occurrence%nparam+model%positive%nparam
        model%aic=2.0_dp*real(model%nparam,dp)-2.0_dp*model%loglik
        model%converged=.true.
    end subroutine alm_fit_occurrence

    subroutine alm_fit_arima_errors(x,y,distribution,orders,model,max_iter,tol,n_outer)
        real(dp), intent(in) :: x(:,:), y(:)
        character(len=*), intent(in) :: distribution
        integer, intent(in) :: orders(3)
        type(alm_dynamic_model), intent(out) :: model
        integer, intent(in), optional :: max_iter, n_outer
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: yd(:), xd(:,:), residual(:), design(:,:), target(:)
        type(alm_model) :: fit
        integer :: pord, dord, qord, lag, n, px, m, i, j, iter, outer

        pord=max(0,orders(1)); dord=max(0,orders(2)); qord=max(0,orders(3))
        if (size(x,1) /= size(y)) error stop 'alm_fit_arima_errors: incompatible x and y'
        call difference_data(x,y,dord,xd,yd)
        n=size(yd); px=size(xd,2); lag=max(pord,qord)
        if (n-lag <= px+pord+qord) then
            model%converged=.false.; return
        end if
        allocate(residual(n)); residual=0.0_dp
        outer=8; if(present(n_outer)) outer=max(1,n_outer)
        do iter=1,outer
            m=n-lag
            if (allocated(design)) deallocate(design,target)
            allocate(design(m,px+pord+qord),target(m))
            do i=1,m
                design(i,1:px)=xd(lag+i,:)
                do j=1,pord
                    design(i,px+j)=yd(lag+i-j)
                end do
                do j=1,qord
                    design(i,px+pord+j)=residual(lag+i-j)
                end do
                target(i)=yd(lag+i)
            end do
            call alm_fit(design,target,distribution,fit,max_iter=max_iter,tol=tol)
            if (.not.fit%converged) then
                model%converged=.false.; return
            end if
            residual=0.0_dp
            residual(lag+1:n)=fit%residuals
        end do
        model%distribution=distribution
        model%orders=orders
        model%conditional=fit
        allocate(model%exogenous_beta(px),model%ar(pord),model%ma(qord), &
            model%fitted(n),model%residuals(n))
        model%exogenous_beta=fit%beta(1:px)
        if(pord>0) model%ar=fit%beta(px+1:px+pord)
        if(qord>0) model%ma=fit%beta(px+pord+1:px+pord+qord)
        model%fitted=0.0_dp
        model%residuals=0.0_dp
        model%fitted(lag+1:n)=fit%fitted
        model%residuals(lag+1:n)=fit%residuals
        model%converged=.true.
    end subroutine alm_fit_arima_errors

    subroutine difference_data(x,y,d,xd,yd)
        real(dp), intent(in) :: x(:,:), y(:)
        integer, intent(in) :: d
        real(dp), allocatable, intent(out) :: xd(:,:), yd(:)
        real(dp), allocatable :: xt(:,:), yt(:), xn(:,:), yn(:)
        integer :: k, n
        allocate(xt(size(x,1),size(x,2)),yt(size(y)))
        xt=x; yt=y
        do k=1,d
            n=size(yt)-1
            allocate(xn(n,size(x,2)),yn(n))
            xn=xt(2:,:)-xt(:size(xt,1)-1,:)
            yn=yt(2:)-yt(:size(yt)-1)
            call move_alloc(xn,xt); call move_alloc(yn,yt)
        end do
        call move_alloc(xt,xd); call move_alloc(yt,yd)
    end subroutine difference_data

end module greybox_dynamic
