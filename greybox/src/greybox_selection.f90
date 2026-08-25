module greybox_selection
    use greybox_kinds, only: dp
    use greybox_regression, only: alm_model, alm_fit, alm_predict, point_aic, point_aicc, point_bic, point_bicc
    implicit none
    private
    public :: calm_model, stepwise_fit, calm_fit, calm_predict
    public :: rolling_origin_alm, recursive_lm
    public :: calm_dynamic_model, lm_dynamic_fit

    type :: calm_model
        character(len=16) :: distribution = 'dnorm'
        character(len=8) :: criterion = 'AICc'
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: scale_beta(:)
        real(dp), allocatable :: inclusion_probability(:)
        real(dp), allocatable :: model_weight(:)
        real(dp) :: best_ic = huge(1.0_dp)
        integer :: nmodels = 0
        logical :: converged = .false.
    end type calm_model

    type :: calm_dynamic_model
        character(len=16) :: distribution = 'dnorm'
        character(len=8) :: criterion = 'AICc'
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: coefficients_dynamic(:,:)
        real(dp), allocatable :: importance(:,:)
        real(dp), allocatable :: weights(:,:)
        real(dp), allocatable :: point_ic(:,:)
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residuals(:)
        real(dp) :: span = 0.67_dp
        integer :: nmodels = 0
        logical :: converged = .false.
    end type calm_dynamic_model

contains

    subroutine stepwise_fit(x,y,distribution,model,selected,criterion,max_iter,tol)
        real(dp), intent(in) :: x(:,:), y(:)
        character(len=*), intent(in) :: distribution
        type(alm_model), intent(out) :: model
        logical, intent(out) :: selected(:)
        character(len=*), intent(in), optional :: criterion
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        type(alm_model) :: current, candidate
        logical, allocatable :: trial(:)
        real(dp), allocatable :: xsub(:,:)
        real(dp) :: current_ic, candidate_ic, best_ic, tolerance
        integer :: p, j, best_j, itmax
        character(len=8) :: icname

        p = size(x,2)
        if (size(selected) /= p) error stop 'stepwise_fit: selected has wrong size'
        icname = 'AICc'
        if (present(criterion)) icname = adjustl(criterion)
        itmax = 800
        if (present(max_iter)) itmax = max_iter
        tolerance = 1.0e-6_dp
        if (present(tol)) tolerance = tol

        selected = .false.
        selected(1) = .true.
        allocate(trial(p))
        call make_subset(x,selected,xsub)
        call alm_fit(xsub,y,distribution,current,max_iter=itmax,tol=tolerance)
        current_ic = model_ic(current,icname)

        do
            best_ic = current_ic
            best_j = 0
            do j=2,p
                if (selected(j)) cycle
                trial = selected
                trial(j) = .true.
                call make_subset(x,trial,xsub)
                call alm_fit(xsub,y,distribution,candidate,max_iter=itmax,tol=tolerance)
                if (.not.candidate%converged) cycle
                candidate_ic = model_ic(candidate,icname)
                if (candidate_ic < best_ic - 1.0e-10_dp) then
                    best_ic = candidate_ic
                    best_j = j
                end if
            end do
            if (best_j == 0) exit
            selected(best_j) = .true.
            call make_subset(x,selected,xsub)
            call alm_fit(xsub,y,distribution,current,max_iter=itmax,tol=tolerance)
            current_ic = model_ic(current,icname)
        end do
        model = current
    end subroutine stepwise_fit

    subroutine calm_fit(x,y,distribution,result,criterion,max_predictors,max_iter,tol)
        real(dp), intent(in) :: x(:,:), y(:)
        character(len=*), intent(in) :: distribution
        type(calm_model), intent(out) :: result
        character(len=*), intent(in), optional :: criterion
        integer, intent(in), optional :: max_predictors, max_iter
        real(dp), intent(in), optional :: tol
        type(alm_model) :: fit
        logical, allocatable :: mask(:)
        real(dp), allocatable :: xsub(:,:), all_ic(:), all_beta(:,:), all_scale_beta(:,:)
        real(dp) :: min_ic, sw, tolerance
        integer :: p, q, m, nmodels, maxp, itmax, j, k
        character(len=8) :: icname

        p = size(x,2)
        q = max(0,p-1)
        maxp = 14
        if (present(max_predictors)) maxp = max_predictors
        if (q > maxp) error stop 'calm_fit: too many predictors for exhaustive averaging'
        nmodels = 2**q
        icname = 'AICc'
        if (present(criterion)) icname = adjustl(criterion)
        itmax = 600
        if (present(max_iter)) itmax = max_iter
        tolerance = 1.0e-6_dp
        if (present(tol)) tolerance = tol

        allocate(mask(p),all_ic(nmodels),all_beta(p,nmodels),all_scale_beta(p,nmodels))
        all_ic = huge(1.0_dp)
        all_beta = 0.0_dp
        all_scale_beta = 0.0_dp
        do m=0,nmodels-1
            mask = .false.
            mask(1) = .true.
            do j=1,q
                mask(j+1) = btest(m,j-1)
            end do
            call make_subset(x,mask,xsub)
            call alm_fit(xsub,y,distribution,fit,max_iter=itmax,tol=tolerance)
            if (.not.fit%converged) cycle
            all_ic(m+1) = model_ic(fit,icname)
            k = 0
            do j=1,p
                if (.not.mask(j)) cycle
                k = k + 1
                all_beta(j,m+1) = fit%beta(k)
                if (trim(distribution) == 'dbeta' .and. allocated(fit%scale_beta)) then
                    all_scale_beta(j,m+1) = fit%scale_beta(k)
                end if
            end do
        end do
        min_ic = minval(all_ic)
        allocate(result%beta(p),result%inclusion_probability(p),result%model_weight(nmodels))
        if (trim(distribution) == 'dbeta') allocate(result%scale_beta(p))
        result%model_weight = 0.0_dp
        where(all_ic < huge(1.0_dp)/10.0_dp)
            result%model_weight = exp(-0.5_dp*(all_ic-min_ic))
        end where
        sw = sum(result%model_weight)
        if (sw <= tiny(1.0_dp)) then
            result%converged = .false.
            result%beta = 0.0_dp
            if (allocated(result%scale_beta)) result%scale_beta = 0.0_dp
            result%inclusion_probability = 0.0_dp
            return
        end if
        result%model_weight = result%model_weight/sw
        result%beta = matmul(all_beta,result%model_weight)
        if (allocated(result%scale_beta)) result%scale_beta = matmul(all_scale_beta,result%model_weight)
        result%inclusion_probability = 0.0_dp
        result%inclusion_probability(1) = 1.0_dp
        do m=0,nmodels-1
            do j=1,q
                if (btest(m,j-1)) then
                    result%inclusion_probability(j+1) = result%inclusion_probability(j+1) + &
                        result%model_weight(m+1)
                end if
            end do
        end do
        result%distribution = distribution
        result%criterion = icname
        result%best_ic = min_ic
        result%nmodels = nmodels
        result%converged = .true.
    end subroutine calm_fit

    subroutine calm_predict(model,xnew,pred)
        type(calm_model), intent(in) :: model
        real(dp), intent(in) :: xnew(:,:)
        real(dp), intent(out) :: pred(:)
        real(dp) :: eta(size(pred))
        eta = matmul(xnew,model%beta)
        select case(trim(model%distribution))
        case('dexp','dgamma','dinvgauss','dpois','dgeom','dnbinom','dbinom')
            pred = exp(min(eta,700.0_dp))
        case('dlnorm','dllaplace','dls','dlgnorm')
            pred = exp(eta)
        case('dbeta')
            if (.not.allocated(model%scale_beta)) then
                pred = 0.5_dp
            else
                pred = exp(min(eta,700.0_dp)) / (exp(min(eta,700.0_dp)) + &
                    exp(min(matmul(xnew,model%scale_beta),700.0_dp)))
            end if
        case('plogis','dlogitnorm')
            where(eta >= 0.0_dp)
                pred = 1.0_dp/(1.0_dp+exp(-eta))
            elsewhere
                pred = exp(eta)/(1.0_dp+exp(eta))
            end where
        case default
            pred = eta
        end select
    end subroutine calm_predict

    subroutine lm_dynamic_fit(x,y,distribution,result,criterion,smooth_weights,span,max_predictors,max_iter,tol)
        real(dp), intent(in) :: x(:,:), y(:)
        character(len=*), intent(in) :: distribution
        type(calm_dynamic_model), intent(out) :: result
        character(len=*), intent(in), optional :: criterion
        logical, intent(in), optional :: smooth_weights
        real(dp), intent(in), optional :: span
        integer, intent(in), optional :: max_predictors, max_iter
        real(dp), intent(in), optional :: tol
        type(alm_model) :: fit
        logical, allocatable :: mask(:)
        real(dp), allocatable :: xsub(:,:), all_beta(:,:), all_pred(:,:), pic(:,:), w(:,:), z(:,:)
        real(dp), allocatable :: trial_w(:,:), sm(:)
        real(dp) :: row_min, sw, fspan, best_span, score, best_score, tolerance
        integer :: n,p,q,nmodels,m,j,k,maxp,itmax,i,ig
        character(len=8) :: icname
        logical :: do_smooth

        n=size(y);p=size(x,2);q=max(0,p-1)
        maxp=12;if(present(max_predictors))maxp=max_predictors
        if(q>maxp)error stop 'lm_dynamic_fit: too many predictors for exhaustive averaging'
        nmodels=2**q
        icname='AICc';if(present(criterion))icname=adjustl(criterion)
        do_smooth=.true.;if(present(smooth_weights))do_smooth=smooth_weights
        itmax=500;if(present(max_iter))itmax=max_iter
        tolerance=1.0e-6_dp;if(present(tol))tolerance=tol
        allocate(mask(p),all_beta(p,nmodels),all_pred(n,nmodels),pic(n,nmodels),w(n,nmodels))
        all_beta=0.0_dp;all_pred=0.0_dp;pic=huge(1.0_dp)
        do m=0,nmodels-1
            mask=.false.;mask(1)=.true.
            do j=1,q
                mask(j+1)=btest(m,j-1)
            end do
            call make_subset(x,mask,xsub)
            call alm_fit(xsub,y,distribution,fit,max_iter=itmax,tol=tolerance)
            if(.not.fit%converged)cycle
            k=0
            do j=1,p
                if(.not.mask(j))cycle
                k=k+1;all_beta(j,m+1)=fit%beta(k)
            end do
            call alm_predict(fit,xsub,all_pred(:,m+1))
            select case(trim(icname))
            case('AIC','aic');pic(:,m+1)=point_aic(fit)
            case('BIC','bic');pic(:,m+1)=point_bic(fit)
            case('BICc','bicc');pic(:,m+1)=point_bicc(fit)
            case default;pic(:,m+1)=point_aicc(fit)
            end select
        end do
        do i=1,n
            row_min=minval(pic(i,:))
            w(i,:)=0.0_dp
            where(pic(i,:)<huge(1.0_dp)/10.0_dp)
                w(i,:)=exp(-0.5_dp*(pic(i,:)-row_min))
            end where
            sw=sum(w(i,:));if(sw>tiny(1.0_dp))w(i,:)=w(i,:)/sw
        end do
        fspan=0.67_dp
        if(present(span))fspan=max(0.05_dp,min(1.0_dp,span))
        if(do_smooth)then
            allocate(z(n,nmodels),trial_w(n,nmodels),sm(n))
            z=log(max(w,1.0e-12_dp)/max(1.0_dp-w,1.0e-12_dp))
            if(.not.present(span))then
                best_score=huge(1.0_dp);best_span=fspan
                do ig=2,10
                    fspan=0.1_dp*real(ig,dp)
                    call smooth_weight_matrix(z,fspan,trial_w)
                    score=sum(trial_w*pic)/real(n,dp)
                    if(score<best_score)then;best_score=score;best_span=fspan;end if
                end do
                fspan=best_span
            end if
            call smooth_weight_matrix(z,fspan,w)
        end if
        allocate(result%beta(p),result%coefficients_dynamic(n,p),result%importance(n,p), &
            result%weights(n,nmodels),result%point_ic(n,nmodels),result%fitted(n),result%residuals(n))
        result%weights=w;result%point_ic=pic;result%coefficients_dynamic=matmul(w,transpose(all_beta))
        result%beta=sum(result%coefficients_dynamic,dim=1)/real(n,dp)
        result%importance=0.0_dp;result%importance(:,1)=1.0_dp
        do m=0,nmodels-1
            do j=1,q
                if(btest(m,j-1))result%importance(:,j+1)=result%importance(:,j+1)+w(:,m+1)
            end do
        end do
        do i=1,n
            result%fitted(i)=sum(w(i,:)*all_pred(i,:))
        end do
        result%residuals=y-result%fitted
        result%distribution=distribution;result%criterion=icname;result%span=fspan
        result%nmodels=nmodels;result%converged=.true.
    end subroutine lm_dynamic_fit

    subroutine smooth_weight_matrix(logit_w,span,wout)
        real(dp),intent(in)::logit_w(:,:),span
        real(dp),intent(out)::wout(:,:)
        real(dp),allocatable::sm(:)
        real(dp)::sw
        integer::j,i
        allocate(sm(size(logit_w,1)))
        do j=1,size(logit_w,2)
            call lowess_numeric(logit_w(:,j),span,sm)
            where(sm>=0.0_dp)
                wout(:,j)=1.0_dp/(1.0_dp+exp(-min(sm,700.0_dp)))
            elsewhere
                wout(:,j)=exp(max(sm,-700.0_dp))/(1.0_dp+exp(max(sm,-700.0_dp)))
            end where
        end do
        do i=1,size(wout,1)
            sw=sum(wout(i,:));if(sw>tiny(1.0_dp))wout(i,:)=wout(i,:)/sw
        end do
    end subroutine smooth_weight_matrix

    subroutine lowess_numeric(y,span,smooth)
        real(dp),intent(in)::y(:),span
        real(dp),intent(out)::smooth(:)
        real(dp),allocatable::robust(:),weights(:),resid(:),absres(:)
        real(dp)::radius,d,u,s0,s1,s2,t0,t1,den,med,r
        integer::n,i,j,m,it
        n=size(y);m=max(2,min(n,ceiling(max(0.05_dp,min(1.0_dp,span))*real(n,dp))))
        allocate(robust(n),weights(n),resid(n),absres(n));robust=1.0_dp
        do it=1,3
            do i=1,n
                radius=real(max(i-1,n-i),dp)
                call kth_distance(i,n,m,radius)
                radius=max(radius,1.0_dp)
                s0=0.0_dp;s1=0.0_dp;s2=0.0_dp;t0=0.0_dp;t1=0.0_dp
                do j=1,n
                    d=abs(real(j-i,dp))/radius
                    if(d>=1.0_dp)then
                        weights(j)=0.0_dp
                    else
                        u=(1.0_dp-d**3)**3;weights(j)=u*robust(j)
                    end if
                    s0=s0+weights(j);s1=s1+weights(j)*real(j-i,dp)
                    s2=s2+weights(j)*real(j-i,dp)**2
                    t0=t0+weights(j)*y(j);t1=t1+weights(j)*real(j-i,dp)*y(j)
                end do
                den=s0*s2-s1*s1
                if(abs(den)>tiny(1.0_dp))then
                    smooth(i)=(t0*s2-t1*s1)/den
                else if(s0>tiny(1.0_dp))then
                    smooth(i)=t0/s0
                else
                    smooth(i)=y(i)
                end if
            end do
            if(it==3)exit
            resid=y-smooth;absres=abs(resid);call sort_local(absres)
            if(mod(n,2)==0)then;med=0.5_dp*(absres(n/2)+absres(n/2+1));else;med=absres((n+1)/2);end if
            if(med<=tiny(1.0_dp))exit
            do j=1,n
                r=abs(resid(j))/(6.0_dp*med)
                if(r>=1.0_dp)then;robust(j)=0.0_dp;else;robust(j)=(1.0_dp-r*r)**2;end if
            end do
        end do
    end subroutine lowess_numeric

    subroutine kth_distance(i,n,m,radius)
        integer,intent(in)::i,n,m
        real(dp),intent(out)::radius
        real(dp)::d(n),t
        integer::j,k
        do j=1,n;d(j)=abs(real(j-i,dp));end do
        do j=2,n;t=d(j);k=j-1;do while(k>=1);if(d(k)<=t)exit;d(k+1)=d(k);k=k-1;end do;d(k+1)=t;end do
        radius=d(m)
    end subroutine kth_distance

    subroutine sort_local(a)
        real(dp),intent(inout)::a(:)
        real(dp)::t;integer::i,j
        do i=2,size(a);t=a(i);j=i-1;do while(j>=1);if(a(j)<=t)exit;a(j+1)=a(j);j=j-1;end do;a(j+1)=t;end do
    end subroutine sort_local

    subroutine rolling_origin_alm(x,y,distribution,initial,forecast,criterion)
        real(dp), intent(in) :: x(:,:), y(:)
        character(len=*), intent(in) :: distribution
        integer, intent(in) :: initial
        real(dp), intent(out) :: forecast(:)
        character(len=*), intent(in), optional :: criterion
        type(alm_model) :: fit
        real(dp) :: one(1)
        integer :: n, i
        character(len=8) :: dummy

        n = size(y)
        if (initial < 2 .or. initial >= n) error stop 'rolling_origin_alm: invalid initial'
        if (size(forecast) /= n-initial) error stop 'rolling_origin_alm: forecast has wrong size'
        dummy = 'AICc'
        if (present(criterion)) dummy = criterion
        do i=initial,n-1
            call alm_fit(x(1:i,:),y(1:i),distribution,fit,max_iter=500,tol=1.0e-6_dp)
            call alm_predict(fit,x(i+1:i+1,:),one)
            forecast(i-initial+1) = one(1)
        end do
        if (len_trim(dummy) < 0) forecast = forecast
    end subroutine rolling_origin_alm

    subroutine recursive_lm(x,y,forgetting,beta_time,fitted)
        real(dp), intent(in) :: x(:,:), y(:), forgetting
        real(dp), intent(out) :: beta_time(:,:), fitted(:)
        real(dp), allocatable :: pmat(:,:), beta(:), xt(:), gain(:), row(:)
        real(dp) :: denom, err
        integer :: n, p, i, j

        n = size(y)
        p = size(x,2)
        if (size(beta_time,1) /= p .or. size(beta_time,2) /= n) then
            error stop 'recursive_lm: beta_time has wrong shape'
        end if
        if (size(fitted) /= n) error stop 'recursive_lm: fitted has wrong size'
        if (forgetting <= 0.0_dp .or. forgetting > 1.0_dp) then
            error stop 'recursive_lm: forgetting must be in (0,1]'
        end if
        allocate(pmat(p,p),beta(p),xt(p),gain(p),row(p))
        pmat = 0.0_dp
        do j=1,p
            pmat(j,j) = 1.0e6_dp
        end do
        beta = 0.0_dp
        do i=1,n
            xt = x(i,:)
            denom = forgetting + dot_product(xt,matmul(pmat,xt))
            gain = matmul(pmat,xt)/max(denom,tiny(1.0_dp))
            fitted(i) = dot_product(xt,beta)
            err = y(i)-fitted(i)
            beta = beta + gain*err
            row = matmul(xt,pmat)
            pmat = (pmat - spread(gain,2,p)*spread(row,1,p))/forgetting
            beta_time(:,i) = beta
            fitted(i) = dot_product(xt,beta)
        end do
    end subroutine recursive_lm

    subroutine make_subset(x,mask,xsub)
        real(dp), intent(in) :: x(:,:)
        logical, intent(in) :: mask(:)
        real(dp), allocatable, intent(out) :: xsub(:,:)
        integer :: j, k
        allocate(xsub(size(x,1),count(mask)))
        k = 0
        do j=1,size(mask)
            if (.not.mask(j)) cycle
            k = k + 1
            xsub(:,k) = x(:,j)
        end do
    end subroutine make_subset

    pure real(dp) function model_ic(model,criterion) result(ic)
        type(alm_model), intent(in) :: model
        character(len=*), intent(in) :: criterion
        select case(trim(criterion))
        case('AIC','aic')
            ic = model%aic
        case('BIC','bic')
            ic = model%bic
        case('BICc','bicc')
            ic = model%bicc
        case default
            ic = model%aicc
        end select
    end function model_ic

end module greybox_selection
