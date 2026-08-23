module greybox_regression
    use greybox_kinds, only: dp, pi
    use greybox_special, only: normal_pdf, normal_cdf, student_t_pdf, logistic_pdf, &
        nan_dp, log1pexp
    use greybox_distributions, only: dlaplace, ds, dgnorm, dalaplace, dbcnorm, dfnorm, drectnorm, dlogitnorm
    use greybox_linalg, only: least_squares, invert_matrix
    use greybox_optimizer, only: pattern_search, numerical_hessian
    implicit none
    private
    public :: alm_model, alm_fit, alm_fit_beta, alm_predict, alm_information, coef_bootstrap, scale_model_fit
    public :: point_aic, point_aicc, point_bic, point_bicc

    type :: alm_model
        character(len=16) :: distribution='dnorm'
        character(len=16) :: loss='likelihood'
        real(dp),allocatable :: beta(:)
        real(dp),allocatable :: scale_beta(:)
        real(dp),allocatable :: fitted(:)
        real(dp),allocatable :: residuals(:)
        real(dp),allocatable :: point_loglik(:)
        real(dp),allocatable :: vcov(:,:)
        real(dp),allocatable :: vcov_scale(:,:), vcov_cross(:,:)
        real(dp) :: scale=1.0_dp
        real(dp) :: other=0.0_dp
        real(dp) :: loglik=0.0_dp
        real(dp) :: aic=0.0_dp
        real(dp) :: aicc=0.0_dp
        real(dp) :: bic=0.0_dp
        real(dp) :: bicc=0.0_dp
        integer :: nobs=0
        integer :: nparam=0
        logical :: converged=.false.
    end type alm_model

contains

    subroutine alm_fit(x,y,distribution,model,loss,binom_size,max_iter,tol,lambda,trim_fraction)
        real(dp),intent(in)::x(:,:),y(:)
        character(len=*),intent(in)::distribution
        type(alm_model),intent(out)::model
        character(len=*),intent(in),optional::loss
        integer,intent(in),optional::binom_size,max_iter
        real(dp),intent(in),optional::tol,lambda,trim_fraction
        real(dp),allocatable::par(:),b0(:),worky(:),hess(:,:),hinv(:,:)
        real(dp)::fval,tolerance,penalty,trim_value,quantile_value
        integer::p,n,nextra,info,itmax,bs
        character(len=16)::loss_name
        n=size(y);p=size(x,2);bs=1;if(present(binom_size))bs=binom_size
        if (trim(distribution) == 'dbeta') then
            call alm_fit_beta(x,y,model,max_iter=max_iter,tol=tol)
            return
        end if
        loss_name='likelihood';if(present(loss))loss_name=lower_ascii(trim(adjustl(loss)))
        tolerance=1.0e-7_dp;if(present(tol))tolerance=tol
        itmax=1200;if(present(max_iter))itmax=max_iter
        penalty=0.0_dp;if(present(lambda))penalty=lambda
        trim_value=0.05_dp;if(present(trim_fraction))trim_value=max(0.0_dp,min(0.49_dp,trim_fraction))
        quantile_value=0.5_dp;if(present(lambda))quantile_value=max(0.0_dp,min(1.0_dp,lambda))
        nextra=extra_count(distribution)
        allocate(par(p+nextra),b0(p),worky(n));par=0.0_dp
        worky=y
        select case(trim(distribution))
        case('dlnorm','dllaplace','dls','dlgnorm')
            if(any(y<=0.0_dp))then;call init_invalid(model,distribution,n,p);return;end if
            worky=log(y)
        case('dlogitnorm')
            if(any(y<=0.0_dp).or.any(y>=1.0_dp))then;call init_invalid(model,distribution,n,p);return;end if
            worky=log(y/(1.0_dp-y))
        case('dexp','dgamma','dinvgauss','dpois','dgeom','dnbinom','dbinom')
            worky=log(max(y,1.0e-3_dp))
        case('plogis','pnorm')
            worky=0.0_dp
        end select
        call least_squares(x,worky,b0,info,ridge=1.0e-10_dp)
        if(info/=0)b0=0.0_dp
        par(1:p)=b0
        if(nextra==1)then
            select case(trim(distribution))
            case('dgnorm','dlgnorm');par(p+1)=log(2.0_dp)
            case('dalaplace');par(p+1)=0.0_dp
            case('dfnorm','drectnorm');par(p+1)=log(max(stddev(y),1.0e-3_dp))
            case('dt');par(p+1)=log(6.0_dp-2.0_dp)
            case('dnbinom');par(p+1)=log(max(sum(y)/real(n,dp),1.0_dp))
            case('dbcnorm');par(p+1)=0.0_dp
            end select
        end if
        if(is_direct_loss(loss_name,distribution))then
            select case(trim(loss_name))
            case('mse')
                call least_squares(x,worky,par(1:p),info,ridge=max(0.0_dp,penalty))
                fval=objective(par)
            case default
                call pattern_search(objective,par,fval,max_iter=itmax,tol=tolerance)
            end select
        else
            call pattern_search(objective,par,fval,max_iter=itmax,tol=tolerance)
        end if
        allocate(model%beta(p),model%fitted(n),model%residuals(n),model%point_loglik(n),model%vcov(p,p))
        model%distribution=trim(distribution);model%loss=loss_name;model%beta=par(1:p)
        call unpack_model(par,model%fitted,model%scale,model%other)
        model%residuals=y-model%fitted
        call fill_point_loglik(par,model%point_loglik)
        select case(trim(loss_name))
        case('role')
            model%loglik=real(n,dp)*trimmed_mean(model%point_loglik,trim_value)
        case('quale')
            model%loglik=real(n,dp)*sample_quantile(model%point_loglik,quantile_value)
        case('lasso','ridge')
            model%loglik=nan_dp()
        case default
            model%loglik=sum(model%point_loglik)
        end select
        model%nobs=n;model%nparam=p+nextra+scale_parameter_count(distribution)
        call set_ic(model)
        model%converged=.true.
        ! Numerical covariance for coefficients from full objective Hessian.
        allocate(hess(size(par),size(par)),hinv(size(par),size(par)))
        call numerical_hessian(objective,par,hess);call invert_matrix(hess,hinv,info)
        if(info==0)then;model%vcov=hinv(1:p,1:p);else;model%vcov=0.0_dp;end if

    contains
        function objective(theta) result(v)
            real(dp),intent(in)::theta(:)
            real(dp)::v,mu(n),sc,oth,e
            integer::i
            select case(trim(loss_name))
            case('mse')
                v=sum((target_values()-linear_values(theta))**2)/real(n,dp)
            case('mae')
                v=sum(abs(target_values()-linear_values(theta)))/real(n,dp)
            case('ham')
                v=sum(sqrt(abs(target_values()-linear_values(theta))))/real(n,dp)
            case('lasso','ridge')
                call unpack_model(theta,mu,sc,oth)
                v=(1.0_dp-penalty)*sum((y-mu)**2)/real(n,dp)
                if (trim(loss_name) == 'lasso') then
                    v=v+penalty*sum(abs(theta(1:p)))
                else if (p > 1) then
                    v=v+penalty*sqrt(sum(theta(2:p)**2))
                end if
            case('role','quale')
                call unpack_latent(theta,mu,sc,oth)
                block
                    real(dp) :: vals(n)
                    do i=1,n
                        vals(i)=point_loglik(y(i),mu(i),sc,oth,distribution,bs)
                        if(.not.(vals(i)>-huge(1.0_dp).and.vals(i)<huge(1.0_dp)))then
                            v=huge(1.0_dp)/100.0_dp;return
                        end if
                    end do
                    if (trim(loss_name) == 'role') then
                        v=-real(n,dp)*trimmed_mean(vals,trim_value)
                    else
                        v=-real(n,dp)*sample_quantile(vals,quantile_value)
                    end if
                end block
            case default
                call unpack_latent(theta,mu,sc,oth)
                v=0.0_dp
                do i=1,n
                    e=point_loglik(y(i),mu(i),sc,oth,distribution,bs)
                    if(.not.(e>-huge(1.0_dp).and.e<huge(1.0_dp)))then;v=huge(1.0_dp)/100.0_dp;return;end if
                    v=v-e
                end do
            end select
        end function objective
        function target_values() result(z)
            real(dp)::z(n);z=worky
        end function target_values
        function linear_values(theta) result(z)
            real(dp),intent(in)::theta(:);real(dp)::z(n);z=matmul(x,theta(1:p))
        end function linear_values
        subroutine unpack_model(theta,mu,sc,oth)
            real(dp),intent(in)::theta(:);real(dp),intent(out)::mu(:),sc,oth
            real(dp)::eta(n),res(n),trans(n),eps
            eta=matmul(x,theta(1:p));oth=0.0_dp;eps=1.0e-12_dp
            select case(trim(distribution))
            case('dexp','dgamma','dinvgauss','dpois','dgeom','dnbinom','dbinom')
                mu=exp(min(eta,700.0_dp))
            case('dlnorm','dllaplace','dls','dlgnorm','dlogitnorm','dbcnorm')
                mu=eta
            case('plogis')
                where(eta>=0.0_dp);mu=1.0_dp/(1.0_dp+exp(-eta));elsewhere;mu=exp(eta)/(1.0_dp+exp(eta));end where
            case('pnorm')
                mu=normal_cdf(eta,0.0_dp,1.0_dp)
            case default
                mu=eta
            end select
            select case(trim(distribution))
            case('dgnorm','dlgnorm');oth=exp(theta(p+1))
            case('dalaplace');oth=1.0_dp/(1.0_dp+exp(-theta(p+1)))
            case('dfnorm','drectnorm');oth=exp(theta(p+1))
            case('dt');oth=2.0_dp+exp(theta(p+1))
            case('dnbinom');oth=exp(theta(p+1))
            case('dbcnorm');oth=max(-5.0_dp,min(5.0_dp,theta(p+1)))
            end select
            select case(trim(distribution))
            case('dnorm');res=y-mu;sc=sqrt(sum(res*res)/real(n,dp))
            case('dlaplace');res=y-mu;sc=sum(abs(res))/real(n,dp)
            case('ds');res=y-mu;sc=sum(sqrt(abs(res)))/(2.0_dp*real(n,dp))
            case('dgnorm');res=y-mu;sc=(oth*sum(abs(res)**oth)/real(n,dp))**(1.0_dp/oth)
            case('dlogis');res=y-mu;sc=sqrt(sum(res*res)/real(n,dp)*3.0_dp/pi**2)
            case('dalaplace');res=y-mu;sc=sum(res*(oth-merge(1.0_dp,0.0_dp,y<=mu)))/real(n,dp);sc=max(abs(sc),eps)
            case('dlnorm','dllaplace','dls','dlgnorm')
                trans=log(max(y,eps));res=trans-mu
                select case(trim(distribution))
                case('dlnorm');sc=sqrt(sum(res*res)/real(n,dp))
                case('dllaplace');sc=sum(abs(res))/real(n,dp)
                case('dls');sc=sum(sqrt(abs(res)))/(2.0_dp*real(n,dp))
                case default;sc=(oth*sum(abs(res)**oth)/real(n,dp))**(1.0_dp/oth)
                end select
            case('dbcnorm')
                trans=boxcox_vec(y,oth);res=trans-mu;sc=sqrt(sum(res*res)/real(n,dp))
            case('dinvgauss');sc=sum((y/mu-1.0_dp)**2/(y/mu))/real(n,dp)
            case('dgamma');sc=sum((y/mu-1.0_dp)**2)/real(n,dp)
            case('dlogitnorm');trans=log(y/(1.0_dp-y));res=trans-mu;sc=sqrt(sum(res*res)/real(n,dp))
            case('dfnorm','drectnorm','dt','dnbinom');sc=oth
            case default;sc=1.0_dp
            end select
            sc=max(sc,eps)
            ! Store means/fitted values on response scale rather than latent eta.
            select case(trim(distribution))
            case('dlnorm','dllaplace','dls','dlgnorm');mu=exp(mu)
            case('dlogitnorm');where(mu>=0.0_dp);mu=1.0_dp/(1.0_dp+exp(-mu));elsewhere;mu=exp(mu)/(1.0_dp+exp(mu));end where
            case('dbcnorm');mu=inv_boxcox_vec(mu,oth)
            case('dfnorm');mu=sqrt(2.0_dp/pi)*sc*exp(-mu*mu/(2.0_dp*sc*sc))+mu*(1.0_dp-2.0_dp*normal_cdf(-mu/sc,0.0_dp,1.0_dp))
            case('drectnorm');mu=mu*(1.0_dp-normal_cdf(0.0_dp,mu,sc))+sc*normal_pdf(0.0_dp,mu,sc)
            end select
        end subroutine unpack_model
        subroutine fill_point_loglik(theta,values)
            real(dp),intent(in)::theta(:)
            real(dp),intent(out)::values(:)
            real(dp)::sc,oth,latent(n)
            integer::i
            call unpack_latent(theta,latent,sc,oth)
            do i=1,n
                values(i)=point_loglik(y(i),latent(i),sc,oth,distribution,bs)
            end do
        end subroutine fill_point_loglik
        subroutine unpack_latent(theta,mu,sc,oth)
            real(dp),intent(in)::theta(:);real(dp),intent(out)::mu(:),sc,oth
            real(dp)::fitted_tmp(n)
            ! Reproduce parameter/scale calculation but return latent mu.
            call unpack_model(theta,fitted_tmp,sc,oth)
            mu=matmul(x,theta(1:p))
            select case(trim(distribution))
            case('dexp','dgamma','dinvgauss','dpois','dgeom','dnbinom','dbinom');mu=exp(min(mu,700.0_dp))
            case('plogis');where(mu>=0.0_dp);mu=1.0_dp/(1.0_dp+exp(-mu));elsewhere;mu=exp(mu)/(1.0_dp+exp(mu));end where
            case('pnorm');mu=normal_cdf(mu,0.0_dp,1.0_dp)
            end select
        end subroutine unpack_latent
    end subroutine alm_fit

    subroutine alm_fit_beta(x,y,model,max_iter,tol)
        real(dp), intent(in) :: x(:,:), y(:)
        type(alm_model), intent(out) :: model
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: par(:), start(:), z(:), hess(:,:), hinv(:,:)
        real(dp) :: fval, tolerance, yy, a, b
        integer :: n, p, i, info, itmax

        n = size(y); p = size(x,2)
        if (size(x,1) /= n .or. n < 2 .or. any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            call init_invalid(model,'dbeta',n,p)
            return
        end if
        tolerance = 1.0e-7_dp
        if (present(tol)) tolerance = tol
        itmax = 1600
        if (present(max_iter)) itmax = max_iter
        allocate(par(2*p),start(p),z(n))
        z = log(max(y,1.0e-10_dp)/max(1.0_dp-y,1.0e-10_dp))
        call least_squares(x,z,start,info,ridge=1.0e-8_dp)
        if (info /= 0) start = 0.0_dp
        par(1:p) = start
        par(p+1:2*p) = -start
        call pattern_search(objective_beta,par,fval,max_iter=itmax,tol=tolerance)

        allocate(model%beta(p),model%scale_beta(p),model%fitted(n),model%residuals(n), &
            model%point_loglik(n),model%vcov(p,p),model%vcov_scale(p,p),model%vcov_cross(p,p))
        model%distribution = 'dbeta'
        model%loss = 'likelihood'
        model%beta = par(1:p)
        model%scale_beta = par(p+1:2*p)
        do i=1,n
            a = exp(min(dot_product(x(i,:),model%beta),700.0_dp))
            b = exp(min(dot_product(x(i,:),model%scale_beta),700.0_dp))
            model%fitted(i) = a/(a+b)
            yy = min(1.0_dp-1.0e-10_dp,max(1.0e-10_dp,y(i)))
            model%point_loglik(i) = (a-1.0_dp)*log(yy) + (b-1.0_dp)*log(1.0_dp-yy) + &
                log_gamma(a+b)-log_gamma(a)-log_gamma(b)
        end do
        model%residuals = y-model%fitted
        model%scale = sum(exp(min(matmul(x,model%beta),700.0_dp))* &
            exp(min(matmul(x,model%scale_beta),700.0_dp)) / &
            ((exp(min(matmul(x,model%beta),700.0_dp))+exp(min(matmul(x,model%scale_beta),700.0_dp)))**2 * &
             (exp(min(matmul(x,model%beta),700.0_dp))+exp(min(matmul(x,model%scale_beta),700.0_dp))+1.0_dp))) / real(n,dp)
        model%loglik = sum(model%point_loglik)
        model%nobs = n
        model%nparam = 2*p
        call set_ic(model)
        model%converged = .true.
        allocate(hess(2*p,2*p),hinv(2*p,2*p))
        call numerical_hessian(objective_beta,par,hess)
        call invert_matrix(hess,hinv,info)
        if (info == 0) then
            model%vcov = hinv(1:p,1:p)
            model%vcov_scale = hinv(p+1:2*p,p+1:2*p)
            model%vcov_cross = hinv(1:p,p+1:2*p)
        else
            model%vcov = 0.0_dp
            model%vcov_scale = 0.0_dp
            model%vcov_cross = 0.0_dp
        end if

    contains
        function objective_beta(theta) result(v)
            real(dp), intent(in) :: theta(:)
            real(dp) :: v, aa, bb, ysafe
            integer :: j
            v = 0.0_dp
            do j=1,n
                aa = exp(min(dot_product(x(j,:),theta(1:p)),700.0_dp))
                bb = exp(min(dot_product(x(j,:),theta(p+1:2*p)),700.0_dp))
                ysafe = min(1.0_dp-1.0e-10_dp,max(1.0e-10_dp,y(j)))
                v = v - ((aa-1.0_dp)*log(ysafe) + (bb-1.0_dp)*log(1.0_dp-ysafe) + &
                    log_gamma(aa+bb)-log_gamma(aa)-log_gamma(bb))
                if (.not.(v < huge(1.0_dp)/100.0_dp)) then
                    v = huge(1.0_dp)/100.0_dp
                    return
                end if
            end do
        end function objective_beta
    end subroutine alm_fit_beta

    subroutine alm_predict(model,xnew,pred)
        type(alm_model),intent(in)::model
        real(dp),intent(in)::xnew(:,:)
        real(dp),intent(out)::pred(:)
        real(dp)::eta(size(xnew,1))
        eta=matmul(xnew,model%beta)
        select case(trim(model%distribution))
        case('dexp','dgamma','dinvgauss','dpois','dgeom','dnbinom','dbinom');pred=exp(min(eta,700.0_dp))
        case('dlnorm','dllaplace','dls','dlgnorm');pred=exp(eta)
        case('dbeta')
            if (.not.allocated(model%scale_beta)) then
                pred = nan_dp()
            else
                pred = exp(min(eta,700.0_dp)) / &
                    (exp(min(eta,700.0_dp)) + exp(min(matmul(xnew,model%scale_beta),700.0_dp)))
            end if
        case('dlogitnorm','plogis')
            where(eta >= 0.0_dp)
                pred = 1.0_dp/(1.0_dp+exp(-eta))
            elsewhere
                pred = exp(eta)/(1.0_dp+exp(eta))
            end where
        case('pnorm');pred=normal_cdf(eta,0.0_dp,1.0_dp)
        case('dbcnorm');pred=inv_boxcox_vec(eta,model%other)
        case('dfnorm')
            pred = sqrt(2.0_dp/pi)*model%scale* &
                exp(-eta*eta/(2.0_dp*model%scale**2)) + eta* &
                (1.0_dp-2.0_dp*normal_cdf(-eta/model%scale,0.0_dp,1.0_dp))
        case('drectnorm');pred=eta*(1.0_dp-normal_cdf(0.0_dp,eta,model%scale))+model%scale*normal_pdf(0.0_dp,eta,model%scale)
        case default;pred=eta
        end select
    end subroutine alm_predict

    pure subroutine alm_information(model,aic,aicc,bic,bicc)
        type(alm_model),intent(in)::model
        real(dp),intent(out),optional::aic,aicc,bic,bicc
        if(present(aic))aic=model%aic;if(present(aicc))aicc=model%aicc
        if(present(bic))bic=model%bic;if(present(bicc))bicc=model%bicc
    end subroutine alm_information

    subroutine coef_bootstrap(x,y,distribution,nsim,coefs)
        real(dp),intent(in)::x(:,:),y(:)
        character(len=*),intent(in)::distribution
        integer,intent(in)::nsim
        real(dp),intent(out)::coefs(:,:)
        real(dp),allocatable::xb(:,:),yb(:)
        type(alm_model)::m
        integer::n,p,s,i,idx
        real(dp)::u
        n=size(y);p=size(x,2);allocate(xb(n,p),yb(n))
        do s=1,nsim
            do i=1,n;call random_number(u);idx=min(n,1+int(u*real(n,dp)));xb(i,:)=x(idx,:);yb(i)=y(idx);end do
            call alm_fit(xb,yb,distribution,m,max_iter=400,tol=1.0e-5_dp);coefs(:,s)=m%beta
        end do
    end subroutine coef_bootstrap

    subroutine scale_model_fit(z,residuals,beta,info)
        real(dp),intent(in)::z(:,:),residuals(:)
        real(dp),intent(out)::beta(:)
        integer,intent(out),optional::info
        real(dp)::target(size(residuals))
        target=log(max(residuals*residuals,1.0e-12_dp))
        call least_squares(z,target,beta,info)
    end subroutine scale_model_fit

    pure integer function extra_count(distribution) result(k)
        character(len=*),intent(in)::distribution
        select case(trim(distribution))
        case('dgnorm','dlgnorm','dalaplace','dfnorm','drectnorm','dt','dnbinom','dbcnorm')
            k = 1
        case default
            k = 0
        end select
    end function extra_count
    pure integer function scale_parameter_count(distribution) result(k)
        character(len=*),intent(in)::distribution
        select case(trim(distribution))
        case('dpois','dgeom','dbinom','plogis','pnorm','dexp','dfnorm','drectnorm','dt','dnbinom')
            k = 0
        case default
            k = 1
        end select
    end function scale_parameter_count
    pure logical function is_direct_loss(loss,distribution) result(v)
        character(len=*),intent(in)::loss,distribution
        v = (trim(loss) /= 'likelihood') .or. (len_trim(distribution) < 0)
    end function is_direct_loss
    subroutine init_invalid(model,distribution,n,p)
        type(alm_model),intent(out)::model;character(len=*),intent(in)::distribution;integer,intent(in)::n,p
        allocate(model%beta(p),model%fitted(n),model%residuals(n),model%point_loglik(n),model%vcov(p,p))
        model%beta = 0.0_dp
        model%fitted = 0.0_dp
        model%residuals = 0.0_dp
        model%point_loglik = 0.0_dp
        model%vcov = 0.0_dp
        model%distribution = distribution
        model%converged = .false.
    end subroutine init_invalid
    pure real(dp) function stddev(y) result(s)
        real(dp),intent(in)::y(:);real(dp)::m;m=sum(y)/real(size(y),dp);s=sqrt(sum((y-m)**2)/real(max(1,size(y)-1),dp))
    end function stddev
    pure function boxcox_vec(y,lambda) result(z)
        real(dp),intent(in)::y(:),lambda
        real(dp)::z(size(y))
        if(abs(lambda)<=epsilon(1.0_dp))then
            z=log(y)
        else
            z=(y**lambda-1.0_dp)/lambda
        end if
    end function boxcox_vec
    pure function inv_boxcox_vec(z,lambda) result(y)
        real(dp),intent(in)::z(:),lambda
        real(dp)::y(size(z))
        if(abs(lambda)<=epsilon(1.0_dp))then
            y=exp(z)
        else
            y=max(0.0_dp,lambda*z+1.0_dp)**(1.0_dp/lambda)
        end if
    end function inv_boxcox_vec

    pure real(dp) function point_loglik(y,mu,scale,other,distribution,binom_size) result(ll)
        real(dp),intent(in)::y,mu,scale,other
        character(len=*),intent(in)::distribution
        integer,intent(in)::binom_size
        real(dp)::prob,shape,sc
        integer::iy
        ll=-huge(1.0_dp)
        select case(trim(distribution))
        case('dnorm');ll=log(max(normal_pdf(y,mu,scale),tiny(1.0_dp)))
        case('dlaplace');ll=dlaplace(y,mu,scale,.true.)
        case('ds');ll=ds(y,mu,scale,.true.)
        case('dgnorm');ll=dgnorm(y,mu,scale,other,.true.)
        case('dlogis');ll=log(max(logistic_pdf(y,mu,scale),tiny(1.0_dp)))
        case('dt');ll=log(max(student_t_pdf(y-mu,scale),tiny(1.0_dp)))
        case('dalaplace');ll=dalaplace(y,mu,scale,other,.true.)
        case('dlnorm')
            if(y>0.0_dp)ll=log(max(normal_pdf(log(y),mu,scale),tiny(1.0_dp)))-log(y)
        case('dllaplace')
            if(y>0.0_dp)ll=dlaplace(log(y),mu,scale,.true.)-log(y)
        case('dls')
            if(y>0.0_dp)ll=ds(log(y),mu,scale,.true.)-log(y)
        case('dlgnorm')
            if(y>0.0_dp)ll=dgnorm(log(y),mu,scale,other,.true.)-log(y)
        case('dbcnorm');ll=dbcnorm(y,mu,scale,other,.true.)
        case('dfnorm');ll=dfnorm(y,mu,scale,.true.)
        case('drectnorm');ll=drectnorm(y,mu,scale,.true.)
        case('dinvgauss')
            if(y>0.0_dp.and.mu>0.0_dp.and.scale>0.0_dp)ll=-0.5_dp*(log(2.0_dp*pi*scale*y**3/mu)+(y-mu)**2/(scale*mu*y))
        case('dgamma')
            if(y>0.0_dp.and.mu>0.0_dp.and.scale>0.0_dp)then
                shape=1.0_dp/scale
                sc=scale*mu
                ll=(shape-1.0_dp)*log(y)-y/sc-log_gamma(shape)-shape*log(sc)
            end if
        case('dexp')
            if(y>=0.0_dp.and.mu>0.0_dp)ll=-log(mu)-y/mu
        case('dpois')
            if(y>=0.0_dp.and.mu>0.0_dp)then;iy=nint(y);if(abs(y-real(iy,dp))<1.0e-8_dp)ll=y*log(mu)-mu-log_gamma(y+1.0_dp);end if
        case('dgeom')
            if(y>=0.0_dp.and.mu>=0.0_dp)then;prob=1.0_dp/(mu+1.0_dp);ll=log(prob)+y*log(1.0_dp-prob);end if
        case('dnbinom')
            if(y>=0.0_dp.and.mu>0.0_dp.and.scale>0.0_dp)then
                ll=log_gamma(y+scale)-log_gamma(scale)-log_gamma(y+1.0_dp) + &
                    scale*log(scale/(scale+mu))+y*log(mu/(scale+mu))
            end if
        case('dbinom')
            iy=nint(y);prob=1.0_dp/(mu+1.0_dp)
            if(iy>=0.and.iy<=binom_size)then
                ll=log_gamma(real(binom_size+1,dp))-log_gamma(real(iy+1,dp)) - &
                    log_gamma(real(binom_size-iy+1,dp))+real(iy,dp)*log(prob) + &
                    real(binom_size-iy,dp)*log(1.0_dp-prob)
            end if
        case('dlogitnorm');ll=dlogitnorm(y,mu,scale,.true.)
        case('plogis','pnorm')
            prob=max(1.0e-15_dp,min(1.0_dp-1.0e-15_dp,mu));ll=y*log(prob)+(1.0_dp-y)*log(1.0_dp-prob)
        end select
    end function point_loglik

    pure function point_aic(model) result(ic)
        type(alm_model),intent(in)::model
        real(dp)::ic(size(model%point_loglik))
        ic=2.0_dp*real(model%nparam,dp)-2.0_dp*real(model%nobs,dp)*model%point_loglik
    end function point_aic

    pure function point_aicc(model) result(ic)
        type(alm_model),intent(in)::model
        real(dp)::ic(size(model%point_loglik)),n,k
        n=real(model%nobs,dp);k=real(model%nparam,dp)
        ic=point_aic(model)
        if(n>k+1.0_dp)ic=ic+2.0_dp*k*(k+1.0_dp)/(n-k-1.0_dp)
    end function point_aicc

    pure function point_bic(model) result(ic)
        type(alm_model),intent(in)::model
        real(dp)::ic(size(model%point_loglik)),n,k
        n=real(model%nobs,dp);k=real(model%nparam,dp)
        ic=log(n)*k-2.0_dp*n*model%point_loglik
    end function point_bic

    pure function point_bicc(model) result(ic)
        type(alm_model),intent(in)::model
        real(dp)::ic(size(model%point_loglik)),n,k
        n=real(model%nobs,dp);k=real(model%nparam,dp)
        if(n>k+1.0_dp)then
            ic=(k*log(n)*n)/(n-k-1.0_dp)-2.0_dp*n*model%point_loglik
        else
            ic=huge(1.0_dp)
        end if
    end function point_bicc

    subroutine set_ic(model)
        type(alm_model),intent(inout)::model
        real(dp)::n,k
        n=real(model%nobs,dp);k=real(model%nparam,dp)
        model%aic=2.0_dp*k-2.0_dp*model%loglik;model%bic=log(n)*k-2.0_dp*model%loglik
        if(n>k+1.0_dp)then
            model%aicc=model%aic+2.0_dp*k*(k+1.0_dp)/(n-k-1.0_dp)
            model%bicc=-2.0_dp*model%loglik+k*log(n)*n/(n-k-1.0_dp)
        else;model%aicc=huge(1.0_dp);model%bicc=huge(1.0_dp);end if
    end subroutine set_ic

    pure function lower_ascii(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, k
        out = s
        do i=1,len(s)
            k = iachar(out(i:i))
            if (k >= iachar('A') .and. k <= iachar('Z')) out(i:i)=achar(k+32)
        end do
    end function lower_ascii

    pure real(dp) function sample_quantile(x,p) result(q)
        real(dp), intent(in) :: x(:), p
        real(dp) :: a(size(x)), t, h
        integer :: i, j, lo
        a=x
        do i=2,size(a)
            t=a(i);j=i-1
            do while(j>=1)
                if(a(j)<=t) exit
                a(j+1)=a(j);j=j-1
            end do
            a(j+1)=t
        end do
        if (size(a)==1) then
            q=a(1); return
        end if
        h=1.0_dp+(real(size(a)-1,dp))*min(1.0_dp,max(0.0_dp,p))
        lo=int(floor(h))
        if(lo>=size(a))then
            q=a(size(a))
        else
            q=a(lo)+(h-real(lo,dp))*(a(lo+1)-a(lo))
        end if
    end function sample_quantile

    pure real(dp) function trimmed_mean(x,trim_fraction) result(m)
        real(dp), intent(in) :: x(:), trim_fraction
        real(dp) :: a(size(x)), t
        integer :: i,j,klo,khi
        a=x
        do i=2,size(a)
            t=a(i);j=i-1
            do while(j>=1)
                if(a(j)<=t) exit
                a(j+1)=a(j);j=j-1
            end do
            a(j+1)=t
        end do
        klo=int(floor(trim_fraction*real(size(a),dp)))+1
        khi=size(a)-int(floor(trim_fraction*real(size(a),dp)))
        if(khi<klo)then
            m=sum(a)/real(size(a),dp)
        else
            m=sum(a(klo:khi))/real(khi-klo+1,dp)
        end if
    end function trimmed_mean

end module greybox_regression
