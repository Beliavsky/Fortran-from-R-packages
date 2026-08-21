module circstats_models
    use circstats_kinds, only: dp, twopi
    use circstats_types, only: vm_fit_result, wrapped_cauchy_fit_result, circ_reg_result, vm_bootstrap_result
    use circstats_utils, only: wrap_2pi, quantile_type7, ols_fit, solve_linear, randu
    use circstats_core, only: circ_mean, est_kappa, a1inv
    implicit none
    private
    public :: vm_ml, wrpcauchy_ml, circ_reg, vm_bootstrap_ci, pp_fit

contains

    pure function vm_ml(x,bias) result(fit)
        real(dp), intent(in) :: x(:)
        logical, intent(in), optional :: bias
        type(vm_fit_result) :: fit
        fit%mu=circ_mean(x)
        fit%kappa=est_kappa(x,bias)
    end function vm_ml

    pure function pp_fit(x) result(fit)
        real(dp), intent(in) :: x(:)
        type(vm_fit_result) :: fit
        fit%mu=wrap_2pi(circ_mean(x))
        fit%kappa=est_kappa(x)
    end function pp_fit

    function wrpcauchy_ml(x,mu0,rho0,acc,max_iter) result(fit)
        real(dp), intent(in) :: x(:),mu0,rho0
        real(dp), intent(in), optional :: acc
        integer, intent(in), optional :: max_iter
        type(wrapped_cauchy_fit_result) :: fit
        real(dp), allocatable :: w(:)
        real(dp) :: m1old,m2old,m1new,m2new,tol,mconst
        integer :: iter,limit
        if (rho0 < 0.0_dp .or. rho0 >= 1.0_dp) error stop "wrpcauchy_ml: rho0 must be in [0,1)"
        tol=1.0e-12_dp
        if (present(acc)) tol=acc
        limit=100000
        if (present(max_iter)) limit=max_iter
        allocate(w(size(x)))
        m1old=2.0_dp*rho0*cos(mu0)/(1.0_dp+rho0*rho0)
        m2old=2.0_dp*rho0*sin(mu0)/(1.0_dp+rho0*rho0)
        m1new=m1old
        m2new=m2old
        do iter=1,limit
            w=1.0_dp/(1.0_dp-m1old*cos(x)-m2old*sin(x))
            m1new=sum(w*cos(x))/sum(w)
            m2new=sum(w*sin(x))/sum(w)
            if (abs(m1new-m1old) < tol .and. abs(m2new-m2old) < tol) then
                fit%converged=.true.
                exit
            end if
            m1old=m1new
            m2old=m2new
        end do
        fit%iterations=min(iter,limit)
        mconst=hypot(m1new,m2new)
        if (mconst <= 100.0_dp*epsilon(1.0_dp)) then
            fit%rho=0.0_dp
            fit%mu=0.0_dp
        else
            mconst=min(mconst,1.0_dp-epsilon(1.0_dp))
            fit%rho=(1.0_dp-sqrt(max(0.0_dp,1.0_dp-mconst*mconst)))/mconst
            fit%mu=wrap_2pi(atan2(m2new,m1new))
        end if
    end function wrpcauchy_ml

    function circ_reg(alpha,theta,order,level) result(fit)
        real(dp), intent(in) :: alpha(:),theta(:)
        integer, intent(in), optional :: order
        real(dp), intent(in), optional :: level
        type(circ_reg_result) :: fit
        real(dp), allocatable :: xmat(:,:),wmat(:,:),bcos(:),bsin(:),bw1(:),bw2(:)
        real(dp), allocatable :: fcos(:),fsin(:),fw1(:),fw2(:),yr1(:),yr2(:),wr(:,:),h(:,:),rhs(:),sol(:)
        real(dp) :: lev,n1,d1,n2,d2,t1,t2
        integer :: n,ord,p,j,info,cdfree
        n=min(size(alpha),size(theta))
        ord=1
        if (present(order)) ord=order
        if (ord < 0) error stop "circ_reg: order must be nonnegative"
        lev=0.05_dp
        if (present(level)) lev=level
        p=2*ord+1
        if (n <= p+2) error stop "circ_reg: insufficient observations"
        allocate(xmat(n,p),wmat(n,2),bcos(p),bsin(p),bw1(p),bw2(p))
        allocate(fcos(n),fsin(n),fw1(n),fw2(n),yr1(n),yr2(n),wr(n,2),h(2,2),rhs(2),sol(2))
        xmat(:,1)=1.0_dp
        do j=1,ord
            xmat(:,1+j)=cos(real(j,dp)*alpha(1:n))
            xmat(:,1+ord+j)=sin(real(j,dp)*alpha(1:n))
        end do
        call ols_fit(xmat,cos(theta(1:n)),bcos,fcos,info)
        if (info /= 0) error stop "circ_reg: singular cosine regression"
        call ols_fit(xmat,sin(theta(1:n)),bsin,fsin,info)
        if (info /= 0) error stop "circ_reg: singular sine regression"
        wmat(:,1)=cos(real(ord+1,dp)*alpha(1:n))
        wmat(:,2)=sin(real(ord+1,dp)*alpha(1:n))
        call ols_fit(xmat,wmat(:,1),bw1,fw1,info)
        if (info /= 0) error stop "circ_reg: singular higher-order projection"
        call ols_fit(xmat,wmat(:,2),bw2,fw2,info)
        if (info /= 0) error stop "circ_reg: singular higher-order projection"
        wr(:,1)=wmat(:,1)-fw1
        wr(:,2)=wmat(:,2)-fw2
        h=matmul(transpose(wr),wr)
        yr1=cos(theta(1:n))-fcos
        yr2=sin(theta(1:n))-fsin
        rhs=matmul(transpose(wr),yr1)
        call solve_linear(h,rhs,sol,info)
        if (info /= 0) error stop "circ_reg: singular higher-order test"
        n1=dot_product(rhs,sol)
        d1=dot_product(yr1,yr1)
        rhs=matmul(transpose(wr),yr2)
        call solve_linear(h,rhs,sol,info)
        if (info /= 0) error stop "circ_reg: singular higher-order test"
        n2=dot_product(rhs,sol)
        d2=dot_product(yr2,yr2)
        cdfree=n-p
        t1=real(cdfree,dp)*n1/d1
        t2=real(cdfree,dp)*n2/d2
        fit%pvalues(1)=exp(-0.5_dp*t1)
        fit%pvalues(2)=exp(-0.5_dp*t2)
        fit%higher_order_significant=fit%pvalues(1)<lev .and. fit%pvalues(2)<lev
        fit%rho=sqrt((dot_product(fcos,fcos)+dot_product(fsin,fsin))/real(n,dp))
        allocate(fit%fitted(n),fit%residuals(n),fit%coef(p,2))
        fit%fitted=wrap_2pi(atan2(fsin,fcos))
        fit%residuals=wrap_2pi(theta(1:n)-atan2(fsin,fcos))
        fit%coef(:,1)=bcos
        fit%coef(:,2)=bsin
        fit%a_k=sum(cos(fit%residuals))/real(n,dp)
        fit%kappa=a1inv(fit%a_k)
    end function circ_reg

    function vm_bootstrap_ci(x,bias,alpha,reps) result(res)
        real(dp), intent(in) :: x(:)
        logical, intent(in), optional :: bias
        real(dp), intent(in), optional :: alpha
        integer, intent(in), optional :: reps
        type(vm_bootstrap_result) :: res
        real(dp), allocatable :: sample(:),sorted_mu(:)
        real(dp) :: aa,offset,maxgap,gap
        integer :: b,n,bb,i,maxidx,idx
        logical :: bc
        n=size(x)
        bb=1000
        if (present(reps)) bb=reps
        aa=0.05_dp
        if (present(alpha)) aa=alpha
        bc=.false.
        if (present(bias)) bc=bias
        allocate(res%mu_reps(bb),res%kappa_reps(bb),sample(n),sorted_mu(bb))
        do b=1,bb
            do i=1,n
                idx=1+int(randu()*real(n,dp))
                idx=min(idx,n)
                sample(i)=x(idx)
            end do
            res%mu_reps(b)=wrap_2pi(circ_mean(sample))
            res%kappa_reps(b)=est_kappa(sample,bc)
        end do
        sorted_mu=res%mu_reps
        call sort_local(sorted_mu)
        maxgap=sorted_mu(1)-sorted_mu(bb)+twopi
        maxidx=bb
        do i=1,bb-1
            gap=sorted_mu(i+1)-sorted_mu(i)
            if (gap > maxgap) then
                maxgap=gap
                maxidx=i
            end if
        end do
        if (maxidx /= bb) then
            offset=twopi-sorted_mu(maxidx+1)
            sorted_mu=wrap_2pi(sorted_mu+offset)
            call sort_local(sorted_mu)
            res%mu_ci(1)=quantile_type7(sorted_mu,aa/2.0_dp)-offset
            res%mu_ci(2)=quantile_type7(sorted_mu,1.0_dp-aa/2.0_dp)-offset
            res%mu_ci=wrap_2pi(res%mu_ci)
        else
            res%mu_ci(1)=quantile_type7(sorted_mu,aa/2.0_dp)
            res%mu_ci(2)=quantile_type7(sorted_mu,1.0_dp-aa/2.0_dp)
        end if
        res%kappa_ci(1)=quantile_type7(res%kappa_reps,aa/2.0_dp)
        res%kappa_ci(2)=quantile_type7(res%kappa_reps,1.0_dp-aa/2.0_dp)
    contains
        subroutine sort_local(v)
            real(dp), intent(inout) :: v(:)
            real(dp) :: key
            integer :: ii,jj
            do ii=2,size(v)
                key=v(ii)
                jj=ii-1
                do while (jj >= 1)
                    if (v(jj) <= key) exit
                    v(jj+1)=v(jj)
                    jj=jj-1
                end do
                v(jj+1)=key
            end do
        end subroutine sort_local
    end function vm_bootstrap_ci
end module circstats_models
