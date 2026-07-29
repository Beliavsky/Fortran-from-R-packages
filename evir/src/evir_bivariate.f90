! SPDX-License-Identifier: GPL-2.0-or-later
module evir_bivariate
    use evir_kinds, only : dp
    use evir_types, only : gpdbiv_fit_result, gpd_fit_result, evir_ok, &
        evir_invalid_input, evir_optimization_failed, evir_singular_hessian, evir_domain_error
    use evir_fitting, only : fit_gpd
    use evir_optimize, only : nelder_mead, numerical_hessian
    use evir_math, only : safe_nan
    implicit none
    private

    type :: gpdbiv_context
        real(dp), allocatable :: data1(:), data2(:)
        real(dp) :: u1 = 0.0_dp, u2 = 0.0_dp
        real(dp) :: lambda1 = 0.0_dp, lambda2 = 0.0_dp
        real(dp) :: marginal(4) = 0.0_dp
        logical :: global_fit = .false.
    end type gpdbiv_context

    public :: fit_gpdbiv, interpret_gpdbiv
    public :: logistic_exponent, bivariate_cdf, bivariate_survivor

contains

    function fit_gpdbiv(data1,data2,u1,u2,ne1,ne2,global_fit) result(out)
        real(dp), intent(in) :: data1(:),data2(:)
        real(dp), intent(in), optional :: u1,u2
        integer, intent(in), optional :: ne1,ne2
        logical, intent(in), optional :: global_fit
        type(gpdbiv_fit_result) :: out
        type(gpd_fit_result) :: margin1,margin2
        type(gpdbiv_context) :: ctx
        real(dp), allocatable :: x0(:),xbest(:),hess(:,:),cov(:,:)
        real(dp) :: fbest
        logical :: global,conv
        integer :: npar,hstatus,i,j

        if(size(data1)/=size(data2).or.size(data1)<3) then
            out%status=evir_invalid_input
            return
        end if
        if(present(u1).eqv.present(ne1)) then
            out%status=evir_invalid_input
            return
        end if
        if(present(u2).eqv.present(ne2)) then
            out%status=evir_invalid_input
            return
        end if
        if(present(u1)) then
            margin1=fit_gpd(data1,threshold=u1)
        else
            margin1=fit_gpd(data1,nextremes=ne1)
        end if
        if(present(u2)) then
            margin2=fit_gpd(data2,threshold=u2)
        else
            margin2=fit_gpd(data2,nextremes=ne2)
        end if
        if(margin1%status/=evir_ok.or.margin2%status/=evir_ok) then
            out%status=merge(margin1%status,margin2%status,margin1%status/=evir_ok)
            return
        end if
        global=.false.
        if(present(global_fit)) global=global_fit
        ctx%data1=data1; ctx%data2=data2
        ctx%u1=margin1%threshold; ctx%u2=margin2%threshold
        ctx%lambda1=real(margin1%n_exceed,dp)/real(size(data1),dp)
        ctx%lambda2=real(margin2%n_exceed,dp)/real(size(data2),dp)
        ctx%marginal=[margin1%xi,margin1%beta,margin2%xi,margin2%beta]
        ctx%global_fit=global
        npar=merge(5,1,global)
        allocate(x0(npar),xbest(npar),hess(npar,npar),cov(npar,npar))
        x0(1)=0.8_dp
        if(global) x0(2:5)=ctx%marginal
        call nelder_mead(gpdbiv_objective,ctx,x0,xbest,fbest,conv,initial_step=0.06_dp,max_iter=8000)
        call numerical_hessian(gpdbiv_objective,ctx,xbest,hess,cov,hstatus)

        out%global_fit=global
        out%converged=conv
        out%u1=ctx%u1; out%u2=ctx%u2
        out%ne1=margin1%n_exceed; out%ne2=margin2%n_exceed
        out%lambda1=ctx%lambda1; out%lambda2=ctx%lambda2
        out%alpha=xbest(1); out%nllh=fbest
        if(global) then
            out%par1=xbest(2:3); out%par2=xbest(4:5)
        else
            out%par1=ctx%marginal(1:2); out%par2=ctx%marginal(3:4)
        end if
        if(hstatus==0) then
            out%alpha_se=sqrt(max(cov(1,1),0.0_dp))
            if(global) then
                out%se1=[sqrt(max(cov(2,2),0.0_dp)),sqrt(max(cov(3,3),0.0_dp))]
                out%se2=[sqrt(max(cov(4,4),0.0_dp)),sqrt(max(cov(5,5),0.0_dp))]
            else
                out%se1=margin1%se; out%se2=margin2%se
            end if
        else
            out%status=evir_singular_hessian
            out%alpha_se=safe_nan(); out%se1=safe_nan(); out%se2=safe_nan()
        end if
        if(.not.conv) out%status=evir_optimization_failed

        allocate(out%data1(margin1%n_exceed),out%joint1(margin1%n_exceed))
        allocate(out%data2(margin2%n_exceed),out%joint2(margin2%n_exceed))
        i=0
        do j=1,size(data1)
            if(data1(j)>ctx%u1) then
                i=i+1; out%data1(i)=data1(j); out%joint1(i)=data2(j)>ctx%u2
            end if
        end do
        i=0
        do j=1,size(data2)
            if(data2(j)>ctx%u2) then
                i=i+1; out%data2(i)=data2(j); out%joint2(i)=data1(j)>ctx%u1
            end if
        end do
    end function fit_gpdbiv

    pure real(dp) function logistic_exponent(x,y,alpha) result(v)
        real(dp), intent(in) :: x,y,alpha
        if(x<=0.0_dp.or.y<=0.0_dp.or.alpha<=0.0_dp.or.alpha>1.0_dp) then
            v=safe_nan()
        else
            v=(x**(-1.0_dp/alpha)+y**(-1.0_dp/alpha))**alpha
        end if
    end function logistic_exponent

    pure real(dp) function bivariate_cdf(out,x,y) result(p)
        type(gpdbiv_fit_result), intent(in) :: out
        real(dp), intent(in) :: x,y
        real(dp) :: z1,z2
        if(x<out%u1.or.y<out%u2) then
            p=safe_nan(); return
        end if
        z1=z_function(x,out%u1,out%lambda1,out%par1(1),out%par1(2))
        z2=z_function(y,out%u2,out%lambda2,out%par2(1),out%par2(2))
        p=1.0_dp-logistic_exponent(z1,z2,out%alpha)
    end function bivariate_cdf

    pure real(dp) function bivariate_survivor(out,x,y) result(p)
        type(gpdbiv_fit_result), intent(in) :: out
        real(dp), intent(in) :: x,y
        real(dp) :: f1,f2,joint
        if(x<out%u1.or.y<out%u2) then
            p=safe_nan(); return
        end if
        f1=marginal_cdf(x,out%u1,out%lambda1,out%par1(1),out%par1(2))
        f2=marginal_cdf(y,out%u2,out%lambda2,out%par2(1),out%par2(2))
        joint=bivariate_cdf(out,x,y)
        p=1.0_dp-f1-f2+joint
    end function bivariate_survivor

    subroutine interpret_gpdbiv(out,x,y,probabilities,status)
        type(gpdbiv_fit_result), intent(in) :: out
        real(dp), intent(in) :: x,y
        real(dp), intent(out) :: probabilities(6)
        integer, intent(out), optional :: status
        real(dp) :: p1,p2,p12
        if(x<out%u1.or.y<out%u2) then
            probabilities=safe_nan()
            if(present(status)) status=evir_invalid_input
            return
        end if
        p1=1.0_dp-marginal_cdf(x,out%u1,out%lambda1,out%par1(1),out%par1(2))
        p2=1.0_dp-marginal_cdf(y,out%u2,out%lambda2,out%par2(1),out%par2(2))
        p12=bivariate_survivor(out,x,y)
        probabilities(1:4)=[p1,p2,p12,p1*p2]
        if(p1>0.0_dp) then
            probabilities(5)=p12/p1
        else
            probabilities(5)=safe_nan()
        end if
        if(p2>0.0_dp) then
            probabilities(6)=p12/p2
        else
            probabilities(6)=safe_nan()
        end if
        if(present(status)) status=evir_ok
    end subroutine interpret_gpdbiv

    function gpdbiv_objective(theta,context) result(f)
        real(dp), intent(in) :: theta(:)
        class(*), intent(in) :: context
        real(dp) :: f
        real(dp) :: alpha,xi1,sigma1,xi2,sigma2
        real(dp) :: z1,z2,k1,k2,v,v1,v2,term
        logical :: d1,d2
        integer :: i
        select type(context)
        type is(gpdbiv_context)
            alpha=theta(1)
            if(context%global_fit) then
                xi1=theta(2); sigma1=theta(3); xi2=theta(4); sigma2=theta(5)
            else
                xi1=context%marginal(1); sigma1=context%marginal(2)
                xi2=context%marginal(3); sigma2=context%marginal(4)
            end if
            if(alpha<=0.0_dp.or.alpha>=1.0_dp.or.sigma1<=0.0_dp.or.sigma2<=0.0_dp) then
                f=huge(1.0_dp)/100.0_dp; return
            end if
            f=0.0_dp
            do i=1,size(context%data1)
                d1=context%data1(i)>context%u1
                d2=context%data2(i)>context%u2
                if(d1) then
                    z1=z_function(context%data1(i),context%u1,context%lambda1,xi1,sigma1)
                    k1=k_function(context%data1(i),context%u1,context%lambda1,xi1,sigma1)
                else
                    z1=1.0_dp/context%lambda1
                    k1=0.0_dp
                end if
                if(d2) then
                    z2=z_function(context%data2(i),context%u2,context%lambda2,xi2,sigma2)
                    k2=k_function(context%data2(i),context%u2,context%lambda2,xi2,sigma2)
                else
                    z2=1.0_dp/context%lambda2
                    k2=0.0_dp
                end if
                if(z1<=0.0_dp.or.z2<=0.0_dp) then
                    f=huge(1.0_dp)/100.0_dp; return
                end if
                if(d1.and.d2) then
                    v2=v_second(z1,z2,alpha)
                    term=k1*k2*v2
                else if(d1) then
                    v1=v_first(z1,1.0_dp/context%lambda2,alpha)
                    term=k1*v1
                else if(d2) then
                    v1=v_first(z2,1.0_dp/context%lambda1,alpha)
                    term=k2*v1
                else
                    v=logistic_exponent(1.0_dp/context%lambda1,1.0_dp/context%lambda2,alpha)
                    term=1.0_dp-v
                end if
                if(term<=0.0_dp) then
                    f=huge(1.0_dp)/100.0_dp; return
                end if
                f=f-log(term)
            end do
        class default
            f=huge(1.0_dp)/100.0_dp
        end select
    end function gpdbiv_objective

    pure real(dp) function z_function(y,u,lambda,xi,sigma) result(z)
        real(dp), intent(in) :: y,u,lambda,xi,sigma
        real(dp) :: excess,t
        excess=max(y-u,0.0_dp)
        if(abs(xi)<=1.0e-8_dp) then
            z=exp(excess/sigma)/lambda
        else
            t=1.0_dp+xi*excess/sigma
            if(t<=0.0_dp) then
                z=-1.0_dp
            else
                z=t**(1.0_dp/xi)/lambda
            end if
        end if
    end function z_function

    pure real(dp) function k_function(y,u,lambda,xi,sigma) result(k)
        real(dp), intent(in) :: y,u,lambda,xi,sigma
        real(dp) :: z
        z=z_function(y,u,lambda,xi,sigma)
        if(z<=0.0_dp.or.sigma<=0.0_dp) then
            k=safe_nan()
        else
            k=-lambda**(-xi)*z**(1.0_dp-xi)/sigma
        end if
    end function k_function

    pure real(dp) function v_first(x,y,alpha) result(v1)
        real(dp), intent(in) :: x,y,alpha
        real(dp) :: s
        s=x**(-1.0_dp/alpha)+y**(-1.0_dp/alpha)
        v1=-x**(-1.0_dp/alpha-1.0_dp)*s**(alpha-1.0_dp)
    end function v_first

    pure real(dp) function v_second(x,y,alpha) result(v2)
        real(dp), intent(in) :: x,y,alpha
        real(dp) :: s
        s=x**(-1.0_dp/alpha)+y**(-1.0_dp/alpha)
        v2=-(alpha-1.0_dp)/alpha*(x*y)**(-1.0_dp/alpha-1.0_dp)*s**(alpha-2.0_dp)
    end function v_second

    pure real(dp) function marginal_cdf(x,u,lambda,xi,sigma) result(f)
        real(dp), intent(in) :: x,u,lambda,xi,sigma
        real(dp) :: t
        if(abs(xi)<=1.0e-8_dp) then
            f=1.0_dp-lambda*exp(-(x-u)/sigma)
        else
            t=1.0_dp+xi*(x-u)/sigma
            if(t<=0.0_dp) then
                f=merge(1.0_dp,0.0_dp,xi<0.0_dp)
            else
                f=1.0_dp-lambda*t**(-1.0_dp/xi)
            end if
        end if
    end function marginal_cdf

end module evir_bivariate
