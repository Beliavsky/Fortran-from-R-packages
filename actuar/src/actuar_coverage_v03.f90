module actuar_coverage_v03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use actuar_kinds, only: dp
    implicit none
    private
    public :: coverage_spec_t, coverage_cdf, coverage_pdf

    abstract interface
        function coverage_model_fun(x, par) result(v)
            import dp
            real(dp), intent(in) :: x
            real(dp), intent(in) :: par(:)
            real(dp) :: v
        end function coverage_model_fun
    end interface

    type :: coverage_spec_t
        real(dp) :: deductible = 0.0_dp
        logical :: franchise = .false.
        real(dp) :: limit = huge(1.0_dp)
        real(dp) :: coinsurance = 1.0_dp
        real(dp) :: inflation = 0.0_dp
        logical :: per_loss = .false.
    contains
        procedure :: valid => coverage_valid
        procedure :: maximum_payment
    end type coverage_spec_t
contains
    pure logical function coverage_valid(self) result(ok)
        class(coverage_spec_t),intent(in)::self
        ok=self%deductible>=0.0_dp .and. self%limit>=0.0_dp .and. &
           self%coinsurance>=0.0_dp .and. self%coinsurance<=1.0_dp .and. &
           self%inflation>=0.0_dp .and. self%limit>self%deductible
    end function coverage_valid

    pure real(dp) function maximum_payment(self) result(v)
        class(coverage_spec_t),intent(in)::self
        if(.not.ieee_is_finite(self%limit)) then
            v=huge(1.0_dp)
        else if(self%franchise) then
            v=self%coinsurance*self%limit
        else
            v=self%coinsurance*(self%limit-self%deductible)
        end if
    end function maximum_payment

    real(dp) function coverage_cdf(y,par,cdf,spec) result(v)
        real(dp),intent(in)::y,par(:)
        procedure(coverage_model_fun)::cdf
        type(coverage_spec_t),intent(in)::spec
        real(dp)::r,d,u,fd,sd,x,b1,b2
        logical::has_limit
        v=0.0_dp;if(.not.spec%valid())return
        if(y<0.0_dp)return
        r=1.0_dp+spec%inflation;d=spec%deductible/r;u=spec%limit/r
        has_limit=ieee_is_finite(spec%limit)
        if(spec%deductible>0.0_dp) then;fd=clip01(cdf(d,par));else;fd=clip01(cdf(0.0_dp,par));end if
        sd=max(0.0_dp,1.0_dp-fd)
        if(spec%franchise) then
            b1=spec%coinsurance*spec%deductible;b2=spec%maximum_payment()
            if(y<=b1) then
                if(spec%per_loss .and. spec%deductible>0.0_dp)v=fd
                return
            end if
            if(has_limit .and. y>=b2) then;v=1.0_dp;return;end if
            if(spec%coinsurance<=0.0_dp) then;v=1.0_dp;return;end if
            x=y/(spec%coinsurance*r)
        else
            b2=spec%maximum_payment()
            if(y==0.0_dp) then
                if(spec%per_loss .and. spec%deductible>0.0_dp)v=fd
                return
            end if
            if(has_limit .and. y>=b2) then;v=1.0_dp;return;end if
            if(spec%coinsurance<=0.0_dp) then;v=1.0_dp;return;end if
            x=(y/spec%coinsurance+spec%deductible)/r
        end if
        v=clip01(cdf(x,par))
        if(.not.spec%per_loss .and. spec%deductible>0.0_dp) then
            if(sd>0.0_dp) then;v=clip01((v-fd)/sd);else;v=1.0_dp;end if
        end if
    end function coverage_cdf

    real(dp) function coverage_pdf(y,par,pdf,cdf,spec) result(v)
        real(dp),intent(in)::y,par(:)
        procedure(coverage_model_fun)::pdf,cdf
        type(coverage_spec_t),intent(in)::spec
        real(dp)::r,d,u,fd,sd,su,x,b1,b2,jac
        logical::has_limit
        v=0.0_dp;if(.not.spec%valid() .or. y<0.0_dp)return
        r=1.0_dp+spec%inflation;d=spec%deductible/r;u=spec%limit/r
        has_limit=ieee_is_finite(spec%limit)
        if(spec%deductible>0.0_dp) then;fd=clip01(cdf(d,par));else;fd=clip01(cdf(0.0_dp,par));end if
        sd=max(0.0_dp,1.0_dp-fd)
        if(spec%franchise) then
            b1=spec%coinsurance*spec%deductible;b2=spec%maximum_payment()
            if(y==0.0_dp) then
                if(spec%per_loss .and. spec%deductible>0.0_dp)v=fd
                return
            end if
            if(y<=b1)return
            if(has_limit .and. y==b2) then
                su=max(0.0_dp,1.0_dp-clip01(cdf(u,par)));v=su
                if(.not.spec%per_loss .and. spec%deductible>0.0_dp .and. sd>0.0_dp)v=v/sd
                return
            end if
            if(has_limit .and. y>b2)return
            if(spec%coinsurance<=0.0_dp)return
            x=y/(spec%coinsurance*r)
        else
            b2=spec%maximum_payment()
            if(y==0.0_dp) then
                if(spec%per_loss .and. spec%deductible>0.0_dp)v=fd
                return
            end if
            if(has_limit .and. y==b2) then
                su=max(0.0_dp,1.0_dp-clip01(cdf(u,par)));v=su
                if(.not.spec%per_loss .and. spec%deductible>0.0_dp .and. sd>0.0_dp)v=v/sd
                return
            end if
            if(has_limit .and. y>b2)return
            if(spec%coinsurance<=0.0_dp)return
            x=(y/spec%coinsurance+spec%deductible)/r
        end if
        jac=spec%coinsurance*r
        if(jac<=0.0_dp)return
        v=max(0.0_dp,pdf(x,par))/jac
        if(.not.spec%per_loss .and. spec%deductible>0.0_dp) then
            if(sd>0.0_dp)then;v=v/sd;else;v=0.0_dp;end if
        end if
    end function coverage_pdf

    pure real(dp) function clip01(x) result(v)
        real(dp),intent(in)::x
        v=min(1.0_dp,max(0.0_dp,x))
    end function clip01
end module actuar_coverage_v03
