! SPDX-License-Identifier: GPL-2.0-only
module nlsr_models
    use nlsr_kinds, only : dp
    use nlsr_linalg, only : solve_linear
    implicit none
    private
    public :: logistic_value, logistic_jacobian, logistic_initial

contains

    pure subroutine logistic_value(input, asym, xmid, scal, value)
        real(dp), intent(in) :: input(:), asym, xmid, scal
        real(dp), intent(out) :: value(:)
        value = asym/(1.0_dp + exp((xmid-input)/scal))
    end subroutine logistic_value

    pure subroutine logistic_jacobian(input, asym, xmid, scal, jac)
        real(dp), intent(in) :: input(:), asym, xmid, scal
        real(dp), intent(out) :: jac(:,:)
        real(dp), allocatable :: e(:), den(:), xm(:)
        allocate(e(size(input)),den(size(input)),xm(size(input)))
        e=exp((xmid-input)/scal)
        den=1.0_dp+e
        jac(:,1)=1.0_dp/den
        xm=asym*e/scal/(den*den)
        jac(:,2)=-xm
        jac(:,3)=xm*((xmid-input)/scal)
    end subroutine logistic_jacobian

    subroutine logistic_initial(input,y,par,ok)
        real(dp), intent(in) :: input(:),y(:)
        real(dp), intent(out) :: par(3)
        logical, intent(out) :: ok
        real(dp), allocatable :: z(:),xx(:,:),rhs(:),coef(:)
        real(dp) :: asym
        integer :: n
        logical :: lsok
        n=size(y); par=0.0_dp; ok=.false.
        if (size(input)/=n .or. n<4) return
        asym=2.0_dp*maxval(max(y,1.0e-9_dp))
        allocate(z(n),xx(2,2),rhs(2),coef(2))
        z=log(asym/max(y,1.0e-9_dp)-1.0_dp)
        xx(1,1)=real(n,dp); xx(1,2)=sum(input); xx(2,1)=xx(1,2); xx(2,2)=dot_product(input,input)
        rhs(1)=sum(z); rhs(2)=dot_product(input,z)
        call solve_linear(xx,rhs,coef,lsok)
        if (.not. lsok .or. abs(coef(2))<=tiny(1.0_dp)) return
        par(1)=asym
        par(3)=-1.0_dp/coef(2)
        par(2)=coef(1)*par(3)
        ok=.true.
    end subroutine logistic_initial
end module nlsr_models
