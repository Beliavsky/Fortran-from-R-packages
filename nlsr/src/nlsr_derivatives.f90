! SPDX-License-Identifier: GPL-2.0-only
module nlsr_derivatives
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use nlsr_kinds, only : dp
    use nlsr_types, only : residual_fn, jac_forward, jac_backward, jac_central, jac_richardson
    implicit none
    private
    public :: numerical_jacobian, residual_gradient, residual_sum_squares, jafwd, jaback, jacentral, jand

contains


    subroutine jafwd(resfn,par,resbest,jac,ndstep,bdmask,feval,ierr)
        procedure(residual_fn) :: resfn
        real(dp),intent(in) :: par(:),resbest(:),ndstep
        real(dp),intent(out) :: jac(:,:)
        integer,intent(in) :: bdmask(:)
        integer,intent(inout) :: feval
        integer,intent(out) :: ierr
        call numerical_jacobian(resfn,par,resbest,jac,jac_forward,ndstep,bdmask,feval,ierr)
    end subroutine jafwd

    subroutine jaback(resfn,par,resbest,jac,ndstep,bdmask,feval,ierr)
        procedure(residual_fn) :: resfn
        real(dp),intent(in) :: par(:),resbest(:),ndstep
        real(dp),intent(out) :: jac(:,:)
        integer,intent(in) :: bdmask(:)
        integer,intent(inout) :: feval
        integer,intent(out) :: ierr
        call numerical_jacobian(resfn,par,resbest,jac,jac_backward,ndstep,bdmask,feval,ierr)
    end subroutine jaback

    subroutine jacentral(resfn,par,resbest,jac,ndstep,bdmask,feval,ierr)
        procedure(residual_fn) :: resfn
        real(dp),intent(in) :: par(:),resbest(:),ndstep
        real(dp),intent(out) :: jac(:,:)
        integer,intent(in) :: bdmask(:)
        integer,intent(inout) :: feval
        integer,intent(out) :: ierr
        call numerical_jacobian(resfn,par,resbest,jac,jac_central,ndstep,bdmask,feval,ierr)
    end subroutine jacentral

    subroutine jand(resfn,par,resbest,jac,ndstep,bdmask,feval,ierr)
        procedure(residual_fn) :: resfn
        real(dp),intent(in) :: par(:),resbest(:),ndstep
        real(dp),intent(out) :: jac(:,:)
        integer,intent(in) :: bdmask(:)
        integer,intent(inout) :: feval
        integer,intent(out) :: ierr
        call numerical_jacobian(resfn,par,resbest,jac,jac_richardson,ndstep,bdmask,feval,ierr)
    end subroutine jand

    subroutine numerical_jacobian(resfn, par, res0, jac, method, ndstep, bdmask, feval, ierr)
        procedure(residual_fn) :: resfn
        real(dp), intent(in) :: par(:), res0(:)
        real(dp), intent(out) :: jac(:,:)
        integer, intent(in) :: method
        real(dp), intent(in) :: ndstep
        integer, intent(in), optional :: bdmask(:)
        integer, intent(inout) :: feval
        integer, intent(out) :: ierr
        integer :: p, n, j
        real(dp) :: h, h2
        real(dp), allocatable :: pp(:), pm(:), rp(:), rm(:), d1(:), d2(:)

        p = size(par); n = size(res0)
        jac = 0.0_dp
        ierr = 0
        if (size(jac,1) /= n .or. size(jac,2) /= p) then
            ierr = 1
            return
        end if
        allocate(pp(p),pm(p),rp(n),rm(n),d1(n),d2(n))
        do j = 1, p
            if (present(bdmask)) then
                if (bdmask(j) == 0) cycle
            end if
            h = ndstep*(abs(par(j))+ndstep)
            if (h <= 0.0_dp) h = sqrt(epsilon(1.0_dp))*(abs(par(j))+1.0_dp)
            select case(method)
            case(jac_forward)
                pp = par; pp(j) = pp(j)+h
                call invoke_residual(resfn,pp,rp,ierr); feval=feval+1
                if (ierr /= 0) return
                jac(:,j) = (rp-res0)/h
            case(jac_backward)
                pm = par; pm(j) = pm(j)-h
                call invoke_residual(resfn,pm,rm,ierr); feval=feval+1
                if (ierr /= 0) return
                jac(:,j) = (res0-rm)/h
            case(jac_richardson)
                ! Two-level centered Richardson extrapolation: O(h^4).
                pp=par; pm=par; pp(j)=pp(j)+h; pm(j)=pm(j)-h
                call invoke_residual(resfn,pp,rp,ierr); feval=feval+1
                if (ierr /= 0) return
                call invoke_residual(resfn,pm,rm,ierr); feval=feval+1
                if (ierr /= 0) return
                d1=(rp-rm)/(2.0_dp*h)
                h2=0.5_dp*h
                pp=par; pm=par; pp(j)=pp(j)+h2; pm(j)=pm(j)-h2
                call invoke_residual(resfn,pp,rp,ierr); feval=feval+1
                if (ierr /= 0) return
                call invoke_residual(resfn,pm,rm,ierr); feval=feval+1
                if (ierr /= 0) return
                d2=(rp-rm)/(2.0_dp*h2)
                jac(:,j)=d2+(d2-d1)/3.0_dp
            case default
                pp=par; pm=par; pp(j)=pp(j)+h; pm(j)=pm(j)-h
                call invoke_residual(resfn,pp,rp,ierr); feval=feval+1
                if (ierr /= 0) return
                call invoke_residual(resfn,pm,rm,ierr); feval=feval+1
                if (ierr /= 0) return
                jac(:,j)=(rp-rm)/(2.0_dp*h)
            end select
            if (any(.not. ieee_is_finite(jac(:,j)))) then
                ierr = 2
                return
            end if
        end do
    end subroutine numerical_jacobian

    subroutine residual_gradient(resfn, par, jac, grad, ss, ierr)
        procedure(residual_fn) :: resfn
        real(dp), intent(in) :: par(:), jac(:,:)
        real(dp), intent(out) :: grad(:), ss
        integer, intent(out) :: ierr
        real(dp), allocatable :: r(:)
        allocate(r(size(jac,1)))
        call invoke_residual(resfn,par,r,ierr)
        if (ierr /= 0) return
        grad = 2.0_dp*matmul(transpose(jac),r)
        ss = dot_product(r,r)
    end subroutine residual_gradient

    pure real(dp) function residual_sum_squares(residual) result(ss)
        real(dp), intent(in) :: residual(:)
        ss = dot_product(residual,residual)
    end function residual_sum_squares

    subroutine invoke_residual(fn,par,res,ierr)
        procedure(residual_fn) :: fn
        real(dp), intent(in) :: par(:)
        real(dp), intent(out) :: res(:)
        integer, intent(out) :: ierr
        call fn(par,res,ierr)
    end subroutine invoke_residual

end module nlsr_derivatives
