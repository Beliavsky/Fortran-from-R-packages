! SPDX-License-Identifier: GPL-2.0-or-later
! Generic delta-method and score-information utilities corresponding to msm's inference helpers.
module msm_inference
    use msm_kinds, only : dp
    implicit none
    private
    public :: delta_method, outer_product_information
    abstract interface
        subroutine vector_function(x,y)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), allocatable, intent(out) :: y(:)
        end subroutine vector_function
    end interface
contains
    subroutine delta_method(fun, theta, covariance, estimate, out_covariance, se, relstep)
        procedure(vector_function) :: fun
        real(dp), intent(in) :: theta(:),covariance(:,:)
        real(dp), allocatable, intent(out) :: estimate(:),out_covariance(:,:),se(:)
        real(dp), intent(in), optional :: relstep
        real(dp), allocatable :: yp(:),ym(:),xp(:),xm(:),jac(:,:)
        real(dp) :: h,rs
        integer :: p,m,j
        p=size(theta); if(any(shape(covariance)/=[p,p])) error stop "delta_method: covariance shape"
        call fun(theta,estimate); m=size(estimate); allocate(jac(m,p),xp(p),xm(p))
        rs=sqrt(epsilon(1.0_dp)); if(present(relstep)) rs=relstep
        do j=1,p
            h=rs*max(1.0_dp,abs(theta(j))); xp=theta; xm=theta; xp(j)=xp(j)+h; xm(j)=xm(j)-h
            call fun(xp,yp); call fun(xm,ym); if(size(yp)/=m.or.size(ym)/=m) error stop "delta_method: output size changed"
            jac(:,j)=(yp-ym)/(2.0_dp*h)
        end do
        out_covariance=matmul(jac,matmul(covariance,transpose(jac)))
        allocate(se(m)); do j=1,m; se(j)=sqrt(max(0.0_dp,out_covariance(j,j))); end do
    end subroutine delta_method

    function outer_product_information(scores) result(info)
        real(dp), intent(in) :: scores(:,:)
        real(dp), allocatable :: info(:,:)
        info=matmul(transpose(scores),scores)
    end function outer_product_information
end module msm_inference
