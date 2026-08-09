! SPDX-License-Identifier: GPL-3.0-only
module roi_qpoases
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
    use qpoases, only : dp, qpoases_options, qpoases_result, solve_qproblem, &
        hst_zero, hst_identity, hst_unknown
    implicit none
    private
    public :: roi_solve_qp, infer_hessian_type

contains

    integer function infer_hessian_type(q) result(ht)
        real(dp), intent(in) :: q(:,:)
        real(dp) :: tol
        integer :: i, j, n
        logical :: ident

        n = size(q,1)
        tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp,maxval(abs(q)))
        if (maxval(abs(q)) <= tol) then
            ht = hst_zero
            return
        end if
        ident = size(q,2) == n
        if (ident) then
            do j = 1, n
                do i = 1, n
                    if (i == j) then
                        if (abs(q(i,j)-1.0_dp) > tol) ident = .false.
                    else
                        if (abs(q(i,j)) > tol) ident = .false.
                    end if
                end do
            end do
        end if
        if (ident) then
            ht = hst_identity
        else
            ht = hst_unknown
        end if
    end function infer_hessian_type

    subroutine roi_solve_qp(q, linear, constraints, direction, rhs, result, &
                            maximum, lower, upper, options, hessian_type, max_nwsr)
        real(dp), intent(in) :: q(:,:), linear(:), constraints(:,:), rhs(:)
        character(len=*), intent(in) :: direction(:)
        type(qpoases_result), intent(out) :: result
        logical, intent(in), optional :: maximum
        real(dp), intent(in), optional :: lower(:), upper(:)
        type(qpoases_options), intent(in), optional :: options
        integer, intent(in), optional :: hessian_type, max_nwsr

        real(dp), allocatable :: h(:,:), g(:), lb(:), ub(:), lba(:), uba(:)
        real(dp) :: pinf, ninf
        logical :: maximize
        integer :: n, m, i, ht, nwsr

        n = size(linear)
        m = size(rhs)
        allocate(h(n,n),g(n),lb(n),ub(n),lba(m),uba(m))
        pinf = ieee_value(1.0_dp,ieee_positive_inf)
        ninf = ieee_value(1.0_dp,ieee_negative_inf)

        h = q
        g = linear
        maximize = .false.
        if (present(maximum)) maximize = maximum
        if (maximize) then
            h = -h
            g = -g
        end if

        ! ROI's default variable lower bound is zero.
        lb = 0.0_dp
        ub = pinf
        if (present(lower)) lb = lower
        if (present(upper)) ub = upper

        lba = ninf
        uba = pinf
        do i = 1, m
            select case (trim(direction(i)))
            case ("==","=")
                lba(i) = rhs(i)
                uba(i) = rhs(i)
            case ("<=","<")
                uba(i) = rhs(i)
            case (">=",">")
                lba(i) = rhs(i)
            case default
                result%status = 3
                return
            end select
        end do

        ht = infer_hessian_type(h)
        if (present(hessian_type)) ht = hessian_type
        nwsr = 2000
        if (present(max_nwsr)) nwsr = max_nwsr

        if (present(options)) then
            call solve_qproblem(h,g,constraints,lb,ub,lba,uba,result,options,ht,nwsr)
        else
            call solve_qproblem(h,g,constraints,lb,ub,lba,uba,result, &
                                 hessian_type=ht,max_nwsr=nwsr)
        end if

        if (allocated(result%x)) then
            result%objval = 0.5_dp * dot_product(result%x,matmul(q,result%x)) + &
                            dot_product(linear,result%x)
        end if
    end subroutine roi_solve_qp
end module roi_qpoases
