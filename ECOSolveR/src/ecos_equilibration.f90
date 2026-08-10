! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_equilibration
    use ecos_types, only : dp, ecos_problem, ecos_result, ecos_settings
    use ecos_sparse, only : csc_to_csr
    implicit none
    private

    type, public :: ecos_scaling
        logical :: active = .false.
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: eq(:)
        real(dp), allocatable :: cone(:)
    end type ecos_scaling

    public :: equilibrate_problem_sparse, unscale_result, scale_warm_vectors

contains

    pure real(dp) function clipped_inverse_sqrt(v, lo, hi) result(s)
        real(dp), intent(in) :: v, lo, hi
        real(dp) :: vv
        vv = max(v, 1.0e-16_dp)
        s = 1.0_dp/sqrt(vv)
        s = min(hi, max(lo, s))
    end function clipped_inverse_sqrt

    subroutine column_max_norm(prob, norm)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(out) :: norm(:)
        integer :: j, k
        norm = 0.0_dp
        do j = 1, prob%g_csc%ncol
            do k = prob%g_csc%colptr(j), prob%g_csc%colptr(j+1)-1
                norm(j) = max(norm(j), abs(prob%g_csc%values(k)))
            end do
        end do
        do j = 1, prob%a_csc%ncol
            do k = prob%a_csc%colptr(j), prob%a_csc%colptr(j+1)-1
                norm(j) = max(norm(j), abs(prob%a_csc%values(k)))
            end do
        end do
    end subroutine column_max_norm

    subroutine scale_columns(prob, s)
        type(ecos_problem), intent(inout) :: prob
        real(dp), intent(in) :: s(:)
        integer :: j, k
        do j = 1, prob%g_csc%ncol
            do k = prob%g_csc%colptr(j), prob%g_csc%colptr(j+1)-1
                prob%g_csc%values(k) = prob%g_csc%values(k)*s(j)
            end do
        end do
        do j = 1, prob%a_csc%ncol
            do k = prob%a_csc%colptr(j), prob%a_csc%colptr(j+1)-1
                prob%a_csc%values(k) = prob%a_csc%values(k)*s(j)
            end do
        end do
        prob%c = prob%c*s
    end subroutine scale_columns

    subroutine equality_row_scale(prob, stg, total)
        type(ecos_problem), intent(inout) :: prob
        type(ecos_settings), intent(in) :: stg
        real(dp), intent(inout) :: total(:)
        real(dp), allocatable :: rn(:), sf(:)
        integer :: i, k
        if (prob%neq() == 0) return
        allocate(rn(prob%neq()), sf(prob%neq()))
        rn = 0.0_dp
        do i = 1, prob%a_csr%nrow
            do k = prob%a_csr%rowptr(i), prob%a_csr%rowptr(i+1)-1
                rn(i) = max(rn(i), abs(prob%a_csr%values(k)))
            end do
            sf(i) = clipped_inverse_sqrt(rn(i), stg%equilibration_min, stg%equilibration_max)
        end do
        do k = 1, size(prob%a_csc%values)
            i = prob%a_csc%rowind(k)
            prob%a_csc%values(k) = prob%a_csc%values(k)*sf(i)
        end do
        prob%b = prob%b*sf
        total = total*sf
    end subroutine equality_row_scale

    subroutine cone_block_scales(prob, stg, sf)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        real(dp), intent(out) :: sf(:)
        real(dp), allocatable :: rn(:)
        real(dp) :: vmax, fac
        integer :: i, k, row, iq, qn, j
        allocate(rn(prob%ncone()))
        rn = 0.0_dp
        do i = 1, prob%g_csr%nrow
            do k = prob%g_csr%rowptr(i), prob%g_csr%rowptr(i+1)-1
                rn(i) = max(rn(i), abs(prob%g_csr%values(k)))
            end do
        end do
        sf = 1.0_dp
        do i = 1, prob%dims%l
            sf(i) = clipped_inverse_sqrt(rn(i), stg%equilibration_min, stg%equilibration_max)
        end do
        row = prob%dims%l + 1
        if (allocated(prob%dims%q)) then
            do iq = 1, size(prob%dims%q)
                qn = prob%dims%q(iq)
                vmax = 0.0_dp
                do j = row, row+qn-1
                    vmax = max(vmax, rn(j))
                end do
                fac = clipped_inverse_sqrt(vmax, stg%equilibration_min, stg%equilibration_max)
                sf(row:row+qn-1) = fac
                row = row + qn
            end do
        end if
        do iq = 1, prob%dims%e
            vmax = max(rn(row), max(rn(row+1), rn(row+2)))
            fac = clipped_inverse_sqrt(vmax, stg%equilibration_min, stg%equilibration_max)
            sf(row:row+2) = fac
            row = row + 3
        end do
    end subroutine cone_block_scales

    subroutine cone_row_scale(prob, stg, total)
        type(ecos_problem), intent(inout) :: prob
        type(ecos_settings), intent(in) :: stg
        real(dp), intent(inout) :: total(:)
        real(dp), allocatable :: sf(:)
        integer :: k, i
        if (prob%ncone() == 0) return
        allocate(sf(prob%ncone()))
        call cone_block_scales(prob, stg, sf)
        do k = 1, size(prob%g_csc%values)
            i = prob%g_csc%rowind(k)
            prob%g_csc%values(k) = prob%g_csc%values(k)*sf(i)
        end do
        prob%h = prob%h*sf
        total = total*sf
    end subroutine cone_row_scale

    subroutine equilibrate_problem_sparse(prob, stg, scaled, scale, info)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        type(ecos_problem), intent(out) :: scaled
        type(ecos_scaling), intent(out) :: scale
        integer, intent(out) :: info
        real(dp), allocatable :: cn(:), cs(:)
        integer :: n, p, m, it

        info = 0
        scaled = prob
        n = prob%nvar(); p = prob%neq(); m = prob%ncone()
        allocate(scale%x(n), scale%eq(p), scale%cone(m))
        scale%x = 1.0_dp; scale%eq = 1.0_dp; scale%cone = 1.0_dp
        if (.not. stg%equilibrate .or. .not. prob%sparse_storage) then
            scale%active = .false.
            return
        end if
        allocate(cn(n), cs(n))
        do it = 1, max(0, stg%equilibration_iters)
            call csc_to_csr(scaled%g_csc, scaled%g_csr)
            call csc_to_csr(scaled%a_csc, scaled%a_csr)
            call column_max_norm(scaled, cn)
            cs = 1.0_dp
            where (cn > 0.0_dp)
                cs = 1.0_dp/sqrt(cn)
            end where
            cs = min(stg%equilibration_max, max(stg%equilibration_min, cs))
            call scale_columns(scaled, cs)
            scale%x = scale%x*cs
            call csc_to_csr(scaled%g_csc, scaled%g_csr)
            call csc_to_csr(scaled%a_csc, scaled%a_csr)
            call equality_row_scale(scaled, stg, scale%eq)
            call cone_row_scale(scaled, stg, scale%cone)
        end do
        call csc_to_csr(scaled%g_csc, scaled%g_csr)
        call csc_to_csr(scaled%a_csc, scaled%a_csr)
        scale%active = .true.
    end subroutine equilibrate_problem_sparse

    subroutine unscale_result(scale, result)
        type(ecos_scaling), intent(in) :: scale
        type(ecos_result), intent(inout) :: result
        if (.not. scale%active) return
        if (allocated(result%x)) result%x = scale%x*result%x
        if (allocated(result%y)) then
            if (size(result%y) == size(scale%eq)) result%y = scale%eq*result%y
        end if
        if (allocated(result%s)) then
            if (size(result%s) == size(scale%cone)) result%s = result%s/max(scale%cone, 1.0e-300_dp)
        end if
        if (allocated(result%z)) then
            if (size(result%z) == size(scale%cone)) result%z = scale%cone*result%z
        end if
    end subroutine unscale_result

    subroutine scale_warm_vectors(scale, x, y)
        type(ecos_scaling), intent(in) :: scale
        real(dp), intent(inout) :: x(:), y(:)
        if (.not. scale%active) return
        if (size(x) == size(scale%x)) x = x/max(scale%x, 1.0e-300_dp)
        if (size(y) == size(scale%eq)) y = y/max(scale%eq, 1.0e-300_dp)
    end subroutine scale_warm_vectors

end module ecos_equilibration
