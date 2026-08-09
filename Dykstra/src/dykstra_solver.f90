! SPDX-License-Identifier: GPL-2.0-or-later
module dykstra_solver
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use dykstra_kinds, only : dp
    use dykstra_linalg, only : symmetric_eigen_jacobi
    implicit none
    private

    type, public :: dykstra_result
        real(dp), allocatable :: solution(:)
        real(dp), allocatable :: unconstrained(:)
        real(dp) :: value = 0.0_dp
        integer :: iterations = 0
        logical :: converged = .false.
        integer :: status = 0
        character(len=:), allocatable :: message
    end type dykstra_result

    public :: dykstra_solve

contains

    subroutine dykstra_solve(dmat, dvec, amat, result, bvec, meq, factorized, maxit, eps)
        real(dp), intent(in) :: dmat(:,:), dvec(:), amat(:,:)
        type(dykstra_result), intent(out) :: result
        real(dp), intent(in), optional :: bvec(:)
        integer, intent(in), optional :: meq, maxit
        logical, intent(in), optional :: factorized
        real(dp), intent(in), optional :: eps

        real(dp), allocatable :: b(:), atrans(:,:), rinv(:,:), eval(:), evec(:,:)
        real(dp), allocatable :: rdiag(:), ddiag(:), gvec(:), acss(:)
        real(dp), allocatable :: beta_change(:,:), beta_solution(:), beta_old(:), beta_work(:)
        real(dp) :: user_eps, tol, ctol, maxb0, ai, shift, spectral_scale
        real(dp) :: symtol, matscale, value_q, value_p
        integer :: n, ncon, neq, max_cycles, i, nvals, einfo, iter
        logical :: fact, diag_flag, passive

        call initialize_result(result)
        n = size(dvec)
        if (size(dmat,1) /= n .or. size(dmat,2) /= n) then
            call fail(result, -1, "Inputs 'Dmat' and 'dvec' are incompatible.")
            return
        end if
        if (size(amat,1) /= n) then
            call fail(result, -2, "Input 'Amat' must have nrow(Amat) == length(dvec).")
            return
        end if
        ncon = size(amat,2)

        allocate(b(ncon))
        if (present(bvec)) then
            if (size(bvec) /= ncon) then
                call fail(result, -3, "Input 'bvec' must have length ncol(Amat).")
                return
            end if
            b = bvec
        else
            b = 0.0_dp
        end if

        neq = 0
        if (present(meq)) neq = meq
        if (neq < 0 .or. neq > ncon) then
            call fail(result, -4, "Input 'meq' must be between 0 and length(bvec).")
            return
        end if

        max_cycles = 30 * max(1, n)
        if (present(maxit)) max_cycles = maxit
        if (max_cycles < 1) then
            call fail(result, -5, "Input 'maxit' must be a positive integer.")
            return
        end if

        user_eps = epsilon(1.0_dp) * real(max(1,n), dp)
        if (present(eps)) user_eps = eps
        if (user_eps < 0.0_dp) then
            call fail(result, -6, "Input 'eps' must be a non-negative scalar.")
            return
        end if
        tol = epsilon(1.0_dp) * real(max(1,n), dp)

        fact = .false.
        if (present(factorized)) fact = factorized
        if (.not. fact) then
            matscale = max(1.0_dp, maxval(abs(dmat)))
            symtol = 100.0_dp * epsilon(1.0_dp) * matscale
            if (maxval(abs(dmat - transpose(dmat))) > symtol) then
                call fail(result, -7, "Input 'Dmat' must be symmetric when factorized is false.")
                return
            end if
        end if

        diag_flag = max_upper_offdiag(dmat) <= tol
        allocate(gvec(n), result%unconstrained(n), atrans(n,ncon))

        if (diag_flag) then
            allocate(rdiag(n))
            if (fact) then
                do i = 1, n
                    rdiag(i) = dmat(i,i)
                end do
            else
                allocate(ddiag(n))
                do i = 1, n
                    ddiag(i) = dmat(i,i)
                end do
                if (any(ddiag < -tol)) then
                    call fail(result, -8, "Input 'Dmat' must be positive definite or semidefinite.")
                    return
                end if
                spectral_scale = max(1.0_dp, maxval(ddiag))
                nvals = count(ddiag > tol * spectral_scale)
                if (nvals < n) then
                    shift = tol * spectral_scale - minval(ddiag)
                    if (shift <= 0.0_dp) shift = tol * spectral_scale
                    ddiag = ddiag + shift
                end if
                if (any(ddiag <= 0.0_dp)) then
                    call fail(result, -9, "Unable to regularize diagonal Dmat to positive definiteness.")
                    return
                end if
                rdiag = 1.0_dp / sqrt(ddiag)
            end if

            gvec = rdiag * dvec
            do i = 1, n
                atrans(i,:) = rdiag(i) * amat(i,:)
            end do
            ! Preserve the upstream Dykstra 1.0-0 output formula exactly.
            result%unconstrained = dvec / (rdiag * rdiag)
        else
            if (fact) then
                allocate(rinv(n,n))
                rinv = dmat
            else
                call symmetric_eigen_jacobi(dmat, eval, evec, einfo)
                if (einfo < 0) then
                    call fail(result, -10, "Symmetric eigendecomposition failed.")
                    return
                end if
                if (any(eval < -tol)) then
                    call fail(result, -8, "Input 'Dmat' must be positive definite or semidefinite.")
                    return
                end if
                spectral_scale = max(1.0_dp, maxval(eval))
                nvals = count(eval > tol * spectral_scale)
                if (nvals < n) then
                    shift = tol * spectral_scale - minval(eval)
                    if (shift <= 0.0_dp) shift = tol * spectral_scale
                    eval = eval + shift
                end if
                if (any(eval <= 0.0_dp)) then
                    call fail(result, -11, "Unable to regularize Dmat to positive definiteness.")
                    return
                end if
                allocate(rinv(n,n))
                do i = 1, n
                    rinv(:,i) = evec(:,i) / sqrt(eval(i))
                end do
            end if
            gvec = matmul(transpose(rinv), dvec)
            atrans = matmul(transpose(rinv), amat)
            result%unconstrained = matmul(rinv, matmul(transpose(rinv), dvec))
        end if

        allocate(acss(ncon))
        if (ncon > 0) acss = sum(atrans * atrans, dim=1)
        do i = 1, ncon
            if (acss(i) <= tiny(1.0_dp)) then
                call fail(result, -12, "Constraint matrix contains a numerically zero column.")
                return
            end if
        end do

        maxb0 = 0.0_dp
        if (n > 0) maxb0 = maxval(abs(result%unconstrained))
        if (maxb0 > 1.0_dp) user_eps = user_eps * maxb0

        allocate(beta_change(n,ncon), beta_solution(n), beta_old(n), beta_work(n))
        beta_change = 0.0_dp
        beta_solution = gvec
        beta_old = gvec
        iter = 0
        ctol = user_eps + 1.0_dp

        do while (ctol > user_eps .and. iter < max_cycles)
            do i = 1, ncon
                beta_work = beta_solution - beta_change(:,i)
                ai = dot_product(atrans(:,i), beta_work)
                if (i <= neq) then
                    passive = ai <= b(i) .and. ai >= b(i)
                else
                    passive = ai >= b(i)
                end if
                if (passive) then
                    beta_change(:,i) = 0.0_dp
                    beta_solution = beta_work
                else
                    beta_change(:,i) = (b(i) - ai) * atrans(:,i) / acss(i)
                    beta_solution = beta_work + beta_change(:,i)
                end if
            end do
            ctol = maxval(abs(beta_solution - beta_old))
            beta_old = beta_solution
            iter = iter + 1
        end do

        allocate(result%solution(n))
        if (diag_flag) then
            result%solution = rdiag * beta_solution
        else
            result%solution = matmul(rinv, beta_solution)
        end if

        result%iterations = iter
        result%converged = ctol <= user_eps
        if (result%converged) then
            result%status = 0
            result%message = "converged"
        else
            result%status = 1
            result%message = "maximum number of cycles reached"
        end if

        if (fact) then
            result%value = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            value_q = 0.5_dp * dot_product(result%solution, matmul(dmat, result%solution))
            value_p = dot_product(result%solution, dvec)
            result%value = value_q - value_p
        end if
    end subroutine dykstra_solve

    subroutine initialize_result(result)
        type(dykstra_result), intent(out) :: result
        result%value = 0.0_dp
        result%iterations = 0
        result%converged = .false.
        result%status = 0
        result%message = ""
    end subroutine initialize_result

    subroutine fail(result, status, message)
        type(dykstra_result), intent(inout) :: result
        integer, intent(in) :: status
        character(len=*), intent(in) :: message
        result%status = status
        result%message = message
        result%converged = .false.
        result%value = ieee_value(0.0_dp, ieee_quiet_nan)
    end subroutine fail

    pure real(dp) function max_upper_offdiag(a) result(ans)
        real(dp), intent(in) :: a(:,:)
        integer :: i, j
        ans = 0.0_dp
        do j = 2, size(a,2)
            do i = 1, min(j - 1, size(a,1))
                ans = max(ans, abs(a(i,j)))
            end do
        end do
    end function max_upper_offdiag

end module dykstra_solver
