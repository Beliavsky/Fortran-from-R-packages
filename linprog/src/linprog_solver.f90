! SPDX-License-Identifier: GPL-2.0-or-later
module linprog_solver
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
        ieee_positive_inf, ieee_negative_inf
    use linprog_types
    use lpsolve, only: lpsolve_solve_lp => solve_lp, lpsolve_control => lp_control, &
        lpsolve_result => lp_result, LP_MIN, LP_MAX, LP_LE, LP_GE, LP_EQ, LP_OPTIMAL
    implicit none
    private

    interface solveLP
        module procedure solveLP_char
    end interface solveLP

    public :: solveLP, solveLP_dirs, constraint_direction

contains

    pure integer function constraint_direction(text) result(dir)
        character(len=*), intent(in) :: text
        select case (trim(adjustl(text)))
        case ('<=', '<')
            dir = LINPROG_LE
        case ('>=', '>')
            dir = LINPROG_GE
        case ('=', '==')
            dir = LINPROG_EQ
        case default
            dir = 99
        end select
    end function constraint_direction

    recursive subroutine solveLP_char(cvec, bvec, amat, result, control, const_dir)
        real(dp), intent(in) :: cvec(:), bvec(:), amat(:,:)
        type(linprog_result), intent(out) :: result
        type(linprog_control), intent(in), optional :: control
        character(len=*), intent(in), optional :: const_dir(:)
        integer, allocatable :: idir(:)
        integer :: i

        allocate(idir(size(bvec)))
        idir = LINPROG_LE
        if (present(const_dir)) then
            if (size(const_dir) /= size(bvec)) error stop "solveLP: const_dir has wrong size"
            do i = 1, size(bvec)
                idir(i) = constraint_direction(const_dir(i))
                if (idir(i) == 99) error stop "solveLP: invalid constraint direction"
            end do
        end if
        call solveLP_dirs(cvec, bvec, amat, result, control, idir)
    end subroutine solveLP_char

    recursive subroutine solveLP_dirs(cvec, bvec, amat, result, control, const_dir)
        real(dp), intent(in) :: cvec(:), bvec(:), amat(:,:)
        type(linprog_result), intent(out) :: result
        type(linprog_control), intent(in), optional :: control
        integer, intent(in), optional :: const_dir(:)
        type(linprog_control) :: ctl, dctl
        integer, allocatable :: idir(:), ddir(:)
        real(dp), allocatable :: dc(:), db(:), da(:,:)
        type(linprog_result) :: dres
        integer :: i, sgn

        ctl = linprog_control()
        if (present(control)) ctl = control
        if (size(amat,1) /= size(bvec) .or. size(amat,2) /= size(cvec)) then
            error stop "solveLP: matrix A has inconsistent dimensions"
        end if
        allocate(idir(size(bvec)))
        idir = LINPROG_LE
        if (present(const_dir)) then
            if (size(const_dir) /= size(bvec)) error stop "solveLP: const_dir has wrong size"
            idir = const_dir
        end if
        if (any(idir /= LINPROG_LE .and. idir /= LINPROG_GE .and. idir /= LINPROG_EQ)) then
            error stop "solveLP: invalid constraint direction"
        end if

        result%maximum = ctl%maximum
        result%used_lpsolve = ctl%use_lpsolve
        result%solved_dual = ctl%solve_dual
        result%maxiter = ctl%maxiter

        if (ctl%use_lpsolve) then
            call solve_with_lpsolve(cvec, bvec, amat, idir, ctl, result)
        else
            call solve_with_tableau(cvec, bvec, amat, idir, ctl, result)
        end if

        if (result%status == LINPROG_SUCCESS .and. allocated(result%con_free)) then
            if (minval(result%con_free) < -ctl%tol) result%status = LINPROG_CONSTRAINT_VIOLATION
        end if

        if (ctl%solve_dual .and. result%status == LINPROG_SUCCESS) then
            if (any(idir == LINPROG_EQ)) then
                error stop "solveLP: dual solve with equality constraints is unsupported"
            end if
            allocate(dc(size(bvec)), db(size(cvec)), da(size(cvec),size(bvec)))
            allocate(ddir(size(cvec)))
            sgn = merge(-1, 1, ctl%maximum)
            dc = bvec * real(idir,dp) * real(sgn,dp)
            db = cvec
            do i = 1, size(bvec)
                da(:,i) = amat(i,:) * real(idir(i),dp) * real(sgn,dp)
            end do
            if (ctl%maximum) then
                ddir = LINPROG_GE
            else
                ddir = LINPROG_LE
            end if
            dctl = ctl
            dctl%maximum = .not. ctl%maximum
            dctl%solve_dual = .false.
            dctl%tol = ctl%dualtol
            call solveLP_dirs(dc, db, da, dres, dctl, ddir)
            result%dual_status = dres%status
            if (dres%status == LINPROG_SUCCESS) then
                if (allocated(result%con_dual)) then
                    allocate(result%con_dual_primal(size(result%con_dual)))
                    result%con_dual_primal = result%con_dual
                end if
                if (allocated(result%con_dual)) deallocate(result%con_dual)
                allocate(result%con_dual(size(bvec)))
                result%con_dual = dres%solution
            else
                result%status = LINPROG_DUAL_FAILED
            end if
        end if
    end subroutine solveLP_dirs

    subroutine solve_with_lpsolve(cvec, bvec, amat, idir, ctl, result)
        real(dp), intent(in) :: cvec(:), bvec(:), amat(:,:)
        integer, intent(in) :: idir(:)
        type(linprog_control), intent(in) :: ctl
        type(linprog_result), intent(inout) :: result
        type(lpsolve_result) :: lr
        type(lpsolve_control) :: lc
        integer, allocatable :: sense(:)
        integer :: i, ncon
        real(dp), allocatable :: actual(:)

        ncon = size(bvec)
        allocate(sense(ncon))
        do i = 1, ncon
            select case (idir(i))
            case (LINPROG_LE)
                sense(i) = LP_LE
            case (LINPROG_GE)
                sense(i) = LP_GE
            case default
                sense(i) = LP_EQ
            end select
        end do
        lc = lpsolve_control()
        lc%max_simplex_iter = max(ctl%maxiter, 1)
        lc%feasibility_tol = min(ctl%tol, 1.0e-9_dp)
        call lpsolve_solve_lp(merge(LP_MAX, LP_MIN, ctl%maximum), cvec, amat, &
            sense, bvec, lr, lc)
        result%lp_status = lr%status
        if (lr%status == LP_OPTIMAL) then
            if (minval(lr%solution) < -ctl%tol) result%lp_status = 7
            allocate(actual(ncon))
            actual = matmul(amat, lr%solution)
            do i = 1, ncon
                if ((bvec(i)-actual(i))*real(idir(i),dp) > ctl%tol) result%lp_status = 3
            end do
        end if
        if (result%lp_status /= 0) then
            result%status = LINPROG_LPSOLVE_FAILED
            return
        end if

        result%status = LINPROG_SUCCESS
        result%opt = lr%objective
        allocate(result%solution(size(cvec)))
        result%solution = lr%solution
        allocate(result%con_actual(ncon), result%con_bvec(ncon), result%con_free(ncon))
        allocate(result%con_dir(ncon))
        result%con_actual = round_vector(matmul(amat, result%solution), ctl%zero)
        result%con_bvec = bvec
        result%con_dir = idir
        result%con_free = round_vector(bvec-result%con_actual, ctl%zero)
        where (idir == LINPROG_GE) result%con_free = -result%con_free
        where (idir == LINPROG_EQ) result%con_free = -abs(result%con_free)
    end subroutine solve_with_lpsolve

    subroutine solve_with_tableau(cvec, bvec, amat, idir, ctl, result)
        real(dp), intent(in) :: cvec(:), bvec(:), amat(:,:)
        integer, intent(in) :: idir(:)
        type(linprog_control), intent(in) :: ctl
        type(linprog_result), intent(inout) :: result
        integer :: nvar, ncon, nbase, nart, i, k, pcol, prow, rhscol
        integer, allocatable :: basis(:), artrow(:)
        real(dp), allocatable :: tab(:,:), tab2(:,:), cvec2(:)
        logical :: ok

        nvar = size(cvec)
        ncon = size(bvec)
        nbase = nvar + ncon
        allocate(cvec2(nbase), basis(ncon))
        cvec2 = 0.0_dp
        cvec2(1:nvar) = cvec
        allocate(tab(ncon+1,nbase+1))
        tab = 0.0_dp
        do i = 1, ncon
            tab(i,1:nvar) = -amat(i,:) * real(idir(i),dp)
            tab(i,nvar+i) = 1.0_dp
            tab(i,nbase+1) = -bvec(i) * real(idir(i),dp)
            basis(i) = nvar+i
        end do
        tab(ncon+1,1:nbase) = cvec2 * merge(-1.0_dp, 1.0_dp, ctl%maximum)

        nart = count(tab(1:ncon,nbase+1) < 0.0_dp)
        result%iter1 = 0
        if (nart > 0) then
            allocate(tab2(ncon+2,nbase+nart+1), artrow(ncon))
            tab2 = 0.0_dp
            artrow = 0
            tab2(1:ncon+1,1:nbase) = tab(:,1:nbase)
            tab2(1:ncon+1,nbase+nart+1) = tab(:,nbase+1)
            k = 0
            do i = 1, ncon
                if (tab(i,nbase+1) < 0.0_dp) then
                    tab2(i,1:nbase) = -tab2(i,1:nbase)
                    tab2(i,nbase+nart+1) = -tab2(i,nbase+nart+1)
                    k = k + 1
                    artrow(i) = nbase+k
                    tab2(i,nbase+k) = 1.0_dp
                    tab2(ncon+2,nbase+k) = 1.0_dp
                end if
            end do
            do i = 1, ncon
                if (artrow(i) > 0) then
                    tab2(ncon+2,1:nbase+nart) = tab2(ncon+2,1:nbase+nart) - &
                        tab2(i,1:nbase+nart)
                end if
            end do
            rhscol = nbase+nart+1
            do i = 1, ncon
                if (artrow(i) > 0) then
                    tab2(ncon+2,rhscol) = tab2(ncon+2,rhscol) - &
                        tab2(i,rhscol)*tab2(ncon+2,basis(i))
                end if
            end do

            do while (minval(tab2(ncon+2,1:nbase+nart)) < -ctl%zero .and. &
                result%iter1 < ctl%maxiter)
                result%iter1 = result%iter1 + 1
                call choose_pivot(tab2, ncon, ncon+2, rhscol, nbase+nart, ctl%zero, &
                    pcol, prow, ok)
                if (.not. ok) exit
                basis(prow) = pcol
                call pivot_tableau(tab2, prow, pcol)
            end do
            tab(:,1:nbase) = tab2(1:ncon+1,1:nbase)
            tab(:,nbase+1) = tab2(1:ncon+1,rhscol)
        end if

        result%iter2 = 0
        do while (minval(tab(ncon+1,1:nbase)) < -ctl%zero .and. result%iter2 < ctl%maxiter)
            result%iter2 = result%iter2 + 1
            where (abs(tab) < ctl%zero) tab = 0.0_dp
            call choose_pivot(tab, ncon, ncon+1, nbase+1, nbase, ctl%zero, pcol, prow, ok)
            if (.not. ok) exit
            basis(prow) = pcol
            call pivot_tableau(tab, prow, pcol)
        end do

        call extract_tableau_results(cvec, bvec, idir, ctl, tab, basis, result)
        if (result%iter1 >= ctl%maxiter) result%status = LINPROG_PHASE1_MAXITER
        if (result%iter2 >= ctl%maxiter) result%status = LINPROG_PHASE2_MAXITER
    end subroutine solve_with_tableau

    subroutine choose_pivot(tab, ncon, objrow, rhscol, ncols, zero, pcol, prow, ok)
        real(dp), intent(in) :: tab(:,:)
        integer, intent(in) :: ncon, objrow, rhscol, ncols
        real(dp), intent(in) :: zero
        integer, intent(out) :: pcol, prow
        logical, intent(out) :: ok
        real(dp), allocatable :: decval(:)
        real(dp) :: ratio, best_ratio, best_dec, best_obj
        integer :: i, j, best_dec_col
        logical :: have_dec, have_ratio

        allocate(decval(ncols))
        decval = linprog_nan()
        have_dec = .false.
        best_dec = huge(1.0_dp)
        best_dec_col = 0
        do j = 1, ncols
            if (tab(objrow,j) < 0.0_dp) then
                best_ratio = huge(1.0_dp)
                have_ratio = .false.
                do i = 1, ncon
                    if (tab(i,j) > 0.0_dp) then
                        ratio = tab(i,rhscol)/tab(i,j)
                        if (ieee_is_finite(ratio)) then
                            if (.not. have_ratio .or. ratio < best_ratio) best_ratio = ratio
                            have_ratio = .true.
                        end if
                    end if
                end do
                if (have_ratio) then
                    decval(j) = tab(objrow,j)*best_ratio
                    if (.not. have_dec .or. decval(j) < best_dec) then
                        best_dec = decval(j)
                        best_dec_col = j
                    end if
                    have_dec = .true.
                end if
            end if
        end do

        if (have_dec .and. best_dec < -zero) then
            pcol = best_dec_col
        else
            pcol = 0
            best_obj = huge(1.0_dp)
            do j = 1, ncols
                if (tab(objrow,j) < best_obj) then
                    best_obj = tab(objrow,j)
                    pcol = j
                end if
            end do
        end if
        if (pcol <= 0) then
            ok = .false.
            prow = 0
            return
        end if

        best_ratio = huge(1.0_dp)
        prow = 0
        do i = 1, ncon
            if (tab(i,pcol) > 0.0_dp) then
                ratio = tab(i,rhscol)/tab(i,pcol)
                if (ieee_is_finite(ratio)) then
                    if (prow == 0 .or. ratio < best_ratio) then
                        best_ratio = ratio
                        prow = i
                    end if
                end if
            end if
        end do
        ok = prow > 0
    end subroutine choose_pivot

    subroutine pivot_tableau(tab, prow, pcol)
        real(dp), intent(inout) :: tab(:,:)
        integer, intent(in) :: prow, pcol
        integer :: i
        real(dp) :: f
        tab(prow,:) = tab(prow,:)/tab(prow,pcol)
        do i = 1, size(tab,1)
            if (i /= prow) then
                f = tab(i,pcol)
                tab(i,:) = tab(i,:) - tab(prow,:)*f
            end if
        end do
    end subroutine pivot_tableau

    subroutine extract_tableau_results(cvec, bvec, idir, ctl, tab, basis, result)
        real(dp), intent(in) :: cvec(:), bvec(:), tab(:,:)
        integer, intent(in) :: idir(:), basis(:)
        type(linprog_control), intent(in) :: ctl
        type(linprog_result), intent(inout) :: result
        integer :: nvar, ncon, nbase, i, row, slack
        real(dp), allocatable :: cvec2(:), quot(:)
        real(dp) :: infp, infn, signobj, qmin, qmax, qposmin, qnegmax
        logical :: has_pos, has_neg, has_any, has_finite

        nvar = size(cvec)
        ncon = size(bvec)
        nbase = nvar+ncon
        infp = ieee_value(0.0_dp, ieee_positive_inf)
        infn = ieee_value(0.0_dp, ieee_negative_inf)
        signobj = merge(-1.0_dp, 1.0_dp, ctl%maximum)
        allocate(cvec2(nbase))
        cvec2 = 0.0_dp
        cvec2(1:nvar) = cvec

        allocate(result%basis(ncon), result%basic_values(ncon))
        result%basis = basis
        do i = 1, ncon
            result%basic_values(i) = tab(i,nbase+1)
        end do
        allocate(result%allvar_opt(nbase), result%allvar_cvec(nbase))
        allocate(result%allvar_min_c(nbase), result%allvar_max_c(nbase))
        allocate(result%allvar_marg(nbase), result%allvar_marg_reg(nbase))
        result%allvar_opt = linprog_nan()
        result%allvar_cvec = cvec2
        result%allvar_min_c = linprog_nan()
        result%allvar_max_c = linprog_nan()
        result%allvar_marg = linprog_nan()
        result%allvar_marg_reg = linprog_nan()
        allocate(quot(nbase))

        do i = 1, nbase
            row = find_basis_row(basis, i)
            if (row > 0) then
                result%allvar_opt(i) = tab(row,nbase+1)
                call safe_quotients(tab(ncon+1,1:nbase), tab(row,1:nbase), quot)
                call quotient_stats(quot, has_any, qmin, qmax, has_pos, qposmin, &
                    has_neg, qnegmax, has_finite)
                if (ctl%maximum) then
                    if (has_pos .and. qmax > 0.0_dp) result%allvar_min_c(i) = cvec2(i)-qposmin
                    if (has_finite .and. qmin < 0.0_dp) then
                        if (has_neg .and. qnegmax > -1.0e14_dp) then
                            result%allvar_max_c(i) = cvec2(i)-qnegmax
                        else
                            result%allvar_max_c(i) = infp
                        end if
                    else
                        result%allvar_max_c(i) = infp
                    end if
                else
                    if (has_pos .and. qmax > 0.0_dp) result%allvar_max_c(i) = cvec2(i)+qposmin
                    if (has_any .and. qmin < 0.0_dp) then
                        if (has_neg .and. qnegmax > -1.0e14_dp) then
                            result%allvar_min_c(i) = -cvec2(i)+qnegmax
                        end if
                    end if
                end if
            else
                result%allvar_opt(i) = 0.0_dp
                if (i <= nvar) then
                    if (ctl%maximum) then
                        result%allvar_min_c(i) = infn
                        result%allvar_max_c(i) = tab(ncon+1,i)+cvec2(i)
                    else
                        result%allvar_min_c(i) = 99.0_dp
                        result%allvar_max_c(i) = 77.0_dp
                    end if
                end if
            end if

            if (.not. (row > 0 .and. i <= nvar)) then
                result%allvar_marg(i) = tab(ncon+1,i)*signobj
                if (row == 0 .and. i > nvar) then
                    if (ctl%maximum) then
                        result%allvar_max_c(i) = tab(ncon+1,i)
                        result%allvar_min_c(i) = infn
                    else
                        result%allvar_min_c(i) = -tab(ncon+1,i)
                        result%allvar_max_c(i) = infp
                    end if
                end if
                if (row == 0) then
                    call safe_quotients(tab(1:ncon,nbase+1), tab(1:ncon,i), quot(1:ncon))
                    call min_positive(quot(1:ncon), has_pos, qposmin)
                    if (has_pos) result%allvar_marg_reg(i) = qposmin
                end if
            end if
        end do

        where (result%allvar_min_c > 1.0e16_dp) result%allvar_min_c = infp
        where (result%allvar_min_c < -1.0e16_dp) result%allvar_min_c = infn

        allocate(result%solution(nvar))
        result%solution = round_vector(result%allvar_opt(1:nvar), ctl%zero)
        allocate(result%con_actual(ncon), result%con_bvec(ncon), result%con_free(ncon))
        allocate(result%con_dual(ncon), result%con_dual_reg(ncon), result%con_dir(ncon))
        result%con_bvec = bvec
        result%con_dir = idir
        do i = 1, ncon
            slack = nvar+i
            row = find_basis_row(basis, slack)
            if (row > 0) then
                result%con_actual(i) = bvec(i) + tab(row,nbase+1)*real(idir(i),dp)
            else
                result%con_actual(i) = bvec(i)
            end if
            result%con_actual(i) = round_scalar(result%con_actual(i), ctl%zero)
            if (abs(result%allvar_opt(slack)) <= tiny(1.0_dp)) then
                result%con_dual(i) = result%allvar_marg(slack)*signobj
                result%con_dual_reg(i) = result%allvar_marg_reg(slack)
            else
                result%con_dual(i) = 0.0_dp
                result%con_dual_reg(i) = result%allvar_opt(slack)
            end if
        end do
        result%con_free = round_vector(bvec-result%con_actual, ctl%zero)
        where (idir == LINPROG_GE) result%con_free = -result%con_free
        where (idir == LINPROG_EQ) result%con_free = -abs(result%con_free)
        result%opt = round_scalar(-tab(ncon+1,nbase+1), ctl%zero)*signobj
        allocate(result%tableau(size(tab,1),size(tab,2)))
        result%tableau = tab
        result%status = LINPROG_SUCCESS
    end subroutine extract_tableau_results

    pure integer function find_basis_row(basis, col) result(row)
        integer, intent(in) :: basis(:), col
        integer :: i
        row = 0
        do i = 1, size(basis)
            if (basis(i) == col) then
                row = i
                return
            end if
        end do
    end function find_basis_row

    subroutine safe_quotients(num, den, q)
        real(dp), intent(in) :: num(:), den(:)
        real(dp), intent(out) :: q(:)
        integer :: i
        real(dp) :: infp
        infp = ieee_value(0.0_dp, ieee_positive_inf)
        do i = 1, size(q)
            if (abs(den(i)) <= tiny(1.0_dp)) then
                if (abs(num(i)) <= tiny(1.0_dp)) then
                    q(i) = linprog_nan()
                else if (num(i) > 0.0_dp) then
                    q(i) = infp
                else
                    q(i) = -infp
                end if
            else
                q(i) = num(i)/den(i)
            end if
        end do
    end subroutine safe_quotients

    subroutine quotient_stats(q, has_any, qmin, qmax, has_pos, qposmin, has_neg, &
        qnegmax, has_finite)
        real(dp), intent(in) :: q(:)
        logical, intent(out) :: has_any, has_pos, has_neg, has_finite
        real(dp), intent(out) :: qmin, qmax, qposmin, qnegmax
        integer :: i
        has_any = .false.
        has_pos = .false.
        has_neg = .false.
        has_finite = .false.
        qmin = huge(1.0_dp)
        qmax = -huge(1.0_dp)
        qposmin = huge(1.0_dp)
        qnegmax = -huge(1.0_dp)
        do i = 1, size(q)
            if (ieee_is_nan(q(i))) cycle
            has_any = .true.
            qmin = min(qmin, q(i))
            qmax = max(qmax, q(i))
            if (ieee_is_finite(q(i))) has_finite = .true.
            if (q(i) > 0.0_dp) then
                has_pos = .true.
                qposmin = min(qposmin, q(i))
            else if (q(i) < 0.0_dp) then
                has_neg = .true.
                qnegmax = max(qnegmax, q(i))
            end if
        end do
    end subroutine quotient_stats

    subroutine min_positive(q, found, value)
        real(dp), intent(in) :: q(:)
        logical, intent(out) :: found
        real(dp), intent(out) :: value
        integer :: i
        found = .false.
        value = huge(1.0_dp)
        do i = 1, size(q)
            if (.not. ieee_is_nan(q(i)) .and. q(i) > 0.0_dp) then
                value = min(value, q(i))
                found = .true.
            end if
        end do
    end subroutine min_positive

    pure real(dp) function round_scalar(x, zero) result(y)
        real(dp), intent(in) :: x, zero
        integer :: digits
        real(dp) :: scale
        if (.not. ieee_is_finite(x)) then
            y = x
            return
        end if
        digits = -nint(log10(zero))
        scale = 10.0_dp**digits
        y = real(nint(x*scale, kind=int64),dp)/scale
    end function round_scalar

    pure function round_vector(x, zero) result(y)
        real(dp), intent(in) :: x(:), zero
        real(dp) :: y(size(x))
        integer :: i
        do i = 1, size(x)
            y(i) = round_scalar(x(i), zero)
        end do
    end function round_vector

end module linprog_solver
