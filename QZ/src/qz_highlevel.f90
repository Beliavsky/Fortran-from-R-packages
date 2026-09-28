module qz_highlevel
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qz_kinds, only: dp
    use qz_types, only: qz_eigen_real_result, qz_eigen_complex_result, &
        qz_geigen_real_result, qz_geigen_complex_result, &
        qz_schur_real_result, qz_schur_complex_result, &
        qz_gschur_real_result, qz_gschur_complex_result, &
        qz_reorder_real_result, qz_reorder_complex_result, &
        qz_greorder_real_result, qz_greorder_complex_result, &
        qz_fda_geigen_result
    use qz_lapack, only: qz_dgees, qz_dgeev, qz_dgges, qz_dggev, &
        qz_dtgsen, qz_dtrsen, qz_zgees, qz_zgeev, qz_zgges, qz_zggev, &
        qz_ztgsen, qz_ztrsen
    implicit none
    private

    public :: qz_decompose, qz_geigen, geigen, ordqz
    public :: fda_geigen, h_transpose

    interface qz_decompose
        module procedure qz_real_standard
        module procedure qz_real_generalized
        module procedure qz_complex_standard
        module procedure qz_complex_generalized
    end interface qz_decompose

    interface qz_geigen
        module procedure qz_geigen_real_standard
        module procedure qz_geigen_real_generalized
        module procedure qz_geigen_complex_standard
        module procedure qz_geigen_complex_generalized
    end interface qz_geigen

    interface geigen
        module procedure qz_geigen_real_standard
        module procedure qz_geigen_real_generalized
        module procedure qz_geigen_complex_standard
        module procedure qz_geigen_complex_generalized
    end interface geigen

    interface ordqz
        module procedure ordqz_real_standard
        module procedure ordqz_real_generalized
        module procedure ordqz_complex_standard
        module procedure ordqz_complex_generalized
    end interface ordqz

    interface h_transpose
        module procedure h_transpose_real
        module procedure h_transpose_complex
    end interface h_transpose

contains

    subroutine qz_real_standard(a, result, select, only_values)
        real(dp), intent(in) :: a(:,:) !! Real square matrix to reduce to Schur form.
        type(qz_reorder_real_result), intent(out) :: result !! Schur decomposition, optional ordering, and status.
        logical, intent(in), optional :: select(:) !! Optional eigenvalue selection mask for reordering the Schur form.
        logical, intent(in), optional :: only_values !! If true and SELECT is absent, omit Schur vectors.
        type(qz_schur_real_result) :: base
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        if (present(select)) vectors = .true.
        call qz_dgees(a, base, vs=vectors)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        if (present(select)) then
            if (size(select) /= size(a, 1)) then
                result%info = -1002
                return
            end if
            call qz_dtrsen(base%t, base%q, select, result, job='B', want_q=.true.)
        else
            call copy_schur_real(base, result)
        end if
    end subroutine qz_real_standard

    subroutine qz_complex_standard(a, result, select, only_values)
        complex(dp), intent(in) :: a(:,:) !! Complex square matrix to reduce to Schur form.
        type(qz_reorder_complex_result), intent(out) :: result !! Schur decomposition, optional ordering, and status.
        logical, intent(in), optional :: select(:) !! Optional eigenvalue selection mask for reordering the Schur form.
        logical, intent(in), optional :: only_values !! If true and SELECT is absent, omit Schur vectors.
        type(qz_schur_complex_result) :: base
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        if (present(select)) vectors = .true.
        call qz_zgees(a, base, vs=vectors)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        if (present(select)) then
            if (size(select) /= size(a, 1)) then
                result%info = -1002
                return
            end if
            call qz_ztrsen(base%t, base%q, select, result, job='B', want_q=.true.)
        else
            call copy_schur_complex(base, result)
        end if
    end subroutine qz_complex_standard

    subroutine qz_real_generalized(a, b, result, select, only_values)
        real(dp), intent(in) :: a(:,:) !! First real matrix in the generalized Schur decomposition.
        real(dp), intent(in) :: b(:,:) !! Second real matrix in the generalized Schur decomposition.
        type(qz_greorder_real_result), intent(out) :: result !! Generalized Schur decomposition, optional ordering, and status.
        logical, intent(in), optional :: select(:) !! Optional generalized-eigenvalue selection mask for reordering.
        logical, intent(in), optional :: only_values !! If true and SELECT is absent, omit left/right Schur vectors.
        type(qz_gschur_real_result) :: base
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        if (present(select)) vectors = .true.
        call qz_dgges(a, b, base, vsl=vectors, vsr=vectors)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        if (present(select)) then
            if (size(select) /= size(a, 1)) then
                result%info = -1002
                return
            end if
            call qz_dtgsen(base%s, base%t, base%q, base%z, select, result, &
                ijob=4, want_q=.true., want_z=.true.)
        else
            call copy_gschur_real(base, result)
        end if
    end subroutine qz_real_generalized

    subroutine qz_complex_generalized(a, b, result, select, only_values)
        complex(dp), intent(in) :: a(:,:) !! First complex matrix in the generalized Schur decomposition.
        complex(dp), intent(in) :: b(:,:) !! Second complex matrix in the generalized Schur decomposition.
        type(qz_greorder_complex_result), intent(out) :: result !! Generalized Schur decomposition, optional ordering, and status.
        logical, intent(in), optional :: select(:) !! Optional generalized-eigenvalue selection mask for reordering.
        logical, intent(in), optional :: only_values !! If true and SELECT is absent, omit left/right Schur vectors.
        type(qz_gschur_complex_result) :: base
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        if (present(select)) vectors = .true.
        call qz_zgges(a, b, base, vsl=vectors, vsr=vectors)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        if (present(select)) then
            if (size(select) /= size(a, 1)) then
                result%info = -1002
                return
            end if
            call qz_ztgsen(base%s, base%t, base%q, base%z, select, result, &
                ijob=4, want_q=.true., want_z=.true.)
        else
            call copy_gschur_complex(base, result)
        end if
    end subroutine qz_complex_generalized

    subroutine qz_geigen_real_standard(a, result, only_values)
        real(dp), intent(in) :: a(:,:) !! Real square matrix whose eigenvalues or eigensystem are required.
        type(qz_eigen_real_result), intent(out) :: result !! Standard eigenvalues, optional left/right eigenvectors, and status.
        logical, intent(in), optional :: only_values !! If true, compute eigenvalues without eigenvectors.
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        call qz_dgeev(a, result, vl=vectors, vr=vectors)
    end subroutine qz_geigen_real_standard

    subroutine qz_geigen_complex_standard(a, result, only_values)
        complex(dp), intent(in) :: a(:,:) !! Complex square matrix whose eigenvalues or eigensystem are required.
        type(qz_eigen_complex_result), intent(out) :: result !! Standard eigenvalues, optional left/right eigenvectors, and status.
        logical, intent(in), optional :: only_values !! If true, compute eigenvalues without eigenvectors.
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        call qz_zgeev(a, result, vl=vectors, vr=vectors)
    end subroutine qz_geigen_complex_standard

    subroutine qz_geigen_real_generalized(a, b, result, only_values)
        real(dp), intent(in) :: a(:,:) !! First real matrix of the generalized eigenproblem.
        real(dp), intent(in) :: b(:,:) !! Second real matrix of the generalized eigenproblem.
        type(qz_geigen_real_result), intent(out) :: result !! Generalized eigenvalue factors, optional vectors, and status.
        logical, intent(in), optional :: only_values !! If true, compute generalized eigenvalues without eigenvectors.
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        call qz_dggev(a, b, result, vl=vectors, vr=vectors)
    end subroutine qz_geigen_real_generalized

    subroutine qz_geigen_complex_generalized(a, b, result, only_values)
        complex(dp), intent(in) :: a(:,:) !! First complex matrix of the generalized eigenproblem.
        complex(dp), intent(in) :: b(:,:) !! Second complex matrix of the generalized eigenproblem.
        type(qz_geigen_complex_result), intent(out) :: result !! Generalized eigenvalue factors, optional vectors, and status.
        logical, intent(in), optional :: only_values !! If true, compute generalized eigenvalues without eigenvectors.
        logical :: vectors

        vectors = .true.
        if (present(only_values)) vectors = .not. only_values
        call qz_zggev(a, b, result, vl=vectors, vr=vectors)
    end subroutine qz_geigen_complex_generalized

    subroutine ordqz_real_standard(a, result, keyword, cluster)
        real(dp), intent(in) :: a(:,:) !! Real square matrix whose Schur form is to be ordered.
        type(qz_reorder_real_result), intent(out) :: result !! Ordered real Schur decomposition and condition estimates.
        character(len=*), intent(in), optional :: keyword !! Ordering rule: lhp, rhp, udi, udo, ref, cef, or finite-only variants.
        logical, intent(in), optional :: cluster(:) !! Explicit selection mask; when present it overrides KEYWORD.
        type(qz_schur_real_result) :: base
        logical, allocatable :: select(:)

        call qz_dgees(a, base, vs=.true.)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        allocate(select(size(a, 1)))
        if (present(cluster)) then
            if (size(cluster) /= size(select)) then
                result%info = -1002
                return
            end if
            select = cluster
        else
            call select_by_keyword(base%w, keyword_or_default(keyword), select)
        end if
        call qz_dtrsen(base%t, base%q, select, result, job='B', want_q=.true.)
    end subroutine ordqz_real_standard

    subroutine ordqz_complex_standard(a, result, keyword, cluster)
        complex(dp), intent(in) :: a(:,:) !! Complex square matrix whose Schur form is to be ordered.
        type(qz_reorder_complex_result), intent(out) :: result !! Ordered complex Schur decomposition and condition estimates.
        character(len=*), intent(in), optional :: keyword !! Ordering rule: lhp, rhp, udi, udo, ref, cef, or finite-only variants.
        logical, intent(in), optional :: cluster(:) !! Explicit selection mask; when present it overrides KEYWORD.
        type(qz_schur_complex_result) :: base
        logical, allocatable :: select(:)

        call qz_zgees(a, base, vs=.true.)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        allocate(select(size(a, 1)))
        if (present(cluster)) then
            if (size(cluster) /= size(select)) then
                result%info = -1002
                return
            end if
            select = cluster
        else
            call select_by_keyword(base%w, keyword_or_default(keyword), select)
        end if
        call qz_ztrsen(base%t, base%q, select, result, job='B', want_q=.true.)
    end subroutine ordqz_complex_standard

    subroutine ordqz_real_generalized(a, b, result, keyword, cluster)
        real(dp), intent(in) :: a(:,:) !! First real matrix whose generalized Schur form is to be ordered.
        real(dp), intent(in) :: b(:,:) !! Second real matrix whose generalized Schur form is to be ordered.
        type(qz_greorder_real_result), intent(out) :: result !! Ordered generalized real Schur decomposition and estimates.
        character(len=*), intent(in), optional :: keyword !! Ordering rule: lhp, rhp, udi, udo, ref, cef, or finite-only variants.
        logical, intent(in), optional :: cluster(:) !! Explicit selection mask; when present it overrides KEYWORD.
        type(qz_gschur_real_result) :: base
        logical, allocatable :: select(:)

        call qz_dgges(a, b, base, vsl=.true., vsr=.true.)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        allocate(select(size(a, 1)))
        if (present(cluster)) then
            if (size(cluster) /= size(select)) then
                result%info = -1002
                return
            end if
            select = cluster
        else
            call select_by_keyword(base%lambda, keyword_or_default(keyword), select)
        end if
        call qz_dtgsen(base%s, base%t, base%q, base%z, select, result, &
            ijob=4, want_q=.true., want_z=.true.)
    end subroutine ordqz_real_generalized

    subroutine ordqz_complex_generalized(a, b, result, keyword, cluster)
        complex(dp), intent(in) :: a(:,:) !! First complex matrix whose generalized Schur form is to be ordered.
        complex(dp), intent(in) :: b(:,:) !! Second complex matrix whose generalized Schur form is to be ordered.
        type(qz_greorder_complex_result), intent(out) :: result !! Ordered generalized complex Schur decomposition and estimates.
        character(len=*), intent(in), optional :: keyword !! Ordering rule: lhp, rhp, udi, udo, ref, cef, or finite-only variants.
        logical, intent(in), optional :: cluster(:) !! Explicit selection mask; when present it overrides KEYWORD.
        type(qz_gschur_complex_result) :: base
        logical, allocatable :: select(:)

        call qz_zgges(a, b, base, vsl=.true., vsr=.true.)
        if (base%info /= 0) then
            result%info = base%info
            return
        end if
        allocate(select(size(a, 1)))
        if (present(cluster)) then
            if (size(cluster) /= size(select)) then
                result%info = -1002
                return
            end if
            select = cluster
        else
            call select_by_keyword(base%lambda, keyword_or_default(keyword), select)
        end if
        call qz_ztgsen(base%s, base%t, base%q, base%z, select, result, &
            ijob=4, want_q=.true., want_z=.true.)
    end subroutine ordqz_complex_generalized

    subroutine fda_geigen(amat, bmat, cmat, result)
        real(dp), intent(in) :: amat(:,:) !! Cross-block matrix of shape p by q.
        real(dp), intent(in) :: bmat(:,:) !! Symmetric positive-definite p by p weighting matrix.
        real(dp), intent(in) :: cmat(:,:) !! Symmetric positive-definite q by q weighting matrix.
        type(qz_fda_geigen_result), intent(out) :: result !! Generalized eigenvalues and normalized L/M coefficient matrices.
        type(qz_geigen_real_result) :: eig
        real(dp), allocatable :: aaug(:,:), baug(:,:), bsym(:,:), csym(:,:)
        integer, allocatable :: order(:)
        integer :: p, qn, sdim, i, j, best, tmp
        real(dp) :: scale_l, scale_m

        p = size(amat, 1)
        qn = size(amat, 2)
        if (size(bmat, 1) /= p .or. size(bmat, 2) /= p .or. &
            size(cmat, 1) /= qn .or. size(cmat, 2) /= qn) then
            result%info = -1001
            return
        end if
        sdim = min(p, qn)
        allocate(bsym(p, p), csym(qn, qn), aaug(p + qn, p + qn), baug(p + qn, p + qn))
        bsym = 0.5_dp * (bmat + transpose(bmat))
        csym = 0.5_dp * (cmat + transpose(cmat))
        aaug = 0.0_dp
        baug = 0.0_dp
        aaug(1:p, p + 1:p + qn) = amat
        aaug(p + 1:p + qn, 1:p) = transpose(amat)
        baug(1:p, 1:p) = bsym
        baug(p + 1:p + qn, p + 1:p + qn) = csym
        call qz_dggev(aaug, baug, eig, vl=.true., vr=.true.)
        result%info = eig%info
        if (result%info /= 0) return

        allocate(order(p + qn))
        order = [(i, i=1, p + qn)]
        do i = 1, p + qn - 1
            best = i
            do j = i + 1, p + qn
                if (real(eig%lambda(order(j)), dp) > real(eig%lambda(order(best)), dp)) best = j
            end do
            if (best /= i) then
                tmp = order(i)
                order(i) = order(best)
                order(best) = tmp
            end if
        end do

        allocate(result%values(sdim), result%lmat(p, sdim), result%mmat(qn, sdim))
        do i = 1, sdim
            result%values(i) = real(eig%lambda(order(i)), dp)
            result%lmat(:, i) = real(eig%u(1:p, order(i)), dp)
            result%mmat(:, i) = real(eig%u(p + 1:p + qn, order(i)), dp)
            scale_l = sqrt(max(0.0_dp, dot_product(result%lmat(:, i), matmul(bsym, result%lmat(:, i)))))
            scale_m = sqrt(max(0.0_dp, dot_product(result%mmat(:, i), matmul(csym, result%mmat(:, i)))))
            if (scale_l > tiny(1.0_dp)) result%lmat(:, i) = result%lmat(:, i) / scale_l
            if (scale_m > tiny(1.0_dp)) result%mmat(:, i) = result%mmat(:, i) / scale_m
        end do
    end subroutine fda_geigen

    pure function h_transpose_real(x) result(y)
        real(dp), intent(in) :: x(:,:) !! Real matrix to transpose; conjugation is a no-op for real values.
        real(dp) :: y(size(x, 2), size(x, 1))

        y = transpose(x)
    end function h_transpose_real

    pure function h_transpose_complex(x) result(y)
        complex(dp), intent(in) :: x(:,:) !! Complex matrix to conjugate-transpose.
        complex(dp) :: y(size(x, 2), size(x, 1))

        y = conjg(transpose(x))
    end function h_transpose_complex

    subroutine select_by_keyword(values, keyword, select)
        complex(dp), intent(in) :: values(:) !! Standard or generalized eigenvalues on which the ordering rule is evaluated.
        character(len=*), intent(in) :: keyword !! Normalized ordering keyword.
        logical, intent(out) :: select(:) !! Logical mask selecting eigenvalues to move to the leading Schur block.
        character(len=:), allocatable :: key
        integer :: i
        real(dp) :: tol
        logical :: finite_value, real_value

        key = lower_ascii(trim(keyword))
        select = .false.
        do i = 1, size(values)
            finite_value = ieee_is_finite(real(values(i), dp)) .and. ieee_is_finite(aimag(values(i)))
            tol = 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(values(i)))
            real_value = abs(aimag(values(i))) <= tol
            select case (key)
            case ('lhp')
                select(i) = real(values(i), dp) < 0.0_dp
            case ('rhp')
                select(i) = real(values(i), dp) >= 0.0_dp
            case ('udi')
                select(i) = abs(values(i)) < 1.0_dp
            case ('udo')
                select(i) = abs(values(i)) >= 1.0_dp
            case ('lhp.fo')
                select(i) = real(values(i), dp) < 0.0_dp .and. finite_value
            case ('rhp.fo')
                select(i) = real(values(i), dp) >= 0.0_dp .and. finite_value
            case ('udi.fo')
                select(i) = abs(values(i)) < 1.0_dp .and. finite_value
            case ('udo.fo')
                select(i) = abs(values(i)) >= 1.0_dp .and. finite_value
            case ('ref')
                select(i) = real_value
            case ('cef')
                select(i) = .not. real_value
            case default
                select(i) = real(values(i), dp) < 0.0_dp
            end select
        end do
    end subroutine select_by_keyword

    function keyword_or_default(keyword) result(value)
        character(len=*), intent(in), optional :: keyword !! Optional ordering keyword; absent means left-half-plane ordering.
        character(len=:), allocatable :: value

        if (present(keyword)) then
            value = trim(keyword)
        else
            value = 'lhp'
        end if
    end function keyword_or_default

    pure function lower_ascii(text) result(lower)
        character(len=*), intent(in) :: text !! ASCII text to convert to lower case.
        character(len=len(text)) :: lower
        integer :: i, code

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
        end do
    end function lower_ascii

    subroutine copy_schur_real(base, result)
        type(qz_schur_real_result), intent(in) :: base !! Unordered real Schur decomposition returned by QZ_DGEES.
        type(qz_reorder_real_result), intent(out) :: result !! High-level result populated from the unordered decomposition.

        result%info = base%info
        result%work_opt = base%work_opt
        allocate(result%t, source=base%t)
        allocate(result%wr, source=base%wr)
        allocate(result%wi, source=base%wi)
        allocate(result%w, source=base%w)
        if (allocated(base%q)) allocate(result%q, source=base%q)
    end subroutine copy_schur_real

    subroutine copy_schur_complex(base, result)
        type(qz_schur_complex_result), intent(in) :: base !! Unordered complex Schur decomposition returned by QZ_ZGEES.
        type(qz_reorder_complex_result), intent(out) :: result !! High-level result populated from the unordered decomposition.

        result%info = base%info
        result%work_opt = base%work_opt
        allocate(result%t, source=base%t)
        allocate(result%w, source=base%w)
        if (allocated(base%q)) allocate(result%q, source=base%q)
    end subroutine copy_schur_complex

    subroutine copy_gschur_real(base, result)
        type(qz_gschur_real_result), intent(in) :: base !! Unordered real generalized Schur decomposition from QZ_DGGES.
        type(qz_greorder_real_result), intent(out) :: result !! High-level result populated from the unordered decomposition.

        result%info = base%info
        result%work_opt = base%work_opt
        allocate(result%s, source=base%s)
        allocate(result%t, source=base%t)
        allocate(result%alphar, source=base%alphar)
        allocate(result%alphai, source=base%alphai)
        allocate(result%beta, source=base%beta)
        allocate(result%alpha, source=base%alpha)
        allocate(result%lambda, source=base%lambda)
        if (allocated(base%q)) allocate(result%q, source=base%q)
        if (allocated(base%z)) allocate(result%z, source=base%z)
    end subroutine copy_gschur_real

    subroutine copy_gschur_complex(base, result)
        type(qz_gschur_complex_result), intent(in) :: base !! Unordered complex generalized Schur decomposition from QZ_ZGGES.
        type(qz_greorder_complex_result), intent(out) :: result !! High-level result populated from the unordered decomposition.

        result%info = base%info
        result%work_opt = base%work_opt
        allocate(result%s, source=base%s)
        allocate(result%t, source=base%t)
        allocate(result%alpha, source=base%alpha)
        allocate(result%beta, source=base%beta)
        allocate(result%lambda, source=base%lambda)
        if (allocated(base%q)) allocate(result%q, source=base%q)
        if (allocated(base%z)) allocate(result%z, source=base%z)
    end subroutine copy_gschur_complex
end module qz_highlevel
