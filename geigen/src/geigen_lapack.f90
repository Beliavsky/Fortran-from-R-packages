module geigen_lapack
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_positive_inf, ieee_value
    use geigen_kinds, only: dp
    use geigen_types, only: geigen_result, gqz_result, gsvd_result
    use la_lapack, only: ggev, gges, sygv, hegv, gesvd, syev, heev, gerqf, ormrq, unmrq
    implicit none
    private

    public :: geigen_real, geigen_complex
    public :: gqz_real, gqz_complex
    public :: gsvd_real, gsvd_complex

contains

    subroutine geigen_real(a, b, result, symmetric, only_values)
        real(dp), intent(in) :: a(:,:) !! First real square matrix in A*x = lambda*B*x.
        real(dp), intent(in) :: b(:,:) !! Second real square matrix, with the same order as A.
        type(geigen_result), intent(out) :: result !! Generalized eigenvalues, optional right eigenvectors, and status.
        logical, intent(in), optional :: symmetric !! Use the symmetric-definite solver; defaults to automatic symmetry detection.
        logical, intent(in), optional :: only_values !! If true, omit right eigenvectors; default is false.
        logical :: use_symmetric, values_only

        use_symmetric = is_symmetric_real(a) .and. is_symmetric_real(b)
        if (present(symmetric)) use_symmetric = symmetric
        values_only = .false.
        if (present(only_values)) values_only = only_values

        if (use_symmetric) then
            call geigen_symmetric_real(a, b, result, values_only)
        else
            call geigen_general_real(a, b, result, values_only)
        end if
    end subroutine geigen_real

    subroutine geigen_complex(a, b, result, symmetric, only_values)
        complex(dp), intent(in) :: a(:,:) !! First complex square matrix in A*x = lambda*B*x.
        complex(dp), intent(in) :: b(:,:) !! Second complex square matrix, with the same order as A.
        type(geigen_result), intent(out) :: result !! Generalized eigenvalues, optional right eigenvectors, and status.
        logical, intent(in), optional :: symmetric !! Use the Hermitian-definite solver; defaults to automatic Hermitian detection.
        logical, intent(in), optional :: only_values !! If true, omit right eigenvectors; default is false.
        logical :: use_symmetric, values_only

        use_symmetric = is_hermitian_complex(a) .and. is_hermitian_complex(b)
        if (present(symmetric)) use_symmetric = symmetric
        values_only = .false.
        if (present(only_values)) values_only = only_values

        if (use_symmetric) then
            call geigen_hermitian_complex(a, b, result, values_only)
        else
            call geigen_general_complex(a, b, result, values_only)
        end if
    end subroutine geigen_complex

    subroutine geigen_symmetric_real(a, b, result, only_values)
        real(dp), intent(in) :: a(:,:) !! Real symmetric matrix A.
        real(dp), intent(in) :: b(:,:) !! Real symmetric positive-definite matrix B.
        type(geigen_result), intent(out) :: result !! Real generalized eigenpairs represented in complex storage.
        logical, intent(in) :: only_values !! If true, compute eigenvalues only.
        integer :: n, lwork
        character :: jobz
        real(dp), allocatable :: acopy(:,:), bcopy(:,:), w(:), work(:)

        n = size(a, 1)
        if (.not. valid_square_pair_real(a, b)) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_real(a) .or. .not. all_finite_real(b)) then
            result%info = -1002
            return
        end if

        jobz = merge('N', 'V', only_values)
        lwork = max(1, 8 * n + 32)
        allocate(acopy(n, n), bcopy(n, n), w(n), work(lwork))
        acopy = a
        bcopy = b
        call sygv(1, jobz, 'L', n, acopy, n, bcopy, n, w, work, lwork, result%info)
        if (result%info /= 0) return

        allocate(result%values(n))
        result%values = cmplx(w, 0.0_dp, kind=dp)
        if (.not. only_values) then
            allocate(result%vectors(n, n))
            result%vectors = cmplx(acopy, 0.0_dp, kind=dp)
        end if
        result%symmetric_used = .true.
    end subroutine geigen_symmetric_real

    subroutine geigen_hermitian_complex(a, b, result, only_values)
        complex(dp), intent(in) :: a(:,:) !! Complex Hermitian matrix A.
        complex(dp), intent(in) :: b(:,:) !! Complex Hermitian positive-definite matrix B.
        type(geigen_result), intent(out) :: result !! Hermitian generalized eigenpairs and status.
        logical, intent(in) :: only_values !! If true, compute eigenvalues only.
        integer :: n, lwork
        character :: jobz
        complex(dp), allocatable :: acopy(:,:), bcopy(:,:), work(:)
        real(dp), allocatable :: w(:), rwork(:)

        n = size(a, 1)
        if (.not. valid_square_pair_complex(a, b)) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_complex(a) .or. .not. all_finite_complex(b)) then
            result%info = -1002
            return
        end if

        jobz = merge('N', 'V', only_values)
        lwork = max(1, 8 * n + 32)
        allocate(acopy(n, n), bcopy(n, n), w(n), work(lwork), rwork(max(1, 3 * n - 2)))
        acopy = a
        bcopy = b
        call hegv(1, jobz, 'L', n, acopy, n, bcopy, n, w, work, lwork, rwork, result%info)
        if (result%info /= 0) return

        allocate(result%values(n))
        result%values = cmplx(w, 0.0_dp, kind=dp)
        if (.not. only_values) then
            allocate(result%vectors(n, n))
            result%vectors = acopy
        end if
        result%symmetric_used = .true.
    end subroutine geigen_hermitian_complex

    subroutine geigen_general_real(a, b, result, only_values)
        real(dp), intent(in) :: a(:,:) !! First real square matrix A.
        real(dp), intent(in) :: b(:,:) !! Second real square matrix B.
        type(geigen_result), intent(out) :: result !! Generalized eigenpairs, alpha/beta factors, and status.
        logical, intent(in) :: only_values !! If true, compute eigenvalues only.
        integer :: n, lwork, ldvr
        character :: jobvr
        real(dp), allocatable :: acopy(:,:), bcopy(:,:), alphar(:), alphai(:), beta(:)
        real(dp), allocatable :: vl(:,:), vr(:,:), work(:)
        complex(dp), allocatable :: alpha_c(:), beta_c(:)

        n = size(a, 1)
        if (.not. valid_square_pair_real(a, b)) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_real(a) .or. .not. all_finite_real(b)) then
            result%info = -1002
            return
        end if

        jobvr = merge('N', 'V', only_values)
        ldvr = merge(1, n, only_values)
        lwork = max(1, 8 * n + 32)
        allocate(acopy(n, n), bcopy(n, n), alphar(n), alphai(n), beta(n), work(lwork))
        allocate(vl(1, max(1, n)), vr(ldvr, max(1, n)))
        acopy = a
        bcopy = b
        vl = 0.0_dp
        vr = 0.0_dp
        call ggev('N', jobvr, n, acopy, n, bcopy, n, alphar, alphai, beta, &
            vl, 1, vr, ldvr, work, lwork, result%info)
        if (result%info /= 0) return

        allocate(alpha_c(n), beta_c(n), result%alpha(n), result%beta(n), result%values(n))
        alpha_c = cmplx(alphar, alphai, kind=dp)
        beta_c = cmplx(beta, 0.0_dp, kind=dp)
        result%alpha = alpha_c
        result%beta = beta_c
        call complex_ratio(alpha_c, beta_c, result%values)
        if (.not. only_values) then
            call unpack_real_eigenvectors(vr, alphai, result%vectors)
        end if
        result%symmetric_used = .false.
    end subroutine geigen_general_real

    subroutine geigen_general_complex(a, b, result, only_values)
        complex(dp), intent(in) :: a(:,:) !! First complex square matrix A.
        complex(dp), intent(in) :: b(:,:) !! Second complex square matrix B.
        type(geigen_result), intent(out) :: result !! Generalized eigenpairs, alpha/beta factors, and status.
        logical, intent(in) :: only_values !! If true, compute eigenvalues only.
        integer :: n, lwork, ldvr
        character :: jobvr
        complex(dp), allocatable :: acopy(:,:), bcopy(:,:), vl(:,:), vr(:,:), work(:)
        real(dp), allocatable :: rwork(:)

        n = size(a, 1)
        if (.not. valid_square_pair_complex(a, b)) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_complex(a) .or. .not. all_finite_complex(b)) then
            result%info = -1002
            return
        end if

        jobvr = merge('N', 'V', only_values)
        ldvr = merge(1, n, only_values)
        lwork = max(1, 8 * n + 32)
        allocate(acopy(n, n), bcopy(n, n), result%alpha(n), result%beta(n), result%values(n))
        allocate(vl(1, max(1, n)), vr(ldvr, max(1, n)), work(lwork), rwork(max(1, 8 * n)))
        acopy = a
        bcopy = b
        vl = cmplx(0.0_dp, 0.0_dp, kind=dp)
        vr = cmplx(0.0_dp, 0.0_dp, kind=dp)
        call ggev('N', jobvr, n, acopy, n, bcopy, n, result%alpha, result%beta, &
            vl, 1, vr, ldvr, work, lwork, rwork, result%info)
        if (result%info /= 0) return

        call complex_ratio(result%alpha, result%beta, result%values)
        if (.not. only_values) then
            allocate(result%vectors(n, n))
            result%vectors = vr(1:n, 1:n)
        end if
        result%symmetric_used = .false.
    end subroutine geigen_general_complex

    subroutine gqz_real(a, b, result, sort)
        real(dp), intent(in) :: a(:,:) !! First real square matrix in the generalized Schur decomposition.
        real(dp), intent(in) :: b(:,:) !! Second real square matrix, with the same order as A.
        type(gqz_result), intent(out) :: result !! Generalized Schur factors, vectors, eigenvalues, selected dimension, and status.
        character(len=*), intent(in), optional :: sort !! Ordering code N, -, +, S, B, or R; default is N.
        character :: sort_use, sort_job
        integer :: n, lwork
        real(dp), allocatable :: s(:,:), t(:,:), alphar(:), alphai(:), beta(:), q(:,:), z(:,:), work(:)
        logical, allocatable :: bwork(:)

        n = size(a, 1)
        if (.not. valid_square_pair_real(a, b)) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_real(a) .or. .not. all_finite_real(b)) then
            result%info = -1002
            return
        end if

        sort_use = 'N'
        if (present(sort)) then
            if (len_trim(sort) > 0) sort_use = sort(1:1)
        end if
        if (index('N-+SBR', sort_use) == 0) then
            result%info = -1003
            return
        end if
        sort_job = merge('N', 'S', sort_use == 'N')
        lwork = max(1, 8 * n + 32)
        allocate(s(n, n), t(n, n), alphar(n), alphai(n), beta(n), q(n, n), z(n, n), work(lwork))
        allocate(bwork(max(1, n)))
        s = a
        t = b
        q = 0.0_dp
        z = 0.0_dp
        bwork = .false.
        select case (sort_use)
        case ('N')
            call gges('V', 'V', sort_job, select_none_real, n, s, n, t, n, result%sdim, &
                alphar, alphai, beta, q, n, z, n, work, lwork, bwork, result%info)
        case ('-')
            call gges('V', 'V', sort_job, select_negative_real, n, s, n, t, n, result%sdim, &
                alphar, alphai, beta, q, n, z, n, work, lwork, bwork, result%info)
        case ('+')
            call gges('V', 'V', sort_job, select_positive_real, n, s, n, t, n, result%sdim, &
                alphar, alphai, beta, q, n, z, n, work, lwork, bwork, result%info)
        case ('S')
            call gges('V', 'V', sort_job, select_inside_real, n, s, n, t, n, result%sdim, &
                alphar, alphai, beta, q, n, z, n, work, lwork, bwork, result%info)
        case ('B')
            call gges('V', 'V', sort_job, select_outside_real, n, s, n, t, n, result%sdim, &
                alphar, alphai, beta, q, n, z, n, work, lwork, bwork, result%info)
        case ('R')
            call gges('V', 'V', sort_job, select_real_real, n, s, n, t, n, result%sdim, &
                alphar, alphai, beta, q, n, z, n, work, lwork, bwork, result%info)
        end select
        if (result%info /= 0) return

        allocate(result%s(n, n), result%t(n, n), result%q(n, n), result%z(n, n))
        allocate(result%alpha(n), result%beta(n), result%values(n))
        result%s = cmplx(s, 0.0_dp, kind=dp)
        result%t = cmplx(t, 0.0_dp, kind=dp)
        result%q = cmplx(q, 0.0_dp, kind=dp)
        result%z = cmplx(z, 0.0_dp, kind=dp)
        result%alpha = cmplx(alphar, alphai, kind=dp)
        result%beta = cmplx(beta, 0.0_dp, kind=dp)
        call complex_ratio(result%alpha, result%beta, result%values)
    end subroutine gqz_real

    subroutine gqz_complex(a, b, result, sort)
        complex(dp), intent(in) :: a(:,:) !! First complex square matrix in the generalized Schur decomposition.
        complex(dp), intent(in) :: b(:,:) !! Second complex square matrix, with the same order as A.
        type(gqz_result), intent(out) :: result !! Generalized Schur factors, vectors, eigenvalues, selected dimension, and status.
        character(len=*), intent(in), optional :: sort !! Ordering code N, -, +, S, B, or R; default is N.
        character :: sort_use, sort_job
        integer :: n, lwork
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)
        logical, allocatable :: bwork(:)

        n = size(a, 1)
        if (.not. valid_square_pair_complex(a, b)) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_complex(a) .or. .not. all_finite_complex(b)) then
            result%info = -1002
            return
        end if

        sort_use = 'N'
        if (present(sort)) then
            if (len_trim(sort) > 0) sort_use = sort(1:1)
        end if
        if (index('N-+SBR', sort_use) == 0) then
            result%info = -1003
            return
        end if
        sort_job = merge('N', 'S', sort_use == 'N')
        lwork = max(1, 8 * n + 32)
        allocate(result%s(n, n), result%t(n, n), result%q(n, n), result%z(n, n))
        allocate(result%alpha(n), result%beta(n), result%values(n), work(lwork))
        allocate(rwork(max(1, 8 * n)), bwork(max(1, n)))
        result%s = a
        result%t = b
        result%q = cmplx(0.0_dp, 0.0_dp, kind=dp)
        result%z = cmplx(0.0_dp, 0.0_dp, kind=dp)
        bwork = .false.
        select case (sort_use)
        case ('N')
            call gges('V', 'V', sort_job, select_none_complex, n, result%s, n, result%t, n, &
                result%sdim, result%alpha, result%beta, result%q, n, result%z, n, work, lwork, &
                rwork, bwork, result%info)
        case ('-')
            call gges('V', 'V', sort_job, select_negative_complex, n, result%s, n, result%t, n, &
                result%sdim, result%alpha, result%beta, result%q, n, result%z, n, work, lwork, &
                rwork, bwork, result%info)
        case ('+')
            call gges('V', 'V', sort_job, select_positive_complex, n, result%s, n, result%t, n, &
                result%sdim, result%alpha, result%beta, result%q, n, result%z, n, work, lwork, &
                rwork, bwork, result%info)
        case ('S')
            call gges('V', 'V', sort_job, select_inside_complex, n, result%s, n, result%t, n, &
                result%sdim, result%alpha, result%beta, result%q, n, result%z, n, work, lwork, &
                rwork, bwork, result%info)
        case ('B')
            call gges('V', 'V', sort_job, select_outside_complex, n, result%s, n, result%t, n, &
                result%sdim, result%alpha, result%beta, result%q, n, result%z, n, work, lwork, &
                rwork, bwork, result%info)
        case ('R')
            call gges('V', 'V', sort_job, select_real_complex, n, result%s, n, result%t, n, &
                result%sdim, result%alpha, result%beta, result%q, n, result%z, n, work, lwork, &
                rwork, bwork, result%info)
        end select
        if (result%info /= 0) return
        call complex_ratio(result%alpha, result%beta, result%values)
    end subroutine gqz_complex

    subroutine gsvd_real(a, b, result)
        real(dp), intent(in) :: a(:,:) !! Real M-by-N matrix in the matrix pair.
        real(dp), intent(in) :: b(:,:) !! Real P-by-N matrix sharing the column count of A.
        type(gsvd_result), intent(out) :: result !! Generalized SVD factors, common triangular factor, and status.
        integer :: m, n, p, rows, mn, r, j, c, lwork, info_eig, info_rq
        real(dp) :: rank_tol, ab_tol, norm_ab
        real(dp), allocatable :: stack(:,:), sval(:), us(:,:), vt(:,:), work(:)
        real(dp), allocatable :: w1(:,:), w2(:,:), gram(:,:), eval(:), zrot(:,:), temp(:,:), xinv(:,:)
        real(dp), allocatable :: rq(:,:), tau(:), qreal(:,:), ureal(:,:), vreal(:,:)
        real(dp), allocatable :: alpha_r(:), beta_r(:), eigwork(:), rqwork(:)

        m = size(a, 1)
        n = size(a, 2)
        p = size(b, 1)
        result%m = m
        result%n = n
        result%p = p
        if (m < 1 .or. p < 1 .or. n < 1 .or. size(b, 2) /= n) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_real(a) .or. .not. all_finite_real(b)) then
            result%info = -1002
            return
        end if

        rows = m + p
        mn = min(rows, n)
        lwork = max(1, 12 * max(rows, n) + 8 * mn + 64)
        allocate(stack(rows, n), sval(mn), us(rows, mn), vt(mn, n), work(lwork))
        stack(1:m, :) = a
        stack(m + 1:rows, :) = b
        call gesvd('S', 'S', rows, n, stack, rows, sval, us, rows, vt, mn, work, lwork, result%info)
        if (result%info /= 0) return

        rank_tol = real(max(rows, n), dp) * epsilon(1.0_dp)
        if (mn > 0) rank_tol = rank_tol * max(1.0_dp, sval(1))
        r = count(sval > rank_tol)
        result%rank = r
        if (r == 0) then
            result%k = 0
            result%l = 0
            allocate(result%alpha(n), result%beta(n), result%u(m, m), result%v(p, p), result%q(n, n), result%r(0, 0))
            result%alpha = 0.0_dp
            result%beta = 0.0_dp
            call identity_complex(result%u)
            call identity_complex(result%v)
            call identity_complex(result%q)
            return
        end if

        allocate(w1(m, r), w2(p, r), gram(r, r), eval(r), zrot(r, r), alpha_r(r), beta_r(r))
        w1 = us(1:m, 1:r)
        w2 = us(m + 1:rows, 1:r)
        gram = matmul(transpose(w1), w1)
        allocate(eigwork(max(1, 8 * r + 32)))
        call syev('V', 'U', r, gram, r, eval, eigwork, size(eigwork), info_eig)
        if (info_eig /= 0) then
            result%info = 2000 + info_eig
            return
        end if
        ab_tol = 1000.0_dp * real(max(m, max(p, n)), dp) * epsilon(1.0_dp)
        do j = 1, r
            zrot(:, j) = gram(:, r - j + 1)
            if (1.0_dp - min(1.0_dp, eval(r - j + 1)) <= ab_tol) then
                alpha_r(j) = 1.0_dp
                beta_r(j) = 0.0_dp
            else if (max(0.0_dp, eval(r - j + 1)) <= ab_tol) then
                alpha_r(j) = 0.0_dp
                beta_r(j) = 1.0_dp
            else
                alpha_r(j) = sqrt(max(0.0_dp, min(1.0_dp, eval(r - j + 1))))
                beta_r(j) = sqrt(max(0.0_dp, 1.0_dp - alpha_r(j) * alpha_r(j)))
                norm_ab = sqrt(alpha_r(j) * alpha_r(j) + beta_r(j) * beta_r(j))
                alpha_r(j) = alpha_r(j) / norm_ab
                beta_r(j) = beta_r(j) / norm_ab
            end if
        end do
        result%k = count(beta_r <= ab_tol)
        result%l = r - result%k

        allocate(ureal(m, m), vreal(p, p))
        ureal = 0.0_dp
        vreal = 0.0_dp
        c = 0
        do j = 1, r
            if (alpha_r(j) > ab_tol) then
                c = c + 1
                ureal(:, c) = matmul(w1, zrot(:, j)) / alpha_r(j)
            end if
        end do
        call complete_real_basis(ureal, c)
        c = 0
        do j = result%k + 1, r
            if (beta_r(j) > ab_tol) then
                c = c + 1
                vreal(:, c) = matmul(w2, zrot(:, j)) / beta_r(j)
            end if
        end do
        call complete_real_basis(vreal, c)

        allocate(temp(r, n), xinv(r, n))
        temp = vt(1:r, :)
        do j = 1, r
            temp(j, :) = sval(j) * temp(j, :)
        end do
        xinv = matmul(transpose(zrot), temp)

        allocate(rq(r, n), tau(r), qreal(n, n), rqwork(max(1, 12 * max(n, r) + 64)))
        rq = xinv
        call gerqf(r, n, rq, r, tau, rqwork, size(rqwork), info_rq)
        if (info_rq /= 0) then
            result%info = 3000 + info_rq
            return
        end if
        qreal = 0.0_dp
        do j = 1, n
            qreal(j, j) = 1.0_dp
        end do
        call ormrq('R', 'T', n, n, r, rq, r, tau, qreal, n, rqwork, size(rqwork), info_rq)
        if (info_rq /= 0) then
            result%info = 3100 + info_rq
            return
        end if

        allocate(result%alpha(n), result%beta(n), result%u(m, m), result%v(p, p), result%q(n, n), result%r(r, r))
        result%alpha = 0.0_dp
        result%beta = 0.0_dp
        result%alpha(1:r) = alpha_r
        result%beta(1:r) = beta_r
        result%u = cmplx(ureal, 0.0_dp, kind=dp)
        result%v = cmplx(vreal, 0.0_dp, kind=dp)
        result%q = cmplx(qreal, 0.0_dp, kind=dp)
        result%r = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 1, r
            result%r(j, j:r) = cmplx(rq(j, n - r + j:n), 0.0_dp, kind=dp)
        end do
    end subroutine gsvd_real

    subroutine gsvd_complex(a, b, result)
        complex(dp), intent(in) :: a(:,:) !! Complex M-by-N matrix in the matrix pair.
        complex(dp), intent(in) :: b(:,:) !! Complex P-by-N matrix sharing the column count of A.
        type(gsvd_result), intent(out) :: result !! Generalized SVD factors, common triangular factor, and status.
        integer :: m, n, p, rows, mn, r, j, c, lwork, info_eig, info_rq
        real(dp) :: rank_tol, ab_tol, norm_ab
        complex(dp), allocatable :: stack(:,:), us(:,:), vt(:,:), work(:), w1(:,:), w2(:,:)
        complex(dp), allocatable :: gram(:,:), zrot(:,:), temp(:,:), xinv(:,:), rq(:,:), tau(:), qcmp(:,:), ucmp(:,:), vcmp(:,:)
        complex(dp), allocatable :: eigwork(:), rqwork(:)
        real(dp), allocatable :: sval(:), eval(:), alpha_r(:), beta_r(:), rwork(:), eigrwork(:)

        m = size(a, 1)
        n = size(a, 2)
        p = size(b, 1)
        result%m = m
        result%n = n
        result%p = p
        if (m < 1 .or. p < 1 .or. n < 1 .or. size(b, 2) /= n) then
            result%info = -1001
            return
        end if
        if (.not. all_finite_complex(a) .or. .not. all_finite_complex(b)) then
            result%info = -1002
            return
        end if

        rows = m + p
        mn = min(rows, n)
        lwork = max(1, 12 * max(rows, n) + 8 * mn + 64)
        allocate(stack(rows, n), sval(mn), us(rows, mn), vt(mn, n), work(lwork), rwork(max(1, 5 * mn)))
        stack(1:m, :) = a
        stack(m + 1:rows, :) = b
        call gesvd('S', 'S', rows, n, stack, rows, sval, us, rows, vt, mn, work, lwork, rwork, result%info)
        if (result%info /= 0) return

        rank_tol = real(max(rows, n), dp) * epsilon(1.0_dp)
        if (mn > 0) rank_tol = rank_tol * max(1.0_dp, sval(1))
        r = count(sval > rank_tol)
        result%rank = r
        if (r == 0) then
            result%k = 0
            result%l = 0
            allocate(result%alpha(n), result%beta(n), result%u(m, m), result%v(p, p), result%q(n, n), result%r(0, 0))
            result%alpha = 0.0_dp
            result%beta = 0.0_dp
            call identity_complex(result%u)
            call identity_complex(result%v)
            call identity_complex(result%q)
            return
        end if

        allocate(w1(m, r), w2(p, r), gram(r, r), eval(r), zrot(r, r), alpha_r(r), beta_r(r))
        w1 = us(1:m, 1:r)
        w2 = us(m + 1:rows, 1:r)
        gram = matmul(conjg(transpose(w1)), w1)
        allocate(eigwork(max(1, 8 * r + 32)), eigrwork(max(1, 3 * r - 2)))
        call heev('V', 'U', r, gram, r, eval, eigwork, size(eigwork), eigrwork, info_eig)
        if (info_eig /= 0) then
            result%info = 2000 + info_eig
            return
        end if
        ab_tol = 1000.0_dp * real(max(m, max(p, n)), dp) * epsilon(1.0_dp)
        do j = 1, r
            zrot(:, j) = gram(:, r - j + 1)
            if (1.0_dp - min(1.0_dp, eval(r - j + 1)) <= ab_tol) then
                alpha_r(j) = 1.0_dp
                beta_r(j) = 0.0_dp
            else if (max(0.0_dp, eval(r - j + 1)) <= ab_tol) then
                alpha_r(j) = 0.0_dp
                beta_r(j) = 1.0_dp
            else
                alpha_r(j) = sqrt(max(0.0_dp, min(1.0_dp, eval(r - j + 1))))
                beta_r(j) = sqrt(max(0.0_dp, 1.0_dp - alpha_r(j) * alpha_r(j)))
                norm_ab = sqrt(alpha_r(j) * alpha_r(j) + beta_r(j) * beta_r(j))
                alpha_r(j) = alpha_r(j) / norm_ab
                beta_r(j) = beta_r(j) / norm_ab
            end if
        end do
        result%k = count(beta_r <= ab_tol)
        result%l = r - result%k

        allocate(ucmp(m, m), vcmp(p, p))
        ucmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
        vcmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
        c = 0
        do j = 1, r
            if (alpha_r(j) > ab_tol) then
                c = c + 1
                ucmp(:, c) = matmul(w1, zrot(:, j)) / alpha_r(j)
            end if
        end do
        call complete_complex_basis(ucmp, c)
        c = 0
        do j = result%k + 1, r
            if (beta_r(j) > ab_tol) then
                c = c + 1
                vcmp(:, c) = matmul(w2, zrot(:, j)) / beta_r(j)
            end if
        end do
        call complete_complex_basis(vcmp, c)

        allocate(temp(r, n), xinv(r, n))
        temp = vt(1:r, :)
        do j = 1, r
            temp(j, :) = sval(j) * temp(j, :)
        end do
        xinv = matmul(conjg(transpose(zrot)), temp)

        allocate(rq(r, n), tau(r), qcmp(n, n), rqwork(max(1, 12 * max(n, r) + 64)))
        rq = xinv
        call gerqf(r, n, rq, r, tau, rqwork, size(rqwork), info_rq)
        if (info_rq /= 0) then
            result%info = 3000 + info_rq
            return
        end if
        qcmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 1, n
            qcmp(j, j) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do
        call unmrq('R', 'C', n, n, r, rq, r, tau, qcmp, n, rqwork, size(rqwork), info_rq)
        if (info_rq /= 0) then
            result%info = 3100 + info_rq
            return
        end if

        allocate(result%alpha(n), result%beta(n), result%u(m, m), result%v(p, p), result%q(n, n), result%r(r, r))
        result%alpha = 0.0_dp
        result%beta = 0.0_dp
        result%alpha(1:r) = alpha_r
        result%beta(1:r) = beta_r
        result%u = ucmp
        result%v = vcmp
        result%q = qcmp
        result%r = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 1, r
            result%r(j, j:r) = rq(j, n - r + j:n)
        end do
    end subroutine gsvd_complex

    pure logical function valid_square_pair_real(a, b) result(ok)
        real(dp), intent(in) :: a(:,:) !! First matrix to validate.
        real(dp), intent(in) :: b(:,:) !! Second matrix to validate.
        integer :: n
        n = size(a, 1)
        ok = n > 0 .and. size(a, 2) == n .and. size(b, 1) == n .and. size(b, 2) == n
    end function valid_square_pair_real

    pure logical function valid_square_pair_complex(a, b) result(ok)
        complex(dp), intent(in) :: a(:,:) !! First matrix to validate.
        complex(dp), intent(in) :: b(:,:) !! Second matrix to validate.
        integer :: n
        n = size(a, 1)
        ok = n > 0 .and. size(a, 2) == n .and. size(b, 1) == n .and. size(b, 2) == n
    end function valid_square_pair_complex

    pure logical function all_finite_real(a) result(ok)
        real(dp), intent(in) :: a(:,:) !! Real matrix whose entries are checked for finiteness.
        ok = all(ieee_is_finite(a))
    end function all_finite_real

    pure logical function all_finite_complex(a) result(ok)
        complex(dp), intent(in) :: a(:,:) !! Complex matrix whose real and imaginary parts are checked for finiteness.
        ok = all(ieee_is_finite(real(a, dp))) .and. all(ieee_is_finite(aimag(a)))
    end function all_finite_complex

    pure logical function is_symmetric_real(a) result(ok)
        real(dp), intent(in) :: a(:,:) !! Real matrix tested for numerical symmetry.
        real(dp) :: scale, tol
        if (size(a, 1) /= size(a, 2)) then
            ok = .false.
            return
        end if
        scale = max(1.0_dp, maxval(abs(a)))
        tol = 100.0_dp * epsilon(1.0_dp) * scale
        ok = maxval(abs(a - transpose(a))) <= tol
    end function is_symmetric_real

    pure logical function is_hermitian_complex(a) result(ok)
        complex(dp), intent(in) :: a(:,:) !! Complex matrix tested for numerical Hermitian symmetry.
        real(dp) :: scale, tol
        if (size(a, 1) /= size(a, 2)) then
            ok = .false.
            return
        end if
        scale = max(1.0_dp, maxval(abs(a)))
        tol = 100.0_dp * epsilon(1.0_dp) * scale
        ok = maxval(abs(a - conjg(transpose(a)))) <= tol
    end function is_hermitian_complex

    subroutine unpack_real_eigenvectors(vr, wi, vectors)
        real(dp), intent(in) :: vr(:,:) !! LAPACK packed real right-eigenvector matrix.
        real(dp), intent(in) :: wi(:) !! Imaginary parts of eigenvalues identifying conjugate pairs.
        complex(dp), allocatable, intent(out) :: vectors(:,:) !! Complex right eigenvectors, one eigenvector per column.
        integer :: n, j

        n = size(wi)
        allocate(vectors(n, n))
        vectors = cmplx(0.0_dp, 0.0_dp, kind=dp)
        j = 1
        do while (j <= n)
            if (wi(j) > 0.0_dp .and. j < n) then
                vectors(:, j) = cmplx(vr(:, j), vr(:, j + 1), kind=dp)
                vectors(:, j + 1) = cmplx(vr(:, j), -vr(:, j + 1), kind=dp)
                j = j + 2
            else
                vectors(:, j) = cmplx(vr(:, j), 0.0_dp, kind=dp)
                j = j + 1
            end if
        end do
    end subroutine unpack_real_eigenvectors

    pure subroutine complex_ratio(alpha, beta, values)
        complex(dp), intent(in) :: alpha(:) !! Generalized-eigenvalue numerators.
        complex(dp), intent(in) :: beta(:) !! Generalized-eigenvalue denominators.
        complex(dp), intent(out) :: values(:) !! Stable componentwise ratios alpha/beta, including infinite roots.
        integer :: j
        real(dp) :: inf

        inf = ieee_value(1.0_dp, ieee_positive_inf)
        do j = 1, size(alpha)
            if (.not. ieee_is_finite(abs(beta(j)))) then
                values(j) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            else if (abs(beta(j)) > 0.0_dp) then
                values(j) = alpha(j) / beta(j)
            else
                values(j) = cmplx(signed_infinity(real(alpha(j), dp), inf), &
                    signed_infinity(aimag(alpha(j)), inf), kind=dp)
            end if
        end do
    end subroutine complex_ratio

    pure real(dp) function signed_infinity(x, inf) result(y)
        real(dp), intent(in) :: x !! Numerator component whose sign is transferred to infinity.
        real(dp), intent(in) :: inf !! Positive IEEE infinity for the active real kind.
        if (x > 0.0_dp) then
            y = inf
        else if (x < 0.0_dp) then
            y = -inf
        else
            y = 0.0_dp
        end if
    end function signed_infinity

    subroutine complete_real_basis(q, n_fixed)
        real(dp), intent(inout) :: q(:,:) !! Square matrix whose leading columns are orthonormalized and completed.
        integer, intent(in) :: n_fixed !! Number of leading columns already containing desired basis vectors.
        integer :: d, j, i, filled
        real(dp) :: nv
        real(dp), allocatable :: v(:)

        d = size(q, 1)
        allocate(v(d))
        filled = 0
        do j = 1, min(n_fixed, d)
            v = q(:, j)
            do i = 1, filled
                v = v - dot_product(q(:, i), v) * q(:, i)
            end do
            nv = sqrt(max(0.0_dp, dot_product(v, v)))
            if (nv > 100.0_dp * epsilon(1.0_dp)) then
                filled = filled + 1
                q(:, filled) = v / nv
            end if
        end do
        do j = 1, d
            if (filled == d) exit
            v = 0.0_dp
            v(j) = 1.0_dp
            do i = 1, filled
                v = v - dot_product(q(:, i), v) * q(:, i)
            end do
            nv = sqrt(max(0.0_dp, dot_product(v, v)))
            if (nv > 100.0_dp * epsilon(1.0_dp)) then
                filled = filled + 1
                q(:, filled) = v / nv
            end if
        end do
    end subroutine complete_real_basis

    subroutine complete_complex_basis(q, n_fixed)
        complex(dp), intent(inout) :: q(:,:) !! Square matrix whose leading columns are orthonormalized and completed.
        integer, intent(in) :: n_fixed !! Number of leading columns already containing desired basis vectors.
        integer :: d, j, i, filled
        real(dp) :: nv
        complex(dp), allocatable :: v(:)

        d = size(q, 1)
        allocate(v(d))
        filled = 0
        do j = 1, min(n_fixed, d)
            v = q(:, j)
            do i = 1, filled
                v = v - dot_product(q(:, i), v) * q(:, i)
            end do
            nv = sqrt(max(0.0_dp, real(dot_product(v, v), dp)))
            if (nv > 100.0_dp * epsilon(1.0_dp)) then
                filled = filled + 1
                q(:, filled) = v / nv
            end if
        end do
        do j = 1, d
            if (filled == d) exit
            v = cmplx(0.0_dp, 0.0_dp, kind=dp)
            v(j) = cmplx(1.0_dp, 0.0_dp, kind=dp)
            do i = 1, filled
                v = v - dot_product(q(:, i), v) * q(:, i)
            end do
            nv = sqrt(max(0.0_dp, real(dot_product(v, v), dp)))
            if (nv > 100.0_dp * epsilon(1.0_dp)) then
                filled = filled + 1
                q(:, filled) = v / nv
            end if
        end do
    end subroutine complete_complex_basis

    subroutine identity_complex(a)
        complex(dp), intent(out) :: a(:,:) !! Square output matrix set to the complex identity.
        integer :: j
        a = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do j = 1, min(size(a, 1), size(a, 2))
            a(j, j) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do
    end subroutine identity_complex

    pure logical function select_none_real(ar, ai, beta) result(sel)
        real(dp), intent(in) :: ar !! Real part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: ai !! Imaginary part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: beta !! Generalized-eigenvalue denominator.
        sel = .false.
        if (ieee_is_nan(ar) .or. ieee_is_nan(ai) .or. ieee_is_nan(beta)) sel = .false.
    end function select_none_real

    pure logical function select_negative_real(ar, ai, beta) result(sel)
        real(dp), intent(in) :: ar !! Real part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: ai !! Imaginary part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: beta !! Generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. ar * beta < 0.0_dp
        if (ieee_is_nan(ai)) sel = .false.
    end function select_negative_real

    pure logical function select_positive_real(ar, ai, beta) result(sel)
        real(dp), intent(in) :: ar !! Real part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: ai !! Imaginary part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: beta !! Generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. ar * beta > 0.0_dp
        if (ieee_is_nan(ai)) sel = .false.
    end function select_positive_real

    pure logical function select_inside_real(ar, ai, beta) result(sel)
        real(dp), intent(in) :: ar !! Real part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: ai !! Imaginary part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: beta !! Generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. abs(cmplx(ar, ai, kind=dp)) < abs(beta)
    end function select_inside_real

    pure logical function select_outside_real(ar, ai, beta) result(sel)
        real(dp), intent(in) :: ar !! Real part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: ai !! Imaginary part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: beta !! Generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. abs(cmplx(ar, ai, kind=dp)) > abs(beta)
    end function select_outside_real

    pure logical function select_real_real(ar, ai, beta) result(sel)
        real(dp), intent(in) :: ar !! Real part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: ai !! Imaginary part of the generalized-eigenvalue numerator.
        real(dp), intent(in) :: beta !! Generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. abs(ai / beta) <= 100.0_dp * epsilon(1.0_dp)
        if (ieee_is_nan(ar)) sel = .false.
    end function select_real_real

    pure logical function select_none_complex(alpha, beta) result(sel)
        complex(dp), intent(in) :: alpha !! Complex generalized-eigenvalue numerator.
        complex(dp), intent(in) :: beta !! Complex generalized-eigenvalue denominator.
        sel = .false.
        if (ieee_is_nan(real(alpha, dp)) .or. ieee_is_nan(aimag(alpha))) sel = .false.
        if (ieee_is_nan(real(beta, dp)) .or. ieee_is_nan(aimag(beta))) sel = .false.
    end function select_none_complex

    pure logical function select_negative_complex(alpha, beta) result(sel)
        complex(dp), intent(in) :: alpha !! Complex generalized-eigenvalue numerator.
        complex(dp), intent(in) :: beta !! Complex generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. real(alpha, dp) * real(beta, dp) < 0.0_dp
    end function select_negative_complex

    pure logical function select_positive_complex(alpha, beta) result(sel)
        complex(dp), intent(in) :: alpha !! Complex generalized-eigenvalue numerator.
        complex(dp), intent(in) :: beta !! Complex generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. real(alpha, dp) * real(beta, dp) > 0.0_dp
    end function select_positive_complex

    pure logical function select_inside_complex(alpha, beta) result(sel)
        complex(dp), intent(in) :: alpha !! Complex generalized-eigenvalue numerator.
        complex(dp), intent(in) :: beta !! Complex generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. abs(alpha) < abs(beta)
    end function select_inside_complex

    pure logical function select_outside_complex(alpha, beta) result(sel)
        complex(dp), intent(in) :: alpha !! Complex generalized-eigenvalue numerator.
        complex(dp), intent(in) :: beta !! Complex generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. abs(alpha) > abs(beta)
    end function select_outside_complex

    pure logical function select_real_complex(alpha, beta) result(sel)
        complex(dp), intent(in) :: alpha !! Complex generalized-eigenvalue numerator.
        complex(dp), intent(in) :: beta !! Complex generalized-eigenvalue denominator.
        sel = abs(beta) > 0.0_dp .and. abs(aimag(alpha / beta)) <= 100.0_dp * epsilon(1.0_dp)
    end function select_real_complex

end module geigen_lapack
