! SPDX-License-Identifier: GPL-2.0-or-later
module tensora_linalg
    use tensora_kinds, only : dp, axis_name_len
    use tensora_types, only : tensor_t
    use tensora_index, only : complement, full_order, positions_by_name, product_shape
    use tensora_core, only : reorder_tensor_pos
    implicit none
    private

    type, public :: tensor_svd_t
        type(tensor_t) :: u
        type(tensor_t) :: d
        type(tensor_t) :: v
    end type tensor_svd_t

    interface inv_tensor
        module procedure inv_tensor_pos
        module procedure inv_tensor_names
    end interface inv_tensor

    interface solve_tensor
        module procedure solve_tensor_names
    end interface solve_tensor

    interface svd_tensor
        module procedure svd_tensor_pos
        module procedure svd_tensor_names
    end interface svd_tensor

    interface chol_tensor
        module procedure chol_tensor_pos
        module procedure chol_tensor_names
    end interface chol_tensor

    interface power_tensor
        module procedure power_tensor_pos
        module procedure power_tensor_names
    end interface power_tensor

    public :: inv_tensor, inv_tensor_pos, inv_tensor_names
    public :: solve_tensor, solve_tensor_names
    public :: svd_tensor, svd_tensor_pos, svd_tensor_names
    public :: chol_tensor, chol_tensor_pos, chol_tensor_names
    public :: power_tensor, power_tensor_pos, power_tensor_names
    public :: opnorm_tensor, opnorm_by_tensor, to_matrix_tensor
    public :: matrix_inverse, matrix_pseudoinverse, svd_matrix, cholesky_upper

contains

    function inv_tensor_pos(x, image, by, allow_singular, eps) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: image(:)
        integer, intent(in), optional :: by(:)
        logical, intent(in), optional :: allow_singular
        real(dp), intent(in), optional :: eps
        type(tensor_t) :: z
        integer, allocatable :: b(:), domain(:), order(:)
        type(tensor_t) :: tmp
        complex(dp), allocatable :: a(:,:), ainv(:,:), out(:,:,:)
        integer :: m, n, nbatch, q, k, p
        logical :: pseudo
        real(dp) :: tol

        pseudo = .false.
        if (present(allow_singular)) pseudo = allow_singular
        tol = 1.0e-10_dp
        if (present(eps)) tol = eps
        if (present(by)) then
            allocate(b(size(by)))
            b = by
        else
            allocate(b(0))
        end if
        if (has_overlap(image, b)) error stop "inv_tensor: image and by overlap"
        domain = complement(x%rank(), [image, b])
        order = [image, domain, b]
        tmp = reorder_tensor_pos(x, order)
        m = product_shape(x%shape(image))
        n = product_shape(x%shape(domain))
        nbatch = product_shape(x%shape(b))
        if (.not. pseudo .and. m /= n) error stop "inv_tensor: non-square mapping"
        allocate(a(m,n), out(n,m,nbatch))
        do q = 1, nbatch
            a = reshape(tmp%data((q-1)*m*n+1:q*m*n), [m,n])
            if (pseudo) then
                call matrix_pseudoinverse(a, ainv, tol)
            else
                call matrix_inverse(a, ainv)
            end if
            out(:,:,q) = ainv
        end do
        allocate(z%shape(size(domain)+size(image)+size(b)), z%axis(size(domain)+size(image)+size(b)))
        p = 0
        do k = 1, size(domain)
            p = p + 1
            z%shape(p) = x%shape(domain(k))
            z%axis(p) = x%axis(domain(k))
        end do
        do k = 1, size(image)
            p = p + 1
            z%shape(p) = x%shape(image(k))
            z%axis(p) = x%axis(image(k))
        end do
        do k = 1, size(b)
            p = p + 1
            z%shape(p) = x%shape(b(k))
            z%axis(p) = x%axis(b(k))
        end do
        allocate(z%data(size(out)))
        z%data = reshape(out, [size(out)])
    end function inv_tensor_pos

    function inv_tensor_names(x, image, by, allow_singular, eps) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: image(:)
        character(len=*), intent(in), optional :: by(:)
        logical, intent(in), optional :: allow_singular
        real(dp), intent(in), optional :: eps
        type(tensor_t) :: z
        integer, allocatable :: ip(:), bp(:)
        logical :: pseudo
        real(dp) :: tol

        ip = positions_by_name(x, image)
        pseudo = .false.
        if (present(allow_singular)) pseudo = allow_singular
        tol = 1.0e-10_dp
        if (present(eps)) tol = eps
        if (present(by)) then
            bp = positions_by_name(x, by)
            z = inv_tensor_pos(x, ip, bp, pseudo, tol)
        else
            z = inv_tensor_pos(x, ip, allow_singular=pseudo, eps=tol)
        end if
    end function inv_tensor_names

    function solve_tensor_names(a, b, ia_names, jb_names, by_names, allow_singular, eps) result(z)
        type(tensor_t), intent(in) :: a, b
        character(len=*), intent(in) :: ia_names(:), jb_names(:)
        character(len=*), intent(in), optional :: by_names(:)
        logical, intent(in), optional :: allow_singular
        real(dp), intent(in), optional :: eps
        type(tensor_t) :: z
        integer, allocatable :: ia(:), jb(:), common_a(:), common_b(:), aonly(:)
        integer, allocatable :: ka(:), lb(:), order_a(:), order_b(:)
        type(tensor_t) :: ar, br
        complex(dp), allocatable :: amat(:,:), bmat(:,:), inva(:,:), xmat(:,:)
        complex(dp), allocatable :: out(:,:,:,:)
        integer :: m, n, qcols, nc, na, ci, ai, p, k, pa, pb, ncommon, naonly
        logical :: pseudo
        real(dp) :: tol

        ia = positions_by_name(a, ia_names)
        jb = positions_by_name(b, jb_names)
        if (size(ia) /= size(jb)) error stop "solve_tensor: equation-axis rank mismatch"
        if (any(a%shape(ia) /= b%shape(jb))) error stop "solve_tensor: equation dimensions do not match"
        pseudo = .false.
        if (present(allow_singular)) pseudo = allow_singular
        tol = 1.0e-10_dp
        if (present(eps)) tol = eps

        if (present(by_names)) then
            allocate(common_a(size(by_names)), common_b(size(by_names)), aonly(size(by_names)))
            ncommon = 0
            naonly = 0
            do k = 1, size(by_names)
                pa = a%axis_pos(by_names(k))
                pb = b%axis_pos(by_names(k))
                if (pa > 0 .and. pb > 0) then
                    if (a%shape(pa) /= b%shape(pb)) error stop "solve_tensor: shared by dimension mismatch"
                    ncommon = ncommon + 1
                    common_a(ncommon) = pa
                    common_b(ncommon) = pb
                else if (pa > 0 .and. pb == 0) then
                    naonly = naonly + 1
                    aonly(naonly) = pa
                end if
            end do
            common_a = common_a(:ncommon)
            common_b = common_b(:ncommon)
            aonly = aonly(:naonly)
        else
            allocate(common_a(0), common_b(0), aonly(0))
        end if
        if (has_overlap(ia, [common_a, aonly]) .or. has_overlap(jb, common_b)) then
            error stop "solve_tensor: equation and by axes overlap"
        end if
        ka = complement(a%rank(), [ia, common_a, aonly])
        lb = complement(b%rank(), [jb, common_b])
        order_a = [ia, ka, common_a, aonly]
        order_b = [jb, lb, common_b]
        ar = reorder_tensor_pos(a, order_a)
        br = reorder_tensor_pos(b, order_b)
        m = product_shape(a%shape(ia))
        n = product_shape(a%shape(ka))
        qcols = product_shape(b%shape(lb))
        nc = product_shape(a%shape(common_a))
        na = product_shape(a%shape(aonly))
        if (.not. pseudo .and. m /= n) error stop "solve_tensor: non-square mapping"
        allocate(amat(m,n), bmat(m,qcols), out(n,qcols,nc,na))
        do ci = 1, nc
            bmat = reshape(br%data((ci-1)*m*qcols+1:ci*m*qcols), [m,qcols])
            do ai = 1, na
                p = ((ai-1)*nc + ci - 1)*m*n
                amat = reshape(ar%data(p+1:p+m*n), [m,n])
                if (pseudo) then
                    call matrix_pseudoinverse(amat, inva, tol)
                else
                    call matrix_inverse(amat, inva)
                end if
                xmat = matmul(inva, bmat)
                out(:,:,ci,ai) = xmat
            end do
        end do
        allocate(z%shape(size(ka)+size(lb)+size(common_a)+size(aonly)))
        allocate(z%axis(size(z%shape)))
        p = 0
        do k = 1, size(ka)
            p = p + 1
            z%shape(p) = a%shape(ka(k))
            z%axis(p) = a%axis(ka(k))
        end do
        do k = 1, size(lb)
            p = p + 1
            z%shape(p) = b%shape(lb(k))
            z%axis(p) = b%axis(lb(k))
        end do
        do k = 1, size(common_a)
            p = p + 1
            z%shape(p) = a%shape(common_a(k))
            z%axis(p) = a%axis(common_a(k))
        end do
        do k = 1, size(aonly)
            p = p + 1
            z%shape(p) = a%shape(aonly(k))
            z%axis(p) = a%axis(aonly(k))
        end do
        allocate(z%data(size(out)))
        z%data = reshape(out, [size(out)])
    end function solve_tensor_names

    function svd_tensor_pos(x, image, domain, by, name) result(res)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: image(:)
        integer, intent(in), optional :: domain(:), by(:)
        character(len=*), intent(in), optional :: name
        type(tensor_svd_t) :: res
        integer, allocatable :: j(:), b(:), order(:)
        type(tensor_t) :: tmp
        complex(dp), allocatable :: a(:,:), u(:,:), v(:,:)
        complex(dp), allocatable :: uall(:,:,:), vall(:,:,:)
        real(dp), allocatable :: s(:), dall(:,:)
        integer :: m, n, r, nbatch, q, p, k
        character(len=axis_name_len) :: lname

        lname = 'lambda'
        if (present(name)) lname = trim(name)
        if (present(by)) then
            allocate(b(size(by)))
            b = by
        else
            allocate(b(0))
        end if
        if (present(domain)) then
            allocate(j(size(domain)))
            j = domain
            if (.not. present(by)) b = complement(x%rank(), [image, j])
        else
            j = complement(x%rank(), [image, b])
        end if
        if (has_overlap(image, [j,b]) .or. has_overlap(j,b)) error stop "svd_tensor: overlapping axis groups"
        order = [image, j, b]
        tmp = reorder_tensor_pos(x, order)
        m = product_shape(x%shape(image))
        n = product_shape(x%shape(j))
        r = min(m,n)
        nbatch = product_shape(x%shape(b))
        allocate(uall(m,r,nbatch), vall(n,r,nbatch), dall(r,nbatch), a(m,n))
        do q = 1, nbatch
            a = reshape(tmp%data((q-1)*m*n+1:q*m*n), [m,n])
            call svd_matrix(a, u, s, v)
            uall(:,:,q) = u
            vall(:,:,q) = v
            dall(:,q) = s
        end do
        allocate(res%u%shape(size(image)+1+size(b)), res%u%axis(size(image)+1+size(b)))
        allocate(res%d%shape(1+size(b)), res%d%axis(1+size(b)))
        allocate(res%v%shape(size(j)+1+size(b)), res%v%axis(size(j)+1+size(b)))
        p = 0
        do k = 1, size(image)
            p = p + 1
            res%u%shape(p) = x%shape(image(k))
            res%u%axis(p) = x%axis(image(k))
        end do
        p = p + 1
        res%u%shape(p) = r
        res%u%axis(p) = lname
        do k = 1, size(b)
            p = p + 1
            res%u%shape(p) = x%shape(b(k))
            res%u%axis(p) = x%axis(b(k))
        end do
        res%d%shape(1) = r
        res%d%axis(1) = lname
        do k = 1, size(b)
            res%d%shape(k+1) = x%shape(b(k))
            res%d%axis(k+1) = x%axis(b(k))
        end do
        p = 0
        do k = 1, size(j)
            p = p + 1
            res%v%shape(p) = x%shape(j(k))
            res%v%axis(p) = x%axis(j(k))
        end do
        p = p + 1
        res%v%shape(p) = r
        res%v%axis(p) = lname
        do k = 1, size(b)
            p = p + 1
            res%v%shape(p) = x%shape(b(k))
            res%v%axis(p) = x%axis(b(k))
        end do
        allocate(res%u%data(size(uall)), res%d%data(size(dall)), res%v%data(size(vall)))
        res%u%data = reshape(uall, [size(uall)])
        res%d%data = cmplx(reshape(dall, [size(dall)]), 0.0_dp, dp)
        res%v%data = reshape(vall, [size(vall)])
    end function svd_tensor_pos

    function svd_tensor_names(x, image, domain, by, name) result(res)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: image(:)
        character(len=*), intent(in), optional :: domain(:), by(:)
        character(len=*), intent(in), optional :: name
        type(tensor_svd_t) :: res
        integer, allocatable :: i(:), j(:), b(:)
        character(len=axis_name_len) :: lname

        i = positions_by_name(x, image)
        lname = 'lambda'
        if (present(name)) lname = trim(name)
        if (present(domain)) then
            j = positions_by_name(x, domain)
            if (present(by)) then
                b = positions_by_name(x, by)
                res = svd_tensor_pos(x, i, j, b, lname)
            else
                res = svd_tensor_pos(x, i, j, name=lname)
            end if
        else
            if (present(by)) then
                b = positions_by_name(x, by)
                res = svd_tensor_pos(x, i, by=b, name=lname)
            else
                res = svd_tensor_pos(x, i, name=lname)
            end if
        end if
    end function svd_tensor_names

    function chol_tensor_pos(x, image, domain, name) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: image(:), domain(:)
        character(len=*), intent(in), optional :: name
        type(tensor_t) :: z
        integer, allocatable :: b(:), order(:)
        type(tensor_t) :: tmp
        complex(dp), allocatable :: a(:,:), rmat(:,:), out(:,:,:)
        integer :: m, n, nbatch, q, p, k
        character(len=axis_name_len) :: lname

        lname = 'lambda'
        if (present(name)) lname = trim(name)
        if (size(image) /= size(domain)) error stop "chol_tensor: axis-rank mismatch"
        if (any(x%shape(image) /= x%shape(domain))) error stop "chol_tensor: dimension mismatch"
        b = complement(x%rank(), [image, domain])
        order = [image, domain, b]
        tmp = reorder_tensor_pos(x, order)
        m = product_shape(x%shape(image))
        n = product_shape(x%shape(domain))
        if (m /= n) error stop "chol_tensor: matrix must be square"
        nbatch = product_shape(x%shape(b))
        allocate(a(m,n), out(n,n,nbatch))
        do q = 1, nbatch
            a = reshape(tmp%data((q-1)*m*n+1:q*m*n), [m,n])
            call cholesky_upper(a, rmat)
            out(:,:,q) = rmat
        end do
        allocate(z%shape(1+size(image)+size(b)), z%axis(1+size(image)+size(b)))
        z%shape(1) = n
        z%axis(1) = lname
        p = 1
        do k = 1, size(image)
            p = p + 1
            z%shape(p) = x%shape(image(k))
            z%axis(p) = x%axis(image(k))
        end do
        do k = 1, size(b)
            p = p + 1
            z%shape(p) = x%shape(b(k))
            z%axis(p) = x%axis(b(k))
        end do
        allocate(z%data(size(out)))
        z%data = reshape(out, [size(out)])
    end function chol_tensor_pos

    function chol_tensor_names(x, image, domain, name) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: image(:), domain(:)
        character(len=*), intent(in), optional :: name
        type(tensor_t) :: z
        integer, allocatable :: i(:), j(:)
        character(len=axis_name_len) :: lname

        i = positions_by_name(x, image)
        j = positions_by_name(x, domain)
        lname = 'lambda'
        if (present(name)) lname = trim(name)
        z = chol_tensor_pos(x, i, j, lname)
    end function chol_tensor_names

    function power_tensor_pos(x, image, domain, pwr, by) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: image(:), domain(:)
        real(dp), intent(in) :: pwr
        integer, intent(in), optional :: by(:)
        type(tensor_t) :: z, tmp
        integer, allocatable :: b(:), order(:), invorder(:)
        complex(dp), allocatable :: a(:,:), u(:,:), v(:,:), recon(:,:), work(:,:,:)
        real(dp), allocatable :: s(:)
        integer :: m, n, r, nbatch, q, k

        allocate(b(0))
        if (present(by)) then
            deallocate(b)
            allocate(b(size(by)))
            b = by
        else
            b = complement(x%rank(), [image, domain])
        end if
        order = [image, domain, b]
        tmp = reorder_tensor_pos(x, order)
        m = product_shape(x%shape(image))
        n = product_shape(x%shape(domain))
        r = min(m,n)
        nbatch = product_shape(x%shape(b))
        allocate(a(m,n), work(m,n,nbatch))
        do q = 1, nbatch
            a = reshape(tmp%data((q-1)*m*n+1:q*m*n), [m,n])
            call svd_matrix(a, u, s, v)
            do k = 1, r
                u(:,k) = u(:,k) * s(k)**pwr
            end do
            recon = matmul(u, transpose(conjg(v)))
            work(:,:,q) = recon
        end do
        tmp%data = reshape(work, [size(work)])
        allocate(invorder(x%rank()))
        do k = 1, x%rank()
            invorder(order(k)) = k
        end do
        z = reorder_tensor_pos(tmp, invorder)
    end function power_tensor_pos

    function power_tensor_names(x, image, domain, pwr, by) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: image(:), domain(:)
        real(dp), intent(in) :: pwr
        character(len=*), intent(in), optional :: by(:)
        type(tensor_t) :: z
        integer, allocatable :: i(:), j(:), b(:)

        i = positions_by_name(x, image)
        j = positions_by_name(x, domain)
        if (present(by)) then
            b = positions_by_name(x, by)
            z = power_tensor_pos(x, i, j, pwr, b)
        else
            z = power_tensor_pos(x, i, j, pwr)
        end if
    end function power_tensor_names

    real(dp) function opnorm_tensor(x, image, domain, by) result(ans)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: image(:)
        character(len=*), intent(in), optional :: domain(:), by(:)
        type(tensor_svd_t) :: sres

        if (present(domain)) then
            if (present(by)) then
                sres = svd_tensor_names(x, image, domain, by)
            else
                sres = svd_tensor_names(x, image, domain)
            end if
        else
            if (present(by)) then
                sres = svd_tensor_names(x, image, by=by)
            else
                sres = svd_tensor_names(x, image)
            end if
        end if
        ans = maxval(real(sres%d%data, dp))
    end function opnorm_tensor

    function opnorm_by_tensor(x, image, by, domain) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: image(:), by(:)
        character(len=*), intent(in), optional :: domain(:)
        type(tensor_t) :: z
        type(tensor_svd_t) :: sres
        integer :: r, nb, q

        if (present(domain)) then
            sres = svd_tensor_names(x, image, domain, by)
        else
            sres = svd_tensor_names(x, image, by=by)
        end if
        r = sres%d%shape(1)
        nb = product_shape(sres%d%shape(2:))
        allocate(z%shape(size(by)), z%axis(size(by)), z%data(nb))
        z%shape = sres%d%shape(2:)
        z%axis = sres%d%axis(2:)
        do q = 1, nb
            z%data(q) = sres%d%data((q-1)*r+1)
        end do
    end function opnorm_by_tensor

    function to_matrix_tensor(x, image, domain, by) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: image(:), domain(:)
        character(len=*), intent(in), optional :: by(:)
        type(tensor_t) :: z, tmp
        integer, allocatable :: i(:), j(:), b(:), order(:)
        integer :: k

        i = positions_by_name(x, image)
        j = positions_by_name(x, domain)
        allocate(b(0))
        if (present(by)) then
            deallocate(b)
            b = positions_by_name(x, by)
        else
            b = complement(x%rank(), [i,j])
        end if
        order = [i,j,b]
        tmp = reorder_tensor_pos(x, order)
        allocate(z%shape(2+size(b)), z%axis(2+size(b)), z%data(x%nelem()))
        z%shape(1) = product_shape(x%shape(i))
        z%shape(2) = product_shape(x%shape(j))
        z%axis(1) = 'i'
        z%axis(2) = 'j'
        do k = 1, size(b)
            z%shape(k+2) = x%shape(b(k))
            z%axis(k+2) = x%axis(b(k))
        end do
        z%data = tmp%data
    end function to_matrix_tensor

    subroutine matrix_inverse(a, ainv)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable, intent(out) :: ainv(:,:)
        complex(dp), allocatable :: aug(:,:)
        complex(dp) :: pivot, factor
        real(dp) :: best
        integer :: n, i, k, piv

        n = size(a,1)
        if (size(a,2) /= n) error stop "matrix_inverse: matrix is not square"
        allocate(aug(n,2*n), ainv(n,n))
        aug(:,1:n) = a
        aug(:,n+1:) = (0.0_dp, 0.0_dp)
        do i = 1, n
            aug(i,n+i) = (1.0_dp, 0.0_dp)
        end do
        do k = 1, n
            piv = k
            best = abs(aug(k,k))
            do i = k + 1, n
                if (abs(aug(i,k)) > best) then
                    best = abs(aug(i,k))
                    piv = i
                end if
            end do
            if (best <= 100.0_dp*epsilon(1.0_dp)) error stop "matrix_inverse: singular matrix"
            if (piv /= k) call swap_rows(aug, piv, k)
            pivot = aug(k,k)
            aug(k,:) = aug(k,:) / pivot
            do i = 1, n
                if (i == k) cycle
                factor = aug(i,k)
                if (abs(factor) > 0.0_dp) aug(i,:) = aug(i,:) - factor*aug(k,:)
            end do
        end do
        ainv = aug(:,n+1:)
    end subroutine matrix_inverse

    subroutine matrix_pseudoinverse(a, pinv, eps)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable, intent(out) :: pinv(:,:)
        real(dp), intent(in) :: eps
        complex(dp), allocatable :: u(:,:), v(:,:)
        real(dp), allocatable :: s(:)
        integer :: k, r
        real(dp) :: cutoff

        call svd_matrix(a, u, s, v)
        r = size(s)
        allocate(pinv(size(a,2),size(a,1)))
        pinv = (0.0_dp, 0.0_dp)
        if (r == 0) return
        cutoff = eps * maxval(s)
        do k = 1, r
            if (s(k) > cutoff) then
                pinv = pinv + outer_complex(v(:,k), conjg(u(:,k))) / s(k)
            end if
        end do
    end subroutine matrix_pseudoinverse

    subroutine svd_matrix(a, u, s, v)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable, intent(out) :: u(:,:), v(:,:)
        real(dp), allocatable, intent(out) :: s(:)
        complex(dp), allocatable :: ah(:,:), u2(:,:), v2(:,:)
        real(dp), allocatable :: s2(:)
        integer :: m, n

        m = size(a,1)
        n = size(a,2)
        if (m >= n) then
            call jacobi_svd_tall(a, u, s, v)
        else
            allocate(ah(n,m))
            ah = transpose(conjg(a))
            call jacobi_svd_tall(ah, u2, s2, v2)
            allocate(u(m,m), s(m), v(n,m))
            u = v2
            s = s2
            v = u2
        end if
    end subroutine svd_matrix

    subroutine jacobi_svd_tall(a, u, s, v)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable, intent(out) :: u(:,:), v(:,:)
        real(dp), allocatable, intent(out) :: s(:)
        complex(dp), allocatable :: b(:,:)
        complex(dp) :: apq, alpha
        complex(dp), allocatable :: bp(:), bq(:), vp(:), vq(:)
        real(dp) :: app, aqq, r, tau, t, c, sn, tol
        integer :: m, n, p, q, sweep, maxsweep, k
        logical :: changed

        m = size(a,1)
        n = size(a,2)
        if (m < n) error stop "jacobi_svd_tall: requires rows >= columns"
        allocate(b(m,n), v(n,n), s(n), u(m,n))
        b = a
        v = (0.0_dp, 0.0_dp)
        do k = 1, n
            v(k,k) = (1.0_dp, 0.0_dp)
        end do
        allocate(bp(m), bq(m), vp(n), vq(n))
        tol = 20.0_dp * epsilon(1.0_dp)
        maxsweep = max(30, 12*n)
        do sweep = 1, maxsweep
            changed = .false.
            do p = 1, n - 1
                do q = p + 1, n
                    app = sum(abs(b(:,p))**2)
                    aqq = sum(abs(b(:,q))**2)
                    if (app <= tiny(1.0_dp) .or. aqq <= tiny(1.0_dp)) cycle
                    apq = sum(conjg(b(:,p)) * b(:,q))
                    r = abs(apq)
                    if (r <= tol*sqrt(app*aqq)) cycle
                    changed = .true.
                    alpha = apq / r
                    tau = (aqq - app) / (2.0_dp*r)
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
                    else
                        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
                    end if
                    c = 1.0_dp / sqrt(1.0_dp + t*t)
                    sn = t*c
                    bp = b(:,p)
                    bq = b(:,q)
                    vp = v(:,p)
                    vq = v(:,q)
                    b(:,p) = c*bp - sn*conjg(alpha)*bq
                    b(:,q) = sn*bp + c*conjg(alpha)*bq
                    v(:,p) = c*vp - sn*conjg(alpha)*vq
                    v(:,q) = sn*vp + c*conjg(alpha)*vq
                end do
            end do
            if (.not. changed) exit
        end do
        do k = 1, n
            s(k) = sqrt(sum(abs(b(:,k))**2))
            if (s(k) > tiny(1.0_dp)) then
                u(:,k) = b(:,k) / s(k)
            else
                u(:,k) = (0.0_dp, 0.0_dp)
            end if
        end do
        call sort_svd(u, s, v)
    end subroutine jacobi_svd_tall

    subroutine sort_svd(u, s, v)
        complex(dp), intent(inout) :: u(:,:), v(:,:)
        real(dp), intent(inout) :: s(:)
        complex(dp), allocatable :: tmpu(:), tmpv(:)
        real(dp) :: tmps
        integer :: i, j, imax, n

        n = size(s)
        allocate(tmpu(size(u,1)), tmpv(size(v,1)))
        do i = 1, n - 1
            imax = i
            do j = i + 1, n
                if (s(j) > s(imax)) imax = j
            end do
            if (imax /= i) then
                tmps = s(i)
                s(i) = s(imax)
                s(imax) = tmps
                tmpu = u(:,i)
                u(:,i) = u(:,imax)
                u(:,imax) = tmpu
                tmpv = v(:,i)
                v(:,i) = v(:,imax)
                v(:,imax) = tmpv
            end if
        end do
    end subroutine sort_svd

    subroutine cholesky_upper(a, r)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable, intent(out) :: r(:,:)
        complex(dp), allocatable :: l(:,:)
        complex(dp) :: sumv
        real(dp) :: diagv
        integer :: n, i, j, k

        n = size(a,1)
        if (size(a,2) /= n) error stop "cholesky_upper: matrix is not square"
        allocate(l(n,n), r(n,n))
        l = (0.0_dp, 0.0_dp)
        do i = 1, n
            do j = 1, i
                sumv = a(i,j)
                do k = 1, j - 1
                    sumv = sumv - l(i,k)*conjg(l(j,k))
                end do
                if (i == j) then
                    if (abs(aimag(sumv)) > 1.0e-10_dp*max(1.0_dp,abs(real(sumv,dp)))) then
                        error stop "cholesky_upper: matrix is not Hermitian"
                    end if
                    diagv = real(sumv,dp)
                    if (diagv <= 0.0_dp) error stop "cholesky_upper: matrix is not positive definite"
                    l(i,j) = cmplx(sqrt(diagv),0.0_dp,dp)
                else
                    l(i,j) = sumv / l(j,j)
                end if
            end do
        end do
        r = transpose(conjg(l))
    end subroutine cholesky_upper

    subroutine swap_rows(a, i, j)
        complex(dp), intent(inout) :: a(:,:)
        integer, intent(in) :: i, j
        complex(dp), allocatable :: tmp(:)

        allocate(tmp(size(a,2)))
        tmp = a(i,:)
        a(i,:) = a(j,:)
        a(j,:) = tmp
    end subroutine swap_rows

    pure function outer_complex(a, b) result(c)
        complex(dp), intent(in) :: a(:), b(:)
        complex(dp) :: c(size(a),size(b))
        integer :: i, j

        do j = 1, size(b)
            do i = 1, size(a)
                c(i,j) = a(i)*b(j)
            end do
        end do
    end function outer_complex

    pure logical function has_overlap(a, b) result(overlap)
        integer, intent(in) :: a(:), b(:)
        integer :: k

        overlap = .false.
        do k = 1, size(a)
            if (any(b == a(k))) then
                overlap = .true.
                return
            end if
        end do
    end function has_overlap

end module tensora_linalg
