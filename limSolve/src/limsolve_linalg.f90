! Upstream license declaration: GPL (version unspecified)
module limsolve_linalg
    use limsolve_kinds, only: dp
    implicit none
    private
    public :: dense_solve, dense_inverse, symmetric_eigen_jacobi
    public :: pseudoinverse, matrix_rank, null_space, resolution_matrix
    public :: least_squares, identity_matrix

contains

    function identity_matrix(n) result(a)
        integer, intent(in) :: n
        real(dp), allocatable :: a(:,:)
        integer :: i
        allocate(a(n,n))
        a = 0.0_dp
        do i = 1, n
            a(i,i) = 1.0_dp
        end do
    end function identity_matrix

    subroutine dense_solve(a, b, x, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: m(:,:), rhs(:), tmp_row(:)
        real(dp) :: pivot, factor, tmp
        integer :: n, i, j, k, p

        n = size(a,1)
        info = 0
        x = 0.0_dp
        if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) then
            info = -1
            return
        end if
        if (n == 0) return
        allocate(m(n,n), rhs(n), tmp_row(n))
        m = a
        rhs = b
        do k = 1, n-1
            p = k
            pivot = abs(m(k,k))
            do i = k+1, n
                if (abs(m(i,k)) > pivot) then
                    pivot = abs(m(i,k))
                    p = i
                end if
            end do
            if (pivot <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(m)))) then
                info = k
                return
            end if
            if (p /= k) then
                tmp_row = m(k,:)
                m(k,:) = m(p,:)
                m(p,:) = tmp_row
                tmp = rhs(k)
                rhs(k) = rhs(p)
                rhs(p) = tmp
            end if
            do i = k+1, n
                factor = m(i,k) / m(k,k)
                m(i,k) = 0.0_dp
                do j = k+1, n
                    m(i,j) = m(i,j) - factor * m(k,j)
                end do
                rhs(i) = rhs(i) - factor * rhs(k)
            end do
        end do
        if (abs(m(n,n)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(m)))) then
            info = n
            return
        end if
        x(n) = rhs(n) / m(n,n)
        do i = n-1, 1, -1
            x(i) = (rhs(i) - dot_product(m(i,i+1:n), x(i+1:n))) / m(i,i)
        end do
    end subroutine dense_solve

    subroutine dense_inverse(a, ainv, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        integer, intent(out) :: info
        real(dp), allocatable :: e(:), x(:)
        integer :: n, j, st
        n = size(a,1)
        info = 0
        ainv = 0.0_dp
        if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) then
            info = -1
            return
        end if
        allocate(e(n), x(n))
        do j = 1, n
            e = 0.0_dp
            e(j) = 1.0_dp
            call dense_solve(a, e, x, st)
            if (st /= 0) then
                info = st
                return
            end if
            ainv(:,j) = x
        end do
    end subroutine dense_inverse

    subroutine symmetric_eigen_jacobi(a, eval, evec, info, tol, max_sweeps)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: eval(:)
        real(dp), intent(out) :: evec(:,:)
        integer, intent(out) :: info
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_sweeps
        real(dp), allocatable :: d(:,:)
        real(dp) :: threshold, app, aqq, apq, tau, t, c, s, dik, diq, vip, viq
        integer :: n, p, q, i, sweep, sweeps, imax
        real(dp) :: eps, tmp

        n = size(a,1)
        info = 0
        if (size(a,2) /= n .or. size(eval) /= n .or. &
            size(evec,1) /= n .or. size(evec,2) /= n) then
            info = -1
            return
        end if
        allocate(d(n,n))
        d = 0.5_dp * (a + transpose(a))
        evec = 0.0_dp
        do i = 1, n
            evec(i,i) = 1.0_dp
        end do
        eps = sqrt(epsilon(1.0_dp))
        if (present(tol)) eps = tol
        sweeps = max(50, 20*n*n)
        if (present(max_sweeps)) sweeps = max_sweeps

        do sweep = 1, sweeps
            threshold = 0.0_dp
            do p = 1, n-1
                do q = p+1, n
                    threshold = max(threshold, abs(d(p,q)))
                end do
            end do
            if (threshold <= eps * max(1.0_dp, maxval(abs(d)))) exit
            do p = 1, n-1
                do q = p+1, n
                    apq = d(p,q)
                    if (abs(apq) <= eps * max(1.0_dp, abs(d(p,p))+abs(d(q,q)))) cycle
                    app = d(p,p)
                    aqq = d(q,q)
                    tau = (aqq-app)/(2.0_dp*apq)
                    if (tau >= 0.0_dp) then
                        t = 1.0_dp/(tau + sqrt(1.0_dp+tau*tau))
                    else
                        t = -1.0_dp/(-tau + sqrt(1.0_dp+tau*tau))
                    end if
                    c = 1.0_dp/sqrt(1.0_dp+t*t)
                    s = t*c
                    do i = 1, n
                        if (i == p .or. i == q) cycle
                        dik = d(i,p)
                        diq = d(i,q)
                        d(i,p) = c*dik - s*diq
                        d(p,i) = d(i,p)
                        d(i,q) = s*dik + c*diq
                        d(q,i) = d(i,q)
                    end do
                    d(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
                    d(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
                    d(p,q) = 0.0_dp
                    d(q,p) = 0.0_dp
                    do i = 1, n
                        vip = evec(i,p)
                        viq = evec(i,q)
                        evec(i,p) = c*vip - s*viq
                        evec(i,q) = s*vip + c*viq
                    end do
                end do
            end do
        end do
        if (sweep > sweeps) info = 1
        do i = 1, n
            eval(i) = d(i,i)
        end do
        ! Sort descending by eigenvalue.
        do i = 1, n-1
            imax = i
            do p = i+1, n
                if (eval(p) > eval(imax)) imax = p
            end do
            if (imax /= i) then
                tmp = eval(i)
                eval(i) = eval(imax)
                eval(imax) = tmp
                do q = 1, n
                    tmp = evec(q,i)
                    evec(q,i) = evec(q,imax)
                    evec(q,imax) = tmp
                end do
            end if
        end do
    end subroutine symmetric_eigen_jacobi

    subroutine pseudoinverse(a, ap, rank, tol, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ap(:,:)
        integer, intent(out) :: rank
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        real(dp), allocatable :: ata(:,:), eval(:), v(:,:), av(:)
        real(dp) :: threshold, t
        integer :: m, n, i, st

        m = size(a,1)
        n = size(a,2)
        ap = 0.0_dp
        rank = 0
        st = 0
        if (size(ap,1) /= n .or. size(ap,2) /= m) then
            st = -1
            if (present(info)) info = st
            return
        end if
        allocate(ata(n,n), eval(n), v(n,n), av(m))
        ata = matmul(transpose(a), a)
        call symmetric_eigen_jacobi(ata, eval, v, st)
        if (st < 0) then
            if (present(info)) info = st
            return
        end if
        threshold = sqrt(epsilon(1.0_dp))
        if (present(tol)) threshold = tol
        if (n > 0) then
            t = max(0.0_dp, eval(1))
            threshold = threshold*threshold * max(1.0_dp, t)
        end if
        do i = 1, n
            if (eval(i) > threshold) then
                av = matmul(a, v(:,i))
                ap = ap + spread(v(:,i)/eval(i), 2, m) * spread(av, 1, n)
                rank = rank + 1
            end if
        end do
        if (present(info)) info = st
    end subroutine pseudoinverse

    integer function matrix_rank(a, tol) result(rank)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: ata(:,:), eval(:), v(:,:)
        real(dp) :: threshold
        integer :: n, st
        n = size(a,2)
        rank = 0
        if (n == 0) return
        allocate(ata(n,n), eval(n), v(n,n))
        ata = matmul(transpose(a),a)
        call symmetric_eigen_jacobi(ata,eval,v,st)
        if (st < 0) return
        threshold = sqrt(epsilon(1.0_dp))
        if (present(tol)) threshold = tol
        threshold = threshold*threshold * max(1.0_dp, maxval(max(0.0_dp,eval)))
        rank = count(eval > threshold)
    end function matrix_rank

    subroutine null_space(a, z, rank, tol, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable, intent(out) :: z(:,:)
        integer, intent(out), optional :: rank
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        real(dp), allocatable :: ata(:,:), eval(:), v(:,:)
        real(dp) :: threshold
        integer :: n, r, st
        n = size(a,2)
        allocate(ata(n,n), eval(n), v(n,n))
        ata = matmul(transpose(a),a)
        call symmetric_eigen_jacobi(ata,eval,v,st)
        threshold = sqrt(epsilon(1.0_dp))
        if (present(tol)) threshold = tol
        if (n > 0) threshold = threshold*threshold * max(1.0_dp,maxval(max(0.0_dp,eval)))
        r = count(eval > threshold)
        allocate(z(n,n-r))
        if (n-r > 0) z = v(:,r+1:n)
        if (present(rank)) rank = r
        if (present(info)) info = st
    end subroutine null_space

    subroutine resolution_matrix(a, rowres, colres, rank, tol, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: rowres(:), colres(:)
        integer, intent(out) :: rank
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        real(dp), allocatable :: ata(:,:), eval(:), v(:,:), av(:)
        real(dp) :: threshold
        integer :: m, n, i, st
        m = size(a,1)
        n = size(a,2)
        rowres = 0.0_dp
        colres = 0.0_dp
        rank = 0
        st = 0
        if (size(rowres) /= m .or. size(colres) /= n) then
            st = -1
            if (present(info)) info = st
            return
        end if
        allocate(ata(n,n),eval(n),v(n,n),av(m))
        ata = matmul(transpose(a),a)
        call symmetric_eigen_jacobi(ata,eval,v,st)
        threshold = sqrt(epsilon(1.0_dp))
        if (present(tol)) threshold = tol
        if (n > 0) threshold = threshold*threshold * max(1.0_dp,maxval(max(0.0_dp,eval)))
        do i = 1, n
            if (eval(i) > threshold) then
                av = matmul(a,v(:,i))/sqrt(eval(i))
                rowres = rowres + av*av
                colres = colres + v(:,i)*v(:,i)
                rank = rank + 1
            end if
        end do
        if (present(info)) info = st
    end subroutine resolution_matrix

    subroutine least_squares(a,b,x,rank,tol,info)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out), optional :: rank
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: info
        real(dp), allocatable :: ap(:,:)
        integer :: r, st
        allocate(ap(size(a,2),size(a,1)))
        if (present(tol)) then
            call pseudoinverse(a,ap,r,tol,st)
        else
            call pseudoinverse(a,ap,r,info=st)
        end if
        x = matmul(ap,b)
        if (present(rank)) rank = r
        if (present(info)) info = st
    end subroutine least_squares

end module limsolve_linalg
