! Upstream license declaration: GPL (version unspecified)
module limsolve_inverse
    use limsolve_kinds, only: dp
    use limsolve_types
    use limsolve_linalg
    use quadprog, only: qp_result, solve_qp, qp_success
    implicit none
    private
    public :: solve_generalized, nnls, ldp, ldei, lsei, resolution

contains

    subroutine solve_generalized(a, b, x, tol, rank, status)
        real(dp), intent(in) :: a(:,:), b(:,:)
        real(dp), intent(out) :: x(:,:)
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: rank, status
        real(dp), allocatable :: ap(:,:)
        integer :: r, st
        if (size(b,1) /= size(a,1) .or. size(x,1) /= size(a,2) .or. &
            size(x,2) /= size(b,2)) then
            x = 0.0_dp
            if (present(rank)) rank = 0
            if (present(status)) status = LS_INVALID
            return
        end if
        allocate(ap(size(a,2),size(a,1)))
        if (present(tol)) then
            call pseudoinverse(a,ap,r,tol,st)
        else
            call pseudoinverse(a,ap,r,info=st)
        end if
        x = matmul(ap,b)
        if (present(rank)) rank = r
        if (present(status)) status = merge(LS_SUCCESS,LS_NUMERICAL,st == 0)
    end subroutine solve_generalized

    subroutine nnls(a,b,result,tol,max_iter)
        real(dp), intent(in) :: a(:,:), b(:)
        type(solve_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_iter
        logical, allocatable :: passive(:)
        real(dp), allocatable :: x(:), z(:), w(:), r(:), ap(:,:), bp(:), zp(:)
        real(dp) :: eps, alpha, best
        integer :: m, n, j, k, t, iter, lim, np, st, rankp

        m = size(a,1); n = size(a,2)
        call init_result(result,n)
        if (size(b) /= m .or. m <= 0 .or. n <= 0) return
        eps = sqrt(epsilon(1.0_dp)); if (present(tol)) eps = tol
        lim = 3*n; if (present(max_iter)) lim = max_iter
        allocate(passive(n),x(n),z(n),w(n),r(m))
        passive = .false.; x = 0.0_dp; iter = 0
        do
            r = b - matmul(a,x)
            w = matmul(transpose(a),r)
            best = eps
            t = 0
            do j = 1, n
                if (.not. passive(j) .and. w(j) > best) then
                    best = w(j); t = j
                end if
            end do
            if (t == 0) exit
            passive(t) = .true.
            do
                np = count(passive)
                allocate(ap(m,np),bp(m),zp(np))
                bp = b
                k = 0
                do j = 1, n
                    if (passive(j)) then
                        k = k + 1
                        ap(:,k) = a(:,j)
                    end if
                end do
                call least_squares(ap,bp,zp,rankp,eps,st)
                z = 0.0_dp
                k = 0
                do j = 1, n
                    if (passive(j)) then
                        k = k + 1
                        z(j) = zp(k)
                    end if
                end do
                deallocate(ap,bp,zp)
                if (all(pack(z,passive) > eps)) then
                    x = z
                    exit
                end if
                alpha = huge(1.0_dp)
                do j = 1, n
                    if (passive(j) .and. z(j) <= eps .and. x(j)-z(j) > 0.0_dp) then
                        alpha = min(alpha,x(j)/(x(j)-z(j)))
                    end if
                end do
                if (alpha >= huge(1.0_dp)/2.0_dp) alpha = 0.0_dp
                x = x + alpha*(z-x)
                do j = 1, n
                    if (passive(j) .and. x(j) <= eps) then
                        x(j) = 0.0_dp
                        passive(j) = .false.
                    end if
                end do
                iter = iter + 1
                if (iter > lim) then
                    result%status = LS_MAXITER
                    result%is_error = .true.
                    result%numiter = iter
                    result%x = x
                    result%residual_norm = -sum(min(x,0.0_dp))
                    result%solution_norm = sum(abs(matmul(a,x)-b))
                    return
                end if
            end do
            iter = iter + 1
            if (iter > lim) exit
        end do
        result%x = x
        result%residual_norm = -sum(min(x,0.0_dp))
        result%solution_norm = sum(abs(matmul(a,x)-b))
        result%numiter = iter
        if (iter > lim) then
            result%status = LS_MAXITER; result%is_error = .true.
        else
            result%status = LS_SUCCESS; result%is_error = .false.
        end if
    end subroutine nnls

    subroutine ldp(g,h,result,tol,lower,upper)
        real(dp), intent(in) :: g(:,:), h(:)
        type(solve_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: lower(:), upper(:)
        real(dp), allocatable :: gg(:,:), hh(:), dmat(:,:), dvec(:), amat(:,:)
        type(qp_result) :: qres
        real(dp) :: eps
        integer :: n
        n = size(g,2)
        call init_result(result,n)
        if (size(h) /= size(g,1) .or. n <= 0) return
        eps = sqrt(epsilon(1.0_dp)); if (present(tol)) eps = tol
        call augment_bounds(g,h,n,gg,hh,lower,upper)
        allocate(dmat(n,n),dvec(n),amat(n,size(gg,1)))
        dmat = identity_matrix(n); dvec = 0.0_dp; amat = transpose(gg)
        qres = solve_qp(dmat,dvec,amat,hh,meq=0)
        if (qres%status == qp_success) then
            result%x = qres%solution
            where(abs(result%x)<eps) result%x = 0.0_dp
            result%residual_norm = -sum(min(matmul(gg,result%x)-hh,0.0_dp))
            result%solution_norm = sum(result%x*result%x)
            result%numiter = qres%iterations(1)
            result%status = LS_SUCCESS
            result%is_error = result%residual_norm > eps
        else
            result%status = LS_INFEASIBLE
            result%is_error = .true.
        end if
    end subroutine ldp

    subroutine ldei(e,f,g,h,result,tol,lower,upper)
        real(dp), intent(in) :: e(:,:), f(:), g(:,:), h(:)
        type(solve_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: lower(:), upper(:)
        real(dp), allocatable :: ep(:,:), x0(:), z(:,:), gg(:,:), hh(:), gz(:,:), hz(:)
        type(solve_result) :: sub
        real(dp) :: eps
        integer :: n, ranke, st
        n = size(e,2)
        call init_result(result,n)
        if (size(f) /= size(e,1) .or. size(g,2) /= n .or. size(h) /= size(g,1)) return
        eps = sqrt(epsilon(1.0_dp)); if (present(tol)) eps = tol
        allocate(ep(n,size(e,1)),x0(n))
        call pseudoinverse(e,ep,ranke,eps,st)
        x0 = matmul(ep,f)
        if (sum(abs(matmul(e,x0)-f)) > eps) then
            result%status = LS_INFEASIBLE; result%is_error = .true.; result%x = x0
            return
        end if
        result%unconstrained_solution = x0
        call augment_bounds(g,h,n,gg,hh,lower,upper)
        if (size(gg,1) == 0 .or. all(matmul(gg,x0)-hh >= -eps)) then
            result%x = x0; result%status = LS_SUCCESS; result%is_error = .false.
            result%solution_norm = sum(x0*x0)
            return
        end if
        call null_space(e,z,ranke,eps,st)
        if (size(z,2) == 0) then
            result%x = x0; result%residual_norm = -sum(min(matmul(gg,x0)-hh,0.0_dp))
            result%status = LS_INFEASIBLE; result%is_error = .true.; return
        end if
        allocate(gz(size(gg,1),size(z,2)),hz(size(gg,1)))
        gz = matmul(gg,z); hz = hh - matmul(gg,x0)
        call ldp(gz,hz,sub,eps)
        if (.not. sub%succeeded()) then
            result%status = sub%status; result%is_error = .true.; result%x = x0; return
        end if
        result%x = x0 + matmul(z,sub%x)
        where(abs(result%x)<eps) result%x = 0.0_dp
        result%residual_norm = sum(abs(matmul(e,result%x)-f)) - &
            sum(min(matmul(gg,result%x)-hh,0.0_dp))
        result%solution_norm = sum(result%x*result%x)
        result%numiter = sub%numiter
        result%status = LS_SUCCESS; result%is_error = result%residual_norm > eps
    end subroutine ldei

    subroutine lsei(a,b,e,f,g,h,result,tol,wx,wa,lower,upper,fulloutput)
        real(dp), intent(in) :: a(:,:), b(:), e(:,:), f(:), g(:,:), h(:)
        type(solve_result), intent(out) :: result
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: wx(:), wa(:), lower(:), upper(:)
        logical, intent(in), optional :: fulloutput
        real(dp), allocatable :: aa(:,:), bb(:), ee(:,:), gg(:,:), hh(:), scale(:)
        real(dp), allocatable :: ep(:,:), y0(:), z(:,:), ar(:,:), br(:), gr(:,:), hr(:)
        real(dp), allocatable :: dmat(:,:), dvec(:), amat(:,:), y(:), q(:), covq(:,:), pinv_d(:,:)
        type(qp_result) :: qres
        real(dp) :: eps, ridge, colnorm
        integer :: n, me, ma, mg, rke, rka, st, j, qdim
        logical :: want_full, eq_contradictory

        n = max(size(a,2),max(size(e,2),size(g,2)))
        call init_result(result,n)
        ma=size(a,1); me=size(e,1); mg=size(g,1)
        if (n <= 0 .or. size(a,2) /= n .or. size(e,2) /= n .or. size(g,2) /= n .or. &
            size(b) /= ma .or. size(f) /= me .or. size(h) /= mg) return
        eps=sqrt(epsilon(1.0_dp)); if (present(tol)) eps=tol
        want_full=.false.; if (present(fulloutput)) want_full=fulloutput
        allocate(aa(ma,n),bb(ma),ee(me,n),scale(n))
        aa=a; bb=b; ee=e; scale=1.0_dp
        if (present(wa)) then
            if (size(wa)==ma) then
                do j=1,ma
                    aa(j,:)=aa(j,:)*wa(j); bb(j)=bb(j)*wa(j)
                end do
            end if
        end if
        call augment_bounds(g,h,n,gg,hh,lower,upper)
        if (present(wx)) then
            if (size(wx)==1) then
                do j=1,n
                    colnorm=sqrt(sum(aa(:,j)**2)+sum(ee(:,j)**2)+sum(gg(:,j)**2))
                    if (colnorm>eps) scale(j)=1.0_dp/colnorm
                end do
            else if (size(wx)==n) then
                scale=wx
            end if
        end if
        do j=1,n
            aa(:,j)=aa(:,j)*scale(j); ee(:,j)=ee(:,j)*scale(j); gg(:,j)=gg(:,j)*scale(j)
        end do

        allocate(y0(n)); y0=0.0_dp
        eq_contradictory=.false.
        if (me>0) then
            allocate(ep(n,me))
            call pseudoinverse(ee,ep,rke,eps,st)
            y0=matmul(ep,f)
            eq_contradictory=sum(abs(matmul(ee,y0)-f))>eps
            call null_space(ee,z,rke,eps,st)
        else
            rke=0; z=identity_matrix(n)
        end if
        qdim=size(z,2)
        if (qdim==0) then
            y=y0
            if (size(gg,1)>0 .and. any(matmul(gg,y)-hh < -eps)) then
                result%status=LS_INFEASIBLE; result%is_error=.true.
            else
                result%status=LS_SUCCESS; result%is_error=.false.
            end if
        else
            allocate(ar(ma,qdim),br(ma),gr(size(gg,1),qdim),hr(size(gg,1)))
            if (ma>0) then
                ar=matmul(aa,z); br=bb-matmul(aa,y0)
            else
                ar=0.0_dp; br=0.0_dp
            end if
            if (size(gg,1)>0) then
                gr=matmul(gg,z); hr=hh-matmul(gg,y0)
            else
                gr=0.0_dp; hr=0.0_dp
            end if
            allocate(dmat(qdim,qdim),dvec(qdim),amat(qdim,size(gr,1)))
            if (ma>0) then
                dmat=matmul(transpose(ar),ar); dvec=matmul(transpose(ar),br)
            else
                dmat=0.0_dp; dvec=0.0_dp
            end if
            ridge=1.0e-10_dp*max(1.0_dp,maxval(abs(dmat)))
            do j=1,qdim
                dmat(j,j)=dmat(j,j)+ridge
            end do
            amat=transpose(gr)
            qres=solve_qp(dmat,dvec,amat,hr,meq=0)
            if (qres%status/=qp_success) then
                result%status=LS_INFEASIBLE; result%is_error=.true.; y=y0
            else
                q=qres%solution; y=y0+matmul(z,q)
                result%status=LS_SUCCESS; result%is_error=.false.; result%numiter=qres%iterations(1)
            end if
            rka=matrix_rank(ar,eps)
            if (want_full) then
                allocate(pinv_d(qdim,qdim),covq(qdim,qdim))
                call pseudoinverse(dmat,pinv_d,rka,eps,st)
                covq=matmul(matmul(z,pinv_d),transpose(z))
                allocate(result%covariance(n,n)); result%covariance=0.0_dp
                do j=1,n
                    result%covariance(j,:)=scale(j)*covq(j,:)*scale
                end do
            end if
        end if
        result%x=scale*y
        where(abs(result%x)<eps) result%x=0.0_dp
        result%residual_norm=0.0_dp
        if (me>0) result%residual_norm=result%residual_norm+sum(abs(matmul(e,result%x)-f))
        if (size(gg,1)>0) then
            ! Recompute unscaled augmented constraints using x through bound helper.
            call augment_bounds(g,h,n,gr,hr,lower,upper)
            result%residual_norm=result%residual_norm-sum(min(matmul(gr,result%x)-hr,0.0_dp))
        end if
        result%solution_norm=0.0_dp
        if (ma>0) result%solution_norm=sum((matmul(aa,y)-bb)**2)
        result%rank_eq=rke
        if (qdim==0) result%rank_app=0
        if (qdim>0) result%rank_app=rka
        if (eq_contradictory) result%is_error=.false. ! Lawson-Hanson mode 1 is meaningful.
        if (result%residual_norm>eps .and. .not. eq_contradictory) result%is_error=.true.
    end subroutine lsei

    function resolution(a,tol) result(res)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tol
        type(resolution_result) :: res
        integer :: st
        allocate(res%row(size(a,1)),res%col(size(a,2)))
        if (present(tol)) then
            call resolution_matrix(a,res%row,res%col,res%nsolvable,tol,st)
        else
            call resolution_matrix(a,res%row,res%col,res%nsolvable,info=st)
        end if
    end function resolution

    subroutine augment_bounds(g,h,n,gg,hh,lower,upper)
        real(dp), intent(in) :: g(:,:),h(:)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: gg(:,:),hh(:)
        real(dp), intent(in), optional :: lower(:),upper(:)
        integer :: m,nl,nu,i,k
        logical :: scalar
        m=size(g,1); nl=0; nu=0
        if (present(lower)) then
            if (size(lower)==1) then
                if (abs(lower(1))<huge(1.0_dp)/10.0_dp) nl=n
            else
                nl=count(abs(lower)<huge(1.0_dp)/10.0_dp)
            end if
        end if
        if (present(upper)) then
            if (size(upper)==1) then
                if (abs(upper(1))<huge(1.0_dp)/10.0_dp) nu=n
            else
                nu=count(abs(upper)<huge(1.0_dp)/10.0_dp)
            end if
        end if
        allocate(gg(m+nl+nu,n),hh(m+nl+nu)); gg=0.0_dp; hh=0.0_dp
        if (m>0) then; gg(1:m,:)=g; hh(1:m)=h; end if
        k=m
        if (present(lower)) then
            scalar=size(lower)==1
            do i=1,n
                if (scalar) then
                    if (abs(lower(1))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=1.0_dp; hh(k)=lower(1)
                else if (i<=size(lower)) then
                    if (abs(lower(i))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=1.0_dp; hh(k)=lower(i)
                end if
            end do
        end if
        if (present(upper)) then
            scalar=size(upper)==1
            do i=1,n
                if (scalar) then
                    if (abs(upper(1))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=-1.0_dp; hh(k)=-upper(1)
                else if (i<=size(upper)) then
                    if (abs(upper(i))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=-1.0_dp; hh(k)=-upper(i)
                end if
            end do
        end if
    end subroutine augment_bounds

    subroutine init_result(result,n)
        type(solve_result), intent(out) :: result
        integer, intent(in) :: n
        allocate(result%x(max(0,n))); result%x=0.0_dp
        result%status=LS_INVALID; result%is_error=.true.; result%numiter=0
        result%residual_norm=0.0_dp; result%solution_norm=0.0_dp
    end subroutine init_result

end module limsolve_inverse
