! Basic self-contained linear algebra used by sna-fortran.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_linalg
    use sna_kinds, only : dp, sna_eps
    implicit none
    private

    public :: solve_linear, inverse_matrix, dominant_eigenvector, jacobi_eigen_symmetric
    public :: least_squares, logistic_irls, matrix_rank

contains

    subroutine solve_linear(a, b, x, info)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: aug(:,:)
        real(dp) :: piv, fac, tmp
        integer :: n, i, j, k, p

        n = size(b)
        if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
            info = -1
            return
        end if
        allocate(aug(n,n+1))
        aug(:,1:n) = a
        aug(:,n+1) = b
        info = 0
        do k = 1, n
            p = k
            piv = abs(aug(k,k))
            do i = k+1, n
                if (abs(aug(i,k)) > piv) then
                    piv = abs(aug(i,k))
                    p = i
                end if
            end do
            if (piv <= sna_eps) then
                info = k
                x = 0.0_dp
                return
            end if
            if (p /= k) then
                do j = k, n+1
                    tmp = aug(k,j)
                    aug(k,j) = aug(p,j)
                    aug(p,j) = tmp
                end do
            end if
            do i = k+1, n
                fac = aug(i,k)/aug(k,k)
                aug(i,k:n+1) = aug(i,k:n+1) - fac*aug(k,k:n+1)
            end do
        end do
        do i = n, 1, -1
            tmp = aug(i,n+1)
            if (i < n) tmp = tmp - dot_product(aug(i,i+1:n), x(i+1:n))
            x(i) = tmp/aug(i,i)
        end do
    end subroutine solve_linear

    subroutine inverse_matrix(a, ainv, info, tol)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        integer, intent(out) :: info
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: aug(:,:)
        real(dp) :: piv, fac, tmp, atol
        integer :: n, i, j, k, p

        n = size(a,1)
        if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) then
            info = -1
            return
        end if
        atol = sna_eps
        if (present(tol)) atol = tol
        allocate(aug(n,2*n))
        aug = 0.0_dp
        aug(:,1:n) = a
        do i = 1, n
            aug(i,n+i) = 1.0_dp
        end do
        info = 0
        do k = 1, n
            p = k
            piv = abs(aug(k,k))
            do i = k+1, n
                if (abs(aug(i,k)) > piv) then
                    piv = abs(aug(i,k))
                    p = i
                end if
            end do
            if (piv <= atol) then
                info = k
                ainv = 0.0_dp
                return
            end if
            if (p /= k) then
                do j = 1, 2*n
                    tmp = aug(k,j)
                    aug(k,j) = aug(p,j)
                    aug(p,j) = tmp
                end do
            end if
            aug(k,:) = aug(k,:)/aug(k,k)
            do i = 1, n
                if (i /= k) then
                    fac = aug(i,k)
                    aug(i,:) = aug(i,:) - fac*aug(k,:)
                end if
            end do
        end do
        ainv = aug(:,n+1:2*n)
    end subroutine inverse_matrix

    subroutine dominant_eigenvector(a, vec, value, tol, maxiter, converged)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: vec(:)
        real(dp), intent(out) :: value
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        logical, intent(out), optional :: converged
        real(dp), allocatable :: vnew(:), av(:)
        real(dp) :: normv, err, atol
        integer :: n, it, mit
        logical :: ok

        n = size(a,1)
        atol = 1.0e-10_dp
        if (present(tol)) atol = tol
        mit = 100000
        if (present(maxiter)) mit = maxiter
        allocate(vnew(n), av(n))
        vec = 1.0_dp/sqrt(real(max(1,n),dp))
        ok = .false.
        do it = 1, mit
            vnew = matmul(a, vec)
            normv = sqrt(sum(vnew*vnew))
            if (normv <= sna_eps) then
                vec = 0.0_dp
                value = 0.0_dp
                ok = .true.
                exit
            end if
            vnew = vnew/normv
            if (dot_product(vnew,vec) < 0.0_dp) vnew = -vnew
            err = maxval(abs(vnew-vec))
            vec = vnew
            if (err <= atol) then
                ok = .true.
                exit
            end if
        end do
        av = matmul(a,vec)
        value = dot_product(vec,av)/max(dot_product(vec,vec),sna_eps)
        if (present(converged)) converged = ok
    end subroutine dominant_eigenvector

    subroutine jacobi_eigen_symmetric(a, values, vectors, info, tol, maxiter)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: values(:), vectors(:,:)
        integer, intent(out) :: info
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        real(dp), allocatable :: d(:,:)
        real(dp) :: app, aqq, apq, phi, c, s, dip, diq, vip, viq, atol, mx
        integer :: n, i, p, q, it, mit

        n = size(a,1)
        if (size(a,2)/=n .or. size(values)/=n .or. size(vectors,1)/=n .or. size(vectors,2)/=n) then
            info = -1
            return
        end if
        atol = 1.0e-12_dp
        if (present(tol)) atol = tol
        mit = max(50*n*n,1000)
        if (present(maxiter)) mit = maxiter
        allocate(d(n,n))
        d = 0.5_dp*(a+transpose(a))
        vectors = 0.0_dp
        do i=1,n
            vectors(i,i)=1.0_dp
        end do
        info = 1
        do it=1,mit
            mx = 0.0_dp
            p=1
            q=min(2,n)
            do i=1,n-1
                do q=i+1,n
                    if (abs(d(i,q))>mx) then
                        mx=abs(d(i,q))
                        p=i
                        ! retain q selected below via helper assignment
                    end if
                end do
            end do
            ! Redo to preserve both indices without relying on loop variable after loop.
            mx=0.0_dp
            p=1
            q=min(2,n)
            block
                integer :: ii,jj
                do ii=1,n-1
                    do jj=ii+1,n
                        if (abs(d(ii,jj))>mx) then
                            mx=abs(d(ii,jj))
                            p=ii
                            q=jj
                        end if
                    end do
                end do
            end block
            if (mx<=atol .or. n<2) then
                info=0
                exit
            end if
            app=d(p,p)
            aqq=d(q,q)
            apq=d(p,q)
            phi=0.5_dp*atan2(2.0_dp*apq,aqq-app)
            c=cos(phi)
            s=sin(phi)
            do i=1,n
                if (i/=p .and. i/=q) then
                    dip=d(i,p)
                    diq=d(i,q)
                    d(i,p)=c*dip-s*diq
                    d(p,i)=d(i,p)
                    d(i,q)=s*dip+c*diq
                    d(q,i)=d(i,q)
                end if
            end do
            d(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
            d(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
            d(p,q)=0.0_dp
            d(q,p)=0.0_dp
            do i=1,n
                vip=vectors(i,p)
                viq=vectors(i,q)
                vectors(i,p)=c*vip-s*viq
                vectors(i,q)=s*vip+c*viq
            end do
        end do
        do i=1,n
            values(i)=d(i,i)
        end do
        call sort_eigen_desc(values,vectors)
    end subroutine jacobi_eigen_symmetric

    subroutine sort_eigen_desc(values,vectors)
        real(dp), intent(inout) :: values(:), vectors(:,:)
        real(dp) :: tv
        real(dp), allocatable :: col(:)
        integer :: i,j,k,n
        n=size(values)
        allocate(col(size(vectors,1)))
        do i=1,n-1
            k=i
            do j=i+1,n
                if(values(j)>values(k)) k=j
            end do
            if(k/=i) then
                tv=values(i)
                values(i)=values(k)
                values(k)=tv
                col=vectors(:,i)
                vectors(:,i)=vectors(:,k)
                vectors(:,k)=col
            end if
        end do
    end subroutine sort_eigen_desc

    integer function matrix_rank(a, tol) result(r)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: b(:,:)
        real(dp) :: atol, fac, mx, tmp
        integer :: m,n,row,col,p,i,j
        m=size(a,1)
        n=size(a,2)
        allocate(b(m,n))
        b=a
        atol=1.0e-10_dp
        if(present(tol)) atol=tol
        r=0
        row=1
        do col=1,n
            if(row>m) exit
            p=row
            mx=abs(b(row,col))
            do i=row+1,m
                if(abs(b(i,col))>mx) then
                mx=abs(b(i,col))
                p=i
                end if
            end do
            if(mx<=atol) cycle
            if(p/=row) then
                do j=col,n
                    tmp=b(row,j)
                    b(row,j)=b(p,j)
                    b(p,j)=tmp
                end do
            end if
            do i=row+1,m
                fac=b(i,col)/b(row,col)
                b(i,col:n)=b(i,col:n)-fac*b(row,col:n)
            end do
            r=r+1
            row=row+1
        end do
    end function matrix_rank

    subroutine least_squares(x, y, beta, cov_beta, sigma2, info)
        real(dp), intent(in) :: x(:,:), y(:)
        real(dp), intent(out) :: beta(:), cov_beta(:,:), sigma2
        integer, intent(out) :: info
        real(dp), allocatable :: xtx(:,:), xty(:), inv(:,:), res(:)
        integer :: n,p,rk
        n=size(x,1)
        p=size(x,2)
        if(size(y)/=n .or. size(beta)/=p) then
        info=-1
        return
        end if
        allocate(xtx(p,p),xty(p),inv(p,p),res(n))
        xtx=matmul(transpose(x),x)
        xty=matmul(transpose(x),y)
        call inverse_matrix(xtx,inv,info,1.0e-12_dp)
        if(info/=0) then
        beta=0.0_dp
        cov_beta=0.0_dp
        sigma2=0.0_dp
        return
        end if
        beta=matmul(inv,xty)
        res=y-matmul(x,beta)
        rk=matrix_rank(x)
        if(n>rk) then
        sigma2=sum(res*res)/real(n-rk,dp)
        else
        sigma2=0.0_dp
        end if
        cov_beta=sigma2*inv
    end subroutine least_squares

    subroutine logistic_irls(x,y,beta,cov_beta,loglik,info,tol,maxiter)
        real(dp), intent(in) :: x(:,:), y(:)
        real(dp), intent(out) :: beta(:),cov_beta(:,:),loglik
        integer,intent(out)::info
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter
        real(dp),allocatable::eta(:),mu(:),w(:),z(:),xtwx(:,:),xtwz(:),inv(:,:),bnew(:)
        real(dp)::atol,den
        integer::n,p,it,mit,i,j,k
        n=size(x,1)
        p=size(x,2)
        atol=1.0e-8_dp
        if(present(tol))atol=tol
        mit=100
        if(present(maxiter))mit=maxiter
        allocate(eta(n),mu(n),w(n),z(n),xtwx(p,p),xtwz(p),inv(p,p),bnew(p))
        beta=0.0_dp
        info=1
        do it=1,mit
            eta=matmul(x,beta)
            do i=1,n
                if(eta(i)>=0.0_dp) then
                    mu(i)=1.0_dp/(1.0_dp+exp(-min(eta(i),700.0_dp)))
                else
                    mu(i)=exp(max(eta(i),-700.0_dp))/(1.0_dp+exp(max(eta(i),-700.0_dp)))
                end if
                mu(i)=min(max(mu(i),1.0e-12_dp),1.0_dp-1.0e-12_dp)
                w(i)=max(mu(i)*(1.0_dp-mu(i)),1.0e-12_dp)
                z(i)=eta(i)+(y(i)-mu(i))/w(i)
            end do
            xtwx=0.0_dp
            xtwz=0.0_dp
            do j=1,p
                do k=1,p
                    xtwx(j,k)=sum(w*x(:,j)*x(:,k))
                end do
                xtwz(j)=sum(w*x(:,j)*z)
            end do
            call inverse_matrix(xtwx,inv,info,1.0e-12_dp)
            if(info/=0)return
            bnew=matmul(inv,xtwz)
            if(maxval(abs(bnew-beta))<=atol) then
                beta=bnew
                info=0
                exit
            end if
            beta=bnew
        end do
        cov_beta=inv
        eta=matmul(x,beta)
        loglik=0.0_dp
        do i=1,n
            if(eta(i)>=0.0_dp) then
                den=eta(i)+log(1.0_dp+exp(-eta(i)))
            else
                den=log(1.0_dp+exp(eta(i)))
            end if
            loglik=loglik+y(i)*eta(i)-den
        end do
    end subroutine logistic_irls

end module sna_linalg
