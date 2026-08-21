module boot_linalg
    use boot_kinds, only : dp
    implicit none
    private
    public :: solve_linear, least_squares, determinant_abs, inverse_matrix
contains
    subroutine solve_linear(a,b,x,info)
        real(dp),intent(in)::a(:,:),b(:)
        real(dp),intent(out)::x(size(b))
        integer,intent(out),optional::info
        real(dp),allocatable::m(:,:),rhs(:)
        real(dp)::piv,tmp
        integer::n,i,j,k,p
        n=size(b)
        if(size(a,1)/=n .or. size(a,2)/=n)error stop "solve_linear: size mismatch"
        allocate(m(n,n),rhs(n))
        m=a
        rhs=b
        if(present(info))info=0
        do k=1,n
            p=k
            do i=k+1,n
                if(abs(m(i,k))>abs(m(p,k)))p=i
            end do
            if(abs(m(p,k))<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(m))))then
                x=0.0_dp
                if(present(info))info=1
                return
            end if
            if(p/=k)then
                do j=1,n
                tmp=m(k,j)
                m(k,j)=m(p,j)
                m(p,j)=tmp
                end do
                tmp=rhs(k)
                rhs(k)=rhs(p)
                rhs(p)=tmp
            end if
            piv=m(k,k)
            m(k,k:n)=m(k,k:n)/piv
            rhs(k)=rhs(k)/piv
            do i=1,n
                if(i==k)cycle
                tmp=m(i,k)
                m(i,k:n)=m(i,k:n)-tmp*m(k,k:n)
                rhs(i)=rhs(i)-tmp*rhs(k)
            end do
        end do
        x=rhs
    end subroutine solve_linear

    subroutine least_squares(xmat,y,beta,info,ridge)
        real(dp),intent(in)::xmat(:,:),y(:)
        real(dp),intent(out)::beta(size(xmat,2))
        integer,intent(out),optional::info
        real(dp),intent(in),optional::ridge
        real(dp),allocatable::xtx(:,:),xty(:)
        real(dp)::rr
        integer::j,p,istat
        if(size(xmat,1)/=size(y))error stop "least_squares: size mismatch"
        p=size(xmat,2)
        allocate(xtx(p,p),xty(p))
        xtx=matmul(transpose(xmat),xmat)
        xty=matmul(transpose(xmat),y)
        rr=0.0_dp
        if(present(ridge))rr=ridge
        do j=1,p
        xtx(j,j)=xtx(j,j)+rr
        end do
        call solve_linear(xtx,xty,beta,istat)
        if(present(info))info=istat
    end subroutine least_squares

    real(dp) function determinant_abs(a) result(det)
        real(dp),intent(in)::a(:,:)
        real(dp),allocatable::m(:,:)
        real(dp)::tmp,piv
        integer::n,i,j,k,p
        n=size(a,1)
        if(size(a,2)/=n)error stop "determinant_abs: square matrix required"
        allocate(m(n,n))
        m=a
        det=1.0_dp
        do k=1,n
            p=k
            do i=k+1,n
            if(abs(m(i,k))>abs(m(p,k)))p=i
            end do
            if(abs(m(p,k))<=tiny(1.0_dp))then
            det=0.0_dp
            return
            end if
            if(p/=k)then
                do j=1,n
                tmp=m(k,j)
                m(k,j)=m(p,j)
                m(p,j)=tmp
                end do
            end if
            piv=m(k,k)
            det=det*abs(piv)
            do i=k+1,n
            m(i,k+1:n)=m(i,k+1:n)-m(i,k)/piv*m(k,k+1:n)
            end do
        end do
    end function determinant_abs

    subroutine inverse_matrix(a,ainv,info)
        real(dp),intent(in)::a(:,:)
        real(dp),intent(out)::ainv(size(a,1),size(a,2))
        integer,intent(out),optional::info
        real(dp),allocatable::e(:),x(:)
        integer::n,j,istat
        n=size(a,1)
        if(size(a,2)/=n)error stop "inverse_matrix: square matrix required"
        allocate(e(n),x(n))
        ainv=0.0_dp
        if(present(info))info=0
        do j=1,n
            e=0.0_dp
            e(j)=1.0_dp
            call solve_linear(a,e,x,istat)
            if(istat/=0)then
                if(present(info))info=istat
                return
            end if
            ainv(:,j)=x
        end do
    end subroutine inverse_matrix
end module boot_linalg
