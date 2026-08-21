module boot_simplex
    use boot_kinds, only : dp
    use boot_linalg, only : solve_linear
    implicit none
    private
    public :: simplex_result, simplex_solve
    type :: simplex_result
        real(dp),allocatable :: x(:)
        real(dp) :: value=0.0_dp
        integer :: status=0 ! 1 optimal, 0 iteration limit, -1 infeasible, -2 unbounded
        integer :: iterations=0
    end type simplex_result
contains
    subroutine simplex_solve(c,result,a_le,b_le,a_ge,b_ge,a_eq,b_eq,maximize,max_iter,tol)
        real(dp),intent(in)::c(:)
        type(simplex_result),intent(out)::result
        real(dp),intent(in),optional::a_le(:,:),b_le(:),a_ge(:,:),b_ge(:),a_eq(:,:),b_eq(:)
        logical,intent(in),optional::maximize
        integer,intent(in),optional::max_iter
        real(dp),intent(in),optional::tol
        integer::n,m1,m2,m3,m,nvar,nart,row,j,sl,art,itmax,status,it
        real(dp)::eps,phase1,sgn,m_penalty
        real(dp),allocatable::a(:,:),b(:),cost1(:),cost2(:),xall(:)
        integer,allocatable::basis(:)
        n=size(c)
        m1=0
        m2=0
        m3=0
        if(present(a_le))m1=size(a_le,1)
        if(present(a_ge))m2=size(a_ge,1)
        if(present(a_eq))m3=size(a_eq,1)
        m=m1+m2+m3
        nart=m2+m3
        nvar=n+m1+m2+nart
        allocate(a(m,nvar),b(m),cost1(nvar),cost2(nvar),xall(nvar),basis(m))
        a=0.0_dp
        b=0.0_dp
        row=0
        sl=n
        art=n+m1+m2
        if(m1>0)then
            if(.not.present(b_le))error stop "simplex_solve: missing b_le"
            do j=1,m1
            row=row+1
            call add_row(a_le(j,:),b_le(j),-1,row,sl,art,a,b,basis)
            end do
        end if
        if(m2>0)then
            if(.not.present(b_ge))error stop "simplex_solve: missing b_ge"
            do j=1,m2
            row=row+1
            call add_row(a_ge(j,:),b_ge(j),1,row,sl,art,a,b,basis)
            end do
        end if
        if(m3>0)then
            if(.not.present(b_eq))error stop "simplex_solve: missing b_eq"
            do j=1,m3
            row=row+1
            call add_row(a_eq(j,:),b_eq(j),0,row,sl,art,a,b,basis)
            end do
        end if
        eps=1.0e-10_dp
        if(present(tol))eps=tol
        itmax=max(100,n+2*m)
        if(present(max_iter))itmax=max_iter
        cost1=0.0_dp
        do j=n+m1+m2+1,nvar
        cost1(j)=1.0_dp
        end do
        call revised_simplex(a,b,cost1,basis,xall,phase1,status,it,itmax,eps)
        if(status<0 .or. phase1>sqrt(eps))then
            allocate(result%x(n))
            result%x=0.0_dp
            result%status=-1
            result%value=0.0_dp
            result%iterations=it
            return
        end if
        cost2=0.0_dp
        sgn=1.0_dp
        if(present(maximize))then
        if(maximize)sgn=-1.0_dp
        end if
        cost2(1:n)=sgn*c
        m_penalty=max(1.0_dp,maxval(abs(c)))*1.0e8_dp
        do j=n+m1+m2+1,nvar
        cost2(j)=m_penalty
        end do
        call revised_simplex(a,b,cost2,basis,xall,result%value,status,it,itmax,eps)
        allocate(result%x(n))
        result%x=xall(1:n)
        result%status=status
        result%iterations=it
        result%value=dot_product(c,result%x)
    end subroutine simplex_solve

    subroutine add_row(orig,rhs,rel,row,sl,art,a,b,basis)
        real(dp),intent(in)::orig(:),rhs
        integer,intent(in)::rel,row
        integer,intent(inout)::sl,art
        real(dp),intent(inout)::a(:,:),b(:)
        integer,intent(inout)::basis(:)
        real(dp)::fac
        integer::rrel
        fac=1.0_dp
        rrel=rel
        if(rhs<0.0_dp)then
        fac=-1.0_dp
        rrel=-rel
        end if
        a(row,1:size(orig))=fac*orig
        b(row)=fac*rhs
        select case(rrel)
        case(-1)
            sl=sl+1
            a(row,sl)=1.0_dp
            basis(row)=sl
        case(1)
            sl=sl+1
            a(row,sl)=-1.0_dp
            art=art+1
            a(row,art)=1.0_dp
            basis(row)=art
        case(0)
            art=art+1
            a(row,art)=1.0_dp
            basis(row)=art
        end select
    end subroutine add_row

    subroutine revised_simplex(a,b,c,basis,x,value,status,iterations,maxit,eps)
        real(dp),intent(in)::a(:,:),b(:),c(:),eps
        integer,intent(inout)::basis(:)
        real(dp),intent(out)::x(:),value
        integer,intent(out)::status,iterations
        integer,intent(in)::maxit
        integer::m,n,it,j,enter,leave,info
        real(dp),allocatable::bmat(:,:),xb(:),cb(:),pi_vec(:),rhs(:),d(:),red(:)
        real(dp)::best,theta,ratio
        logical,allocatable::is_basic(:)
        m=size(a,1)
        n=size(a,2)
        allocate(bmat(m,m),xb(m),cb(m),pi_vec(m),rhs(m),d(m),red(n),is_basic(n))
        status=0
        x=0.0_dp
        do it=1,maxit
            do j=1,m
            bmat(:,j)=a(:,basis(j))
            cb(j)=c(basis(j))
            end do
            call solve_linear(bmat,b,xb,info)
            if(info/=0)then
            status=-1
            exit
            end if
            call solve_linear(transpose(bmat),cb,pi_vec,info)
            if(info/=0)then
            status=-1
            exit
            end if
            red=c-matmul(transpose(a),pi_vec)
            is_basic=.false.
            is_basic(basis)=.true.
            best=-eps
            enter=0
            do j=1,n
            if(.not.is_basic(j) .and. red(j)<best)then
            best=red(j)
            enter=j
            end if
            end do
            if(enter==0)then
            status=1
            exit
            end if
            rhs=a(:,enter)
            call solve_linear(bmat,rhs,d,info)
            if(info/=0)then
            status=-1
            exit
            end if
            theta=huge(1.0_dp)
            leave=0
            do j=1,m
                if(d(j)>eps)then
                ratio=xb(j)/d(j)
                if(ratio<theta-eps)then
                theta=ratio
                leave=j
                end if
                end if
            end do
            if(leave==0)then
            status=-2
            exit
            end if
            basis(leave)=enter
        end do
        iterations=min(it,maxit)
        do j=1,m
        bmat(:,j)=a(:,basis(j))
        end do
        call solve_linear(bmat,b,xb,info)
        x=0.0_dp
        if(info==0)x(basis)=xb
        value=dot_product(c,x)
    end subroutine revised_simplex
end module boot_simplex
