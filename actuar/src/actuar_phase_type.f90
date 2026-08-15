module actuar_phase_type
    use actuar_kinds, only: dp
    implicit none
    private
    public :: dphtype, pphtype, mphtype, mgfphtype, rphtype

contains

    pure subroutine eye_matrix(a)
        real(dp),intent(out)::a(:,:)
        integer::i
        a=0.0_dp
        do i=1,min(size(a,1),size(a,2));a(i,i)=1.0_dp;end do
    end subroutine eye_matrix

    pure function matrix_exp(a) result(ea)
        real(dp),intent(in)::a(:,:)
        real(dp),allocatable::ea(:,:)
        real(dp),allocatable::b(:,:),term(:,:),ident(:,:)
        real(dp)::norm1
        integer::n,s,k
        n=size(a,1)
        allocate(ea(n,n),b(n,n),term(n,n),ident(n,n))
        call eye_matrix(ident)
        norm1=maxval(sum(abs(a),dim=1))
        if(norm1<=0.5_dp) then
            s=0
        else
            s=max(0,ceiling(log(norm1/0.5_dp)/log(2.0_dp)))
        end if
        b=a/(2.0_dp**s)
        ea=ident;term=ident
        do k=1,80
            term=matmul(term,b)/real(k,dp)
            ea=ea+term
            if(maxval(abs(term))<1.0e-16_dp*max(1.0_dp,maxval(abs(ea)))) exit
        end do
        do k=1,s;ea=matmul(ea,ea);end do
    end function matrix_exp

    pure subroutine solve_linear(a,b,x,ok)
        real(dp),intent(in)::a(:,:),b(:)
        real(dp),intent(out)::x(:)
        logical,intent(out)::ok
        real(dp),allocatable::m(:,:),rhs(:),row(:)
        real(dp)::fac,piv
        integer::n,i,j,k,imax
        n=size(b);allocate(m(n,n),rhs(n),row(n));m=a;rhs=b;ok=.true.
        do k=1,n
            imax=k
            do i=k+1,n;if(abs(m(i,k))>abs(m(imax,k))) imax=i;end do
            if(abs(m(imax,k))<1.0e-14_dp) then;ok=.false.;x=0.0_dp;return;end if
            if(imax/=k) then
                row=m(k,:);m(k,:)=m(imax,:);m(imax,:)=row
                piv=rhs(k);rhs(k)=rhs(imax);rhs(imax)=piv
            end if
            do i=k+1,n
                fac=m(i,k)/m(k,k);m(i,k:n)=m(i,k:n)-fac*m(k,k:n);rhs(i)=rhs(i)-fac*rhs(k)
            end do
        end do
        x=0.0_dp
        do i=n,1,-1
            x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
        end do
    end subroutine solve_linear

    pure real(dp) function dphtype(x,prob,tmat) result(f)
        real(dp),intent(in)::x,prob(:),tmat(:,:)
        real(dp),allocatable::e(:,:),tv(:)
        integer::m
        m=size(prob)
        if(x<0.0_dp .or. size(tmat,1)/=m .or. size(tmat,2)/=m) then;f=0.0_dp;return;end if
        if(x==0.0_dp) then;f=max(0.0_dp,1.0_dp-sum(prob));return;end if
        allocate(tv(m));tv=-sum(tmat,dim=2);e=matrix_exp(x*tmat)
        f=max(0.0_dp,dot_product(prob,matmul(e,tv)))
    end function dphtype

    pure real(dp) function pphtype(x,prob,tmat) result(p)
        real(dp),intent(in)::x,prob(:),tmat(:,:)
        real(dp),allocatable::e(:,:),ones(:)
        integer::m
        m=size(prob)
        if(x<0.0_dp) then;p=0.0_dp;return;end if
        allocate(ones(m));ones=1.0_dp;e=matrix_exp(x*tmat)
        p=1.0_dp-dot_product(prob,matmul(e,ones));p=max(0.0_dp,min(1.0_dp,p))
    end function pphtype

    pure real(dp) function mphtype(order,prob,tmat) result(mom)
        integer,intent(in)::order
        real(dp),intent(in)::prob(:),tmat(:,:)
        real(dp),allocatable::v(:),rhs(:)
        logical::ok
        integer::j,n
        n=size(prob);allocate(v(n),rhs(n));v=1.0_dp
        do j=1,order
            rhs=v;call solve_linear(-tmat,rhs,v,ok)
            if(.not.ok) then;mom=huge(1.0_dp);return;end if
        end do
        mom=gamma(real(order+1,dp))*dot_product(prob,v)
    end function mphtype

    pure real(dp) function mgfphtype(s,prob,tmat) result(mgf)
        real(dp),intent(in)::s,prob(:),tmat(:,:)
        real(dp),allocatable::a(:,:),tv(:),v(:)
        logical::ok
        integer::n,i
        n=size(prob);allocate(a(n,n),tv(n),v(n));a=-tmat;tv=-sum(tmat,dim=2)
        do i=1,n;a(i,i)=a(i,i)-s;end do
        call solve_linear(a,tv,v,ok)
        if(.not.ok) then;mgf=huge(1.0_dp);return;end if
        mgf=max(0.0_dp,1.0_dp-sum(prob)+dot_product(prob,v))
    end function mgfphtype

    real(dp) function rphtype(prob,tmat) result(x)
        real(dp),intent(in)::prob(:),tmat(:,:)
        integer::m,state,j
        real(dp)::u,rate,cum
        real(dp),allocatable::rowp(:)
        m=size(prob);allocate(rowp(m+1));call random_number(u);cum=0.0_dp;state=m+1
        do j=1,m;cum=cum+prob(j);if(u<=cum) then;state=j;exit;end if;end do
        x=0.0_dp
        do while(state<=m)
            rate=-tmat(state,state)
            call random_number(u);x=x-log(max(u,tiny(1.0_dp)))/rate
            rowp=0.0_dp
            do j=1,m
                if(j/=state) rowp(j)=max(0.0_dp,tmat(state,j)/rate)
            end do
            rowp(m+1)=max(0.0_dp,1.0_dp-sum(rowp(1:m)))
            call random_number(u);cum=0.0_dp
            do j=1,m+1;cum=cum+rowp(j);if(u<=cum) then;state=j;exit;end if;end do
        end do
    end function rphtype

end module actuar_phase_type
