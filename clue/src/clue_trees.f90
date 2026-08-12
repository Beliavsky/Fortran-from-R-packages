! SPDX-License-Identifier: GPL-2.0-only
module clue_trees
    use clue_kinds, only: dp
    implicit none
    private
    public :: non_ultrametricity, non_additivity
    public :: ultrametricity_gradient, additivity_gradient
    public :: fit_ultrametric_ip, fit_ultrametric_ir, fit_addtree_ip, fit_addtree_ir
    public :: ultrametrify, centroid_addtree, is_ultrametric, is_additive
contains
    function non_ultrametricity(d,use_max) result(v)
        real(dp),intent(in)::d(:,:)
        logical,intent(in),optional::use_max
        real(dp)::v
        real(dp)::a,b,c,delta
        integer::i,j,k,n
        logical::mx
        mx=.false.
        if(present(use_max))mx=use_max
        n=size(d,1)
        v=0
        do i=1,n-2
        do j=i+1,n-1
        a=d(i,j)
        do k=j+1,n
        b=d(i,k)
        c=d(j,k)
            if(a<=b.and.a<=c)then
            delta=c-b
            else if(b<=c)then
            delta=a-c
            else
            delta=b-a
            end if
            if(mx)then
            v=max(v,abs(delta))
            else
            v=v+delta*delta
            end if
        end do
        end do
        end do
    end function
    function ultrametricity_gradient(d) result(g)
        real(dp),intent(in)::d(:,:)
        real(dp),allocatable::g(:,:)
        real(dp)::a,b,c,delta
        integer::i,j,k,n
        n=size(d,1)
        allocate(g(n,n))
        g=0
        do i=1,n-2
        do j=i+1,n-1
        a=d(i,j)
        do k=j+1,n
        b=d(i,k)
        c=d(j,k)
            if(a<=b.and.a<=c)then
            delta=2*(b-c)
            g(i,k)=g(i,k)+delta
            g(j,k)=g(j,k)-delta
            else if(b<=c)then
            delta=2*(c-a)
            g(j,k)=g(j,k)+delta
            g(i,j)=g(i,j)-delta
            else
            delta=2*(a-b)
            g(i,j)=g(i,j)+delta
            g(i,k)=g(i,k)-delta
            end if
        end do
        end do
        end do
        g=g+transpose(g)
    end function
    function non_additivity(d,use_max) result(v)
        real(dp),intent(in)::d(:,:)
        logical,intent(in),optional::use_max
        real(dp)::v
        real(dp)::a,b,c,delta
        integer::i,j,k,l,n
        logical::mx
        mx=.false.
        if(present(use_max))mx=use_max
        n=size(d,1)
        v=0
        do i=1,n-3
        do j=i+1,n-2
        do k=j+1,n-1
        do l=k+1,n
            a=d(i,j)+d(k,l)
            b=d(i,k)+d(j,l)
            c=d(i,l)+d(j,k)
            if(a<=b.and.a<=c)then
            delta=c-b
            else if(b<=c)then
            delta=a-c
            else
            delta=b-a
            end if
            if(mx)then
            v=max(v,abs(delta))
            else
            v=v+delta*delta
            end if
        end do
        end do
        end do
        end do
    end function
    function additivity_gradient(d) result(g)
        real(dp),intent(in)::d(:,:)
        real(dp),allocatable::g(:,:)
        real(dp)::a,b,c,delta
        integer::i,j,k,l,n
        n=size(d,1)
        allocate(g(n,n))
        g=0
        do i=1,n-3
        do j=i+1,n-2
        do k=j+1,n-1
        do l=k+1,n
            a=d(i,j)+d(k,l)
            b=d(i,k)+d(j,l)
            c=d(i,l)+d(j,k)
            if(a<=b.and.a<=c)then
                delta=2*(b-c)
                g(i,l)=g(i,l)-delta
                g(j,k)=g(j,k)-delta
                g(i,k)=g(i,k)+delta
                g(j,l)=g(j,l)+delta
            else if(b<=c)then
                delta=2*(c-a)
                g(i,l)=g(i,l)+delta
                g(j,k)=g(j,k)+delta
                g(i,j)=g(i,j)-delta
                g(k,l)=g(k,l)-delta
            else
                delta=2*(a-b)
                g(i,k)=g(i,k)-delta
                g(j,l)=g(j,l)-delta
                g(i,j)=g(i,j)+delta
                g(k,l)=g(k,l)+delta
            end if
        end do
        end do
        end do
        end do
        g=g+transpose(g)
    end function
    logical function is_ultrametric(d,tol) result(ok)
        real(dp),intent(in)::d(:,:)
        real(dp),intent(in),optional::tol
        real(dp)::t
        t=0
        if(present(tol))t=tol
        ok=non_ultrametricity(d,.true.)<=t
    end function
    logical function is_additive(d,tol) result(ok)
        real(dp),intent(in)::d(:,:)
        real(dp),intent(in),optional::tol
        real(dp)::t
        t=0
        if(present(tol))t=tol
        ok=non_additivity(d,.true.)<=t
    end function

    subroutine sort3(i,j,k)
        integer,intent(inout)::i,j,k
        integer::a(3),p,q,tmp
        a=[i,j,k]
        do p=1,2
        do q=p+1,3
        if(a(q)<a(p))then
        tmp=a(p)
        a(p)=a(q)
        a(q)=tmp
        end if
        end do
        end do
        i=a(1)
        j=a(2)
        k=a(3)
    end subroutine
    subroutine sort4(i,j,k,l)
        integer,intent(inout)::i,j,k,l
        integer::a(4),p,q,tmp
        a=[i,j,k,l]
        do p=1,3
        do q=p+1,4
        if(a(q)<a(p))then
        tmp=a(p)
        a(p)=a(q)
        a(q)=tmp
        end if
        end do
        end do
        i=a(1)
        j=a(2)
        k=a(3)
        l=a(4)
    end subroutine

    subroutine fit_ultrametric_ip(d,order,maxiter,tol,iterations)
        real(dp),intent(inout)::d(:,:)
        integer,intent(in),optional::order(:),maxiter
        real(dp),intent(in),optional::tol
        integer,intent(out),optional::iterations
        integer,allocatable::o(:)
        integer::mi,it,n,i1,j1,k1,i,j,k
        real(dp)::t,a,b,c,delta
        n=size(d,1)
        allocate(o(n))
        if(present(order))then
        o=order
        else
        o=[(i,i=1,n)]
        end if
        mi=10000
        if(present(maxiter))mi=maxiter
        t=1e-8_dp
        if(present(tol))t=tol
        do it=0,mi-1
        delta=0
            do i1=1,n-2
            do j1=i1+1,n-1
            do k1=j1+1,n
            i=o(i1)
            j=o(j1)
            k=o(k1)
            call sort3(i,j,k)
            a=d(i,j)
            b=d(i,k)
            c=d(j,k)
                if(a<=b.and.a<=c)then
                d(i,k)=(b+c)/2
                d(j,k)=d(i,k)
                delta=delta+abs(b-c)
                else if(b<=c)then
                d(i,j)=(c+a)/2
                d(j,k)=d(i,j)
                delta=delta+abs(c-a)
                else
                d(i,j)=(a+b)/2
                d(i,k)=d(i,j)
                delta=delta+abs(a-b)
                end if
            end do
            end do
            end do
            if(delta<t)exit
        end do
        do i=1,n
        d(i,i)=0
        do j=i+1,n
        d(j,i)=d(i,j)
        end do
        end do
        if(present(iterations))iterations=it+1
    end subroutine

    subroutine fit_ultrametric_ir(d,order,maxiter,tol,iterations)
        real(dp),intent(inout)::d(:,:)
        integer,intent(in),optional::order(:),maxiter
        real(dp),intent(in),optional::tol
        integer,intent(out),optional::iterations
        integer,allocatable::o(:)
        integer::mi,it,n,i1,j1,k1,i,j,k,n3
        real(dp)::t,a,b,c,dq,delta,tmp
        n=size(d,1)
        allocate(o(n))
        if(present(order))then
        o=order
        else
        o=[(i,i=1,n)]
        end if
        mi=10000
        if(present(maxiter))mi=maxiter
        t=1e-8_dp
        if(present(tol))t=tol
        do i=1,n-1
        do j=i+1,n
        d(i,j)=0
        end do
        end do
        n3=max(1,n-2)
        do it=0,mi-1
            do i1=1,n-2
            do j1=i1+1,n-1
            do k1=j1+1,n
            i=o(i1)
            j=o(j1)
            k=o(k1)
            call sort3(i,j,k)
            a=d(j,i)
            b=d(k,i)
            c=d(k,j)
                if(a<=b.and.a<=c)then
                dq=(c-b)/2
                d(i,k)=d(i,k)+dq
                d(j,k)=d(j,k)-dq
                else if(b<=c)then
                dq=(c-a)/2
                d(i,j)=d(i,j)+dq
                d(j,k)=d(j,k)-dq
                else
                dq=(b-a)/2
                d(i,j)=d(i,j)+dq
                d(i,k)=d(i,k)-dq
                end if
            end do
            end do
            end do
            delta=0
            do i=1,n-1
            do j=i+1,n
            tmp=d(i,j)/real(n3,dp)
            d(j,i)=d(j,i)+tmp
            d(i,j)=0
            delta=delta+abs(tmp)
            end do
            end do
            if(delta<t)exit
        end do
        do i=1,n
        d(i,i)=0
        do j=i+1,n
        d(i,j)=d(j,i)
        end do
        end do
        if(present(iterations))iterations=it+1
    end subroutine

    subroutine fit_addtree_ip(d,order,maxiter,tol,iterations)
        real(dp),intent(inout)::d(:,:)
        integer,intent(in),optional::order(:),maxiter
        real(dp),intent(in),optional::tol
        integer,intent(out),optional::iterations
        integer,allocatable::o(:)
        integer::mi,it,n,i1,j1,k1,l1,i,j,k,l
        real(dp)::t,a,b,c,dq,delta
        n=size(d,1)
        allocate(o(n))
        if(present(order))then
        o=order
        else
        o=[(i,i=1,n)]
        end if
        mi=10000
        if(present(maxiter))mi=maxiter
        t=1e-8_dp
        if(present(tol))t=tol
        do it=0,mi-1
        delta=0
            do i1=1,n-3
            do j1=i1+1,n-2
            do k1=j1+1,n-1
            do l1=k1+1,n
            i=o(i1)
            j=o(j1)
            k=o(k1)
            l=o(l1)
            call sort4(i,j,k,l)
                a=d(i,j)+d(k,l)
                b=d(i,k)+d(j,l)
                c=d(i,l)+d(j,k)
                if(a<=b.and.a<=c)then
                dq=(c-b)/4
                d(i,l)=d(i,l)-dq
                d(j,k)=d(j,k)-dq
                d(i,k)=d(i,k)+dq
                d(j,l)=d(j,l)+dq
                delta=delta+abs(c-b)
                else if(b<=c)then
                dq=(a-c)/4
                d(i,l)=d(i,l)+dq
                d(j,k)=d(j,k)+dq
                d(i,j)=d(i,j)-dq
                d(k,l)=d(k,l)-dq
                delta=delta+abs(a-c)
                else
                dq=(b-a)/4
                d(i,k)=d(i,k)-dq
                d(j,l)=d(j,l)-dq
                d(i,j)=d(i,j)+dq
                d(k,l)=d(k,l)+dq
                delta=delta+abs(b-a)
                end if
            end do
            end do
            end do
            end do
            if(delta<t)exit
        end do
        do i=1,n
        d(i,i)=0
        do j=i+1,n
        d(j,i)=d(i,j)
        end do
        end do
        if(present(iterations))iterations=it+1
    end subroutine

    subroutine fit_addtree_ir(d,order,maxiter,tol,iterations)
        real(dp),intent(inout)::d(:,:)
        integer,intent(in),optional::order(:),maxiter
        real(dp),intent(in),optional::tol
        integer,intent(out),optional::iterations
        integer,allocatable::o(:)
        integer::mi,it,n,i1,j1,k1,l1,i,j,k,l
        real(dp)::t,a,b,c,dq,delta,tmp,n3
        n=size(d,1)
        allocate(o(n))
        if(present(order))then
        o=order
        else
        o=[(i,i=1,n)]
        end if
        mi=10000
        if(present(maxiter))mi=maxiter
        t=1e-8_dp
        if(present(tol))t=tol
        do i=1,n-1
        do j=i+1,n
        d(i,j)=0
        end do
        end do
        n3=max(1.0_dp,real((n-2)*(n-3),dp)/2)
        do it=0,mi-1
            do i1=1,n-3
            do j1=i1+1,n-2
            do k1=j1+1,n-1
            do l1=k1+1,n
            i=o(i1)
            j=o(j1)
            k=o(k1)
            l=o(l1)
            call sort4(i,j,k,l)
                a=d(j,i)+d(l,k)
                b=d(k,i)+d(l,j)
                c=d(l,i)+d(k,j)
                if(a<=b.and.a<=c)then
                dq=(c-b)/4
                d(i,l)=d(i,l)-dq
                d(j,k)=d(j,k)-dq
                d(i,k)=d(i,k)+dq
                d(j,l)=d(j,l)+dq
                else if(b<=c)then
                dq=(a-c)/4
                d(i,l)=d(i,l)+dq
                d(j,k)=d(j,k)+dq
                d(i,j)=d(i,j)-dq
                d(k,l)=d(k,l)-dq
                else
                dq=(b-a)/4
                d(i,k)=d(i,k)-dq
                d(j,l)=d(j,l)-dq
                d(i,j)=d(i,j)+dq
                d(k,l)=d(k,l)+dq
                end if
            end do
            end do
            end do
            end do
            delta=0
            do i=1,n-1
            do j=i+1,n
            tmp=d(i,j)/n3
            d(j,i)=d(j,i)+tmp
            d(i,j)=0
            delta=delta+abs(tmp)
            end do
            end do
            if(delta<t)exit
        end do
        do i=1,n
        d(i,i)=0
        do j=i+1,n
        d(i,j)=d(j,i)
        end do
        end do
        if(present(iterations))iterations=it+1
    end subroutine

    function ultrametrify(d) result(u)
        ! Subdominant ultrametric: all-pairs minimax closure, equivalent to single linkage.
        real(dp),intent(in)::d(:,:)
        real(dp),allocatable::u(:,:)
        integer::i,j,k,n
        n=size(d,1)
        allocate(u(n,n))
        u=d
        do k=1,n
        do i=1,n
        do j=1,n
        u(i,j)=min(u(i,j),max(u(i,k),u(k,j)))
        end do
        end do
        end do
        do i=1,n
        u(i,i)=0
        end do
    end function

    function centroid_addtree(d) result(c)
        real(dp),intent(in)::d(:,:)
        real(dp),allocatable::c(:,:),g(:)
        integer::n,i,j
        n=size(d,1)
        allocate(c(n,n))
        if(n<=2)then
        c=0
        return
        end if
        allocate(g(n))
        g=sum(d,dim=2)/real(n-2,dp)-sum(d)/(2.0_dp*real((n-1)*(n-2),dp))
        do i=1,n
        do j=1,n
        c(i,j)=g(i)+g(j)
        end do
        c(i,i)=0
        end do
    end function
end module clue_trees
