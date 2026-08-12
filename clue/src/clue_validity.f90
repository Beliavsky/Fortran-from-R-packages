! SPDX-License-Identifier: GPL-2.0-only
module clue_validity
    use clue_kinds, only: dp
    implicit none
    private
    public :: dissimilarity_accounted_for, variance_accounted_for, deviance_accounted_for
    public :: silhouette_widths
contains
    function dissimilarity_accounted_for(m,d) result(v)
        real(dp),intent(in)::m(:,:),d(:,:)
        real(dp)::v,num,den,mean_d
        real(dp),allocatable::z(:),w(:,:)
        integer::k
        num=0.0_dp
        den=0.0_dp
        do k=1,size(m,2)
            z=m(:,k)
            w=spread(z,2,size(z))*spread(z,1,size(z))
            num=num+sum(w*d)
            den=den+sum(w)
        end do
        mean_d=sum(d)/real(max(1,size(d)),dp)
        if(den>0.0_dp .and. abs(mean_d)>tiny(1.0_dp))then
        v=1.0_dp-(num/den)/mean_d
        else
        v=0.0_dp
        end if
    end function
    function variance_accounted_for(u,d) result(v)
        real(dp),intent(in)::u(:,:),d(:,:)
        real(dp)::v,md,den
        md=sum(d)/real(max(1,size(d)),dp)
        den=sum((d-md)**2)
        if(den>0)then
        v=max(1.0_dp-sum((d-u)**2)/den,0.0_dp)
        else
        v=0
        end if
    end function
    function deviance_accounted_for(u,d) result(v)
        real(dp),intent(in)::u(:,:),d(:,:)
        real(dp)::v,med,den
        real(dp),allocatable::x(:)
        integer::i,j,k,n
        n=size(d,1)
        allocate(x(n*(n-1)/2))
        k=0
        do i=2,n
        do j=1,i-1
        k=k+1
        x(k)=d(i,j)
        end do
        end do
        call sort_real(x)
        med=median_sorted(x)
        den=sum(abs(d-med))
        if(den>0)then
        v=max(1.0_dp-sum(abs(d-u))/den,0.0_dp)
        else
        v=0
        end if
    end function
    function silhouette_widths(ids,d) result(s)
        integer,intent(in)::ids(:)
        real(dp),intent(in)::d(:,:)
        real(dp),allocatable::s(:)
        integer::n,k,i,j,c,ki
        real(dp)::a,b,tmp
        integer::na,nb
        n=size(ids)
        k=maxval(ids)
        allocate(s(n))
        s=0
        do i=1,n
            ki=ids(i)
            a=0
            na=0
            do j=1,n
            if(j/=i.and.ids(j)==ki)then
            a=a+d(i,j)
            na=na+1
            end if
            end do
            if(na>0)a=a/na
            b=huge(1.0_dp)
            do c=1,k
            if(c==ki)cycle
            tmp=0
            nb=0
            do j=1,n
            if(ids(j)==c)then
            tmp=tmp+d(i,j)
            nb=nb+1
            end if
            end do
            if(nb>0)b=min(b,tmp/nb)
            end do
            if(na==0)then
            s(i)=0
            else if(max(a,b)>0.and.b<huge(1.0_dp))then
            s(i)=(b-a)/max(a,b)
            else
            s(i)=0
            end if
        end do
    end function
    subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::t
    do i=1,size(x)-1
    do j=i+1,size(x)
    if(x(j)<x(i))then
    t=x(i)
    x(i)=x(j)
    x(j)=t
    end if
    end do
    end do
    end subroutine
    pure real(dp) function median_sorted(x) result(v)
    real(dp),intent(in)::x(:)
    integer::n
    n=size(x)
    if(n==0)then
    v=0
    else if(mod(n,2)==1)then
    v=x((n+1)/2)
    else
    v=0.5_dp*(x(n/2)+x(n/2+1))
    end if
    end function
end module clue_validity
