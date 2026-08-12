! SPDX-License-Identifier: GPL-2.0-only
module clue_pava
    use clue_kinds, only: dp
    implicit none
    private
    public :: pava_mean, pava_median, weighted_median
contains
    function pava_mean(x,w) result(y)
        real(dp),intent(in)::x(:)
        real(dp),intent(in),optional::w(:)
        real(dp),allocatable::y(:),ww(:),val(:),wt(:)
        integer,allocatable::lo(:),hi(:)
        integer::n,m,i,j
        n=size(x)
        allocate(y(n),ww(n),val(n),wt(n),lo(n),hi(n))
        if(present(w))then
        ww=w
        else
        ww=1.0_dp
        end if
        m=0
        do i=1,n
            m=m+1
            lo(m)=i
            hi(m)=i
            wt(m)=ww(i)
            val(m)=x(i)
            do
                if(m<=1) exit
                if(val(m-1)<=val(m)) exit
                if(wt(m-1)+wt(m)>0.0_dp)then
                    val(m-1)=(wt(m-1)*val(m-1)+wt(m)*val(m))/(wt(m-1)+wt(m))
                else
                    val(m-1)=0.5_dp*(val(m-1)+val(m))
                end if
                wt(m-1)=wt(m-1)+wt(m)
                hi(m-1)=hi(m)
                m=m-1
            end do
        end do
        do i=1,m
            do j=lo(i),hi(i)
            y(j)=val(i)
            end do
        end do
    end function

    function weighted_median(x,w) result(v)
        real(dp),intent(in)::x(:),w(:)
        real(dp)::v
        real(dp),allocatable::xx(:),ww(:)
        real(dp)::half,s
        integer::i,j,n
        real(dp)::tx,tw
        n=size(x)
        allocate(xx(n),ww(n))
        xx=x
        ww=w
        do i=1,n-1
        do j=i+1,n
        if(xx(j)<xx(i))then
        tx=xx(i)
        xx(i)=xx(j)
        xx(j)=tx
        tw=ww(i)
        ww(i)=ww(j)
        ww(j)=tw
        end if
        end do
        end do
        half=0.5_dp*sum(ww)
        s=0
        v=xx(max(1,n))
        do i=1,n
        s=s+ww(i)
        if(s>=half)then
        v=xx(i)
        return
        end if
        end do
    end function

    function pava_median(x,w) result(y)
        real(dp),intent(in)::x(:)
        real(dp),intent(in),optional::w(:)
        real(dp),allocatable::y(:),ww(:)
        integer,allocatable::lo(:),hi(:)
        real(dp),allocatable::val(:)
        integer::n,m,i,j,k,nblk
        n=size(x)
        allocate(y(n),ww(n),lo(n),hi(n),val(n))
        if(present(w))then
        ww=w
        else
        ww=1
        end if
        m=0
        do i=1,n
            m=m+1
            lo(m)=i
            hi(m)=i
            val(m)=x(i)
            do
                if(m<=1) exit
                if(val(m-1)<=val(m)) exit
                hi(m-1)=hi(m)
                nblk=hi(m-1)-lo(m-1)+1
                val(m-1)=weighted_median(x(lo(m-1):hi(m-1)),ww(lo(m-1):hi(m-1)))
                m=m-1
            end do
        end do
        do k=1,m
        do j=lo(k),hi(k)
        y(j)=val(k)
        end do
        end do
    end function
end module clue_pava
