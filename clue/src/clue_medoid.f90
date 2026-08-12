! SPDX-License-Identifier: GPL-2.0-only
module clue_medoid
    use clue_kinds, only: dp
    use lpsolve, only: solve_lp, lp_result, lp_control, LP_MIN, LP_EQ, LP_LE, LP_OPTIMAL
    implicit none
    private
    public :: kmedoids_result, kmedoids, pam_kmedoids, medoid_id
    type :: kmedoids_result
        integer, allocatable :: class_ids(:)
        integer, allocatable :: medoid_ids(:)
        real(dp) :: criterion=huge(1.0_dp)
        integer :: status=1
    end type
contains
    integer function medoid_id(d) result(id)
        real(dp),intent(in)::d(:,:)
        real(dp),allocatable::s(:)
        s=sum(d,dim=2)
        id=minloc(s,dim=1)
    end function

    function kmedoids(d,k,control) result(out)
        real(dp),intent(in)::d(:,:)
        integer,intent(in)::k
        type(lp_control),intent(in),optional::control
        type(kmedoids_result)::out
        type(lp_result)::res
        integer::n,n2,nvar,ncon,i,j,q,row
        real(dp),allocatable::c(:),a(:,:),rhs(:)
        integer,allocatable::sense(:),ints(:)
        integer::mcount
        n=size(d,1)
        if(size(d,2)/=n .or. k<1 .or. k>n)return
        n2=n*n
        nvar=n2+n
        ncon=n+n2+1+n
        allocate(c(nvar),a(ncon,nvar),rhs(ncon),sense(ncon),ints(nvar))
        c=0
        a=0
        rhs=0
        sense=LP_LE
        ints=[(i,i=1,nvar)]
        q=0
        do j=1,n
        do i=1,n
        q=q+1
        c(q)=d(i,j)
        end do
        end do
        row=0
        do i=1,n
        row=row+1
        sense(row)=LP_EQ
        rhs(row)=1
        do j=1,n
        q=(j-1)*n+i
        a(row,q)=1
        end do
        end do
        do j=1,n
        do i=1,n
        row=row+1
        q=(j-1)*n+i
        a(row,q)=1
        a(row,n2+j)=-1
        sense(row)=LP_LE
        rhs(row)=0
        end do
        end do
        row=row+1
        a(row,n2+1:n2+n)=1
        sense(row)=LP_EQ
        rhs(row)=real(k,dp)
        do j=1,n
        row=row+1
        a(row,n2+j)=1
        sense(row)=LP_LE
        rhs(row)=1
        end do
        if(present(control))then
        call solve_lp(LP_MIN,c,a,sense,rhs,res,control,integer_variables=ints)
        else
        call solve_lp(LP_MIN,c,a,sense,rhs,res,integer_variables=ints)
        end if
        out%status=res%status
        if(res%status/=LP_OPTIMAL .or. .not.allocated(res%solution))return
        allocate(out%medoid_ids(k),out%class_ids(n))
        mcount=0
        do j=1,n
        if(res%solution(n2+j)>0.5_dp)then
        mcount=mcount+1
        if(mcount<=k)out%medoid_ids(mcount)=j
        end if
        end do
        do i=1,n
            out%class_ids(i)=1
            do j=1,n
            q=(j-1)*n+i
            if(res%solution(q)>0.5_dp)then
            do mcount=1,k
            if(out%medoid_ids(mcount)==j)out%class_ids(i)=mcount
            end do
            exit
            end if
            end do
        end do
        out%criterion=res%objective
    end function

    function pam_kmedoids(d,k,maxiter) result(out)
        real(dp),intent(in)::d(:,:)
        integer,intent(in)::k
        integer,intent(in),optional::maxiter
        type(kmedoids_result)::out
        integer::n,mi,i,j,it,bestj,oldj
        integer,allocatable::med(:),cls(:)
        real(dp)::best,cost,oldcost
        n=size(d,1)
        mi=100
        if(present(maxiter))mi=maxiter
        allocate(med(k),cls(n))
        med(1)=medoid_id(d)
        do j=2,k
        best=-1
        bestj=1
        do i=1,n
        if(any(med(1:j-1)==i))cycle
        cost=minval(d(:,i))
        if(cost>best)then
        best=cost
        bestj=i
        end if
        end do
        med(j)=bestj
        end do
        oldcost=huge(1.0_dp)
        cost=0.0_dp
        do it=1,mi
            cost=0
            do i=1,n
            cls(i)=minloc(d(i,med),dim=1)
            cost=cost+d(i,med(cls(i)))
            end do
            if(abs(oldcost-cost)<=1e-12_dp*(1+abs(oldcost)))exit
            oldcost=cost
            do j=1,k
            best=huge(1.0_dp)
            bestj=med(j)
            do i=1,n
            if(cls(i)/=j)cycle
            if(sum(d(pack([(oldj,oldj=1,n)],cls==j),i))<best)then
            best=sum(d(pack([(oldj,oldj=1,n)],cls==j),i))
            bestj=i
            end if
            end do
            med(j)=bestj
            end do
        end do
        allocate(out%class_ids(size(cls)),out%medoid_ids(size(med)))
        out%class_ids=cls
        out%medoid_ids=med
        out%criterion=cost
        out%status=0
    end function
end module clue_medoid
