! SPDX-License-Identifier: GPL-2.0-only
module clue_pclust
    use clue_kinds, only: dp
    use clue_partition, only: memberships_from_cross_dissimilarities, class_ids_from_membership
    implicit none
    private
    public :: pclust_result, pclust_euclidean
    type :: pclust_result
        integer, allocatable :: class_ids(:)
        real(dp), allocatable :: membership(:,:)
        real(dp), allocatable :: prototypes(:,:)
        real(dp) :: objective=huge(1.0_dp)
        integer :: iterations=0
        logical :: converged=.false.
    end type
contains
    function pclust_euclidean(x,k,m,weights,start,maxiter,reltol) result(out)
        ! Matrix-specialized version of clue::pclust using squared Euclidean D and weighted centroids.
        real(dp),intent(in)::x(:,:)
        integer,intent(in)::k
        real(dp),intent(in),optional::m,weights(:),start(:,:)
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::reltol
        type(pclust_result)::out
        real(dp)::mm,tol,oldv,newv,den
        integer::mi,it,i,j,n,p
        real(dp),allocatable::w(:),proto(:,:),d(:,:),u(:,:),pw(:)
        n=size(x,1)
        p=size(x,2)
        mm=1.0_dp
        if(present(m))mm=m
        mi=100
        if(present(maxiter))mi=maxiter
        tol=sqrt(epsilon(1.0_dp))
        if(present(reltol))tol=reltol
        allocate(w(n))
        if(present(weights))then
        w=weights
        else
        w=1
        end if
        allocate(proto(k,p))
        if(present(start))then
        proto=start
        else
        do j=1,k
        proto(j,:)=x(1+mod(j-1,n),:)
        end do
        end if
        allocate(d(n,k))
        call distances()
        newv=0.0_dp
        if(abs(mm-1.0_dp)<=epsilon(1.0_dp))then
            allocate(out%class_ids(n))
            out%class_ids=minloc(d,dim=2)
            oldv=objective_hard(out%class_ids)
            do it=1,mi
                do j=1,k
                pw=merge(w,0.0_dp,out%class_ids==j)
                den=sum(pw)
                if(den>0)proto(j,:)=matmul(pw,x)/den
                end do
                call distances()
                out%class_ids=minloc(d,dim=2)
                newv=objective_hard(out%class_ids)
                if(abs(oldv-newv)<tol*(abs(oldv)+tol))exit
                oldv=newv
            end do
            allocate(out%membership(n,k))
            out%membership=0
            do i=1,n
            out%membership(i,out%class_ids(i))=1
            end do
        else
            u=memberships_from_cross_dissimilarities(d,[mm])
            oldv=sum(spread(w,2,k)*(u**mm)*d)
            do it=1,mi
                do j=1,k
                pw=w*u(:,j)**mm
                den=sum(pw)
                if(den>0)proto(j,:)=matmul(pw,x)/den
                end do
                call distances()
                u=memberships_from_cross_dissimilarities(d,[mm])
                newv=sum(spread(w,2,k)*(u**mm)*d)
                if(abs(oldv-newv)<tol*(abs(oldv)+tol))exit
                oldv=newv
            end do
            out%membership=u
            out%class_ids=class_ids_from_membership(u)
        end if
        out%prototypes=proto
        out%objective=newv
        out%iterations=it
        out%converged=(it<=mi)
    contains
        subroutine distances()
        integer::ii,jj
        do ii=1,n
        do jj=1,k
        d(ii,jj)=sum((x(ii,:)-proto(jj,:))**2)
        end do
        end do
        end subroutine
        function objective_hard(ids) result(v)
        integer,intent(in)::ids(:)
        real(dp)::v
        integer::ii
        v=0
        do ii=1,n
        v=v+w(ii)*d(ii,ids(ii))
        end do
        end function
    end function
end module clue_pclust
