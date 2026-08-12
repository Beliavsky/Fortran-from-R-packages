! SPDX-License-Identifier: GPL-2.0-only
module clue_consensus
    use clue_kinds, only: dp
    use clue_lsap, only: solve_lsap
    use clue_partition, only: class_ids_from_membership
    use clue_trees, only: fit_ultrametric_ip
    use lpsolve, only: solve_lp, lp_result, LP_MIN, LP_EQ, LP_OPTIMAL, LP_SUBOPTIMAL
    implicit none
    private
    public :: consensus_dwh, consensus_soft_euclidean, consensus_hard_euclidean
    public :: consensus_soft_manhattan, consensus_hard_manhattan
    public :: consensus_hierarchy_euclidean, ensemble_medoid_index
contains
    subroutine match_euclidean(reference, member, matched)
        real(dp), intent(in) :: reference(:,:), member(:,:)
        real(dp), allocatable, intent(out) :: matched(:,:)
        real(dp), allocatable :: c(:,:)
        integer, allocatable :: p(:)
        c = matmul(transpose(reference), member)
        call solve_lsap(c, p, .true.)
        matched = member(:,p)
    end subroutine

    subroutine match_manhattan(reference, member, matched)
        real(dp), intent(in) :: reference(:,:), member(:,:)
        real(dp), allocatable, intent(out) :: matched(:,:)
        real(dp), allocatable :: c(:,:)
        integer, allocatable :: p(:)
        integer :: i, j
        allocate(c(size(reference,2), size(member,2)))
        do i = 1, size(c,1)
            do j = 1, size(c,2)
                c(i,j) = sum(abs(reference(:,i) - member(:,j)))
            end do
        end do
        call solve_lsap(c, p)
        matched = member(:,p)
    end subroutine

    subroutine stochastify(m)
        real(dp), intent(inout) :: m(:,:)
        integer :: i
        real(dp) :: s
        m = max(m, 0.0_dp)
        do i = 1, size(m,1)
            s = sum(m(i,:))
            if (s > 0.0_dp) then
                m(i,:) = m(i,:) / s
            else if (size(m,2) > 0) then
                m(i,:) = 1.0_dp / real(size(m,2),dp)
            end if
        end do
    end subroutine

    function consensus_dwh(memberships, weights, k) result(m)
        ! memberships(n,kmax,b).  This follows clue's sequential DWH update.
        real(dp), intent(in) :: memberships(:,:,:)
        real(dp), intent(in), optional :: weights(:)
        integer, intent(in), optional :: k
        real(dp), allocatable :: m(:,:), mb(:,:), w(:)
        real(dp) :: s, cumulative
        integer :: b, kk, kmax, n
        n = size(memberships,1)
        kmax = size(memberships,2)
        kk = kmax
        if (present(k)) kk = min(max(1,k),kmax)
        allocate(w(size(memberships,3)))
        if (present(weights)) then
            w = weights
        else
            w = 1.0_dp
        end if
        allocate(m(n,kmax))
        m = memberships(:,:,1)
        cumulative = w(1)
        do b = 2, size(memberships,3)
            call match_euclidean(m, memberships(:,:,b), mb)
            cumulative = cumulative + w(b)
            if (abs(cumulative)>tiny(1.0_dp)) then
                s = w(b) / cumulative
            else
                s = 0.0_dp
            end if
            m = (1.0_dp-s)*m + s*mb
            if (kk < kmax) call project_leading(m,kk)
        end do
        m = m(:,1:kk)
        call stochastify(m)
    end function

    subroutine project_leading(m,k)
        real(dp), intent(inout) :: m(:,:)
        integer, intent(in) :: k
        integer :: i,j,jbest
        if (k >= size(m,2)) return
        do i=1,size(m,1)
            do j=k+1,size(m,2)
                jbest=maxloc(m(i,1:k),dim=1)
                m(i,jbest)=m(i,jbest)+m(i,j)
                m(i,j)=0.0_dp
            end do
        end do
    end subroutine

    function consensus_soft_euclidean(memberships,weights,k,maxiter,reltol) result(m)
        real(dp), intent(in) :: memberships(:,:,:)
        real(dp), intent(in), optional :: weights(:)
        integer, intent(in), optional :: k,maxiter
        real(dp), intent(in), optional :: reltol
        real(dp), allocatable :: m(:,:), aligned(:,:,:), mb(:,:), w(:)
        real(dp) :: oldv,newv,tol
        integer :: kk,mi,it,b,kmax,n
        n=size(memberships,1)
        kmax=size(memberships,2)
        kk=kmax
        if(present(k))kk=min(max(1,k),kmax)
        mi=100
        if(present(maxiter))mi=maxiter
        tol=sqrt(epsilon(1.0_dp))
        if(present(reltol))tol=reltol
        allocate(w(size(memberships,3)))
        if(present(weights))then
        w=weights
        else
        w=1
        end if
        if(abs(sum(w))>tiny(1.0_dp))w=w/sum(w)
        m=consensus_dwh(memberships,w,kk)
        allocate(aligned(n,kmax,size(memberships,3)))
        if(size(m,2)<kmax)then
            mb=m
            deallocate(m)
            allocate(m(n,kmax))
            m=0.0_dp
            m(:,1:kk)=mb
        end if
        do b=1,size(memberships,3)
        call match_euclidean(m,memberships(:,:,b),mb)
        aligned(:,:,b)=mb
        end do
        oldv=objective_e2(m,aligned,w)
        do it=1,mi
            m=0
            do b=1,size(aligned,3)
            m=m+w(b)*aligned(:,:,b)
            end do
            if(kk<kmax)call project_leading(m,kk)
            do b=1,size(memberships,3)
            call match_euclidean(m,memberships(:,:,b),mb)
            aligned(:,:,b)=mb
            end do
            newv=objective_e2(m,aligned,w)
            if(abs(oldv-newv)<tol*(abs(oldv)+tol))exit
            oldv=newv
        end do
        m=m(:,1:kk)
        call stochastify(m)
    end function

    pure function objective_e2(m,a,w) result(v)
        real(dp),intent(in)::m(:,:),a(:,:,:),w(:)
        real(dp)::v
        integer::b
        v=0
        do b=1,size(a,3)
        v=v+w(b)*sum((a(:,:,b)-m)**2)
        end do
    end function

    function consensus_hard_euclidean(memberships,weights,k,maxiter,reltol) result(ids)
        real(dp),intent(in)::memberships(:,:,:)
        real(dp),intent(in),optional::weights(:)
        integer,intent(in),optional::k,maxiter
        real(dp),intent(in),optional::reltol
        integer,allocatable::ids(:)
        call consensus_hard_aos(memberships,weights,k,maxiter,reltol,.false.,ids)
    end function

    function consensus_soft_manhattan(memberships,weights,k,maxiter,reltol) result(m)
        real(dp), intent(in) :: memberships(:,:,:)
        real(dp), intent(in), optional :: weights(:)
        integer, intent(in), optional :: k, maxiter
        real(dp), intent(in), optional :: reltol
        real(dp), allocatable :: m(:,:), aligned(:,:,:), mb(:,:), w(:)
        real(dp) :: oldv, newv, tol
        integer :: b,it,mi,kk,kmax,n
        n=size(memberships,1)
        kmax=size(memberships,2)
        kk=kmax
        if(present(k)) kk=min(max(1,k),kmax)
        mi=100
        if(present(maxiter)) mi=maxiter
        tol=sqrt(epsilon(1.0_dp))
        if(present(reltol)) tol=reltol
        allocate(w(size(memberships,3)))
        if(present(weights)) then
        w=weights
        else
        w=1.0_dp
        end if
        if(abs(sum(w))>tiny(1.0_dp)) w=w/sum(w)
        allocate(m(n,kmax),aligned(n,kmax,size(memberships,3)))
        m=memberships(:,:,1)
        do b=1,size(memberships,3)
            call match_manhattan(m,memberships(:,:,b),mb)
            aligned(:,:,b)=mb
        end do
        oldv=objective_l1(m,aligned,w)
        do it=1,mi
            call fit_l1_stochastic_rows(aligned, w, kk, m)
            if(kk<kmax) m(:,kk+1:kmax)=0.0_dp
            do b=1,size(memberships,3)
                call match_manhattan(m,memberships(:,:,b),mb)
                aligned(:,:,b)=mb
            end do
            newv=objective_l1(m,aligned,w)
            if(abs(oldv-newv)<tol*(abs(oldv)+tol)) exit
            oldv=newv
        end do
        m=m(:,1:kk)
        call stochastify(m)
    end function

    pure function objective_l1(m,a,w) result(v)
        real(dp),intent(in)::m(:,:),a(:,:,:),w(:)
        real(dp)::v
        integer::b
        v=0.0_dp
        do b=1,size(a,3)
            v=v+w(b)*sum(abs(a(:,:,b)-m))
        end do
    end function

    function consensus_hard_manhattan(memberships,weights,k,maxiter,reltol) result(ids)
        real(dp),intent(in)::memberships(:,:,:)
        real(dp),intent(in),optional::weights(:)
        integer,intent(in),optional::k,maxiter
        real(dp),intent(in),optional::reltol
        integer,allocatable::ids(:)
        call consensus_hard_aos(memberships,weights,k,maxiter,reltol,.true.,ids)
    end function

    subroutine fit_l1_stochastic_rows(aligned,w,k,m)
        real(dp), intent(in) :: aligned(:,:,:), w(:)
        integer, intent(in) :: k
        real(dp), intent(out) :: m(:,:)
        real(dp), allocatable :: obj(:), a(:,:), rhs(:)
        integer, allocatable :: sense(:)
        type(lp_result) :: res
        integer :: n,bcount,l,ncon,i,b,j,base,ix
        n=size(aligned,1)
        bcount=size(aligned,3)
        l=(2*bcount+1)*k
        ncon=bcount*k+1
        allocate(obj(l),a(ncon,l),rhs(ncon),sense(ncon))
        obj=0.0_dp
        do b=1,bcount
            base=2*(b-1)*k
            obj(base+1:base+2*k)=w(b)
        end do
        sense=LP_EQ
        m=0.0_dp
        do i=1,n
            a=0.0_dp
            rhs=0.0_dp
            do b=1,bcount
                base=2*(b-1)*k
                do j=1,k
                    ix=(b-1)*k+j
                    a(ix,base+j)=1.0_dp
                    a(ix,base+k+j)=-1.0_dp
                    a(ix,2*bcount*k+j)=1.0_dp
                    rhs(ix)=aligned(i,j,b)
                end do
            end do
            a(ncon,2*bcount*k+1:(2*bcount+1)*k)=1.0_dp
            rhs(ncon)=1.0_dp
            call solve_lp(LP_MIN,obj,a,sense,rhs,res)
            if(res%status==LP_OPTIMAL .or. res%status==LP_SUBOPTIMAL) then
                m(i,1:k)=res%solution(2*bcount*k+1:(2*bcount+1)*k)
            else
                m(i,1:k)=sum(aligned(i,1:k,:)*reshape(w,[1,bcount]),dim=2)
                if(sum(m(i,1:k))>0.0_dp) m(i,1:k)=m(i,1:k)/sum(m(i,1:k))
            end if
        end do
    end subroutine

    subroutine consensus_hard_aos(memberships,weights,k,maxiter,reltol,use_manhattan,ids)
        real(dp),intent(in)::memberships(:,:,:)
        real(dp),intent(in),optional::weights(:)
        integer,intent(in),optional::k,maxiter
        real(dp),intent(in),optional::reltol
        logical,intent(in)::use_manhattan
        integer,allocatable,intent(out)::ids(:)
        real(dp),allocatable::m(:,:),aligned(:,:,:),mb(:,:),w(:),meanm(:,:)
        real(dp)::oldv,newv,tol
        integer::n,kmax,kk,mi,it,b,i,jbest
        n=size(memberships,1)
        kmax=size(memberships,2)
        kk=kmax
        if(present(k))kk=min(max(1,k),kmax)
        mi=100
        if(present(maxiter))mi=maxiter
        tol=sqrt(epsilon(1.0_dp))
        if(present(reltol))tol=reltol
        allocate(w(size(memberships,3)))
        if(present(weights))then
            w=weights
        else
            w=1.0_dp
        end if
        if(sum(w)>0.0_dp)w=w/sum(w)
        allocate(m(n,kmax),aligned(n,kmax,size(memberships,3)),meanm(n,kmax))
        m=0.0_dp
        mb=consensus_dwh(memberships,w,kk)
        m(:,1:kk)=mb
        do b=1,size(memberships,3)
            if(use_manhattan)then
                call match_manhattan(m,memberships(:,:,b),mb)
            else
                call match_euclidean(m,memberships(:,:,b),mb)
            end if
            aligned(:,:,b)=mb
        end do
        if(use_manhattan)then
            oldv=objective_l1(m,aligned,w)
        else
            oldv=objective_e2(m,aligned,w)
        end if
        do it=1,mi
            meanm=0.0_dp
            do b=1,size(aligned,3)
                meanm=meanm+w(b)*aligned(:,:,b)
            end do
            m=0.0_dp
            do i=1,n
                jbest=maxloc(meanm(i,1:kk),dim=1)
                m(i,jbest)=1.0_dp
            end do
            do b=1,size(memberships,3)
                if(use_manhattan)then
                    call match_manhattan(m,memberships(:,:,b),mb)
                else
                    call match_euclidean(m,memberships(:,:,b),mb)
                end if
                aligned(:,:,b)=mb
            end do
            if(use_manhattan)then
                newv=objective_l1(m,aligned,w)
            else
                newv=objective_e2(m,aligned,w)
            end if
            if(abs(oldv-newv)<tol*(abs(oldv)+tol))exit
            oldv=newv
        end do
        ids=class_ids_from_membership(m(:,1:kk))
    end subroutine

    function consensus_hierarchy_euclidean(ultrametrics,weights) result(u)
        real(dp),intent(in)::ultrametrics(:,:,:)
        real(dp),intent(in),optional::weights(:)
        real(dp),allocatable::u(:,:),w(:)
        integer::b,n
        n=size(ultrametrics,1)
        allocate(u(n,n),w(size(ultrametrics,3)))
        if(present(weights))then
        w=weights
        else
        w=1
        end if
        if(abs(sum(w))>tiny(1.0_dp))w=w/sum(w)
        u=0
        do b=1,size(ultrametrics,3)
        u=u+w(b)*ultrametrics(:,:,b)
        end do
        call fit_ultrametric_ip(u,maxiter=1000,tol=1e-10_dp)
    end function

    integer function ensemble_medoid_index(dissimilarities) result(id)
        real(dp),intent(in)::dissimilarities(:,:)
        id=minloc(sum(dissimilarities,dim=2),dim=1)
    end function
end module clue_consensus
