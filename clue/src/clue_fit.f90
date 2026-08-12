! SPDX-License-Identifier: GPL-2.0-only
module clue_fit
    use clue_kinds, only: dp
    use clue_sumt, only: sumt_result, sumt_optimize
    use clue_trees, only: non_ultrametricity, ultrametricity_gradient, &
        non_additivity, additivity_gradient, ultrametrify, fit_ultrametric_ip, &
        fit_ultrametric_ir, fit_addtree_ip, fit_addtree_ir
    implicit none
    private

    type, public :: tree_fit_result
        real(dp), allocatable :: distance(:,:)
        real(dp) :: objective = huge(1.0_dp)
        integer :: status = 0
        integer :: iterations = 0
    end type tree_fit_result

    public :: fit_ultrametric_sumt, fit_addtree_sumt
    public :: fit_l1_ultrametric_sumt, fit_l1_ultrametric_irip
    public :: fit_sum_of_ultrametrics

contains

    function fit_ultrametric_sumt(d,weights,eps,q,max_outer,max_inner) result(out)
        real(dp), intent(in) :: d(:,:)
        real(dp), intent(in), optional :: weights(:,:)
        real(dp), intent(in), optional :: eps,q
        integer, intent(in), optional :: max_outer,max_inner
        type(tree_fit_result) :: out
        real(dp), allocatable :: x(:), w(:), fitted(:,:)
        type(sumt_result) :: sr
        integer :: n

        n=size(d,1)
        call matrix_to_veclh(d,x)
        call make_pair_weights(n,weights,w)
        w=w/sum(w)
        sr=sumt_optimize(x,loss,penalty,grad_loss,grad_penalty,eps,q,max_outer,max_inner)
        call veclh_to_matrix(sr%x,n,fitted)
        fitted=ultrametrify(fitted)
        allocate(out%distance(size(fitted,1),size(fitted,2)))
        out%distance=fitted
        out%objective=loss_matrix(fitted,d,w)
        out%status=merge(1,0,sr%penalty>1.0e-8_dp)
        out%iterations=sr%outer_iterations
    contains
        function loss(z) result(v)
            real(dp),intent(in)::z(:)
            real(dp)::v
            v=sum(w*(z-x)**2)
        end function
        function penalty(z) result(v)
            real(dp),intent(in)::z(:)
            real(dp)::v
            real(dp),allocatable::m(:,:)
            call veclh_to_matrix(z,n,m)
            v=non_ultrametricity(m)+sum(min(z,0.0_dp)**2)
        end function
        subroutine grad_loss(z,g)
            real(dp),intent(in)::z(:)
            real(dp),intent(out)::g(:)
            g=2.0_dp*w*(z-x)
        end subroutine
        subroutine grad_penalty(z,g)
            real(dp),intent(in)::z(:)
            real(dp),intent(out)::g(:)
            real(dp),allocatable::m(:,:),gm(:,:),gv(:)
            call veclh_to_matrix(z,n,m)
            gm=ultrametricity_gradient(m)
            call matrix_to_veclh(gm,gv)
            g=gv
            ! Preserve clue's R-level penalty-gradient expression.
            g=g+2.0_dp*sum(min(z,0.0_dp))
        end subroutine
    end function fit_ultrametric_sumt

    function fit_addtree_sumt(d,weights,eps,q,max_outer,max_inner) result(out)
        real(dp), intent(in) :: d(:,:)
        real(dp), intent(in), optional :: weights(:,:)
        real(dp), intent(in), optional :: eps,q
        integer, intent(in), optional :: max_outer,max_inner
        type(tree_fit_result) :: out
        real(dp), allocatable :: x(:), w(:), fitted(:,:)
        type(sumt_result) :: sr
        integer :: n

        n=size(d,1)
        call matrix_to_veclh(d,x)
        call make_pair_weights(n,weights,w)
        w=w/sum(w)
        sr=sumt_optimize(x,loss,penalty,grad_loss,grad_penalty,eps,q,max_outer,max_inner)
        call veclh_to_matrix(sr%x,n,fitted)
        fitted=max(fitted,0.0_dp)
        allocate(out%distance(size(fitted,1),size(fitted,2)))
        out%distance=fitted
        out%objective=loss_matrix(fitted,d,w)
        out%status=merge(1,0,sr%penalty>1.0e-8_dp)
        out%iterations=sr%outer_iterations
    contains
        function loss(z) result(v)
            real(dp),intent(in)::z(:)
            real(dp)::v
            v=sum(w*(z-x)**2)
        end function
        function penalty(z) result(v)
            real(dp),intent(in)::z(:)
            real(dp)::v
            real(dp),allocatable::m(:,:)
            call veclh_to_matrix(z,n,m)
            v=non_additivity(m)+sum(min(z,0.0_dp)**2)
        end function
        subroutine grad_loss(z,g)
            real(dp),intent(in)::z(:)
            real(dp),intent(out)::g(:)
            g=2.0_dp*w*(z-x)
        end subroutine
        subroutine grad_penalty(z,g)
            real(dp),intent(in)::z(:)
            real(dp),intent(out)::g(:)
            real(dp),allocatable::m(:,:),gm(:,:),gv(:)
            call veclh_to_matrix(z,n,m)
            gm=additivity_gradient(m)
            call matrix_to_veclh(gm,gv)
            g=gv
            ! Preserve clue's R-level penalty-gradient expression.
            g=g+2.0_dp*sum(min(z,0.0_dp))
        end subroutine
    end function fit_addtree_sumt

    function fit_l1_ultrametric_sumt(d,weights,eps,q,max_outer,max_inner) result(out)
        real(dp), intent(in) :: d(:,:)
        real(dp), intent(in), optional :: weights(:,:)
        real(dp), intent(in), optional :: eps,q
        integer, intent(in), optional :: max_outer,max_inner
        type(tree_fit_result) :: out
        real(dp), allocatable :: x(:),w(:),fitted(:,:)
        type(sumt_result) :: sr
        integer :: n

        n=size(d,1)
        call matrix_to_veclh(d,x)
        call make_pair_weights(n,weights,w)
        w=w/sum(w)
        sr=sumt_optimize(x,loss,penalty,grad_loss,grad_penalty,eps,q,max_outer,max_inner)
        call veclh_to_matrix(sr%x,n,fitted)
        fitted=ultrametrify(fitted)
        allocate(out%distance(size(fitted,1),size(fitted,2)))
        out%distance=fitted
        out%objective=l1_loss_matrix(fitted,d,w)
        out%status=merge(1,0,sr%penalty>1.0e-8_dp)
        out%iterations=sr%outer_iterations
    contains
        function loss(z) result(v)
            real(dp),intent(in)::z(:)
            real(dp)::v
            v=sum(w*abs(z-x))
        end function
        function penalty(z) result(v)
            real(dp),intent(in)::z(:)
            real(dp)::v
            real(dp),allocatable::m(:,:)
            call veclh_to_matrix(z,n,m)
            v=non_ultrametricity(m)+sum(min(z,0.0_dp)**2)
        end function
        subroutine grad_loss(z,g)
            real(dp),intent(in)::z(:)
            real(dp),intent(out)::g(:)
            g=w*sign(1.0_dp,z-x)
            where(abs(z-x)<=epsilon(1.0_dp)) g=0.0_dp
        end subroutine
        subroutine grad_penalty(z,g)
            real(dp),intent(in)::z(:)
            real(dp),intent(out)::g(:)
            real(dp),allocatable::m(:,:),gm(:,:),gv(:)
            call veclh_to_matrix(z,n,m)
            gm=ultrametricity_gradient(m)
            call matrix_to_veclh(gm,gv)
            g=gv
            g=g+2.0_dp*sum(min(z,0.0_dp))
        end subroutine
    end function fit_l1_ultrametric_sumt

    function fit_l1_ultrametric_irip(d,weights,min_delta,eps,maxiter,reltol) result(out)
        real(dp), intent(in) :: d(:,:)
        real(dp), intent(in), optional :: weights(:,:)
        real(dp), intent(in), optional :: min_delta,eps,reltol
        integer, intent(in), optional :: maxiter
        type(tree_fit_result) :: out
        real(dp),allocatable::x(:),w(:),u(:),rw(:),umat(:,:),wm(:,:)
        real(dp)::mind,tol,rt,lold,lnew,du,dl
        integer::n,it,mi
        type(tree_fit_result)::ls

        n=size(d,1)
        call matrix_to_veclh(d,x)
        call make_pair_weights(n,weights,w)
        w=w/sum(w)
        mind=1.0e-3_dp
        if(present(min_delta))mind=min_delta
        tol=1.0e-6_dp
        if(present(eps))tol=eps
        rt=1.0e-6_dp
        if(present(reltol))rt=reltol
        mi=100
        if(present(maxiter))mi=maxiter
        u=x
        lnew=sum(w*abs(x-u))
        do it=1,mi
            lold=lnew
            rw=w/max(abs(u-x),mind)
            call vector_weights_to_matrix(rw,n,wm)
            call veclh_to_matrix(x,n,umat)
            ls=fit_ultrametric_sumt(umat,wm,max_outer=50,max_inner=300)
            call matrix_to_veclh(ls%distance,rw)
            du=maxval(abs(u-rw))
            u=rw
            lnew=sum(w*abs(x-u))
            dl=lold-lnew
            if(du<tol)exit
            if(dl>=0.0_dp .and. dl<=rt*(abs(lold)+rt))exit
        end do
        call veclh_to_matrix(u,n,umat)
        umat=ultrametrify(umat)
        allocate(out%distance(size(umat,1),size(umat,2)))
        out%distance=umat
        out%objective=lnew
        out%iterations=it
        out%status=merge(1,0,it>=mi)
    end function fit_l1_ultrametric_irip

    function fit_sum_of_ultrametrics(d,nterms,method,eps,maxiter,reltol) result(terms)
        real(dp), intent(in) :: d(:,:)
        integer, intent(in) :: nterms
        character(*), intent(in), optional :: method
        real(dp), intent(in), optional :: eps,reltol
        integer, intent(in), optional :: maxiter
        real(dp), allocatable :: terms(:,:,:)
        real(dp),allocatable::resid(:,:),old(:,:),total(:,:)
        real(dp)::tol,rt,lold,lnew,du,dl
        integer::n,nt,mi,it,i,j
        character(8)::meth
        type(tree_fit_result)::fr

        n=size(d,1)
        nt=max(1,nterms)
        allocate(terms(n,n,nt),resid(n,n),old(n,n),total(n,n))
        terms=0.0_dp
        tol=1.0e-6_dp
        if(present(eps))tol=eps
        rt=1.0e-6_dp
        if(present(reltol))rt=reltol
        mi=100
        if(present(maxiter))mi=maxiter
        meth='SUMT'
        if(present(method))meth=adjustl(method)
        lnew=sum(lower_squared(d))
        do it=1,mi
            lold=lnew
            du=0.0_dp
            do i=1,nt
                total=0.0_dp
                do j=1,nt
                    if(j/=i)total=total+terms(:,:,j)
                end do
                resid=d-total
                old=terms(:,:,i)
                select case(trim(meth))
                case('IP','ip')
                    terms(:,:,i)=resid
                    call fit_ultrametric_ip(terms(:,:,i),maxiter=10000,tol=1.0e-8_dp)
                    terms(:,:,i)=max(terms(:,:,i),0.0_dp)
                case('IR','ir')
                    terms(:,:,i)=resid
                    call fit_ultrametric_ir(terms(:,:,i),maxiter=10000,tol=1.0e-8_dp)
                    terms(:,:,i)=max(terms(:,:,i),0.0_dp)
                case default
                    fr=fit_ultrametric_sumt(resid,max_outer=50,max_inner=300)
                    terms(:,:,i)=fr%distance
                end select
                du=max(du,maxval(abs(terms(:,:,i)-old)))
            end do
            total=sum(terms,dim=3)
            lnew=sum(lower_squared(d-total))
            dl=lold-lnew
            if(du<tol)exit
            if(dl>=0.0_dp .and. dl<=rt*(abs(lold)+rt))exit
        end do
    end function fit_sum_of_ultrametrics

    subroutine matrix_to_veclh(a,v)
        real(dp),intent(in)::a(:,:)
        real(dp),allocatable,intent(out)::v(:)
        integer::n,i,j,p
        n=size(a,1)
        allocate(v(n*(n-1)/2))
        p=0
        do j=1,n-1
            do i=j+1,n
                p=p+1
                v(p)=a(i,j)
            end do
        end do
    end subroutine matrix_to_veclh

    subroutine veclh_to_matrix(v,n,a)
        real(dp),intent(in)::v(:)
        integer,intent(in)::n
        real(dp),allocatable,intent(out)::a(:,:)
        integer::i,j,p
        allocate(a(n,n))
        a=0.0_dp
        p=0
        do j=1,n-1
            do i=j+1,n
                p=p+1
                a(i,j)=v(p)
                a(j,i)=v(p)
            end do
        end do
    end subroutine veclh_to_matrix

    subroutine make_pair_weights(n,wm,w)
        integer,intent(in)::n
        real(dp),intent(in),optional::wm(:,:)
        real(dp),allocatable,intent(out)::w(:)
        integer::i,j,p
        allocate(w(n*(n-1)/2))
        w=1.0_dp
        if(.not.present(wm))return
        p=0
        do j=1,n-1
            do i=j+1,n
                p=p+1
                w(p)=wm(i,j)
            end do
        end do
        if(.not.any(w>0.0_dp))w=1.0_dp
    end subroutine make_pair_weights

    subroutine vector_weights_to_matrix(w,n,wm)
        real(dp),intent(in)::w(:)
        integer,intent(in)::n
        real(dp),allocatable,intent(out)::wm(:,:)
        integer::i,j,p
        allocate(wm(n,n))
        wm=0.0_dp
        p=0
        do j=1,n-1
            do i=j+1,n
                p=p+1
                wm(i,j)=w(p)
                wm(j,i)=w(p)
            end do
        end do
    end subroutine vector_weights_to_matrix

    function loss_matrix(a,b,w) result(v)
        real(dp),intent(in)::a(:,:),b(:,:),w(:)
        real(dp)::v
        real(dp),allocatable::av(:),bv(:)
        call matrix_to_veclh(a,av)
        call matrix_to_veclh(b,bv)
        v=sum(w*(av-bv)**2)
    end function loss_matrix

    function l1_loss_matrix(a,b,w) result(v)
        real(dp),intent(in)::a(:,:),b(:,:),w(:)
        real(dp)::v
        real(dp),allocatable::av(:),bv(:)
        call matrix_to_veclh(a,av)
        call matrix_to_veclh(b,bv)
        v=sum(w*abs(av-bv))
    end function l1_loss_matrix

    function lower_squared(a) result(v)
        real(dp),intent(in)::a(:,:)
        real(dp),allocatable::v(:)
        call matrix_to_veclh(a,v)
        v=v*v
    end function lower_squared

end module clue_fit
