module actuar_ruin_v02
    use actuar_kinds, only: dp
    use actuar_phase_type, only: pphtype
    implicit none
    private
    public :: ruin_result_t, ruin_phase_type, make_exponential_phase, make_erlang_phase

    type :: ruin_result_t
        real(dp), allocatable :: prob(:)
        real(dp), allocatable :: rates(:,:)
        logical :: converged = .false.
        integer :: iterations = 0
    contains
        procedure :: probability => ruin_probability
    end type ruin_result_t
contains
    pure real(dp) function ruin_probability(self,u) result(psi)
        class(ruin_result_t),intent(in)::self
        real(dp),intent(in)::u
        if(.not.allocated(self%prob) .or. .not.allocated(self%rates)) then
            psi=0.0_dp
        else
            psi=max(0.0_dp,min(1.0_dp,1.0_dp-pphtype(max(0.0_dp,u),self%prob,self%rates)))
        end if
    end function ruin_probability

    subroutine make_exponential_phase(rates,weights,prob,tmat)
        real(dp),intent(in)::rates(:)
        real(dp),intent(in),optional::weights(:)
        real(dp),allocatable,intent(out)::prob(:),tmat(:,:)
        integer::i,n
        n=size(rates);allocate(prob(n),tmat(n,n));tmat=0.0_dp
        if(present(weights)) then
            prob=weights/sum(weights)
        else if(n==1) then
            prob=1.0_dp
        else
            prob=1.0_dp/real(n,dp)
        end if
        do i=1,n;tmat(i,i)=-rates(i);end do
    end subroutine make_exponential_phase

    subroutine make_erlang_phase(shapes,rates,weights,prob,tmat)
        integer,intent(in)::shapes(:)
        real(dp),intent(in)::rates(:)
        real(dp),intent(in),optional::weights(:)
        real(dp),allocatable,intent(out)::prob(:),tmat(:,:)
        integer::ncomp,nstate,c,i,j,first,last
        real(dp)::sw
        ncomp=size(shapes);nstate=sum(shapes)
        allocate(prob(nstate),tmat(nstate,nstate));prob=0.0_dp;tmat=0.0_dp
        sw=1.0_dp
        if(present(weights)) sw=sum(weights)
        first=1
        do c=1,ncomp
            last=first+shapes(c)-1
            if(present(weights)) then
                prob(first)=weights(c)/sw
            else if(ncomp==1) then
                prob(first)=1.0_dp
            else
                prob(first)=1.0_dp/real(ncomp,dp)
            end if
            do i=first,last
                tmat(i,i)=-rates(min(c,size(rates)))
            end do
            do j=first,last-1
                tmat(j,j+1)=rates(min(c,size(rates)))
            end do
            first=last+1
        end do
    end subroutine make_erlang_phase

    function ruin_phase_type(prob_claim,t_claim,prob_wait,t_wait,premium_rate,tol,maxit) result(res)
        real(dp),intent(in)::prob_claim(:),t_claim(:,:),prob_wait(:),t_wait(:,:),premium_rate
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(ruin_result_t)::res
        integer::n,m,i,j,it,mi
        real(dp)::eps,denom
        real(dp),allocatable::q(:,:),qold(:,:),rs_claim(:),rs_wait(:),rhs(:,:),sol(:,:), &
            identn(:,:),identm(:,:),a(:,:),b(:,:),c(:,:),kq(:,:),matsolve(:,:),t0pi(:,:),tmp(:,:),y(:)
        logical::ok
        n=size(prob_claim);m=size(prob_wait)
        if(n==0 .or. m==0 .or. premium_rate<=0.0_dp) return
        if(size(t_claim,1)/=n .or. size(t_claim,2)/=n) return
        if(size(t_wait,1)/=m .or. size(t_wait,2)/=m) return
        eps=sqrt(epsilon(1.0_dp));if(present(tol)) eps=tol
        mi=200;if(present(maxit)) mi=maxit
        allocate(q(n,n),qold(n,n),rs_claim(n),rs_wait(m))
        rs_claim=sum(t_claim,dim=2);rs_wait=sum(t_wait,dim=2)
        q=t_claim
        if(m==1) then
            allocate(y(n),rhs(n,1),sol(n,1))
            rhs(:,1)=prob_claim
            call solve_matrix(transpose(t_claim),rhs,sol,ok)
            if(.not.ok) return
            y=sol(:,1)
            allocate(res%prob(n),res%rates(n,n))
            res%prob=t_wait(1,1)*y/premium_rate
            do i=1,n
                do j=1,n
                    q(i,j)=t_claim(i,j)-rs_claim(i)*res%prob(j)
                end do
            end do
            res%rates=q;res%converged=.true.;res%iterations=0
            return
        end if

        allocate(identn(n,n),identm(m,m));call identity_matrix(identn);call identity_matrix(identm)
        allocate(a(n,n*m),b(n*m,n*m),c(n*m,n),t0pi(n,n))
        a=kron(identn,reshape(prob_wait,[1,m]))
        b=kron(identn,t_wait)
        c=kron(identn,reshape(-rs_wait,[m,1]))
        do i=1,n
            do j=1,n
                t0pi(i,j)=-rs_claim(i)*prob_claim(j)
            end do
        end do
        do it=1,mi
            qold=q
            kq=kron(q,identm)
            matsolve=kq+b
            call solve_matrix(matsolve,c,sol,ok)
            if(.not.ok) exit
            tmp=matmul(a,sol)
            q=t_claim-matmul(t0pi,tmp)
            if(maxval(sum(abs(q-qold),dim=2))<eps) then
                res%converged=.true.;exit
            end if
        end do
        res%iterations=min(it,mi)
        denom=-sum(t_claim)*premium_rate
        if(abs(denom)<=tiny(1.0_dp)) return
        allocate(res%prob(n),res%rates(n,n));res%rates=q
        do j=1,n
            res%prob(j)=sum(q(:,j)-t_claim(:,j))/denom
        end do
        where(abs(res%prob)<1.0e-14_dp) res%prob=0.0_dp
    end function ruin_phase_type

    pure subroutine identity_matrix(a)
        real(dp),intent(out)::a(:,:)
        integer::i
        a=0.0_dp
        do i=1,min(size(a,1),size(a,2));a(i,i)=1.0_dp;end do
    end subroutine identity_matrix

    pure function kron(a,b) result(c)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp),allocatable::c(:,:)
        integer::i,j,ib,jb,br,bc
        br=size(b,1);bc=size(b,2)
        allocate(c(size(a,1)*br,size(a,2)*bc));c=0.0_dp
        do i=1,size(a,1)
            do j=1,size(a,2)
                do ib=1,br
                    do jb=1,bc
                        c((i-1)*br+ib,(j-1)*bc+jb)=a(i,j)*b(ib,jb)
                    end do
                end do
            end do
        end do
    end function kron

    pure subroutine solve_matrix(a,b,x,ok)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp),allocatable,intent(out)::x(:,:)
        logical,intent(out)::ok
        real(dp),allocatable::m(:,:),rhs(:,:),row(:),rrow(:)
        real(dp)::fac,piv
        integer::n,nrhs,i,k,imax
        n=size(a,1);nrhs=size(b,2)
        allocate(m(n,n),rhs(n,nrhs),row(n),rrow(nrhs),x(n,nrhs));m=a;rhs=b;ok=.true.
        do k=1,n
            imax=k
            do i=k+1,n;if(abs(m(i,k))>abs(m(imax,k))) imax=i;end do
            if(abs(m(imax,k))<1.0e-13_dp) then;ok=.false.;x=0.0_dp;return;end if
            if(imax/=k) then
                row=m(k,:);m(k,:)=m(imax,:);m(imax,:)=row
                rrow=rhs(k,:);rhs(k,:)=rhs(imax,:);rhs(imax,:)=rrow
            end if
            do i=k+1,n
                fac=m(i,k)/m(k,k);m(i,k:n)=m(i,k:n)-fac*m(k,k:n);rhs(i,:)=rhs(i,:)-fac*rhs(k,:)
            end do
        end do
        x=0.0_dp
        do i=n,1,-1
            if(i<n) then
                x(i,:)=(rhs(i,:)-matmul(m(i,i+1:n),x(i+1:n,:)))/m(i,i)
            else
                x(i,:)=rhs(i,:)/m(i,i)
            end if
        end do
    end subroutine solve_matrix
end module actuar_ruin_v02
