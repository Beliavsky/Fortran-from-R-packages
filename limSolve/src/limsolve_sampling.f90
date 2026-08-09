! Upstream license declaration: GPL (version unspecified)
module limsolve_sampling
    use limsolve_kinds, only: dp
    use limsolve_types
    use limsolve_linalg, only: null_space, identity_matrix
    use limsolve_inverse, only: lsei
    use limsolve_ranges, only: linp
    implicit none
    private
    public :: xsample, seed_rng

contains

    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n,i
        integer, allocatable :: s(:)
        call random_seed(size=n); allocate(s(n))
        do i=1,n
            s(i)=mod(abs(seed)+104729*i,2147483646)+1
        end do
        call random_seed(put=s)
    end subroutine seed_rng

    subroutine xsample(a,b,e,f,g,h,result,iter,outputlength,burnin,type,jump,tol,xstart, &
        sd_b,weights,ispos,lower,upper,seed,fulloutput)
        real(dp), intent(in) :: a(:,:),b(:),e(:,:),f(:),g(:,:),h(:)
        type(sample_result), intent(out) :: result
        integer, intent(in), optional :: iter,outputlength,burnin,seed
        character(len=*), intent(in), optional :: type
        real(dp), intent(in), optional :: jump(:),tol,xstart(:),sd_b(:),weights(:)
        logical, intent(in), optional :: ispos,fulloutput
        real(dp), intent(in), optional :: lower(:),upper(:)
        real(dp), allocatable :: gg(:,:),hh(:),z(:,:),x0(:),gq(:,:),hq(:),aq(:,:),bq(:)
        real(dp), allocatable :: q1(:),q2(:),jmp(:),xcur(:),aw(:,:),bw(:)
        type(solve_result) :: start,lpstart
        real(dp) :: eps,p1,p2,u
        integer :: n,k,niter,nout,nburn,ou,i,ii,accepted,me,ma,st,ranke
        logical :: pos,wantq
        character(len=12) :: method

        n=max(size(a,2),max(size(e,2),size(g,2))); me=size(e,1); ma=size(a,1)
        result%status=LS_INVALID
        if(n<=0 .or. size(a,2)/=n .or. size(e,2)/=n .or. size(g,2)/=n) return
        if(size(b)/=ma .or. size(f)/=me .or. size(h)/=size(g,1)) return
        eps=sqrt(epsilon(1.0_dp)); if(present(tol)) eps=tol
        niter=3000; if(present(iter)) niter=max(1,iter)
        nout=niter; if(present(outputlength)) nout=max(1,min(outputlength,niter))
        nburn=0; if(present(burnin)) nburn=max(0,burnin)
        method='mirror'; if(present(type)) method=adjustl(type)
        pos=.false.; if(present(ispos)) pos=ispos
        wantq=.false.; if(present(fulloutput)) wantq=fulloutput
        if(present(seed)) call seed_rng(seed)
        call augment_bounds_local(g,h,n,gg,hh,lower,upper)

        allocate(x0(n))
        if(present(xstart)) then
            if(size(xstart)/=n) return
            x0=xstart
        else
            call lsei_bounds_dispatch(a,b,e,f,g,h,start,eps,lower,upper)
            if(start%status==LS_SUCCESS .and. start%residual_norm<=1.0e-6_dp) then
                x0=start%x
            else
                call linp_bounds_dispatch(e,f,g,h,[(1.0_dp,i=1,n)],lpstart,pos,lower,upper)
                if(lpstart%status/=LS_SUCCESS .or. lpstart%residual_norm>1.0e-6_dp) then
                    result%status=LS_INFEASIBLE; return
                end if
                x0=lpstart%x
            end if
        end if
        if(me>0) then
            call null_space(e,z,ranke,eps,st)
        else
            z=identity_matrix(n)
        end if
        k=size(z,2)
        if(k==0) then
            allocate(result%x(1,n),result%p(1),result%jump(0))
            result%x(1,:)=x0; result%p=1.0_dp; result%accepted_ratio=1.0_dp
            if(wantq) allocate(result%q(1,0))
            result%status=LS_SUCCESS; return
        end if
        allocate(gq(size(gg,1),k),hq(size(gg,1)),aq(ma,k),bq(ma))
        if(size(gg,1)>0) then
            gq=matmul(gg,z); hq=hh-matmul(gg,x0)
        end if
        if(ma>0) then
            aq=matmul(a,z); bq=b-matmul(a,x0)
            allocate(aw(ma,k),bw(ma)); aw=aq; bw=bq
            if(present(weights)) then
                if(size(weights)==1) then
                    aw=aw*weights(1); bw=bw*weights(1)
                else if(size(weights)==ma) then
                    do i=1,ma; aw(i,:)=aw(i,:)*weights(i); bw(i)=bw(i)*weights(i); end do
                end if
            end if
            if(present(sd_b)) then
                if(size(sd_b)==1) then
                    aw=aw/sd_b(1); bw=bw/sd_b(1)
                else if(size(sd_b)==ma) then
                    do i=1,ma; aw(i,:)=aw(i,:)/sd_b(i); bw(i)=bw(i)/sd_b(i); end do
                end if
            end if
            aq=aw; bq=bw
        end if
        allocate(jmp(k)); jmp=1.0_dp
        if(present(jump)) then
            if(size(jump)==1) then; jmp=jump(1)
            else if(size(jump)==k) then; jmp=jump
            end if
        else if(size(gq,1)>0) then
            call estimate_jump(gq,hq,jmp)
        end if
        allocate(q1(k),q2(k),xcur(n)); q1=0.0_dp
        p1=probability(q1,aq,bq)
        do i=1,nburn
            call propose(method,q1,gq,hq,jmp,q2)
            p2=probability(q2,aq,bq); call random_number(u)
            if(p2>=p1 .or. u<min(1.0_dp,p2/max(p1,tiny(1.0_dp)))) then; q1=q2; p1=p2; end if
        end do
        allocate(result%x(nout,n),result%p(nout),result%jump(k))
        if(wantq) allocate(result%q(nout,k))
        result%jump=jmp; result%x=0.0_dp; result%p=0.0_dp
        result%x(1,:)=x0+matmul(z,q1); result%p(1)=p1
        if(wantq) result%q(1,:)=q1
        accepted=1; ou=ceiling(real(niter,dp)/real(nout,dp))
        do i=2,nout
            do ii=1,ou
                call propose(method,q1,gq,hq,jmp,q2)
                p2=probability(q2,aq,bq); call random_number(u)
                if(p2>=p1 .or. u<min(1.0_dp,p2/max(p1,tiny(1.0_dp)))) then
                    q1=q2; p1=p2; accepted=accepted+1
                end if
            end do
            xcur=x0+matmul(z,q1); result%x(i,:)=xcur; result%p(i)=p1
            if(wantq) result%q(i,:)=q1
        end do
        result%accepted_ratio=real(accepted,dp)/real(max(1,niter),dp)
        result%status=LS_SUCCESS
    end subroutine xsample

    subroutine propose(method,q,g,h,jmp,q2)
        character(len=*),intent(in)::method
        real(dp),intent(in)::q(:),g(:,:),h(:),jmp(:)
        real(dp),intent(out)::q2(:)
        select case(trim(method))
        case('rda'); call propose_rda(q,g,h,q2)
        case('cda'); call propose_cda(q,g,h,q2)
        case default; call propose_mirror(q,g,h,jmp,q2)
        end select
    end subroutine propose

    subroutine propose_mirror(q,g,h,jmp,q2)
        real(dp),intent(in)::q(:),g(:,:),h(:),jmp(:); real(dp),intent(out)::q2(:)
        real(dp),allocatable::q10(:),epsv(:),res(:),alpha(:)
        real(dp)::d,best,den; integer::i,hit,loops
        q2=q+jmp*randn_vec(size(q)); if(size(g,1)==0) return
        allocate(q10(size(q)),epsv(size(q)),res(size(g,1)),alpha(size(g,1)))
        q10=q; res=matmul(g,q2)-h; loops=0
        do while(any(res<0.0_dp) .and. loops<1000)
            loops=loops+1; epsv=q2-q10; alpha=huge(1.0_dp)
            do i=1,size(g,1)
                den=dot_product(g(i,:),epsv)
                if(res(i)<0.0_dp .and. abs(den)>tiny(1.0_dp)) alpha(i)=(h(i)-dot_product(g(i,:),q10))/den
            end do
            hit=0; best=huge(1.0_dp)
            do i=1,size(g,1)
                if(res(i)<0.0_dp .and. alpha(i)>=0.0_dp .and. alpha(i)<best) then; best=alpha(i); hit=i; end if
            end do
            if(hit==0) exit
            den=sum(g(hit,:)**2); if(den<=tiny(1.0_dp)) exit
            d=-res(hit)/den; q2=q2+2.0_dp*d*g(hit,:)
            q10=q10+best*epsv; res=matmul(g,q2)-h
        end do
    end subroutine propose_mirror

    subroutine propose_rda(q,g,h,q2)
        real(dp),intent(in)::q(:),g(:,:),h(:); real(dp),intent(out)::q2(:)
        real(dp),allocatable::d(:); real(dp)::lo,hi,den,a,u; integer::i
        d=randn_vec(size(q)); den=sqrt(sum(d*d)); if(den<=tiny(1.0_dp)) then; q2=q; return; end if
        d=d/den; lo=-huge(1.0_dp); hi=huge(1.0_dp)
        do i=1,size(g,1)
            den=dot_product(g(i,:),d); a=h(i)-dot_product(g(i,:),q)
            if(abs(den)<=tiny(1.0_dp)) then
                if(a>0.0_dp) then; q2=q; return; end if
            else if(den>0.0_dp) then; lo=max(lo,a/den)
            else; hi=min(hi,a/den); end if
        end do
        if(lo>hi .or. abs(lo)>=huge(1.0_dp)/10.0_dp .or. abs(hi)>=huge(1.0_dp)/10.0_dp) then; q2=q; return; end if
        call random_number(u); q2=q+(lo+u*(hi-lo))*d
    end subroutine propose_rda

    subroutine propose_cda(q,g,h,q2)
        real(dp),intent(in)::q(:),g(:,:),h(:); real(dp),intent(out)::q2(:)
        real(dp)::lo,hi,a,coef,u; integer::i,j
        call random_number(u); j=1+int(u*real(size(q),dp)); j=min(j,size(q)); q2=q
        lo=-huge(1.0_dp); hi=huge(1.0_dp)
        do i=1,size(g,1)
            coef=g(i,j); a=h(i)-dot_product(g(i,:),q)+coef*q(j)
            if(abs(coef)<=tiny(1.0_dp)) then
                if(a>0.0_dp) return
            else if(coef>0.0_dp) then; lo=max(lo,a/coef)
            else; hi=min(hi,a/coef); end if
        end do
        if(lo>hi .or. abs(lo)>=huge(1.0_dp)/10.0_dp .or. abs(hi)>=huge(1.0_dp)/10.0_dp) return
        call random_number(u); q2(j)=lo+u*(hi-lo)
    end subroutine propose_cda

    function probability(q,a,b) result(p)
        real(dp),intent(in)::q(:),a(:,:),b(:); real(dp)::p,ss
        if(size(a,1)==0) then; p=1.0_dp
        else; ss=sum((matmul(a,q)-b)**2); p=exp(-0.5_dp*min(ss,1400.0_dp)); end if
    end function probability

    function randn_vec(n) result(z)
        integer,intent(in)::n; real(dp),allocatable::z(:)
        real(dp)::u1,u2,rad,theta; integer::i
        allocate(z(n)); i=1
        do while(i<=n)
            call random_number(u1); call random_number(u2); u1=max(u1,tiny(1.0_dp))
            rad=sqrt(-2.0_dp*log(u1)); theta=2.0_dp*acos(-1.0_dp)*u2
            z(i)=rad*cos(theta); if(i+1<=n) z(i+1)=rad*sin(theta); i=i+2
        end do
    end function randn_vec

    subroutine estimate_jump(g,h,jmp)
        real(dp),intent(in)::g(:,:),h(:); real(dp),intent(out)::jmp(:)
        real(dp)::lo,hi,a,coef; integer::j,i
        do j=1,size(jmp)
            lo=-huge(1.0_dp); hi=huge(1.0_dp)
            do i=1,size(g,1)
                coef=g(i,j); a=h(i)
                if(abs(coef)<=tiny(1.0_dp)) cycle
                if(coef>0.0_dp) then; lo=max(lo,a/coef)
                else; hi=min(hi,a/coef); end if
            end do
            if(abs(lo)<huge(1.0_dp)/10.0_dp .and. abs(hi)<huge(1.0_dp)/10.0_dp .and. hi>lo) then
                jmp(j)=max((hi-lo)/5.0_dp,1.0e-6_dp)
            else; jmp(j)=1.0_dp; end if
        end do
    end subroutine estimate_jump

    subroutine lsei_bounds_dispatch(a,b,e,f,g,h,res,tol,lower,upper)
        real(dp),intent(in)::a(:,:),b(:),e(:,:),f(:),g(:,:),h(:),tol
        type(solve_result),intent(out)::res; real(dp),intent(in),optional::lower(:),upper(:)
        if(present(lower).and.present(upper)) then; call lsei(a,b,e,f,g,h,res,tol,lower=lower,upper=upper)
        else if(present(lower)) then; call lsei(a,b,e,f,g,h,res,tol,lower=lower)
        else if(present(upper)) then; call lsei(a,b,e,f,g,h,res,tol,upper=upper)
        else; call lsei(a,b,e,f,g,h,res,tol); end if
    end subroutine lsei_bounds_dispatch

    subroutine linp_bounds_dispatch(e,f,g,h,c,res,pos,lower,upper)
        real(dp),intent(in)::e(:,:),f(:),g(:,:),h(:),c(:); type(solve_result),intent(out)::res
        logical,intent(in)::pos; real(dp),intent(in),optional::lower(:),upper(:)
        if(present(lower).and.present(upper)) then; call linp(e,f,g,h,c,res,pos,lower=lower,upper=upper)
        else if(present(lower)) then; call linp(e,f,g,h,c,res,pos,lower=lower)
        else if(present(upper)) then; call linp(e,f,g,h,c,res,pos,upper=upper)
        else; call linp(e,f,g,h,c,res,pos); end if
    end subroutine linp_bounds_dispatch

    subroutine augment_bounds_local(g,h,n,gg,hh,lower,upper)
        real(dp),intent(in)::g(:,:),h(:); integer,intent(in)::n
        real(dp),allocatable,intent(out)::gg(:,:),hh(:); real(dp),intent(in),optional::lower(:),upper(:)
        integer::m,nl,nu,i,k; logical::scalar
        m=size(g,1); nl=0; nu=0
        if(present(lower)) then
            if(size(lower)==1) then; if(abs(lower(1))<huge(1.0_dp)/10.0_dp) nl=n
            else; nl=count(abs(lower)<huge(1.0_dp)/10.0_dp); end if
        end if
        if(present(upper)) then
            if(size(upper)==1) then; if(abs(upper(1))<huge(1.0_dp)/10.0_dp) nu=n
            else; nu=count(abs(upper)<huge(1.0_dp)/10.0_dp); end if
        end if
        allocate(gg(m+nl+nu,n),hh(m+nl+nu)); gg=0.0_dp; hh=0.0_dp; if(m>0) then; gg(1:m,:)=g; hh(1:m)=h; end if; k=m
        if (present(lower)) then
            scalar = size(lower) == 1
            do i = 1, n
                if (scalar) then
                    if (abs(lower(1)) >= huge(1.0_dp)/10.0_dp) cycle
                    k = k + 1
                    gg(k,i) = 1.0_dp
                    hh(k) = lower(1)
                else if (i <= size(lower)) then
                    if (abs(lower(i)) >= huge(1.0_dp)/10.0_dp) cycle
                    k = k + 1
                    gg(k,i) = 1.0_dp
                    hh(k) = lower(i)
                end if
            end do
        end if
        if (present(upper)) then
            scalar = size(upper) == 1
            do i = 1, n
                if (scalar) then
                    if (abs(upper(1)) >= huge(1.0_dp)/10.0_dp) cycle
                    k = k + 1
                    gg(k,i) = -1.0_dp
                    hh(k) = -upper(1)
                else if (i <= size(upper)) then
                    if (abs(upper(i)) >= huge(1.0_dp)/10.0_dp) cycle
                    k = k + 1
                    gg(k,i) = -1.0_dp
                    hh(k) = -upper(i)
                end if
            end do
        end if
    end subroutine augment_bounds_local

end module limsolve_sampling
