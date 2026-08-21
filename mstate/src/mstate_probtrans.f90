module mstate_probtrans
    use mstate_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use mstate_types, only : transition_map, hazard_type, probtrans_type, relative_bootstrap_type
    implicit none
    private
    public :: probtrans, probtrans_bootstrap, expected_length_of_stay

contains

    subroutine probtrans(hz, tr, predt, result, direction, method, variance, info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        real(dp),intent(in)::predt
        type(probtrans_type),intent(out)::result
        character(len=*),intent(in),optional::direction,method
        logical,intent(in),optional::variance
        integer,intent(out),optional::info
        character(len=16)::dir,meth
        logical::dovar
        integer::s,k,ntsel,i,j,r,idx,n2,oldidx
        integer,allocatable::sel(:)
        real(dp),allocatable::p(:,:),ia(:,:),dah(:),b(:,:),dc(:,:),varda(:,:),varp(:,:),tmp1(:,:),tmp2(:,:)
        real(dp),allocatable::left(:,:),right(:,:),id(:,:),sevec(:),cprev(:,:),ccur(:,:)

        if(present(info)) info=0
        dir='forward'; if(present(direction)) dir=trim(adjustl(direction))
        meth='aalen'; if(present(method)) meth=trim(adjustl(method))
        dovar=.true.; if(present(variance)) dovar=variance
        s=tr%nstate; k=tr%ntrans; n2=s*s
        if(hz%ntrans/=k) then; if(present(info))info=1; return; end if
        if(.not.allocated(hz%haz)) then; if(present(info))info=2; return; end if
        if(dovar .and. .not.allocated(hz%varhaz)) dovar=.false.
        if(dir=='forward') then
            ntsel=count(hz%time>predt); allocate(sel(ntsel)); j=0
            do i=1,hz%nt; if(hz%time(i)>predt) then; j=j+1;sel(j)=i;end if;end do
            result%nt=ntsel+1
        else
            ntsel=count(hz%time<=predt); allocate(sel(ntsel)); j=0
            do i=1,hz%nt; if(hz%time(i)<=predt) then;j=j+1;sel(j)=i;end if;end do
            result%nt=ntsel+1
        end if
        result%nstate=s
        allocate(result%time(result%nt),result%p(result%nt,s,s),result%se(result%nt,s,s)); result%se=0.0_dp
        allocate(p(s,s),ia(s,s),id(s,s),dah(k),b(n2,k),dc(k,k),varda(n2,n2),varp(n2,n2), &
                 tmp1(n2,n2),tmp2(n2,n2),left(n2,n2),right(n2,n2),sevec(n2),cprev(k,k),ccur(k,k))
        id=0.0_dp; do i=1,s;id(i,i)=1.0_dp;end do
        b=0.0_dp
        do r=1,k
            b(tr%from(r)+(tr%to(r)-1)*s,r)=1.0_dp
            b(tr%from(r)+(tr%from(r)-1)*s,r)=-1.0_dp
        end do
        p=id; varp=0.0_dp
        if(dir=='forward') then
            result%time(1)=predt; result%p(1,:,:)=p
            cprev=0.0_dp
            oldidx=0
            do i=1,hz%nt
                if(hz%time(i)<=predt) then
                    if(dovar) cprev=hz%varhaz(i,:,:)
                    oldidx=i
                end if
            end do
            do j=1,ntsel
                idx=sel(j); ia=id
                do r=1,k
                    if(idx==1) then; dah(r)=hz%haz(idx,r); else; dah(r)=hz%haz(idx,r)-hz%haz(idx-1,r); end if
                    if(oldidx>0 .and. idx==oldidx+1) dah(r)=hz%haz(idx,r)-hz%haz(oldidx,r)
                    ia(tr%from(r),tr%to(r))=ia(tr%from(r),tr%to(r))+dah(r)
                    ia(tr%from(r),tr%from(r))=ia(tr%from(r),tr%from(r))-dah(r)
                end do
                if(dovar) then
                    ccur=hz%varhaz(idx,:,:); dc=ccur-cprev; cprev=ccur
                    varda=matmul(b,matmul(dc,transpose(b)))
                    if(meth=='aalen') then
                        p=matmul(p,ia)
                        call kron(transpose(ia),id,left); call kron(ia,id,right)
                        tmp1=matmul(left,matmul(varp,right))
                        call kron(id,p,left); call kron(id,transpose(p),right)
                        tmp2=matmul(left,matmul(varda,right)); varp=tmp1+tmp2
                    else
                        call kron(transpose(ia),id,left); call kron(ia,id,right)
                        tmp1=matmul(left,matmul(varp,right))
                        call kron(id,p,left); call kron(id,transpose(p),right)
                        tmp2=matmul(left,matmul(varda,right)); varp=tmp1+tmp2
                        p=matmul(p,ia)
                    end if
                    do r=1,n2; sevec(r)=sqrt(max(0.0_dp,varp(r,r))); end do
                    result%se(j+1,:,:)=reshape(sevec,[s,s])
                else
                    p=matmul(p,ia)
                end if
                result%time(j+1)=hz%time(idx); result%p(j+1,:,:)=p
                oldidx=idx
            end do
        else
            result%time(result%nt)=predt; result%p(result%nt,:,:)=p
            ccur=0.0_dp
            if(ntsel>0.and.dovar) ccur=hz%varhaz(sel(ntsel),:,:)
            do j=ntsel,1,-1
                idx=sel(j); ia=id
                do r=1,k
                    if(idx==1) then; dah(r)=hz%haz(idx,r); else; dah(r)=hz%haz(idx,r)-hz%haz(idx-1,r); end if
                    ia(tr%from(r),tr%to(r))=ia(tr%from(r),tr%to(r))+dah(r)
                    ia(tr%from(r),tr%from(r))=ia(tr%from(r),tr%from(r))-dah(r)
                end do
                if(dovar) then
                    if(j>1) then;cprev=hz%varhaz(sel(j-1),:,:);else;cprev=0.0_dp;end if
                    dc=ccur-cprev; ccur=cprev; varda=matmul(b,matmul(dc,transpose(b)))
                    if(meth=='aalen') then
                        call kron(id,ia,left); call kron(id,transpose(ia),right)
                        tmp1=matmul(left,matmul(varp,right))
                        call kron(transpose(p),ia,left); call kron(p,transpose(ia),right)
                        tmp2=matmul(left,matmul(varda,right)); varp=tmp1+tmp2
                    else
                        call kron(id,ia,left); call kron(id,transpose(ia),right)
                        tmp1=matmul(left,matmul(varp,right))
                        call kron(transpose(p),id,left); call kron(p,id,right)
                        tmp2=matmul(left,matmul(varda,right)); varp=tmp1+tmp2
                    end if
                    p=matmul(ia,p)
                    do r=1,n2;sevec(r)=sqrt(max(0.0_dp,varp(r,r)));end do
                    result%se(j,:,:)=reshape(sevec,[s,s])
                else
                    p=matmul(ia,p)
                end if
                result%time(j)=hz%time(idx); result%p(j,:,:)=p
            end do
        end if
    end subroutine probtrans


    subroutine probtrans_bootstrap(hz,tr,boot,predt,result,direction,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        type(relative_bootstrap_type),intent(in)::boot
        real(dp),intent(in)::predt
        type(probtrans_type),intent(out)::result
        character(len=*),intent(in),optional::direction
        integer,intent(out),optional::info
        type(hazard_type)::bhz
        type(probtrans_type)::bpt
        real(dp),allocatable::psum(:,:,:),psq(:,:,:)
        character(len=16)::dir
        integer::ib,i,j,k,nv,ierr

        if(present(info))info=0
        dir='forward';if(present(direction))dir=trim(adjustl(direction))
        if(boot%nt/=hz%nt.or.boot%ntrans/=tr%ntrans.or.boot%b<1)then
            if(present(info))info=1
            return
        end if
        call probtrans(hz,tr,predt,result,direction=dir,variance=.false.,info=ierr)
        if(ierr/=0)then;if(present(info))info=10+ierr;return;end if
        allocate(psum(result%nt,tr%nstate,tr%nstate),psq(result%nt,tr%nstate,tr%nstate))
        psum=0.0_dp;psq=0.0_dp;nv=0
        bhz%nt=hz%nt;bhz%ntrans=hz%ntrans
        allocate(bhz%time(hz%nt),bhz%haz(hz%nt,hz%ntrans),bhz%varhaz(hz%nt,hz%ntrans,hz%ntrans))
        bhz%time=hz%time;bhz%varhaz=0.0_dp
        do ib=1,boot%b
            if(allocated(boot%valid_rep))then
                if(.not.boot%valid_rep(ib))cycle
            end if
            bhz%haz=0.0_dp
            do i=1,hz%nt
                do j=1,hz%ntrans
                    if(ieee_is_finite(boot%haz(i,j,ib)))bhz%haz(i,j)=boot%haz(i,j,ib)
                end do
            end do
            call probtrans(bhz,tr,predt,bpt,direction=dir,variance=.false.,info=ierr)
            if(ierr/=0.or.bpt%nt/=result%nt)cycle
            nv=nv+1;psum=psum+bpt%p;psq=psq+bpt%p*bpt%p
        end do
        if(nv>=2)then
            do i=1,result%nt
                do j=1,tr%nstate
                    do k=1,tr%nstate
                        result%se(i,j,k)=sqrt(max(0.0_dp,(psq(i,j,k)-psum(i,j,k)**2/real(nv,dp))/real(nv-1,dp)))
                    end do
                end do
            end do
        else if(present(info))then
            info=2
        end if
    end subroutine probtrans_bootstrap

    subroutine kron(a,b,c)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp),intent(out)::c(:,:)
        integer::i,j,m,n
        m=size(b,1);n=size(b,2)
        do i=1,size(a,1);do j=1,size(a,2)
            c((i-1)*m+1:i*m,(j-1)*n+1:j*n)=a(i,j)*b
        end do;end do
    end subroutine kron

    subroutine expected_length_of_stay(pt,tau,elos,info)
        type(probtrans_type),intent(in)::pt
        real(dp),intent(in)::tau
        real(dp),allocatable,intent(out)::elos(:,:)
        integer,intent(out),optional::info
        integer::s,g,h,i
        real(dp)::t1,t2
        if(present(info))info=0
        s=pt%nstate; allocate(elos(s,s)); elos=0.0_dp
        if(pt%nt<1 .or. tau<pt%time(1)) then;if(present(info))info=1;return;end if
        do i=1,pt%nt
            t1=pt%time(i); if(t1>=tau) exit
            if(i<pt%nt) then;t2=min(tau,pt%time(i+1));else;t2=tau;end if
            if(t2<=t1) cycle
            do g=1,s;do h=1,s;elos(g,h)=elos(g,h)+(t2-t1)*pt%p(i,g,h);end do;end do
        end do
    end subroutine expected_length_of_stay
end module mstate_probtrans
