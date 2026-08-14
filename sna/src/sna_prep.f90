! Data preparation and graph-matrix utilities translated from R/sna.
! Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_prep
    use sna_kinds, only : dp, sna_nan, is_missing, sna_eps
    implicit none
    private

    public :: diag_remove, lower_tri_remove, upper_tri_remove, graph_transpose
    public :: symmetrize, event2dichot, make_stochastic, nties, gvectorize
    public :: interval_graph, add_isolates, ego_extract_mask, sr2css, stackcount
    public :: log_sum, log_mean, log_sub, sample_quantile

    interface stackcount
        module procedure stackcount_matrix, stackcount_stack
    end interface stackcount

contains


    pure integer function stackcount_matrix(d) result(nstack)
        real(dp), intent(in) :: d(:,:)
        nstack=1
    end function stackcount_matrix

    pure integer function stackcount_stack(d) result(nstack)
        real(dp), intent(in) :: d(:,:,:)
        nstack=size(d,1)
    end function stackcount_stack

    pure function diag_remove(dat, remove_val) result(out)
        real(dp), intent(in) :: dat(:,:)
        real(dp), intent(in), optional :: remove_val
        real(dp), allocatable :: out(:,:)
        real(dp) :: rv
        integer :: i
        rv=sna_nan()
        if(present(remove_val))rv=remove_val
        allocate(out(size(dat,1),size(dat,2)))
        out=dat
        do i=1,min(size(dat,1),size(dat,2))
        out(i,i)=rv
        end do
    end function diag_remove

    pure function lower_tri_remove(dat, remove_val) result(out)
        real(dp), intent(in)::dat(:,:)
        real(dp),intent(in),optional::remove_val
        real(dp),allocatable::out(:,:)
        real(dp)::rv
        integer::i,j
        rv=sna_nan()
        if(present(remove_val))rv=remove_val
        allocate(out(size(dat,1),size(dat,2)))
        out=dat
        do j=1,size(dat,2)
        do i=j+1,size(dat,1)
        out(i,j)=rv
        end do
        end do
    end function lower_tri_remove

    pure function upper_tri_remove(dat, remove_val) result(out)
        real(dp), intent(in)::dat(:,:)
        real(dp),intent(in),optional::remove_val
        real(dp),allocatable::out(:,:)
        real(dp)::rv
        integer::i,j
        rv=sna_nan()
        if(present(remove_val))rv=remove_val
        allocate(out(size(dat,1),size(dat,2)))
        out=dat
        do j=1,size(dat,2)
        do i=1,min(j-1,size(dat,1))
        out(i,j)=rv
        end do
        end do
    end function upper_tri_remove

    pure function graph_transpose(dat) result(out)
        real(dp),intent(in)::dat(:,:)
        real(dp),allocatable::out(:,:)
        allocate(out(size(dat,2),size(dat,1)))
        out=transpose(dat)
    end function graph_transpose

    pure function add_isolates(dat,nadd) result(out)
        real(dp),intent(in)::dat(:,:)
        integer,intent(in)::nadd
        real(dp),allocatable::out(:,:)
        integer::n,m
        n=size(dat,1)
        m=size(dat,2)
        allocate(out(n+max(0,nadd),m+max(0,nadd)))
        out=0.0_dp
        out(1:n,1:m)=dat
    end function add_isolates

    pure function symmetrize(dat,rule) result(out)
        real(dp),intent(in)::dat(:,:)
        character(*),intent(in),optional::rule
        real(dp),allocatable::out(:,:)
        character(len=16)::r
        integer::i,j,n
        n=size(dat,1)
        allocate(out(n,size(dat,2)))
        out=dat
        r='weak'
        if(present(rule))r=adjustl(rule)
        do i=1,n
            do j=i+1,min(n,size(dat,2))
                select case(trim(r))
                case('upper')
                    out(j,i)=out(i,j)
                case('lower')
                    out(i,j)=out(j,i)
                case('strong')
                    if(.not.is_missing(dat(i,j)) .and. .not.is_missing(dat(j,i))) then
                        if(dat(i,j)/=0.0_dp .and. dat(j,i)/=0.0_dp) then
                            out(i,j)=1.0_dp
                            out(j,i)=1.0_dp
                        else
                            out(i,j)=0.0_dp
                            out(j,i)=0.0_dp
                        end if
                    else
                        out(i,j)=sna_nan()
                        out(j,i)=sna_nan()
                    end if
                case default ! weak
                    if((.not.is_missing(dat(i,j)) .and. dat(i,j)/=0.0_dp) .or. &
                       (.not.is_missing(dat(j,i)) .and. dat(j,i)/=0.0_dp)) then
                        out(i,j)=1.0_dp
                        out(j,i)=1.0_dp
                    else if(is_missing(dat(i,j)) .or. is_missing(dat(j,i))) then
                        out(i,j)=sna_nan()
                        out(j,i)=sna_nan()
                    else
                        out(i,j)=0.0_dp
                        out(j,i)=0.0_dp
                    end if
                end select
            end do
        end do
    end function symmetrize

    function event2dichot(m,method,thresh,leq) result(out)
        real(dp),intent(in)::m(:,:)
        character(*),intent(in),optional::method
        real(dp),intent(in),optional::thresh
        logical,intent(in),optional::leq
        real(dp),allocatable::out(:,:)
        character(len=20)::meth
        real(dp)::th,q
        logical::flip
        integer::i,j
        real(dp),allocatable::v(:)
        meth='quantile'
        if(present(method))meth=adjustl(method)
        th=0.5_dp
        if(present(thresh))th=thresh
        flip=.false.
        if(present(leq))flip=leq
        allocate(out(size(m,1),size(m,2)))
        out=0.0_dp
        select case(trim(meth))
        case('quantile')
            v=pack(m,.not.is_missing(m))
            q=sample_quantile(v,th)
            where(.not.is_missing(m)) out=merge(1.0_dp,0.0_dp,m>q)
            where(is_missing(m)) out=sna_nan()
        case('rquantile')
            do i=1,size(m,1)
                v=pack(m(i,:),.not.is_missing(m(i,:)))
                q=sample_quantile(v,th)
                do j=1,size(m,2)
                    if(is_missing(m(i,j)))then
                    out(i,j)=sna_nan()
                    else
                    out(i,j)=merge(1.0_dp,0.0_dp,m(i,j)>q)
                    end if
                end do
            end do
        case('cquantile')
            do j=1,size(m,2)
                v=pack(m(:,j),.not.is_missing(m(:,j)))
                q=sample_quantile(v,th)
                do i=1,size(m,1)
                    if(is_missing(m(i,j)))then
                    out(i,j)=sna_nan()
                    else
                    out(i,j)=merge(1.0_dp,0.0_dp,m(i,j)>q)
                    end if
                end do
            end do
        case('mean')
            v=pack(m,.not.is_missing(m))
            q=sum(v)/real(max(1,size(v)),dp)
            where(.not.is_missing(m)) out=merge(1.0_dp,0.0_dp,m>q)
            where(is_missing(m))out=sna_nan()
        case('rmean')
            do i=1,size(m,1)
                v=pack(m(i,:),.not.is_missing(m(i,:)))
                q=sum(v)/real(max(1,size(v)),dp)
                do j=1,size(m,2)
                    if(is_missing(m(i,j)))then
                    out(i,j)=sna_nan()
                    else
                    out(i,j)=merge(1.0_dp,0.0_dp,m(i,j)>q)
                    end if
                end do
            end do
        case('cmean')
            do j=1,size(m,2)
                v=pack(m(:,j),.not.is_missing(m(:,j)))
                q=sum(v)/real(max(1,size(v)),dp)
                do i=1,size(m,1)
                    if(is_missing(m(i,j)))then
                    out(i,j)=sna_nan()
                    else
                    out(i,j)=merge(1.0_dp,0.0_dp,m(i,j)>q)
                    end if
                end do
            end do
        case('absolute')
            where(.not.is_missing(m))out=merge(1.0_dp,0.0_dp,m>th)
            where(is_missing(m))out=sna_nan()
        case('rank')
            call threshold_rank_matrix(m,th,out)
        case('rrank')
            do i=1,size(m,1)
            call threshold_rank_vector(m(i,:),th,out(i,:))
            end do
        case('crank')
            do j=1,size(m,2)
            call threshold_rank_vector(m(:,j),th,out(:,j))
            end do
        case default
            error stop 'event2dichot: unknown method'
        end select
        if(flip) where(.not.is_missing(out)) out=1.0_dp-out
    end function event2dichot

    subroutine threshold_rank_matrix(m,th,out)
        real(dp),intent(in)::m(:,:),th
        real(dp),intent(out)::out(:,:)
        real(dp),allocatable::v(:),o(:)
        integer::i,j,k
        allocate(v(size(m)),o(size(m)))
        v=reshape(m,[size(m)])
        call threshold_rank_vector(v,th,o)
        out=reshape(o,shape(out))
    end subroutine threshold_rank_matrix

    subroutine threshold_rank_vector(v,th,out)
        real(dp),intent(in)::v(:),th
        real(dp),intent(out)::out(:)
        integer,allocatable::idx(:)
        integer::i,j,t,n
        n=size(v)
        allocate(idx(n))
        idx=[(i,i=1,n)]
        do i=1,n-1
            do j=i+1,n
                if((is_missing(v(idx(i))) .and. .not.is_missing(v(idx(j)))) .or. &
                   (.not.is_missing(v(idx(i))) .and. .not.is_missing(v(idx(j))) .and. v(idx(i))>v(idx(j)))) then
                    t=idx(i)
                    idx(i)=idx(j)
                    idx(j)=t
                end if
            end do
        end do
        out=0.0_dp
        do i=1,n
            if(is_missing(v(i))) then
            out(i)=sna_nan()
            cycle
            end if
            do j=1,n
                if(idx(j)==i) exit
            end do
            out(i)=merge(1.0_dp,0.0_dp,real(n-j+1,dp)<th)
        end do
    end subroutine threshold_rank_vector

    function make_stochastic(dat,mode,tol,maxiter) result(out)
        real(dp),intent(in)::dat(:,:)
        character(*),intent(in),optional::mode
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter
        real(dp),allocatable::out(:,:)
        character(len=12)::md
        real(dp)::s,atol,err
        integer::i,j,it,mit,n,m
        n=size(dat,1)
        m=size(dat,2)
        allocate(out(n,m))
        out=dat
        md='rowcol'
        if(present(mode))md=adjustl(mode)
        atol=0.005_dp
        if(present(tol))atol=tol
        mit=n*m*100
        if(present(maxiter))mit=maxiter
        select case(trim(md))
        case('row')
            do i=1,n
            s=sum(out(i,:),mask=.not.is_missing(out(i,:)))
            if(abs(s)>sna_eps)out(i,:)=out(i,:)/s
            end do
        case('col')
            do j=1,m
            s=sum(out(:,j),mask=.not.is_missing(out(:,j)))
            if(abs(s)>sna_eps)out(:,j)=out(:,j)/s
            end do
        case('total')
            s=sum(out,mask=.not.is_missing(out))
            if(abs(s)>sna_eps)out=out/s
        case('rowcol')
            ! Iterative proportional fitting is deterministic and converges to the same
            ! row/column stochastic target when the support admits one.  The R routine
            ! uses a stochastic annealer after row/column seeding.
            do it=1,mit
                do i=1,n
                    s=sum(out(i,:),mask=.not.is_missing(out(i,:)))
                    if(abs(s)>sna_eps)out(i,:)=out(i,:)/s
                end do
                do j=1,m
                    s=sum(out(:,j),mask=.not.is_missing(out(:,j)))
                    if(abs(s)>sna_eps)out(:,j)=out(:,j)/s
                end do
                err=0.0_dp
                do i=1,n
                err=err+abs(sum(out(i,:),mask=.not.is_missing(out(i,:)))-1.0_dp)
                end do
                do j=1,m
                err=err+abs(sum(out(:,j),mask=.not.is_missing(out(:,j)))-1.0_dp)
                end do
                if(err<=real(n+m,dp)*atol)exit
            end do
        case default
            error stop 'make_stochastic: unknown mode'
        end select
        where(is_missing(out)) out=0.0_dp
    end function make_stochastic

    pure integer function nties(dat,mode,diag) result(count)
        real(dp),intent(in)::dat(:,:)
        character(*),intent(in),optional::mode
        logical,intent(in),optional::diag
        character(len=12)::md
        logical::dg
        integer::n,m
        n=size(dat,1)
        m=size(dat,2)
        md='digraph'
        if(present(mode))md=adjustl(mode)
        dg=.false.
        if(present(diag))dg=diag
        select case(trim(md))
        case('graph')
        count=(n*n-n)/2+n
        case('hgraph','twomode')
        count=n*m
        dg=.true.
        case default
        count=n*n
        end select
        if(.not.dg)count=count-n
    end function nties

    function gvectorize(m,mode,diag,censor_as_na) result(v)
        real(dp),intent(in)::m(:,:)
        character(*),intent(in),optional::mode
        logical,intent(in),optional::diag,censor_as_na
        real(dp),allocatable::v(:)
        character(len=12)::md
        logical::dg,ca
        integer::i,j,k,n,nn
        md='digraph'
        if(present(mode))md=adjustl(mode)
        dg=.false.
        if(present(diag))dg=diag
        ca=.true.
        if(present(censor_as_na))ca=censor_as_na
        n=size(m,1)
        if(ca)then
            allocate(v(size(m)))
            v=reshape(m,[size(m)])
            do j=1,size(m,2)
            do i=1,size(m,1)
            k=i+(j-1)*size(m,1)
                if((.not.dg .and. i==j).or.(trim(md)=='graph'.and.i<j))v(k)=sna_nan()
            end do
            end do
        else
            nn=nties(m,md,dg)
            allocate(v(nn))
            k=0
            do j=1,size(m,2)
            do i=1,size(m,1)
                if(.not.dg.and.i==j)cycle
                if(trim(md)=='graph'.and.i<j)cycle
                k=k+1
                v(k)=m(i,j)
            end do
            end do
        end if
    end function gvectorize

    function interval_graph(spells,type,diag) result(g)
        ! spells(:,1)=type, spells(:,2)=onset, spells(:,3)=termination
        real(dp),intent(in)::spells(:,:)
        character(*),intent(in),optional::type
        logical,intent(in),optional::diag
        real(dp),allocatable::g(:,:)
        character(len=16)::tp
        logical::dg
        integer::i,j,n
        real(dp)::ov,li,lj
        n=size(spells,1)
        allocate(g(n,n))
        tp='simple'
        if(present(type))tp=adjustl(type)
        dg=.false.
        if(present(diag))dg=diag
        do i=1,n
        do j=1,n
            ov=max(min(spells(i,3),spells(j,3))-max(spells(i,2),spells(j,2)),0.0_dp)
            li=spells(i,3)-spells(i,2)
            lj=spells(j,3)-spells(j,2)
            select case(trim(tp))
            case('simple')
            g(i,j)=merge(1.0_dp,0.0_dp,spells(i,2)<=spells(j,3).and.spells(i,3)>=spells(j,2))
            case('overlap')
            g(i,j)=ov
            case('fracxy')
            g(i,j)=merge(ov/li,0.0_dp,abs(li)>sna_eps)
            case('fracyx')
            g(i,j)=merge(ov/lj,0.0_dp,abs(lj)>sna_eps)
            case('jntfrac')
            g(i,j)=merge(2.0_dp*ov/(li+lj),0.0_dp,abs(li+lj)>sna_eps)
            case default
            error stop 'interval_graph: unknown type'
            end select
        end do
        end do
        if(.not.dg) then
            do i=1,n
                g(i,i)=0.0_dp
            end do
        end if
    end function interval_graph

    function ego_extract_mask(dat,ego,neighborhood) result(mask)
        real(dp),intent(in)::dat(:,:)
        integer,intent(in)::ego
        character(*),intent(in),optional::neighborhood
        logical,allocatable::mask(:)
        character(len=16)::nh
        integer::i
        allocate(mask(size(dat,1)))
        mask=.false.
        nh='combined'
        if(present(neighborhood))nh=adjustl(neighborhood)
        mask(ego)=.true.
        do i=1,size(dat,1)
            select case(trim(nh))
            case('in')
            mask(i)=mask(i).or.(.not.is_missing(dat(i,ego)).and.dat(i,ego)>0.0_dp)
            case('out')
            mask(i)=mask(i).or.(.not.is_missing(dat(ego,i)).and.dat(ego,i)>0.0_dp)
            case default
            mask(i)=mask(i).or.(.not.is_missing(dat(i,ego)).and.dat(i,ego)>0.0_dp).or. &
                (.not.is_missing(dat(ego,i)).and.dat(ego,i)>0.0_dp)
            end select
        end do
    end function ego_extract_mask


    function sr2css(net) result(css)
        ! Convert row-wise self reports to a CSS stack, matching R sna::sr2css.
        real(dp), intent(in) :: net(:,:)
        real(dp), allocatable :: css(:,:,:)
        integer :: n, i
        n=size(net,1)
        if(size(net,2)/=n) error stop 'sr2css: square matrix required'
        allocate(css(n,n,n))
        css=sna_nan()
        do i=1,n
            css(i,i,:)=net(i,:)
        end do
    end function sr2css

    pure real(dp) function log_sum(x) result(v)
        real(dp),intent(in)::x(:)
        real(dp)::mx
        if(size(x)==0)then
        v=-huge(1.0_dp)
        return
        end if
        mx=maxval(x)
        if(mx<=-huge(1.0_dp)/2)then
        v=mx
        else
        v=mx+log(sum(exp(x-mx)))
        end if
    end function log_sum

    pure real(dp) function log_mean(x) result(v)
        real(dp),intent(in)::x(:)
        if(size(x)==0)then
        v=sna_nan()
        else
        v=log_sum(x)-log(real(size(x),dp))
        end if
    end function log_mean

    pure function log_sub(x,y) result(v)
        real(dp),intent(in)::x(:),y(:)
        real(dp),allocatable::v(:)
        integer::i
        allocate(v(size(x)))
        do i=1,size(x)
            if(y(i)>=x(i))then
            v(i)=sna_nan()
            else
            v(i)=x(i)+log(1.0_dp-exp(y(i)-x(i)))
            end if
        end do
    end function log_sub

    function sample_quantile(x,p) result(q)
        real(dp),intent(in)::x(:),p
        real(dp)::q,h,frac,tmp
        real(dp),allocatable::z(:)
        integer::i,j,n,k
        n=size(x)
        if(n==0)then
        q=sna_nan()
        return
        end if
        allocate(z(n))
        z=x
        do i=1,n-1
        do j=i+1,n
        if(z(j)<z(i))then
        tmp=z(i)
        z(i)=z(j)
        z(j)=tmp
        end if
        end do
        end do
        if(n==1)then
        q=z(1)
        return
        end if
        ! R quantile type 7: h=(n-1)*p+1
        h=(real(n-1,dp)*min(max(p,0.0_dp),1.0_dp))+1.0_dp
        k=floor(h)
        frac=h-real(k,dp)
        if(k>=n)then
        q=z(n)
        else
        q=(1.0_dp-frac)*z(k)+frac*z(k+1)
        end if
    end function sample_quantile

end module sna_prep
