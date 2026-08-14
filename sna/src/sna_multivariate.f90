! Multivariate graph analysis translated from R/sna gmultiv.R.
! Upstream copyright (C) 2004-2024 Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_multivariate
    use sna_kinds, only : dp, sna_nan, is_missing, sna_eps
    use sna_prep, only : nties
    use sna_permutation, only : lab_optimize_mc, lab_optimize_hillclimb, lab_optimize_anneal, lab_optimize_exhaustive
    implicit none
    private
    public :: centralgraph, graph_covariance, graph_correlation, hdist
    public :: structdist, sdmat, gscov, gscor
contains

    function centralgraph(dat,normalize) result(out)
        real(dp), intent(in) :: dat(:,:,:)
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: out(:,:)
        integer :: i,j,k,n,m,c
        real(dp) :: s
        logical :: norm
        n=size(dat,2)
        m=size(dat,1)
        norm=.false.
        if(present(normalize)) norm=normalize
        allocate(out(n,n))
        do i=1,n
            do j=1,n
                s=0.0_dp
                c=0
                do k=1,m
                    if(.not.is_missing(dat(k,i,j))) then
                    s=s+dat(k,i,j)
                    c=c+1
                    end if
                end do
                if(c==0) then
                    out(i,j)=sna_nan()
                else if(norm) then
                    out(i,j)=s/real(c,dp)
                else
                    out(i,j)=merge(1.0_dp,0.0_dp,s/real(c,dp)>=0.5_dp)
                end if
            end do
        end do
    end function centralgraph

    function graph_covariance(a,b,diag,mode) result(v)
        real(dp), intent(in) :: a(:,:),b(:,:)
        logical,intent(in),optional :: diag
        character(len=*),intent(in),optional :: mode
        real(dp)::v
        real(dp),allocatable::aa(:,:),bb(:,:)
        aa=masked(a,diag,mode)
        bb=masked(b,diag,mode)
        v=matrix_covariance(aa,bb)
    end function graph_covariance

    function graph_correlation(a,b,diag,mode) result(v)
        real(dp), intent(in) :: a(:,:),b(:,:)
        logical,intent(in),optional :: diag
        character(len=*),intent(in),optional :: mode
        real(dp)::v
        real(dp),allocatable::aa(:,:),bb(:,:)
        aa=masked(a,diag,mode)
        bb=masked(b,diag,mode)
        v=matrix_correlation(aa,bb)
    end function graph_correlation

    function hdist(a,b,normalize,diag,mode) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        logical,intent(in),optional::normalize,diag
        character(len=*),intent(in),optional::mode
        real(dp)::v
        real(dp),allocatable::aa(:,:),bb(:,:)
        logical::norm,dg
        character(len=16)::md
        integer::denom
        aa=masked(a,diag,mode)
        bb=masked(b,diag,mode)
        v=matrix_hamming(aa,bb)
        norm=.false.
        if(present(normalize))norm=normalize
        if(norm)then
            dg=.false.
            if(present(diag))dg=diag
            md='digraph'
            if(present(mode))md=trim(mode)
            denom=nties(a,md,dg)
            if(denom>0)v=v/real(denom,dp)
        end if
    end function hdist

    function gscov(a,b,method,exchange_list,diag,mode,reps) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        character(len=*),intent(in),optional::method,mode
        integer,intent(in),optional::exchange_list(:),reps
        logical,intent(in),optional::diag
        real(dp)::v
        real(dp),allocatable::aa(:,:),bb(:,:)
        aa=masked(a,diag,mode)
        bb=masked(b,diag,mode)
        v=optimize_pair(aa,bb,.true.,.false.,method,exchange_list,reps)
    end function gscov

    function gscor(a,b,method,exchange_list,diag,mode,reps) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        character(len=*),intent(in),optional::method,mode
        integer,intent(in),optional::exchange_list(:),reps
        logical,intent(in),optional::diag
        real(dp)::v
        real(dp),allocatable::aa(:,:),bb(:,:)
        aa=masked(a,diag,mode)
        bb=masked(b,diag,mode)
        v=optimize_pair(aa,bb,.true.,.true.,method,exchange_list,reps)
    end function gscor

    function structdist(a,b,normalize,diag,mode,method,exchange_list,reps) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        logical,intent(in),optional::normalize,diag
        character(len=*),intent(in),optional::mode,method
        integer,intent(in),optional::exchange_list(:),reps
        real(dp)::v
        real(dp),allocatable::aa(:,:),bb(:,:)
        logical::norm,dg
        character(len=16)::md
        integer::denom
        aa=masked(a,diag,mode)
        bb=masked(b,diag,mode)
        v=optimize_pair(aa,bb,.false.,.false.,method,exchange_list,reps)
        norm=.false.
        if(present(normalize))norm=normalize
        if(norm)then
            dg=.false.
            if(present(diag))dg=diag
            md='digraph'
            if(present(mode))md=trim(mode)
            denom=nties(a,md,dg)
            if(denom>0)v=v/real(denom,dp)
        end if
    end function structdist

    function sdmat(dat,normalize,diag,mode,method,exchange_list,reps) result(d)
        real(dp),intent(in)::dat(:,:,:)
        logical,intent(in),optional::normalize,diag
        character(len=*),intent(in),optional::mode,method
        integer,intent(in),optional::exchange_list(:),reps
        real(dp),allocatable::d(:,:)
        integer::i,j,m
        m=size(dat,1)
        allocate(d(m,m))
        d=0.0_dp
        do i=1,m
            do j=i+1,m
                d(i,j)=structdist(dat(i,:,:),dat(j,:,:),normalize,diag,mode,method,exchange_list,reps)
                d(j,i)=d(i,j)
            end do
        end do
    end function sdmat

    function optimize_pair(a,b,seek_max,is_cor,method,exchange_list,reps) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        logical,intent(in)::seek_max,is_cor
        character(len=*),intent(in),optional::method
        integer,intent(in),optional::exchange_list(:),reps
        real(dp)::v
        character(len=24)::meth
        integer::nr
        meth='anneal'
        if(present(method))meth=trim(method)
        nr=1000
        if(present(reps))nr=reps
        if(trim(meth)=='none')then
            if(seek_max)then
                if(is_cor)v=matrix_correlation(a,b)
                if(.not.is_cor)v=matrix_covariance(a,b)
            else
                v=matrix_hamming(a,b)
            end if
            return
        end if
        if(seek_max.and.is_cor)then
            select case(trim(meth))
            case('mc')
            v=lab_optimize_mc(a,b,matrix_correlation,exchange_list,.true.,nr)
            case('hillclimb')
            v=lab_optimize_hillclimb(a,b,matrix_correlation,exchange_list,.true.)
            case('exhaustive')
            v=lab_optimize_exhaustive(a,b,matrix_correlation,exchange_list,.true.)
            case default
            v=lab_optimize_anneal(a,b,matrix_correlation,exchange_list,.true.,0.9_dp,0.85_dp,25,.true.)
            end select
        else if(seek_max)then
            select case(trim(meth))
            case('mc')
            v=lab_optimize_mc(a,b,matrix_covariance,exchange_list,.true.,nr)
            case('hillclimb')
            v=lab_optimize_hillclimb(a,b,matrix_covariance,exchange_list,.true.)
            case('exhaustive')
            v=lab_optimize_exhaustive(a,b,matrix_covariance,exchange_list,.true.)
            case default
            v=lab_optimize_anneal(a,b,matrix_covariance,exchange_list,.true.,0.9_dp,0.85_dp,25,.true.)
            end select
        else
            select case(trim(meth))
            case('mc')
            v=lab_optimize_mc(a,b,matrix_hamming,exchange_list,.false.,nr)
            case('hillclimb')
            v=lab_optimize_hillclimb(a,b,matrix_hamming,exchange_list,.false.)
            case('exhaustive')
            v=lab_optimize_exhaustive(a,b,matrix_hamming,exchange_list,.false.)
            case default
            v=lab_optimize_anneal(a,b,matrix_hamming,exchange_list,.false.,0.9_dp,0.85_dp,25,.true.)
            end select
        end if
    end function optimize_pair

    function masked(a,diag,mode) result(x)
        real(dp),intent(in)::a(:,:)
        logical,intent(in),optional::diag
        character(len=*),intent(in),optional::mode
        real(dp),allocatable::x(:,:)
        logical::dg
        character(len=16)::md
        integer::i,j,n
        n=size(a,1)
        allocate(x(n,n))
        x=a
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        if(.not.dg) then
            do i=1,n
                x(i,i)=sna_nan()
            end do
        end if
        if(trim(md)=='graph')then
            do i=1,n
            do j=i+1,n
            x(i,j)=sna_nan()
            end do
            end do
        end if
    end function masked

    function matrix_covariance(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v,sx,sy
        integer::i,j,c
        sx=0.0_dp
        sy=0.0_dp
        c=0
        do i=1,size(a,1)
        do j=1,size(a,2)
            if(.not.is_missing(a(i,j)).and..not.is_missing(b(i,j)))then
            sx=sx+a(i,j)
            sy=sy+b(i,j)
            c=c+1
            end if
        end do
        end do
        if(c<2)then
        v=sna_nan()
        return
        end if
        sx=sx/real(c,dp)
        sy=sy/real(c,dp)
        v=0.0_dp
        do i=1,size(a,1)
        do j=1,size(a,2)
            if(.not.is_missing(a(i,j)).and..not.is_missing(b(i,j)))v=v+(a(i,j)-sx)*(b(i,j)-sy)
        end do
        end do
        v=v/real(c-1,dp)
    end function matrix_covariance

    function matrix_correlation(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v,sx,sy,sxx,syy,sxy
        integer::i,j,c
        sx=0
        sy=0
        c=0
        do i=1,size(a,1)
        do j=1,size(a,2)
            if(.not.is_missing(a(i,j)).and..not.is_missing(b(i,j)))then
            sx=sx+a(i,j)
            sy=sy+b(i,j)
            c=c+1
            end if
        end do
        end do
        if(c<2)then
        v=sna_nan()
        return
        end if
        sx=sx/real(c,dp)
        sy=sy/real(c,dp)
        sxx=0
        syy=0
        sxy=0
        do i=1,size(a,1)
        do j=1,size(a,2)
            if(.not.is_missing(a(i,j)).and..not.is_missing(b(i,j)))then
                sxx=sxx+(a(i,j)-sx)**2
                syy=syy+(b(i,j)-sy)**2
                sxy=sxy+(a(i,j)-sx)*(b(i,j)-sy)
            end if
        end do
        end do
        if(sxx<=sna_eps.or.syy<=sna_eps)then
        v=sna_nan()
        else
        v=sxy/sqrt(sxx*syy)
        end if
    end function matrix_correlation

    function matrix_hamming(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v
        integer::i,j
        v=0.0_dp
        do i=1,size(a,1)
        do j=1,size(a,2)
            if(.not.is_missing(a(i,j)).and..not.is_missing(b(i,j)))v=v+abs(a(i,j)-b(i,j))
        end do
        end do
    end function matrix_hamming
end module sna_multivariate
