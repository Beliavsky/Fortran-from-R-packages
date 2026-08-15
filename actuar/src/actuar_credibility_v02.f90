module actuar_credibility_v02
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use actuar_kinds, only: dp
    implicit none
    private
    public :: credibility_level_t, hierarchical_credibility_result_t, hierarchical_credibility

    type :: credibility_level_t
        real(dp), allocatable :: weight(:)
        real(dp), allocatable :: mean(:)
        real(dp), allocatable :: credibility(:)
        real(dp), allocatable :: premium(:)
        integer, allocatable :: parent(:)
    end type credibility_level_t

    type :: hierarchical_credibility_result_t
        integer :: nlevels = 0
        real(dp) :: process_variance = 0.0_dp
        real(dp), allocatable :: between_variance(:)
        type(credibility_level_t), allocatable :: level(:)
    end type hierarchical_credibility_result_t
contains
    function hierarchical_credibility(ratios,weights,classification,ohlsson) result(res)
        real(dp),intent(in)::ratios(:,:),weights(:,:)
        integer,intent(in)::classification(:,:)
        logical,intent(in),optional::ohlsson
        type(hierarchical_credibility_result_t)::res
        integer::nr,np,l,k,i,j,g,p,c,neff,nobs,nent_eff
        integer,allocatable::gid(:,:),gcount(:),eff(:),child_parent(:)
        real(dp),allocatable::roww(:),rowm(:),bi(:),ci(:)
        real(dp)::s2,between,num,den
        logical::use_ohlsson
        nr=size(ratios,1);np=size(ratios,2);l=size(classification,2)
        if(size(weights,1)/=nr .or. size(weights,2)/=np .or. size(classification,1)/=nr .or. l<1) return
        use_ohlsson=.false.;if(present(ohlsson)) use_ohlsson=ohlsson
        allocate(gid(l,nr),gcount(l));call build_group_ids(classification,gid,gcount)
        res%nlevels=l;allocate(res%between_variance(l),res%level(0:l));res%between_variance=0.0_dp
        allocate(roww(nr),rowm(nr));roww=0.0_dp;rowm=0.0_dp;nobs=0;nent_eff=0
        do i=1,nr
            do j=1,np
                if(weights(i,j)>0.0_dp .and. .not.ieee_is_nan(ratios(i,j)) .and. .not.ieee_is_nan(weights(i,j))) then
                    roww(i)=roww(i)+weights(i,j);rowm(i)=rowm(i)+weights(i,j)*ratios(i,j);nobs=nobs+1
                end if
            end do
            if(roww(i)>0.0_dp) then;rowm(i)=rowm(i)/roww(i);nent_eff=nent_eff+1;end if
        end do
        s2=0.0_dp
        do i=1,nr
            do j=1,np
                if(weights(i,j)>0.0_dp .and. .not.ieee_is_nan(ratios(i,j)) .and. .not.ieee_is_nan(weights(i,j))) then
                    s2=s2+weights(i,j)*(ratios(i,j)-rowm(i))**2
                end if
            end do
        end do
        if(nobs>nent_eff) s2=s2/real(nobs-nent_eff,dp)
        res%process_variance=s2

        allocate(res%level(l)%weight(gcount(l)),res%level(l)%mean(gcount(l)))
        res%level(l)%weight=0.0_dp;res%level(l)%mean=0.0_dp
        do i=1,nr
            g=gid(l,i)
            res%level(l)%weight(g)=res%level(l)%weight(g)+roww(i)
            res%level(l)%mean(g)=res%level(l)%mean(g)+roww(i)*rowm(i)
        end do
        do g=1,gcount(l)
            if(res%level(l)%weight(g)>0.0_dp) res%level(l)%mean(g)=res%level(l)%mean(g)/res%level(l)%weight(g)
        end do

        do k=l,1,-1
            if(k==1) then
                p=1
            else
                p=gcount(k-1)
            end if
            c=gcount(k)
            allocate(res%level(k)%parent(c));call map_child_to_parent(k,gid,gcount,res%level(k)%parent)
            allocate(res%level(k-1)%weight(p),res%level(k-1)%mean(p))
            res%level(k-1)%weight=0.0_dp;res%level(k-1)%mean=0.0_dp
            do i=1,c
                g=res%level(k)%parent(i)
                res%level(k-1)%weight(g)=res%level(k-1)%weight(g)+res%level(k)%weight(i)
                res%level(k-1)%mean(g)=res%level(k-1)%mean(g)+res%level(k)%weight(i)*res%level(k)%mean(i)
            end do
            do g=1,p
                if(res%level(k-1)%weight(g)>0.0_dp) then
                    res%level(k-1)%mean(g)=res%level(k-1)%mean(g)/res%level(k-1)%weight(g)
                end if
            end do
            between=s2
            if(k<l) then
                do j=k+1,l
                    if(res%between_variance(j)/=0.0_dp) then;between=res%between_variance(j);exit;end if
                end do
            end if
            allocate(bi(p),ci(p),eff(p));bi=0.0_dp;ci=0.0_dp;eff=0
            do i=1,c
                g=res%level(k)%parent(i)
                if(res%level(k)%weight(i)>0.0_dp) eff(g)=eff(g)+1
                bi(g)=bi(g)+res%level(k)%weight(i)*(res%level(k)%mean(i)-res%level(k-1)%mean(g))**2
                ci(g)=ci(g)+res%level(k)%weight(i)**2
            end do
            do g=1,p
                bi(g)=bi(g)-real(max(0,eff(g)-1),dp)*between
                if(res%level(k-1)%weight(g)>0.0_dp) then
                    ci(g)=res%level(k-1)%weight(g)-ci(g)/res%level(k-1)%weight(g)
                else
                    ci(g)=0.0_dp
                end if
            end do
            if(use_ohlsson) then
                num=0.0_dp;den=0.0_dp
                do g=1,p
                    if(ci(g)>0.0_dp) then;num=num+bi(g);den=den+ci(g);end if
                end do
                if(den>0.0_dp) res%between_variance(k)=num/den
            else
                num=0.0_dp;neff=0
                do g=1,p
                    if(ci(g)>0.0_dp) then;num=num+max(bi(g)/ci(g),0.0_dp);neff=neff+1;end if
                end do
                if(neff>0) res%between_variance(k)=num/real(neff,dp)
            end if
            allocate(res%level(k)%credibility(c));res%level(k)%credibility=0.0_dp
            if(res%between_variance(k)>0.0_dp) then
                do i=1,c
                    if(res%level(k)%weight(i)>0.0_dp) then
                        res%level(k)%credibility(i)=1.0_dp/(1.0_dp+between/(res%between_variance(k)*res%level(k)%weight(i)))
                    end if
                end do
                res%level(k-1)%weight=0.0_dp;res%level(k-1)%mean=0.0_dp
                do i=1,c
                    g=res%level(k)%parent(i)
                    res%level(k-1)%weight(g)=res%level(k-1)%weight(g)+res%level(k)%credibility(i)
                    res%level(k-1)%mean(g)=res%level(k-1)%mean(g)+res%level(k)%credibility(i)*res%level(k)%mean(i)
                end do
                do g=1,p
                    if(res%level(k-1)%weight(g)>0.0_dp) then
                        res%level(k-1)%mean(g)=res%level(k-1)%mean(g)/res%level(k-1)%weight(g)
                    end if
                end do
            end if
            deallocate(bi,ci,eff)
        end do

        allocate(res%level(0)%premium(1));res%level(0)%premium=res%level(0)%mean
        do k=1,l
            c=gcount(k);allocate(res%level(k)%premium(c))
            do i=1,c
                p=res%level(k)%parent(i)
                res%level(k)%premium(i)=res%level(k-1)%premium(p)+res%level(k)%credibility(i) &
                    *(res%level(k)%mean(i)-res%level(k-1)%premium(p))
            end do
        end do
    end function hierarchical_credibility

    subroutine build_group_ids(classification,gid,gcount)
        integer,intent(in)::classification(:,:)
        integer,intent(out)::gid(:,:),gcount(:)
        integer::nr,l,k,i,j,g
        logical::same
        nr=size(classification,1);l=size(classification,2);gid=0;gcount=0
        do k=1,l
            do i=1,nr
                g=0
                do j=1,i-1
                    same=all(classification(i,1:k)==classification(j,1:k))
                    if(same) then;g=gid(k,j);exit;end if
                end do
                if(g==0) then;gcount(k)=gcount(k)+1;g=gcount(k);end if
                gid(k,i)=g
            end do
        end do
    end subroutine build_group_ids

    subroutine map_child_to_parent(k,gid,gcount,parent)
        integer,intent(in)::k,gid(:,:),gcount(:)
        integer,intent(out)::parent(:)
        integer::row,c
        parent=0
        do row=1,size(gid,2)
            c=gid(k,row)
            if(parent(c)==0) then
                if(k==1) then
                    parent(c)=1
                else
                    parent(c)=gid(k-1,row)
                end if
            end if
        end do
    end subroutine map_child_to_parent
end module actuar_credibility_v02
