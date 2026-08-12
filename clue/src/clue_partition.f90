! SPDX-License-Identifier: GPL-2.0-only
module clue_partition
    use clue_kinds, only: dp
    implicit none
    private
    public :: hard_partition, soft_partition, make_hard_partition, make_soft_partition
    public :: membership_from_ids, class_ids_from_membership, canonicalize_ids
    public :: co_membership, partition_margin, fuzziness_pc, fuzziness_pe
    public :: partition_meet, partition_join, memberships_from_cross_dissimilarities
    public :: n_of_classes, n_of_objects

    type :: hard_partition
        integer, allocatable :: id(:)
    contains
        procedure :: nobjects => hard_nobjects
        procedure :: nclasses => hard_nclasses
        procedure :: membership => hard_membership
    end type hard_partition

    type :: soft_partition
        real(dp), allocatable :: u(:,:)
    contains
        procedure :: nobjects => soft_nobjects
        procedure :: nclasses => soft_nclasses
        procedure :: class_ids => soft_class_ids
    end type soft_partition

    interface n_of_classes
        module procedure nclasses_ids, nclasses_membership
    end interface
    interface n_of_objects
        module procedure nobjects_ids, nobjects_membership
    end interface
contains
    function make_hard_partition(ids) result(p)
        integer, intent(in) :: ids(:)
        type(hard_partition) :: p
        allocate(p%id(size(ids)))
        p%id=canonicalize_ids(ids)
    end function

    function make_soft_partition(u, normalize) result(p)
        real(dp), intent(in) :: u(:,:)
        logical, intent(in), optional :: normalize
        type(soft_partition) :: p
        logical :: norm
        integer :: i
        real(dp) :: s
        norm=.true.
        if(present(normalize)) norm=normalize
        allocate(p%u(size(u,1),size(u,2)))
        p%u=max(u,0.0_dp)
        if(norm) then
            do i=1,size(u,1)
                s=sum(p%u(i,:))
                if(s>0.0_dp) then
                    p%u(i,:)=p%u(i,:)/s
                else if(size(u,2)>0) then
                    p%u(i,:)=1.0_dp/real(size(u,2),dp)
                end if
            end do
        end if
    end function

    function canonicalize_ids(ids) result(out)
        integer, intent(in) :: ids(:)
        integer, allocatable :: out(:), seen(:)
        integer :: i,j,k,nseen
        allocate(out(size(ids)),seen(size(ids)))
        seen=0
        nseen=0
        do i=1,size(ids)
            k=0
            do j=1,nseen
                if(seen(j)==ids(i)) then
                k=j
                exit
                end if
            end do
            if(k==0) then
                nseen=nseen+1
                seen(nseen)=ids(i)
                k=nseen
            end if
            out(i)=k
        end do
    end function

    function membership_from_ids(ids,k) result(m)
        integer, intent(in) :: ids(:)
        integer, intent(in), optional :: k
        real(dp), allocatable :: m(:,:)
        integer, allocatable :: cid(:)
        integer :: kk,i
        cid=canonicalize_ids(ids)
        if(size(ids)==0) then
            kk=0
        else
            kk=maxval(cid)
        end if
        if(present(k)) kk=max(kk,k)
        allocate(m(size(ids),kk))
        m=0.0_dp
        do i=1,size(ids)
            if(cid(i)>=1 .and. cid(i)<=kk) m(i,cid(i))=1.0_dp
        end do
    end function

    function class_ids_from_membership(u) result(ids)
        real(dp), intent(in) :: u(:,:)
        integer, allocatable :: ids(:)
        integer :: i
        allocate(ids(size(u,1)))
        ids=0
        do i=1,size(u,1)
            if(size(u,2)>0) ids(i)=maxloc(u(i,:),dim=1)
        end do
    end function

    function co_membership(u) result(c)
        real(dp), intent(in) :: u(:,:)
        real(dp), allocatable :: c(:,:)
        allocate(c(size(u,1),size(u,1)))
        c=matmul(u,transpose(u))
    end function

    function partition_margin(u) result(m)
        real(dp), intent(in) :: u(:,:)
        real(dp), allocatable :: m(:)
        real(dp) :: a,b
        integer :: i,j,jmax
        allocate(m(size(u,1)))
        do i=1,size(u,1)
            if(size(u,2)<2) then
            m(i)=1.0_dp
            cycle
            end if
            jmax=maxloc(u(i,:),dim=1)
            a=u(i,jmax)
            b=-huge(1.0_dp)
            do j=1,size(u,2)
                if(j/=jmax) b=max(b,u(i,j))
            end do
            m(i)=a-b
        end do
    end function

    function fuzziness_pc(u,normalize) result(v)
        real(dp), intent(in) :: u(:,:)
        logical, intent(in), optional :: normalize
        real(dp) :: v
        logical :: norm
        integer :: k
        norm=.true.
        if(present(normalize)) norm=normalize
        if(size(u,1)==0 .or. size(u,2)==0) then
        v=0.0_dp
        return
        end if
        v=sum(u*u)/real(size(u,1),dp)
        k=size(u,2)
        if(norm .and. k>1) v=(1.0_dp-v)/(1.0_dp-1.0_dp/real(k,dp))
    end function

    function fuzziness_pe(u,normalize) result(v)
        real(dp), intent(in) :: u(:,:)
        logical, intent(in), optional :: normalize
        real(dp) :: v
        logical :: norm
        integer :: i,j,k
        norm=.true.
        if(present(normalize)) norm=normalize
        v=0.0_dp
        do j=1,size(u,2)
        do i=1,size(u,1)
            if(u(i,j)>0.0_dp) v=v-u(i,j)*log(u(i,j))
        end do
        end do
        if(size(u,1)>0) v=v/real(size(u,1),dp)
        k=size(u,2)
        if(norm .and. k>1) v=v/log(real(k,dp))
    end function

    function partition_meet(a,b) result(ids)
        integer, intent(in) :: a(:),b(:)
        integer, allocatable :: ids(:),ca(:),cb(:)
        integer :: i,j,k,n
        n=size(a)
        if(size(b)/=n) then
        allocate(ids(0))
        return
        end if
        ca=canonicalize_ids(a)
        cb=canonicalize_ids(b)
        allocate(ids(n))
        ids=0
        k=0
        do i=1,n
            do j=1,i-1
                if(ca(j)==ca(i) .and. cb(j)==cb(i)) then
                ids(i)=ids(j)
                exit
                end if
            end do
            if(i==1 .or. ids(i)==0) then
            k=k+1
            ids(i)=k
            end if
        end do
    end function

    function partition_join(a,b) result(ids)
        integer, intent(in) :: a(:),b(:)
        integer, allocatable :: ids(:),ca(:),cb(:),parent(:),rootmap(:)
        integer :: n,i,j,ra,k
        n=size(a)
        if(size(b)/=n) then
        allocate(ids(0))
        return
        end if
        ca=canonicalize_ids(a)
        cb=canonicalize_ids(b)
        allocate(parent(n),rootmap(n),ids(n))
        parent=[(i,i=1,n)]
        rootmap=0
        do i=1,n-1
            do j=i+1,n
                if(ca(i)==ca(j) .or. cb(i)==cb(j)) call unite(i,j,parent)
            end do
        end do
        k=0
        do i=1,n
            ra=find_root(i,parent)
            if(rootmap(ra)==0) then
            k=k+1
            rootmap(ra)=k
            end if
            ids(i)=rootmap(ra)
        end do
    contains
        recursive integer function find_root(x,p) result(r)
            integer,intent(in)::x
            integer,intent(inout)::p(:)
            if(p(x)==x) then
            r=x
            else
            p(x)=find_root(p(x),p)
            r=p(x)
            end if
        end function
        subroutine unite(x,y,p)
            integer,intent(in)::x,y
            integer,intent(inout)::p(:)
            integer::rx,ry
            rx=find_root(x,p)
            ry=find_root(y,p)
            if(rx/=ry) p(ry)=rx
        end subroutine
    end function

    function memberships_from_cross_dissimilarities(d,power) result(u)
        real(dp),intent(in)::d(:,:)
        real(dp),intent(in),optional::power(:)
        real(dp),allocatable::u(:,:)
        real(dp)::expo,mx,s
        integer::i,j,nz
        logical,allocatable::z(:)
        if(present(power)) then
            if(size(power)==1) then
            expo=1.0_dp/(1.0_dp-power(1))
            else
            expo=power(2)/(1.0_dp-power(1))
            end if
        else
        expo=-1.0_dp
        end if
        allocate(u(size(d,1),size(d,2)))
        u=0.0_dp
        allocate(z(size(d,2)))
        do i=1,size(d,1)
            z=.not.(d(i,:)>0.0_dp)
            nz=count(z)
            if(nz>0) then
                do j=1,size(d,2)
                if(z(j)) u(i,j)=1.0_dp/real(nz,dp)
                end do
            else
                mx=maxval(expo*log(d(i,:)))
                u(i,:)=exp(expo*log(d(i,:))-mx)
                s=sum(u(i,:))
                if(s>0) u(i,:)=u(i,:)/s
            end if
        end do
    end function

    integer function hard_nobjects(self) result(n)
        class(hard_partition),intent(in)::self
        n=size(self%id)
    end function
    integer function hard_nclasses(self) result(n)
        class(hard_partition),intent(in)::self
        n=nclasses_ids(self%id)
    end function
    function hard_membership(self,k) result(m)
        class(hard_partition),intent(in)::self
        integer,intent(in),optional::k
        real(dp),allocatable::m(:,:)
        if(present(k)) then
        m=membership_from_ids(self%id,k)
        else
        m=membership_from_ids(self%id)
        end if
    end function
    integer function soft_nobjects(self) result(n)
        class(soft_partition),intent(in)::self
        n=size(self%u,1)
    end function
    integer function soft_nclasses(self) result(n)
        class(soft_partition),intent(in)::self
        n=size(self%u,2)
    end function
    function soft_class_ids(self) result(ids)
        class(soft_partition),intent(in)::self
        integer,allocatable::ids(:)
        ids=class_ids_from_membership(self%u)
    end function
    integer function nclasses_ids(ids) result(n)
        integer,intent(in)::ids(:)
        integer,allocatable::c(:)
        if(size(ids)==0) then
        n=0
        else
        c=canonicalize_ids(ids)
        n=maxval(c)
        end if
    end function
    integer function nclasses_membership(u) result(n)
        real(dp),intent(in)::u(:,:)
        n=size(u,2)
    end function
    integer function nobjects_ids(ids) result(n)
        integer,intent(in)::ids(:)
        n=size(ids)
    end function
    integer function nobjects_membership(u) result(n)
        real(dp),intent(in)::u(:,:)
        n=size(u,1)
    end function
end module clue_partition
