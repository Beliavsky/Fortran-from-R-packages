module actuar_hier_sim_v03
    use actuar_kinds, only: dp
    implicit none
    private
    public :: hierarchy_node_counts_t, hierarchy_real_model_t, hierarchy_count_model_t
    public :: hierarchical_portfolio_t, rcomphierarc_simulate

    abstract interface
        subroutine hierarchy_real_draw(ancestors, values)
            import dp
            real(dp),intent(in)::ancestors(:,:)
            real(dp),intent(out)::values(:)
        end subroutine hierarchy_real_draw
        subroutine hierarchy_count_draw(ancestors, values)
            import dp
            real(dp),intent(in)::ancestors(:,:)
            integer,intent(out)::values(:)
        end subroutine hierarchy_count_draw
    end interface

    type :: hierarchy_node_counts_t
        integer,allocatable::count(:)
    end type hierarchy_node_counts_t
    type :: hierarchy_real_model_t
        procedure(hierarchy_real_draw),pointer,nopass::draw=>null()
    end type hierarchy_real_model_t
    type :: hierarchy_count_model_t
        procedure(hierarchy_count_draw),pointer,nopass::draw=>null()
    end type hierarchy_count_model_t
    type :: real_vector_t
        real(dp),allocatable::value(:)
    end type real_vector_t
    type :: int_vector_t
        integer,allocatable::value(:)
    end type int_vector_t
    type :: int_matrix_t
        integer,allocatable::value(:,:)
    end type int_matrix_t
    type :: hierarchical_portfolio_t
        integer::nlevels=0
        integer,allocatable::node_count(:)
        integer,allocatable::frequency(:),claim_node(:),terminal_path(:,:)
        real(dp),allocatable::claims(:),aggregate(:)
    end type hierarchical_portfolio_t
contains
    function rcomphierarc_simulate(nodes,freq_mix,freq_final,sev_mix,sev_final) result(res)
        type(hierarchy_node_counts_t),intent(in)::nodes(:)
        type(hierarchy_real_model_t),intent(in),optional::freq_mix(:),sev_mix(:)
        type(hierarchy_count_model_t),intent(in),optional::freq_final
        type(hierarchy_real_model_t),intent(in),optional::sev_final
        type(hierarchical_portfolio_t)::res
        type(int_vector_t),allocatable::parent(:)
        type(int_matrix_t),allocatable::path(:)
        type(real_vector_t),allocatable::fpar(:),spar(:)
        real(dp),allocatable::anc(:,:),claimanc(:,:)
        integer::l,k,npar,nnode,i,j,nclaims,pos,cnt,idx
        l=size(nodes);if(l<1)return
        allocate(parent(l),path(l),fpar(max(1,l-1)),spar(max(1,l-1)),res%node_count(l))
        npar=1
        do k=1,l
            call expand_level(nodes(k)%count,npar,parent(k)%value)
            nnode=size(parent(k)%value);res%node_count(k)=nnode
            allocate(path(k)%value(nnode,k))
            if(k==1)then
                do i=1,nnode;path(k)%value(i,1)=i;end do
            else
                do i=1,nnode
                    path(k)%value(i,:k-1)=path(k-1)%value(parent(k)%value(i),:)
                    path(k)%value(i,k)=i
                end do
            end if
            if(k<l)then
                allocate(fpar(k)%value(nnode),spar(k)%value(nnode));fpar(k)%value=0.0_dp;spar(k)%value=0.0_dp
                call make_ancestors(k,path(k)%value,fpar,anc)
                if(present(freq_mix))then
                    if(k<=size(freq_mix))then
                        if(associated(freq_mix(k)%draw))call freq_mix(k)%draw(anc,fpar(k)%value)
                    end if
                end if
                call make_ancestors(k,path(k)%value,spar,anc)
                if(present(sev_mix))then
                    if(k<=size(sev_mix))then
                        if(associated(sev_mix(k)%draw))call sev_mix(k)%draw(anc,spar(k)%value)
                    end if
                end if
            end if
            npar=nnode
        end do
        nnode=res%node_count(l);allocate(res%frequency(nnode),res%aggregate(nnode),res%terminal_path(nnode,l))
        res%terminal_path=path(l)%value;res%aggregate=0.0_dp
        call make_ancestors(l,path(l)%value,fpar,anc)
        if(present(freq_final))then
            if(associated(freq_final%draw))then;call freq_final%draw(anc,res%frequency)
            else;res%frequency=1;end if
        else;res%frequency=1;end if
        res%frequency=max(res%frequency,0);nclaims=sum(res%frequency)
        allocate(res%claims(nclaims),res%claim_node(nclaims));if(nclaims==0)return
        allocate(claimanc(nclaims,max(0,l-1)));pos=0
        do i=1,nnode
            cnt=res%frequency(i)
            do j=1,cnt
                pos=pos+1;res%claim_node(pos)=i
                do k=1,l-1
                    idx=path(l)%value(i,k);claimanc(pos,k)=spar(k)%value(idx)
                end do
            end do
        end do
        if(present(sev_final))then
            if(associated(sev_final%draw))then;call sev_final%draw(claimanc,res%claims)
            else;res%claims=1.0_dp;end if
        else;res%claims=1.0_dp;end if
        do i=1,nclaims;res%aggregate(res%claim_node(i))=res%aggregate(res%claim_node(i))+res%claims(i);end do
    end function rcomphierarc_simulate

    subroutine expand_level(raw,nparent,parent)
        integer,intent(in)::raw(:),nparent
        integer,allocatable,intent(out)::parent(:)
        integer,allocatable::counts(:)
        integer::i,j,n,pos
        if(size(raw)<1 .or. nparent<1)then;allocate(parent(0));return;end if
        allocate(counts(nparent));do i=1,nparent;counts(i)=max(0,raw(1+mod(i-1,size(raw))));end do
        n=sum(counts);allocate(parent(n));pos=0
        do i=1,nparent;do j=1,counts(i);pos=pos+1;parent(pos)=i;end do;end do
    end subroutine expand_level

    subroutine make_ancestors(k,path,fpar,anc)
        integer,intent(in)::k,path(:,:)
        type(real_vector_t),intent(in)::fpar(:)
        real(dp),allocatable,intent(out)::anc(:,:)
        integer::i,j
        allocate(anc(size(path,1),max(0,k-1)))
        if(k<=1)return
        do i=1,size(path,1);do j=1,k-1;anc(i,j)=fpar(j)%value(path(i,j));end do;end do
    end subroutine make_ancestors
end module actuar_hier_sim_v03
