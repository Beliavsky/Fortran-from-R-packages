module mstate_transitions
    use mstate_types, only : transition_map
    implicit none
    private
    public :: transition_from_matrix, trans_comprisk, trans_illdeath
    public :: trans2q, absorbing_states, is_circular, enumerate_paths
    public :: trans2tra, tra2trans

contains

    subroutine trans2tra(tr,tra)
        type(transition_map),intent(in)::tr
        logical,allocatable,intent(out)::tra(:,:)
        allocate(tra(tr%nstate,tr%nstate));tra=tr%trans>0
    end subroutine trans2tra

    subroutine tra2trans(tra,tr,info)
        logical,intent(in)::tra(:,:)
        type(transition_map),intent(out)::tr
        integer,intent(out),optional::info
        integer,allocatable::mat(:,:)
        integer::i,j,val,inf0
        if(present(info))info=0
        if(size(tra,1)/=size(tra,2))then;if(present(info))info=1;return;end if
        allocate(mat(size(tra,1),size(tra,2)));mat=0;val=0
        do i=1,size(tra,1);do j=1,size(tra,2)
            if(tra(i,j))then;val=val+1;mat(i,j)=val;end if
        end do;end do
        call transition_from_matrix(mat,tr,inf0)
        if(inf0/=0.and.present(info))info=10+inf0
    end subroutine tra2trans

    subroutine transition_from_matrix(mat, tr, info)
        integer, intent(in) :: mat(:, :)
        type(transition_map), intent(out) :: tr
        integer, intent(out), optional :: info
        integer :: i, j, k, mx, n
        logical, allocatable :: seen(:)
        if (present(info)) info = 0
        n = size(mat, 1)
        if (size(mat, 2) /= n) then
            if (present(info)) info = 1
            return
        end if
        mx = maxval(mat)
        if (mx < 0) mx = 0
        tr%nstate = n
        tr%ntrans = mx
        allocate(tr%trans(n,n)); tr%trans = mat
        allocate(tr%from(mx), tr%to(mx)); tr%from = 0; tr%to = 0
        if (mx > 0) then
            allocate(seen(mx)); seen = .false.
            do i = 1, n
                do j = 1, n
                    k = mat(i,j)
                    if (k > 0) then
                        if (k > mx .or. seen(k)) then
                            if (present(info)) info = 2
                            return
                        end if
                        seen(k) = .true.; tr%from(k) = i; tr%to(k) = j
                    end if
                end do
            end do
            if (.not. all(seen)) then
                if (present(info)) info = 3
                return
            end if
        end if
    end subroutine transition_from_matrix

    subroutine trans_comprisk(k, tr)
        integer, intent(in) :: k
        type(transition_map), intent(out) :: tr
        integer, allocatable :: mat(:, :)
        integer :: j, info
        allocate(mat(k+1,k+1)); mat = 0
        do j=1,k
            mat(1,j+1)=j
        end do
        call transition_from_matrix(mat,tr,info)
    end subroutine trans_comprisk

    subroutine trans_illdeath(tr)
        type(transition_map), intent(out) :: tr
        integer :: mat(3,3), info
        mat = 0
        mat(1,2)=1; mat(1,3)=2; mat(2,3)=3
        call transition_from_matrix(mat,tr,info)
    end subroutine trans_illdeath

    subroutine trans2q(tr, q)
        type(transition_map), intent(in) :: tr
        integer, allocatable, intent(out) :: q(:, :)
        integer, allocatable :: p(:, :), pk(:, :), prev(:, :), tmp(:, :)
        integer :: n, k
        n=tr%nstate
        allocate(p(n,n),pk(n,n),prev(n,n),tmp(n,n),q(n,n))
        p=merge(1,0,tr%trans>0)
        do k=1,n; p(k,k)=1; end do
        pk=p
        do k=1,n; pk(k,k)=0; end do
        prev=pk; q=pk
        do k=2,n
            tmp=matmul(pk,p)
            where(tmp>1) tmp=1
            q=q+k*(tmp-prev)
            prev=tmp; pk=tmp
        end do
    end subroutine trans2q

    subroutine absorbing_states(tr, states)
        type(transition_map), intent(in) :: tr
        integer, allocatable, intent(out) :: states(:)
        integer :: i, n
        n=count([(all(tr%trans(i,:)==0), i=1,tr%nstate)])
        allocate(states(n)); n=0
        do i=1,tr%nstate
            if (all(tr%trans(i,:)==0)) then
                n=n+1; states(n)=i
            end if
        end do
    end subroutine absorbing_states

    logical function is_circular(tr) result(ans)
        type(transition_map), intent(in) :: tr
        integer, allocatable :: q(:, :)
        integer :: i
        call trans2q(tr,q)
        ans=.false.
        do i=1,tr%nstate
            if (q(i,i)>0) ans=.true.
        end do
    end function is_circular

    subroutine enumerate_paths(tr, start, paths, lengths, info)
        type(transition_map), intent(in) :: tr
        integer, intent(in) :: start
        integer, allocatable, intent(out) :: paths(:, :), lengths(:)
        integer, intent(out), optional :: info
        integer, allocatable :: work(:, :), lens(:), cur(:)
        integer :: maxpaths, npath
        if (present(info)) info=0
        if (start<1 .or. start>tr%nstate .or. is_circular(tr)) then
            allocate(paths(0,0),lengths(0))
            if (present(info)) info=1
            return
        end if
        maxpaths=max(1,2**min(tr%nstate,20))
        allocate(work(maxpaths,tr%nstate),lens(maxpaths),cur(tr%nstate))
        work=0; lens=0; cur=0; npath=0
        call dfs(start,1)
        allocate(paths(npath,tr%nstate),lengths(npath))
        if (npath>0) then
            paths=work(1:npath,:); lengths=lens(1:npath)
        end if
    contains
        recursive subroutine dfs(state, depth)
            integer,intent(in)::state,depth
            integer::j
            cur(depth)=state
            ! mstate::paths() includes the current prefix itself before recursively
            ! appending all extensions from each reachable next state.
            npath=npath+1
            if (npath>maxpaths) return
            work(npath,:)=0; work(npath,1:depth)=cur(1:depth); lens(npath)=depth
            do j=1,tr%nstate
                if (tr%trans(state,j)>0) call dfs(j,depth+1)
            end do
        end subroutine dfs
    end subroutine enumerate_paths
end module mstate_transitions
