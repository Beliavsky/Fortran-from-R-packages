! SPDX-License-Identifier: GPL-3.0-or-later
!
! Sparse matrix utilities and LDL^T factorization used by the ECOS sparse KKT
! backend.  The LDL symbolic/numeric algorithm follows the algorithmic form of
! SuiteSparse LDL (Timothy A. Davis), as used by upstream ECOS.  See NOTICE.md.
module ecos_sparse
    use, intrinsic :: iso_fortran_env, only : int64
    use ecos_types, only : dp, ecos_csc_matrix, ecos_csr_matrix, ecos_sparse_cache
    implicit none
    private

    type, public :: sparse_triplet_builder
        integer :: nrow = 0
        integer :: ncol = 0
        integer :: nnz = 0
        integer :: capacity = 0
        integer, allocatable :: row(:)
        integer, allocatable :: col(:)
        real(dp), allocatable :: val(:)
    contains
        procedure :: init => triplet_init
        procedure :: add => triplet_add
        procedure :: clear => triplet_clear
    end type sparse_triplet_builder


    type :: integer_neighbor_list
        integer, allocatable :: item(:)
        integer :: n = 0
    end type integer_neighbor_list

    type :: degree_heap
        integer, allocatable :: vertex(:)
        integer, allocatable :: key(:)
        integer :: n = 0
    end type degree_heap

    type, public :: sparse_ldl_factor
        integer :: n = 0
        integer, allocatable :: perm(:)   ! new index -> old index
        integer, allocatable :: iperm(:)  ! old index -> new index
        integer, allocatable :: parent(:)
        integer, allocatable :: lp(:)
        integer, allocatable :: li(:)
        integer, allocatable :: lnz(:)
        integer, allocatable :: flag(:)
        integer, allocatable :: pattern(:)
        real(dp), allocatable :: lx(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: y(:)
        integer :: symbolic_nnz = 0
        logical :: analyzed = .false.
    contains
        procedure :: analyze => ldl_analyze
        procedure :: factorize => ldl_factorize
        procedure :: solve => ldl_solve
        procedure :: solve_refined => ldl_solve_refined
        procedure :: save_symbolic => ldl_save_symbolic
        procedure :: load_symbolic => ldl_load_symbolic
    end type sparse_ldl_factor

    public :: csc_to_csr, csr_to_csc, csr_matvec, csr_tmatvec
    public :: csc_matvec, csc_tmatvec, csr_row_dot
    public :: triplet_to_csc, triplet_to_csr, symmetric_csc_matvec
    public :: reverse_cuthill_mckee, approximate_minimum_degree, sym_permute_upper
    public :: sparse_structure_equal, sparse_pattern_hash

contains

    subroutine triplet_init(this,nrow,ncol,capacity)
        class(sparse_triplet_builder), intent(inout) :: this
        integer, intent(in) :: nrow,ncol
        integer, intent(in), optional :: capacity
        integer :: cap
        cap=64
        if (present(capacity)) cap=max(1,capacity)
        this%nrow=nrow; this%ncol=ncol; this%nnz=0; this%capacity=cap
        if (allocated(this%row)) deallocate(this%row,this%col,this%val)
        allocate(this%row(cap),this%col(cap),this%val(cap))
    end subroutine triplet_init

    subroutine triplet_add(this,i,j,v)
        class(sparse_triplet_builder), intent(inout) :: this
        integer, intent(in) :: i,j
        real(dp), intent(in) :: v
        integer, allocatable :: ir(:),jc(:)
        real(dp), allocatable :: vv(:)
        integer :: newcap
        if (.not.allocated(this%row)) call this%init(max(i,1),max(j,1),64)
        if (i<1 .or. i>this%nrow .or. j<1 .or. j>this%ncol) return
        if (this%nnz>=this%capacity) then
            newcap=max(this%capacity+1,2*this%capacity)
            allocate(ir(newcap),jc(newcap),vv(newcap))
            if (this%nnz>0) then
                ir(1:this%nnz)=this%row(1:this%nnz)
                jc(1:this%nnz)=this%col(1:this%nnz)
                vv(1:this%nnz)=this%val(1:this%nnz)
            end if
            call move_alloc(ir,this%row)
            call move_alloc(jc,this%col)
            call move_alloc(vv,this%val)
            this%capacity=newcap
        end if
        this%nnz=this%nnz+1
        this%row(this%nnz)=i
        this%col(this%nnz)=j
        this%val(this%nnz)=v
    end subroutine triplet_add

    subroutine triplet_clear(this)
        class(sparse_triplet_builder), intent(inout) :: this
        this%nnz=0
    end subroutine triplet_clear

    recursive subroutine quicksort_key(key,row,col,val,lo,hi)
        integer(int64), intent(inout) :: key(:)
        integer, intent(inout) :: row(:),col(:)
        real(dp), intent(inout) :: val(:)
        integer, intent(in) :: lo,hi
        integer :: i,j,tr,tc
        integer(int64) :: pivot,tk
        real(dp) :: tv
        if (lo>=hi) return
        i=lo; j=hi; pivot=key((lo+hi)/2)
        do
            do while(key(i)<pivot); i=i+1; end do
            do while(key(j)>pivot); j=j-1; end do
            if (i<=j) then
                tk=key(i); key(i)=key(j); key(j)=tk
                tr=row(i); row(i)=row(j); row(j)=tr
                tc=col(i); col(i)=col(j); col(j)=tc
                tv=val(i); val(i)=val(j); val(j)=tv
                i=i+1; j=j-1
            end if
            if (i>j) exit
        end do
        if (lo<j) call quicksort_key(key,row,col,val,lo,j)
        if (i<hi) call quicksort_key(key,row,col,val,i,hi)
    end subroutine quicksort_key

    subroutine triplet_to_csc(builder,a,upper_only)
        type(sparse_triplet_builder), intent(in) :: builder
        type(ecos_csc_matrix), intent(out) :: a
        logical, intent(in), optional :: upper_only
        logical :: upper
        integer(int64), allocatable :: key(:)
        integer, allocatable :: rr(:),cc(:),ur(:),uc(:),counts(:)
        real(dp), allocatable :: vv(:),uv(:)
        integer :: k,m,u,j,pos
        upper=.false.; if(present(upper_only)) upper=upper_only
        m=builder%nnz
        allocate(rr(m),cc(m),vv(m),key(m))
        if (m>0) then
            rr=builder%row(1:m); cc=builder%col(1:m); vv=builder%val(1:m)
            if (upper) then
                do k=1,m
                    if (rr(k)>cc(k)) then
                        j=rr(k); rr(k)=cc(k); cc(k)=j
                    end if
                end do
            end if
            do k=1,m
                key(k)=int(cc(k)-1,int64)*int(builder%nrow,int64)+int(rr(k),int64)
            end do
            call quicksort_key(key,rr,cc,vv,1,m)
        end if
        allocate(ur(max(1,m)),uc(max(1,m)),uv(max(1,m)))
        u=0
        k=1
        do while(k<=m)
            u=u+1; ur(u)=rr(k); uc(u)=cc(k); uv(u)=vv(k)
            pos=k+1
            do while(pos<=m)
                if(key(pos)/=key(k)) exit
                uv(u)=uv(u)+vv(pos)
                pos=pos+1
            end do
            k=pos
        end do
        a%nrow=builder%nrow; a%ncol=builder%ncol
        allocate(a%colptr(a%ncol+1),a%rowind(u),a%values(u),counts(a%ncol))
        counts=0
        do k=1,u
            counts(uc(k))=counts(uc(k))+1
        end do
        a%colptr(1)=1
        do j=1,a%ncol
            a%colptr(j+1)=a%colptr(j)+counts(j)
        end do
        if (u>0) then
            a%rowind=ur(1:u); a%values=uv(1:u)
        end if
    end subroutine triplet_to_csc

    subroutine triplet_to_csr(builder,a)
        type(sparse_triplet_builder), intent(in) :: builder
        type(ecos_csr_matrix), intent(out) :: a
        type(sparse_triplet_builder) :: bt
        type(ecos_csc_matrix) :: at
        integer :: k
        call bt%init(builder%ncol,builder%nrow,max(1,builder%nnz))
        do k=1,builder%nnz
            call bt%add(builder%col(k),builder%row(k),builder%val(k))
        end do
        call triplet_to_csc(bt,at)
        a%nrow=builder%nrow; a%ncol=builder%ncol
        allocate(a%rowptr(a%nrow+1),a%colind(size(at%rowind)),a%values(size(at%values)))
        a%rowptr=at%colptr
        a%colind=at%rowind
        a%values=at%values
    end subroutine triplet_to_csr

    subroutine csc_to_csr(a,b)
        type(ecos_csc_matrix), intent(in) :: a
        type(ecos_csr_matrix), intent(out) :: b
        integer, allocatable :: count(:),next(:)
        integer :: j,k,r,p
        b%nrow=a%nrow; b%ncol=a%ncol
        allocate(count(b%nrow),next(b%nrow))
        count=0
        do k=1,size(a%rowind)
            count(a%rowind(k))=count(a%rowind(k))+1
        end do
        allocate(b%rowptr(b%nrow+1))
        b%rowptr(1)=1
        do r=1,b%nrow
            b%rowptr(r+1)=b%rowptr(r)+count(r)
        end do
        allocate(b%colind(size(a%rowind)),b%values(size(a%values)))
        next=b%rowptr(1:b%nrow)
        do j=1,a%ncol
            do k=a%colptr(j),a%colptr(j+1)-1
                r=a%rowind(k); p=next(r)
                b%colind(p)=j; b%values(p)=a%values(k); next(r)=p+1
            end do
        end do
    end subroutine csc_to_csr

    subroutine csr_to_csc(a,b)
        type(ecos_csr_matrix), intent(in) :: a
        type(ecos_csc_matrix), intent(out) :: b
        integer, allocatable :: count(:),next(:)
        integer :: i,k,c,p
        b%nrow=a%nrow; b%ncol=a%ncol
        allocate(count(b%ncol),next(b%ncol)); count=0
        do k=1,size(a%colind); count(a%colind(k))=count(a%colind(k))+1; end do
        allocate(b%colptr(b%ncol+1)); b%colptr(1)=1
        do c=1,b%ncol; b%colptr(c+1)=b%colptr(c)+count(c); end do
        allocate(b%rowind(size(a%colind)),b%values(size(a%values)))
        next=b%colptr(1:b%ncol)
        do i=1,a%nrow
            do k=a%rowptr(i),a%rowptr(i+1)-1
                c=a%colind(k); p=next(c)
                b%rowind(p)=i; b%values(p)=a%values(k); next(c)=p+1
            end do
        end do
    end subroutine csr_to_csc

    subroutine csr_matvec(a,x,y)
        type(ecos_csr_matrix), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: i,k
        y=0.0_dp
        do i=1,a%nrow
            do k=a%rowptr(i),a%rowptr(i+1)-1
                y(i)=y(i)+a%values(k)*x(a%colind(k))
            end do
        end do
    end subroutine csr_matvec

    subroutine csr_tmatvec(a,x,y)
        type(ecos_csr_matrix), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: i,k
        y=0.0_dp
        do i=1,a%nrow
            do k=a%rowptr(i),a%rowptr(i+1)-1
                y(a%colind(k))=y(a%colind(k))+a%values(k)*x(i)
            end do
        end do
    end subroutine csr_tmatvec

    subroutine csc_matvec(a,x,y)
        type(ecos_csc_matrix), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: j,k
        y=0.0_dp
        do j=1,a%ncol
            do k=a%colptr(j),a%colptr(j+1)-1
                y(a%rowind(k))=y(a%rowind(k))+a%values(k)*x(j)
            end do
        end do
    end subroutine csc_matvec

    subroutine csc_tmatvec(a,x,y)
        type(ecos_csc_matrix), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: j,k
        y=0.0_dp
        do j=1,a%ncol
            do k=a%colptr(j),a%colptr(j+1)-1
                y(j)=y(j)+a%values(k)*x(a%rowind(k))
            end do
        end do
    end subroutine csc_tmatvec

    pure real(dp) function csr_row_dot(a,i,b,j) result(v)
        type(ecos_csr_matrix), intent(in) :: a,b
        integer, intent(in) :: i,j
        integer :: pa,pb,ea,eb,ca,cb
        v=0.0_dp
        pa=a%rowptr(i); ea=a%rowptr(i+1)-1
        pb=b%rowptr(j); eb=b%rowptr(j+1)-1
        do while(pa<=ea .and. pb<=eb)
            ca=a%colind(pa); cb=b%colind(pb)
            if(ca==cb) then
                v=v+a%values(pa)*b%values(pb); pa=pa+1; pb=pb+1
            else if(ca<cb) then
                pa=pa+1
            else
                pb=pb+1
            end if
        end do
    end function csr_row_dot

    subroutine symmetric_csc_matvec(a,x,y)
        type(ecos_csc_matrix), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: j,k,i
        y=0.0_dp
        do j=1,a%ncol
            do k=a%colptr(j),a%colptr(j+1)-1
                i=a%rowind(k)
                y(i)=y(i)+a%values(k)*x(j)
                if(i/=j) y(j)=y(j)+a%values(k)*x(i)
            end do
        end do
    end subroutine symmetric_csc_matvec

    logical function sparse_structure_equal(a,b) result(eq)
        type(ecos_csc_matrix), intent(in) :: a,b
        eq=a%nrow==b%nrow .and. a%ncol==b%ncol
        if(.not.eq) return
        if(.not.allocated(a%colptr) .or. .not.allocated(b%colptr)) then
            eq=.false.; return
        end if
        eq=size(a%colptr)==size(b%colptr) .and. size(a%rowind)==size(b%rowind)
        if(.not.eq) return
        eq=all(a%colptr==b%colptr) .and. all(a%rowind==b%rowind)
    end function sparse_structure_equal

    subroutine reverse_cuthill_mckee(a,perm)
        type(ecos_csc_matrix), intent(in) :: a
        integer, allocatable, intent(out) :: perm(:)
        integer, allocatable :: deg(:),ptr(:),adj(:),next(:),queue(:),order(:),nbr(:)
        logical, allocatable :: seen(:)
        integer :: n,j,k,i,p,nnzadj,start,head,tail,no,nn,ii,jj,t
        n=a%ncol
        allocate(deg(n)); deg=0
        do j=1,n
            do k=a%colptr(j),a%colptr(j+1)-1
                i=a%rowind(k)
                if(i/=j) then
                    deg(i)=deg(i)+1; deg(j)=deg(j)+1
                end if
            end do
        end do
        allocate(ptr(n+1)); ptr(1)=1
        do i=1,n; ptr(i+1)=ptr(i)+deg(i); end do
        nnzadj=ptr(n+1)-1
        allocate(adj(max(1,nnzadj)),next(n)); next=ptr(1:n)
        do j=1,n
            do k=a%colptr(j),a%colptr(j+1)-1
                i=a%rowind(k)
                if(i/=j) then
                    adj(next(i))=j; next(i)=next(i)+1
                    adj(next(j))=i; next(j)=next(j)+1
                end if
            end do
        end do
        allocate(seen(n),queue(n),order(n),nbr(max(1,n)),perm(n)); seen=.false.; no=0
        do while(no<n)
            start=0
            do i=1,n
                if(.not.seen(i)) then
                    if(start==0) then
                        start=i
                    else if(deg(i)<deg(start)) then
                        start=i
                    end if
                end if
            end do
            head=1; tail=1; queue(1)=start; seen(start)=.true.
            do while(head<=tail)
                i=queue(head); head=head+1; no=no+1; order(no)=i
                nn=0
                do p=ptr(i),ptr(i+1)-1
                    j=adj(p)
                    if(.not.seen(j)) then
                        nn=nn+1; nbr(nn)=j; seen(j)=.true.
                    end if
                end do
                ! Stable insertion sort of neighbors by degree, then index.
                do ii=2,nn
                    t=nbr(ii); jj=ii-1
                    do while(jj>=1)
                        if(deg(nbr(jj))<deg(t)) exit
                        if(deg(nbr(jj))==deg(t) .and. nbr(jj)<t) exit
                        nbr(jj+1)=nbr(jj); jj=jj-1
                    end do
                    nbr(jj+1)=t
                end do
                if(nn>0) then
                    queue(tail+1:tail+nn)=nbr(1:nn); tail=tail+nn
                end if
            end do
        end do
        do i=1,n
            perm(i)=order(n-i+1)
        end do
    end subroutine reverse_cuthill_mckee


    logical function list_contains(list, value) result(found)
        type(integer_neighbor_list), intent(in) :: list
        integer, intent(in) :: value
        integer :: k
        found = .false.
        do k = 1, list%n
            if (list%item(k) == value) then
                found = .true.
                return
            end if
        end do
    end function list_contains

    subroutine list_append_unique(list, value)
        type(integer_neighbor_list), intent(inout) :: list
        integer, intent(in) :: value
        integer, allocatable :: tmp(:)
        integer :: cap
        if (list_contains(list, value)) return
        if (.not. allocated(list%item)) then
            allocate(list%item(8))
        else if (list%n >= size(list%item)) then
            cap = max(8, 2*size(list%item))
            allocate(tmp(cap))
            if (list%n > 0) tmp(1:list%n) = list%item(1:list%n)
            call move_alloc(tmp, list%item)
        end if
        list%n = list%n + 1
        list%item(list%n) = value
    end subroutine list_append_unique

    subroutine heap_reserve(heap, need)
        type(degree_heap), intent(inout) :: heap
        integer, intent(in) :: need
        integer, allocatable :: iv(:), ik(:)
        integer :: cap
        if (.not. allocated(heap%vertex)) then
            cap = max(64, need)
            allocate(heap%vertex(cap), heap%key(cap))
        else if (need > size(heap%vertex)) then
            cap = max(need, 2*size(heap%vertex))
            allocate(iv(cap), ik(cap))
            if (heap%n > 0) then
                iv(1:heap%n) = heap%vertex(1:heap%n)
                ik(1:heap%n) = heap%key(1:heap%n)
            end if
            call move_alloc(iv, heap%vertex)
            call move_alloc(ik, heap%key)
        end if
    end subroutine heap_reserve

    subroutine heap_push(heap, vertex, key)
        type(degree_heap), intent(inout) :: heap
        integer, intent(in) :: vertex, key
        integer :: i, parent, tv, tk
        call heap_reserve(heap, heap%n + 1)
        heap%n = heap%n + 1
        i = heap%n
        heap%vertex(i) = vertex
        heap%key(i) = key
        do while (i > 1)
            parent = i/2
            if (heap%key(parent) < heap%key(i)) exit
            if (heap%key(parent) == heap%key(i) .and. heap%vertex(parent) <= heap%vertex(i)) exit
            tv = heap%vertex(parent); heap%vertex(parent) = heap%vertex(i); heap%vertex(i) = tv
            tk = heap%key(parent); heap%key(parent) = heap%key(i); heap%key(i) = tk
            i = parent
        end do
    end subroutine heap_push

    subroutine heap_pop(heap, vertex, key, ok)
        type(degree_heap), intent(inout) :: heap
        integer, intent(out) :: vertex, key
        logical, intent(out) :: ok
        integer :: i, left, right, child, tv, tk
        if (heap%n <= 0) then
            vertex = 0; key = 0; ok = .false.
            return
        end if
        vertex = heap%vertex(1); key = heap%key(1); ok = .true.
        heap%vertex(1) = heap%vertex(heap%n)
        heap%key(1) = heap%key(heap%n)
        heap%n = heap%n - 1
        i = 1
        do
            left = 2*i; right = left + 1
            if (left > heap%n) exit
            child = left
            if (right <= heap%n) then
                if (heap%key(right) < heap%key(left)) child = right
                if (heap%key(right) == heap%key(left) .and. &
                    heap%vertex(right) < heap%vertex(left)) child = right
            end if
            if (heap%key(i) < heap%key(child)) exit
            if (heap%key(i) == heap%key(child) .and. heap%vertex(i) <= heap%vertex(child)) exit
            tv = heap%vertex(i); heap%vertex(i) = heap%vertex(child); heap%vertex(child) = tv
            tk = heap%key(i); heap%key(i) = heap%key(child); heap%key(child) = tk
            i = child
        end do
    end subroutine heap_pop

    subroutine approximate_minimum_degree(a, perm)
        ! Native approximate minimum-degree ordering.  The implementation uses
        ! explicit sparse adjacency lists and lazy heap updates.  Small fronts
        ! receive exact clique fill updates; very large fronts use degree-only
        ! approximation to prevent O(d^2) work from dominating ordering time.
        type(ecos_csc_matrix), intent(in) :: a
        integer, allocatable, intent(out) :: perm(:)
        type(integer_neighbor_list), allocatable :: adj(:)
        type(degree_heap) :: heap
        logical, allocatable :: alive(:), touched(:)
        integer, allocatable :: degree(:), nbr(:)
        integer :: n, i, j, k, v, kv, step, nn, p, q, u, w
        logical :: ok
        integer, parameter :: exact_front_limit = 128

        n = a%ncol
        allocate(perm(n), adj(n), alive(n), touched(n), degree(n), nbr(max(1,n)))
        alive = .true.; touched = .false.; degree = 0
        do j = 1, n
            do k = a%colptr(j), a%colptr(j+1)-1
                i = a%rowind(k)
                if (i == j) cycle
                call list_append_unique(adj(i), j)
                call list_append_unique(adj(j), i)
            end do
        end do
        do i = 1, n
            degree(i) = adj(i)%n
            call heap_push(heap, i, degree(i))
        end do

        do step = 1, n
            do
                call heap_pop(heap, v, kv, ok)
                if (.not. ok) then
                    v = 0
                    do i = 1, n
                        if (alive(i)) then
                            v = i
                            exit
                        end if
                    end do
                    exit
                end if
                if (alive(v) .and. kv == degree(v)) exit
            end do
            if (v == 0) exit
            perm(step) = v
            nn = 0
            do k = 1, adj(v)%n
                u = adj(v)%item(k)
                if (alive(u)) then
                    nn = nn + 1
                    nbr(nn) = u
                end if
            end do
            alive(v) = .false.
            touched = .false.
            do p = 1, nn
                u = nbr(p)
                degree(u) = max(0, degree(u)-1)
                touched(u) = .true.
            end do
            if (nn <= exact_front_limit) then
                do p = 1, nn-1
                    u = nbr(p)
                    do q = p+1, nn
                        w = nbr(q)
                        if (.not. list_contains(adj(u), w)) then
                            call list_append_unique(adj(u), w)
                            call list_append_unique(adj(w), u)
                            degree(u) = degree(u) + 1
                            degree(w) = degree(w) + 1
                            touched(u) = .true.; touched(w) = .true.
                        end if
                    end do
                end do
            end if
            do p = 1, nn
                u = nbr(p)
                if (alive(u) .and. touched(u)) call heap_push(heap, u, degree(u))
            end do
        end do
    end subroutine approximate_minimum_degree

    subroutine sym_permute_upper(a,perm,b)
        type(ecos_csc_matrix), intent(in) :: a
        integer, intent(in) :: perm(:)
        type(ecos_csc_matrix), intent(out) :: b
        integer, allocatable :: ip(:)
        type(sparse_triplet_builder) :: tb
        integer :: n,j,k,i,ii,jj
        n=a%ncol; allocate(ip(n))
        do i=1,n; ip(perm(i))=i; end do
        call tb%init(n,n,max(1,size(a%values)))
        do j=1,n
            do k=a%colptr(j),a%colptr(j+1)-1
                i=a%rowind(k); ii=ip(i); jj=ip(j)
                call tb%add(min(ii,jj),max(ii,jj),a%values(k))
            end do
        end do
        call triplet_to_csc(tb,b,.true.)
    end subroutine sym_permute_upper

    integer function sparse_pattern_hash(a) result(h)
        type(ecos_csc_matrix), intent(in) :: a
        integer(int64) :: z
        integer :: k
        z = int(1469598103934665603_int64, int64)
        z = ieor(z, int(a%nrow, int64)); z = z*1099511628211_int64
        z = ieor(z, int(a%ncol, int64)); z = z*1099511628211_int64
        do k = 1, size(a%colptr)
            z = ieor(z, int(a%colptr(k), int64)); z = z*1099511628211_int64
        end do
        do k = 1, size(a%rowind)
            z = ieor(z, int(a%rowind(k), int64)); z = z*1099511628211_int64
        end do
        h = int(iand(z, int(z'7fffffff', int64)))
    end function sparse_pattern_hash

    subroutine ldl_save_symbolic(this, cache, pattern_hash)
        class(sparse_ldl_factor), intent(in) :: this
        type(ecos_sparse_cache), intent(inout) :: cache
        integer, intent(in) :: pattern_hash
        if (.not. this%analyzed) return
        cache%n = this%n
        cache%structure_hash = pattern_hash
        cache%symbolic_nnz = this%symbolic_nnz
        if (allocated(cache%perm)) deallocate(cache%perm)
        if (allocated(cache%parent)) deallocate(cache%parent)
        if (allocated(cache%lp)) deallocate(cache%lp)
        allocate(cache%perm(this%n), cache%parent(this%n), cache%lp(this%n+1))
        cache%perm = this%perm
        cache%parent = this%parent
        cache%lp = this%lp
        cache%symbolic_valid = .true.
    end subroutine ldl_save_symbolic

    subroutine ldl_load_symbolic(this, cache, info)
        class(sparse_ldl_factor), intent(inout) :: this
        type(ecos_sparse_cache), intent(in) :: cache
        integer, intent(out) :: info
        integer :: n, i, total
        info = 0
        if (.not. cache%symbolic_valid) then
            info = 1
            return
        end if
        n = cache%n
        if (n < 1 .or. .not. allocated(cache%perm) .or. .not. allocated(cache%parent) .or. &
            .not. allocated(cache%lp)) then
            info = 2
            return
        end if
        this%n = n
        if (allocated(this%perm)) deallocate(this%perm,this%iperm,this%parent,this%lp,this%lnz, &
            this%flag,this%pattern,this%y,this%d)
        if (allocated(this%li)) deallocate(this%li,this%lx)
        allocate(this%perm(n),this%iperm(n),this%parent(n),this%lp(n+1),this%lnz(n), &
            this%flag(n),this%pattern(n),this%y(n),this%d(n))
        this%perm = cache%perm
        do i = 1, n
            this%iperm(this%perm(i)) = i
        end do
        this%parent = cache%parent
        this%lp = cache%lp
        total = cache%symbolic_nnz
        allocate(this%li(max(1,total)),this%lx(max(1,total)))
        this%symbolic_nnz = total
        this%analyzed = .true.
    end subroutine ldl_load_symbolic

    subroutine ldl_analyze(this,a,use_rcm,info,use_amd)
        class(sparse_ldl_factor), intent(inout) :: this
        type(ecos_csc_matrix), intent(in) :: a
        logical, intent(in), optional :: use_rcm
        integer, intent(out) :: info
        logical, intent(in), optional :: use_amd
        logical :: do_rcm, do_amd
        type(ecos_csc_matrix) :: ap
        integer :: n,k,p,i,total
        do_rcm=.false.; if(present(use_rcm)) do_rcm=use_rcm
        do_amd=.true.; if(present(use_amd)) do_amd=use_amd
        info=0
        if(a%nrow/=a%ncol .or. a%ncol<1) then; info=1; return; end if
        n=a%ncol; this%n=n
        if(allocated(this%perm)) deallocate(this%perm,this%iperm,this%parent,this%lp,this%lnz, &
            this%flag,this%pattern,this%y,this%d)
        if(allocated(this%li)) deallocate(this%li,this%lx)
        allocate(this%perm(n),this%iperm(n),this%parent(n),this%lp(n+1),this%lnz(n), &
            this%flag(n),this%pattern(n),this%y(n),this%d(n))
        if(do_amd) then
            call approximate_minimum_degree(a,this%perm)
        else if(do_rcm) then
            call reverse_cuthill_mckee(a,this%perm)
        else
            this%perm=[(i,i=1,n)]
        end if
        do i=1,n; this%iperm(this%perm(i))=i; end do
        call sym_permute_upper(a,this%perm,ap)
        this%parent=0; this%flag=0; this%lnz=0
        do k=1,n
            this%flag(k)=k
            do p=ap%colptr(k),ap%colptr(k+1)-1
                i=ap%rowind(p)
                if(i>=k) cycle
                do while(this%flag(i)/=k)
                    if(this%parent(i)==0) this%parent(i)=k
                    this%lnz(i)=this%lnz(i)+1
                    this%flag(i)=k
                    i=this%parent(i)
                end do
            end do
        end do
        this%lp(1)=1
        do k=1,n; this%lp(k+1)=this%lp(k)+this%lnz(k); end do
        total=this%lp(n+1)-1
        allocate(this%li(max(1,total)),this%lx(max(1,total)))
        this%symbolic_nnz=total
        this%analyzed=.true.
    end subroutine ldl_analyze

    subroutine ldl_factorize(this,a,sign,delta,eps,info)
        class(sparse_ldl_factor), intent(inout) :: this
        type(ecos_csc_matrix), intent(in) :: a
        real(dp), intent(in) :: sign(:),delta,eps
        integer, intent(out) :: info
        type(ecos_csc_matrix) :: ap
        integer :: n,k,p,i,len,top,p2,pos
        integer, allocatable :: path(:)
        real(dp) :: yi,lki,sgn
        info=0
        if(.not.this%analyzed) then; info=1; return; end if
        n=this%n
        if(size(sign)/=n) then; info=2; return; end if
        call sym_permute_upper(a,this%perm,ap)
        allocate(path(n))
        this%y=0.0_dp; this%flag=0; this%lnz=0; this%lx=0.0_dp; this%li=0; this%d=0.0_dp
        do k=1,n
            this%y(k)=0.0_dp; top=n+1; this%flag(k)=k
            do p=ap%colptr(k),ap%colptr(k+1)-1
                i=ap%rowind(p)
                if(i>k) cycle
                this%y(i)=this%y(i)+ap%values(p)
                len=0
                do while(this%flag(i)/=k)
                    len=len+1
                    path(len)=i
                    this%flag(i)=k
                    i=this%parent(i)
                    if(i==0) exit
                end do
                do while(len>0)
                    top=top-1
                    this%pattern(top)=path(len)
                    len=len-1
                end do
            end do
            this%d(k)=this%y(k); this%y(k)=0.0_dp
            do while(top<=n)
                i=this%pattern(top); top=top+1
                yi=this%y(i); this%y(i)=0.0_dp
                p2=this%lp(i)+this%lnz(i)-1
                do p=this%lp(i),p2
                    this%y(this%li(p))=this%y(this%li(p))-this%lx(p)*yi
                end do
                if(abs(this%d(i))<=tiny(1.0_dp)) then; info=3; return; end if
                lki=yi/this%d(i)
                this%d(k)=this%d(k)-lki*yi
                pos=this%lp(i)+this%lnz(i)
                if(pos>this%lp(i+1)-1) then; info=4; return; end if
                this%li(pos)=k; this%lx(pos)=lki; this%lnz(i)=this%lnz(i)+1
            end do
            sgn=sign(this%perm(k))
            if(sgn*this%d(k)<=eps) this%d(k)=sgn*max(delta,eps)
        end do
    end subroutine ldl_factorize

    subroutine ldl_solve(this,b,x,info)
        class(sparse_ldl_factor), intent(in) :: this
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: xp(:)
        integer :: n,i,j,p
        info=0; n=this%n
        if(.not.this%analyzed .or. size(b)/=n .or. size(x)/=n) then; info=1; return; end if
        allocate(xp(n))
        do i=1,n; xp(i)=b(this%perm(i)); end do
        do j=1,n
            do p=this%lp(j),this%lp(j+1)-1
                xp(this%li(p))=xp(this%li(p))-this%lx(p)*xp(j)
            end do
        end do
        do i=1,n
            if(abs(this%d(i))<=tiny(1.0_dp)) then; info=2; return; end if
            xp(i)=xp(i)/this%d(i)
        end do
        do j=n,1,-1
            do p=this%lp(j),this%lp(j+1)-1
                xp(j)=xp(j)-this%lx(p)*xp(this%li(p))
            end do
        end do
        x=0.0_dp
        do i=1,n; x(this%perm(i))=xp(i); end do
    end subroutine ldl_solve

    subroutine ldl_solve_refined(this,a,b,x,max_refine,tol,nref,info)
        class(sparse_ldl_factor), intent(in) :: this
        type(ecos_csc_matrix), intent(in) :: a
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: max_refine
        real(dp), intent(in) :: tol
        integer, intent(out) :: nref,info
        real(dp), allocatable :: ax(:),r(:),dx(:)
        real(dp) :: nr,nprev,nb
        integer :: k,ie
        allocate(ax(size(b)),r(size(b)),dx(size(b)))
        call this%solve(b,x,info); if(info/=0) return
        nb=max(1.0_dp,sqrt(dot_product(b,b))); nprev=huge(1.0_dp); nref=0
        do k=1,max_refine
            call symmetric_csc_matvec(a,x,ax)
            r=b-ax; nr=sqrt(dot_product(r,r))/nb
            if(nr<=tol .or. nr>=nprev) exit
            nprev=nr
            call this%solve(r,dx,ie); if(ie/=0) exit
            x=x+dx; nref=nref+1
        end do
    end subroutine ldl_solve_refined

end module ecos_sparse
