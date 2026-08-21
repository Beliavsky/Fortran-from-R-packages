module boot_resampling
    use boot_kinds, only : dp
    implicit none
    private
    public :: ordinary_array, balanced_array, permutation_array, antithetic_array
    public :: importance_array, balanced_importance_array, frequency_array
    public :: normalize_weights, random_permutation, weighted_sample

contains

    subroutine random_permutation(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, tmp
        real(dp) :: u
        do i = size(x), 2, -1
            call random_number(u)
            j = 1 + int(u*real(i,dp))
            if (j > i) j = i
            tmp = x(i)
            x(i) = x(j)
            x(j) = tmp
        end do
    end subroutine random_permutation

    integer function weighted_sample(items, prob) result(item)
        integer, intent(in) :: items(:)
        real(dp), intent(in) :: prob(:)
        real(dp) :: u, s, total
        integer :: i
        total = sum(max(prob,0.0_dp))
        if (total <= 0.0_dp) error stop "weighted_sample: nonpositive probability sum"
        call random_number(u)
        u = u*total
        s = 0.0_dp
        do i=1,size(items)
            s = s + max(prob(i),0.0_dp)
            if (u <= s) then
                item = items(i)
                return
            end if
        end do
        item = items(size(items))
    end function weighted_sample

    subroutine normalize_weights(w, strata, out)
        real(dp), intent(in) :: w(:)
        integer, intent(in) :: strata(:)
        real(dp), intent(out) :: out(:)
        integer :: i, j
        real(dp) :: s
        if (size(w) /= size(strata) .or. size(out) /= size(w)) error stop "normalize_weights: size mismatch"
        do i=1,size(w)
            s = 0.0_dp
            do j=1,size(w)
                if (strata(j) == strata(i)) s = s + w(j)
            end do
            if (s <= 0.0_dp) error stop "normalize_weights: nonpositive stratum sum"
            out(i) = w(i)/s
        end do
    end subroutine normalize_weights

    subroutine ordinary_array(n, r, strata, idx)
        integer, intent(in) :: n, r
        integer, intent(in) :: strata(n)
        integer, intent(out) :: idx(r,n)
        integer :: i, j, k, ng, g(n)
        real(dp) :: u
        do i=1,n
            ng = 0
            do j=1,n
                if (strata(j) == strata(i)) then
                    ng = ng + 1
                    g(ng)=j
                end if
            end do
            do k=1,r
                if (ng == 1) then
                    idx(k,i)=g(1)
                else
                    call random_number(u)
                    j = 1 + int(u*real(ng,dp))
                    if (j>ng) j=ng
                    idx(k,i)=g(j)
                end if
            end do
        end do
    end subroutine ordinary_array

    subroutine permutation_array(n, r, strata, idx)
        integer, intent(in) :: n, r, strata(n)
        integer, intent(out) :: idx(r,n)
        integer :: i, j, k, ng, g(n), p(n)
        logical :: done(n)
        done=.false.
        do i=1,n
            if (done(i)) cycle
            ng=0
            do j=1,n
                if (strata(j)==strata(i)) then
                    ng=ng+1
                    g(ng)=j
                    done(j)=.true.
                end if
            end do
            do k=1,r
                p(1:ng)=g(1:ng)
                call random_permutation(p(1:ng))
                do j=1,ng
                    idx(k,g(j))=p(j)
                end do
            end do
        end do
    end subroutine permutation_array

    subroutine balanced_array(n, r, strata, idx)
        integer, intent(in) :: n, r, strata(n)
        integer, intent(out) :: idx(r,n)
        integer :: i,j,k,ng,g(n), pool(n*r), pos
        logical :: done(n)
        done=.false.
        idx=0
        do i=1,n
            if(done(i)) cycle
            ng=0
            do j=1,n
                if(strata(j)==strata(i)) then
                    ng=ng+1
                    g(ng)=j
                    done(j)=.true.
                end if
            end do
            pos=0
            do k=1,r
                do j=1,ng
                    pos=pos+1
                    pool(pos)=g(j)
                end do
            end do
            call random_permutation(pool(1:pos))
            pos=0
            do k=1,r
                do j=1,ng
                    pos=pos+1
                    idx(k,g(j))=pool(pos)
                end do
            end do
        end do
    end subroutine balanced_array

    subroutine antithetic_array(n, r, influence, strata, idx)
        integer, intent(in) :: n,r,strata(n)
        real(dp), intent(in) :: influence(n)
        integer, intent(out) :: idx(r,n)
        integer :: i,j,k,ng,g(n), rank(n), rev(n), half, sample(n), pick, ii
        logical :: done(n)
        real(dp) :: u
        done=.false.
        idx=0
        half=r/2
        do i=1,n
            if(done(i)) cycle
            ng=0
            do j=1,n
                if(strata(j)==strata(i)) then
                    ng=ng+1
                    g(ng)=j
                    done(j)=.true.
                end if
            end do
            call ranks_unique(influence(g(1:ng)), rank(1:ng))
            do j=1,ng
                do k=1,ng
                    if(rank(k)==ng+1-rank(j)) then
                        rev(j)=g(k)
                        exit
                    end if
                end do
            end do
            do k=1,half
                do j=1,ng
                    call random_number(u)
                    pick = 1 + int(u*real(ng,dp))
                    if(pick>ng) pick=ng
                    sample(j)=g(pick)
                    idx(k,g(j))=sample(j)
                end do
                do j=1,ng
                    do ii=1,ng
                        if(g(ii)==sample(j)) then
                            idx(k+half,g(j))=rev(ii)
                            exit
                        end if
                    end do
                end do
            end do
            if (2*half < r) then
                do j=1,ng
                    call random_number(u)
                    k=1+int(u*real(ng,dp))
                    if(k>ng)k=ng
                    idx(r,g(j))=g(k)
                end do
            end if
        end do
    end subroutine antithetic_array

    subroutine ranks_unique(x, rank)
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: rank(size(x))
        integer :: i,j,n,ord(size(x)),tmp
        real(dp) :: u
        n=size(x)
        ord=[(i,i=1,n)]
        do i=1,n-1
            do j=i+1,n
                if(x(ord(j)) < x(ord(i))) then
                    tmp=ord(i)
                    ord(i)=ord(j)
                    ord(j)=tmp
                else if(abs(x(ord(j))-x(ord(i))) <= epsilon(1.0_dp)*max(1.0_dp,abs(x(ord(i))))) then
                    call random_number(u)
                    if(u<0.5_dp) then
                        tmp=ord(i)
                        ord(i)=ord(j)
                        ord(j)=tmp
                    end if
                end if
            end do
        end do
        do i=1,n
        rank(ord(i))=i
        end do
    end subroutine ranks_unique

    subroutine importance_array(n, r, weights, strata, idx)
        integer, intent(in) :: n,r,strata(n)
        real(dp), intent(in) :: weights(n)
        integer, intent(out) :: idx(r,n)
        integer :: i,j,k,ng,g(n)
        real(dp) :: p(n)
        do i=1,n
            ng=0
            do j=1,n
                if(strata(j)==strata(i)) then
                    ng=ng+1
                    g(ng)=j
                    p(ng)=weights(j)
                end if
            end do
            do k=1,r
                idx(k,i)=weighted_sample(g(1:ng),p(1:ng))
            end do
        end do
    end subroutine importance_array

    subroutine balanced_importance_array(n,r,weights,strata,idx)
        integer,intent(in)::n,r,strata(n)
        real(dp),intent(in)::weights(n)
        integer,intent(out)::idx(r,n)
        integer::i,j,k,ng,g(n),cnt(n),pool(n*r),pos,need,chosen
        real(dp)::p(n), frac(n), s
        logical::done(n)
        done=.false.
        idx=0
        do i=1,n
            if(done(i))cycle
            ng=0
            s=0.0_dp
            do j=1,n
                if(strata(j)==strata(i))then
                    ng=ng+1
                    g(ng)=j
                    p(ng)=max(weights(j),0.0_dp)
                    s=s+p(ng)
                    done(j)=.true.
                end if
            end do
            p(1:ng)=p(1:ng)/s
            cnt(1:ng)=floor(real(ng*r,dp)*p(1:ng))
            frac(1:ng)=real(ng*r,dp)*p(1:ng)-real(cnt(1:ng),dp)
            need=ng*r-sum(cnt(1:ng))
            do k=1,need
                chosen=weighted_sample([(j,j=1,ng)],frac(1:ng))
                cnt(chosen)=cnt(chosen)+1
                frac(chosen)=0.0_dp
            end do
            pos=0
            do j=1,ng
                do k=1,cnt(j)
                pos=pos+1
                pool(pos)=g(j)
                end do
            end do
            call random_permutation(pool(1:pos))
            pos=0
            do k=1,r
                do j=1,ng
                pos=pos+1
                idx(k,g(j))=pool(pos)
                end do
            end do
        end do
    end subroutine balanced_importance_array

    subroutine frequency_array(idx, freq)
        integer,intent(in)::idx(:,:)
        integer,intent(out)::freq(size(idx,1),size(idx,2))
        integer::r,j,n
        n=size(idx,2)
        freq=0
        do r=1,size(idx,1)
            do j=1,n
                if(idx(r,j)>=1 .and. idx(r,j)<=n) freq(r,idx(r,j))=freq(r,idx(r,j))+1
            end do
        end do
    end subroutine frequency_array
end module boot_resampling
