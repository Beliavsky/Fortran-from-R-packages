module matchingmarkets_utils
   use matchingmarkets_kinds, only : dp
   use matchingmarkets_types, only : assignment_result_t
   implicit none
   private
   public :: rank_matrix, rank_of, pairs_from_assignment, singles_from_assignment
   public :: college_lists_from_assignment, sort_unique_pairs, normal_cdf, normal_pdf
contains
   function rank_matrix(pref, missing_value) result(rank)
      integer, intent(in) :: pref(:,:)
      integer, intent(in), optional :: missing_value
      integer, allocatable :: rank(:,:)
      integer :: i, j, v, miss, mx
      miss = 0
      if (present(missing_value)) miss = missing_value
      mx = max(maxval(pref), size(pref,1), size(pref,2))
      allocate(rank(mx,size(pref,2)))
      rank = huge(1)
      do j = 1, size(pref,2)
         do i = 1, size(pref,1)
            v = pref(i,j)
            if (v == miss .or. v < 1 .or. v > mx) cycle
            rank(v,j) = i
         end do
      end do
   end function rank_matrix

   integer function rank_of(pref_col, value) result(r)
      integer, intent(in) :: pref_col(:)
      integer, intent(in) :: value
      integer :: i
      r = huge(1)
      do i = 1, size(pref_col)
         if (pref_col(i) == value) then
            r = i
            return
         end if
      end do
   end function rank_of

   function pairs_from_assignment(assign) result(pairs)
      integer, intent(in) :: assign(:)
      integer, allocatable :: pairs(:,:)
      integer :: i, k
      allocate(pairs(2,count(assign>0)))
      k = 0
      do i = 1, size(assign)
         if (assign(i)>0) then
            k = k + 1
            pairs(:,k) = [i, assign(i)]
         end if
      end do
   end function pairs_from_assignment

   function singles_from_assignment(assign) result(singles)
      integer, intent(in) :: assign(:)
      integer, allocatable :: singles(:)
      integer :: i,k
      allocate(singles(count(assign==0)))
      k=0
      do i=1,size(assign)
         if(assign(i)==0) then
            k=k+1; singles(k)=i
         end if
      end do
   end function singles_from_assignment

   subroutine college_lists_from_assignment(assign,nc,slots,lists)
      integer,intent(in)::assign(:),nc,slots(:)
      integer,allocatable,intent(out)::lists(:,:)
      integer :: c,s,k
      allocate(lists(max(1,maxval(slots)),nc)); lists=0
      do s=1,size(assign)
         c=assign(s)
         if(c<1 .or. c>nc) cycle
         k=count(lists(:,c)>0)+1
         if(k<=slots(c)) lists(k,c)=s
      end do
   end subroutine college_lists_from_assignment

   subroutine sort_unique_pairs(pairs)
      integer, allocatable, intent(inout) :: pairs(:,:)
      integer, allocatable :: tmp(:,:)
      integer :: i,j,k,n,a,b
      if(.not.allocated(pairs)) return
      n=size(pairs,2)
      do i=1,n-1
         do j=i+1,n
            if(pairs(1,j)<pairs(1,i) .or. &
               (pairs(1,j)==pairs(1,i) .and. pairs(2,j)<pairs(2,i))) then
               a=pairs(1,i); b=pairs(2,i)
               pairs(:,i)=pairs(:,j); pairs(:,j)=[a,b]
            end if
         end do
      end do
      if(n==0) return
      allocate(tmp(2,n)); k=1; tmp(:,1)=pairs(:,1)
      do i=2,n
         if(any(pairs(:,i)/=tmp(:,k))) then
            k=k+1;tmp(:,k)=pairs(:,i)
         end if
      end do
      pairs=tmp(:,:k)
   end subroutine sort_unique_pairs

   pure real(dp) function normal_cdf(x) result(p)
      real(dp),intent(in)::x
      p=0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function normal_pdf(x) result(p)
      real(dp),intent(in)::x
      p=exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp))
   end function normal_pdf
end module matchingmarkets_utils
