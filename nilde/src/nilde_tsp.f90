! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_tsp
   use nilde_kinds, only : i8, dp
   use nilde_types, only : tsp_result_t
   use nilde_collect, only : append_int_column
   implicit none
   private
   public :: tsp_solver, assignment_lower_bound
contains

   function tsp_solver(cost, upper_bound, lower_bound, no_go, method) result(res)
      integer(i8), intent(in) :: cost(:,:)
      integer(i8), intent(in), optional :: upper_bound, lower_bound, no_go
      character(len=*), intent(in), optional :: method
      type(tsp_result_t) :: res
      integer :: n, i, j, nt
      integer(i8) :: ng, lb, ub, target, maxc
      integer(i8), allocatable :: c(:,:)
      character(len=:), allocatable :: meth
      integer, allocatable :: succ(:), tour(:), tours_tmp(:,:)
      logical, allocatable :: used_in(:)

      n=size(cost,1)
      if (n < 2 .or. size(cost,2) /= n) error stop "tsp_solver: square cost matrix required"
      maxc=0_i8
      do j=1,n; do i=1,n
         if (i/=j .and. cost(i,j)>=0_i8) maxc=max(maxc,cost(i,j))
      end do; end do
      if (maxc <= 0_i8) return
      ng=maxc*100000_i8; if (present(no_go)) ng=no_go
      allocate(c(n,n)); c=cost
      do i=1,n
         c(i,i)=ng
         do j=1,n
            if (i/=j .and. c(i,j)<0_i8) c(i,j)=ng
         end do
      end do
      if (present(lower_bound)) then
         lb=lower_bound
      else
         lb=assignment_lower_bound(c)
      end if
      if (present(upper_bound)) then
         ub=upper_bound
      else
         meth='cheapest_insertion'; if (present(method)) meth=trim(method)
         select case(meth)
         case('nearest_neighbor')
            ub=nearest_neighbor_upper_bound(c,ng)
         case default
            ub=cheapest_insertion_upper_bound(c,ng)
         end select
      end if
      if (ub < lb .or. ub >= ng) return
      res%initial_lower_bound=lb; res%final_lower_bound=lb; res%upper_bound=ub
      allocate(succ(n),used_in(n),tour(n)); succ=0; used_in=.false.; nt=0

      do target=lb,ub
         res%iterations=res%iterations+1
         call search_assign(1,0_i8,target)
         if (nt > 0) then
            res%tour_length=target; res%final_lower_bound=target; exit
         end if
      end do
      if (nt > 0) then
         res%ntours=nt; allocate(res%tours(n,nt)); res%tours=tours_tmp(:,1:nt)
      end if

   contains
      recursive subroutine search_assign(row, partial, target_cost)
         integer, intent(in) :: row
         integer(i8), intent(in) :: partial, target_cost
         integer :: col, krow
         integer(i8) :: newcost, minrem
         if (partial > target_cost) return
         if (row > n) then
            res%edge_subsets_tested=res%edge_subsets_tested+1_i8
            if (partial /= target_cost) return
            res%degree_feasible=res%degree_feasible+1_i8
            if (one_cycle(succ,tour)) call append_unique_tour(tour)
            return
         end if
         ! Cheap optimistic bound from the cheapest still-available incoming edge per remaining row.
         minrem=0_i8
         do krow=row,n
            minrem=minrem+row_min_available(krow)
            if (minrem >= ng) return
         end do
         if (partial+minrem > target_cost) return
         do col=1,n
            if (col==row .or. used_in(col) .or. c(row,col)>=ng) cycle
            newcost=partial+c(row,col)
            if (newcost > target_cost) cycle
            succ(row)=col; used_in(col)=.true.
            call search_assign(row+1,newcost,target_cost)
            used_in(col)=.false.; succ(row)=0
         end do
      end subroutine search_assign

      integer(i8) function row_min_available(row) result(v)
         integer, intent(in) :: row
         integer :: q
         v=ng
         do q=1,n
            if (q/=row .and. .not. used_in(q)) v=min(v,c(row,q))
         end do
      end function row_min_available

      logical function one_cycle(s, t) result(ok)
         integer, intent(in) :: s(:)
         integer, intent(out) :: t(:)
         logical :: seen(size(s))
         integer :: k, v
         seen=.false.; v=1
         do k=1,size(s)
            if (v<1 .or. v>size(s) .or. seen(v)) then; ok=.false.; return; end if
            t(k)=v; seen(v)=.true.; v=s(v)
         end do
         ok=(v==1 .and. all(seen))
      end function one_cycle

      subroutine append_unique_tour(t)
         integer, intent(in) :: t(:)
         integer :: q
         do q=1,nt
            if (all(tours_tmp(:,q)==t)) return
         end do
         call append_int_column(tours_tmp,nt,t)
      end subroutine append_unique_tour
   end function tsp_solver

   function assignment_lower_bound(cost) result(obj)
      integer(i8), intent(in) :: cost(:,:)
      integer(i8) :: obj
      integer :: n, i, j, j0, j1, i0
      real(dp), allocatable :: u(:), v(:), minv(:)
      integer, allocatable :: p(:), way(:)
      logical, allocatable :: used(:)
      real(dp) :: cur, delta

      n=size(cost,1)
      if (size(cost,2)/=n) error stop "assignment_lower_bound: square matrix required"
      allocate(u(0:n),v(0:n),p(0:n),way(0:n),minv(0:n),used(0:n))
      u=0.0_dp; v=0.0_dp; p=0; way=0
      do i=1,n
         p(0)=i; j0=0; minv=huge(1.0_dp); used=.false.
         do
            used(j0)=.true.; i0=p(j0); delta=huge(1.0_dp); j1=0
            do j=1,n
               if (used(j)) cycle
               cur=real(cost(i0,j),dp)-u(i0)-v(j)
               if (cur<minv(j)) then; minv(j)=cur; way(j)=j0; end if
               if (minv(j)<delta) then; delta=minv(j); j1=j; end if
            end do
            do j=0,n
               if (used(j)) then; u(p(j))=u(p(j))+delta; v(j)=v(j)-delta
               else if (j>0) then; minv(j)=minv(j)-delta; end if
            end do
            j0=j1
            if (p(j0)==0) exit
         end do
         do
            j1=way(j0); p(j0)=p(j1); j0=j1
            if (j0==0) exit
         end do
      end do
      obj=nint(-v(0),kind=i8)
   end function assignment_lower_bound

   function nearest_neighbor_upper_bound(cost, no_go) result(best)
      integer(i8), intent(in) :: cost(:,:), no_go
      integer(i8) :: best
      integer :: n, start, cur, next, step, j
      integer(i8) :: total, wbest
      logical, allocatable :: used(:)
      n=size(cost,1); allocate(used(n)); best=no_go
      ! Best nearest-neighbor tour over all starting cities. This replaces TSP::solve_TSP
      ! while retaining the role of a deterministic heuristic upper bound.
      do start=1,n
         used=.false.; used(start)=.true.; cur=start; total=0_i8
         do step=2,n
            next=0; wbest=no_go
            do j=1,n
               if (.not. used(j) .and. cost(cur,j)<wbest) then
                  next=j; wbest=cost(cur,j)
               end if
            end do
            if (next==0 .or. wbest>=no_go) then; total=no_go; exit; end if
            total=total+wbest; used(next)=.true.; cur=next
         end do
         if (total<no_go .and. cost(cur,start)<no_go) total=total+cost(cur,start)
         if (total<best) best=total
      end do
   end function nearest_neighbor_upper_bound

   function cheapest_insertion_upper_bound(cost, no_go) result(best)
      integer(i8), intent(in) :: cost(:,:), no_go
      integer(i8) :: best
      integer :: n, a, b, m, k, p, q, insk, insp, t
      integer, allocatable :: cyc(:), work(:)
      logical, allocatable :: used(:)
      integer(i8) :: total, delta, bestdelta

      n=size(cost,1); best=no_go
      allocate(cyc(n),work(n),used(n))
      ! Try every feasible two-city seed. This is deterministic and robust for sparse directed graphs.
      do a=1,n
         do b=1,n
            if (a==b) cycle
            if (cost(a,b)>=no_go .or. cost(b,a)>=no_go) cycle
            cyc=0; cyc(1)=a; cyc(2)=b; m=2
            used=.false.; used(a)=.true.; used(b)=.true.
            do while (m<n)
               bestdelta=no_go; insk=0; insp=0
               do k=1,n
                  if (used(k)) cycle
                  do p=1,m
                     q=merge(p+1,1,p<m)
                     if (cost(cyc(p),k)>=no_go .or. cost(k,cyc(q))>=no_go) cycle
                     if (cost(cyc(p),cyc(q))>=no_go) cycle
                     delta=cost(cyc(p),k)+cost(k,cyc(q))-cost(cyc(p),cyc(q))
                     if (delta<bestdelta) then
                        bestdelta=delta; insk=k; insp=p
                     end if
                  end do
               end do
               if (insk==0) exit
               work=cyc
               cyc(1:insp)=work(1:insp)
               cyc(insp+1)=insk
               if (insp<m) cyc(insp+2:m+1)=work(insp+1:m)
               m=m+1; used(insk)=.true.
            end do
            if (m/=n) cycle
            total=0_i8
            do t=1,n
               q=merge(t+1,1,t<n)
               if (cost(cyc(t),cyc(q))>=no_go) then
                  total=no_go; exit
               end if
               total=total+cost(cyc(t),cyc(q))
            end do
            best=min(best,total)
         end do
      end do
      if (best>=no_go) best=nearest_neighbor_upper_bound(cost,no_go)
   end function cheapest_insertion_upper_bound

end module nilde_tsp
