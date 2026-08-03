! SPDX-License-Identifier: GPL-3.0-only
module portvine_ordering
   use portvine_kinds, only : dp
   use rugarch, only : normal_quantile, invert_matrix
   implicit none
   private
   public :: greedy_dvine_order, normal_scores, partial_correlation_value

contains

   subroutine normal_scores(u, z)
      real(dp), intent(in) :: u(:,:)
      real(dp), intent(out) :: z(size(u,1),size(u,2))
      real(dp), allocatable :: values(:)
      integer, allocatable :: index(:)
      integer :: i, j, k, n

      n = size(u,2)
      allocate(values(n),index(n))
      do i = 1, size(u,1)
         values = u(i,:)
         index = [(j,j=1,n)]
         call stable_index_sort(values,index)
         do k = 1, n
            z(i,index(k)) = normal_quantile((real(k,dp)-0.5_dp)/real(n,dp))
         end do
      end do
   end subroutine normal_scores

   subroutine stable_index_sort(values,index)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: index(:)
      integer :: i, j, key
      do i = 2, size(index)
         key = index(i)
         j = i-1
         do while (j >= 1)
            if (values(index(j)) <= values(key)) exit
            index(j+1) = index(j)
            j = j-1
         end do
         index(j+1) = key
      end do
   end subroutine stable_index_sort

   real(dp) function partial_correlation_value(x, first, second, conditioning) result(value)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: first, second
      integer, intent(in), optional :: conditioning(:)
      real(dp), allocatable :: r(:,:), rinv(:,:), y(:,:)
      real(dp) :: mi, si, sj
      integer :: p, n, i, j, info
      integer, allocatable :: ind(:)

      if (present(conditioning)) then
         p = 2+size(conditioning)
         allocate(ind(p))
         ind(1:2) = [first,second]
         if (p > 2) ind(3:) = conditioning
      else
         p = 2
         allocate(ind(2))
         ind = [first,second]
      end if
      n = size(x,2)
      allocate(y(p,n),r(p,p),rinv(p,p))
      do i = 1, p
         y(i,:) = x(ind(i),:)
      end do
      do i = 1, p
         mi = sum(y(i,:))/real(n,dp)
         y(i,:) = y(i,:)-mi
      end do
      do i = 1, p
         si = sqrt(max(tiny(1.0_dp),sum(y(i,:)**2)))
         do j = 1, p
            sj = sqrt(max(tiny(1.0_dp),sum(y(j,:)**2)))
            r(i,j) = sum(y(i,:)*y(j,:))/(si*sj)
         end do
      end do
      call invert_matrix(r,rinv,info)
      if (info /= 0) then
         do i = 1, p
            r(i,i) = r(i,i)+1.0e-8_dp
         end do
         call invert_matrix(r,rinv,info)
      end if
      if (info /= 0) then
         value = r(1,2)
      else
         value = -rinv(1,2)/sqrt(max(tiny(1.0_dp),rinv(1,1)*rinv(2,2)))
      end if
      value = max(-1.0_dp,min(1.0_dp,value))
   end function partial_correlation_value

   subroutine greedy_dvine_order(copula_data, order, cond_indices, cutoff_depth, status)
      real(dp), intent(in) :: copula_data(:,:)
      integer, intent(out) :: order(size(copula_data,1))
      integer, intent(in), optional :: cond_indices(:), cutoff_depth
      integer, intent(out), optional :: status
      real(dp), allocatable :: z(:,:)
      integer, allocatable :: free(:), suffix(:)
      logical, allocatable :: used(:)
      integer :: d, i, j, n_order, best, depth, max_depth, ncond
      real(dp) :: score, best_score, corr

      d = size(copula_data,1)
      order = 0
      if (present(status)) status = 0
      if (d < 2 .or. size(copula_data,2) < 3) then
         if (present(status)) status = 1
         return
      end if
      ncond = 0
      if (present(cond_indices)) ncond = size(cond_indices)
      if (ncond > 2) then
         if (present(status)) status = 1
         return
      end if
      allocate(z(d,size(copula_data,2)),used(d),free(d))
      call normal_scores(copula_data,z)
      used = .false.
      n_order = 0
      if (ncond > 0) then
         do i = 1, ncond
            if (cond_indices(i) < 1 .or. cond_indices(i) > d .or. used(cond_indices(i))) then
               if (present(status)) status = 1
               return
            end if
            n_order = n_order+1
            order(n_order) = cond_indices(i)
            used(cond_indices(i)) = .true.
         end do
      else
         best_score = -1.0_dp
         best = 1
         do i = 1, d-1
            do j = i+1, d
               corr = abs(partial_correlation_value(z,i,j))
               if (corr > best_score) then
                  best_score = corr
                  best = i
               end if
            end do
         end do
         n_order = 1
         order(1) = best
         used(best) = .true.
      end if

      do while (n_order < d)
         best = 0
         best_score = -huge(1.0_dp)
         do i = 1, d
            if (used(i)) cycle
            score = 0.0_dp
            max_depth = n_order
            if (present(cutoff_depth)) max_depth = min(max_depth,max(1,cutoff_depth))
            do depth = 1, max_depth
               j = n_order-depth+1
               if (depth == 1) then
                  corr = partial_correlation_value(z,i,order(j))
               else
                  allocate(suffix(depth-1))
                  suffix = order(j+1:n_order)
                  corr = partial_correlation_value(z,i,order(j),suffix)
                  deallocate(suffix)
               end if
               score = score+abs(corr)
            end do
            if (score > best_score) then
               best_score = score
               best = i
            end if
         end do
         if (best == 0) exit
         n_order = n_order+1
         order(n_order) = best
         used(best) = .true.
      end do
   end subroutine greedy_dvine_order

end module portvine_ordering
