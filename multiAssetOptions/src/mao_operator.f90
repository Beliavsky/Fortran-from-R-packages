! This file is part of multiAssetOptions-fortran.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module mao_operator
   use mao_kinds, only: dp
   use mao_status, only: status_type, clear_status, set_status, &
      mao_invalid_argument, mao_allocation_error
   use mao_types, only: grid_set
   use mao_grid, only: decode_index
   use mao_sparse, only: csr_matrix
   implicit none
   private

   public :: build_fdm_operator

contains

   subroutine build_fdm_operator(grid, rf, q, vol, rho, operator, status)
      type(grid_set), intent(in) :: grid
      real(dp), intent(in) :: rf
      real(dp), intent(in) :: q(:), vol(:), rho(:,:)
      type(csr_matrix), intent(out) :: operator
      type(status_type), intent(out) :: status

      integer, allocatable :: indices(:), local_col(:)
      integer, allocatable :: all_col(:), row_ptr(:), tmp_col(:)
      real(dp), allocatable :: local_value(:), all_value(:), tmp_value(:)
      integer :: n, row, i, j, k, l, n1, n2, nd2, nlocal
      integer :: pos1(2), pos2(2), posd2(3)
      real(dp) :: coef1(2), coef2(2), coefd2(3)
      real(dp) :: s_i, drift, diffusion, cross_scale
      integer :: max_per_row, max_nnz, nnz, stat, target_col
      real(dp), parameter :: drop_tolerance = 100.0_dp * tiny(1.0_dp)

      call clear_status(status)
      n = size(grid%asset)
      if (n < 1 .or. size(q) /= n .or. size(vol) /= n .or. &
          size(rho,1) /= n .or. size(rho,2) /= n) then
         call set_status(status, mao_invalid_argument, &
            'invalid FDM operator dimensions')
         return
      end if
      do i = 1, n
         if (grid%dims(i) < 3) then
            call set_status(status, mao_invalid_argument, &
               'each asset grid requires at least three nodes')
            return
         end if
         do j = 1, grid%dims(i)-1
            if (grid%asset(i)%x(j+1) <= grid%asset(i)%x(j)) then
               call set_status(status, mao_invalid_argument, &
                  'asset grids must be strictly increasing')
               return
            end if
         end do
      end do

      max_per_row = 2*n*n + 1
      if (grid%n_nodes > huge(max_nnz) / max_per_row) then
         call set_status(status, mao_invalid_argument, &
            'FDM sparse matrix is too large for default integer indexing')
         return
      end if
      max_nnz = grid%n_nodes * max_per_row
      allocate(indices(n), local_col(max_per_row), local_value(max_per_row), &
         all_col(max_nnz), all_value(max_nnz), row_ptr(grid%n_nodes+1), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate FDM operator workspace')
         return
      end if

      nnz = 0
      row_ptr(1) = 1
      do row = 1, grid%n_nodes
         call decode_index(row,grid%dims,indices)
         nlocal = 0
         call add_local(row,-rf,local_col,local_value,nlocal)

         do i = 1, n
            s_i = grid%asset(i)%x(indices(i))
            drift = (rf-q(i))*s_i
            diffusion = 0.5_dp*vol(i)**2*s_i**2

            call first_stencil(grid%asset(i)%x,indices(i),pos1,coef1,n1)
            do k = 1, n1
               target_col = row + pos1(k)*grid%strides(i)
               call add_local(target_col,drift*coef1(k), &
                  local_col,local_value,nlocal)
            end do

            call second_stencil(grid%asset(i)%x,indices(i),posd2,coefd2,nd2)
            do k = 1, nd2
               target_col = row + posd2(k)*grid%strides(i)
               call add_local(target_col,diffusion*coefd2(k), &
                  local_col,local_value,nlocal)
            end do
         end do

         do i = 2, n
            call first_stencil(grid%asset(i)%x,indices(i),pos1,coef1,n1)
            do j = 1, i-1
               call first_stencil(grid%asset(j)%x,indices(j),pos2,coef2,n2)
               cross_scale = rho(i,j)*vol(i)*vol(j) * &
                  grid%asset(i)%x(indices(i))*grid%asset(j)%x(indices(j))
               do k = 1, n1
                  do l = 1, n2
                     target_col = row + pos1(k)*grid%strides(i) + &
                        pos2(l)*grid%strides(j)
                     call add_local(target_col,cross_scale*coef1(k)*coef2(l), &
                        local_col,local_value,nlocal)
                  end do
               end do
            end do
         end do

         call sort_local(local_col,local_value,nlocal)
         do k = 1, nlocal
            if (abs(local_value(k)) > drop_tolerance .or. local_col(k) == row) then
               nnz = nnz + 1
               all_col(nnz) = local_col(k)
               all_value(nnz) = local_value(k)
            end if
         end do
         row_ptr(row+1) = nnz + 1
      end do

      allocate(tmp_col(nnz), tmp_value(nnz), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to finalize FDM operator')
         return
      end if
      tmp_col = all_col(1:nnz)
      tmp_value = all_value(1:nnz)

      operator%nrow = grid%n_nodes
      operator%ncol = grid%n_nodes
      call move_alloc(row_ptr,operator%row_ptr)
      call move_alloc(tmp_col,operator%col_ind)
      call move_alloc(tmp_value,operator%value)
   end subroutine build_fdm_operator

   subroutine first_stencil(x, index_value, positions, coefficients, count)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: index_value
      integer, intent(out) :: positions(2), count
      real(dp), intent(out) :: coefficients(2)
      real(dp) :: hminus, hplus

      if (index_value == 1) then
         hplus = x(2)-x(1)
         count = 2
         positions = [0,1]
         coefficients = [-1.0_dp/hplus,1.0_dp/hplus]
      else if (index_value == size(x)) then
         hminus = x(index_value)-x(index_value-1)
         count = 2
         positions = [-1,0]
         coefficients = [-1.0_dp/hminus,1.0_dp/hminus]
      else
         hminus = x(index_value)-x(index_value-1)
         hplus = x(index_value+1)-x(index_value)
         count = 2
         positions = [-1,1]
         coefficients = [-1.0_dp/(hminus+hplus), &
            1.0_dp/(hminus+hplus)]
      end if
   end subroutine first_stencil

   subroutine second_stencil(x, index_value, positions, coefficients, count)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: index_value
      integer, intent(out) :: positions(3), count
      real(dp), intent(out) :: coefficients(3)
      real(dp) :: hminus, hplus

      positions = 0
      coefficients = 0.0_dp
      if (index_value == 1 .or. index_value == size(x)) then
         count = 0
      else
         hminus = x(index_value)-x(index_value-1)
         hplus = x(index_value+1)-x(index_value)
         count = 3
         positions = [-1,0,1]
         coefficients(1) = 2.0_dp/(hminus*(hminus+hplus))
         coefficients(2) = -2.0_dp/(hminus*hplus)
         coefficients(3) = 2.0_dp/(hplus*(hminus+hplus))
      end if
   end subroutine second_stencil

   subroutine add_local(column, value, columns, values, count)
      integer, intent(in) :: column
      real(dp), intent(in) :: value
      integer, intent(inout) :: columns(:), count
      real(dp), intent(inout) :: values(:)
      integer :: i

      do i = 1, count
         if (columns(i) == column) then
            values(i) = values(i) + value
            return
         end if
      end do
      count = count + 1
      columns(count) = column
      values(count) = value
   end subroutine add_local

   subroutine sort_local(columns, values, count)
      integer, intent(inout) :: columns(:)
      real(dp), intent(inout) :: values(:)
      integer, intent(in) :: count
      integer :: i, j, column_key
      real(dp) :: value_key

      do i = 2, count
         column_key = columns(i)
         value_key = values(i)
         j = i-1
         do while (j >= 1)
            if (columns(j) <= column_key) exit
            columns(j+1) = columns(j)
            values(j+1) = values(j)
            j = j-1
         end do
         columns(j+1) = column_key
         values(j+1) = value_key
      end do
   end subroutine sort_local

end module mao_operator
