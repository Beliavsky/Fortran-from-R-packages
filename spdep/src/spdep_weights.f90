! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_weights
   use spdep_kinds, only : dp
   use spdep_types, only : int_vector, real_vector, neighbor_list, spatial_weights, weights_constants
   use spdep_graph, only : nb_adjacency_matrix
   implicit none
   private

   public :: card
   public :: nb2listw
   public :: nb2mat
   public :: mat2listw
   public :: listw2mat
   public :: lag_listw
   public :: lag_listw_matrix
   public :: spweights_constants
   public :: szero
   public :: listw2U

contains

   pure function card(nb) result(counts)
      type(neighbor_list), intent(in) :: nb !! Neighbor list whose per-region cardinalities are returned.
      integer, allocatable :: counts(:)
      integer :: i

      allocate(counts(nb%size()))
      do i = 1, nb%size()
         if (allocated(nb%neighbors(i)%values)) then
            counts(i) = size(nb%neighbors(i)%values)
         else
            counts(i) = 0
         end if
      end do
   end function card

   pure function nb2listw(nb, style, glist, zero_policy) result(listw)
      type(neighbor_list), intent(in) :: nb !! Neighbor topology defining which region pairs receive weights.
      character(len=*), intent(in), optional :: style !! Coding style B, W, C, U, S, or minmax; default is W.
      type(real_vector), intent(in), optional :: glist(:) !! General link weights conformable with nb; default is one per link.
      logical, intent(in), optional :: zero_policy !! Whether zero-neighbor rows are allowed; default is true in this Fortran API.
      type(spatial_weights) :: listw
      character(len=8) :: use_style
      logical :: allow_zero
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      real(dp) :: d
      real(dp) :: q
      real(dp) :: row_sum
      real(dp) :: scale
      real(dp) :: max_row
      real(dp) :: max_col
      real(dp), allocatable :: col_sum(:)
      type(real_vector), allocatable :: base(:)

      n = nb%size()
      use_style = "W"
      if (present(style)) use_style = adjustl(style)
      allow_zero = .true.
      if (present(zero_policy)) allow_zero = zero_policy
      allocate(base(n), listw%weights(n))
      listw%nb = nb
      listw%style = trim(use_style)
      listw%zero_policy = allow_zero
      do i = 1, n
         allocate(base(i)%values(size(nb%neighbors(i)%values)))
         if (present(glist)) then
            if (size(glist) /= n) then
               base(i)%values = 0.0_dp
            else if (size(glist(i)%values) /= size(base(i)%values)) then
               base(i)%values = 0.0_dp
            else
               base(i)%values = glist(i)%values
            end if
         else
            base(i)%values = 1.0_dp
         end if
      end do
      select case (trim(use_style))
      case ("B", "b")
         do i = 1, n
            allocate(listw%weights(i)%values(size(base(i)%values)))
            listw%weights(i)%values = base(i)%values
         end do
      case ("W", "w")
         do i = 1, n
            allocate(listw%weights(i)%values(size(base(i)%values)))
            row_sum = sum(base(i)%values)
            if (row_sum == 0.0_dp) then
               listw%weights(i)%values = 0.0_dp
            else
               listw%weights(i)%values = base(i)%values / row_sum
            end if
         end do
      case ("C", "c", "U", "u")
         d = 0.0_dp
         do i = 1, n
            d = d + sum(base(i)%values)
         end do
         scale = 0.0_dp
         if (d /= 0.0_dp) scale = real(n, dp) / d
         if (trim(use_style) == "U" .or. trim(use_style) == "u") then
            scale = scale / real(max(1, n), dp)
         end if
         do i = 1, n
            allocate(listw%weights(i)%values(size(base(i)%values)))
            listw%weights(i)%values = scale * base(i)%values
         end do
      case ("S", "s")
         q = 0.0_dp
         do i = 1, n
            row_sum = sqrt(sum(base(i)%values ** 2))
            allocate(listw%weights(i)%values(size(base(i)%values)))
            if (row_sum == 0.0_dp) then
               listw%weights(i)%values = 0.0_dp
            else
               listw%weights(i)%values = base(i)%values / row_sum
            end if
            q = q + sum(listw%weights(i)%values)
         end do
         scale = 0.0_dp
         if (q /= 0.0_dp) scale = real(n, dp) / q
         do i = 1, n
            listw%weights(i)%values = scale * listw%weights(i)%values
         end do
      case ("minmax", "MINMAX", "Minmax")
         allocate(col_sum(n))
         col_sum = 0.0_dp
         max_row = 0.0_dp
         do i = 1, n
            row_sum = sum(base(i)%values)
            max_row = max(max_row, row_sum)
            do k = 1, size(nb%neighbors(i)%values)
               j = nb%neighbors(i)%values(k)
               if (j >= 1 .and. j <= n) col_sum(j) = col_sum(j) + base(i)%values(k)
            end do
         end do
         max_col = 0.0_dp
         if (n > 0) max_col = maxval(col_sum)
         scale = min(max_row, max_col)
         do i = 1, n
            allocate(listw%weights(i)%values(size(base(i)%values)))
            if (scale == 0.0_dp) then
               listw%weights(i)%values = 0.0_dp
            else
               listw%weights(i)%values = base(i)%values / scale
            end if
         end do
      case default
         do i = 1, n
            allocate(listw%weights(i)%values(size(base(i)%values)))
            listw%weights(i)%values = base(i)%values
         end do
      end select
   end function nb2listw

   pure function listw2mat(listw) result(w)
      type(spatial_weights), intent(in) :: listw !! Spatial-weights object converted to a dense square weight matrix.
      real(dp), allocatable :: w(:, :)
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = listw%size()
      allocate(w(n, n))
      w = 0.0_dp
      do i = 1, n
         do k = 1, size(listw%nb%neighbors(i)%values)
            j = listw%nb%neighbors(i)%values(k)
            if (j >= 1 .and. j <= n) w(i, j) = listw%weights(i)%values(k)
         end do
      end do
   end function listw2mat

   pure function nb2mat(nb, style, glist, zero_policy) result(w)
      type(neighbor_list), intent(in) :: nb !! Neighbor topology converted to a dense numeric matrix.
      character(len=*), intent(in), optional :: style !! Weight coding style B, W, C, U, S, or minmax; default is W.
      type(real_vector), intent(in), optional :: glist(:) !! Optional link weights conformable with nb.
      logical, intent(in), optional :: zero_policy !! Whether zero-neighbor rows are permitted; default is true.
      real(dp), allocatable :: w(:, :)
      type(spatial_weights) :: listw

      listw = nb2listw(nb, style, glist, zero_policy)
      w = listw2mat(listw)
   end function nb2mat

   pure function mat2listw(w, style, zero_policy) result(listw)
      real(dp), intent(in) :: w(:, :) !! Square matrix; each nonzero element defines a directed link and general weight.
      character(len=*), intent(in), optional :: style !! Optional recoding style; default B retains matrix entries exactly.
      logical, intent(in), optional :: zero_policy !! Whether zero-neighbor rows are permitted; default is true.
      type(spatial_weights) :: listw
      type(neighbor_list) :: nb
      type(real_vector), allocatable :: glist(:)
      character(len=8) :: use_style
      integer :: n
      integer :: i
      integer :: j
      integer :: k
      integer :: m

      n = size(w, 1)
      if (size(w, 2) /= n) then
         allocate(nb%neighbors(0))
         listw = nb2listw(nb, "B")
         return
      end if
      allocate(nb%neighbors(n), glist(n))
      nb%self_included = .false.
      do i = 1, n
         m = count(w(i, :) /= 0.0_dp)
         allocate(nb%neighbors(i)%values(m), glist(i)%values(m))
         k = 0
         do j = 1, n
            if (w(i, j) /= 0.0_dp) then
               k = k + 1
               nb%neighbors(i)%values(k) = j
               glist(i)%values(k) = w(i, j)
               if (i == j) nb%self_included = .true.
            end if
         end do
      end do
      use_style = "B"
      if (present(style)) use_style = adjustl(style)
      listw = nb2listw(nb, trim(use_style), glist, zero_policy)
   end function mat2listw

   pure function lag_listw(listw, x) result(lagged)
      type(spatial_weights), intent(in) :: listw !! Spatial weights used to form each region's weighted neighbor sum.
      real(dp), intent(in) :: x(:) !! Numeric vector with one value per region.
      real(dp), allocatable :: lagged(:)
      integer :: n
      integer :: i
      integer :: j
      integer :: k

      n = listw%size()
      allocate(lagged(n))
      lagged = 0.0_dp
      if (size(x) /= n) return
      do i = 1, n
         do k = 1, size(listw%nb%neighbors(i)%values)
            j = listw%nb%neighbors(i)%values(k)
            if (j >= 1 .and. j <= n) then
               lagged(i) = lagged(i) + listw%weights(i)%values(k) * x(j)
            end if
         end do
      end do
   end function lag_listw

   pure function lag_listw_matrix(listw, x) result(lagged)
      type(spatial_weights), intent(in) :: listw !! Spatial weights applied independently to every matrix column.
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with one region per row and arbitrary response columns.
      real(dp), allocatable :: lagged(:, :)
      integer :: n
      integer :: j

      n = listw%size()
      allocate(lagged(n, size(x, 2)))
      lagged = 0.0_dp
      if (size(x, 1) /= n) return
      do j = 1, size(x, 2)
         lagged(:, j) = lag_listw(listw, x(:, j))
      end do
   end function lag_listw_matrix

   pure function spweights_constants(listw) result(constants)
      type(spatial_weights), intent(in) :: listw !! Spatial weights whose S0, S1, and S2 constants are required.
      type(weights_constants) :: constants
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: rs(:)
      real(dp), allocatable :: cs(:)
      integer :: i
      integer :: j
      integer :: n

      n = listw%size()
      constants%n = n
      w = listw2mat(listw)
      constants%s0 = sum(w)
      constants%s1 = 0.0_dp
      do i = 1, n
         do j = 1, n
            constants%s1 = constants%s1 + w(i, j) ** 2 + w(i, j) * w(j, i)
         end do
      end do
      allocate(rs(n), cs(n))
      rs = sum(w, dim = 2)
      cs = sum(w, dim = 1)
      constants%s2 = sum((rs + cs) ** 2)
   end function spweights_constants

   pure real(dp) function sZero(listw) result(s0)
      type(spatial_weights), intent(in) :: listw !! Spatial weights whose total weight S0 is required.
      integer :: i

      s0 = 0.0_dp
      do i = 1, listw%size()
         s0 = s0 + sum(listw%weights(i)%values)
      end do
   end function sZero

   pure function listw2U(listw) result(out)
      type(spatial_weights), intent(in) :: listw !! Spatial weights symmetrized as (W + transpose(W))/2.
      type(spatial_weights) :: out
      real(dp), allocatable :: w(:, :)

      w = listw2mat(listw)
      w = 0.5_dp * (w + transpose(w))
      out = mat2listw(w, "B", listw%zero_policy)
      out%style = "U"
   end function listw2U

end module spdep_weights
