! SPDX-License-Identifier: GPL-2.0-or-later
! Multiple-wavelet pre/postfilter kernels translated from wavethresh 4.7.3.
module wavethresh_multiwavelet
   use wavethresh_types, only : dp, mwd_t, wt_matrix_t, multiple_filter_t
   use wavethresh_filters, only : mfilter_select
   implicit none
   private

   public :: mprefilter, mpostfilter, mwd, mwr

contains

   function mwd(data, prefilter_type, filter_type, boundary) result(object)
      real(dp), intent(in) :: data(:) !! Scalar observations for the periodic discrete multiple-wavelet transform.
      character(len=*), intent(in), optional :: prefilter_type !! Prefilter name; default selects upstream interpolation.
      character(len=*), intent(in), optional :: filter_type !! Multiple-wavelet family; default is Geronimo.
      character(len=*), intent(in), optional :: boundary !! Boundary handling; only periodic is currently supported.
      type(mwd_t) :: object
      type(wt_matrix_t), allocatable :: temp_detail(:)
      type(wt_matrix_t), allocatable :: temp_scaling(:)
      type(multiple_filter_t) :: filter
      real(dp), allocatable :: detail(:,:)
      real(dp), allocatable :: smooth(:,:)
      real(dp), allocatable :: work(:,:)
      character(len=16) :: bc
      character(len=16) :: pf
      character(len=24) :: ft
      integer :: base_levels
      integer :: levels
      integer :: level
      integer :: starting_columns
      integer :: step

      pf = "default"
      if (present(prefilter_type)) pf = prefilter_type
      ft = "Geronimo"
      if (present(filter_type)) ft = filter_type
      bc = "periodic"
      if (present(boundary)) bc = boundary
      object%n_original = size(data)
      object%prefilter = pf
      object%boundary = bc
      if (trim(bc) /= "periodic") then
         object%message = "mwd currently supports periodic boundary handling only"
         return
      end if
      filter = mfilter_select(trim(ft))
      object%filter = filter
      if (.not. filter%ok) then
         object%message = filter%message
         return
      end if
      if (filter%nphi /= filter%npsi) then
         object%message = "mwd currently requires equal multiple scaling and wavelet ranks"
         return
      end if
      base_levels = exact_power_levels(size(data) / filter%nphi, filter%ndecim)
      if (mod(size(data), filter%nphi) /= 0 .or. base_levels < 1) then
         object%message = "data length is incompatible with the selected multiple wavelet"
         return
      end if
      levels = base_levels
      if (trim(pf) == "Repeat") levels = levels + 1
      if (trim(pf) == "Repeat" .and. trim(ft) /= "Geronimo") then
         object%message = "Repeat prefiltering is implemented only for Geronimo"
         return
      end if
      starting_columns = filter%ndecim**levels
      work = mprefilter(data, trim(pf), trim(ft), levels, starting_columns, filter%nphi, filter%npsi, filter%ndecim)
      if (size(work, 1) /= filter%nphi .or. size(work, 2) /= starting_columns) then
         object%message = "requested multiple-wavelet prefilter is not available"
         return
      end if

      object%nlevels = levels
      allocate(object%detail(0:levels - 1))
      allocate(object%scaling(0:levels))
      object%scaling(levels)%values = work
      allocate(temp_detail(levels))
      allocate(temp_scaling(levels))
      do step = 1, levels
         call multiple_analysis_step(work, filter, detail, smooth)
         temp_detail(step)%values = detail
         temp_scaling(step)%values = smooth
         work = smooth
      end do
      do step = 1, levels
         level = levels - step
         object%detail(level)%values = temp_detail(step)%values
         object%scaling(level)%values = temp_scaling(step)%values
      end do
      object%ok = .true.
      object%message = "ok"
   end function mwd

   function mwr(object, start_level) result(data)
      type(mwd_t), intent(in) :: object !! Periodic multiple-wavelet decomposition to reconstruct.
      integer, intent(in), optional :: start_level !! Coarsest R-style level used for synthesis; default is zero.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: next(:,:)
      real(dp), allocatable :: work(:,:)
      integer :: first
      integer :: level

      if (.not. object%ok .or. object%nlevels < 1 .or. trim(object%boundary) /= "periodic") then
         allocate(data(0))
         return
      end if
      first = 0
      if (present(start_level)) first = start_level
      if (first < 0 .or. first >= object%nlevels) then
         allocate(data(0))
         return
      end if
      if (.not. allocated(object%scaling(first)%values)) then
         allocate(data(0))
         return
      end if
      work = object%scaling(first)%values
      do level = first, object%nlevels - 1
         if (.not. allocated(object%detail(level)%values)) then
            allocate(data(0))
            return
         end if
         call multiple_synthesis_step(work, object%detail(level)%values, object%filter, next)
         work = next
      end do
      data = mpostfilter(work, trim(object%prefilter), trim(object%filter%filter_type), object%filter%nphi, &
         object%filter%npsi, object%filter%ndecim, object%nlevels)
   end function mwr

   function mprefilter(data, prefilter_type, filter_type, nlevels, nvecs_c, nphi, npsi, ndecim) result(c)
      real(dp), intent(in) :: data(:) !! Original scalar observations to convert to vector-valued scaling coefficients.
      character(len=*), intent(in) :: prefilter_type !! Upstream prefilter name, such as Identity, Interp, or Xia.
      character(len=*), intent(in) :: filter_type !! Multiple-wavelet family; Geronimo and Donovan3 are supported.
      integer, intent(in) :: nlevels !! Number of multiple-wavelet decomposition levels; must be positive.
      integer, intent(in) :: nvecs_c !! Number of scaling-vector columns reserved by the first/last database.
      integer, intent(in) :: nphi !! Number of scaling-function components for the selected multiple wavelet.
      integer, intent(in) :: npsi !! Number of wavelet components for the selected multiple wavelet.
      integer, intent(in) :: ndecim !! Decimation factor of the selected multiple wavelet.
      real(dp), allocatable :: c(:,:)
      real(dp) :: a
      real(dp) :: b
      real(dp) :: cc
      real(dp) :: d
      real(dp) :: epsilon1
      real(dp) :: epsilon2
      real(dp) :: r
      real(dp) :: root2
      real(dp) :: s
      real(dp) :: w
      real(dp) :: x
      real(dp) :: qb(2, 2)
      real(dp) :: qa(2, 2)
      real(dp) :: qz(2, 2)
      real(dp) :: partition(2, size(data) / 2)
      integer :: j
      integer :: left
      integer :: lc
      integer :: ncolumns
      integer :: right

      if (nlevels < 1 .or. nvecs_c < 1 .or. nphi < 1 .or. npsi < 1 .or. ndecim < 2) then
         allocate(c(0, 0))
         return
      end if
      if (size(data) < 1) then
         allocate(c(0, 0))
         return
      end if
      allocate(c(nphi, nvecs_c), source=0.0_dp)

      select case (trim(filter_type))
      case ("Geronimo")
         if (nphi /= 2 .or. npsi /= 2 .or. ndecim /= 2 .or. mod(size(data), 2) /= 0) then
            deallocate(c)
            allocate(c(0, 0))
            return
         end if
         select case (trim(prefilter_type))
         case ("Minimal")
            ncolumns = size(data) / 2
            if (nvecs_c < ncolumns) then
               deallocate(c)
               allocate(c(0, 0))
               return
            end if
            do j = 1, ncolumns
               c(1, j) = sqrt(2.0_dp) * data(2 * j - 1) - data(2 * j) / sqrt(2.0_dp)
               c(2, j) = 0.5_dp * data(2 * j - 1)
            end do
         case ("Identity")
            ncolumns = size(data) / 2
            if (nvecs_c < ncolumns) then
               deallocate(c)
               allocate(c(0, 0))
               return
            end if
            do j = 1, ncolumns
               c(1, j) = data(2 * j - 1)
               c(2, j) = data(2 * j)
            end do
         case ("Repeat")
            ncolumns = size(data)
            if (nvecs_c < ncolumns) then
               deallocate(c)
               allocate(c(0, 0))
               return
            end if
            c(1, 1:ncolumns) = sqrt(2.0_dp) * data
            c(2, 1:ncolumns) = data
         case ("Interp", "default")
            ncolumns = size(data) / 2
            if (nvecs_c < ncolumns) then
               deallocate(c)
               allocate(c(0, 0))
               return
            end if
            r = sqrt(25.0_dp / 96.0_dp)
            s = sqrt(1.0_dp / 3.0_dp)
            a = -0.3_dp
            do j = 1, ncolumns
               c(2, j) = s * data(2 * j)
            end do
            c(1, 1) = r * (data(1) - a * (data(size(data)) + data(2)))
            do j = 2, ncolumns
               c(1, j) = r * (data(2 * j - 1) - a * (data(2 * j - 2) + data(2 * j)))
            end do
         case ("Xia")
            ncolumns = size(data) / 2
            if (nvecs_c < ncolumns) then
               deallocate(c)
               allocate(c(0, 0))
               return
            end if
            epsilon1 = 0.0_dp
            epsilon2 = 0.1_dp
            root2 = sqrt(2.0_dp)
            x = (2.0_dp * root2) / (5.0_dp * (root2 * epsilon2 - epsilon1))
            a = (x - epsilon1 + 2.0_dp * root2 * epsilon2) / 2.0_dp
            b = (x + epsilon1 - 2.0_dp * root2 * epsilon2) / 2.0_dp
            cc = (x + 4.0_dp * epsilon1 - 3.0_dp * root2 * epsilon2) / (2.0_dp * root2)
            d = (x - 4.0_dp * epsilon1 + 3.0_dp * root2 * epsilon2) / (2.0_dp * root2)
            do j = 1, ncolumns
               c(1, j) = a * data(2 * j) + b * data(2 * j - 1)
               c(2, j) = cc * data(2 * j) + d * data(2 * j - 1)
            end do
         case ("Roach1")
            ncolumns = size(data) / 2
            if (nvecs_c < ncolumns) then
               deallocate(c)
               allocate(c(0, 0))
               return
            end if
            partition = reshape(data, shape(partition))
            qb(1, :) = [0.2318485184_dp, 0.3298205429_dp]
            qb(2, :) = [-0.1629787369_dp, -0.2318485184_dp]
            qa(1, :) = [-0.2945950581_dp, 0.8187567536_dp]
            qa(2, :) = [0.8187567536_dp, 0.2945950581_dp]
            qz(1, :) = [0.2318485184_dp, -0.1629787369_dp]
            qz(2, :) = [0.3298205429_dp, -0.2318485184_dp]
            do j = 1, ncolumns
               left = modulo(j - 2, ncolumns) + 1
               right = modulo(j, ncolumns) + 1
               c(:, j) = matmul(qb, partition(:, left)) + matmul(qa, partition(:, j)) + &
                  matmul(qz, partition(:, right))
            end do
         case ("Roach3")
            ncolumns = size(data) / 2
            if (nvecs_c < ncolumns) then
               deallocate(c)
               allocate(c(0, 0))
               return
            end if
            partition = reshape(data, shape(partition))
            qb(:, 1) = [-0.003600312909_dp, 0.0001535859223_dp]
            qb(:, 2) = [0.08439740344_dp, -0.003600312909_dp]
            qa(:, 1) = [0.9927991855_dp, -0.08485816121_dp]
            qa(:, 2) = [0.08485816121_dp, 0.9927991855_dp]
            qz(:, 1) = [-0.003600312909_dp, -0.08439740344_dp]
            qz(:, 2) = [-0.0001535859223_dp, -0.003600312909_dp]
            do j = 1, ncolumns
               left = modulo(j - 2, ncolumns) + 1
               right = modulo(j, ncolumns) + 1
               c(:, j) = matmul(qb, partition(:, left)) + matmul(qa, partition(:, j)) + &
                  matmul(qz, partition(:, right))
            end do
         case default
            deallocate(c)
            allocate(c(0, 0))
         end select

      case ("Donovan3")
         if (nphi /= 3 .or. npsi /= 3 .or. ndecim /= 2 .or. mod(size(data), 3) /= 0) then
            deallocate(c)
            allocate(c(0, 0))
            return
         end if
         lc = size(data) / 3
         if (nvecs_c < lc) then
            deallocate(c)
            allocate(c(0, 0))
            return
         end if
         select case (trim(prefilter_type))
         case ("Identity")
            do j = 1, lc
               c(1, j) = data(3 * j - 2)
               c(2, j) = data(3 * j - 1)
               c(3, j) = data(3 * j)
            end do
         case ("Linear")
            do j = 1, lc
               c(1, j) = -0.76885512_dp * data(3 * j - 2) + data(3 * j - 1)
               c(2, j) = -0.56536683_dp * data(3 * j - 2) + data(3 * j - 1)
               c(3, j) = 0.09767654_dp * data(3 * j - 2) - data(3 * j - 1) + data(3 * j)
            end do
         case ("Interp", "default")
            w = sqrt(5.0_dp)
            do j = 1, lc
               c(1, j) = data(3 * j - 2) * sqrt(7.0_dp / 11.0_dp)
            end do
            c(3, 1) = donovan_c3(data(2), data(3), c(1, lc), c(1, 1), w)
            do j = 2, lc
               c(3, j) = donovan_c3(data(3 * j - 1), data(3 * j), c(1, j - 1), c(1, j), w)
            end do
            c(2, 1) = donovan_c2(data(2), c(1, lc), c(1, 1), c(3, 1), w)
            do j = 2, lc
               c(2, j) = donovan_c2(data(3 * j - 1), c(1, j - 1), c(1, j), c(3, j), w)
            end do
         case default
            deallocate(c)
            allocate(c(0, 0))
         end select

      case default
         deallocate(c)
         allocate(c(0, 0))
      end select
   end function mprefilter

   function mpostfilter(c, prefilter_type, filter_type, nphi, npsi, ndecim, nlevels) result(data)
      real(dp), intent(in) :: c(:,:) !! Reconstructed vector-valued scaling coefficients at the original data resolution.
      character(len=*), intent(in) :: prefilter_type !! Upstream prefilter name whose matching postfilter is applied.
      character(len=*), intent(in) :: filter_type !! Multiple-wavelet family; Geronimo and Donovan3 are supported.
      integer, intent(in) :: nphi !! Number of scaling-function components for the selected multiple wavelet.
      integer, intent(in) :: npsi !! Number of wavelet components for the selected multiple wavelet.
      integer, intent(in) :: ndecim !! Decimation factor of the selected multiple wavelet.
      integer, intent(in) :: nlevels !! Number of multiple-wavelet decomposition levels used by reconstruction.
      real(dp), allocatable :: data(:)
      real(dp) :: a
      real(dp) :: b
      real(dp) :: cc
      real(dp) :: d
      real(dp) :: epsilon1
      real(dp) :: epsilon2
      real(dp) :: root2
      real(dp) :: t
      real(dp) :: u
      real(dp) :: w
      real(dp) :: x
      real(dp) :: qb(2, 2)
      real(dp) :: qa(2, 2)
      real(dp) :: qz(2, 2)
      real(dp), allocatable :: partition(:,:)
      integer :: j
      integer :: left
      integer :: lc
      integer :: ncolumns
      integer :: ndata
      integer :: right

      if (nlevels < 1 .or. nphi < 1 .or. npsi < 1 .or. ndecim < 2) then
         allocate(data(0))
         return
      end if
      ndata = ndecim**nlevels * nphi
      if (trim(prefilter_type) == "Repeat") ndata = ndecim**(nlevels - 1) * nphi
      if (ndata < 1) then
         allocate(data(0))
         return
      end if
      allocate(data(ndata), source=0.0_dp)

      select case (trim(filter_type))
      case ("Geronimo")
         if (nphi /= 2 .or. npsi /= 2 .or. ndecim /= 2) then
            deallocate(data)
            allocate(data(0))
            return
         end if
         select case (trim(prefilter_type))
         case ("Minimal")
            ncolumns = ndata / 2
            if (size(c, 1) < 2 .or. size(c, 2) < ncolumns) then
               deallocate(data)
               allocate(data(0))
               return
            end if
            do j = 1, ncolumns
               data(2 * j - 1) = 2.0_dp * c(2, j)
               data(2 * j) = -sqrt(2.0_dp) * c(1, j) + 4.0_dp * c(2, j)
            end do
         case ("Identity")
            ncolumns = ndata / 2
            if (size(c, 1) < 2 .or. size(c, 2) < ncolumns) then
               deallocate(data)
               allocate(data(0))
               return
            end if
            do j = 1, ncolumns
               data(2 * j - 1) = c(1, j)
               data(2 * j) = c(2, j)
            end do
         case ("Repeat")
            ncolumns = ndata
            if (size(c, 1) < 2 .or. size(c, 2) < ncolumns) then
               deallocate(data)
               allocate(data(0))
               return
            end if
            do j = 1, ncolumns
               data(j) = (c(2, j) + c(1, j) / sqrt(2.0_dp)) / 2.0_dp
            end do
         case ("Interp", "default")
            ncolumns = ndata / 2
            if (size(c, 1) < 2 .or. size(c, 2) < ncolumns) then
               deallocate(data)
               allocate(data(0))
               return
            end if
            t = sqrt(96.0_dp / 25.0_dp)
            u = sqrt(3.0_dp)
            do j = 1, ncolumns
               data(2 * j) = u * c(2, j)
            end do
            data(1) = t * c(1, 1) - 0.3_dp * (data(ndata) + data(2))
            do j = 2, ncolumns
               data(2 * j - 1) = t * c(1, j) - 0.3_dp * (data(2 * j - 2) + data(2 * j))
            end do
         case ("Xia")
            ncolumns = ndata / 2
            if (size(c, 1) < 2 .or. size(c, 2) < ncolumns) then
               deallocate(data)
               allocate(data(0))
               return
            end if
            epsilon1 = 0.0_dp
            epsilon2 = 0.1_dp
            root2 = sqrt(2.0_dp)
            x = (2.0_dp * root2) / (5.0_dp * (root2 * epsilon2 - epsilon1))
            a = (x - epsilon1 + 2.0_dp * root2 * epsilon2) / 2.0_dp
            b = (x + epsilon1 - 2.0_dp * root2 * epsilon2) / 2.0_dp
            cc = (x + 4.0_dp * epsilon1 - 3.0_dp * root2 * epsilon2) / (2.0_dp * root2)
            d = (x - 4.0_dp * epsilon1 + 3.0_dp * root2 * epsilon2) / (2.0_dp * root2)
            do j = 1, ncolumns
               data(2 * j) = d * c(1, j) - b * c(2, j)
               data(2 * j - 1) = a * c(2, j) - cc * c(1, j)
            end do
         case ("Roach1")
            ncolumns = ndata / 2
            if (size(c, 1) < 2 .or. size(c, 2) < ncolumns) then
               deallocate(data)
               allocate(data(0))
               return
            end if
            qb(1, :) = [0.2318485184_dp, 0.3298205429_dp]
            qb(2, :) = [-0.1629787369_dp, -0.2318485184_dp]
            qa(1, :) = [-0.2945950581_dp, 0.8187567536_dp]
            qa(2, :) = [0.8187567536_dp, 0.2945950581_dp]
            qz(1, :) = [0.2318485184_dp, -0.1629787369_dp]
            qz(2, :) = [0.3298205429_dp, -0.2318485184_dp]
            allocate(partition(2, ncolumns), source=0.0_dp)
            do j = 1, ncolumns
               left = modulo(j - 2, ncolumns) + 1
               right = modulo(j, ncolumns) + 1
               partition(:, j) = matmul(qb, c(:, left)) + matmul(qa, c(:, j)) + matmul(qz, c(:, right))
            end do
            data = reshape(partition, [ndata])
         case ("Roach3")
            ncolumns = ndata / 2
            if (size(c, 1) < 2 .or. size(c, 2) < ncolumns) then
               deallocate(data)
               allocate(data(0))
               return
            end if
            qz(1, :) = [-0.003600312909_dp, 0.0001535859223_dp]
            qz(2, :) = [0.08439740344_dp, -0.003600312909_dp]
            qa(1, :) = [0.9927991855_dp, -0.08485816121_dp]
            qa(2, :) = [0.08485816121_dp, 0.9927991855_dp]
            qb(1, :) = [-0.003600312909_dp, -0.08439740344_dp]
            qb(2, :) = [-0.0001535859223_dp, -0.003600312909_dp]
            allocate(partition(2, ncolumns), source=0.0_dp)
            do j = 1, ncolumns
               left = modulo(j - 2, ncolumns) + 1
               right = modulo(j, ncolumns) + 1
               partition(:, j) = matmul(qb, c(:, left)) + matmul(qa, c(:, j)) + matmul(qz, c(:, right))
            end do
            data = reshape(partition, [ndata])
         case default
            deallocate(data)
            allocate(data(0))
         end select

      case ("Donovan3")
         if (nphi /= 3 .or. npsi /= 3 .or. ndecim /= 2 .or. mod(ndata, 3) /= 0) then
            deallocate(data)
            allocate(data(0))
            return
         end if
         lc = ndata / 3
         if (size(c, 1) < 3 .or. size(c, 2) < lc) then
            deallocate(data)
            allocate(data(0))
            return
         end if
         select case (trim(prefilter_type))
         case ("Identity")
            do j = 1, lc
               data(3 * j - 2) = c(1, j)
               data(3 * j - 1) = c(2, j)
               data(3 * j) = c(3, j)
            end do
         case ("Linear")
            do j = 1, lc
               data(3 * j - 2) = -4.914288_dp * c(1, j) + 4.914288_dp * c(2, j)
               data(3 * j - 1) = -2.778375_dp * c(1, j) + 3.778375_dp * c(2, j)
               data(3 * j) = -2.298365_dp * c(1, j) + 3.298365_dp * c(2, j) + c(3, j)
            end do
         case ("Interp", "default")
            w = sqrt(5.0_dp)
            do j = 1, lc
               data(3 * j - 2) = c(1, j) * sqrt(11.0_dp / 7.0_dp)
            end do
            data(2) = donovan_x2(c(1, lc), c(1, 1), c(2, 1), c(3, 1), w)
            do j = 2, lc
               data(3 * j - 1) = donovan_x2(c(1, j - 1), c(1, j), c(2, j), c(3, j), w)
            end do
            data(3) = donovan_x3(c(1, lc), c(1, 1), c(2, 1), c(3, 1), w)
            do j = 2, lc
               data(3 * j) = donovan_x3(c(1, j - 1), c(1, j), c(2, j), c(3, j), w)
            end do
         case default
            deallocate(data)
            allocate(data(0))
         end select

      case default
         deallocate(data)
         allocate(data(0))
      end select
   end function mpostfilter

   subroutine multiple_analysis_step(c_in, filter, d_out, c_out)
      real(dp), intent(in) :: c_in(:,:) !! Fine-scale multiple scaling coefficients, components by translation.
      type(multiple_filter_t), intent(in) :: filter !! Multiple-wavelet analysis filters and dimensions.
      real(dp), allocatable, intent(out) :: d_out(:,:) !! Coarse detail coefficient matrix produced by this step.
      real(dp), allocatable, intent(out) :: c_out(:,:) !! Coarse scaling coefficient matrix produced by this step.
      integer :: input_index
      integer :: input_component
      integer :: output_component
      integer :: output_index
      integer :: output_count
      integer :: tap
      integer :: h_index
      integer :: g_index

      if (.not. filter%ok .or. size(c_in, 1) /= filter%nphi .or. &
         mod(size(c_in, 2), filter%ndecim) /= 0) then
         allocate(d_out(0, 0), c_out(0, 0))
         return
      end if
      output_count = size(c_in, 2) / filter%ndecim
      allocate(c_out(filter%nphi, output_count), source=0.0_dp)
      allocate(d_out(filter%npsi, output_count), source=0.0_dp)
      do output_index = 0, output_count - 1
         do tap = 0, filter%nh - 1
            input_index = modulo(filter%ndecim * output_index + tap, size(c_in, 2))
            do output_component = 0, filter%nphi - 1
               do input_component = 0, filter%nphi - 1
                  h_index = ((tap * filter%nphi + output_component) * filter%nphi + input_component) + 1
                  c_out(output_component + 1, output_index + 1) = c_out(output_component + 1, output_index + 1) + &
                     filter%h(h_index) * c_in(input_component + 1, input_index + 1)
               end do
            end do
            do output_component = 0, filter%npsi - 1
               do input_component = 0, filter%nphi - 1
                  g_index = ((tap * filter%npsi + output_component) * filter%nphi + input_component) + 1
                  d_out(output_component + 1, output_index + 1) = d_out(output_component + 1, output_index + 1) + &
                     filter%g(g_index) * c_in(input_component + 1, input_index + 1)
               end do
            end do
         end do
      end do
   end subroutine multiple_analysis_step

   subroutine multiple_synthesis_step(c_in, d_in, filter, c_out)
      real(dp), intent(in) :: c_in(:,:) !! Coarse multiple scaling coefficients, components by translation.
      real(dp), intent(in) :: d_in(:,:) !! Coarse multiple detail coefficients paired with c_in.
      type(multiple_filter_t), intent(in) :: filter !! Multiple-wavelet synthesis filters and dimensions.
      real(dp), allocatable, intent(out) :: c_out(:,:) !! Reconstructed next-finer scaling coefficient matrix.
      integer :: coarse_index
      integer :: detail_component
      integer :: fine_component
      integer :: fine_index
      integer :: coarse_count
      integer :: h_index
      integer :: g_index
      integer :: scaling_component
      integer :: tap

      if (.not. filter%ok .or. size(c_in, 1) /= filter%nphi .or. size(d_in, 1) /= filter%npsi .or. &
         size(c_in, 2) /= size(d_in, 2) .or. filter%nphi /= filter%npsi) then
         allocate(c_out(0, 0))
         return
      end if
      coarse_count = size(c_in, 2)
      allocate(c_out(filter%nphi, coarse_count * filter%ndecim), source=0.0_dp)
      do coarse_index = 0, coarse_count - 1
         do tap = 0, filter%nh - 1
            fine_index = modulo(filter%ndecim * coarse_index + tap, size(c_out, 2))
            do scaling_component = 0, filter%nphi - 1
               do fine_component = 0, filter%nphi - 1
                  h_index = ((tap * filter%nphi + scaling_component) * filter%nphi + fine_component) + 1
                  c_out(fine_component + 1, fine_index + 1) = c_out(fine_component + 1, fine_index + 1) + &
                     filter%h(h_index) * c_in(scaling_component + 1, coarse_index + 1)
               end do
            end do
            do detail_component = 0, filter%npsi - 1
               do fine_component = 0, filter%nphi - 1
                  g_index = ((tap * filter%npsi + detail_component) * filter%nphi + fine_component) + 1
                  c_out(fine_component + 1, fine_index + 1) = c_out(fine_component + 1, fine_index + 1) + &
                     filter%g(g_index) * d_in(detail_component + 1, coarse_index + 1)
               end do
            end do
         end do
      end do
   end subroutine multiple_synthesis_step

   pure function exact_power_levels(value, base) result(levels)
      integer, intent(in) :: value !! Positive candidate exact power of base.
      integer, intent(in) :: base !! Integer radix greater than one.
      integer :: levels
      integer :: work

      if (value < 1 .or. base < 2) then
         levels = -1
         return
      end if
      levels = 0
      work = value
      do while (work > 1 .and. mod(work, base) == 0)
         work = work / base
         levels = levels + 1
      end do
      if (work /= 1) levels = -1
   end function exact_power_levels

   pure function donovan_c3(x2, x3, c1_left, c1_here, w) result(value)
      real(dp), intent(in) :: x2 !! Second scalar sample in the current Donovan3 interpolation block.
      real(dp), intent(in) :: x3 !! Third scalar sample in the current Donovan3 interpolation block.
      real(dp), intent(in) :: c1_left !! First scaling component from the preceding periodic block.
      real(dp), intent(in) :: c1_here !! First scaling component from the current block.
      real(dp), intent(in) :: w !! Square root of five used by the Donovan3 interpolation formula.
      real(dp) :: value

      value = (sqrt(3.0_dp) * (x2 - x3) + c1_left * (-1.0_dp + 8.0_dp * w) / (3.0_dp * sqrt(231.0_dp)) + &
         c1_here * (1.0_dp + 8.0_dp * w) / (3.0_dp * sqrt(231.0_dp))) * &
         3.0_dp * sqrt(33.0_dp) * (16.0_dp - 5.0_dp * w) / (-203.0_dp + 88.0_dp * w)
   end function donovan_c3

   pure function donovan_c2(x2, c1_left, c1_here, c3_here, w) result(value)
      real(dp), intent(in) :: x2 !! Second scalar sample in the current Donovan3 interpolation block.
      real(dp), intent(in) :: c1_left !! First scaling component from the preceding periodic block.
      real(dp), intent(in) :: c1_here !! First scaling component from the current block.
      real(dp), intent(in) :: c3_here !! Third scaling component already computed for the current block.
      real(dp), intent(in) :: w !! Square root of five used by the Donovan3 interpolation formula.
      real(dp) :: value

      value = (sqrt(3.0_dp) * x2 + c1_left * (2.0_dp + 6.0_dp * w) / (3.0_dp * sqrt(231.0_dp)) + &
         c1_here * (3.0_dp + 2.0_dp * w) / (3.0_dp * sqrt(231.0_dp)) - &
         c3_here * (103.0_dp - 24.0_dp * w) / (3.0_dp * sqrt(33.0_dp) * (16.0_dp - 5.0_dp * w))) * &
         sqrt(3.0_dp) / 2.0_dp
   end function donovan_c2

   pure function donovan_x2(c1_left, c1_here, c2_here, c3_here, w) result(value)
      real(dp), intent(in) :: c1_left !! First scaling component from the preceding periodic Donovan3 block.
      real(dp), intent(in) :: c1_here !! First scaling component from the current Donovan3 block.
      real(dp), intent(in) :: c2_here !! Second scaling component from the current Donovan3 block.
      real(dp), intent(in) :: c3_here !! Third scaling component from the current Donovan3 block.
      real(dp), intent(in) :: w !! Square root of five used by the Donovan3 interpolation formula.
      real(dp) :: value

      value = (-(2.0_dp + 6.0_dp * w) * c1_left - (3.0_dp + 2.0_dp * w) * c1_here + &
         6.0_dp * sqrt(77.0_dp) * c2_here + (103.0_dp - 24.0_dp * w) * sqrt(7.0_dp) * c3_here / &
         (16.0_dp - 5.0_dp * w)) / (9.0_dp * sqrt(77.0_dp))
   end function donovan_x2

   pure function donovan_x3(c1_left, c1_here, c2_here, c3_here, w) result(value)
      real(dp), intent(in) :: c1_left !! First scaling component from the preceding periodic Donovan3 block.
      real(dp), intent(in) :: c1_here !! First scaling component from the current Donovan3 block.
      real(dp), intent(in) :: c2_here !! Second scaling component from the current Donovan3 block.
      real(dp), intent(in) :: c3_here !! Third scaling component from the current Donovan3 block.
      real(dp), intent(in) :: w !! Square root of five used by the Donovan3 interpolation formula.
      real(dp) :: value

      value = ((-3.0_dp + 2.0_dp * w) * c1_left / (3.0_dp * sqrt(231.0_dp)) + &
         (-2.0_dp + 6.0_dp * w) * c1_here / (3.0_dp * sqrt(231.0_dp)) + &
         2.0_dp * c2_here / sqrt(3.0_dp) + (306.0_dp - 112.0_dp * w) * c3_here / &
         ((16.0_dp - 5.0_dp * w) * 3.0_dp * sqrt(33.0_dp))) / sqrt(3.0_dp)
   end function donovan_x3

end module wavethresh_multiwavelet
