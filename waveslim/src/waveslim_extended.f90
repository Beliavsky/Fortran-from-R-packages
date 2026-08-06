! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_extended
  use waveslim_kinds, only : dp, i8
  use waveslim_status, only : status_type, clear_status, set_status, &
    waveslim_invalid_input, waveslim_invalid_level
  use waveslim_types, only : wavelet_transform, wavelet_transform_2d, &
    packet_transform
  use waveslim_filters, only : wave_filter
  use waveslim_transform_1d, only : dwt, modwt, brick_wall, phase_shift
  use waveslim_transform_nd, only : dwt2_step, idwt2_step
  use waveslim_packet, only : dwpt, idwpt
  use waveslim_math, only : seed_rng
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan
  implicit none
  private

  type, public :: real_matrix
    real(dp), allocatable :: values(:,:)
  end type real_matrix

  type, public :: packet_level_2d
    type(real_matrix), allocatable :: node(:)
  end type packet_level_2d

  type, public :: packet_transform_2d
    type(packet_level_2d), allocatable :: level(:)
    character(len=16) :: wavelet = ''
    character(len=16) :: boundary = 'periodic'
    integer :: original_shape(2) = 0
    type(status_type) :: status
  contains
    procedure :: levels => packet_2d_levels
  end type packet_transform_2d

  type, public :: variance_change_result
    integer, allocatable :: level(:)
    real(dp), allocatable :: statistic(:)
    integer, allocatable :: location_dwt(:)
    integer, allocatable :: location_modwt(:)
    type(status_type) :: status
  contains
    procedure :: count => change_count
  end type variance_change_result

  public :: dwpt_2d, idwpt_2d, dwpt_boot
  public :: brick_wall_2d, phase_shift_packet
  public :: testing_hov, mult_loc

contains

  integer function packet_2d_levels(self)
    class(packet_transform_2d), intent(in) :: self
    if (allocated(self%level)) then
      packet_2d_levels = size(self%level)-1
    else
      packet_2d_levels = 0
    end if
  end function packet_2d_levels

  integer function change_count(self)
    class(variance_change_result), intent(in) :: self
    if (allocated(self%level)) then
      change_count = size(self%level)
    else
      change_count = 0
    end if
  end function change_count

  function dwpt_2d(x, wf, n_levels, boundary) result(tree)
    real(dp), intent(in) :: x(:,:)
    character(len=*), intent(in), optional :: wf, boundary
    integer, intent(in), optional :: n_levels
    type(packet_transform_2d) :: tree
    character(len=16) :: wname, bname
    integer :: levels, j, parent, side, ix, iy
    integer :: xlow, xhigh, ylow, yhigh
    real(dp), allocatable :: ll(:,:), lh(:,:), hl(:,:), hh(:,:)
    type(status_type) :: status

    wname = 'la8'
    if (present(wf)) wname = trim(wf)
    bname = 'periodic'
    if (present(boundary)) bname = trim(boundary)
    levels = 4
    if (present(n_levels)) levels = n_levels
    tree%wavelet = wname
    tree%boundary = bname
    tree%original_shape = shape(x)
    if (bname /= 'periodic') then
      call set_status(tree%status,waveslim_invalid_input, &
        '2D packet transform currently supports periodic boundaries')
      return
    end if
    if (levels < 1 .or. any(mod(shape(x),2**levels) /= 0)) then
      call set_status(tree%status,waveslim_invalid_level, &
        '2D dimensions must be divisible by 2**levels')
      return
    end if
    block
      use waveslim_types, only : wavelet_filter_type
      type(wavelet_filter_type) :: filter
      filter = wave_filter(wname,status)
      if (.not. status%ok()) then
        tree%status = status
        return
      end if
      allocate(tree%level(0:levels))
      allocate(tree%level(0)%node(1))
      tree%level(0)%node(1)%values = x
      do j = 1, levels
        side = 2**j
        allocate(tree%level(j)%node(4**j))
        do parent = 0, 4**(j-1)-1
          ix = parent/2**(j-1)
          iy = modulo(parent,2**(j-1))
          if (mod(ix,2) == 0) then
            xlow = 2*ix
            xhigh = 2*ix+1
          else
            xlow = 2*ix+1
            xhigh = 2*ix
          end if
          if (mod(iy,2) == 0) then
            ylow = 2*iy
            yhigh = 2*iy+1
          else
            ylow = 2*iy+1
            yhigh = 2*iy
          end if
          call dwt2_step(tree%level(j-1)%node(parent+1)%values, &
            filter%hpf,filter%lpf,ll,lh,hl,hh)
          tree%level(j)%node(xlow*side+ylow+1)%values = ll
          tree%level(j)%node(xlow*side+yhigh+1)%values = lh
          tree%level(j)%node(xhigh*side+ylow+1)%values = hl
          tree%level(j)%node(xhigh*side+yhigh+1)%values = hh
        end do
      end do
    end block
    call clear_status(tree%status)
  end function dwpt_2d

  function idwpt_2d(tree) result(x)
    type(packet_transform_2d), intent(in) :: tree
    real(dp), allocatable :: x(:,:)
    type(real_matrix), allocatable :: nodes(:), parents(:)
    integer :: j, parent, side, ix, iy
    integer :: xlow, xhigh, ylow, yhigh
    type(status_type) :: status

    if (.not. allocated(tree%level)) then
      allocate(x(0,0))
      return
    end if
    block
      use waveslim_types, only : wavelet_filter_type
      type(wavelet_filter_type) :: filter
      filter = wave_filter(tree%wavelet,status)
      if (.not. status%ok()) then
        allocate(x(0,0))
        return
      end if
      j = tree%levels()
      allocate(nodes(4**j))
      do parent = 1, size(nodes)
        nodes(parent)%values = tree%level(j)%node(parent)%values
      end do
      do j = tree%levels(), 1, -1
        side = 2**j
        allocate(parents(4**(j-1)))
        do parent = 0, 4**(j-1)-1
          ix = parent/2**(j-1)
          iy = modulo(parent,2**(j-1))
          if (mod(ix,2) == 0) then
            xlow = 2*ix
            xhigh = 2*ix+1
          else
            xlow = 2*ix+1
            xhigh = 2*ix
          end if
          if (mod(iy,2) == 0) then
            ylow = 2*iy
            yhigh = 2*iy+1
          else
            ylow = 2*iy+1
            yhigh = 2*iy
          end if
          call idwt2_step( &
            nodes(xlow*side+ylow+1)%values, &
            nodes(xlow*side+yhigh+1)%values, &
            nodes(xhigh*side+ylow+1)%values, &
            nodes(xhigh*side+yhigh+1)%values, &
            filter%hpf,filter%lpf,parents(parent+1)%values)
        end do
        call move_alloc(parents,nodes)
      end do
      x = nodes(1)%values
    end block
  end function idwpt_2d

  function dwpt_boot(y, wf, n_levels, seed) result(sample)
    real(dp), intent(in) :: y(:)
    character(len=*), intent(in), optional :: wf
    integer, intent(in), optional :: n_levels
    integer(i8), intent(in), optional :: seed
    real(dp), allocatable :: sample(:)
    character(len=16) :: wname
    integer :: levels, node, i, n, pick
    real(dp) :: u
    type(packet_transform) :: tree
    real(dp), allocatable :: source(:)

    wname = 'la8'
    if (present(wf)) wname = trim(wf)
    levels = max(1,int(log(real(size(y),dp))/log(2.0_dp))-1)
    if (present(n_levels)) levels = n_levels
    tree = dwpt(y,wname,levels)
    if (.not. tree%status%ok()) then
      allocate(sample(0))
      return
    end if
    if (present(seed)) call seed_rng(seed)
    do node = 1, size(tree%level(levels)%node)
      source = tree%level(levels)%node(node)%values
      n = size(source)
      do i = 1, n
        call random_number(u)
        pick = min(n,1+int(u*real(n,dp)))
        tree%level(levels)%node(node)%values(i) = source(pick)
      end do
    end do
    sample = idwpt(tree)
  end function dwpt_boot

  subroutine brick_wall_2d(wt, method)
    type(wavelet_transform_2d), intent(inout) :: wt
    character(len=*), intent(in), optional :: method
    character(len=16) :: transform_method
    integer :: j, filter_length, count
    real(dp) :: nan
    type(status_type) :: status
    block
      use waveslim_types, only : wavelet_filter_type
      type(wavelet_filter_type) :: filter
      filter = wave_filter(wt%wavelet,status)
      if (.not. status%ok()) return
      filter_length = filter%length()
    end block
    transform_method = wt%method
    if (present(method)) transform_method = trim(method)
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    do j = 1, size(wt%level)
      if (transform_method == 'dwt') then
        count = ceiling(real(filter_length-2,dp)*(1.0_dp-1.0_dp/ &
          real(2**j,dp)))
      else
        count = (2**j-1)*(filter_length-1)
      end if
      count = min(count,size(wt%level(j)%lh,1),size(wt%level(j)%lh,2))
      if (count <= 0) cycle
      call blank_edges(wt%level(j)%lh,count,nan)
      call blank_edges(wt%level(j)%hl,count,nan)
      call blank_edges(wt%level(j)%hh,count,nan)
    end do
  end subroutine brick_wall_2d

  subroutine blank_edges(x, count, value)
    real(dp), intent(inout) :: x(:,:)
    integer, intent(in) :: count
    real(dp), intent(in) :: value
    x(1:count,:) = value
    x(:,1:count) = value
  end subroutine blank_edges

  subroutine phase_shift_packet(tree)
    type(packet_transform), intent(inout) :: tree
    integer :: j, node, shift_amount
    type(wavelet_transform) :: temporary

    do j = 1, tree%levels()
      do node = 1, size(tree%level(j)%node)
        allocate(temporary%detail(1))
        temporary%detail(1)%values = tree%level(j)%node(node)%values
        temporary%smooth = tree%level(j)%node(node)%values
        temporary%wavelet = tree%wavelet
        temporary%method = tree%method
        temporary%original_length = size(tree%level(j)%node(node)%values)
        shift_amount = modulo(node-1,max(1,size(tree%level(j)%node(node)%values)))
        call phase_shift(temporary)
        tree%level(j)%node(node)%values = cshift_simple( &
          temporary%detail(1)%values,shift_amount)
        deallocate(temporary%detail,temporary%smooth)
      end do
    end do
  end subroutine phase_shift_packet

  function cshift_simple(x, amount) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: amount
    real(dp), allocatable :: y(:)
    integer :: i
    allocate(y(size(x)))
    do i = 1, size(x)
      y(i) = x(modulo(i-1+amount,size(x))+1)
    end do
  end function cshift_simple

  function testing_hov(x, wf, n_levels, min_coefficients) result(changes)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in), optional :: wf
    integer, intent(in), optional :: n_levels, min_coefficients
    type(variance_change_result) :: changes
    character(len=16) :: wname
    integer :: levels, minimum, j, n
    type(wavelet_transform) :: decimated, maximal
    real(dp), allocatable :: dcoef(:), mcoef(:)

    wname = 'la8'
    if (present(wf)) wname = trim(wf)
    levels = min(4,int(log(real(size(x),dp))/log(2.0_dp)))
    if (present(n_levels)) levels = n_levels
    minimum = 128
    if (present(min_coefficients)) minimum = min_coefficients
    decimated = dwt(x,wname,levels)
    maximal = modwt(x,wname,levels)
    if (.not. decimated%status%ok() .or. .not. maximal%status%ok()) then
      call set_status(changes%status,waveslim_invalid_input, &
        'wavelet transform failed in variance homogeneity test')
      return
    end if
    call brick_wall(decimated,wname,'dwt')
    call brick_wall(maximal,wname,'modwt')
    do j = 1, levels
      dcoef = pack(decimated%detail(j)%values, &
        ieee_is_finite(decimated%detail(j)%values))
      mcoef = pack(maximal%detail(j)%values, &
        ieee_is_finite(maximal%detail(j)%values))
      n = min(size(dcoef),size(mcoef))
      if (n > 0) call mult_loc(dcoef,mcoef,j,minimum,1,1,changes)
    end do
    call clear_status(changes%status)
  end function testing_hov

  recursive subroutine mult_loc(dwt_coefficients, modwt_coefficients, level, &
      min_coefficients, left_dwt, left_modwt, changes)
    real(dp), intent(in) :: dwt_coefficients(:), modwt_coefficients(:)
    integer, intent(in) :: level, min_coefficients, left_dwt, left_modwt
    type(variance_change_result), intent(inout) :: changes
    real(dp), allocatable :: statistic_dwt(:), statistic_modwt(:)
    real(dp) :: total, cumulative, critical, maximum
    integer :: i, n_dwt, n_modwt, location_dwt, location_modwt

    n_dwt = size(dwt_coefficients)
    n_modwt = size(modwt_coefficients)
    if (n_dwt <= min_coefficients .or. n_dwt < 3 .or. n_modwt < 3) return
    allocate(statistic_dwt(n_dwt),statistic_modwt(n_modwt))
    total = sum(dwt_coefficients*dwt_coefficients)
    if (total <= tiny(1.0_dp)) return
    cumulative = 0.0_dp
    do i = 1, n_dwt
      cumulative = cumulative+dwt_coefficients(i)**2
      statistic_dwt(i) = max( &
        real(i,dp)/real(n_dwt-1,dp)-cumulative/total, &
        cumulative/total-real(i-1,dp)/real(n_dwt-1,dp))
    end do
    maximum = maxval(statistic_dwt)
    location_dwt = maxloc(statistic_dwt,dim=1)
    total = sum(modwt_coefficients*modwt_coefficients)
    if (total <= tiny(1.0_dp)) return
    cumulative = 0.0_dp
    do i = 1, n_modwt
      cumulative = cumulative+modwt_coefficients(i)**2
      statistic_modwt(i) = max( &
        real(i,dp)/real(n_modwt-1,dp)-cumulative/total, &
        cumulative/total-real(i-1,dp)/real(n_modwt-1,dp))
    end do
    location_modwt = maxloc(statistic_modwt,dim=1)
    critical = sqrt(2.0_dp)*1.358_dp/sqrt(real(n_dwt,dp))
    if (maximum <= critical) return
    call append_change(changes,level,maximum,left_dwt+location_dwt-1, &
      left_modwt+location_modwt-1)
    if (location_dwt > 2 .and. location_modwt > 2) then
      call mult_loc(dwt_coefficients(:location_dwt-1), &
        modwt_coefficients(:location_modwt-1),level,min_coefficients, &
        left_dwt,left_modwt,changes)
    end if
    if (location_dwt < n_dwt-1 .and. location_modwt < n_modwt-1) then
      call mult_loc(dwt_coefficients(location_dwt+1:), &
        modwt_coefficients(location_modwt+1:),level,min_coefficients, &
        left_dwt+location_dwt,left_modwt+location_modwt,changes)
    end if
  end subroutine mult_loc

  subroutine append_change(changes, level, statistic, dwt_location, &
      modwt_location)
    type(variance_change_result), intent(inout) :: changes
    integer, intent(in) :: level, dwt_location, modwt_location
    real(dp), intent(in) :: statistic
    integer, allocatable :: integer_work(:)
    real(dp), allocatable :: real_work(:)
    integer :: n
    n = changes%count()
    allocate(integer_work(n+1))
    if (n > 0) integer_work(:n) = changes%level
    integer_work(n+1) = level
    call move_alloc(integer_work,changes%level)
    allocate(real_work(n+1))
    if (n > 0) real_work(:n) = changes%statistic
    real_work(n+1) = statistic
    call move_alloc(real_work,changes%statistic)
    allocate(integer_work(n+1))
    if (n > 0) integer_work(:n) = changes%location_dwt
    integer_work(n+1) = dwt_location
    call move_alloc(integer_work,changes%location_dwt)
    allocate(integer_work(n+1))
    if (n > 0) integer_work(:n) = changes%location_modwt
    integer_work(n+1) = modwt_location
    call move_alloc(integer_work,changes%location_modwt)
  end subroutine append_change

end module waveslim_extended
