module dice_design_stochastic
  use iso_fortran_env, only : int64
  use dice_design_kinds, only : dp
  use dice_design_rng, only : rng_state
  use dice_design_utils, only : determinant, euclidean_distance
  implicit none
  private

  type, public :: dmax_result
    real(dp), allocatable :: design_init(:, :)
    real(dp), allocatable :: design(:, :)
    real(dp) :: range = 0.0_dp
    real(dp) :: det_init = 0.0_dp
    real(dp) :: det_end = 0.0_dp
    integer :: niter = 0
  end type dmax_result

  type, public :: strauss_result
    real(dp), allocatable :: design_init(:, :)
    real(dp), allocatable :: design(:, :)
    real(dp) :: radius = 0.0_dp
    real(dp) :: alpha = 0.5_dp
    real(dp) :: repulsion = 0.001_dp
    real(dp) :: repulsion1d = 0.0001_dp
    integer :: nmc = 0
    logical :: constraints1d = .false.
  end type strauss_result

  type, public :: wsp_result
    real(dp), allocatable :: design(:, :)
    real(dp), allocatable :: residual_design(:, :)
    integer, allocatable :: selected_indices(:)
    real(dp) :: dmin = 0.0_dp
  end type wsp_result

  public :: dmax_design, strauss_design, wsp_design

contains

  pure function spherical_variogram(h, a) result(gamma)
    real(dp), intent(in) :: h, a
    real(dp) :: gamma, r
    if (h <= a) then
      r = h/a
      gamma = 1.5_dp*r - 0.5_dp*r**3
    else
      gamma = 1.0_dp
    end if
  end function spherical_variogram

  subroutine fill_correlation(p, range, c)
    real(dp), intent(in) :: p(:, :), range
    real(dp), intent(out) :: c(:, :)
    integer :: n, i, j
    real(dp) :: d
    n = size(p,1)
    c = 0.0_dp
    do i = 1, n
      c(i,i) = 1.0_dp
    end do
    do j = 2, n
      do i = 1, j-1
        d = euclidean_distance(p(i,:),p(j,:))
        if (d <= range) c(i,j) = 1.0_dp-spherical_variogram(d,range)
        c(j,i) = c(i,j)
      end do
    end do
  end subroutine fill_correlation

  subroutine dmax_design(n, dimension, range_in, niter_max, result, seed)
    integer, intent(in) :: n, dimension, niter_max
    real(dp), intent(in) :: range_in
    type(dmax_result), intent(out) :: result
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    real(dp), allocatable :: p(:, :), pu(:, :), c(:, :), cu(:, :), u(:)
    real(dp) :: range, dinit, du, dist
    integer :: i, j, iter, chosen

    if (n < dimension) error stop 'dmax_design: number of points is lower than dimension'
    if (n <= 0 .or. dimension <= 0 .or. niter_max < 1) error stop 'dmax_design: invalid dimensions'
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(p(n,dimension),pu(n,dimension),c(n,n),cu(n,n),u(dimension))
    do j = 1, dimension
      do i = 1, n
        p(i,j) = rng%uniform()
      end do
    end do
    allocate(result%design_init(n,dimension))
    result%design_init = p
    range = range_in
    if (abs(range) <= 16.0_dp*epsilon(1.0_dp)) range = sqrt(real(dimension,dp))/2.0_dp
    if (range <= 0.0_dp) error stop 'dmax_design: range must be nonnegative'
    call fill_correlation(p,range,c)
    dinit = determinant(c)
    result%det_init = dinit
    do iter = 2, niter_max
      chosen = rng%integer(1,n)
      do j = 1, dimension
        u(j) = rng%uniform()
      end do
      pu = p
      pu(chosen,:) = u
      cu = c
      cu(chosen,chosen) = 1.0_dp
      do i = 1, n
        if (i == chosen) cycle
        dist = euclidean_distance(u,pu(i,:))
        if (dist <= range) then
          cu(chosen,i) = 1.0_dp-spherical_variogram(dist,range)
        else
          cu(chosen,i) = 0.0_dp
        end if
        cu(i,chosen) = cu(chosen,i)
      end do
      du = determinant(cu)
      if (du >= dinit) then
        p = pu
        c = cu
        dinit = du
      end if
    end do
    allocate(result%design(n,dimension))
    result%design = p
    result%range = range
    result%det_end = dinit
    result%niter = niter_max
  end subroutine dmax_design

  subroutine strauss_design(n, dimension, rnd, result, alpha, repulsion, nmc, constraints1d, &
                            repulsion1d, seed)
    integer, intent(in) :: n, dimension
    real(dp), intent(in) :: rnd
    type(strauss_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha, repulsion, repulsion1d
    integer, intent(in), optional :: nmc
    logical, intent(in), optional :: constraints1d
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    real(dp), allocatable :: v(:, :), u(:)
    real(dp) :: a, gamma_nd, gamma_1d, beta_nd, r1d, p1d, pnd, prob
    real(dp) :: nxr, nur, px, pu, dx, du
    integer :: sweeps, sweep, k, chosen, dim, point, nx, nu
    logical :: use1d

    if (n <= 0 .or. dimension <= 0) error stop 'strauss_design: invalid dimensions'
    if (rnd <= 0.0_dp) error stop 'strauss_design: radius must be positive'
    a = 0.5_dp
    if (present(alpha)) a = alpha
    gamma_nd = 0.001_dp
    if (present(repulsion)) gamma_nd = repulsion
    gamma_1d = 0.0001_dp
    if (present(repulsion1d)) gamma_1d = repulsion1d
    sweeps = 1000
    if (present(nmc)) sweeps = nmc
    use1d = .false.
    if (present(constraints1d)) use1d = constraints1d
    if (gamma_nd <= 0.0_dp .or. gamma_nd > 1.0_dp) error stop 'strauss_design: repulsion must be in (0,1]'
    if (gamma_1d < 0.0_dp .or. gamma_1d > 1.0_dp) error stop 'strauss_design: repulsion1d must be in [0,1]'
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(v(n,dimension),u(dimension))
    do dim = 1, dimension
      do point = 1, n
        v(point,dim) = rng%uniform()
      end do
    end do
    allocate(result%design_init(n,dimension))
    result%design_init = v
    r1d = 0.75_dp/real(n,dp)
    if (abs(a) > 16.0_dp*epsilon(1.0_dp)) beta_nd = -log(gamma_nd)

    do sweep = 1, sweeps
      do k = 1, n
        chosen = rng%integer(1,n)
        p1d = 1.0_dp
        do dim = 1, dimension
          u(dim) = rng%uniform()
          if (use1d) then
            nx = 0
            nu = 0
            do point = 1, n
              if (point == chosen) cycle
              if (abs(v(point,dim)-v(chosen,dim)) < r1d) nx = nx+1
              if (abs(v(point,dim)-u(dim)) < r1d) nu = nu+1
            end do
            p1d = p1d*gamma_1d**max(0,nu-nx)
          end if
        end do

        if (abs(a) > 16.0_dp*epsilon(1.0_dp)) then
          px = 0.0_dp
          pu = 0.0_dp
          do point = 1, n
            if (point == chosen) cycle
            nxr = euclidean_distance(v(point,:),v(chosen,:))
            nur = euclidean_distance(v(point,:),u)
            if (nxr < rnd) px = px+(1.0_dp-nxr/rnd)**a
            if (nur < rnd) pu = pu+(1.0_dp-nur/rnd)**a
          end do
          pnd = min(1.0_dp,exp(-beta_nd*(pu-px)))
          prob = pnd*p1d
        else
          nx = 0
          nu = 0
          do point = 1, n
            if (point == chosen) cycle
            dx = sum((v(point,:)-v(chosen,:))**2)
            du = sum((v(point,:)-u)**2)
            if (dx < rnd*rnd) nx = nx+1
            if (du < rnd*rnd) nu = nu+1
          end do
          prob = gamma_nd**max(0,nu-nx)*p1d
        end if
        if (rng%uniform() < prob) v(chosen,:) = u
      end do
    end do

    allocate(result%design(n,dimension))
    result%design = v
    result%radius = rnd
    result%alpha = a
    result%repulsion = gamma_nd
    result%repulsion1d = gamma_1d
    result%nmc = sweeps
    result%constraints1d = use1d
  end subroutine strauss_design

  subroutine wsp_design(initial_design, dmin, result, random_init, seed)
    real(dp), intent(in) :: initial_design(:, :)
    real(dp), intent(in) :: dmin
    type(wsp_result), intent(out) :: result
    logical, intent(in), optional :: random_init
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    logical, allocatable :: active(:), selected(:)
    integer, allocatable :: idx(:), ridx(:)
    real(dp), allocatable :: anchor(:)
    real(dp) :: dist, best
    integer :: n, d, i, base, next_base, ns, nr
    logical :: random_anchor

    n = size(initial_design,1)
    d = size(initial_design,2)
    if (n == 0 .or. d == 0) error stop 'wsp_design: empty design'
    if (dmin < 0.0_dp) error stop 'wsp_design: dmin must be nonnegative'
    random_anchor = .false.
    if (present(random_init)) random_anchor = random_init
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(anchor(d),active(n),selected(n))
    if (random_anchor) then
      do i = 1, d
        anchor(i) = rng%uniform()
      end do
    else
      anchor = 0.5_dp
    end if
    best = huge(1.0_dp)
    base = 1
    do i = 1, n
      dist = euclidean_distance(initial_design(i,:),anchor)
      if (dist < best) then
        best = dist
        base = i
      end if
    end do
    active = .true.
    selected = .false.

    do
      if (.not. active(base)) exit
      selected(base) = .true.
      active(base) = .false.
      do i = 1, n
        if (.not. active(i)) cycle
        if (euclidean_distance(initial_design(base,:),initial_design(i,:)) < dmin) active(i) = .false.
      end do
      if (.not. any(active)) exit
      best = huge(1.0_dp)
      next_base = 0
      do i = 1, n
        if (.not. active(i)) cycle
        dist = euclidean_distance(initial_design(base,:),initial_design(i,:))
        if (dist < best) then
          best = dist
          next_base = i
        end if
      end do
      if (next_base == 0) exit
      base = next_base
    end do

    ns = count(selected)
    nr = n-ns
    allocate(idx(ns),ridx(nr),result%design(ns,d),result%residual_design(nr,d),result%selected_indices(ns))
    ns = 0
    nr = 0
    do i = 1, n
      if (selected(i)) then
        ns = ns+1
        idx(ns) = i
      else
        nr = nr+1
        ridx(nr) = i
      end if
    end do
    do i = 1, size(idx)
      result%design(i,:) = initial_design(idx(i),:)
    end do
    do i = 1, size(ridx)
      result%residual_design(i,:) = initial_design(ridx(i),:)
    end do
    result%selected_indices = idx
    result%dmin = dmin
  end subroutine wsp_design

end module dice_design_stochastic
