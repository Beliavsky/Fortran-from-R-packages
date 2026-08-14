module dice_design_classic
  use iso_fortran_env, only : int64
  use dice_design_kinds, only : dp
  use dice_design_rng, only : rng_state
  use dice_design_utils, only : empirical_cdf_column, empirical_quantile
  use dice_design_nolh_data, only : load_nolh_base, load_nolhdr
  implicit none
  private

  public :: fact_design, lhs_design, olh_design, nolh_design, nolhdr_design
  public :: runif_faure, faureprime_design, scale_design, unscale_design
  public :: xdrdn_transform

contains

  pure logical function near_value(x, y) result(ok)
    real(dp), intent(in) :: x, y
    real(dp) :: scale
    scale = max(1.0_dp, abs(x), abs(y))
    ok = abs(x-y) <= 16.0_dp*epsilon(1.0_dp)*scale
  end function near_value

  subroutine fact_design(dimension, levels, design)
    integer, intent(in) :: dimension
    integer, intent(in) :: levels(:)
    real(dp), allocatable, intent(out) :: design(:, :)
    integer, allocatable :: lev(:)
    integer :: n, i, j, stride, idx

    if (dimension <= 0) error stop 'fact_design: dimension must be positive'
    allocate(lev(dimension))
    if (size(levels) == 1) then
      lev = levels(1)
    else if (size(levels) == dimension) then
      lev = levels
    else
      error stop 'fact_design: levels must have length 1 or dimension'
    end if
    if (any(lev < 1)) error stop 'fact_design: levels must be positive'
    n = product(lev)
    allocate(design(n, dimension))
    do j = 1, dimension
      stride = 1
      if (j < dimension) stride = product(lev(j+1:dimension))
      do i = 1, n
        idx = modulo((i-1)/stride, lev(j))
        if (lev(j) == 1) then
          design(i,j) = 0.0_dp
        else
          design(i,j) = real(idx,dp)/real(lev(j)-1,dp)
        end if
      end do
    end do
  end subroutine fact_design

  subroutine lhs_design(n, dimension, design, randomized, seed)
    integer, intent(in) :: n, dimension
    real(dp), allocatable, intent(out) :: design(:, :)
    logical, intent(in), optional :: randomized
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    integer, allocatable :: idx(:)
    logical :: random_points
    integer :: i, j
    real(dp) :: ran

    if (n <= 0 .or. dimension <= 0) error stop 'lhs_design: n and dimension must be positive'
    random_points = .true.
    if (present(randomized)) random_points = randomized
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(design(n,dimension), idx(n))
    do i = 1, n
      idx(i) = i
    end do
    do j = 1, dimension
      do i = 1, n
        idx(i) = i
      end do
      call rng%shuffle(idx)
      do i = 1, n
        if (random_points) then
          ran = rng%uniform()
        else
          ran = 0.5_dp
        end if
        design(i,j) = (real(idx(i),dp) - ran)/real(n,dp)
      end do
    end do
  end subroutine lhs_design

  recursive subroutine make_m(power, m)
    integer, intent(in) :: power
    integer, allocatable, intent(out) :: m(:, :)
    integer, allocatable :: old(:, :)
    integer :: nr, nc

    if (power == 1) then
      allocate(m(2,2))
      m = reshape([1,2,2,1],[2,2])
      return
    end if
    call make_m(power-1, old)
    nr = size(old,1)
    nc = size(old,2)
    allocate(m(2*nr,2*nc))
    m(1:nr,1:nc) = old
    m(1:nr,nc+1:2*nc) = old + nc
    m(nr+1:2*nr,1:nc) = old + nc
    m(nr+1:2*nr,nc+1:2*nc) = old
  end subroutine make_m

  recursive subroutine make_s(power, s)
    integer, intent(in) :: power
    integer, allocatable, intent(out) :: s(:, :)
    integer, allocatable :: old(:, :), p(:, :), q(:, :)
    integer :: nr, nc, half

    if (power == 1) then
      allocate(s(2,2))
      s = reshape([1,1,1,-1],[2,2])
      return
    end if
    call make_s(power-1, old)
    nr = size(old,1)
    nc = size(old,2)
    half = nr/2
    allocate(p(half,nc), q(half,nc))
    p = old(1:half,:)
    q = old(half+1:nr,:)
    allocate(s(4*half,2*nc))
    s(1:half,1:nc) = p
    s(1:half,nc+1:2*nc) = p
    s(half+1:2*half,1:nc) = q
    s(half+1:2*half,nc+1:2*nc) = -q
    s(2*half+1:3*half,1:nc) = p
    s(2*half+1:3*half,nc+1:2*nc) = -p
    s(3*half+1:4*half,1:nc) = q
    s(3*half+1:4*half,nc+1:2*nc) = q
  end subroutine make_s

  subroutine olh_design(dimension, design, target_range)
    integer, intent(in) :: dimension
    real(dp), allocatable, intent(out) :: design(:, :)
    real(dp), intent(in), optional :: target_range(2)
    integer, allocatable :: m(:, :), s(:, :), t(:, :), full(:, :)
    real(dp) :: r(2)
    integer :: p, nr, nc

    if (dimension <= 0) error stop 'olh_design: dimension must be positive'
    p = max(1, ceiling(log(real(dimension,dp))/log(2.0_dp)))
    call make_m(p,m)
    call make_s(p,s)
    t = m*s
    nr = size(t,1)
    nc = size(t,2)
    allocate(full(2*nr+1,nc))
    full(1:nr,:) = t
    full(nr+1,:) = 0
    full(nr+2:2*nr+1,:) = -t
    allocate(design(size(full,1),dimension))
    design = real(full(:,1:dimension),dp)
    r = [0.0_dp,1.0_dp]
    if (present(target_range)) r = target_range
    if (near_value(r(1),1.0_dp) .and. near_value(r(2),1.0_dp)) then
      continue
    else if (near_value(r(1),0.0_dp) .and. near_value(r(2),0.0_dp)) then
      design = design + real(size(full,1)-1,dp)/2.0_dp
    else
      design = (design/real(size(full,1)-1,dp)+0.5_dp)*(r(2)-r(1))+r(1)
    end if
  end subroutine olh_design

  subroutine nolh_design(dimension, design, target_range)
    integer, intent(in) :: dimension
    real(dp), allocatable, intent(out) :: design(:, :)
    real(dp), intent(in), optional :: target_range(2)
    real(dp), allocatable :: base(:, :)
    real(dp) :: r(2)
    logical :: ok
    integer :: group

    if (dimension < 2 .or. dimension > 29) error stop 'nolh_design: dimension must be in 2:29'
    if (dimension <= 7) then
      group = 1
    else if (dimension <= 11) then
      group = 2
    else if (dimension <= 16) then
      group = 3
    else if (dimension <= 22) then
      group = 4
    else
      group = 5
    end if
    call load_nolh_base(group,base,ok)
    if (.not. ok) error stop 'nolh_design: embedded design unavailable'
    allocate(design(size(base,1),dimension))
    design = base(:,1:dimension)
    r = [0.0_dp,1.0_dp]
    if (present(target_range)) r = target_range
    call apply_nolh_range(design,r)
  end subroutine nolh_design

  subroutine nolhdr_design(dimension, design, target_range)
    integer, intent(in) :: dimension
    real(dp), allocatable, intent(out) :: design(:, :)
    real(dp), intent(in), optional :: target_range(2)
    real(dp), allocatable :: base(:, :)
    real(dp) :: r(2)
    logical :: ok

    if (dimension < 2 .or. dimension > 29) error stop 'nolhdr_design: dimension must be in 2:29'
    if (dimension <= 7) then
      call load_nolh_base(1,base,ok)
    else
      call load_nolhdr(dimension,base,ok)
    end if
    if (.not. ok) error stop 'nolhdr_design: embedded design unavailable'
    allocate(design(size(base,1),dimension))
    design = base(:,1:dimension)
    r = [0.0_dp,1.0_dp]
    if (present(target_range)) r = target_range
    call apply_nolh_range(design,r)
  end subroutine nolhdr_design

  subroutine apply_nolh_range(design,r)
    real(dp), intent(inout) :: design(:, :)
    real(dp), intent(in) :: r(2)
    integer :: n
    n = size(design,1)
    if (near_value(r(1),1.0_dp) .and. near_value(r(2),1.0_dp)) then
      return
    else if (near_value(r(1),0.0_dp) .and. near_value(r(2),0.0_dp)) then
      design = design + real(n-1,dp)/2.0_dp
    else
      design = (design/real(n-1,dp)+0.5_dp)*(r(2)-r(1))+r(1)
    end if
  end subroutine apply_nolh_range

  subroutine runif_faure(n, dimension, design)
    integer, intent(in) :: n, dimension
    real(dp), allocatable, intent(out) :: design(:, :)
    integer, parameter :: primes(168) = [ &
      2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101, &
      103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199, &
      211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317, &
      331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443, &
      449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577, &
      587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701, &
      709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839, &
      853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997]
    integer(int64), allocatable :: c(:, :), nr(:), next(:)
    real(dp), allocatable :: dg(:)
    integer :: r, m, ip, i, j, k
    integer(int64) :: number, quo, res

    if (n <= 0 .or. dimension <= 0) error stop 'runif_faure: n and dimension must be positive'
    r = 0
    do ip = 1, size(primes)
      if (primes(ip) >= dimension) then
        r = primes(ip)
        exit
      end if
    end do
    if (r == 0) error stop 'runif_faure: dimension exceeds supported prime table'
    m = max(1, ceiling(log(real(n,dp))/log(real(r,dp))))
    allocate(c(m,m),nr(m),next(m),dg(m),design(n,dimension))
    c = 0_int64
    c(1,:) = 1_int64
    do j = 2, m
      c(j,j) = 1_int64
      if (j-1 >= 2) then
        do i = 2, j-1
          c(i,j) = modulo(c(i,j-1)+c(i-1,j-1),int(r,int64))
        end do
      end if
    end do
    do k = 1, m
      dg(k) = real(r,dp)**(-k)
    end do
    do i = 1, n
      nr = 0_int64
      number = int(i,int64)
      do k = 1, m
        quo = number/int(r,int64)
        res = number - quo*int(r,int64)
        number = quo
        nr(k) = res
      end do
      design(i,1) = sum(dg*real(nr,dp))
      do j = 2, dimension
        next = 0_int64
        do k = 1, m
          next(k) = modulo(sum(c(k,:)*nr),int(r,int64))
        end do
        nr = next
        design(i,j) = sum(dg*real(nr,dp))
      end do
    end do
  end subroutine runif_faure

  subroutine faureprime_design(dimension, u, design, prime, target_range)
    integer, intent(in) :: dimension, u
    real(dp), allocatable, intent(out) :: design(:, :)
    integer, intent(out), optional :: prime
    real(dp), intent(in), optional :: target_range(2)
    integer, parameter :: ptab(46) = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67, &
      71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199]
    real(dp), allocatable :: f(:, :)
    real(dp) :: r(2)
    integer :: p, i, pu, n

    if (u < 2) error stop 'faureprime_design: u must be at least 2'
    if (dimension < 3 .or. dimension > 199) error stop 'faureprime_design: dimension must be in 3:199'
    p = 0
    do i = 1, size(ptab)
      if (ptab(i) >= dimension) then
        p = ptab(i)
        exit
      end if
    end do
    if (int(p,int64)**(u+1) > 8000000_int64) error stop 'faureprime_design: prime^(u+1) exceeds 200^3'
    pu = p**u
    n = pu - 1
    call runif_faure(n,dimension,f)
    allocate(design(n,dimension))
    design = real(pu,dp)*f
    r = [0.0_dp,-1.0_dp]
    if (present(target_range)) r = target_range
    if (near_value(r(1),0.0_dp) .and. near_value(r(2),0.0_dp)) then
      continue
    else if (near_value(r(1),1.0_dp) .and. near_value(r(2),1.0_dp)) then
      design = 2.0_dp*design-real(n+1,dp)
    else if (near_value(r(1),0.0_dp) .and. near_value(r(2),1.0_dp)) then
      design = (design-1.0_dp)/real(n-1,dp)
    else if (near_value(r(1),0.0_dp) .and. near_value(r(2),-1.0_dp)) then
      design = design/real(n+1,dp)
    else if (near_value(r(1),-1.0_dp) .and. near_value(r(2),-1.0_dp)) then
      design = (2.0_dp*design-real(n+1,dp))/real(n+1,dp)
    else if (near_value(r(1),-1.0_dp) .and. near_value(r(2),1.0_dp)) then
      design = (2.0_dp*design-real(n+1,dp))/real(n-1,dp)
    else
      if (r(2) <= r(1)) error stop 'faureprime_design: invalid target range'
      design = (design-1.0_dp)/real(n-1,dp)*(r(2)-r(1))+r(1)
    end if
    if (present(prime)) prime = p
  end subroutine faureprime_design

  subroutine scale_design(design, scaled, lower, upper, uniformize)
    real(dp), intent(in) :: design(:, :)
    real(dp), allocatable, intent(out) :: scaled(:, :)
    real(dp), intent(in), optional :: lower(:), upper(:)
    logical, intent(in), optional :: uniformize
    logical :: unif
    real(dp), allocatable :: lo(:), hi(:)
    integer :: n, d, j

    n = size(design,1)
    d = size(design,2)
    allocate(scaled(n,d),lo(d),hi(d))
    unif = .false.
    if (present(uniformize)) unif = uniformize
    if (unif) then
      do j = 1, d
        call empirical_cdf_column(design(:,j),scaled(:,j))
      end do
      return
    end if
    do j = 1, d
      lo(j) = minval(design(:,j))
      hi(j) = maxval(design(:,j))
    end do
    if (present(lower)) then
      if (size(lower) /= d) error stop 'scale_design: lower size mismatch'
      lo = lower
    end if
    if (present(upper)) then
      if (size(upper) /= d) error stop 'scale_design: upper size mismatch'
      hi = upper
    end if
    do j = 1, d
      if (near_value(hi(j),lo(j))) then
        scaled(:,j) = 0.0_dp
      else
        scaled(:,j) = (design(:,j)-lo(j))/(hi(j)-lo(j))
      end if
    end do
  end subroutine scale_design

  subroutine unscale_design(scaled, design, lower, upper, uniformize, initial_design)
    real(dp), intent(in) :: scaled(:, :)
    real(dp), allocatable, intent(out) :: design(:, :)
    real(dp), intent(in), optional :: lower(:), upper(:)
    logical, intent(in), optional :: uniformize
    real(dp), intent(in), optional :: initial_design(:, :)
    logical :: unif
    integer :: n, d, i, j

    n = size(scaled,1)
    d = size(scaled,2)
    allocate(design(n,d))
    unif = .false.
    if (present(uniformize)) unif = uniformize
    if (unif) then
      if (.not. present(initial_design)) error stop 'unscale_design: initial_design is required for uniformize'
      if (size(initial_design,2) /= d) error stop 'unscale_design: initial_design dimension mismatch'
      do j = 1, d
        do i = 1, n
          design(i,j) = empirical_quantile(initial_design(:,j),scaled(i,j))
        end do
      end do
    else
      if (.not. present(lower) .or. .not. present(upper)) error stop 'unscale_design: lower and upper are required'
      if (size(lower) /= d .or. size(upper) /= d) error stop 'unscale_design: bound size mismatch'
      do j = 1, d
        design(:,j) = scaled(:,j)*(upper(j)-lower(j))+lower(j)
      end do
    end if
  end subroutine unscale_design

  subroutine xdrdn_transform(input, output, digits, target_range)
    real(dp), intent(in) :: input(:, :)
    real(dp), allocatable, intent(out) :: output(:, :)
    integer, intent(in), optional :: digits
    real(dp), intent(in), optional :: target_range(2)
    real(dp) :: lo, hi, fac
    integer :: dg

    allocate(output(size(input,1),size(input,2)))
    output = input
    if (present(target_range)) then
      lo = minval(input)
      hi = maxval(input)
      if (hi > lo) output = (input-lo)/(hi-lo)*(target_range(2)-target_range(1))+target_range(1)
    end if
    if (present(digits)) then
      dg = digits
      fac = 10.0_dp**dg
      output = anint(output*fac)/fac
    end if
  end subroutine xdrdn_transform

end module dice_design_classic
